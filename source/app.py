import argparse
import asyncio
import collections
import httpx
import logging
import random
import uvicorn

from crypto import *
from fastapi import FastAPI, Request

# Server for listening to our peers in the cluster (and ourselves).
peer = FastAPI()

# Server for listening to clients. Proxies etcd and transparently handles the
# encryption/decryption of values sent/received by the client.
proxy = FastAPI()

# Using the peer server state as the state for our entire app.
state = peer.state

# Piggybacking off of Uvicorn's logger.
logger = logging.getLogger("uvicorn.error")

@peer.post("/status")
async def status():
    """
    Return a status message indicating readiness to service requests.
    """
    return {"message": "Peer server is ready."}

@peer.post("/keygen/{id}")
async def keygen(id: int):
    """
    Return a share of the server's secret key and the server's public key.
    """
    xⱼ = id
    k  = state.secret_key
    Q  = state.public_key
    kⱼ = k[xⱼ]

    return {"secretKeyShare": kⱼ, "publicKey": Q.to_dict()}

@peer.post("/decrypt")
async def decrypt(request: Request):
    """
    Return the server's share of the decryption of a message.
    """
    body = await request.json()

    kᵢ = state.joint_secret_key_share
    R  = Point.from_dict(body)
    Sᵢ = kᵢ * R

    return Sᵢ.to_dict()

@proxy.get("/status")
async def status():
    """
    Return a status message indicating readiness to service requests.
    """
    return {"message": "Proxy server is ready."}

@proxy.put("/keys/{key}")
async def write(key: str, request: Request):
    """
    Proxy PUT requests to the upstream etcd API and perform encryption.
    """
    body = await request.json()

    # Elliptic-curve cryptography parameters.
    E = state.curve
    q = E.order
    G = E.generator

    # Generate an ephemeral key pair for this message.
    d = random.randint(1, q)
    R = d * G

    # Perform the encryption phase of the ECIES protocol.
    Q = state.joint_public_key
    M = body["value"]
    S = d * Q
    K = symmetric_derive_key(S)
    C = symmetric_encrypt(M, K)

    # Replace the message with the public key and ciphertext.
    body["value"] = json.dumps({"publicKey": R.to_dict(), "ciphertext": base64.urlsafe_b64encode(C).decode("utf-8")})

    # Carry on with the write request to etcd.
    async with httpx.AsyncClient() as client:
        response = await client.put(f"{state.etcd}/v2/keys/{key}", data=body)

    return response.json()

@proxy.get("/keys/{key}")
async def read(key: str):
    """
    Proxy GET requests to the upstream etcd API and perform decryption.
    """
    async with httpx.AsyncClient() as client:
        response = await client.get(f"{state.etcd}/v2/keys/{key}")
        body = response.json()

        # Extract the public key and ciphertext from the etcd value.
        value = json.loads(body["node"]["value"])
        R = Point.from_dict(value["publicKey"])
        C = base64.urlsafe_b64decode(value["ciphertext"].encode("utf-8"))

        # We need to collect decryption shares from our peers.
        t = state.threshold
        S = {}

        # Put the peers in a queue and visit them in a round-robin manner.
        queue = collections.deque(enumerate(state.cluster, 1))

        while len(S) != t:
            xⱼ, host = queue.popleft()
            try:
                # Timeout is more strict here because the client is waiting.
                response = await client.post(f"{host}/decrypt", json=R.to_dict(), timeout=2)
            except:
                # Peer cannot be reached. Put them at the end of the queue and
                # come back later if we still require shares.
                queue.append((xⱼ, host))
                logger.info(f"Peer {host} is unresponsive, moving on.")
            else:
                # Record the decryption share. Do not put the peer back in the
                # queue, we don't want to contact them again.
                S[xⱼ] = Point.from_dict(response.json())

    # Elliptic-curve cryptography parameters.
    E = state.curve
    q = E.order
    I = E.identity

    # Interpolate the decryption shares and reveal the message.
    S = interpolate(S, q, I)
    K = symmetric_derive_key(S)
    M = symmetric_decrypt(C, K)

    # Replace the etcd value with the message.
    body["node"]["value"] = M

    return body

async def run_peer():
    """
    Run the peer server API.
    """
    config = uvicorn.Config(peer, host="0.0.0.0", port=2380, log_level="info")
    server = uvicorn.Server(config)
    await server.serve()

async def run_proxy():
    """
    Run the proxy server API.
    """
    config = uvicorn.Config(proxy, host="0.0.0.0", port=2379, log_level="info")
    server = uvicorn.Server(config)
    await server.serve()

async def run():
    """
    Run the entire program, with initialization.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--id", type=int)
    parser.add_argument("--etcd", type=str)
    parser.add_argument("--cluster", type=str)
    parser.add_argument("--threshold", type=int)

    args = parser.parse_args()

    state.id = args.id
    state.etcd = args.etcd
    state.cluster = args.cluster.split(",")
    state.threshold = args.threshold

    # Threshold cryptography parameters.
    t = state.threshold
    n = len(state.cluster)
    assert t <= n

    # Elliptic-curve cryptography parameters.
    E = Curve("NIST P-256")
    q = E.order
    G = E.generator
    I = E.identity

    # Generate individual secret key shares and public key.
    k = random.randint(1, q)
    Q = k * G
    f = Polynomial.shamir(k, t, q)
    k = {xⱼ: f(xⱼ) for xⱼ in range(1, n + 1)}

    state.curve = E
    state.secret_key = k
    state.public_key = Q

    # Start the peer server.
    asyncio.create_task(run_peer())

    # Perform the distributed key generation algorithm.
    kᵢ = 0
    Q  = I

    async with httpx.AsyncClient() as client:
        for host in state.cluster:
            while True:
                try:
                    response = await client.post(f"{host}/keygen/{state.id}")
                except:
                    logger.info(f"Peer {host} is unresponsive, trying again.")
                    await asyncio.sleep(2)
                else:
                    break
            body = response.json()

            kᵢ += body["secretKeyShare"]
            Q  += Point.from_dict(body["publicKey"])

    state.joint_secret_key_share = kᵢ
    state.joint_public_key = Q

    # Start the proxy server. Clients can now make requests.
    await run_proxy()

if __name__ == "__main__":
    asyncio.run(run())

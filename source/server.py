import argparse
import asyncio
import httpx
import random
import uvicorn

from crypto import *
from fastapi import FastAPI, Request

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
etcd = "127.0.0.1:53667"
peers = ["127.0.0.1:9001", "127.0.0.1:9002", "127.0.0.1:9003"]
nodeid = None

# Elliptic-curve cryptography parameters.
E = Curve('P-256')
q = E.order
G = E.generator
I = E.identity

# Threshold cryptography parameters.
n = len(peers)
t = 2
assert t <= n

# ----------------------------------------------------------------------------
# Persistent State
# ----------------------------------------------------------------------------
# Individual secret key (gets reassigned below).
k_prime = random.randint(1, q)

# Individual public key.
Q_prime = k_prime * G

# Individual polynomial for Shamir's secret sharing.
f_prime = Polynomial.shamir(k_prime, t, q)

# Individual secret key shares.
k_prime = {xⱼ: f_prime(xⱼ) for xⱼ in range(1, n + 1)}

# Joint secret key share (calculated at startup).
kᵢ = 0

# Joint public key (calculated at startup).
Q = I

# ----------------------------------------------------------------------------
# Peer Server API
# ----------------------------------------------------------------------------
peer = FastAPI()

@peer.post("/status")
async def status():
    """
    Return a status message indicating readiness to service requests.
    """
    return {"message": "Peer server is ready"}

@peer.post("/keygen/{nodeid}")
async def keygen(nodeid: int):
    """
    Return a share of the server's secret key and the server's public key.
    """
    xⱼ = nodeid

    return {"secretKeyShare": k_prime[xⱼ], "publicKey": Q_prime.to_dict()}

@peer.post("/decrypt")
async def decrypt(request: Request):
    """
    Return the server's share of the decryption of a message.
    """
    body = await request.json()
    R = Point.from_dict(body)

    return {"decryptionShare": (kᵢ * R).to_dict()}

# ----------------------------------------------------------------------------
# Proxy Server API
# ----------------------------------------------------------------------------
proxy = FastAPI()

@proxy.get("/status")
async def status():
    """
    Return a status message indicating readiness to service requests.
    """
    return {"message": "Proxy server is ready"}

@proxy.put("/keys/{key}")
async def write(key: str, request: Request):
    """
    Proxy PUT requests to the upstream etcd API and perform encryption.
    """
    body = await request.json()

    # Generate an ephemeral key pair for this message.
    d = random.randint(1, q)
    R = d * G

    # Perform the encryption phase of the ECIES protocol.
    M = body["value"]
    S = d * Q
    K = symmetric_derive_key(S)
    C = symmetric_encrypt(M, K)

    # Replace the message with the public key and ciphertext.
    body["value"] = json.dumps({"publicKey": R.to_dict(), "ciphertext": base64.urlsafe_b64encode(C).decode("utf-8")})

    # Proxy the write request to etcd.
    async with httpx.AsyncClient() as client:
        response = await client.put(f"http://{etcd}/v2/keys/{key}", data=body)

    return response.json()

@proxy.get("/keys/{key}")
async def read(key: str):
    """
    Proxy GET requests to the upstream etcd API and perform decryption.
    """
    async with httpx.AsyncClient() as client:
        response = await client.get(f"http://{etcd}/v2/keys/{key}")
        body = response.json()

        # Extract the public key and ciphertext from the etcd value.
        value = json.loads(body["node"]["value"])
        R = Point.from_dict(value["publicKey"])
        C = base64.urlsafe_b64decode(value["ciphertext"].encode("utf-8"))

        # Collect decryption shares from peers.
        S = {}

        for xⱼ, peer in random.sample(sorted(enumerate(peers, 1)), t):
            response = await client.post(f"http://{peer}/decrypt", json=R.to_dict())
            temp = response.json()

            S[xⱼ] = Point.from_dict(temp["decryptionShare"])

    # Interpolate the decryption shares and reveal the message.
    S = interpolate(S, q, I)
    K = symmetric_derive_key(S)
    M = symmetric_decrypt(C, K)

    # Replace the etcd value with the message.
    body["node"]["value"] = M

    return body

# ----------------------------------------------------------------------------
# Main Execution
# ----------------------------------------------------------------------------
async def run_peer(port: int):
    """
    Run the peer server API.
    """
    config = uvicorn.Config(peer, port=port, log_level="info")
    server = uvicorn.Server(config)
    await server.serve()

async def run_proxy(port: int):
    """
    Run the proxy server API.
    """
    config = uvicorn.Config(proxy, port=port, log_level="info")
    server = uvicorn.Server(config)
    await server.serve()

async def run():
    """
    Run the entire program, with initialization.
    """
    global nodeid

    print("Starting peer server...")
    asyncio.create_task(run_peer(9000 + nodeid))

    print("Performing DKG algorithm...")
    global kᵢ, Q

    async with httpx.AsyncClient() as client:
        for peer in peers:
            while True:
                try:
                    response = await client.post(f"http://{peer}/keygen/{nodeid}")
                except:
                    print("Peer unresponsive, trying again...")
                    await asyncio.sleep(1)
                else:
                    break
            temp = response.json()

            kᵢ += temp["secretKeyShare"]
            Q  += Point.from_dict(temp["publicKey"])

    print("Starting proxy server...")
    await run_proxy(8000 + nodeid)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int)

    args = parser.parse_args()

    nodeid = args.index

    asyncio.run(run())

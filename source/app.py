import ssl
import argparse
import asyncio
import base64
import collections
import httpx
import json
import logging
import random
import uvicorn

from crypto import *
from fastapi import FastAPI, Request
from starlette.datastructures import State

# Piggybacking off of Uvicorn's logger.
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("uvicorn.error")

# Server for listening to our peers in the cluster (and ourselves).
peer = FastAPI()

# Using the peer server state as the state for our entire application.
state = peer.state

# Server for listening to clients. Proxies etcd and transparently handles
# encryption/decryption of values sent/received by the client.
proxy = FastAPI()

@peer.post("/status")
async def status():
    """
    Return a status message indicating readiness to service requests.
    """
    return {"message": "Peer server is ready."}

@peer.post("/keygen/{id}")
async def keygen(id: str):
    """
    Return a share of the server's secret key and the server's public key.
    """
    async with httpx.AsyncClient(verify=ssl_context) as client:
        xⱼ = id
        k = state.secret_key
        Q = state.public_key
        kⱼ = k[xⱼ]

        return {"secretKeyShare": kⱼ, "publicKey": Q.to_dict()}
    '''
    xⱼ = id
    k  = state.secret_key
    Q  = state.public_key
    kⱼ = k[xⱼ]

    return {"secretKeyShare": kⱼ, "publicKey": Q.to_dict()}
    '''

@peer.post("/decrypt")
async def decrypt(request: Request):
    """
    Return the server's share of the decryption of a message.
    
    body = await request.json()

    kᵢ = state.joint_secret_key_share
    R  = Point.from_dict(body)
    Sᵢ = kᵢ * R

    return Sᵢ.to_dict()
    """
    try:
        body = await request.json()
        logger.info("Received decryption request.")

        kᵢ = state.joint_secret_key_share
        R = Point.from_dict(body)
        Sᵢ = kᵢ * R

        logger.info("Successfully processed decryption request.")
        return Sᵢ.to_dict()
    except Exception as e:
        logger.error(f"Error during decryption: {e}")
        raise


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
    C = base64.urlsafe_b64encode(C).decode("utf-8")

    # Replace the message with the public key and ciphertext.
    body["value"] = json.dumps({"publicKey": R.to_dict(), "ciphertext": C})
    '''
    # Carry on with the write request to etcd.
    async with httpx.AsyncClient(verify=ssl_context) as client:
        response = await client.put(f"{state.etcd}/v2/keys/{key}", data=body)

    return response.json()
    '''
    try:
        async with httpx.AsyncClient(verify=ssl_context) as client:
            response = await client.put(f"{state.etcd}/v2/keys/{key}", data=body)
            response.raise_for_status()
            logger.info(f"Successfully wrote key {key} to {state.etcd}")
            return response.json()
    except httpx.HTTPStatusError as e:
        logger.error(f"HTTP error occurred while writing key {key}: {e}")
        raise
    except httpx.RequestError as e:
        logger.error(f"Request error occurred while writing key {key}: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error while writing key {key}: {e}")
        raise
        
@proxy.get("/keys/{key}")
async def read(key: str):
    """
    Proxy GET requests to the upstream etcd API and perform decryption.
    """
    try:
        async with httpx.AsyncClient(verify=ssl_context) as client:
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
                    response = await client.post(f"{host}/decrypt", json=R.to_dict(), timeout=1)
                except:
                    # Peer cannot be reached. Put them at the end of the queue and
                    # come back later if we still require shares.
                    queue.append((xⱼ, host))
                    logger.info(f"{host} is unresponsive, moving on.")
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
    except httpx.HTTPStatusError as e:
        logger.error(f"HTTP error occurred while reading key {key}: {e}")
        raise
    except httpx.RequestError as e:
        logger.error(f"Request error occurred while reading key {key}: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error while reading key {key}: {e}")
        raise

class StateEncoder(json.JSONEncoder):
    """
    Custom JSON encoder for saving application state to disk.
    """
    def default(self, obj):
        if isinstance(obj, State):
            return {"__type__": "State", "value": obj._state}

        if isinstance(obj, Point):
            return {"__type__": "Point", "value": obj.to_dict()}

        if isinstance(obj, Curve):
            return {"__type__": "Curve", "value": obj.to_dict()}

        return super().default(obj)

class StateDecoder(json.JSONDecoder):
    """
    Custom JSON decoder for loading application state from disk.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, data):
        if data.get("__type__") == "State":
            return State(data["value"])

        if data.get("__type__") == "Point":
            return Point.from_dict(data["value"])

        if data.get("__type__") == "Curve":
            return Curve.from_dict(data["value"])

        return data

def save_state():
    """
    Save the application state to a file on the disk.
    """
    with open("/var/app/data/state.json", "w") as file:
        json.dump(state, file, cls=StateEncoder)

def load_state():
    """
    Load the application state from a file on the disk.
    """
    global state

    with open("/var/app/data/state.json", "r") as file:
        state = json.load(file, cls=StateDecoder)

def init_state():
    """
    Initialize a new application state.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--id", type=str)
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

    # Generate individual secret key shares and public key.
    k = random.randint(1, q)
    Q = k * G
    f = Polynomial.shamir(k, t, q)

    # TODO: Fix issue caused by JSON converting integer keys to strings.
    k = {str(xⱼ): f(xⱼ) for xⱼ in range(1, n + 1)}

    state.curve = E
    state.secret_key = k
    state.public_key = Q

    # Persist application state on the disk.
    save_state()

def init():
    """
    Initialize the application by loading its state from a file, if it exists,
    or else creating a new state.
    """
    try:
        load_state()
    except FileNotFoundError:
        init_state()

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
    Run the entire application.
    """
    asyncio.create_task(run_peer())

    # Elliptic-curve cryptography parameters.
    E = state.curve
    I = E.identity

    # Perform the distributed key generation algorithm.
    kᵢ = 0
    Q  = I

    async with httpx.AsyncClient() as client:
        for host in state.cluster:
            while True:
                try:
                    response = await client.post(f"{host}/keygen/{state.id}")
                except:
                    logger.info(f"{host} is unresponsive, trying again.")
                    await asyncio.sleep(1)
                else:
                    break
            body = response.json()

            kᵢ += body["secretKeyShare"]
            Q  += Point.from_dict(body["publicKey"])

    state.joint_secret_key_share = kᵢ
    state.joint_public_key = Q

    # Clients can now make requests.
    await run_proxy()

if __name__ == "__main__":
    init()

    try:
        logger.info("Initializing SSL context...")
        
        # Use state.id to dynamically select the peer-specific certificate and key
        certfile = f"/etc/certs/peer{state.id}.crt"
        keyfile = f"/etc/certs/peer{state.id}.key"
        
        ssl_context = ssl.create_default_context(cafile="/etc/certs/ca.crt")
        ssl_context.load_cert_chain(certfile=certfile, keyfile=keyfile)
        ssl_context.verify_mode = ssl.CERT_REQUIRED  # Require certificate verification
        logger.info(f"SSL context successfully configured for peer {state.id}.")
    except Exception as e:
        logger.error(f"Error configuring SSL context for peer {state.id}: {e}")
        raise
        
    asyncio.run(run())

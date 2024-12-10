import httpx
import random

from crypto import *
from fastapi import FastAPI, HTTPException, Request

etcd0 = "http://127.0.0.1:53667"
etcd1 = "http://127.0.0.1:53664"
etcd2 = "http://127.0.0.1:53661"

E = Curve('P-256')
q = E.order
G = E.generator
t = 3
n = 5
k = random.randint(1, q)
Q = k * G
f = Polynomial.Shamir(k, t, q)
k = {i: f(i) for i in range(1, n + 1)}

app = FastAPI()

@app.get("/status")
async def status():
    return {"message": "Server is ready"}

@app.put("/keys/{key}")
async def proxy_write(key: str, req: Request):
    """
    Proxy PUT requests to the upstream etcd API and perform encryption.
    """
    body = await req.json()

    # Generate an ephemeral key pair for this message.
    d = random.randint(1, q)
    R = d * G

    # Complete the ECIES protocol.
    M = body["value"]
    S = d * Q
    K = derive_key(S)
    C = encrypt(M, K)

    # Replace the message with the public key and ciphertext.
    body["value"] = EncryptedValue(R, C).to_json()

    # Proxy the write request to etcd.
    url = f"{etcd0}/v2/keys/{key}"

    async with httpx.AsyncClient() as client:
        res = await client.put(url, data=body)

    if res.status_code >= 400:
        raise HTTPException(status_code=res.status_code, detail=res.text)

    return res.json()

@app.get("/keys/{key}")
async def proxy_read(key: str):
    """
    Proxy GET requests to the upstream etcd API and perform decryption.
    """
    url = f"{etcd0}/v2/keys/{key}"

    async with httpx.AsyncClient() as client:
        res = await client.get(url)

    if res.status_code != 200:
        raise HTTPException(status_code=res.status_code, detail=res.text)

    body = res.json()

    # Extract the public key and ciphertext from the etcd value.
    encrypted = EncryptedValue.from_json(body["node"]["value"])
    R = encrypted.public_key
    C = encrypted.ciphertext

    # Interpolate the decryption shares and reveal the message.
    subset = random.sample(range(1, n + 1), t)
    shares = {i: k[i] * R for i in subset}

    S = interpolate_ecc(shares, q, R.point_at_infinity())
    K = derive_key(S)
    M = decrypt(C, K)

    # Replace the etcd value with the message.
    body["node"]["value"] = M

    return body

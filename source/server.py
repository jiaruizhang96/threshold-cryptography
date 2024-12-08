import httpx

from contextlib import asynccontextmanager
from fastapi import FastAPI
from pydantic import BaseModel

class Value(BaseModel):
    value: str

etcd0 = "http://127.0.0.1:57641"
etcd1 = "http://127.0.0.1:57647"
etcd2 = "http://127.0.0.1:57644"

@asynccontextmanager
async def lifespan(app: FastAPI):
    print("Starting up...")
    yield
    print("Shutting down...")

app = FastAPI(lifespan=lifespan)

@app.get("/status")
async def status():
    return {"message": "Server is ready"}

@app.put("/keys/{key}")
async def write(key: str, value: Value):
    return httpx.put(f"{etcd0}/v2/keys/{key}", data=value.model_dump()).json()

@app.get("/keys/{key}")
async def read(key: str):
    return httpx.get(f"{etcd0}/v2/keys/{key}").json()

@app.delete("/keys/{key}")
async def delete(key: str):
    return httpx.delete(f"{etcd0}/v2/keys/{key}").json()

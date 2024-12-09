from deepdiff import DeepDiff
from fastapi.testclient import TestClient
from server import app

client = TestClient(app)

ignore = [
    "root['node']['createdIndex']",
    "root['node']['modifiedIndex']",
    "root['prevNode']",
]

def test_status():
    response = client.get("/status")
    assert response.status_code == 200
    assert response.json() == {"message": "Server is ready"}

def test_write_then_read():
    response = client.put("/keys/message", json={"value": "Hello world"})
    assert response.status_code == 200

    response = client.get("/keys/message")
    assert response.status_code == 200

    expected = {
        "action": "get",
        "node": {
            "key": "/message",
            "value": "Hello world",
        },
    }
    assert DeepDiff(response.json(), expected, exclude_paths=ignore) == {}

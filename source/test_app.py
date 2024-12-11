import httpx

from deepdiff import DeepDiff

def test_status():
    response = httpx.get("http://127.0.0.1:8080/status")
    assert response.status_code == 200
    assert response.json() == {"message": "Proxy server is ready"}

def test_write_then_read():
    response = httpx.put("http://127.0.0.1:8080/keys/message", json={"value": "Hello world"})
    assert response.status_code == 200

    response = httpx.get("http://127.0.0.1:8080/keys/message")
    assert response.status_code == 200

    expected = {
        "action": "get",
        "node": {
            "key": "/message",
            "value": "Hello world",
        },
    }
    ignore = [
        "root['node']['createdIndex']",
        "root['node']['modifiedIndex']",
        "root['prevNode']",
    ]
    assert DeepDiff(response.json(), expected, exclude_paths=ignore) == {}

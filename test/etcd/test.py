import requests
import time

# Define the etcd server base URL
ETCD_SERVER_URL = "http://localhost:2379/v2/keys"

def put_key(key, value):
    """
    Puts a key-value pair into the etcd server and returns the time taken.
    """
    url = f"{ETCD_SERVER_URL}/{key}"
    data = {'value': value}
    try:
        start_time = time.time()
        response = requests.put(url, data=data)
        response.raise_for_status()  # Raise an error for bad HTTP responses
        end_time = time.time()
        return (end_time - start_time) * 1000  # Convert to milliseconds
    except requests.exceptions.RequestException as e:
        print(f"Error during PUT request: {e}")
        return None

def get_key(key):
    """
    Gets a key-value pair from the etcd server and returns the time taken.
    """
    url = f"{ETCD_SERVER_URL}/{key}"
    try:
        start_time = time.time()
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad HTTP responses
        end_time = time.time()
        return (end_time - start_time) * 1000  # Convert to milliseconds
    except requests.exceptions.RequestException as e:
        print(f"Error during GET request: {e}")
        return None

if __name__ == "__main__":
    # Common value for all keys
    test_value = "123456789"
    
    # Store unique keys for PUT and GET
    keys = [f"test_key_{i}" for i in range(100)]
    
    # Run 100 iterations for PUT
    put_times = []
    print("Performing 100 PUT operations...")
    for key in keys:
        time_taken = put_key(key, test_value)
        if time_taken is not None:
            put_times.append(time_taken)
    avg_put_time = sum(put_times) / len(put_times) if put_times else 0
    print(f"Average PUT time: {avg_put_time:.2f} ms")

    # Run 100 iterations for GET
    get_times = []
    print("\nPerforming 100 GET operations...")
    for key in keys:
        time_taken = get_key(key)
        if time_taken is not None:
            get_times.append(time_taken)
    avg_get_time = sum(get_times) / len(get_times) if get_times else 0
    print(f"Average GET time: {avg_get_time:.2f} ms")

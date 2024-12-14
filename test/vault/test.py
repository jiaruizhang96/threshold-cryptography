import time
import requests
import statistics

# Define the Vault server configuration
base_url = "http://127.0.0.1:8200/v1/secret/data/"
vault_token = "" # replace it with actual token
headers = {
    "X-Vault-Token": vault_token,
    "Content-Type": "application/json"
}
value = "123456789"

# Number of iterations for the benchmark
iterations = 100

# Generate unique keys for the test
keys = [f"benchmark_key_{i}" for i in range(iterations)]

# Measure PUT latency
put_latencies = []
for key in keys:
    url = base_url + key
    data = {"data": {"password": value}}
    start_time = time.time()
    response = requests.post(url, json=data, headers=headers)
    end_time = time.time()
    if response.status_code in [200, 204]:  # Vault API can return 204 for successful writes
        put_latencies.append((end_time - start_time) * 1000)  # Convert to ms
    else:
        print(f"PUT request failed for {key} with status code: {response.status_code} - {response.text}")

# Measure GET latency
get_latencies = []
for key in keys:
    url = base_url + key
    start_time = time.time()
    response = requests.get(url, headers=headers)
    end_time = time.time()
    if response.status_code == 200:
        try:
            response_data = response.json()
            if "data" in response_data and "password" in response_data["data"]["data"]:
                get_latencies.append((end_time - start_time) * 1000)  # Convert to ms
        except ValueError:
            print(f"GET request failed to parse JSON for {key} with response: {response.text}")
    else:
        print(f"GET request failed for {key} with status code: {response.status_code} - {response.text}")

# Calculate average latencies
average_put_latency = statistics.mean(put_latencies) if put_latencies else 0
average_get_latency = statistics.mean(get_latencies) if get_latencies else 0

# Output results
print(f"Average PUT latency: {average_put_latency:.2f} ms")
print(f"Average GET latency: {average_get_latency:.2f} ms")

# Optionally, print additional statistics
if put_latencies:
    print(f"PUT latency - Min: {min(put_latencies):.2f} ms, Max: {max(put_latencies):.2f} ms")
if get_latencies:
    print(f"GET latency - Min: {min(get_latencies):.2f} ms, Max: {max(get_latencies):.2f} ms")

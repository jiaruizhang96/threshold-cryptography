
import time
import requests
import statistics

# Define the endpoints and payload
base_url = "http://localhost:8080/keys/"
headers = {"Content-Type": "application/json"}
value = "123456789"

# Number of iterations
iterations = 100

# Generate unique keys
keys = [f"benchmark_key_{i}" for i in range(iterations)]

# Measure PUT latency
put_latencies = []
for key in keys:
    url = base_url + key
    data = {"value": value}
    start_time = time.time()
    response = requests.put(url, json=data, headers=headers)
    end_time = time.time()
    if response.status_code == 200:
        put_latencies.append((end_time - start_time) * 1000)  # Convert to ms
    else:
        print(f"PUT request failed for {key} with status code: {response.status_code}")

# Measure GET latency
get_latencies = []
for key in keys:
    url = base_url + key
    start_time = time.time()
    response = requests.get(url)
    end_time = time.time()
    if response.status_code == 200:
        get_latencies.append((end_time - start_time) * 1000)  # Convert to ms
    else:
        print(f"GET request failed for {key} with status code: {response.status_code}")

# Calculate average latencies
average_put_latency = statistics.mean(put_latencies) if put_latencies else 0
average_get_latency = statistics.mean(get_latencies) if get_latencies else 0

print(f"Average PUT latency: {average_put_latency:.2f} ms")
print(f"Average GET latency: {average_get_latency:.2f} ms")
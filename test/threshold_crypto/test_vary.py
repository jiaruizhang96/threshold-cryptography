import time
import json
import requests
import random
import matplotlib.pyplot as plt
import numpy as np
import statistics

# Define endpoints and payload
base_url = "http://localhost:8080/keys/"
headers = {"Content-Type": "application/json"}
value = "123456789"

# Workload definitions
workloads = {
    "A": {"read": 0, "write": 100},
    "B": {"read": 10, "write": 90},
    "C": {"read": 20, "write": 80},
    "D": {"read": 30, "write": 70},
    "E": {"read": 40, "write": 60},
    "F": {"read": 50, "write": 50},
    "G": {"read": 60, "write": 40},
    "H": {"read": 70, "write": 30},
    "I": {"read": 80, "write": 20},
    "J": {"read": 90, "write": 10},
    "K": {"read": 100, "write": 0},
}

# Generate keys
iterations = 100
keys = [f"benchmark_key_{i}" for i in range(iterations)]

# Function to perform PUT and GET requests
def perform_request(read_ratio, write_ratio, total_operations):
    latencies = []
    read_operations = int(total_operations * (read_ratio / 100))
    write_operations = int(total_operations * (write_ratio / 100))

    # Counters for operations
    remaining_reads = read_operations
    remaining_writes = write_operations
    cnt = 0 
    while remaining_reads > 0 or remaining_writes > 0:
        if remaining_writes > 0:
            operation = "write"
            remaining_writes -= 1
        elif remaining_reads > 0:
            operation = "read"
            remaining_reads -= 1

        key = keys[cnt]
        cnt += 1
        url = base_url + key

        start_time = time.time()
        if operation == "write":
            data = {"value": value}
            response = requests.put(url, json=data, headers=headers)
        else:  # read
            response = requests.get(url)
        end_time = time.time()

        if response.status_code == 200:
            latencies.append((end_time - start_time) * 1000)  # Convert to ms
        else:
            print(key)
            print(f"{operation} request failed: {response.status_code}")

    return statistics.mean(latencies) if latencies else 0.0


# Collect all latencies across workloads
results = {}
for workload, config in workloads.items():
    print(f"Running workload {workload}...")
    mean_latency = perform_request(read_ratio=config.get("read", 0), write_ratio=config.get("write", 0), total_operations=100)
    results[workload] = mean_latency

# Save latencies to disk
output_filename = "../../logs/vary/n=3/threshold_crypto.json"
with open(output_filename, "w") as f:
    json.dump(results, f, indent=4)
print(f"Latency data saved to {output_filename}")

'''
# Plot Combined CDF
plt.figure(figsize=(8, 6))
sorted_latencies = np.sort(all_latencies)
cdf = np.arange(1, len(sorted_latencies) + 1) / len(sorted_latencies)
plt.plot(sorted_latencies, cdf, label="Combined Workloads")

# Customize plot
plt.xscale('log')  # Use logarithmic scale for latency
plt.xlabel("Latency (ms)")
plt.ylabel("Proportion")
plt.title("Combined Latency CDF for All Workloads")
plt.legend()
plt.grid()
plt.show()
'''
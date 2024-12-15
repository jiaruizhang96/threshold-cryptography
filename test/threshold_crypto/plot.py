#import matplotlib.pyplot as plt
#import seaborn as sns
import pandas as pd
import numpy as np
import os 

# Data for the bar plots
systems = ['Our Solution', 'Vault', 'etcd']
configs = ['n=3', 'n=3', 'n=3', 'n=5', 'n=5', 'n=5']
put_times = [14.85, 27.35, 5.21, 20.52, 40.82, 6.72]  # Combined data for n=3 and n=5
get_times = [1.42, 0.297, 0.890, 2.16, 0.181, 0.729]  # Combined data for n=3 and n=5

# Create DataFrame for PUT and GET operations
put_data = pd.DataFrame({
    'System': systems * 2,  # Repeat systems for each configuration
    'Time (ms)': put_times,
    'Configuration': configs
})

get_data = pd.DataFrame({
    'System': systems * 2,  # Repeat systems for each configuration
    'Time (ms)': get_times,
    'Configuration': configs
})

# Function to create bar plots and add annotations
def create_plot(data, title, filename):

    # Define custom gradient palettes
    gradient_palette = {
        'Our Solution': '#bae4bc',
        'Vault': '#7bccc4',
        'etcd': '#2b8cbe',
    }

    plt.figure(figsize=(8, 6))
    barplot = sns.barplot(x='Configuration', y='Time (ms)', hue='System', data=data)
    plt.title(title)
    plt.ylabel('Time (ms)')
    plt.legend(title='System')

    # Add annotations for bar heights
    for p in barplot.patches:
        height = p.get_height()
        plt.annotate(f'{height:.2f}', (p.get_x() + p.get_width() / 2., height),
                     ha='center', va='center', xytext=(0, 10), textcoords='offset points')

    plt.savefig(f'{filename}.png')
    plt.show()


# Create and save plots
#create_plot(put_data, 'Benchmark Results for PUT Operations', 'put_benchmark_bar')
#create_plot(get_data, 'Benchmark Results for GET Operations', 'get_benchmark_bar')


def gen_avg_time(log_dir, n):
    """
    Compute the average execution time for secret split and restore from log files.
    
    Args:
        log_dir (str): The directory where the log files are stored.
        n (int): The number of servers to process.
        
    Returns:
        dict: A dictionary containing average secret split and restore times.
    """
    put_times = []
    get_times = []
    
    # Process PUT log files
    for i in range(1, n + 1):
        put_file = os.path.join(log_dir, f"put_{i}.txt")
        if os.path.exists(put_file):
            with open(put_file, "r") as file:
                times = [float(line.strip())*1000 for line in file if line.strip()]
                put_times.extend(times)
    
    # Process GET log files
    for i in range(1, n + 1):
        get_file = os.path.join(log_dir, f"get_{i}.txt")
        if os.path.exists(get_file):
            with open(get_file, "r") as file:
                times = [float(line.strip())*1000 for line in file if line.strip()]
                get_times.extend(times)
    
    # Compute averages
    avg_put_time = sum(put_times) / len(put_times) if put_times else 0
    avg_get_time = sum(get_times) / len(get_times) if get_times else 0
    
    # Return as a dictionary
    result = {
        "average_split_time": avg_put_time,
        "average_restore_time": avg_get_time
    }
    
    return result

# Example usage:
log_directory = "../../logs/n=3"  # Replace with the path to your log directory
num_servers = 3
averages = gen_avg_time(log_directory, num_servers)
print("Average Split Time (ms):", averages["average_split_time"])
print("Average Restore Time (ms):", averages["average_restore_time"])
    
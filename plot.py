import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

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
create_plot(put_data, 'Benchmark Results for PUT Operations', 'put_benchmark_bar')
create_plot(get_data, 'Benchmark Results for GET Operations', 'get_benchmark_bar')
'''

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Sample data
data = {
    'model1': {
        'system1': {'ner': [0.92, 0.91, 0.90, 0.89, 0.88, 0.90], 'cer': [0.82, 0.81, 0.80, 0.79, 0.78, 0.80]},
        'system2': {'ner': [0.92, 0.91, 0.90, 0.89, 0.88, 0.88], 'cer': [0.82, 0.81, 0.80, 0.79, 0.78, 0.78]},
        'system3': {'ner': [0.92, 0.91, 0.90, 0.89, 0.88, 0.89], 'cer': [0.82, 0.81, 0.80, 0.79, 0.78, 0.79]},
        'system4': {'ner': [0.92, 0.91, 0.90, 0.89, 0.88, 0.87], 'cer': [0.82, 0.81, 0.80, 0.79, 0.78, 0.77]}
    },
    'model2': {
        'system1': {'ner': [0.82, 0.81, 0.80, 0.79, 0.78, 0.80], 'cer': [0.72, 0.71, 0.70, 0.69, 0.68, 0.70]},
        'system2': {'ner': [0.82, 0.81, 0.80, 0.79, 0.78, 0.78], 'cer': [0.72, 0.71, 0.70, 0.69, 0.68, 0.68]},
        'system3': {'ner': [0.82, 0.81, 0.80, 0.79, 0.78, 0.79], 'cer': [0.72, 0.71, 0.70, 0.69, 0.68, 0.69]},
        'system4': {'ner': [0.82, 0.81, 0.80, 0.79, 0.78, 0.77], 'cer': [0.72, 0.71, 0.70, 0.69, 0.68, 0.67]}
    },
    'model3': {
        'system1': {'ner': [0.92, 0.91, 0.90, 0.89, 0.88, 0.89], 'cer': [0.82, 0.81, 0.80, 0.79, 0.78, 0.79]},
        'system2': {'ner': [0.92, 0.91, 0.90, 0.89, 0.88, 0.88], 'cer': [0.82, 0.81, 0.80, 0.79, 0.78, 0.78]},
        'system3': {'ner': [0.92, 0.91, 0.90, 0.89, 0.88, 0.90], 'cer': [0.82, 0.81, 0.80, 0.79, 0.78, 0.80]},
        'system4': {'ner': [0.92, 0.91, 0.90, 0.89, 0.88, 0.87], 'cer': [0.82, 0.81, 0.80, 0.79, 0.78, 0.77]}
    }
}

# Transforming data into DataFrame format
def prepare_data(data, metric):
    results = []
    for model, systems in data.items():
        for system, scores in systems.items():
            results.append({
                'Model': model,
                'System': system,
                'Value': scores[metric][-1]
            })
    return pd.DataFrame(results)

# Prepare NER and CER data
df_ner_avg = prepare_data(data, 'ner')
df_cer_avg = prepare_data(data, 'cer')

# Define custom gradient palettes
gradient_palette = {
    'system1': '#f0f9e8',
    'system2': '#bae4bc',
    'system3': '#7bccc4',
    'system4': '#2b8cbe'
}

# Plotting NER averages with gradient hues
plt.figure(figsize=(12, 6))
sns.barplot(x='Model', y='Value', hue='System', data=df_ner_avg, palette=gradient_palette)
plt.title('NER Averages Across Models and Systems')
plt.ylabel('NER Average')
plt.show()

# Plotting CER averages with gradient hues
plt.figure(figsize=(12, 6))
sns.barplot(x='Model', y='Value', hue='System', data=df_cer_avg, palette=gradient_palette)
plt.title('CER Averages Across Models and Systems')
plt.ylabel('CER Average')
plt.show()

'''
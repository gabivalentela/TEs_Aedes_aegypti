import os

variation_value = ['1000']
script = 'making_plot_mean_FST.py'

for variation in variation_value:
    command = f"python {script} {variation}"
    print(f"Running command: {command}")
    os.system(command)
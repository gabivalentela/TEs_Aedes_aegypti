import pandas as pd

df = pd.read_csv('read_lengths.csv')

# get max and min of read average_read_length column
max_read_length = df['average_read_length'].max()
min_read_length = df['average_read_length'].min()

print(f'Maximum read length: {max_read_length}')
print(f'Minimum read length: {min_read_length}')
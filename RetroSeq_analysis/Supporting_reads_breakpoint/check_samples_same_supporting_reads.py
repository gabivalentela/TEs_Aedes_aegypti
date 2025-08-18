import re

total_count = 0
count_eq = 0

with open("slurm-9929543.out", 'r') as f:
    # remove the last line from the file
    for line in f:
        if line.startswith('{'):
            total_count += 1
            count1 = re.search(r'{.+:(.+),.+}', line).group(1)
            count2 = re.search(r'{.+:.+,.+:(.+)}', line).group(1)
            if count1 == count2:
                count_eq += 1

print(f"Total count of lines starting with '{{': {total_count}")
print(f"Total count of equal lines: {count_eq}")

print(f'the ratio of equal lines to total lines is {count_eq / total_count:.2f}')
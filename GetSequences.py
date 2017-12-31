import csv
import os
import re
from collections import defaultdict

def get_sequences(filename):
    files = defaultdict(list)
    with open(filename, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            sample = row[1].split('-')[0]
            file = os.path.join(row[4], row[0])
            files[sample].append(file)
    for sample in files.keys():
        files[sample].sort() # need to make sure that R reads follow F reads
    return files


if __name__ == '__main__':
    files = get_sequences("Data/sequences.txt")
    for sample in files.keys():
        print(sample, ','.join(files[sample][::2]), ','.join(files[sample][1::2]))
        
#    for sample in set([string.split('-')[0] for string in files.keys()]):
#        print(sample, reduce(lambda x,y: x+y, [files[key] for key in filter(lambda key: key.startswith(sample), files.keys())]))

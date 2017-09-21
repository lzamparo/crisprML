import os
import sys
from __future__ import print_function

infile = sys.argv[1]
partitions = sys.argv[2]



# read the lines
with open(infile, 'r') as f:
    lines = f.readlines()

# calculate the load for each file
lines_per_partition = len(lines) - 1 // int(partitions)
header = lines[0]

# open a new partition file:
def new_part_file(infile, partition_number, header):
    infile_prefix, extension = os.path.splitext(infile)
    partition_name = infile_prefix + "_" + str(partition_number) + extension
    f = open(partition_name,'w')
    print(header, file=f)
    return f
    

# split the infile into each piece file, write out the guides
partition_num = 0
part_file = new_part_file(infile, partition_num, header)
print(header, file=part_file)

for count, guide in lines[1:]:
    print(guide, file=part_file)
    if count > lines_per_partition:
        part_file.close()
        partition_num = partition_num + 1
        part_file = new_part_file(infile, partition_number, header)

# close the remainind partition file
part_file.close()
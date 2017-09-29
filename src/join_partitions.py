from __future__ import print_function
import os
import sys


indir = sys.argv[1]

os.chdir(indir)
partition_dirs = [f for f in os.listdir('.')]
outfile = 'raw_features_computed_Cas9_sequences_fitting_annotation.csv'


def process_file(filename, outhandle, skip_header=True):
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('seq') and skip_header:
                continue
            print(line, file=outhandle)


def make_pf_filename(num):
    return 'raw_features_computed_Cas9_sequences_fitting_annotation' + '_' + str(num) + '.csv'


# read the lines
with open(outfile, 'w') as out:
    print("Processing partition")
    # process the 0th partition
    partfile = os.path.join('partition_0', make_pf_filename(0))
    process_file(partfile, out, skip_header=False)

    # process the rest
    for p in range(1:100):
        print('.', end='', flush=True)
        partdir = 'partition_' + str(p)
        partfile = os.path.join(partdir, make_pf_filename(p))
        process_file(partfile, out, skip_header=True)

print("aaand")
print("done!")
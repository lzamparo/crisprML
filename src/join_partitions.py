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


def get_pf_filename(pf_dir):
    csvs = [f for f in os.listdir(pf_dir) if f.endswith('.csv')]
    return csvs[0]
    #return 'raw_features_computed_Cas9_sequences_fitting_annotation' + '_' + str(num) + '.csv'


# read the lines
with open(outfile, 'w') as out:
    print("Processing partition")
    # process the 0th partition
    partfile = os.path.join('partition_0', get_pf_filename('partition_0'))
    process_file(partfile, out, skip_header=False)

    # process the rest
    part_dirs = [f for f in os.listdir('.') if not f.endswith('.csv')]
    part_dirs.sort()
    partindex = int(part_dirs[-1].split('_')[1]) + 1

    for p in range(1,partindex):
        print('.', end='', flush=True)
        partdir = 'partition_' + str(p)
        partfile = os.path.join(partdir, get_pf_filename(partdir))
        process_file(partfile, out, skip_header=True)

print("aaand")
print("done!")

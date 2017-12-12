#!/bin/bash

### identify experiments that did not produce data for all input guides

# Go to directory with outputs, count up number of scored guides
curr=`pwd`
cd $1

output_totals=$(find . -name "*.csv" -exec wc -l {} \; | cut -f1 -d' ' | awk '{ sum += $1 } END { print sum }')

# Go to directory with inputs, count up number of scored guides
prefix=$(dirname $1)
suffix=$(basename $1 | sed -e 's/input/output/')
cd $prefix/$suffix

input_totals=$(find . -name "*.txt" -exec wc -l {} \; | cut -f1 -d' ' | awk '{ sum += $1 } END { print sum }')

echo "Found $input_totals guides as input, and $output_totals guides as output"

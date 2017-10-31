#!/bin/bash

# find failed jobs

for dir in $(ls -1); do
    [ -z "`find $dir -type f`" ] && echo "$dir" >> empties.txt
done


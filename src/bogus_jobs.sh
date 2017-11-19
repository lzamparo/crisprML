#!/bin/bash

# find failed jobs
curr=`pwd`
cd $1

for dir in $(ls -1); do
    [ -z "`find $dir -type f`" ] && echo "$dir" >> empties.txt
done

mv empties.txt $curr

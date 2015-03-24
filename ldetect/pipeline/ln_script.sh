#!/bin/bash

# This script will just create links to original vectors with nice names
# Assuming the following directory structure:
# top-level-dir
# - population1
# - population2
# - ...
#
# It should be run from the top level dir and take the population name as an argument
# It also assumes the file name structure as shown below...

cd $1

for a in {1..22}
do
ln -s vector-orig_data_$1-chr$a-* vector-chr$a.txt.gz
done

ls -alh

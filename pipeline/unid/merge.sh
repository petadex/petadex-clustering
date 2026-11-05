#!/bin/bash
export LC_ALL=C
outer=$1
inner=$2

join -t$'\t' -1 1 -2 2 -o 2.1 1.2 <(sort.sh -k1,1 $outer) <(sort.sh -k2,2 $inner)

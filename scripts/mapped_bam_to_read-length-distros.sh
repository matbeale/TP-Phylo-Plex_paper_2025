#!/bin/bash

echo
echo "Take a bam file of mapped reads and outputs counts of the read length distributions"
if [ $# -ne 1 ]
then
 echo "Usage: $0 <mapped.sorted.bam>"
 echo
 exit 1
fi




samtools view $1 | awk '{print length($10)}' | sort | uniq -c | sort -n -k 2 > $1\.read-lengths.tsv

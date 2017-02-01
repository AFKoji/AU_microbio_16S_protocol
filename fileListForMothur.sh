#!/bin/bash

# making list of input files for Mothur 16S analysis
# run this script while you are in the folder containing your raw sequence files
# should be run with one input = short name of project (a.k.a. "basename")
# e.g. MetteMB@dnaseq:~/spider_microbiome/raw$ ./fileListForMothur.sh spider16s

# if you get the message: permission denied,
# first run: chmod u+x fileListForMothur.sh 

BASENAME="$1"

# sample names
# awk sets _ as separator and selects the first field
ls *R1_001.fastq.gz | awk -F'[_]' '{print $1}' > names.txt

# forward seq filenames
ls *R1_001.fastq.gz > forward.txt

# reverse seq filenames
ls *R2_001.fastq.gz > reverse.txt

# put them all together
paste names.txt forward.txt reverse.txt > $BASENAME.files

# remove the text files that are no longer needed rm forward.txt
rm forward.txt reverse.txt names.txt
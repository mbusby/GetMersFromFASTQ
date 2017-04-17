GetMersFromFASTQ
================

This is a script that makes a matrix of KMers present in a FASTQ file. 

Q: What does it do?

Parse the files to see if any Kmers are overrepresented in the data.


What would I want to do that?

To see what your reads are that are not mapping to see if you have contamination from barcodes, etc.  FASTQC does something like this but I couldn't find where they labeled the Y axis.  This also lets you do KMers of different sizes which may be helpful.


Q: How does it work? 

From the installation directory type:

./GetMersFastq -fastq first_read.fastq-out outputFile -mer 10

where 

outputFile is the name of the outputFile

first_read.fastq in the name of the fastq file you want to parse

10 is the number of nucleotides that are parsed


By default only 20% of the reads are processed to save time.  To change this alter the sampling_percent parameter, e.g.:

-sampling_percent 40


This utility should work for any fastq file though it may crash if you don't run it somewhere with enough memory.  


INSTALLATION

Compiling

Go to the folder where everything is installed.

Type:
g++ Main.cpp Handy.cpp -o GetMersFastq



License
Anyone can use this.

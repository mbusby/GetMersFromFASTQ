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

This is likely to be a pain in the neck.

Installation: The file GetMersFromFASTQ is a 64 bit Unix binary. This is easiest if it works on your server. Otherwise you will need to compile it.

How to get it to compile:

Compiling the file may or may not be a big hassle because you need bamtools.

First, though, copy all the files in here to a single directory on your unix server. The server should have g++ is installed (it probably is).

Bamtools

To get this to compile you will need to download Derek Barnett's bamtools API and install it in a folder somewhere. It is here: https://github.com/pezmaster31/bamtools You need the API, not the command line program though is quite useful to have. 

Edit the Makefile

Then, after you install it, YOU NEED TO EDIT THE Makefile in this folder so that everywhere it says "FolderWhereBamToolsIs" you put in the folder where bamtools is located on these lines:

LIBS= -L/FolderWhereBamToolsIs/bamtools/lib -lbamtools -lz

LDFLAGS = -Wl,-rpath /FolderWhereBamToolsIs/bamtools/lib

INCLUDES = -I/FolderWhereBamToolsIs/bamtools/include

Compiling

Go to the folder where everything is installed in type "make". Ignore the warnings about Handy doing stupid stuff. If nothing says ERROR and it makes an executable called ComplexityByStartPos you should be all set.

License
Anyone can use this.

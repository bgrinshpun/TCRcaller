# TCRcaller
Python tool to extract V,J, and CDR3 information from T cells for both A and B receptor chains.

Author: Boris Grinshpun, 2017

## To run demo

`sh TCRcaller.sh -i ../demo.bam -o demo`

Results are stored in a newly created demo directory

FOR HELP WITH OPTIONS
`sh TCRcaller.sh -h`

## Dependencies:
SAMtools: http://samtools.sourceforge.net/

Burrows Wheel Aligner (BWA): http://bio-bwa.sourceforge.net/

If the input is a fastq file and must first be mapped, a properly indexed fasta reference file is required.

# TCRcaller
Python tool to extract V,J, and CDR3 information from T cells for both A and B receptor chains.

Author: Boris Grinshpun, 2017

## To run demo

`sh TCRcaller.sh -i ../demo.bam -o demo`

Results are stored in a newly created demo directory

FOR HELP WITH OPTIONS
`sh TCRcaller.sh -h`

## File Outputs
* TCRX.final.tsv  &ndash;  Final output file with number of occurences of a CDR3 defined by nucleotide sequence, amino acid sequence, Vgene and Jgene.
* VJ.called.TCRX.tsv  &ndash;  List of nucleotide sequence, amino acid sequence, Vgene and Jgene by individual read.
* VorJonly.TCRX.tsv  &ndash;  List of reads where only one cassette (V or J) mapped to the reference.
* discarded.TCRX.tsv  &ndash;  Reads that were discarded for a variety of reasons -- no mapping, more than one V or J cassette mapping, out of frame, etc.

## Dependencies
SAMtools: http://samtools.sourceforge.net/

Burrows Wheel Aligner (BWA): http://bio-bwa.sourceforge.net/

If the input is a fastq file and must first be mapped, a properly indexed fasta reference file is required.

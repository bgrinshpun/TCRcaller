sample=$1

outdir=~/HC42_Nina/data/processed/
perl ReadSam.CDR3DNAseqTRB_8.pl $outdir $sample
perl ReadCDR3DNA_OutputCDR3prot.TRBbyscore.VJ_3.pl $outdir `basename $sample `.CDR3.VJ.seq
perl CategorizeCDR3_TRB_3.pl $outdir `basename $sample`.CDR3.VJ.seq.prot.txt 

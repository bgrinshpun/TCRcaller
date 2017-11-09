fqfile=$1
outdir=$2
ref=$3
name=`basename $fqfile`

bwa mem -t 8 -T 10 -r 1 -O 4 -Y $ref $fqfile | samtools sort -n -@8 -O BAM -o $outdir/$name.bam -


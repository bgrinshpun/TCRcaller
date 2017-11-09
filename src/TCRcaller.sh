#!/bin/sh

# Run TCRcaller pipeline to extract CDR3 and VJ information
# show help if no arguments are passed
#
#TODO:
#      -- better documentation for discard pile
#      -- cython implementation
#      -- suppress line timing / give only total time
#      -- paired end fastqs
#      -- k-mer error correction (fastqs only)
#      -- V,J nucleotide sequence output with variants
#      -- D cassette
#      -- keep track of sequences that have zero chromosome match

DIR=`pwd`
SCRIPT_PATH=$(dirname $0)
scriptname=$(basename $0)

function usage() {
    echo >&2
    echo >&2 "$scriptname: Identify T cell CDR3 sequences and VJ cassettes from raw fastq or mapped bwa files"
    echo >&2 
    echo USAGE: >&2
    echo >&2 "  " $scriptname [-i path/to/input_file] [-o output_name] [-c chain] [-r path/to/genome_reference]
    echo >&2
    echo >&2 "INPUTS"
    echo >&2 "\t-i: Specify an input file either in fastq, sam, or bam format."
    echo >&2 "\t    Presently only single-end reads are allowed." 
    echo >&2
    echo >&2 "\t-o: Output name. A directory with this name will be created."
    echo >&2 "\t    Default [output]"
    echo >&2
    echo >&2 "\t-c: TCR chain to process. Must be either A (TCRA) or B (TCRB)."
    echo >&2 "\t    If this option is not included both chains are processed."
    echo >&2
    echo >&2 "\t-r: Path to reference genome which has been properly for use with bwa software"
    echo >&2 "\t    Can be excluded if input file is already aligned"
    echo >&2
    echo >&2
    echo >&2 "RUN DEMO"
    echo >&2 "\t sh TCRcaller.sh -i ../demo.bam -o demo"
    echo >&2
    echo >&2
    echo >&2 "Author: Boris Grinshpun (bg2178@columbia.edu), 2017"
    echo >&2
}

if [ $# -lt 1 ] 
    then
    usage
    exit 1
fi

showAll=0
# process options
while getopts ahi:o:c:r: opt
  do
  case "$opt" in
      i) file=$OPTARG ;;
      o) output=$OPTARG ;;
      c) chain=$OPTARG ;;
      r) referencePath=$OPTARG ;;
      h) usage; exit 1 ;;
      \?) usage; exit 1 ;;
        
  esac
done

shift $((OPTIND -1))

# Paths to dependencies ----------------------------
BWA=`which bwa`
SAMTOOLS=`which samtools`

if  [ -z "$BWA" ]; then
	echo >&2 ;
	echo >&2 "BWA software installation not found";
	echo >&2 "If no mapping is to be done, please remove line 75 of this script"
	echo >&2 ;
	echo >&2 "EXITING PROGRAM (except if manually turned off)"
	exit 1; # possible removal
fi

if  [ -z "$SAMTOOLS" ]; then
        echo >&2 ;
        echo >&2 "SAMTOOLS software installation not found";
        echo >&2 ;
        echo >&2 "EXITING PROGRAM"
        exit 1;
fi

# included reference files ----------------------------
if [ $SCRIPT_PATH == "." ]; then
        refdir="../ref"
else
	refdir="$(dirname $SCRIPT_PATH)/ref";

fi

codonmap=$refdir/codonmap.txt
TCRAcoords=$refdir/TRA.v37.Coordinates.txt
TCRBcoords=$refdir/TRB.v37p8.Coordinates.txt
TCRA_Vmotifs=$refdir/TRAV.motifs.IMGT.txt
TCRA_Jmotifs=$refdir/TRAJ.motifs.IMGT.txt
TCRB_Vmotifs=$refdir/TRBV.motifs.IMGT.txt
TCRB_Jmotifs=$refdir/TRBJ.motifs.IMGT.txt

# determine input file type ----------------------------

filetype=''
if [ ! -f $file ]; then
       echo >&2 ;
       echo >&2 "Input file not found"
       echo >&2 ;
       echo >&2 "EXITING PROGRAM"
       exit 1;
elif [[ $file == *.sam ]] || [[ $file == *.bam ]]; then
	filetype='mapped'
else
	filetype="raw" # TODO: check explicitly for fastq
fi


# check for reference files ----------------------------
if [ "$filetype" == 'raw' ]; then
	if [ ! -f $referencePath ]; then
		echo >&2 ;
		echo >&2 "Path to reference genome does not exist." 
		echo >&2 "Run TCRcaller.sh -h to see Usage";
		echo >&2 ;
		echo >&2 "EXITING PROGRAM"
		exit 1;
	elif [ ! -f $referencePath.bwt ]; then
		echo >&2 ;
		echo >&2 "Reference exists but bwa compatible file not found";
		echo >&2 "Please check that a \"`basename $referencePath`.bwt\" file exists in the reference directory";
		echo >&2 ;
		echo >&2 "EXITING PROGRAM";
		exit 1;
	fi
fi


# prepare output folder ----------------------------
[ -z $output ] && output="output";
mkdir -p $output
echo
echo "INPUT FILE `basename $file`"
echo
echo "RESULTS WILL BE PLACED IN \"$DIR/$output/\""
echo "------------"
echo

# map fastq  ----------------------------
if [ "$filetype" == 'raw' ]; then
     echo 
     echo \"RUNNING BWA MAPPING\"
     sh bwamap.sh $file $output $referencePath
     echo "------------"
     echo
     mappedfile=$output/
else
     mappedfile=$file
fi

# where the magic happens [TCR identification]  ----------------------------
if [ "$chain" == "A" ]; then
	echo "Calling alpha chain only";
	echo;
	refchr="14";
	samtools view $mappedfile | python findTCR.py $output "A" $TCRAcoords $TCRA_Vmotifs $TCRA_Jmotifs $codonmap $refchr
	echo "alpha chain processed";
	echo
	echo "compiling sequence counts"
	echo
	echo "count\tnucleotide\tcdr3\tVgene\tJgene" > $output/TCRA.final.tsv
        sed 1d $output/VJcalled.TCRA.tsv | cut -f2,3,4,5 | sort | uniq -c | awk 'OFS="\t"{print $1,$2,$3,$4,$5}' | sort -rnk1  >> $output/TCRA.final.tsv
	echo
	echo "FINISHED"

elif [ "$chain" == "B" ]; then
	echo "Calling beta chain only";
	refchr="gi|114841177:91557-667340";
        samtools view $mappedfile | python findTCR.py $output "B" $TCRBcoords $TCRB_Vmotifs $TCRB_Jmotifs $codonmap $refchr
	echo;
        echo "beta chain processed";
        echo
        echo "compiling sequence counts"
        echo
	echo $output
        echo "count\tnucleotide\tcdr3\tVgene\tJgene" > $output/TCRB.final.tsv
        sed 1d $output/VJcalled.TCRB.tsv | cut -f2,3,4,5 | sort | uniq -c | awk 'OFS="\t"{print $1,$2,$3,$4,$5}' | sort -rnk1  >> $output/TCRB.final.tsv
        echo
	echo "FINISHED"

else 
	echo "Calling both chains";
	echo
	echo "* alpha chain in progress *";
        refchr="14";
        samtools view $mappedfile | python $SCRIPT_PATH/findTCR.py $output "A" $TCRAcoords $TCRA_Vmotifs $TCRA_Jmotifs $codonmap $refchr
	echo
        echo "alpha chain processed";
        echo
        echo "compiling sequence counts"
        echo
        echo "count\tnucleotide\tcdr3\tVgene\tJgene" > $output/TCRA.final.tsv
        sed 1d $output/VJcalled.TCRA.tsv | cut -f2,3,4,5 | sort | uniq -c | awk 'OFS="\t"{print $1,$2,$3,$4,$5}' | sort -rnk1  >> $output/TCRA.final.tsv
	echo
	echo
	echo "* beta chain in progress *";
        refchr="gi|114841177:91557-667340";
        samtools view $mappedfile | python $SCRIPT_PATH/findTCR.py $output "B" $TCRBcoords $TCRB_Vmotifs $TCRB_Jmotifs $codonmap $refchr
        echo "beta chain processed";
	echo
        echo "compiling sequence counts"
        echo
        echo "count\tnucleotide\tcdr3\tVgene\tJgene" > $output/TCRB.final.tsv
        sed 1d $output/VJcalled.TCRB.tsv | cut -f2,3,4,5 | sort | uniq -c | awk 'OFS="\t"{print $1,$2,$3,$4,$5}' | sort -rnk1  >> $output/TCRB.final.tsv
        echo
        echo "FINISHED"

fi


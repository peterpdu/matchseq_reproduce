#!/bin/bash

#SBATCH -c 32
#SBATCH -t 36:00:00
#SBATCH --partition=normal,satpathy
#SBATCH --mem 222G
#SBATCH -o proc_4barcodes."%j".out
#SBATCH -e proc_4barcodes."%j".err

source ~/miniconda3/bin/activate /home/groups/bassik/ppdu/miniconda3/envs/py3
# requires regex package

ml biology samtools bowtie

# usage: proc_4_10X_barcodes.sh R1.fastq R2.fastq cell_whitelist.txt sgRNA.chip.txt
# fastq files must be unzipped
# cell_whitelist does not have header (can be unzipped output from cellranger, barcodes.tsv)
# sgRNA.chip.txt has header

# args
CB_MM_ALLOWED=1
SG_MM_ALLOWED=2
NPROC=32
# must be multiple of 4
CHUNK_SIZE=4000000

# used in step 1a
R1_UMI_PATTERN='^.{16}(?P<umi_1>.{12})(?P<discard_1>.*)'

# used in step 2b
R1_PATTERN='CCCCCCCCCCCCCCCCNNNNNNNNNNNN'

# this will take the first G, and remove the last nt from guides that do not start with G in order to always have a 20nt spacer
R2_PATTERN='.+(?:ACACC){e<=1}(?P<sgRNA>.{20})(?:GTTTA){e<=1}.+(?:GTGCTTTTTTT){e<=1}(?P<UMI>.{16})'
R2_UMI_PATTERN='^.{20}(?P<umi_1>.{16})'  # (20 sgRNA, then 16 sgUMI)


# ------------------------

# R1, R2 must be unzipped fastq
R1="$1"
R2="$2"
CELL_WHITELIST="$3"
SGRNA_CHIP="$4"

CHUNK_DIR=`dirname $R1`/tmp
mkdir -p $CHUNK_DIR

echo '0a. building cell barcode bowtie index'
CELL_WHITELIST_FASTA=${CELL_WHITELIST%.txt}.fasta
CELL_IDX=`dirname $CELL_WHITELIST`/cell_bowtie-idx/cell

awk -F"\t" '{ {print ">"$1"\n"$1}}' $CELL_WHITELIST > $CELL_WHITELIST_FASTA
mkdir -p `dirname $CELL_IDX`
bowtie-build $CELL_WHITELIST_FASTA $CELL_IDX > `dirname $CELL_IDX`/build.log

echo '0b. building sgRNA barcode bowtie index'
SGRNA_FASTA=${SGRNA_CHIP%.txt}.fasta
SGRNA_IDX=`dirname $SGRNA_CHIP`/sgRNA_bowtie-idx/sgRNA

if [ ! -f $SGRNA_FASTA ]; then
        # add leading G
        awk -F"\t" '{ if (NR>1) {if ($1 ~ /^G/) {print ">"$2"_"$1"\n"$1} else {print ">"$2"_"$1"\nG"$1}}}' $SGRNA_CHIP > $SGRNA_FASTA
fi


if [ ! -d `dirname $SGRNA_IDX` ]; then
        mkdir -p `dirname $SGRNA_IDX`
        bowtie-build $SGRNA_FASTA $SGRNA_IDX > `dirname $SGRNA_IDX`/build.log
fi

echo
echo '1a. extracting R1 UMI to readname'
R1_UMI=${R1%.fastq}.UMI.fastq

R1_CHUNK_DIR=$CHUNK_DIR/R1_UMI
mkdir -p $R1_CHUNK_DIR
split -l $CHUNK_SIZE $R1 $R1_CHUNK_DIR/chunk_
parallel -I ,, umi_tools extract --extract-method=regex --bc-pattern \'$R1_UMI_PATTERN\' --stdin ,, --stdout ,,.extract.fastq --log ,,.log ::: $R1_CHUNK_DIR/chunk_*
cat $R1_CHUNK_DIR/*.extract.fastq > $R1_UMI

# umi_tools extract --extract-method=regex --bc-pattern=$R1_UMI_PATTERN --stdin "$R1" --stdout $R1_UMI

echo
echo '1b. aligning 10X cell barcode to whitelist'
R1_MAP="${R1%.fastq}".bowtie.map.txt
R1_MAP_FASTQ="${R1_MAP%.txt}".fastq
R1_MAP_RN="${R1_MAP%.txt}".readname.txt

R1_UNMAP="${R1%.fastq}".bowtie.unmap.fastq
CB_LOG="${R1%.fastq}".align.log

CB_ALIGN_DIR=`dirname $R1`/cell_bowtie-align
mkdir -p $CB_ALIGN_DIR

# -n is the number of tolerated mismatches. (seed length is 28 by default, which covers the entire read)
# -p is the number of processors to use.
# -q is the input file and indicates that it is fastq.
# -k 1 specifies to only report 1 alignment per read
# --best ensures the reported alignment is the best one
# mapFile stores the output alignments.
# --un stores the reads with no alignments

echo "bowtie -n $CB_MM_ALLOWED -k 1 --best -p $NPROC $CELL_IDX -q $R1_UMI $R1_MAP --un $R1_UNMAP" >> $CB_LOG
bowtie -n $CB_MM_ALLOWED -k 1 --best -p $NPROC "$CELL_IDX" -q "$R1_UMI" "$R1_MAP" --un "$R1_UNMAP" >> "$CB_LOG" 2>&1
cat "$CB_LOG"

# write corrected to fastq
# awk -F '\t' '{print "@"$1"\n"$3"\n+\n"$6}' $MAP > $MAP_FASTQ

# sets phred score of all UMI bases to "I", doesn't matter - not used
# removes "-1" if it exists in cell barcodes file
awk -F '\t' 'BEGIN{s=sprintf("%12s","");gsub(/ /,"I",s)} {split($1,a," "); split(a[1],b,"_"); sub(/-1/, "", $3); print "@"b[1],a[2]"\n"$3b[2]"\n+\n"$6s}' $R1_MAP > $R1_MAP_FASTQ
awk '(NR%4==1){print substr($1,2)}' $R1_MAP_FASTQ > $R1_MAP_RN

echo
echo '2a. filtering R2 based on R1'
R2_MATCHED=${R2%.fastq}.matched.fastq
seqtk subseq $R2 $R1_MAP_RN > $R2_MATCHED

R2_MATCHED_SORTED=${R2_MATCHED%.fastq}.sort.fastq
R1_MATCHED_SORTED=${R1_MAP_FASTQ%.fastq}.sort.fastq

samtools import -0 $R2_MATCHED | samtools sort -n -@ $NPROC -u | samtools fastq > $R2_MATCHED_SORTED
samtools import -0 $R1_MAP_FASTQ | samtools sort -n -@ $NPROC -u | samtools fastq > $R1_MATCHED_SORTED

echo
echo '2b. transferring R1 barcodes'
R2_UMI=${R2_MATCHED%.fastq}.UMI.fastq

umi_tools extract \
        --bc-pattern=$R1_PATTERN \
        --stdin "$R1_MATCHED_SORTED" --read2-in "$R2_MATCHED_SORTED" \
        --read2-stdout --stdout $R2_UMI

echo
echo '2c. extracting sgRNA and sgUMI from R2'

R2_UMI_EXTRACT=${R2_UMI%.fastq}.extract.fastq
R2_CHUNK_DIR=$CHUNK_DIR/R2_regex
mkdir -p $R2_CHUNK_DIR

# split to parallelize, make sure line number is divisible by 4
split -l $CHUNK_SIZE $R2_UMI $R2_CHUNK_DIR/chunk_

parallel -I ,, python ~/scripts/regex_fastq_v2.py -fq ,, --regex \'$R2_PATTERN\' -q ::: $R2_CHUNK_DIR/chunk_*
# python ~/scripts/regex_fastq_v2.py -fq $R2_UMI --regex $R2_PATTERN --out $R2_UMI_EXTRACT

cat $R2_CHUNK_DIR/*.fastq > $R2_UMI_EXTRACT

echo
echo '2d. extracting R2 sgUMI to readname'
R2_sgUMI=${R2_UMI_EXTRACT%.fastq}.sgUMI.fastq

R2sg_CHUNK_DIR=$CHUNK_DIR/R2_sgUMI
mkdir -p $R2sg_CHUNK_DIR
split -l $CHUNK_SIZE $R2_UMI_EXTRACT $R2sg_CHUNK_DIR/chunk_
parallel -I ,, umi_tools extract --extract-method=regex --bc-pattern \'$R2_UMI_PATTERN\' --stdin ,, --stdout ,,.extract.fastq --log ,,.log ::: $R2sg_CHUNK_DIR/chunk_*
cat $R2sg_CHUNK_DIR/*.extract.fastq > $R2_sgUMI

# umi_tools extract --extract-method=regex --bc-pattern=$R2_UMI_PATTERN --stdin "$R2_UMI_EXTRACT" --stdout $R2_sgUMI

echo
echo '3a. aligning sgRNA to library'
R2_MAP="${R2_sgUMI%.fastq}".bowtie.map.sam
R2_UNMAP="${R2_sgUMI%.fastq}".bowtie.unmap.fastq
SG_LOG="${R2_sgUMI%.fastq}".align.log

SG_ALIGN_DIR=`dirname $R2_sgUMI`/sgRNA_bowtie-align
mkdir -p $SG_ALIGN_DIR


echo "bowtie -n $SG_MM_ALLOWED -k 1 --best -p $NPROC $SGRNA_IDX -q $R2_sgUMI -S $R2_MAP --un $R2_UNMAP" >> $SG_LOG
bowtie -n $SG_MM_ALLOWED -k 1 --best -p $NPROC "$SGRNA_IDX" -q "$R2_sgUMI" -S "$R2_MAP" --un "$R2_UNMAP" >> "$SG_LOG" 2>&1
cat "$SG_LOG"


echo
echo '3b. moving barcodes to tags'
echo '10X CB to XC tag'
echo '10X UMI to XU tag'
echo 'sgUMI to XG tag'
R2_MAP_TAG="${R2_MAP%.sam}".tag.sam
# this is the input to error correction python script
R2_MAP_TABLE="${R2_MAP%.sam}".barcode_table.txt

samtools view -H $R2_MAP > $R2_MAP_TAG
samtools view --no-header $R2_MAP | awk -F '\t' 'BEGIN{OFS="\t"}{split($1,a,"_"); $1=a[1]; print $0"\tXC:Z:"a[2]"\tXU:Z:"a[3]"\tXG:Z:"a[4]}' >> $R2_MAP_TAG
samtools view -F 4 --no-header $R2_MAP | awk -F'\t' 'BEGIN{OFS="\t"; print "cell\tcell_UMI\tsgRNA\tsgUMI"} {split($1,a,"_"); print a[2]"\t"a[3]"\t"$3"\t"a[4]}' > $R2_MAP_TABLE



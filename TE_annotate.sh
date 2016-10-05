#!/bin/bash/

# Transpoable Element Annotation Pipeline

# -----------------------------------------------------------------------------
# AUTHORS
# -----------------------------------------------------------------------------
# Michelle Hwang
# Regina Baucom
# Baucom Lab
# University of Michigan, Ann Arbor 

# -----------------------------------------------------------------------------
# LATEST UPDATE
# -----------------------------------------------------------------------------
# July 2016

# -----------------------------------------------------------------------------
# REQUIRED SOFTWARE
# -----------------------------------------------------------------------------
# Python
# MCL Cluster Algorithm
# BLAST+

# -----------------------------------------------------------------------------
# INCLUDED FILES
# -----------------------------------------------------------------------------
# cnv_mcl2na.pl 
# annotate.py
# sum_bp.py

# -----------------------------------------------------------------------------
# HOW TO USE 
# -----------------------------------------------------------------------------
# sh TE_annotate.sh reads_file db_file output_name eval length bit pident annotate_thresh annotate2_thresh fam_thresh

# EXAMPLE:
# sh TE_annotate.sh TEST_reads.fasta TE_DB.fasta test 1 0 0 0 .2 0.5 5

# All output files will be stored in the directory that this script is run.

# -----------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# -----------------------------------------------------------------------------

RAW_READS=$1 # File of your sequenced reads in FASTA format
TE_DB=$2 # File of your transposable element database in FASTA format
NAME=$3 # Will be used to name all files
EVAL=$4 
LEN=$5
BIT=$6
PID=$7
ANN=$8
ANN2=$9
FAM=${10}


echo -e ">All by all BLAST"
makeblastdb -in $RAW_READS -dbtype nucl -out $NAME.db
blastn -task blastn -num_threads 8 -db $NAME.db -query $RAW_READS -out BLAST_ALL_$NAME.txt -evalue 1e-10 -outfmt "6 qseqid sacc bitscore"

echo -e ">BLAST to NCBI"
blastn -task blastn -num_threads 8 -query $RAW_READS -out BLAST_NCBI_$NAME.txt -db nt -evalue 1e-1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -max_target_seqs 1 -max_hsps 1

echo -e ">BLAST to TE DB"
makeblastdb -in $TE_DB -dbtype nucl -out $TE_DB.db
blastn -task blastn -num_threads 8 -db $TE_DB.db -query $RAW_READS -out BLAST_TE_$NAME.txt -evalue 1e-1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -max_target_seqs 1 -max_hsps 1

echo -e ">Run MCL Clust"
mcl BLAST_ALL_$NAME.txt  --abc

echo -e ">Transpose MCL Output"
perl cnv_mcl2na.pl -i out.BLAST_ALL_$NAME.txt.I20 -o MCL_$NAME.txt
sed -i -e "1d" MCL_$NAME.txt

echo -e ">Run annotation"
python annotate-v2.py MCL_$NAME.txt BLAST_TE_$NAME.txt $RAW_READS -n BLAST_NCBI_$NAME.txt -pre $NAME -e $EVAL -l $LEN -b $BIT -p $PID -t $ANN -st $ANN2 -ft $FAM

#### OPTIONAL
echo -e ">Generate table"
# Summary of all seqs, one per line
grep '^>' ANNOTATE_$NAME.fasta | tr -d ">" | tr '#' \\t > TABLE_$NAME.txt

echo -e ">Getting family bp"
# Summary of all families, one per line, and the total BP in each
join FAM_INFO_$NAME.txt <(python sum_bp.py ANNOTATE_$NAME.fasta) > FAM_SUMS_$NAME.txt



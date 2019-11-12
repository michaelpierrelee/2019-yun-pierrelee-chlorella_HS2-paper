#!/bin/bash

cd /media/hdd/pierrelee/korea_hs2;

DIR="degs_annotation/"
FASTA="longest_orfs.pep"
CPU=16
TEMP_DIR=$DIR"temp_fasta_subsets/"
OUT="longest_orfs"

# number of lines per CPU

N=`more ${DIR}${FASTA} | wc -l` # number of lines
X=$(( $N/$CPU - ($N/$CPU)%2 )) # number of lines in the first $CPU-1 CPUs
Y=$(( $N - $X * ( $CPU - 1 ) )) # number of lines in the last CPU

# split fasta file into temporary smaller files

LIST_SUBSETS=()
TEST=0
for i in $( seq 1 $CPU ); do
	# for the first CPUs
	if (( $i < $CPU )); then
		SEQ=`head -n $(( $i * $X )) ${DIR}${FASTA} | tail -n $(( $X ))`
	# for the last CPU
	else
		SEQ=`more ${DIR}${FASTA} | tail -n $(( $Y ))`
	fi
	# create sub fasta files
	NAME="${FASTA}_subset_${i}"
	echo "$SEQ" > ${TEMP_DIR}${NAME}.fasta.temp
	TEST=$(( $TEST + `more ${TEMP_DIR}${NAME}.fasta.temp | wc -l` ))
	# run hmmscan (output is throw away)
	hmmscan -o /dev/null --cpu 1 --domtblout ${TEMP_DIR}${NAME}.pfam.domtblout.temp \
		../databases/Pfam-A.hmm-30AUG2018/data \
		${TEMP_DIR}${NAME}.fasta.temp 2> ${TEMP_DIR}${NAME}.log &
done

echo $(( ($CPU-1)*$X+$Y )) = $N = $TEST

echo "waiting to finish..."

wait

# concatenate result files

RESULTS=""
i=0
for FILE in `ls ${TEMP_DIR}*.pfam.domtblout.temp`; do
	((i++))
	if (( $i == 1 )); then RESULTS=`more $FILE | head -n -10`
	else RESULTS=${RESULTS}$'\n'`tail -n +4 $FILE | head -n -10`
	fi
done

echo "$RESULTS" > ${DIR}${OUT}.pfam.domtblout

echo "finished."

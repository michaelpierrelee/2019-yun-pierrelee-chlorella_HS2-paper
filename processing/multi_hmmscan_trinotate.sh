#!/bin/bash

cd /media/hdd/pierrelee/korea_hs2;

DIR="degs_annotation/"
FASTA="degs_to_annotate.fasta.transdecoder.pep"
FASTA_FILE=$DIR$FASTA
CPU=16
TEMP_DIR=$DIR"temp_transdecoder_subsets/"
OUT="trinotate"

# number of lines per CPU

N=`more $FASTA_FILE | wc -l` # number of lines
X=$(( $N / $CPU )) # theorical repartition of lines across CPUs

# split fasta file (with sequences splitted over several lines) into temporary smaller files

POINTER_TOP=1
POINTER_BOT=$X

TEST=0
RES=""
for i in $( seq 1 $CPU ); do
	# for the first CPUs
	if (( $i < $CPU )); then
		# we search the closest head toward the bottom
		LAST_LINE="`sed "${POINTER_BOT}q;d" $FASTA_FILE`"
		while [[ $LAST_LINE != \>* ]] && (( $POINTER_BOT + 1 <= $N )); do
			((POINTER_BOT++))
			LAST_LINE="`sed "${POINTER_BOT}q;d" $FASTA_FILE`"
		done
		# we find the head, so we split just before it
		if [[ $LAST_LINE == \>* ]]; then
			SEQ=`head -n $(( $POINTER_BOT - 1 )) $FASTA_FILE | tail -n $(( $POINTER_BOT - $POINTER_TOP  ))`
		# or we already reached the end
		else
			echo "X"
			SEQ=`tail -n $(( $N - $POINTER_TOP + 1 )) $FASTA_FILE`
		fi
		# define the pointers
		POINTER_TOP=$POINTER_BOT
		POINTER_BOT=$(( $POINTER_BOT + $X ))
	# for the last CPU, we take everything from the end
	else
		SEQ=`tail -n $(( $N - $POINTER_TOP + 1  )) $FASTA_FILE`
	fi
	# create sub fasta files
	NAME="${FASTA}_subset_${i}"
	echo "$SEQ" > ${TEMP_DIR}${NAME}.fasta.temp
	TEST=$(( $TEST + `more ${TEMP_DIR}${NAME}.fasta.temp | wc -l` ))
	# hmmscan
	hmmscan -o /dev/null --cpu 1 --domtblout ${TEMP_DIR}${NAME}.pfam.domtblout.temp \
		../databases/Pfam-A.hmm-30AUG2018/data \
		${TEMP_DIR}${NAME}.fasta.temp 2> ${TEMP_DIR}${NAME}.log &
	# progress
	echo "job sent to CPU" $i
	# check
	if (( $i==1 )); then RES=$SEQ
	else RES=$RES$'\n'$SEQ; fi
	# stop if we sent all sequences to CPUs while there are still unused CPUs
	if (( $POINTER_TOP >= $N )); then break; fi
done

echo $N = $TEST = `echo "$RES$'\n'" | wc -l`

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

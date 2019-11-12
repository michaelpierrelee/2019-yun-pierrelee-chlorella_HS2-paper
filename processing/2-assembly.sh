#!/bin/bash

cd /media/hdd/pierrelee/korea_hs2;

# get the list of samples (with two pairs)

P1=()
P2=()

for FILE in `ls fastq/*.fq`; do
	# get sequence number
	PAIR=${FILE%.*}
	PAIR=${PAIR##*_}
	# add element to array
	if [ $PAIR = 1 ]; then P1+=("$FILE")
	else P2+=("$FILE")
	fi
done

P1=$( IFS=$','; echo "${P1[*]}" ) #https://superuser.com/a/462400
P2=$( IFS=$','; echo "${P2[*]}" )

# assembly

Trinity --seqType fq  \
	--left $P1 \
	--right $P2 \
	--output trinity_assembly \
	--CPU 16 --max_memory 28G --monitoring

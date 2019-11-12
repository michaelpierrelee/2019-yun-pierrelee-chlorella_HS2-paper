#!/bin/bash

cd /media/hdd/pierrelee/korea_hs2;

# read representation

bowtie2-build trinity_assembly/Trinity.fasta trinity_assembly/Trinity.fasta --threads 16

for FILE in `ls fastq/*.fq`
	do
	PAIR=${FILE%.*}
	PAIR=${PAIR##*_}
	if [ $PAIR = 1 ]; then
		SAMPLE=${FILE%_*}
		SAMPLE=${SAMPLE##*/}

		# read alignment statistics
		bowtie2 -p 16 -q --no-unal -k 20 -x trinity_assembly/Trinity.fasta \
			-1 fastq/${SAMPLE}_1.fq \
			-2 fastq/${SAMPLE}_2.fq \
			2>bowtie2_out/${SAMPLE}_align_stats.txt| samtools view -@ 15 -m 2G -Sb -o bowtie2_out/${SAMPLE}_bowtie2.bam

		cat 2>&1 bowtie2_out/${SAMPLE}_align_stats.txt

	fi
done

wait

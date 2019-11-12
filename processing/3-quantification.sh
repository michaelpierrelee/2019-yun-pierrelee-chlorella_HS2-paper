#!/bin/bash

cd /media/hdd/pierrelee/korea_hs2;

# quantification for each file

for FILE in `ls fastq/*.fq`
	do
	PAIR=${FILE%.*}
	PAIR=${PAIR##*_}
	if [ $PAIR = 1 ]; then
		SAMPLE=${FILE%_*}
		SAMPLE=${SAMPLE##*/}
		~/miniconda3/envs/hs2_transcriptome_env/opt/TRINITY_HOME/util/align_and_estimate_abundance.pl \
			--seqType fq \
			--left fastq/${SAMPLE}_1.fq \
			--right fastq/${SAMPLE}_2.fq \
			--transcripts trinity_assembly/Trinity.fasta  \
			--est_method RSEM \
			--aln_method bowtie \
			--gene_trans_map trinity_assembly/Trinity.fasta.gene_trans_map \
			--prep_reference \
			--coordsort_bam \
			--output_dir rsem_quantification/${SAMPLE} \
			--thread_count 16
	fi
done

wait

# gene count matrix

~/miniconda3/envs/hs2_transcriptome_env/opt/TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method RSEM --cross_sample_norm none --gene_trans_map trinity_assembly/Trinity.fasta.gene_trans_map --out_prefix rsem_quantification/Trinity_genes `ls rsem_quantification/*/RSEM.genes.results` --name_sample_by_basedir

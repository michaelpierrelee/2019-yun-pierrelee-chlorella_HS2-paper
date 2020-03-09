# Methods for the data analysis

---
Michaël Pierrelée, Aix Marseille Univ, CNRS, IBDM, UMR7288, FRANCE - <michael.pierrelee@univ-amu.fr>

*Apache License 2.0*

---

This file describes the steps to reproduce the data analysis. **All scripts and log files are available on: <https://gitlab.com/habermann_lab/yun2019_hs2>.**

More information can be found at the beginning of each notebook, where the header describes precisely the workflow. The reader is invited to access them through the GitLab repository.

[TOC]

## Setup working environment

Data processing, as well as analysis with Python and R, was performed within conda environments.

### Install conda

* Download latest conda version https://conda.io/miniconda.html.
* Start conda installation `bash Miniconda3-latest-Linux-x86_64.sh`.
* Configure conda:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set auto_activate_base False
```

### Install Python environment

```bash
conda create -n hs2_python_env
conda activate hs2_python_env
conda install python ipykernel scipy numpy pandas matplotlib seaborn xlrd openpyxl scikit-learn pydot biopython matplotlib-venn gseapy networkx
```

### Install R environment

DESeq2 v.1.20.0.

```bash
conda create -n hs2_r_env
conda activate hs2_r_env
conda install r r-irkernel r-devtools bioconductor-deseq2
```

Then, to setup SVCD (`cdnormbio` package v0.1.0) on R:

```
devtools::install_github("https://github.com/carlosproca/cdnormbio.git", dependencies = FALSE)
```

### Install data processing environment

```bash
conda create -n hs2_transcriptome_env
conda activate hs2_transcriptome_env
conda install fastqc=0.11.8 multiqc=1.7 trinity=2.8.5 rsem=1.3.2 busco=3.0.2 bowtie2=2.3.5 blast=2.9.0 samtools=1.9 bowtie=1.2.2
conda install trinotate=3.2.0 sqlite transdecoder=5.5.0 hmmer=3.2.1 diamond=0.8.36
```

## Perform de novo assembly

Assembly was performed following recommendations from:
* https://github.com/trinityrnaseq/trinityrnaseq/wiki
* https://github.com/trinityrnaseq/EMBOtrinityWorkshopSept2016/wiki/De-novo-Assembly,-Quantitation,-and-Differential-Expression

BUSCO dataset (pre-release Chlorophyta odb10) was retrieved from:

* https://busco.ezlab.org/datasets/prerelease/chlorophyta_odb10.tar.gz.

```bash
cd /media/hdd/pierrelee/korea_hs2
conda activate hs2_transcriptome_env

# quality control
bash 1-raw_qc.sh

# assembly
nohup bash 2-assembly.sh >& trinity.log &

# trinity statistics
~/miniconda3/envs/hs2_transcriptome_env/opt/TRINITY_HOME/util/TrinityStats.pl trinity_assembly/Trinity.fasta > trinity_stats.txt

# gene count matrix
nohup bash 3-quantification.sh >& rsem.log &

# BUSCO results
nohup bash 4-busco.sh >& busco.log &

# bowtie2 mapping statistics
nohup bash 5-read_representation.sh >& read_representation.log &
```

Please see the bash files on gitlab repository.

## Perform differential expression analysis

1. **Quality control of raw counts:** execute the *IPython* notebook `1-raw_count_quality.py.ipynb` within the conda environment `hs2_python_env`.
2. **Differential expression analysis:** `2-differential_expression.r.ipynb` within `hs2_r_env`.
3. **Analysis of DEGs and SVCD-DESeq2 workflow:** `3-DE_genes.py.ipynb` within `hs2_python_env`.

The latter notebook will generate `degs_to_annotate.fasta`, a subset of transcripts sequences used in the downstream steps.

## Perform functional annotation

Functional annotation was performed following recommendations from:

*  https://github.com/Trinotate/Trinotate.github.io/wiki/Software-installation-and-data-required
* https://littlebioinformatician.wordpress.com/ngs-protocols/trinotate/
* https://github.com/TransDecoder/TransDecoder/wiki

Annotations were retrieved from:

* <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz> (release 31-JUL-2019)
* <ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz> (release 30-AUG-2018)

### Prepare databases

```bash
cd /media/hdd/pierrelee/databases
conda activate hs2_transcriptome_env

# prepare BLAST+ database
makeblastdb -in uniprot_sprot-31JUL2019.faa -title uniprot_sprot -dbtype prot -parse_seqids -out db_uniprot_sprot-31JUL2019/uniprot_sprot

# prepare DIAMOND database
diamond makedb --in uniprot_sprot-31JUL2019.faa -d uniprot_sprot-31JUL2019

# TransDecoder - prepare pfam db
hmmpress Pfam-A.hmm

# Trinotate - Download SQLite database
cd /media/hdd/pierrelee/korea_hs2
Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
mv *sqlite* degs_annotation/trinotate_db/
rm Trinotate.sqlite
rm uniprot_sprot.pep
```

### Predict protein coding regions

```bash
cd /media/hdd/pierrelee/korea_hs2
conda activate hs2_transcriptome_env

# TransDecoder - extract the long open reading frames
nohup ~/miniconda3/envs/hs2_transcriptome_env/opt/transdecoder/TransDecoder.LongOrfs \
	-t degs_annotation/degs_to_annotate.fasta \
	-O degs_annotation/ \
	--gene_trans_map trinity_assembly/Trinity.fasta.gene_trans_map \
	-m 100 >& degs_annotation/transdecoder.longorfs.log &

# TransDecoder - BLAST homology search
nohup diamond blastp -d ../databases/uniprot_sprot-31JUL2019 \
	-q degs_annotation/longest_orfs.pep \
	-o degs_annotation/longest_orfs.blastp.outfmt6 \
	--max-target-seqs 1 --outfmt 6 --evalue 1e-5 --threads 16 \
	-t /dev/shm >& degs_annotation/diamond_blastp_longest_orfs.log &

# TransDecoder - PFAM homology search
bash multi_hmmscan_longestorfs.sh

# TransDecoder - predict the likely coding regions
nohup ~/miniconda3/envs/hs2_transcriptome_env/opt/transdecoder/TransDecoder.Predict \
	-t degs_annotation/degs_to_annotate.fasta \
	-O degs_annotation/ \
	--retain_pfam_hits degs_annotation/longest_orfs.pfam.domtblout \
	--retain_blastp_hits degs_annotation/longest_orfs.blastp.outfmt6 >& degs_annotation/transdecoder.predict.log &
```

### Annotate sequences

```bash
cd /media/hdd/pierrelee/korea_hs2
conda activate hs2_transcriptome_env

# SignalP v5.0b (batch 200,000 for 30 Gb memory)
nohup signalp -fasta degs_annotation/degs_to_annotate.fasta.transdecoder.pep \
	-batch 200000 -format short -org euk \
	-prefix degs_annotation/signalp.out >& degs_annotation/signalp.log &

# Trinotate - Capturing BLASTX Homologies
nohup diamond blastx -d ../databases/uniprot_sprot-31JUL2019 \
	-q degs_annotation/degs_to_annotate.fasta \
	-o degs_annotation/degs_to_annotate.blastx.outfmt6 \
	--threads 16 --max-target-seqs 1 --outfmt 6 --evalue 1e-5 >& degs_annotation/diamond_blastx_to_annotate.fasta.log &

# Trinotate - Capturing BLASTP Homologies
nohup diamond blastp -d ../databases/uniprot_sprot-31JUL2019 \
	-q degs_annotation/degs_to_annotate.fasta.transdecoder.pep \
	-o degs_annotation/transdecoder.pep.blastp.outfmt6 \
	--threads 16 --max-target-seqs 1 --outfmt 6 --evalue 1e-5 >& degs_annotation/diamond_blastp_transdecoder.pep.log &

# Trinotate - Running HMMER to identify protein domains 
bash multi_hmmscan_trinotate.sh
```

### Summarize annotations

```bash
cd /media/hdd/pierrelee/korea_hs2
conda activate hs2_transcriptome_env

# Trinotate - Load transcripts and coding regions
nohup Trinotate degs_annotation/trinotate_db/Trinotate.sqlite init \
	--gene_trans_map trinity_assembly/Trinity.fasta.gene_trans_map \
	--transcript_fasta degs_annotation/degs_to_annotate.fasta \
	--transdecoder_pep degs_annotation/degs_to_annotate.fasta.transdecoder.pep >& degs_annotation/trinotate_load_transcripts.log &

# Trinotate - Load BLASTX homologies
Trinotate degs_annotation/trinotate_db/Trinotate.sqlite LOAD_swissprot_blastx degs_annotation/degs_to_annotate.blastx.outfmt6

# Trinotate - Load BLASTP homologies
Trinotate degs_annotation/trinotate_db/Trinotate.sqlite LOAD_swissprot_blastp degs_annotation/transdecoder.pep.blastp.outfmt6 

# Trinotate - Load Pfam domain entries
Trinotate degs_annotation/trinotate_db/Trinotate.sqlite LOAD_pfam degs_annotation/trinotate.pfam.domtblout

# Trinotate - Load signal peptide predictions
Trinotate degs_annotation/trinotate_db/Trinotate.sqlite LOAD_signalp degs_annotation/signalp.out_summary.signalp5

# Trinotate - execution
Trinotate degs_annotation/trinotate_db/Trinotate.sqlite report -E 1e-10 > degs_annotation/trinotate_annotation_report.xls
```

Only BLAST hits with e-value lower than 1e-10 were kept in the downstream steps.

## Data exploration

4. **Aggregation of annotations to get dysregulated KEGG orthologies:** `4-aggregation.py.ipynb` within `hs2_python_env`.
5. **Functional enrichment of KEGG PATHWAY and BRITE databases:** `5-KEGG_enrichment.py.ipynb` within `hs2_python_env`.
   1. GSEAPY v0.9.15. KEGG database release 92.0 (01-OCT 2019).
6. **Extraction of Bowtie2 and MultiQC results and comparison of KEGG and Blast2GO pipelines results:** `6-pipeline_comparison.py.ipynb` within `hs2_python_env`.

## Blast2GO pipeline

We searched to compare the results of this above pipeline to one based on *SuperTranscript*, to aggregate transcript sequences at the gene-level, and on *Blast2GO*, to perform Blast alignment and get the GO annotation for all genes. Contrary to the above pipeline, annotation was performed on all genes with a Blast hit, and not just on a subset of differentially expressed genes. It is thus possible to have an order of magnitude on the number of DEGs used for the data exploration among all genes.

### Execute *SuperTranscript*

```bash
nohup ~/miniconda3/envs/hs2_transcriptome_env/opt/TRINITY_HOME/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
	--trinity_fasta trinity_assembly/Trinity.fasta \
	--out_prefix supertranscripts/trinity_genes >& supertranscripts/supertranscripts.log &
```

### Execute BLAST2GO

1. Download the software from https://www.blast2go.com/blast2go-pro/download-b2g (free v5.2.5).
2. Load `supertranscripts/trinity_genes.fasta`.
3. Execute Blast against `db_uniprot_sprot-31JUL2019/uniprot_sprot`. E-value threshold = 1e.10.
4. Execute GO mapping against OBO database (vOCT2019).
5. Export the result table into the file `blast2go_go_table.txt`.

### Comparison of pipeline results

Then, execute the following *IPython* notebook: `6-pipeline_comparison.py.ipynb` within `hs2_python_env`.

### Functional enrichment

See `5-GO_enrichment.py.ipynb` within `hs2_python_env`. The OBO file of gene ontologies has to be downloaded from http://release.geneontology.org/2019-10-07/ontology/go-basic.obo or, for the latest release, http://purl.obolibrary.org/obo/go/go-basic.obo.

* GSEAPY v0.9.15. Gene-ontology OBO file (release 07-OCT-2019).
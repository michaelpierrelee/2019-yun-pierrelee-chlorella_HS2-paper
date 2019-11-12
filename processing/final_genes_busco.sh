#!/bin/bash

cd /media/hdd/pierrelee/korea_hs2;

run_busco -i final_genes_busco/final_used_genes.fasta -o final_genes_busco -l ../databases/busco-prerelease_chlorophyta_odb10 -m transcriptome -c 16 -e 1e-10 --long -f;

generate_plot -wd final_genes_busco/;

#!/bin/bash

cd /media/hdd/pierrelee/korea_hs2;

run_busco -i trinity_assembly/Trinity.fasta -o busco_hs2_trinity -l ../databases/busco-prerelease_chlorophyta_odb10 -m transcriptome -c 16 -e 1e-10 --long -f;

generate_plot -wd run_busco_hs2_trinity/;

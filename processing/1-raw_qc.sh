#!/bin/bash

cd /media/hdd/pierrelee/korea_hs2;

# 1st fastqc quality control

ls fastq/*.fq | parallel -j0 --progress --workdir . 'fastqc {} -o quality_control/before_trimming'

wait

# 1st multiqc

multiqc quality_control/before_trimming/ -o quality_control/before_trimming;



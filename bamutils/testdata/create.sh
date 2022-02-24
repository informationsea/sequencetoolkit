#!/bin/bash

set -xeu -o pipefail

BWA=./bin/bwa-mac

"$BWA" index ref/MT.fa
"$BWA" mem ref/MT.fa demo1.fastq > demo1.sam
samtools sort --reference ref/MT.fa -o demo1.cram -O CRAM demo1.sam
samtools index demo1.cram
samtools sort --reference ref/MT.fa -o demo1.bam -O BAM demo1.sam
samtools index demo1.bam

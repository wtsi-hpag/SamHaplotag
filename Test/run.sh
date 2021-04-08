#!/usr/bin/env bash

set -eo pipefail

samtools view -h ../Test/test.cram | ./SamHaplotag prefix | samtools view -o test_out.cram
samtools fastq -T BX test_out.cram | ./10xSpoof prefix_SamHaplotag_Clear_BC prefix >test.fq

for file in prefix_SamHaplotag_Clear_BC prefix_SamHaplotag_Missing_BC_QT_tags prefix_SamHaplotag_UnClear_BC test.fq; do diff $file ../Test/Expected/$file; done 
diff <(samtools view test_out.cram) <(samtools view ../Test/Expected/test_out.cram)

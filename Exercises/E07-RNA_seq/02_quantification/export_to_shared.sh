#!/usr/bin/env bash

mkdir /data/shared/AGE2020/Exercises/E07-RNA_seq/02_quantification

cp -r _rmd_images 02_quantification.Rmd 02_quantification.html /data/shared/AGE2020/Exercises/E07-RNA_seq/02_quantification
cp -r ../_bin /data/shared/AGE2020/Exercises/E07-RNA_seq
cp -r ../_data/{sample_sheet_airway.csv,genome,rna_sequences,transcriptome} /data/shared/AGE2020/Exercises/E07-RNA_seq/_data

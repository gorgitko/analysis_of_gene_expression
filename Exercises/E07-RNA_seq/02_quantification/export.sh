#!/usr/bin/env bash

##-- Export to directory specified in the first argument.

cd "$(dirname "$0")"

TO_EXPORT=`echo _rmd_images/ 02_quantification.{Rmd,html} *_out multiqc_quantification sortmerna_*`

if [ -z "$1" ]; then
  echo "Destination path is not set!"
  exit 1
fi

echo "To export: $TO_EXPORT"
mkdir $1
cp -r $TO_EXPORT $1

TO_EXPORT=`echo ../data/{experiment_trimmed_filtered/,genome/,gsnap_index_hg37_chr20/,rna_sequences/,salmon_index_hg37_chr20/,sample_sheet_airway.csv,transcriptome/}`
echo "To export: $TO_EXPORT"
mkdir $1/../data
cp -r $TO_EXPORT $1/../data

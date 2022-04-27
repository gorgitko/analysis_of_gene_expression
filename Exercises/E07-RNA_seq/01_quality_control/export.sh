#!/usr/bin/env bash

##-- Export to directory specified in the first argument.

cd "$(dirname "$0")"

TO_EXPORT=`echo 01_quality_control.{Rmd,html} adapters.* fasta_to_table.py _rmd_images/`

if [ -z "$1" ]; then
  echo "Destination path is not set!"
  exit 1
fi

echo "To export: $TO_EXPORT"
mkdir $1
cp -r $TO_EXPORT $1

mkdir -p $1/../data/experiment
TO_EXPORT=`echo ../data/experiment/*_subset.fastq`
echo "To export: $TO_EXPORT"
cp -r $TO_EXPORT $1/../data/experiment

TO_EXPORT=`echo ../data/experiment_trimmed`
echo "To export: $TO_EXPORT"
cp -r $TO_EXPORT $1/../data

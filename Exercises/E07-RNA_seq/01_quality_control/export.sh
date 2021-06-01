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
mkdir $1/../data
cp -r ../data/experiment/ ../data/experiment_trimmed/ $1/../data

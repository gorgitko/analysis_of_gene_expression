#!/usr/bin/env bash

##-- Export to directory specified in the first argument.

cd "$(dirname "$0")"

TO_EXPORT="MXRA8.gff3 README.txt SRR1039508.bam SRR1039508.bam.bai"

if [ -z "$1" ]; then
  echo "Destination path is not set!"
  exit 1
fi

echo "To export: $TO_EXPORT"
mkdir $1
cp -r $TO_EXPORT $1

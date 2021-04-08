#!/usr/bin/env bash

cd "$(dirname "$0")"

TO_EXPORT=`echo data/ fit.Rds html_reports/ images/ microarrays.{Rmd,html} norm_data.Rds`

if [ -z "$1" ]; then
  echo "Destination path is not set!"
  exit 1
fi

echo "To export: $TO_EXPORT"
mkdir $1
cp -r $TO_EXPORT $1

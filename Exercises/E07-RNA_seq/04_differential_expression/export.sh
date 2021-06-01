#!/usr/bin/env bash

##-- Export to directory specified in the first argument.

cd "$(dirname "$0")"

TO_EXPORT="_rmd_images/ 04_differential_expression.html 04_differential_expression.Rmd dex_trt_vs_untrt_degs.csv"

if [ -z "$1" ]; then
  echo "Destination path is not set!"
  exit 1
fi

echo "To export: $TO_EXPORT"
mkdir $1
cp -r $TO_EXPORT $1

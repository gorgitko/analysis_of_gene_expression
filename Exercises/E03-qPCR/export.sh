#!/usr/bin/env bash

##-- Export to directory specified in the first argument.

cd "$(dirname "$0")"

TO_EXPORT="data/ _rmd_images/ qPCR.Rmd qPCR.html"

if [ -z "$1" ]; then
  echo "Destination path is not set!"
  exit 1
fi

echo "To export: $TO_EXPORT"
mkdir $1
cp -r --parents $TO_EXPORT $1
cp ../age_library_empty.R $1/../age_library.R

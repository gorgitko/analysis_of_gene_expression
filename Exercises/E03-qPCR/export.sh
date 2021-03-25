#!/usr/bin/env bash

cd "$(dirname "$0")"

TO_EXPORT=`echo data/ images/ qPCR.{Rmd,html}`

if [ -z "$1" ]; then
  echo "Destination path is not set!"
  exit 1
fi

echo "To export: $TO_EXPORT"
mkdir $1
cp -r --parents $TO_EXPORT $1
cp ../age_library_empty.R $1/../age_library.R

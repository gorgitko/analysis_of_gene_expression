#!/usr/bin/env bash

##-- Export to directory specified in the first argument.

cd "$(dirname "$0")"

TO_EXPORT=`echo assignment/E02-intro_to_R_assignment.{html,Rmd} assignment/{p_heatmap.png,p2.png} ComplexHeatmap-intro _rmd_images/ ggplot2_example.* intro_to_R.{html,Rmd} intro_to_Rmd.{html,knit.md,Rmd,utf8.md} intro_to_Rmd_tufte.html p_plotly.html tab_design.*`

if [ -z "$1" ]; then
  echo "Destination path is not set!"
  exit 1
fi

echo "To export: $TO_EXPORT"
cp -r --parents $TO_EXPORT $1

#!/usr/bin/env bash

SOURCE_FILE=$1

if [ -z ${SOURCE_FILE+x} ]; then
  echo "Parameter 1 (SOURCE_FILE) is not set."
  exit 1
fi

if [ ! -f $SOURCE_FILE ]; then
  echo "File $SOURCE_FILE doesn't exist."
  exit 1
fi

ASSIGNMENT_FILE=${SOURCE_FILE/_source\.Rmd/.Rmd}
sed "/^#S$/,/^#E$/c\#S" $SOURCE_FILE > $ASSIGNMENT_FILE
sed -i "s/^#S$//g" $ASSIGNMENT_FILE
sed -i '/^#SC$/,/^#EC$/d' $ASSIGNMENT_FILE

SOLUTION_FILE=${SOURCE_FILE/_source\.Rmd/_solutions.Rmd}
sed -E -e "s/^#(S|E)$//g" -e "/^#(SC|EC)$/d" $SOURCE_FILE > $SOLUTION_FILE

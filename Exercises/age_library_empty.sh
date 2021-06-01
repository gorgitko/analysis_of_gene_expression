#!/usr/bin/env bash

# Remove all content between lines starting with #SC and #EC (including those lines).
sed '/^#SC$/,/^#EC$/d' age_library.R > age_library_empty.R

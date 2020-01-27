# Analysis of Gene Expression @ University of Chemistry and Technology in Prague

## [E01](Exercises/E01-intro/intro.html) - Intro ([Rmd](Exercises/E01-intro/intro.Rmd))

- `sshfs` - mount directory on remote server
- Conda package manager
- Programming in R
  - Links to base R tutorials and other useful stuff.
  - Debugging.
  - Introduction to RMarkdown.
    - [Tab design](Exercises/E01-intro/tab_design.html) ([Rmd](Exercises/E01-intro/tab_design.Rmd)).

## [E02](Exercises/E02-intro_to_advanced_R/intro_to_advanced_R.html) - Intro to advanced R ([Rmd](Exercises/E02-intro_to_advanced_R/intro_to_advanced_R.Rmd))

- Introduction to tidyverse
  - `magrittr` - pipe operator
  - `tibble` - enhanced data.frame
  - `dplyr` - data manipulation
  - `tidyr` - tools for tidy data
  - `stringr` - consistent wrappers for common string operations
- `ggplot2`
  - Basic philosophy.
  - Libraries extending the `ggplot2`.
  - Additional themes.
- Other useful libraries
  - `janitor` - table summaries (RIP `table()`)
  - `plotly` - interactive HTML plots
  - `heatmaply` - interactive HTML heatmaps
  - `pheatmap` - pretty heatmaps in base R
- Parallelization

## [E03](Exercises/E03-qPCR/qPCR.html) - qPCR ([Rmd](Exercises/E03-qPCR/qPCR.Rmd))

- Main purpose of this exercise is to practice basic R on a small dataset and to
  implement a basic set of (mainly visualization) functions, which will be used
  later for microarray and RNA-seq data.
- Implemented functions are located in [age_library.R](Exercises/age_library.R),
  skeletons are in [age_library_empty.R](Exercises/age_library_empty.R).

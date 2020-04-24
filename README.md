# Analysis of Gene Expression @ University of Chemistry and Technology in Prague

## [E01](Exercises/E01-intro/intro.html) - Intro ([Rmd](Exercises/E01-intro/intro.Rmd))

- `sshfs` - mount directory on remote server
- `tmux` - termimal multiplexer
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

## [E04](Exercises/E04-microarrays/microarrays.html) - microarrays ([Rmd](Exercises/E04-microarrays/microarrays.Rmd))

- Exercise on Affymetrix microarray analysis.
- Reading in data, technical and biological quality control, normalization, differential expression, reporting.

## [E05](Exercises/E05-multiple_testing_issue/multiple_testing_issue.html) - multiple testing issue ([Rmd](Exercises/E05-multiple_testing_issue/multiple_testing_issue.Rmd))

- Demonstration of multiple testing issue correction methods on fair/skewed coins.

## [E06](Exercises/E06-IGV) - IGV browser

- Files for practising [IGV](http://software.broadinstitute.org/software/igv/) usage.

## [E07](Exercises/E07-RNA_seq) - RNA-seq

### [01](Exercises/E07-RNA_seq/01_quality_control/01_quality_control.html) - technical quality control and trimming ([Rmd](Exercises/E07-RNA_seq/01_quality_control/01_quality_control.Rmd))

- Downloading from SRA (`fasterq-dump`), `FastQC`, `MultiQC`, `Trimmomatic`.

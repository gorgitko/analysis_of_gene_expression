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

- This exercise is using [experimental data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778) from human airway smooth muscle cells treatment,
and is largely based on a great tutorial
[RNA-Seq workflow: gene-level exploratory analysis and differential expression](https://f1000research.com/articles/4-1070/v2),
from which preprocessed R data are used.

### [01](Exercises/E07-RNA_seq/01_quality_control/01_quality_control.html) - technical quality control and trimming ([Rmd](Exercises/E07-RNA_seq/01_quality_control/01_quality_control.Rmd))

- Downloading from SRA (`fasterq-dump`).
- Technical quality control (`FastQC`, `MultiQC`).
- Read trimming (`Trimmomatic`).

### [02](Exercises/E07-RNA_seq/02_quantification/02_quantification.html) - quantification ([Rmd](Exercises/E07-RNA_seq/02_quantification/02_quantification.Rmd))

- Downloading reference files (genome, annotation, etc.).
- Filtering out rRNA and tRNA (`SortMeRNA`).
- Two quantification pipelines:
  - Aligning to genome (`GSNAP`), quality control of the alignment (`RSeQC`, `preseq`) and counting overlaps (`featureCounts`).
  - Mapping to transcriptome (`Salmon`).
- Importing count matrix to R (`tximport`, `DESeq2`).
- Using `DESeqDataSet`.

### [03](Exercises/E07-RNA_seq/03_exploratory_analysis/03_exploratory_analysis.html) - exploratory analysis ([Rmd](Exercises/E07-RNA_seq/03_exploratory_analysis/03_exploratory_analysis.Rmd))

- Running `DESeq2`.
- Gene annotation.
- Count transformations, TPM calculation.
- PCA, hierarchical clustering, boxplots.

### [04](Exercises/E07-RNA_seq/04_differential_expression/04_differential_expression.html) - differential expression ([Rmd](Exercises/E07-RNA_seq/04_differential_expression/04_differential_expression.Rmd))

- Using `DESeq2` - contrasts, interactions, independent filtering, LFC shrinkage.
- Reporting results: MA plot, volcano plot, boxplots, `ReportingTools`.

### [05](Exercises/E07-RNA_seq/05_gene_set_analysis/05_gene_set_analysis.html) - Gene Set Enrichment Analysis (ORA, GSEA, SPIA) ([Rmd](Exercises/E07-RNA_seq/05_gene_set_analysis/05_gene_set_analysis.Rmd))

- Gene set databases.
- Data preparation.
- ORA (`goseq`).
- GSEA by Subramanian (`clusterProfiler`) + visualization.
- Signaling pathway impact analysis (`SPIA`).
- Viewing data in KEGG (`pathview`).
- Online tools.

## [E08](Exercises/E08-scRNA_seq/scRNA_seq.html) - single-cell RNA-seq ([Rmd](Exercises/E08-scRNA_seq/scRNA_seq.Rmd))

- Introduction, software overview, and links to tutorials, lists and other readings.

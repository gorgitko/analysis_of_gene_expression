# Analysis of Gene Expression @ University of Chemistry and Technology in Prague

These materials are for the course [Analysis of Gene Expression](https://student.vscht.cz/eng/predmety/index.php?do=predmet&kod=M143004)
taught at the [University of Chemistry and Technology in Prague](https://www.vscht.cz//?jazyk=en),
and guaranteed by [Department of Informatics and Chemistry](https://lich.vscht.cz/?jazyk=en)
in the study programme [Bioinformatics](http://studuj.bioinformatiku.cz/) (available for the bachelor, master, and PhD. levels).

The authors are from the [Laboratory of Genomics and Bioinformatics](https://www.img.cas.cz/research/michal-kolar/) at the
Institute of Molecular Genetics of the Czech Academy of Sciences:

- Michal Kolar \<kolarmi@img.cas.cz\> (guarantor, theoretical lectures)
- Jiri Novotny \<jiri.novotny@img.cas.cz\> (exercises)

In case of suggestions or problems, create a new issue.
We will be happy to answer your questions, integrate new ideas, or resolve any problems :blush:

# Lectures

Recordings and materials for theoretical lectures are stored at school MS Teams,
and currently available only to course participants.

# Exercises

## Prerequisites

We expect all participants to have a basic knowledge of base R and Linux shell (bash).
Links to relevant materials can be found in E01 - Intro.

## Software prerequisites

We are using virtual machines (VMs) with images based on Debian 10 and including all the necessary software
(R 4.0, RStudio Server, `conda`, and various tools).
We gratefully thank to the [Metacentrum Cloud](https://www.metacentrum.cz/en/Sluzby/Cloud/index.html)
team for a great assistance with virtual machines :heart:

However, it is possible to install all the stuff in order to have the same environment as our VMs offer (or be very close to it).
Generally, we recommend to work on Linux-based system (our tip: [Linux Mint](https://linuxmint.com/)).

### Getting the exercise files

Just download and unzip this repository.
Additional data files for E07 - RNA-seq must be downloaded, see the relevant section below.

### R dependencies

You need [R](https://www.r-project.org/) 4.0+ and [Bioconductor](http://www.bioconductor.org/install/) 3.12+ installed.
We recommend to use [RStudio IDE](https://www.rstudio.com/) for programming.

A lockfile for [renv](https://rstudio.github.io/renv/articles/renv.html) is included -
it captures all packages needed to run the exercises. Moreover, `renv` ensures all packages
are installed to a local R library, and thus, the installation doesn't pollute the system library.

To start the installation of required packages:

1. Create a new RStudio project in `Exercises/` directory. If you are not using RStudio, just change R's working directory to `Exercises/`.
2. Start R.
3. Run `renv::init()`. This will create a new project-specific library and install packages from `renv.lock`.
   If `renv` is not available, install it first by `install.packages("renv")`.

### Other tools

Other tools could be installed through your OS package manager or the `conda` tool (see E01 - Intro).
The latter is recommended for bioinformatics tools, which are mainly used during RNA-seq exercises.

***

## [E01](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E01-intro/intro.html) - Intro ([Rmd](Exercises/E01-intro/intro.Rmd)) - _Jiri Novotny_

- Some information about our virtual machines and files.
- `sshfs` - mount directory on a remote server
- `tmux` - termimal multiplexer
- `fish` - a friendly, interactive shell
- `conda` - package and virtual environment manager
- Links to beginner base R tutorials and other useful stuff.

## [E02](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E02-intro_to_R/intro_to_R.html) - Intro to R ([Rmd](Exercises/E02-intro_to_R/intro_to_R.Rmd)) - _Jiri Novotny_

- [Introduction to RMarkdown](https://raw.githubusercontent.com/gorgitko/analysis_of_gene_expression/master/Exercises/E02-intro_to_R/intro_to_Rmd.html)
  ([Rmd](Exercises/E02-intro_to_R/intro_to_R.Rmd)).
- Reproducible R (project-oriented workflow, consistent paths using [here()](https://here.r-lib.org/),
  namespace conflicts, [renv](https://rstudio.github.io/renv/articles/renv.html), etc.).
- Installing R packages.
- Debugging R.
- Writing your own functions.
- Vectorized operations, avoiding for loops, parallelization.
- Introduction to tidyverse
  - Overview of tidy data and non-standard/tidy evaluation.
  - `magrittr` - pipe operator
  - `tibble` - enhanced data.frame
  - `dplyr` - data manipulation
  - `tidyr` - tools for tidy data
  - `stringr` - consistent wrappers for common string operations
  - `glue` - string interpolation
  - `purrr` - functional programming tools
  - `ggplot2`
    - Basic philosophy and usage.
    - Libraries extending the `ggplot2`.
    - Additional themes.
- Other useful libraries
  - `janitor` - table summaries
  - `plotly` - interactive HTML plots
  - `heatmaply` - interactive HTML heatmaps
  - `pheatmap` - pretty heatmaps in base R
  - `ComplexHeatmap` -
    [introduction](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E02-intro_to_R/ComplexHeatmap-intro/ComplexHeatmap.html)
    ([Rmd](Exercises/E02-intro_to_R/ComplexHeatmap-intro/ComplexHeatmap.Rmd))
  - `BiocParallel` - parallelized `lapply()` and others

## [E03](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E03-qPCR/qPCR.html) - qPCR ([Rmd](Exercises/E03-qPCR/qPCR.Rmd)) - _Jiri Novotny_

- Main purpose of this exercise is to practice basic R on a small dataset and to
  implement a basic set of (mainly visualization) functions, which will be used
  later for microarray and RNA-seq data.
- Implemented functions are located in [age_library.R](Exercises/age_library.R),
  skeletons are in [age_library_empty.R](Exercises/age_library_empty.R).

## [E04](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E04-microarrays/microarrays.html) - microarrays ([Rmd](Exercises/E04-microarrays/microarrays.Rmd)) - _Jiri Novotny_

- Exercise on Affymetrix microarray analysis.
- Reading in data, technical and biological quality control, normalization, differential expression, reporting.

## [E05](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E05-multiple_testing_issue/multiple_testing_issue.html) - multiple testing issue ([Rmd](Exercises/E05-multiple_testing_issue/multiple_testing_issue.Rmd)) - _Michal Kolar_

- Demonstration of multiple testing issue correction methods on fair/skewed coins.

## [E06](Exercises/E06-IGV) - IGV browser - _Michal Kolar_

- Files for practising [IGV](http://software.broadinstitute.org/software/igv/) usage.

## [E07](Exercises/E07-RNA_seq) - RNA-seq - _Jiri Novotny_

- This exercise is using [experimental data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)
from human airway smooth muscle cells treatment, and is largely based on a great tutorial
[RNA-Seq workflow: gene-level exploratory analysis and differential expression](https://f1000research.com/articles/4-1070/v2),
from which preprocessed R data are later used (starting from `03 - exploratory analysis` part).

Additional data files must be downloaded prior from [here](https://onco.img.cas.cz/novotnyj/age/AGE2021_data.tar.gz).
If you are working on a remote server, you can use `wget` for downloading: `wget https://onco.img.cas.cz/novotnyj/age/AGE2021_data.tar`.
Then decompress the downloaded archive to `Exercises/` directory, e.g. `tar xzf AGE2021_data.tar -C /path/to/Exercises`.

### [01](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E07-RNA_seq/01_quality_control/01_quality_control.html) - technical quality control and trimming ([Rmd](Exercises/E07-RNA_seq/01_quality_control/01_quality_control.Rmd))

- Downloading from SRA (`fasterq-dump`).
- Technical quality control (`FastQC`, `MultiQC`).
- Read trimming (`Trimmomatic`).

### [02](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E07-RNA_seq/02_quantification/02_quantification.html) - quantification ([Rmd](Exercises/E07-RNA_seq/02_quantification/02_quantification.Rmd))

- Downloading reference files (genome, annotation, etc.).
- Filtering out rRNA and tRNA (`SortMeRNA`).
- Two quantification pipelines:
  - Aligning to genome (`GSNAP`), quality control of the alignment (`RSeQC`, `preseq`) and counting overlaps (`featureCounts`).
  - Mapping to transcriptome (`Salmon`).
- Importing count matrix to R (`tximport`, `DESeq2`).
- Using `DESeqDataSet`.

### [03](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E07-RNA_seq/03_exploratory_analysis/03_exploratory_analysis.html) - exploratory analysis ([Rmd](Exercises/E07-RNA_seq/03_exploratory_analysis/03_exploratory_analysis.Rmd))

- Running `DESeq2`.
- Gene annotation.
- Count transformations, TPM calculation.
- PCA, hierarchical clustering, boxplots.

### [04](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E07-RNA_seq/04_differential_expression/04_differential_expression.html) - differential expression ([Rmd](Exercises/E07-RNA_seq/04_differential_expression/04_differential_expression.Rmd))

- Using `DESeq2` - contrasts, interactions, independent filtering, LFC shrinkage.
- Reporting results: MA plot, volcano plot, boxplots, `ReportingTools`.

### [05](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E07-RNA_seq/05_gene_set_analysis/05_gene_set_analysis.html) - Gene Set Enrichment Analysis (ORA, GSEA, SPIA) ([Rmd](Exercises/E07-RNA_seq/05_gene_set_analysis/05_gene_set_analysis.Rmd))

- Gene set databases.
- Data preparation.
- ORA (`goseq`).
- GSEA by Subramanian (`clusterProfiler`) + visualization.
- Signaling pathway impact analysis (`SPIA`).
- Viewing data in KEGG (`pathview`).
- Online tools.

## [E08](https://gitcdn.link/repo/gorgitko/analysis_of_gene_expression/master/Exercises/E08-scRNA_seq/scRNA_seq.html) - single-cell RNA-seq ([Rmd](Exercises/E08-scRNA_seq/scRNA_seq.Rmd)) - _Jiri Novotny_

- Introduction, software overview, and links to tutorials, lists and other readings.

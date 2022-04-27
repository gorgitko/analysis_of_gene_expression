library(tidyverse)

render_callr <- function(rmd_file,
                         output_file = NULL,
                         envir = new.env(),
                         knit_root_dir = here::here(),
                         output_format = NULL,
                         warning = FALSE,
                         message = FALSE,
                         eval = TRUE,
                         clean = TRUE,
                         ...) {
  cat(glue::glue("Processing {rmd_file} -> {output_file} ..."))
  knitr::opts_chunk$set(warning = warning, message = message, eval = eval, ...)
  rmd_file <- here::here(rmd_file)
  if (is.null(output_file)) {
    output_file <- fs::path(rmd_file)
    fs::path_ext(output_file) <- ".html"
  }
  output_file <- here::here(output_file)

  res <- tryCatch(
    expr = {
      knitr_msg <- capture.output(rmarkdown::render(
        rmd_file,
        output_file = output_file,
        envir = envir,
        knit_root_dir = knit_root_dir,
        output_format = output_format,
        intermediates_dir = tempfile(),
        clean = clean
      ))
      cat("Done.")
      list(error = FALSE, error_msg = "", knitr_msg = paste0(knitr_msg, collapse = "\n"))
    },
    error = function(e) {
      return(list(error = TRUE, error_msg = glue::glue("Error in {rmd_file}:\n{e$message}", .trim = FALSE), knitr_msg = ""))
    }
  )

  return(tibble::tibble(rmd_file = rmd_file, output_file = output_file, error = res$error, error_msg = res$error_msg, knitr_msg = res$knitr_msg))
}

file_defs <- list(
  ## -- E01 - intro
  list(
    rmd_file = "E01-intro/intro.Rmd",
    eval = FALSE
  ),

  ## -- E02 - intro to R and Rmd
  list(
    rmd_file = "E02-intro_to_R/intro_to_R.Rmd"
  ),
  list(
    rmd_file = "E02-intro_to_R/intro_to_Rmd.Rmd",
    clean = FALSE
  ),
  list(
    rmd_file = "E02-intro_to_R/intro_to_Rmd.Rmd",
    output_file = "E02-intro_to_R/intro_to_Rmd_tufte.html",
    output_format = "tufte::tufte_html"
  ),

  ## -- E03 - qPCR
  list(
    rmd_file = "E03-qPCR/qPCR.Rmd"
  ),

  ## -- E04 - microarrays
  list(
    rmd_file = "E04-microarrays/microarrays.Rmd"
  ),

  ## -- E05 - multiple testing issue
  list(
    rmd_file = "E05-multiple_testing_issue/multiple_testing_issue.Rmd"
  ),

  ## -- E07 - RNA-seq
  list(
    rmd_file = "E07-RNA_seq/01_quality_control/01_quality_control.Rmd",
    eval = FALSE
  ),
  list(
    rmd_file = "E07-RNA_seq/02_quantification/02_quantification.Rmd",
    eval = FALSE
  ),
  list(
    rmd_file = "E07-RNA_seq/03_exploratory_analysis/03_exploratory_analysis.Rmd"
  ),
  list(
    rmd_file = "E07-RNA_seq/04_differential_expression/04_differential_expression.Rmd"
  ),
  list(
    rmd_file = "E07-RNA_seq/05_gene_set_analysis/05_gene_set_analysis.Rmd"
  ),

  ## -- E08 - scRNA-seq
  list(
    rmd_file = "E08-scRNA_seq/scRNA_seq.Rmd"
  )
)

res <- lapply(file_defs, function(file_def) {
  message(file_def$rmd_file)

  stdout_file <- tempfile()
  stderr_file <- tempfile()
  res_render <- callr::r(
    render_callr,
    args = file_def,
    stdout = stdout_file,
    stderr = stderr_file,
    env = c(callr::rcmd_safe_env(), CALLR = "TRUE")
  )

  if (res_render$error) {
    message("  ERROR")
  }

  res_render %>%
    dplyr::mutate(
      stdout = readr::read_file(stdout_file),
      stderr = readr::read_file(stderr_file)
    )
}) %>% dplyr::bind_rows()

res_error <- res %>%
  dplyr::filter(error)

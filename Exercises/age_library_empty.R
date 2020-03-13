### These functions have the same arguments as functions implemented by me and Michal.
### Be free to use whatever arguments you want, but the output should be the same as in exercise HTML.
### Functions with matrix input should be universal, as we will use similar biological matrices later:
### rows are features (genes, probes, etc.), columns are samples, values are measurements (probe intensity, read counts, etc.).

### TASK 1: Implement your own function which takes CP matrix and returns a dendrogram with coloring by a chosen variable.
### Hints:
### - Use dendrapply() to iterate dendrogram nodes.
### - Use attributes() to get node attributes and attr() to set "nodePar" attribute.

#' Calculate hiearchical clustering and plot dendrogram.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   method_distance: Distance metric used for clustering. See ?dist.
#'                    Can be also correlation metrix ("pearson", "kendall", "spearman"). See ?cor.
#'   method_clustering: Clustering method. See ?hclust.
#'   color_by: Vector of discrete values to colorize samples.
#'             Length must be the same as number of columns in 'm'.
#'   color_by_lab: Name of the color legend.
#'   main: Main plot title.
#'
#' Returns:
#'   Object of class 'dendrogram'.
plot_hc <- function(
  m,
  method_distance = "euclidean",
  method_clustering = "complete",
  color_by = NULL,
  color_by_lab = "Group",
  main = "Hierarchical Clustering"
) {}

### TASK 2: Implement a function for PCA visualization with point coloring by chosen variable.

### We will be using this function for microarray and RNA-seq data, where we want to calculate
### PCA using top most variable genes (from few tens of thousand ones),
### that is, genes with the highest variance across samples.
### You can implement this function later.

#' Select features (rows) with the highest variance across the samples (columns).
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   n_top_features: Number of top features.
#'
#' Returns:
#'   Object 'm' subset to 'n_top_features' with the highest variance across the samples.
select_var_features <- function(m, n_top_features) {}

# You can use ggplot2, but be prepared to use data in long format.

#' Using ggplot2, plot first three PCs of samples in expression matrix.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   sample_data: Dataframe describing samples.
#'   n_top_features: Number of top features with the highest variance across the samples.
#'   color_by: Column name in sample_data to use for point coloring.
#'   shape_by: Column name in sample_data to use for point shape.
#'   label_by: Column name in sample_data to use for point labels.
#'   point_size: Point size (numeric).
#'   text_size: Label text size (numeric).
#'   center: Whether to center PCA. See ?prcomp.
#'   scale.: Whether to scale PCA. See ?prcomp.
#'
#' Returns:
#'   list(pca = prcomp object, pca_df = combined dataframe of sample_data and PCs, plot = ggplot2 object)
plot_pca_ggplot2 <- function(
  m,
  sample_data,
  n_top_features = Inf,
  color_by = NULL,
  shape_by = NULL,
  label_by = NULL,
  point_size = 2,
  text_size = 2.5,
  center = TRUE,
  scale. = TRUE
) {}

# You can also use GGally::ggpairs() to plot a grid of PCA plots (PC1 vs. PC2, PC1 vs. PC3, PC2 vs. PC3, etc).

#' Using GGally::ggpairs(), plot grid of PCA plots (PC1 vs. PC2, PC1 vs. PC3, PC2 vs. PC3, etc).
#' When n_components == 2, use normal ggplot2.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   sample_data: Dataframe describing samples.
#'   output_file: File to save plot in.
#'   n_components: Number of PCs to plot.
#'   n_top_features: Number of top features with the highest variance across the samples.
#'   color_by: Column name in sample_data to use for point coloring.
#'   color_legend_lab: Name of the color legend.
#'   shape_by: Column name in sample_data to use for point shape.
#'   shape_legend_lab: Name of the shape legend.
#'   label_by: Column name in sample_data to use for point labels.
#'   point_size: Point size (numeric).
#'   text_size: Label text size (numeric).
#'   title: Plot title.
#'   subtitle: Plot subtitle.
#'   center: Whether to center PCA. See ?prcomp.
#'   scale.: Whether to scale PCA. See ?prcomp.
#'
#' Returns:
#'   list(pca = prcomp object, pca_df = combined dataframe of sample_data and PCs, plot = ggplot2 object)
plot_pca_ggpairs <- function(
  m,
  sample_info,
  output_file = NULL,
  n_components = 5,
  n_top_features = Inf,
  color_by = NULL,
  color_legend_lab = NULL,
  shape_by = NULL,
  shape_legend_lab = NULL,
  label_by = NULL,
  point_size = 2,
  text_size = 2.5,
  title = NULL,
  subtitle = NULL,
  center = TRUE,
  scale. = TRUE
) {}

### TASK 3: Implement your own function for heatmap.

#' Plot heatmap using heatmap.2()
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   color_by: Vector of discrete values to colorize samples.
#'             Length must be the same as number of columns in 'm'.
#'   color_legend_lab: Name of the color legend.
#'   main: Plot title.
#'
#' Returns:
#'   logical
plot_heatmap <- function(
  m,
  color_by = NULL,
  color_legend_lab = "Group",
  main = ""
) {}

# pheatmap is a nicer and easier to use version of heatmap.2()

#' Plot heatmap using pheatmap()
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   sample_data: Dataframe describing samples.
#'   feature_data: Dataframe describing features
#'   column_color_by: Vector of column names of sample_data to color columns in heatmap.
#'                    If NULL and sample_data is provided, all columns will be used for coloring.
#'   row_color_by: Vector of column names of feature_data to color rows in heatmap.
#'                 If NULL and feature_data is provided, all columns will be used for coloring.
#'   main: Plot title.
#'
#' Returns:
#'   pheatmap object
plot_pheatmap <- function(
  m,
  sample_data = NULL,
  feature_data = NULL,
  column_color_by = NULL,
  row_color_by = NULL,
  main = ""
) {}

# You can use a nice library heatmaply which produces an interactive HTML heatmap.

#' Create heatmap using heatmaply.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   sample_data: Dataframe describing samples.
#'   feature_data: Dataframe describing features
#'   column_color_by: Vector of column names of sample_data to color columns in heatmap.
#'                    If NULL and sample_data is provided, all columns will be used for coloring.
#'   row_color_by: Vector of column names of feature_data to color rows in heatmap.
#'                 If NULL and feature_data is provided, all columns will be used for coloring.
#'   main: Plot title.
#'   key.title: Main color legend title.
#'
#' Returns:
#'   plotly object
plot_heatmaply <- function(
  m,
  sample_data = NULL,
  feature_data = NULL,
  column_color_by = NULL,
  row_color_by = NULL,
  main = NULL,
  key.title = NULL
) {}

### TASK 4: Implement a function which computes the M-value.
### Hints:
### - Remember that R is vectorized by default. In this way you can directly write most of equations as you see them on paper.
### - Very useful is `apply()` function. It applies a function to each of rows or columns of a matrix.
###   Whether to apply function to rows or columns can be specified by `MARGIN` parameter: `1` for rows and `2` for columns.
###   That is, when `MARGIN = 1`, you obtain N values from matrix with N rows. Similarly for `MARGIN = 2`, but for M columns.

#' Compute the M value of CP values.
#'
#' Args:
#'  gene: Name of gene to compute the M value.
#'  cp_m: Dataframe of CP values. Rows are genes and columns are samples.
#'
#' Returns:
#'  M value (numeric)
compute_m <- function(gene, cp_m) {}

### TASK 5: Implement a function to produce a boxplot of CP values.
### Hints:
### - ggpubr::ggboxplot() is pretty good.
### - Again, data in long format are required.

#' Plot boxplot of gene expression.
#'
#' Args:
#'   plot_data: data.frame (long format)
#'   x: Column to divide x-axis values to (e.g. sample groups).
#'   y: Column to compute boxplots on y-axis.
#'   feature_col: Column for facetting.
#'   feature: Name of feature from which boxplots will be made.
#'            If NULL, facet by feature_col (facet_by will be ignored).
#'   facet_by: One or two columns used for facetting.
#'   color_by: Column to use for boxplot and point coloring.
#'   x_lab: Name of x-axe.
#'   y_lab: Name of y-axe.
#'   main: Main plot title.
#'   add: Add something more to boxplots.
#'        Allowed values are one or the combination of:
#'        "none", "dotplot", "jitter", "boxplot", "point", "mean",
#'        "mean_se", "mean_sd", "mean_ci", "mean_range", "median",
#'        "median_iqr", "median_mad", "median_range".
#'        See ?ggpubr::ggboxplot
#'   point_size: Size of points inside boxplots.
#'   outlier_shape: Which point shape to use for outliers.
#'   do_t_test: Whether to do the t-test and display a p-value inside the plot.
#'
#' Returns:
#'   ggplot2 object
plot_boxplot_ggplot2 <- function(
  plot_data,
  x,
  y,
  feature_col,
  feature = NULL,
  facet_by = NULL,
  color_by = NULL,
  x_lab = x,
  y_lab = y,
  main = NULL,
  add = "jitter",
  point_size = 2,
  outlier_shape = 0,
  do_t_test = TRUE
) {}

### TASK 6: Implement a function which do the t-test of gene from two groups.
### Hint:
### - Look at the "formula" parameter of "t.test" function.

#' For a single gene, test for statistical significance of difference in group means.
#'
#' Args:
#'   gene: Name of gene to test.
#'   gene_data: Dataframe in long format.
#'   gene_col: Column with genes.
#'   value_col: Column with values to test.
#'   group_col: Column with sample groups. There must be exactly two different groups.
#'   test: Statistical test to perform. It must have the same interface as t.test()
#'   verbose: Whether to print test results.
#'
#' Returns:
#'   htest object
test_gene <- function(gene, gene_data, gene_col, value_col, group_col, test = t.test, verbose = TRUE) {}

### You can also produce a table for all genes.
### Hint:
### - Make use of test_gene() and just put results to dataframe.

#' For all genes in the input dataframe, test for statistical significance of difference in group means.
#'
#' Args:
#'   gene_data: Dataframe in long format.
#'   gene_col: Column with genes.
#'   value_col: Column with values to test.
#'   group_col: Column with sample groups. There must be exactly two different groups.
#'   test: Statistical test to perform. It must have the same interface as t.test()
#'
#' Returns:
#'   tibble object
test_gene_table <- function(gene_data, gene_col, value_col, group_col, test = t.test) {}

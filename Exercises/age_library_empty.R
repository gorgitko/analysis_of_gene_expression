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
#'   title: Main plot title.
#'
#' Returns:
#'   dendrogram object
plot_hc <- function(
  m,
  method_distance = "euclidean",
  method_clustering = "complete",
  color_by = NULL,
  color_by_lab = "Group",
  title = "Hierarchical Clustering"
) {
}

#' Select features (rows) with the highest variance across the samples (columns).
#' This will be useful for visualization (mainly heatmaps) of large expression matrices.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   n_top_features: Number of top features.
#'
#' Returns:
#'   Subset of 'm' with 'n_top_features' with the highest variance across the samples.
select_var_features <- function(m, n_top_features) {
}

#' Using the ggplot2 package, plot the first PC or first to three PCs of samples in expression matrix.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   sample_data: Dataframe describing samples.
#'   plot_type: "single" for PC1 vs. PC2, "multi" for combinations of PC1-3 and their cumulative explained variance.
#'   n_top_features: Number of top features with the highest variance across the samples.
#'   color_by: Column name in sample_data to use for point coloring.
#'   shape_by: Column noneame in sample_data to use for point shape.
#'   label_by: Column name in sample_data to use for point labels.
#'   point_size: Point size (numeric).
#'   text_size: Label text size (numeric).
#'   center: Whether to center PCA. See ?prcomp.
#'   scale.: Whether to scale PCA. See ?prcomp.
#'
#' Returns:
#'   list(pca = <prcomp object>, pca_df = <combined dataframe of sample_data and PCs>, plot = <ggplot2 or patchwork object (depends on plot_type)>)
plot_pca <- function(
  m,
  sample_data,
  plot_type = c("single", "multi"),
  n_top_features = Inf,
  color_by = NULL,
  shape_by = NULL,
  label_by = NULL,
  point_size = 2,
  text_size = 2.5,
  center = TRUE,
  scale. = TRUE
) {
}

#' Using the GGally::ggpairs() function, plot grid of PCA plots (PC1 vs. PC2, PC1 vs. PC3, PC2 vs. PC3, etc).
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
#'   list(pca = <prcomp object>, pca_df = <combined dataframe of sample_data and PCs>, plot = <ggplot2 object>)
plot_pca_ggpairs <- function(
  m,
  sample_data,
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
) {
}

#' Plot heatmap using the ComplexHeatmap package.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   z_score: If TRUE, calculate row z-score.
#'   sample_annotation: Dataframe used for annotation of columns.
#'   feature_annotation: Dataframe used for annotation of rows.
#'   title: Heatmap title.
#'   legend_title: Heatmap color legend title.
#'   show_row_names: If TRUE, show rownames in the heatmap.
#'   show_col_names: If TRUE, show colnames in the heatmap.
#'   color_palette: Function to generate colors for annotations.
#'   color_mapping: Named list of named vectors to map colors to variable levels.
#'                  See https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#simple-annotation
#'
#' Returns:
#'   ComplexHeatmap object
plot_heatmap <- function(
  m,
  z_score = FALSE,
  sample_annotation = NULL,
  feature_annotation = NULL,
  title = "",
  legend_title = "Values",
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  color_palette = scales::hue_pal(),
  color_mapping = NULL
) {
}

#' Create a heatmap using the heatmaply package.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   z_score: If TRUE, calculate row z-score.
#'   sample_annotation: Dataframe used for annotation of columns.
#'   feature_annotation: Dataframe used for annotation of rows.
#'   title: Heatmap title.
#'   legend_title: Heatmap color legend title.
#'
#' Returns:
#'   plotly object
plot_heatmaply <- function(
  m,
  z_score = FALSE,
  sample_annotation = NULL,
  feature_annotation = NULL,
  main = NULL,
  legend_title = NULL,
  showticklabels = c(TRUE, TRUE)
) {
}

#' Using the ggpubr::ggboxplot() function, plot boxplots of gene expression.
#'
#' Args:
#'   plot_data: data.frame (long format)
#'   x: Column to divide x-axis values to (e.g. sample groups).
#'   y: Column to compute boxplots on y-axis.
#'   facet_by: One or two columns used for facetting.
#'   feature: Name of a feature from which boxplots will be made.
#'            Data will be filtered based on facet_by.
#'            E.g. if facet_by = "gene" and feature = "CD24", only boxplots for "CD24" will be made.
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
#'   do_t_test: Whether to do the t-test and display p-values inside the plot.
#'
#' Returns:
#'   ggplot2 object
plot_boxplots <- function(
  plot_data,
  x,
  y,
  facet_by,
  feature = NULL,
  color_by = NULL,
  x_lab = x,
  y_lab = y,
  main = NULL,
  add = "jitter",
  point_size = 2,
  outlier_shape = 0,
  do_t_test = TRUE
) {
}

#' Compute the M value of CP values.
#'
#' Args:
#'  gene: Name of gene to compute the M value for.
#'  cp: Matrix or dataframe of CP values. Rows are genes and columns are samples.
#'
#' Returns:
#'  M value (numeric)
compute_m <- function(gene, cp) {
}

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
test_gene <- function(gene, gene_data, gene_col, value_col, group_col, test = t.test, verbose = TRUE) {
}

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
test_gene_table <- function(gene_data, gene_col, value_col, group_col, test = t.test) {
}

#' Return asterisks according to p-values.
#'
#' Args:
#'   p_value: Vector of p-values.
#'
#' Returns:
#'  character vector
asterisk <- function(p_value) {
}

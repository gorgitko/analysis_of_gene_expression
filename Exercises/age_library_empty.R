### These functions have the same arguments as functions implemented by me and Michal.
### Be free to use whatever arguments you want, but the output should be the same as in exercise HTML.

### TASK 1: Implement your own function which takes CP matrix and returns a dendrogram with coloring by chosen variable.
### Hints:
### - Use dendrapply() to iterate dendrogram nodes.
### - Use attributes() to get node attributes and attr() to set "nodePar" attribute.
plot_hc <- function(
  m,
  method_distance = "pearson",
  method_clustering = "complete",
  group = NULL,
  group_name = "Group",
  main = "Hierarchical Clustering"
) {}

### TASK 2: Implement a function for PCA visualization with point coloring by chosen variable.

### We will be using this function for microarray and RNA-seq data, where we want to calculate
### PCA using top most variable genes (from few tens of thousand ones),
### that is, genes with the highest variance across samples.
### You can implement this function later.
select_var_features <- function(m, n_top_features) {}

plot_pca <- function(
  x,
  group,
  center = TRUE,
  scale. = TRUE,
  group_name = "Group"
) {}

# You can use ggplot2, but be prepared to use data in long format.
# Hint
plot_pca_ggplot2 <- function(
  m,
  pheno_data,
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
plot_pca_ggpairs <- function() {}

### TASK 3: Implement your own function for heatmap.
plot_heatmap <- function() {}
# You can use a nice library heatmaply (https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html) which produces interactive HTML heatmaps.
plot_heatmaply <- function() {}

### TASK 4: Implement a function which computes the M-value.
compute_m <- function() {}

### TASK 5: Implement a function to produce a boxplot of CP values.
plot_boxplot_ggplot <- function() {}

### TASK 6: Implement a function which do the t-test of gene from two groups.
test_gene <- function() {}
# You can also produce a table for all genes.
test_gene_table <- function() {}

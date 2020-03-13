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
) {
  library(RColorBrewer)
  library(glue)
  library(magrittr)

  if (method_distance %in% c("pearson", "kendall", "spearman")) {
    di <- as.dist(1 - cor(m, method = method_distance))
    attributes(di) <- c(attributes(di), method = method_distance)
    xlab <- ""
  } else {
    di <- dist(t(m), method = method_distance)
    xlab <- ""
  }

  hc <- hclust(di, method = method_clustering)
  hc <- as.dendrogram(hc, hang = 0.1)

  if (!is.null(color_by)) {
    if (is.null(names(color_by))) {
      color_by <- setNames(color_by, colnames(m))
    }

    color_by <- factor(color_by)
    if (length(levels(color_by)) > 9) {
      stop(glue("The '{group_name}' factor can have at most 9 distinct levels."))
    } else {
      pal <- levels(color_by) %>% length() %>% brewer.pal("Set1")
    }
    col <- pal[color_by]
    names(col) <- names(color_by)

    set_leaf_color <- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        attr(n, "nodePar") <-
          c(a$nodePar, list(lab.cex = 0.8,
                            lab.col = col[a$label],
                            pch = NA))
      }
      return(n)
    }

    hc <- dendrapply(hc, set_leaf_color)
    plot(hc, main = main, xlab = xlab)
    legend("topright", title = color_by_lab, legend = levels(color_by), pch = 19,
           col = pal, bty = "n")
  } else {
    plot(hc, main = main, xlab = xlab)
  }

  return(invisible(hc))
}

#' Select features (rows) with the highest variance across the samples (columns).
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   n_top_features: Number of top features.
#'
#' Returns:
#'   Object 'm' subset to 'n_top_features' with the highest variance across the samples.
select_var_features <- function(m, n_top_features) {
  if (!is.infinite(n_top_features))
  {
    # Calculate the variance for each gene.
    rv <- matrixStats::rowVars(m)

    # Select the n top genes by variance.
    select <- order(rv, decreasing = TRUE)[min(n_top_features, length(rv)) %>% seq_len()]
    m <- m[select, ]
  }

  return(m)
}

#' Plots PCA of first three PC's, colorized by 'group'. Uses the base R graphics.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   color_by: Vector of discrete values to colorize samples.
#'             Length must be the same as number of columns in 'm'.
#'   n_top_features: Number of top features with the highest variance across the samples.
#'   center: Whether to center PCA. See ?prcomp.
#'   scale.: Whether to scale PCA. See ?prcomp.
#'   color_by_lab: Name of the color legend.
#'
#' Returns:
#'   Object of class 'prcomp'.
plot_pca <- function(
  m,
  color_by,
  n_top_features = Inf,
  center = TRUE,
  scale. = TRUE,
  color_by_lab = "Group"
) {
  library(RColorBrewer)
  library(glue)
  library(magrittr)

  m <- select_var_features(m, n_top_features = n_top_features)

  pca <- prcomp(t(m), center = center, scale. = scale.)

  color_by <- factor(color_by)
  if (length(levels(color_by)) > 9) {
    stop(glue("The '{group_name}' factor can have at most 9 distinct levels."))
  } else {
    pal <- levels(color_by) %>% length() %>% brewer.pal("Set1")
  }
  col <- pal[color_by]
  names(col) <- names(color_by)

  op <- par(mfrow = c(2, 2))
  plot(PC2 ~ PC1, data = pca$x, col = col, pch = 19)
  plot(PC2 ~ PC3, data = pca$x, col = col, pch = 19)
  plot(PC3 ~ PC1, data = pca$x, col = col, pch = 19)

  len <- length(pca$sdev)
  cme <- summary(pca)$importance["Cumulative Proportion", ]
  ymax <- c(0, 100)

  barplot(100 * cme, col = "grey", xlab = "PC rank",
          ylab = "Cumulative variance explained (%)", ylim = 1.5 * ymax)

  legend("topleft", title = color_by_lab, legend = levels(color_by), pch = 19,
         col = pal, bty = "n", ncol = 3, cex = 0.5)
  par(op)
  return(invisible(pca))
}

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
) {
  library(rlist)
  library(ggplot2)

  x <- select_var_features(m, n_top_features = n_top_features)

  pca <- prcomp(t(m), center = center, scale. = scale.)
  pca_df <- data.frame(sample_data, pca$x)

  cme <- summary(pca)$importance["Cumulative Proportion", ]
  cme <- as.data.frame(cme)
  cme$PC <- rownames(cme)

  pc_combinations <- list(c("PC2", "PC1"), c("PC3", "PC1"), c("PC2", "PC3"))
  p_list <- lapply(pc_combinations, function(pc) {
    p <- ggplot(pca_df, aes_string(x = pc[1], y = pc[2], color = color_by, shape = shape_by)) +
      geom_point(size = point_size) +
      ggthemes::theme_few()

    if (!is.null(label_by)) {
      p <- p + ggrepel::geom_text_repel(aes_string(label = label_by), color = "black", size = text_size)
    }

    return(p)
  })

  # Get the legend and align it properly.
  legend <- cowplot::get_legend(
    p_list[[1]] +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.margin = margin(-50, 0, 0, 200)
    )
  )

  # To use a shared legend, we need to hide other legends.
  p_list <- lapply(p_list, function(p) {
    p + theme(legend.position = "none")
  })

  p_list <- list.append(
    p_list,
    ggplot(cme, aes(x = PC, y = cme)) + geom_bar(stat = "identity") +
      ggthemes::theme_few(),
    legend
  )

  p <- cowplot::plot_grid(plotlist = p_list, nrow = 3)

  return(list(pca = pca, pca_df = pca_df, plot = p))
}

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
) {
  library(ggplot2)
  library(magrittr)

  x <- select_var_features(m, n_top_features = n_top_features)

  pca <- prcomp(t(m), center = center, scale. = scale.)
  pca_df <- data.frame(sample_info, pca$x)

  cme <- summary(pca)$importance["Cumulative Proportion", ]
  cme <- as.data.frame(cme)
  cme$PC <- rownames(cme)

  percent_var <- pca$sdev^2 / sum(pca$sdev^2)
  percent_var <- glue("PC{1:ncol(pca$x)} ({round(percent_var * 100)}%)")[1:n_components]

  aes_params <- list()

  if (!is.null(color_by)) {
    aes_params[["color"]] <- color_by
    color_legend <- c(2, 1)
    if (is.null(color_legend_lab)) color_legend_lab <- color_by
  } else {
    color_legend <- NULL
  }

  if (!is.null(shape_by)) {
    aes_params[["shape"]] <- shape_by
    shape_legend <- c(3, 1)
    if (is.null(shape_legend_lab)) shape_legend_lab <- shape_by
  } else {
    shape_legend <- NULL
  }

  aes_points <- do.call(aes_string, aes_params)

  if (n_components == 2) {
    p <- ggplot(
      pca_df,
      aes(x = PC1, y = PC2)
    ) +
      geom_point(size = point_size) +
      xlab(percent_var[1]) +
      ylab(percent_var[2]) +
      theme_bw() +
      aes_points

    if (!is.null(label_by)) {
      p <- p +
        ggrepel::geom_text_repel(aes_string(label = label_by), color = "black", size = text_size)
    }
  } else if (n_components > 2) {
    pc_columns <- str_c("PC", 1:n_components)

    p <- GGally::ggpairs(
      data = pca_df,
      mapping = aes_points,
      columns = pc_columns,
      legend = color_legend,
      upper = list(continuous = "points"),
      diag = list(continuous = "blankDiag"),
      columnLabels = percent_var,
      size = point_size,
      progress = FALSE
    ) +
      theme_bw()
  } else {
    stop(glue("Bad number of principal components (n_components = {n_components}). Must be >= 2."))
  }

  p <- p +
    labs(
      title = title,
      subtitle = subtitle,
      color = color_legend_lab,
      shape = shape_legend_lab
    )

  if (!is.null(output_file)) {
    cowplot::ggsave2(output_file, p, width = 6 + n_components, height = 4 + n_components)
  }

  return(list(pca = pca, pca_df = pca_df, plot = p))
}


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
) {
  library(RColorBrewer)
  library(gplots)
  library(magrittr)

  pal <- colorRampPalette(c("white", "darkblue"))(256)

  if (is.null(color_by)) {
    heatmap.2(m, main = main, col = pal, mar = c(10, 10),
              scale = "none", trace = "none")

  } else {
    color_by <- factor(color_by)
    cpal  <- levels(color_by) %>% length() %>% brewer.pal("Set1")
    csc   <- cpal[color_by]

    heatmap.2(m, ColSideColors = csc, main = main, col = pal, mar = c(10, 10),
              scale = "none", trace = "none")
    legend("topright", title = color_legend_lab, legend = levels(color_by), pch = 19,
           bty = "n", col = cpal)
  }

  return(invisible(TRUE))
}

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
) {
  if (!is.null(sample_data) && !is.null(column_color_by)) {
    sample_data <- sample_data[, column_color_by, drop = FALSE]
  }

  if (!is.null(feature_data) && !is.null(row_color_by)) {
    feature_data <- feature_data[, row_color_by, drop = FALSE]
  }

  pheatmap::pheatmap(m, annotation_col = sample_data, annotation_row = feature_data, main = main)
}

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
) {
  # heatmaply cannot handle NULL parameters.
  # Could be this implemented better?
  params <- list(m, key.title = key.title)

  if (!is.null(sample_data) && !is.null(column_color_by)) {
    params[["col_side_colors"]] <- sample_data[, column_color_by, drop = FALSE]
  }

  if (!is.null(feature_data) && !is.null(row_color_by)) {
    params[["row_side_colors"]] <- feature_data[, row_color_by, drop = FALSE]
  }

  if (!is.null(main)) {
    params[["main"]] <- main
  }

  do.call(heatmaply::heatmaply, params)
}

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
) {
  library(ggplot2)
  library(glue)
  library(friendlyeval)

  if (is.null(color_by)) {
    color <- "black"
  } else {
    color <- color_by
  }

  if (is.null(main)) {
    main <- feature
  }

  subtitle <- ""

  if (is.null(feature))
  {
    facet_by <- feature_col

    if (do_t_test) {
      subtitle <- "p-val: t-test between groups"
    }
  } else {
    subtitle <- feature
    plot_data <- plot_data %>%
      dplyr::filter(!!treat_string_as_col(feature_col) == !!feature)

    if (do_t_test) {
      subtitle <- glue("{subtitle}\np-val: t-test between groups")
    }
  }

  p <- ggpubr::ggboxplot(
    data = plot_data,
    x = x,
    y = y,
    facet.by = facet_by,
    color = color,
    xlab = x_lab,
    ylab = y_lab,
    title = main,
    subtitle = subtitle,
    outlier.shape = outlier_shape,
    add = add,
    repel = TRUE,
    add.params = list(size = point_size),
    ggtheme = theme_bw()
  ) +
    theme(legend.position = "top")

  if (do_t_test) {
    p <- p + ggpubr::stat_compare_means(aes(group = !!treat_string_as_col(x)), inherit.aes = TRUE, label = "p.format", method = "t.test", paired = FALSE)
  }

  return(p)
}

#' Compute the M value of CP values.
#'
#' Args:
#'  gene: Name of gene to compute the M value.
#'  cp_m: Dataframe of CP values. Rows are genes and columns are samples.
#'
#' Returns:
#'  M value (numeric)
compute_m <- function(gene, cp_m) {
  # Save the names of other reference genes.
  other_genes <- setdiff(rownames(cp_m), gene)

  # For each of the other reference genes, we apply the following function to its corresponding row in CP matrix.
  # MARGIN = 1 applies FUN to each row of input matrix (dataframe). MARGIN = 2 does so for columns instead of rows.
  V <- apply(cp_m[other_genes, ], MARGIN = 1, FUN = function(x) {
    # For each sample, we calculate expression ratio of the input reference gene and other reference gene.
    expression_ratios <- cp_m[gene, ] - x

    # We calculate standard deviation of these expression ratios.
    sd(expression_ratios)
  })

  # Final M value is mean of standard deviations of expression ratios.
  M <- mean(V)

  return(M)
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
  library(glue)
  library(magrittr)
  library(friendlyeval)

  test_name <- quote(t.test) %>% as.character()
  test_formula <- glue("{value_col} ~ {group_col}") %>% as.formula()

  sample_groups <- levels(gene_data[, group_col, drop = TRUE])

  if (length(sample_groups) != 2)
    stop("More than two groups to compare.")

  if (verbose)
    glue("{gene}: {sample_groups[2]} vs. {sample_groups[1]}") %>% cat(sep = "\n")

  res <- test(
    test_formula,
    data = dplyr::filter(gene_data, !!treat_string_as_col(gene_col) == !!gene)
  )

  p_value <- res$p.value

  if (verbose)
    glue("{test_name} p-value: {signif(p_value, 3)} {asterisk(p_value)}") %>% cat(sep = "\n")

  return(res)
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
  library(magrittr)

  genes <- unique(gene_data[, gene_col, drop = TRUE])

  test_results <- lapply(genes, function(gene) {
    test_gene(gene, gene_data, gene_col, value_col, group_col, test = test, verbose = FALSE)
  })

  p_values <- lapply(test_results, "[[", "p.value") %>% unlist()
  significance <- lapply(p_values, asterisk) %>% unlist()

  return(tibble::tibble(
    gene = genes,
    p_value = p_values,
    significance = significance,
    test_result = test_results
  ))
}

#' Return asterisks according to p-values.
#'
#' Args:
#'   p_value: Vector of p-values.
#'
#' Returns:
#'  character vector
asterisk <- function(p_value) {
  dplyr::case_when(
    is.na(p_value) ~ NA_character_,
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "NS"
  )
}

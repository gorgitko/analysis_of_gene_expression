#' Calculate hiearchical clustering and plot dendrogram.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   method_distance: Distance metric used for clustering. See ?dist.
#'                    Can be also correlation metrix ("pearson", "kendall", "spearman"). See ?cor.
#'   method_clustering: Clustering method. See ?hclust.
#'   color_by: Vector of discrete values to color samples by.
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
#SC
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

    assertthat::assert_that(length(levels(color_by)) <= 9, msg = glue("The '{group_name}' factor can have at most 9 distinct levels."))

    pal <- levels(color_by) %>% length() %>% brewer.pal("Set1")
    col <- pal[color_by]
    names(col) <- names(color_by)
    set_leaf_color <- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        attr(n, "nodePar") <- c(
          a$nodePar,
          list(lab.cex = 0.8, lab.col = col[a$label], pch = NA)
        )
      }
      return(n)
    }

    hc <- dendrapply(hc, set_leaf_color)
    plot(hc, main = title, xlab = xlab)
    legend(
      "topright", title = color_by_lab, legend = levels(color_by),
      pch = 19, col = pal, bty = "n"
    )
  } else {
    plot(hc, main = title, xlab = xlab)
  }

  return(invisible(hc))
#EC
}

#' Calculate hiearchical clustering and plot dendrogram using the dendextend package.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   method_distance: Distance metric used for clustering. See ?dist.
#'                    Can be also correlation metrix ("pearson", "kendall", "spearman"). See ?cor.
#'   method_clustering: Clustering method. See ?hclust.
#'   color_by: Vector of discrete values to color samples by.
#'             Length must be the same as number of columns in 'm'.
#'   color_by_lab: Name of the color legend.
#'   title: Main plot title.
#'
#' Returns:
#'   dendrogram object
plot_hc2 <- function(
  m,
  method_distance = "euclidean",
  method_clustering = "complete",
  color_by = NULL,
  color_by_lab = "Group",
  title = "Hierarchical Clustering"
) {
#SC
  library(glue)
  library(magrittr)
  library(dendextend)

  if (method_distance %in% c("pearson", "kendall", "spearman")) {
    di <- as.dist(1 - cor(m, method = method_distance))
    xlab <- ""
  } else {
    di <- dist(t(m), method = method_distance)
    xlab <- ""
  }

  hc <- hclust(di, method = method_clustering)
  hc <- as.dendrogram(hc, hang = 0.1)

  if (!is.null(color_by)) {
    if (!is.factor(color_by)) {
      color_by <- factor(color_by)
    }

    color_by <- tibble::tibble(sample_name = colnames(m), group = color_by)
    group_colors <- tibble::tibble(group = levels(color_by$group) %>% factor()) %>%
      dplyr::mutate(color = scales::hue_pal()(dplyr::n()))
    color_by <- dplyr::left_join(color_by, group_colors, by = "group") %>%
      dplyr::arrange(match(sample_name, labels(hc)))

    hc <- color_labels(hc, col = color_by$color)

    plot(hc, main = title, xlab = xlab)
    legend(
      "topright", title = color_by_lab, legend = group_colors$group,
      pch = 19, col = group_colors$color, bty = "n"
    )
  } else {
    plot(hc, main = title, xlab = xlab)
  }

  return(invisible(hc))
#EC
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
#SC
  if (!is.infinite(n_top_features))
  {
    # Calculate the variance for each gene.
    rv <- matrixStats::rowVars(m)

    # Select the n top genes by variance.
    select <- order(rv, decreasing = TRUE)[min(n_top_features, length(rv)) %>% seq_len()]
    m <- m[select, ]
  }

  return(m)
#EC
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
#SC
  library(ggplot2)
  library(patchwork)

  plot_type <- rlang::arg_match(plot_type)

  x <- select_var_features(m, n_top_features = n_top_features)

  pca <- prcomp(t(m), center = center, scale. = scale.)
  percent_var <- pca$sdev^2 / sum(pca$sdev^2)
  percent_var <- glue("PC{1:ncol(pca$x)} ({round(percent_var * 100)}%)")
  names(percent_var) <- glue("PC{1:ncol(pca$x)}")

  pca_df <- data.frame(sample_data, pca$x)

  if (plot_type == "single") {
    p_final <- ggplot(pca_df, aes_string(x = "PC1", y = "PC2", color = color_by, shape = shape_by)) +
      geom_point(size = point_size) +
      labs(x = percent_var[1], y = percent_var[2]) +
      theme_bw()

    if (!is.null(label_by)) {
      p_final <- p_final + ggrepel::geom_text_repel(aes_string(label = label_by), color = "black", size = text_size)
    }
  } else {
    cme <- summary(pca)$importance["Cumulative Proportion", ]
    cme <- as.data.frame(cme) %>%
      dplyr::mutate(PC = 1:nrow(.))

    pc_combinations <- list(c("PC2", "PC1"), c("PC3", "PC1"), c("PC2", "PC3"))
    p_list <- lapply(pc_combinations, function(pc) {
      p <- ggplot(pca_df, aes_string(x = pc[1], y = pc[2], color = color_by, shape = shape_by)) +
        geom_point(size = point_size) +
        labs(x = percent_var[pc[1]], y = percent_var[pc[2]]) +
        theme_bw()

      if (!is.null(label_by)) {
        p <- p + ggrepel::geom_text_repel(aes_string(label = label_by), color = "black", size = text_size)
      }

      return(p)
    })

    p_variance <- ggplot(cme, aes(x = PC, y = cme)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      scale_x_continuous(breaks = 1:nrow(cme)) +
      labs(y = "Cumulative % of var. explained")

    p_final <- wrap_plots(c(p_list, list(p_variance)), nrow = 2, guides = "collect") & theme(legend.position = "bottom")

  }

  return(list(pca = pca, pca_df = pca_df, plot = p_final))
#EC
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
#SC
  library(ggplot2)
  library(magrittr)

  assertthat::assert_that(n_components >= 2, msg = "Invalid n_components: must be >= 2")

  if (n_components == 2) {
    return(
      plot_pca(
        m, sample_data, plot_type = "single", n_top_features = n_top_features,
        color_by = color_by, shape_by = shape_by, label_by = label_by, point_size = point_size,
        text_size = text_size, center = center, scale. = scale.
      )
    )
  }

  x <- select_var_features(m, n_top_features = n_top_features)

  pca <- prcomp(t(m), center = center, scale. = scale.)
  pca_df <- data.frame(sample_data, pca$x)

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

  p <- p +
    labs(
      title = title,
      subtitle = subtitle,
      color = color_legend_lab,
      shape = shape_legend_lab
    )

  return(list(pca = pca, pca_df = pca_df, plot = p))
#EC
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
#SC
  library(ComplexHeatmap)

  if (z_score) {
    m <- t(m) %>% scale() %>% t()
  }

  if (is.null(color_mapping)) {
    ##-- To make annotation colors consistent, we generate a proper number of colors from the color_palette.
    ##-- We determine what is the maximum number of levels (or unique values) in both sample_annotation and feature_annotation,
    ##-- and generate this number of colors from color_palette.
    max_n_levels_s <- 0
    max_n_levels_f <- 0

    if (!is.null(sample_annotation)) {
      sample_annotation <- dplyr::select(sample_annotation, where(is.character), where(is.factor))
      max_n_levels_s <- purrr::map_int(sample_annotation, ~ length(unique(.)))
    }

    if (!is.null(feature_annotation)) {
      feature_annotation <- dplyr::select(feature_annotation, where(is.character), where(is.factor))
      max_n_levels_f <- purrr::map_int(feature_annotation, ~ length(unique(.)))
    }

    max_n_levels <- max(max_n_levels_s, max_n_levels_f)

    if (max_n_levels > 0) {
      annot_colors <- color_palette(max_n_levels)
      color_mapping <- lapply(c(as.list(sample_annotation), as.list(feature_annotation)), function(variable) {
        lvls <- as.factor(variable) %>% levels()
        colors <- annot_colors[1:length(lvls)]
        names(colors) <- lvls
        return(colors)
      })
    }
  }

  if (!is.null(sample_annotation)) {
    top_annotation <- HeatmapAnnotation(df = sample_annotation, col = color_mapping)
  } else {
    top_annotation <- NULL
  }

  if (!is.null(feature_annotation)) {
    right_annotation <- rowAnnotation(df = feature_annotation, col = color_mapping)
  } else {
    right_annotation <- NULL
  }

  p <- Heatmap(
    m,
    name = legend_title,
    column_title = title,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    top_annotation = top_annotation,
    right_annotation = right_annotation
  )

  return(p)
#EC
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
#SC
  if (z_score) {
    m <- t(m) %>% scale() %>% t()
  }

  ##-- heatmaply cannot handle NULL parameters
  ##-- Could be this implemented better?
  params <- list(
    m, key.title = legend_title, showticklabels = showticklabels, main = main,
    col_side_colors = sample_annotation, row_side_colors = feature_annotation
  )

  params <- purrr::discard(params, is.null)

  do.call(heatmaply::heatmaply, params)
#EC
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
#SC
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
    if (do_t_test) {
      subtitle <- "p-val: t-test between groups"
    }
  } else {
    subtitle <- feature
    plot_data <- plot_data %>%
      dplyr::filter(!!treat_string_as_col(facet_by) == !!feature)

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
#EC
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
#SC
  ##-- Save the names of other reference genes.
  other_genes <- setdiff(rownames(cp), gene)

  ##-- For each of the other reference genes, we apply the following function to its corresponding row in CP matrix.
  V <- apply(cp[other_genes, ], MARGIN = 1, FUN = function(x) {
    ##-- For each sample, we calculate expression ratio of the input reference gene and other reference gene.
    ##-- Note that CP values are log2-distributed and log-ratio can be written as difference of log arguments.
    expression_ratios <- cp[gene, ] - x

    ##-- We calculate standard deviation of these expression ratios.
    sd(expression_ratios)
  })

  ##-- Final M value is mean of standard deviations of expression ratios.
  M <- mean(V)

  return(M)
#EC
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
#SC
  library(glue)
  library(magrittr)
  library(friendlyeval)

  test_name <- quote(t.test) %>% as.character()
  test_formula <- glue("{value_col} ~ {group_col}") %>% as.formula()

  sample_groups <- levels(gene_data[, group_col, drop = TRUE])

  assertthat::assert_that(length(sample_groups) == 2, msg = "More than two groups to compare.")

  if (verbose)
    message(glue("{gene}: {sample_groups[2]} vs. {sample_groups[1]}"))

  res <- test(
    test_formula,
    data = dplyr::filter(gene_data, !!treat_string_as_col(gene_col) == !!gene)
  )

  p_value <- res$p.value

  if (verbose)
    message(glue("{test_name} p-value: {signif(p_value, 3)} {asterisk(p_value)}"))

  return(res)
#EC
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
#SC
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
#EC
}

#' Return asterisks according to p-values.
#'
#' Args:
#'   p_value: Vector of p-values.
#'
#' Returns:
#'  character vector
asterisk <- function(p_value) {
#SC
  dplyr::case_when(
    is.na(p_value) ~ NA_character_,
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "NS"
  )
#EC
}

plot_hc <- function(
  m,
  method_distance = "euclidean",
  method_clustering = "complete",
  color_by = NULL,
  color_by_lab = "Group",
  main = "Hierarchical Clustering"
) {
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

select_var_features <- function(m, n_top_features) {
  #' Select features (rows) with the highest variance across the samples (columns).
  #'
  #' Args:
  #'   m: Expression matrix (rows are features, columns are samples).
  #'   n_top_features: Number of top features.
  #'
  #' Returns:
  #'   Object 'm' subset to 'n_top_features' with the highest variance across the samples.

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

plot_pca <- function(
  m,
  color_by,
  n_top_features = Inf,
  center = TRUE,
  scale. = TRUE,
  color_by_lab = "Group"
) {
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
) {
  #' Using ggplot2, plot first three PCs of samples in expression matrix.
  #'
  #' Args:
  #'   m: Expression matrix (rows are features, columns are samples).
  #'   pheno_data: Dataframe describing samples.
  #'   n_top_features: Number of top features with the highest variance across the samples.
  #'   color_by: Column name in pheno_data to use for point coloring.
  #'   shape_by: Column name in pheno_data to use for point shape.
  #'   label_by: Column name in pheno_data to use for point labels.
  #'   point_size: Point size (numeric).
  #'   text_size: Label text size (numeric).
  #'   center: Whether to center PCA. See ?prcomp.
  #'   scale.: Whether to scale PCA. See ?prcomp.
  #'
  #' Returns:
  #'   list(pca = object of class 'prcomp', pca_df = combined dataframe of pheno_data and PCs, plot = ggplot2 object)

  library(rlist)
  library(ggplot2)
  library(ggthemes)
  library(cowplot)
  library(ggrepel)

  x <- select_var_features(m, n_top_features = n_top_features)

  pca <- prcomp(t(m), center = center, scale. = scale.)
  pca_df <- data.frame(pheno_data, pca$x)

  cme <- summary(pca)$importance["Cumulative Proportion", ]
  cme <- as.data.frame(cme)
  cme$PC <- rownames(cme)

  pc_combinations <- list(c("PC2", "PC1"), c("PC3", "PC1"), c("PC2", "PC3"))
  p_list <- lapply(pc_combinations, function(pc) {
    p <- ggplot(pca_df, aes_string(x = pc[1], y = pc[2], color = color_by, shape = shape_by)) +
      geom_point(size = point_size) +
      theme_few()

    if (!is.null(label_by)) {
      p <- p + geom_text_repel(aes_string(label = label_by), color = "black", size = text_size)
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
    ggplot(cme, aes(x = PC, y = cme)) + geom_bar(stat = "identity") + theme_few(),
    legend
  )

  p <- plot_grid(plotlist = p_list, nrow = 3)

  return(list(pca = pca, pca_df = pca_df, plot = p))
}


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
  #' Using GGally::ggpairs(), plot grid of PCA plots (PC1 vs. PC2, PC1 vs. PC3, PC2 vs. PC3, etc).
  #' When n_components == 2, use normal ggplot2.
  #'
  #' Args:
  #'   m: Expression matrix (rows are features, columns are samples).
  #'   pheno_data: Dataframe describing samples.
  #'   output_file: File to save plot in.
  #'   n_components: Number of PCs to plot.
  #'   n_top_features: Number of top features with the highest variance across the samples.
  #'   color_by: Column name in pheno_data to use for point coloring.
  #'   color_legend_lab: Name of the color legend.
  #'   shape_by: Column name in pheno_data to use for point shape.
  #'   shape_legend_lab: Name of the shape legend.
  #'   label_by: Column name in pheno_data to use for point labels.
  #'   point_size: Point size (numeric).
  #'   text_size: Label text size (numeric).
  #'   title: Plot title.
  #'   subtitle: Plot subtitle.
  #'   center: Whether to center PCA. See ?prcomp.
  #'   scale.: Whether to scale PCA. See ?prcomp.

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


plot_heatmap <- function(x, group = NULL, group_name = "Group", main = "") {
  library(RColorBrewer)
  library(gplots)
  library(magrittr)

  pal <- colorRampPalette(c("white", "darkblue"))(256)

  if (is.null(group)) {
    heatmap.2(x, main = main, col = pal, mar = c(10, 10),
              scale = "none", trace = "none")

  } else {
    group <- factor(group)
    cpal  <- levels(group) %>% length() %>% brewer.pal("Set1")
    csc   <- cpal[group]

    heatmap.2(x, ColSideColors = csc, main = main, col = pal, mar = c(10, 10),
              scale = "none", trace = "none")
    legend("topright", title = group_name, legend = levels(group), pch = 19,
           bty = "n", col = cpal)
  }

  return(invisible(TRUE))
}


ComputeM <- function(g, CP) {
  og <- setdiff(rownames(CP), g)

  V <- apply(CP[og, ], 1, function(x) {sd(CP[g, ] - x)})
  M <- mean(V)
  return(M)
}


PlotBW <- function(plotData, gene, ticks = -10:10) {
  library(lattice)

  trellis.device(device = "pdf", file = paste(gene, ".pdf", sep = ""),
                 width = 4, height = 4)

  labels <- 2^ticks
  labels[ticks < 0] <- paste("1/", 2^(-ticks[ticks < 0]), sep = "")

  print(bwplot(as.formula(paste(gene, "~ Sample_Group")),
               groups = Confluence,
               data = plotData,
               main = gene,
               ylab = "Relative expression\n(normalised to controls)",
               xlab = "Treatment\nConfluence",
               horizontal = FALSE,
               scales = list(y = list(at = ticks, labels = labels)),
               panel = function(...) {
                 panel.abline(h = ticks, lwd = 0.25, lty = 1,
                              col = "#eeeeee")
                 panel.stripplot(jitter.data = TRUE, ...)
                 panel.bwplot(do.out = FALSE, ...)
               }, auto.key = TRUE))
  dev.off()
  return(invisible(TRUE))
}


TestGene <- function(gene, data, test = t.test) {
  library(glue)
  library(dplyr)

  print(glue("{gene}: Controls vs. Spirulina"))

  res <- test(as.formula(delta_cp_centered ~ Sample_Group),
              data = data %>% filter(Symbol == gene))$p.value

  print(glue("p-value: {signif(res, 3)} {Asterisk(res)}"))
}


TestGeneTable <- function(data, test = t.test) {
  library(dplyr)

  genes <- unique(data$Symbol)
  p_values <- lapply(genes, function(gene) {
    test(as.formula(delta_cp_centered ~ Sample_Group),
         data = data %>% filter(Symbol == gene))$p.value
  }) %>% unlist()

  return(data.frame(
    Symbol = genes,
    p_value = p_values,
    significance = lapply(p_values, Asterisk) %>% unlist()
  ))
}


Asterisk <- function(p_value) {
  if (is.na(p_value)) {
    return("")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}


PlotBW_ggplot <- function(plot_data, gene, cp_col, facet_by = NULL, point_size = 2, do_t_test = TRUE) {
  #' Plot boxplot of gene expression.
  #'
  #' Args:
  #'   plot_data: data.frame (long format)
  #'   gene: symbol of gene to plot
  #'   cp_col: name of column with CP values
  #'   facet_by: name of column by which to make subplots (panels)
  #'   point_size: size of points inside boxplots
  #'   do_t_test: do t-test and display p-value inside the plot
  #'
  #' Returns:
  #'   ggplot2 object

  library(ggplot2)
  library(ggpubr)
  library(ggthemes)

  if (do_t_test) {
    title = sprintf("%s (p-val: t-test between groups)", gene)
  } else {
    title = gene
  }

  plot_data <- plot_data %>%
    dplyr::filter(Symbol == gene)

  p <- ggboxplot(data = plot_data,
                 x = "Sample_Group", y = cp_col,
                 facet.by = facet_by,
                 color = "Sample_Group",
                 xlab = "Sample Group", ylab = "Relative expression\n(normalised to controls)",
                 title = title,
                 outlier.shape = NA, add = "jitter", repel = TRUE, add.params = list(size = point_size),
                 ggtheme = theme_bw()) +
    theme(legend.position = "top")

  if (do_t_test)
    p <- p + stat_compare_means(aes(group = Sample_Group), inherit.aes = TRUE, label = "p.format", method = "t.test", paired = FALSE)

  return(p)
}


PlotHeatmap_heatmaply <- function(x, col_side_colors = NULL, showticklabels = c(TRUE, FALSE)) {
  #' Create heatmap with heatmaply.
  #'
  #' Args:
  #'   x: matrix or dataframe (wide format)
  #'   col_side_colors: numeric vector of factors to use for column colors; length(col_side_colors) == ncol(x)
  #'
  #' Returns:
  #'   plotly object

  library(heatmaply)

  heatmaply(x,
            col_side_colors = col_side_colors,
            showticklabels = showticklabels)
}

#' @title TfuncViz
#' @description A package for T cell functionality visualization and analysis
#' @import Seurat
#' @import tidyverse
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @import ggrepel
#' @import ggforce
#' @import RColorBrewer
#' @import ComplexHeatmap
#' @import circlize
#' @import GGally
#' @importFrom data.table fread
#' @importFrom stats prcomp scale na.omit
#' @importFrom utils str
#' @importFrom grDevices pdf dev.off colorRamp2
#' @importFrom grid gpar unit
NULL

#' Clean gene list by removing genes not in Seurat object
#'
#' @param gene_list A vector of gene names
#' @param seurat_obj Seurat object
#' @return A vector of genes present in the Seurat object
#' @export
clean_genes <- function(gene_list, seurat_obj = NULL) {
  if(is.null(seurat_obj)) {
    return(gene_list)
  }
  return(intersect(gene_list, rownames(seurat_obj)))
}

#' Calculate T cell functionality scores
#'
#' @param seurat_obj A Seurat object
#' @param gene_sets A list of gene sets for scoring
#' @param cluster_col Column name for clusters in metadata
#' @param group_col Column name for group in metadata
#' @param ctrl Number of control genes for AddModuleScore
#' @return A Seurat object with scores and a list of analysis results
#' @export
calculate_tfunc_scores <- function(seurat_obj,
                                   gene_sets = NULL,
                                   cluster_col = "seurat_clusters",
                                   group_col = "group",
                                   ctrl = 100) {

  # Use built-in gene sets if none provided
  if(is.null(gene_sets)) {
    # load T list
    gene_sets <- list(
      Naive = c("CCR7","SELL","IL7R","TCF7","LEF1","LTB","FOXO1","GZMK"),
      Activation = c("SELL", "CD40LG", "ANXA1", "IL2RA", "CD69"),
      Memory = c("CCR7", "TCF7", "CD69", "NR4A1", "MYADM", "GATA3", "TBX21"),
      Effector = c("BATF3", "IL2RA", "IFNG", "IRF4", "MYC", "SLC7A5", "SLC7A1", "XCL1", "GZMB", "CCL3", "CCL4", "IL2", "PRF1", "NKG7", "GNLY"),
      Tolerant = c("EGR2", "IKZF2", "IZUMO1R", "CD200", "DGKZ", "BTLA", "TOX", "CTLA4"),
      Exhausted = c("PDCD1", "LAG3", "HAVCR2", "CD244", "CD160", "IL10R", "EOMES", "NR4A2", "PTGER4", "TOX", "TOX2", "TIGIT", "CTLA4", "ENTPD1"),
      Proliferation = c("MKI67", "TK1", "STMN1"),
      Cytotoxicity = c("GZMK", "GZMA", "GZMB", "NKG7", "PRF1", "IFNG", "GNLY"),
      Regulatory = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "TNFRSF18", "IL10", "IKZF2", "CCR8")
    )
  }

  # Clean gene sets
  gene_sets <- lapply(gene_sets, clean_genes, seurat_obj = seurat_obj)

  # Calculate module scores
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = gene_sets,
    name = names(gene_sets),
    ctrl = ctrl
  )

  # Extract metadata and scores
  metadata <- seurat_obj@meta.data
  cluster_group <- metadata %>%
    dplyr::select(!!sym(cluster_col), !!sym(group_col))

  # 获取评分列名 - Seurat的AddModuleScore函数会在每个特征集名称后加1
  score_cols <- paste0(names(gene_sets), 1)

  # 检查评分列是否存在
  missing_cols <- setdiff(score_cols, colnames(metadata))
  if(length(missing_cols) > 0) {
    # 尝试查找实际的列名
    available_cols <- grep(paste(names(gene_sets), collapse="|"), colnames(metadata), value=TRUE)
    if(length(available_cols) == 0) {
      stop("无法找到由AddModuleScore函数生成的评分列。请检查Seurat对象。")
    }
    score_cols <- available_cols
  }

  # 提取评分
  scores <- metadata[, score_cols, drop = FALSE]

  # 重命名列,去掉末尾的数字
  names(scores) <- names(gene_sets)

  # 按集群和分组聚合
  score_matrix <- cluster_group %>%
    bind_cols(scores) %>%
    group_by(!!sym(cluster_col), !!sym(group_col)) %>%
    summarise(across(everything(), mean), .groups = "drop")

  # 标准化评分
  scaled_scores <- score_matrix %>%
    dplyr::select(-!!sym(cluster_col), -!!sym(group_col)) %>%
    scale() %>%
    as.data.frame()

  rownames(scaled_scores) <- paste0("C", score_matrix[[cluster_col]], "_", score_matrix[[group_col]])
  colnames(scaled_scores) <- names(gene_sets)

  # Run PCA
  pca_res <- prcomp(scaled_scores, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  pca_df$cluster_group <- rownames(pca_df)
  pca_df <- pca_df %>%
    mutate(
      cluster = paste0("C", str_extract(cluster_group, "\\d+")),
      group = str_extract(cluster_group, paste(unique(score_matrix[[group_col]]), collapse="|"))
    )

  # Return analysis results
  return(list(
    seurat_obj = seurat_obj,
    score_matrix = score_matrix,
    scaled_scores = scaled_scores,
    pca_res = pca_res,
    pca_df = pca_df,
    gene_sets = gene_sets,
    cluster_col = cluster_col,
    group_col = group_col
  ))
}
#' Create bubble plot with connection lines
#'
#' @param tfunc_results Results from calculate_tfunc_scores
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @return A ggplot object
#' @export
plot_tfunc_bubble <- function(tfunc_results,
                              title = "Spatial Distribution of T Cell Functional States",
                              subtitle = "PCA Reveals Cluster Functional Features and Disease-Health Differences") {

  pca_df <- tfunc_results$pca_df

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster, shape = group)) +
    # Add connection lines
    geom_segment(
      aes(xend = c(tail(PC1, -1), NA),
          yend = c(tail(PC2, -1), NA),
          group = rep(seq_along(unique(pca_df$cluster)), each = 2)),
      color = "grey70", linetype = "dashed", linewidth = 0.6,
      data = pca_df %>% group_by(cluster) %>% filter(n() == 2)
    ) +
    # Add points
    geom_point(aes(size = ifelse(group == unique(pca_df$group)[1], 1.2, 1)), alpha = 0.9) +
    # Add cluster labels
    geom_label_repel(
      aes(label = cluster),
      fill = alpha("white", 0.8),
      size = 4.5,
      fontface = "bold",
      box.padding = 0.3,
      label.padding = 0.2,
      show.legend = FALSE
    ) +
    # Add group labels
    geom_text_repel(
      aes(label = substr(group, 1, 1)),
      size = 3.5,
      color = "black",
      fontface = "bold",
      box.padding = 0.2,
      point.padding = 0.3
    ) +
    # Styling
    scale_color_viridis_d() +
    scale_shape_manual(values = c(16, 17)) +
    scale_size_identity() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "PC1 (main)",
      y = "PC2",
      color = "Cluster",
      shape = "Group",
      caption = paste(unique(pca_df$group)[1], "|", unique(pca_df$group)[2])
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray40", margin = margin(b = 20)),
      plot.caption = element_text(hjust = 0.5, size = 12, color = "gray50", margin = margin(t = 15)),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "gray30")
    ) +
    # Add ellipses
    geom_mark_ellipse(
      aes(group = cluster, fill = cluster),
      alpha = 0.08,
      expand = unit(1, "mm"),
      show.legend = FALSE
    ) +
    scale_fill_viridis_d()

  return(p)
}

#' Create vector plot showing direction of change
#'
#' @param tfunc_results Results from calculate_tfunc_scores
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @return A ggplot object
#' @export
plot_tfunc_vector <- function(tfunc_results,
                              title = "Disease-Health Group Functional State Change Vectors",
                              subtitle = "Arrows indicate the direction of functional state changes in the disease group relative to the healthy group") {

  pca_df <- tfunc_results$pca_df

  # Calculate differences between groups
  pca_diff <- pca_df %>%
    group_by(cluster) %>%
    filter(n() == 2) %>%
    summarise(
      delta_PC1 = PC1[group == unique(pca_df$group)[1]] - PC1[group == unique(pca_df$group)[2]],
      delta_PC2 = PC2[group == unique(pca_df$group)[1]] - PC2[group == unique(pca_df$group)[2]],
      centroid_PC1 = mean(PC1),
      centroid_PC2 = mean(PC2)
    )

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
    # Add change vectors
    geom_segment(
      aes(x = centroid_PC1, y = centroid_PC2,
          xend = centroid_PC1 + delta_PC1/3,
          yend = centroid_PC2 + delta_PC2/3),
      data = pca_diff,
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      size = 1.2,
      color = "black"
    ) +
    # Add points
    geom_point(aes(shape = group), size = 4) +
    # Add cluster labels
    geom_label_repel(
      aes(label = cluster),
      fill = alpha("white", 0.8),
      size = 4.5,
      fontface = "bold",
      box.padding = 0.5
    ) +
    # Add vector labels
    geom_text(
      aes(x = centroid_PC1 + delta_PC1/2.5,
          y = centroid_PC2 + delta_PC2/2.5,
          label = sprintf("Δ%.2f", sqrt(delta_PC1^2 + delta_PC2^2))),
      data = pca_diff,
      size = 3.5,
      color = "black",
      fontface = "bold",
      hjust = 0.5,
      vjust = 0.5
    ) +
    # Styling
    scale_color_viridis_d() +
    scale_shape_manual(values = c(16, 17)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "PC1",
      y = "PC2",
      color = "Cluster",
      shape = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    # Add ellipses
    geom_mark_ellipse(
      aes(group = cluster, fill = cluster),
      alpha = 0.08,
      expand = unit(1, "mm"),
      show.legend = FALSE
    ) +
    scale_fill_viridis_d()

  return(p)
}

#' Create faceted group comparison plot
#'
#' @param tfunc_results Results from calculate_tfunc_scores
#' @param title Plot title
#' @return A ggplot object
#' @export
plot_tfunc_facet <- function(tfunc_results,
                             title = "Comparison of Functional State Distribution Between Groups") {

  pca_df <- tfunc_results$pca_df

  # Create group labels for facets
  group_labels <- setNames(
    paste(c("Disease Group", "Health Group"),
          "(",
          unique(pca_df$group),
          ")",
          sep = " "),
    unique(pca_df$group)
  )

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(aes(shape = group), size = 5) +
    geom_text_repel(
      aes(label = paste0(cluster, " (", substr(group, 1, 1), ")")),
      size = 4,
      fontface = "bold",
      box.padding = 0.5
    ) +
    geom_mark_ellipse(
      aes(group = cluster, fill = cluster),
      alpha = 0.08,
      show.legend = FALSE
    ) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    scale_shape_manual(values = c(16, 17)) +
    facet_wrap(~ group, labeller = labeller(group = group_labels)) +
    labs(
      title = title,
      x = "PC1",
      y = "PC2",
      color = "Cluster",
      shape = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey90", color = NA),
      panel.spacing = unit(1.5, "lines")
    )

  return(p)
}

#' Create a heatmap of T cell functionality scores
#'
#' @param tfunc_results Results from calculate_tfunc_scores
#' @param filename Optional filename to save PDF
#' @param width PDF width
#' @param height PDF height
#' @return A ComplexHeatmap object
#' @export
plot_tfunc_heatmap <- function(tfunc_results, filename = NULL, width = 12, height = 8) {

  heatmap_data <- as.matrix(tfunc_results$scaled_scores)
  score_matrix <- tfunc_results$score_matrix

  # Create row annotations
  row_ha <- rowAnnotation(
    Group = score_matrix[[tfunc_results$group_col]],
    Cluster = paste0("C", score_matrix[[tfunc_results$cluster_col]]),
    col = list(
      Group = setNames(brewer.pal(length(unique(score_matrix[[tfunc_results$group_col]])), "Set1"),
                       unique(score_matrix[[tfunc_results$group_col]])),
      Cluster = setNames(brewer.pal(min(9, length(unique(score_matrix[[tfunc_results$cluster_col]]))), "Set1"),
                         paste0("C", unique(score_matrix[[tfunc_results$cluster_col]])))
    )
  )

  # Create heatmap
  hm <- Heatmap(
    heatmap_data,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("#1e90ff", "#ffffff", "#ff4500")),
    right_annotation = row_ha,
    row_split = factor(paste0("C", score_matrix[[tfunc_results$cluster_col]]),
                       levels = paste0("C", sort(unique(score_matrix[[tfunc_results$cluster_col]])))),
    column_split = rep(1, ncol(heatmap_data)),
    row_title = "Cluster_Group",
    column_title = "T cell function score",
    show_row_names = TRUE,
    row_names_side = "left",
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 9),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  )

  # Save to file if filename provided
  if(!is.null(filename)) {
    pdf(filename, width = width, height = height)
    draw(hm)
    dev.off()
  }

  return(hm)
}

#' Create a polar bubble plot of functionality differences
#'
#' @param tfunc_results Results from calculate_tfunc_scores
#' @return A ggplot object
#' @export
plot_tfunc_bubble_diff <- function(tfunc_results) {

  pca_df <- tfunc_results$pca_df

  # Calculate differences
  diff_data <- pca_df %>%
    group_by(cluster) %>%
    filter(n() == 2) %>%
    summarise(
      delta_PC1 = PC1[group == unique(pca_df$group)[1]] - PC1[group == unique(pca_df$group)[2]],
      delta_PC2 = PC2[group == unique(pca_df$group)[1]] - PC2[group == unique(pca_df$group)[2]],
      magnitude = sqrt(delta_PC1^2 + delta_PC2^2),
      direction = atan2(delta_PC2, delta_PC1)
    )

  # Create bubble plot
  p <- ggplot(diff_data, aes(x = cluster, y = magnitude, color = direction)) +
    geom_point(aes(size = magnitude), alpha = 0.8) +
    geom_segment(
      aes(xend = cluster, yend = 0),
      size = 1,
      alpha = 0.6
    ) +
    scale_color_gradient2(
      low = "#1f77b4",
      mid = "#ffffff",
      high = "#ff7f0e",
      midpoint = 0,
      limits = c(-pi, pi),
      breaks = c(-pi, 0, pi),
      labels = c(paste0("← ", unique(pca_df$group)[2], " Advantage"),
                 "No Difference",
                 paste0(unique(pca_df$group)[1], " Advantage →"))
    ) +
    scale_size(range = c(5, 15)) +
    coord_polar(theta = "x", start = -pi/2, direction = -1) +
    labs(
      title = "T-cell Functional States Intergroup Differences",
      subtitle = "Bubble Size Indicates Difference Magnitude, Color Indicates Direction",
      x = "",
      y = "Difference Magnitude",
      color = "Difference Direction",
      size = "Difference Magnitude"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.major = element_line(color = "grey90"),
      axis.text.x = element_text(size = 12, face = "bold"),
      legend.position = "right"
    )

  return(p)
}

#' Create a parallel coordinates plot of T cell functionality
#'
#' @param tfunc_results Results from calculate_tfunc_scores
#' @return A ggplot object
#' @export
plot_tfunc_parallel <- function(tfunc_results) {

  # Prepare data
  parcoord_data <- tfunc_results$scaled_scores %>%
    as.data.frame() %>%
    rownames_to_column("cluster_group") %>%
    mutate(
      cluster = str_extract(cluster_group, "C\\d+"),
      group = str_extract(cluster_group, paste(unique(tfunc_results$score_matrix[[tfunc_results$group_col]]), collapse="|"))
    ) %>%
    dplyr::select(-cluster_group)

  # Create parallel coordinates plot
  p <- ggparcoord(
    parcoord_data,
    columns = 1:ncol(tfunc_results$scaled_scores),
    groupColumn = "cluster",
    scale = "globalminmax",
    alphaLines = 0.6,
    showPoints = TRUE,
    title = "Parallel Coordinates Plot of T Cell Functional States"
  ) +
    facet_grid(~ group) +
    scale_color_brewer(palette = "Set1") +
    labs(
      x = "Functional State",
      y = "Normalized Score",
      color = "Cluster"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  return(p)
}

#' Run the complete T cell functionality analysis pipeline
#'
#' @param seurat_obj A Seurat object
#' @param gene_sets A list of gene sets for scoring (optional)
#' @param cluster_col Column name for clusters in metadata
#' @param group_col Column name for group in metadata
#' @param output_dir Directory to save plots (optional)
#' @param prefix Prefix for output files
#' @return A list with all results and plots
#' @export
analyze_tcell_functionality <- function(seurat_obj,
                                        gene_sets = NULL,
                                        cluster_col = "seurat_clusters",
                                        group_col = "group",
                                        output_dir = NULL,
                                        prefix = "tfunc") {

  # Calculate scores
  results <- calculate_tfunc_scores(seurat_obj, gene_sets, cluster_col, group_col)

  # Create plots
  plots <- list(
    bubble = plot_tfunc_bubble(results),
    vector = plot_tfunc_vector(results),
    facet = plot_tfunc_facet(results),
    bubble_diff = plot_tfunc_bubble_diff(results),
    parallel = plot_tfunc_parallel(results)
  )

  # Save plots if output_dir is provided
  if(!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # Save individual plots
    for(name in names(plots)) {
      filename <- file.path(output_dir, paste0(prefix, "_", name, ".pdf"))
      pdf(filename, width = 10, height = 8)
      print(plots[[name]])
      dev.off()
    }

    # Save heatmap separately
    heatmap_file <- file.path(output_dir, paste0(prefix, "_heatmap.pdf"))
    plot_tfunc_heatmap(results, filename = heatmap_file)
  }

  # Return all results
  return(c(
    results,
    list(plots = plots)
  ))
}

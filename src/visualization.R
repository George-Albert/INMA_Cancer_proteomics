build_pca_plot <- function(combined_data, variance_pc1, variance_pc2) {
  fill_base <- combined_data$Fill[!duplicated(combined_data$short_setup)]
  color_base <- combined_data$Color[!duplicated(combined_data$short_setup)]
  labels <- factor(unique(combined_data$short_setup), levels = unique(combined_data$short_setup))

  combined_data$short_setup <- factor(combined_data$short_setup, levels = unique(combined_data$short_setup))

  ggplot2::ggplot(combined_data, ggplot2::aes(x = PC1, y = PC2, fill = short_setup, color = short_setup)) +
    ggplot2::geom_point(size = 3, shape = 21, stroke = 1.5) +
    ggrepel::geom_label_repel(
      ggplot2::aes(label = short_setup),
      color = "black",
      fill = "white",
      max.overlaps = 21,
      size = 2,
      nudge_x = 0.1,
      nudge_y = 0.1,
      show.legend = FALSE
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = fill_base, name = "Samples", labels = labels) +
    ggplot2::scale_color_manual(values = color_base, name = "Samples", labels = labels) +
    ggplot2::xlab(paste0("PC1:", round(100 * variance_pc1, 2), "% variance explained")) +
    ggplot2::ylab(paste0("PC2:", round(100 * variance_pc2, 2), "% variance explained")) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.line = ggplot2::element_line(colour = "black"),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1, linetype = "solid"),
      legend.text = ggplot2::element_text(size = 14),
      legend.key.size = grid::unit(1, "lines")
    )
}

normalize_expression <- function(reads, rm_ext = "both") {
  exp <- log2(reads + 1)
  preprocessCore::normalize.quantiles.robust(
    as.matrix(exp),
    copy = FALSE,
    remove.extreme = rm_ext,
    n.remove = 1,
    use.median = FALSE,
    use.log2 = FALSE
  )
}

build_design <- function(metadata) {
  metadata$short_setup <- factor(metadata$short_setup, levels = unique(metadata$short_setup))
  metadata$Time <- factor(metadata$Time, levels = c("4h", "24h"))
  metadata$Particle <- factor(metadata$Particle, levels = c("Au", "PEG"))
  group <- factor(metadata$short_setup)
  design <- stats::model.matrix(~0 + group, data = metadata$short_setup)
  colnames(design) <- levels(group)
  design
}

build_default_contrasts <- function(design) {
  au_24_vs_au_4 <- limma::makeContrasts(Au_24h - Au_4h, levels = design)
  peg_24_vs_peg_4 <- limma::makeContrasts(PEG_24h - PEG_4h, levels = design)
  peg_4_vs_au_4 <- limma::makeContrasts(PEG_4h - Au_4h, levels = design)
  peg_24_vs_au_24 <- limma::makeContrasts(PEG_24h - Au_24h, levels = design)
  peg_vs_au <- limma::makeContrasts((PEG_24h - PEG_4h) - (Au_24h - Au_4h), levels = design)

  out <- cbind(au_24_vs_au_4, peg_24_vs_peg_4, peg_4_vs_au_4, peg_24_vs_au_24, peg_vs_au)
  colnames(out) <- c("Au_24hvsAu_4h", "PEG_24hvsPEG_4h", "PEG_4hvsAu_4h", "PEG_24hvsAu_24h", "PEGvsAu")
  out
}

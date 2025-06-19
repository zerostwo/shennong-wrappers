#' Calculate ROGUE Score(s) from a Seurat Object
#'
#' Computes entropy-based **ROGUE** scores to assess cluster purity from a Seurat object.
#' This function supports three usage modes:
#'
#' - **Full dataset ROGUE** calculation (no cluster provided)
#' - **Per-cluster ROGUE** (when both `cluster` and `sample` are provided)
#' - **Multi-resolution ROGUE comparison**, using multiple `cluster` columns
#'
#' @param x A Seurat object containing a `counts` layer.
#' @param cluster A character vector of one or more clustering columns in the metadata (`x@meta.data`).
#'   If multiple, calculates mean ROGUE per column and identifies the elbow/knee point.
#' @param sample *(Optional)* Column name indicating sample identity; required for per-cluster ROGUE.
#' @param span Smoothing parameter for [ROGUE::rogue()] when `sample` is used. Default: `0.9`.
#' @param min_cells Minimum number of cells a gene must be expressed in. Default: `10`.
#' @param min_genes Minimum number of genes a cell must express. Default: `10`.
#' @param platform Either `"UMI"` (droplet-based) or `"full-length"` (SMART-seq, etc.).
#' @param cutoff Adjusted p-value cutoff for [ROGUE::CalculateRogue()]. Default: `0.05`.
#' @param features *(Optional)* Character vector of gene names to use. If `NULL`, all filtered genes are used.
#' @param verbose Whether to show progress messages. Default: `TRUE`.
#'
#' @return A `data.frame` with ROGUE results:
#'
#' - For **multi-resolution comparison**:
#'   Returns columns: `cluster`, `average_rogue`, `is_knee_point`.
#'
#' - For **per-cluster analysis**:
#'   Returns per-cluster ROGUE scores.
#'
#' - For **global** (no cluster):
#'   Returns ROGUE score vector.
#'
#' @details
#' The **ROGUE** (Relative entropy Of Gene Usage Entropy) score quantifies the homogeneity of a cell population
#' using Shannon entropy. High ROGUE scores indicate purer (i.e., less heterogeneous) cell groups.
#'
#' This function wraps around:
#'
#' - [ROGUE::SE_fun()]
#' - [ROGUE::CalculateRogue()]
#' - [ROGUE::rogue()]
#'
#' @references
#' Liu, B., Li, C., Li, Z., Wang, D., Ren, X., & Zhang, Z. (2020).
#' An entropy-based metric for assessing the purity of single cell populations.
#' *Nature Communications*, **11**(1), 1–13.
#' [https://doi.org/10.1038/s41467-020-16904-3](https://doi.org/10.1038/s41467-020-16904-3)
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' seurat_obj <- CreateSeuratObject(counts = example_matrix)
#' seurat_obj <- FindClusters(seurat_obj, resolution = c(0.2, 0.4, 0.6))
#'
#' # Multi-resolution comparison
#' sn_calculate_rogue(seurat_obj, cluster = c("RNA_snn_res.0.2", "RNA_snn_res.0.4", "RNA_snn_res.0.6"))
#'
#' # Per-cluster ROGUE
#' sn_calculate_rogue(seurat_obj, cluster = "seurat_clusters", sample = "sample_id")
#'
#' # Global ROGUE
#' sn_calculate_rogue(seurat_obj)
#' }
#'
#' @export
sn_calculate_rogue <- function(x,
                               cluster = NULL,
                               sample = NULL,
                               span = 0.9,
                               min_cells = 10,
                               min_genes = 10,
                               platform = "UMI",
                               cutoff = 0.05,
                               features = NULL,
                               verbose = TRUE) {
  check_installed_github("ROGUE", "PaulingLiu/ROGUE")
  check_installed("SeuratObject", reason = "to extract counts from Seurat")

  if (!inherits(x, "Seurat")) {
    cli_abort("Input {.arg x} must be a Seurat object.")
  }

  if (!"counts" %in% SeuratObject::Layers(x)) {
    cli_abort("The {.code 'counts'} layer is not found in the Seurat object.")
  }
  platform <- arg_match(platform, c("UMI", "full-length"))

  metadata <- x[[]]
  counts <- SeuratObject::LayerData(x, layer = "counts") |>
    Matrix::as.matrix()

  if (verbose) cli_alert_info("Filtering matrix (min.cells = {min_cells}, min.genes = {min_genes})...")
  counts <- ROGUE::matr.filter(counts, min.cells = min_cells, min.genes = min_genes)

  if (verbose) cli_alert_info("Calculating entropy...")
  entropy <- ROGUE::SE_fun(counts)

  # Multi-cluster average ROGUE mode
  if (length(cluster) > 1) {
    cluster <- as_character(cluster)

    if (!all(cluster %in% colnames(metadata))) {
      missing <- setdiff(cluster, colnames(metadata))
      cli_abort("Column(s) {.val {missing}} not found in metadata.")
    }

    if (verbose) cli_alert_info("Calculating average ROGUE score for {length(cluster)} clustering resolutions...")

    rogue_df <- purrr::map_dfr(cluster, function(clust) {
      labels <- as.character(metadata[[clust]])
      scores <- ROGUE::CalculateRogue(
        entropy,
        platform = platform,
        cutoff = cutoff,
        features = features
      )
      tibble::tibble(cluster = clust, average_rogue = mean(scores, na.rm = TRUE))
    })

    # Elbow/knee point detection (largest jump in rogue score)
    rogue_df <- rogue_df |> dplyr::arrange(cluster)
    d_diff <- diff(rogue_df$average_rogue)
    knee_idx <- which.max(d_diff) + 1
    rogue_df$is_knee_point <- seq_len(nrow(rogue_df)) == knee_idx

    if (verbose) {
      knee_label <- rogue_df$cluster[knee_idx]
      cli_alert_success("Best resolution candidate: {.val {knee_label}} (detected by elbow method)")
    }

    return(rogue_df)
  }

  # Per-cluster ROGUE (requires sample + single cluster)
  if (!is.null(cluster) && !is.null(sample)) {
    cluster <- as_string(cluster)
    sample <- as_string(sample)

    if (!cluster %in% colnames(metadata)) {
      cli_abort("Column {.val {cluster}} not found in metadata.")
    }
    if (!sample %in% colnames(metadata)) {
      cli_abort("Column {.val {sample}} not found in metadata.")
    }

    if (verbose) cli_alert_info("Calculating ROGUE per cluster (by sample)...")

    rogue_result <- ROGUE::rogue(
      expr = counts,
      labels = as.character(metadata[[cluster]]),
      samples = as.character(metadata[[sample]]),
      platform = platform,
      span = span
    )

    rogue_result <- as.data.frame(rogue_result)
    rogue_result$cluster <- rownames(rogue_result)
    rownames(rogue_result) <- NULL

    if (verbose) cli_alert_success("ROGUE score calculation completed.")
    return(rogue_result)
  }

  # Fallback: full dataset rogue score
  if (verbose) cli_alert_info("Calculating ROGUE score for all cells (no cluster provided)...")

  rogue_result <- ROGUE::CalculateRogue(entropy, platform = platform)
  return(rogue_result)
}

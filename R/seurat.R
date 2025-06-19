#' Calculate ROGUE Score(s) from a Seurat Object
#'
#' Computes **ROGUE** (entropy-based) scores to evaluate cell population purity
#' from a Seurat object. Supports both global and per-cluster ROGUE computation.
#'
#' ## Modes:
#' - If `clusters` and `label` are both `NULL`, computes global ROGUE score across all cells.
#' - If both are provided, computes per-cluster ROGUE scores across samples.
#'
#' @param x A Seurat object with a `counts` layer.
#' @param clusters *(optional)* Character vector of metadata column names containing clustering info.
#'   If provided, `label` must also be provided.
#' @param label *(optional)* A metadata column name specifying sample identity. Required when `clusters` is used.
#' @param span Smoothing parameter used in `ROGUE::rogue()` when `clusters` and `label` are used. Default is `0.9`.
#' @param min_cells Minimum number of cells a gene must be detected in. Default: `10`.
#' @param min_genes Minimum number of genes each cell must express. Default: `10`.
#' @param platform Character string, either `"UMI"` or `"full-length"` to match sequencing platform. Default: `"UMI"`.
#' @param cutoff Adjusted p-value threshold for `CalculateRogue`. Default: `0.05`.
#' @param features Optional gene list to use for scoring.
#' @param verbose Logical, whether to print messages. Default: `TRUE`.
#'
#' @return A `data.frame` containing:
#' - `label`: Sample identity or "All"
#' - `cluster`: Cluster identity
#' - `score`: ROGUE score
#' - `group`: Source of the score (for faceting)
#'
#' @references
#' Liu, B., Li, C., Li, Z., Wang, D., Ren, X., & Zhang, Z. (2020).
#' An entropy-based metric for assessing the purity of single cell populations.
#' *Nature Communications*, **11**(1), 1–13.
#' [https://doi.org/10.1038/s41467-020-16904-3](https://doi.org/10.1038/s41467-020-16904-3)
#'
#' @examples
#' \dontrun{
#' # Global ROGUE score
#' sn_calculate_rogue(seurat_obj)
#'
#' # Per-cluster ROGUE score
#' sn_calculate_rogue(seurat_obj, clusters = "seurat_clusters", label = "sample")
#'
#' # Multiple clustering resolutions
#' sn_calculate_rogue(seurat_obj, clusters = c("res_0.4", "res_1.0"), label = "sample")
#' }
#'
#' @export
sn_calculate_rogue <- function(x,
                               clusters = NULL,
                               label = NULL,
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
    stop("Input `x` must be a Seurat object.")
  }

  if (!"counts" %in% SeuratObject::Layers(x)) {
    stop("The 'counts' layer is not found in the Seurat object.")
  }

  platform <- match.arg(platform, c("UMI", "full-length"))

  metadata <- x[[]]

  # Check metadata columns
  if (!is.null(clusters) && any(!clusters %in% colnames(metadata))) {
    stop("Some elements in `clusters` are not found in metadata columns: ",
         paste(setdiff(clusters, colnames(metadata)), collapse = ", "))
  }
  if (xor(is.null(clusters), is.null(label))) {
    stop("`clusters` and `label` must either both be provided or both be NULL.")
  }

  counts <- SeuratObject::LayerData(x, layer = "counts")
  counts <- Matrix::as.matrix(counts)


  counts <- ROGUE::matr.filter(counts, min.cells = min_cells, min.genes = min_genes)
  entropy <- ROGUE::SE_fun(counts)
  all_score <- ROGUE::CalculateRogue(entropy, platform = platform, cutoff = cutoff, features = features)

  result_list <- list(
    data.frame(
      label = "All",
      name = "All",
      value = all_score,
      group = "All",
      stringsAsFactors = FALSE
    )
  )

  if (!is.null(clusters)) {
    for (cluster in clusters) {
      if (verbose) message("Calculating ROGUE for cluster: ", cluster)
      rogue_result <- ROGUE::rogue(
        expr = counts,
        labels = as.character(metadata[[cluster]]),
        samples = as.character(metadata[[label]]),
        platform = platform,
        span = span
      )

      df <- stack(as.data.frame(rogue_result))
      df$label <- rownames(rogue_result)
      df$group <- cluster
      names(df)[1:2] <- c("value", "name")

      result_list[[length(result_list) + 1]] <- df
    }
  }

  result_df <- do.call(rbind, result_list)
  result_df$group <- factor(result_df$group, levels = c("All", clusters))
  colnames(result_df) <- c("label", "cluster", "score", "group")
  return(result_df)
}

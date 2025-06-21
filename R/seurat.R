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
  if (!is.null(clusters) &&
    any(!clusters %in% colnames(metadata))) {
    stop(
      "Some elements in `clusters` are not found in metadata columns: ",
      paste(setdiff(clusters, colnames(metadata)), collapse = ", ")
    )
  }
  if (xor(is.null(clusters), is.null(label))) {
    stop("`clusters` and `label` must either both be provided or both be NULL.")
  }

  counts <- SeuratObject::LayerData(x, layer = "counts")
  counts <- Matrix::as.matrix(counts)


  counts <- ROGUE::matr.filter(counts, min.cells = min_cells, min.genes = min_genes)
  entropy <- ROGUE::SE_fun(counts)
  all_score <- ROGUE::CalculateRogue(entropy,
    platform = platform,
    cutoff = cutoff,
    features = features
  )

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
      if (verbose) {
        message("Calculating ROGUE for cluster: ", cluster)
      }
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

#' Run PySCENIC on a Seurat Object
#'
#' @param x A Seurat object.
#' @param assay The assay to use for the analysis. Default: `"RNA"`.
#' @param layer The layer to use for the analysis. Default: `"counts"`.
#' @param outdir The directory to save the results. Default: `"./pyscenic"`.
#' @param species The species to use for the analysis. Default: `"human"`.
#' @param features The features to use for the analysis. Default: `NULL`.
#' @param tfs_fname The path to the TF list file. Default: `NULL`.
#' @param database_fname The path to the database file. Default: `NULL`.
#' @param annotations_fname The path to the annotations file. Default: `NULL`.
#' @param mode The mode to use for the analysis. Default: `"custom_multiprocessing"`.
#' @param ncores The number of cores to use for the analysis. Default: `8`.
#' @param seed The seed to use for the analysis. Default: `717`.
#' @param overwrite Whether to overwrite the results. Default: `FALSE`.
#'
#' @return A list of results. The results are saved in the `outdir` directory.
#'
#' @examples
#' \dontrun{
#' sn_run_pyscenic(seurat_obj, species = "human")
#' }
#'
#' @export
sn_run_pyscenic <- function(x,
                            assay = "RNA",
                            layer = "counts",
                            outdir = "./pyscenic",
                            species = "human",
                            features = NULL,
                            tfs_fname = NULL,
                            database_fname = NULL,
                            annotations_fname = NULL,
                            mode = "custom_multiprocessing",
                            ncores = 8,
                            seed = 717,
                            overwrite = FALSE) {
  check_installed_github(pkg = "ShennongTools", repo = "zerostwo/shennong-tools")

  if (species == "human") {
    tfs_fname <- tfs_fname %||% "/mnt/resources/cistarget/tf_lists/allTFs_hg38.txt"
    database_fname <- database_fname %||% c(
      "/mnt/resources/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
      "/mnt/resources/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
    )
    annotations_fname <- annotations_fname %||% "/mnt/resources/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
  } else if (species == "mouse") {
    tfs_fname <- tfs_fname %||% "/mnt/resources/cistarget/tf_lists/allTFs_mm.txt"
    database_fname <- database_fname %||% c(
      "/mnt/resources/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
      "/mnt/resources/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
    )
    annotations_fname <- annotations_fname %||% "/mnt/resources/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
  } else {
    stop("Please check if the species is human or mouse!")
  }

  expression_mtx_fname <- file.path(outdir, "expression_mtx.csv")
  module_fname <- file.path(outdir, "expression_mtx.adjacencies.tsv")
  signatures_fname <- file.path(outdir, "regulons.gmt")
  auc_mtx_fname <- file.path(outdir, "auc_mtx.csv")
  tfs_target_fname <- file.path(outdir, "tfs_targer.tsv")
  auc_g_mtx_fname <- file.path(outdir, "auc_g_mtx.csv")
  dir_create(outdir)
  if (!is_null(features)) {
    tfs <- read.table(tfs_fname)[, 1]
    filtered_tfs <- intersect(tfs, features)
    tfs_fname <- file.path(outdir, "tfs.txt")
    write.table(
      x = filtered_tfs,
      file = tfs_fname,
      row.names = FALSE,
      quote = FALSE,
      col.names = FALSE
    )
    counts <- as.data.frame(Matrix::as.matrix(
      SeuratObject::LayerData(
        object = x,
        assay = assay,
        layer = layer
      )[features, ]
    ))
  } else {
    counts <- as.data.frame(Matrix::as.matrix(
      SeuratObject::LayerData(
        object = x,
        assay = assay,
        layer = layer
      )
    ))
  }
  if (!file.exists(expression_mtx_fname)) {
    cli_inform("Write counts to a csv file...")
    data.table::fwrite(
      x = counts,
      file = expression_mtx_fname,
      quote = FALSE,
      row.names = TRUE
    )
  }

  grn <- ShennongTools::sn_run(
    tool_name = "pyscenic",
    command = "grn",
    expression = expression_mtx_fname,
    tf_list = tfs_fname,
    adjacencies = module_fname,
    transpose = TRUE,
    seed = seed,
    method = "grnboost2",
    threads = ncores,
    overwrite = overwrite
  )

  ctx <- ShennongTools::sn_run(
    tool_name = "pyscenic",
    command = "ctx",
    module = grn@outputs$adjacencies,
    databases = database_fname,
    annotations = annotations_fname,
    expression = grn@inputs$expression,
    transpose = TRUE,
    regulons = signatures_fname,
    ncores = ncores,
    overwrite = overwrite
  )

  aucell <- ShennongTools::sn_run(
    tool_name = "pyscenic",
    command = "aucell",
    regulons = ctx@outputs$regulons,
    expression = grn@inputs$expression,
    auc_matrix = auc_mtx_fname,
    transpose = TRUE,
    ncores = ncores,
    overwrite = overwrite
  )

  if (!file.exists(tfs_target_fname)) {
    signatures <- clusterProfiler::read.gmt(signatures_fname)
    signatures_count <- as.data.frame(table(signatures$term))
    colnames(signatures_count) <- c("term", "n")

    signatures <- merge(signatures, signatures_count, by = "term", all.x = TRUE)

    split_term <- strsplit(as.character(signatures$term), "\\(")
    simplified_term <- sapply(split_term, function(x) x[[1]])
    simplified_term <- trimws(simplified_term)
    signatures$symbol <- paste0(simplified_term, " (", signatures$n, "g)")
    signatures$tf <- simplified_term
    signatures$target_gene <- signatures$gene
    signatures <- signatures[, c("symbol", "tf", "target_gene")]
    write.csv(signatures, tfs_target_fname, row.names = FALSE)
  } else {
    signatures <- read.csv(tfs_target_fname)
  }
  invisible(list(
    grn = grn,
    ctx = ctx,
    aucell = aucell,
    signatures = signatures
  ))
}

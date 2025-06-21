#' Fetch Cistarget resources
#'
#' @param species The species to use for the analysis. Default: `"human"`.
#' @param outdir The directory to save the results. Default: `NULL`.
#' @param overwrite Whether to overwrite the results. Default: `FALSE`.
#'
#' @return A list of local paths to the downloaded files.
#'
#' @examples
#' \dontrun{
#' sn_fetch_cistarget_resources(species = "human")
#' }
#'
#' @export
sn_fetch_cistarget_resources <- function(species = c("human", "mouse"),
                                         outdir = NULL,
                                         overwrite = FALSE) {
  species <- match.arg(species)
  if (is.null(outdir)) {
    outdir <- tools::R_user_dir("shennong-wrappers", which = "data")
    outdir <- file.path(outdir, "cistarget", species)
  }

  dir_create(outdir)

  resources <- switch(species,
    human = list(
      tf_list = "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
      motif2tf = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
      databases = c(
        "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
      )
    ),
    mouse = list(
      tf_list = "https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt",
      motif2tf = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
      databases = c(
        "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
      )
    )
  )

  tf_file <- .download_resource(resources$tf_list, outdir, overwrite = overwrite)
  motif_file <- .download_resource(resources$motif2tf, outdir, overwrite = overwrite)
  db_files <- vapply(resources$databases, .download_resource, character(1), outdir = outdir, overwrite = overwrite)
  names(db_files) <- NULL

  invisible(list(
    tf_list = tf_file,
    motif2tf = motif_file,
    databases = db_files
  ))
}


.download_resource <- function(url, outdir, overwrite = FALSE) {
  dest_file <- file.path(outdir, basename(url))
  if (!file.exists(dest_file) || overwrite) {
    message(sprintf("Downloading %s...", basename(url)))
    utils::download.file(url, destfile = dest_file, mode = "wb", quiet = TRUE)
  } else {
    message(sprintf("Using cached file: %s", basename(url)))
  }
  return(dest_file)
}

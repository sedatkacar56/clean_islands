#' Remove Small Islands From a Metadata Label on an Embedding
#'
#' Cleans a categorical metadata column by removing small, isolated groups of
#' cells (islands) and stray cells that sit apart from the main body of each
#' label on a 2D embedding such as UMAP. Cells are grouped with DBSCAN, and only
#' groups with at least `min_cells` cells are kept. The result is written to a
#' new metadata column, so the original column is never changed.
#'
#' @param seurat_obj A Seurat object with the reduction named in `reduction`.
#' @param col Character. Metadata column to clean, for example a region, a
#'   cluster or a cell type column.
#' @param values Character vector or `NULL`. Labels of `col` to clean. Default
#'   `NULL`: every non-`NA` label is cleaned, each one separately. Labels not
#'   listed here are copied unchanged.
#' @param reduction Character. Reduction used for the cell positions.
#'   Default `"umap"`.
#' @param dims Integer vector. Dimensions of the reduction to use.
#'   Default `1:2`.
#' @param min_cells Integer. Smallest group of connected cells that is kept.
#'   Smaller groups are removed. Default `300`.
#' @param eps Numeric or `NULL`. DBSCAN neighborhood radius in embedding units:
#'   two cells closer than this count as neighbors. Default `0.2`, which suits a
#'   typical UMAP. Set `NULL` to use 1 percent of the average width of the embedding.
#'   This is a heuristic, so check the result visually.
#' @param min_pts Integer. Number of neighbors within `eps` that a cell needs to
#'   count as part of a dense area. Default `30`.
#' @param suffix Character. Added to `col` to make the new column name.
#'   Default `"_clean"`.
#' @param new_col Character or `NULL`. Full name of the new column. Overrides
#'   `suffix` when given.
#' @param verbose Logical. Print how many cells were kept and removed per
#'   label. Default `TRUE`.
#'
#' @return The Seurat object with one new metadata column. Cleaned cells keep
#'   their label, removed cells are set to `NA`, and cells that were `NA` stay
#'   `NA`.
#'
#' @details
#' A cell is removed if it is DBSCAN noise (fewer than `min_pts` cells within
#' `eps`) or if it belongs to a connected group smaller than `min_cells`. If
#' `eps` is too small, one real region breaks into pieces and some are lost. If
#' it is too large, separate regions merge together, so check the result with
#' `DimPlot(obj, group.by = "<new column>")`.
#'
#' The function only changes labels. It does not test anything, and it works
#' with the positions on the embedding, not with expression. Requires the
#' dbscan and Seurat packages.
#'
#' @examples
#' \dontrun{
#' # clean a region label, new column "region_clean"
#' obj <- clean_islands(obj, col = "region")
#'
#' # larger minimum size and a chosen radius
#' obj <- clean_islands(obj, "region", min_cells = 500, eps = 0.25)
#'
#' # only clean one label of a cell type column, custom column name
#' obj <- clean_islands(obj, "celltype", values = "CD8_T", new_col = "celltype_pure")
#'
#' Seurat::DimPlot(obj, group.by = "region_clean")
#' }
#'
#' @export
clean_islands <- function(seurat_obj, col, values = NULL,
                          reduction = "umap", dims = 1:2,
                          min_cells = 300, eps = 0.2, min_pts = 30,
                          suffix = "_clean", new_col = NULL, verbose = TRUE) {

  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is needed. Install it with install.packages('dbscan').")
  }
  meta <- seurat_obj@meta.data
  if (!col %in% colnames(meta)) stop("Column '", col, "' not found in the metadata.")
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Reduction '", reduction, "' not found in the object.")
  }

  emb <- Seurat::Embeddings(seurat_obj, reduction)[, dims, drop = FALSE]
  lab <- as.character(meta[[col]])

  if (is.null(eps)) {
    eps <- 0.01 * mean(apply(emb, 2, function(z) diff(range(z))))
  }
  if (verbose) message("eps = ", signif(eps, 3), ", min_pts = ", min_pts,
                       ", min_cells = ", min_cells)

  targets <- if (is.null(values)) sort(unique(lab[!is.na(lab)])) else values
  missing_vals <- setdiff(targets, unique(lab))
  if (length(missing_vals)) {
    stop("Label(s) not found in '", col, "': ", paste(missing_vals, collapse = ", "))
  }

  out <- lab
  rep <- data.frame(label = targets, before = NA_integer_, after = NA_integer_)

  for (k in seq_along(targets)) {
    idx <- which(lab == targets[k])
    rep$before[k] <- length(idx)

    if (length(idx) < min_pts) {
      out[idx] <- NA
      rep$after[k] <- 0L
      next
    }

    db   <- dbscan::dbscan(emb[idx, , drop = FALSE], eps = eps, minPts = min_pts)
    sz   <- table(db$cluster[db$cluster > 0])
    keep <- db$cluster %in% as.integer(names(sz)[sz >= min_cells])

    out[idx[!keep]] <- NA
    rep$after[k] <- sum(keep)
  }

  if (verbose) {
    rep$removed <- rep$before - rep$after
    message("Cells kept and removed per label:")
    print(rep, row.names = FALSE)
  }

  if (is.null(new_col)) new_col <- paste0(col, suffix)
  names(out) <- rownames(meta)
  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = out, col.name = new_col)
  seurat_obj
}

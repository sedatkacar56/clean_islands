# cleanislands

Remove small, isolated groups of cells ("islands") from a categorical metadata label on a 2D embedding (for example UMAP) of a Seurat object.

`clean_islands()` groups the cells of each label with DBSCAN and keeps only groups with at least `min_cells` cells. The result goes into a **new** metadata column. Your original column is never changed.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sedatkacar56/clean_islands")
```

Dependencies: [Seurat](https://satijalab.org/seurat/) and [dbscan](https://cran.r-project.org/package=dbscan).

## Usage

```r
library(cleanislands)

# clean a region label, new column "region_clean"
obj <- clean_islands(obj, col = "region")

# larger minimum size and a chosen radius
obj <- clean_islands(obj, "region", min_cells = 500, eps = 0.25)

# clean only one label of a cell type column, custom column name
obj <- clean_islands(obj, "celltype", values = "CD8_T", new_col = "celltype_pure")

# always check the result visually
Seurat::DimPlot(obj, group.by = "region_clean", raster = TRUE)
```

## Arguments

| Argument | Default | Meaning |
|---|---|---|
| `seurat_obj` | | Seurat object with the chosen reduction |
| `col` | | Metadata column to clean |
| `values` | `NULL` | Labels to clean. `NULL` cleans every non-`NA` label, each separately. Other labels are copied unchanged |
| `reduction` | `"umap"` | Reduction used for cell positions |
| `dims` | `1:2` | Dimensions of the reduction |
| `min_cells` | `300` | Smallest connected group that is kept |
| `eps` | `NULL` | DBSCAN radius. `NULL` = 1% of the average embedding width (about 0.15 to 0.2 for a typical UMAP) |
| `min_pts` | `30` | Neighbors within `eps` a cell needs to count as dense |
| `suffix` | `"_clean"` | Added to `col` to name the new column |
| `new_col` | `NULL` | Full name of the new column (overrides `suffix`) |
| `verbose` | `TRUE` | Print cells kept and removed per label |

## Behavior

- The original column is never overwritten. The output goes to a new column.
- Removed cells become `NA`. Cells that were `NA` stay `NA`.
- A cell is removed if it is DBSCAN noise (fewer than `min_pts` cells within `eps`) or belongs to a connected group smaller than `min_cells`.
- Only labels change. Nothing is tested, and expression is not used, only positions on the embedding.

## Tuning

The automatic `eps` is a rule of thumb. **Always check the result with `DimPlot`.**

- `eps` too small: one real region breaks into pieces and some are lost.
- `eps` too large: separate regions merge together.
- A thin streak disappears: lower `min_pts` (for example to 15) or raise `eps`.

## Performance

DBSCAN runs on every cell of a label. This is fine for 100k cells, but slow for several hundred thousand.

## Status

Early release (0.1.0). Run it on your own label first: the main blobs should stay and the tiny fragments should be gone.

## Acknowledgements

Alhamdulillah. Developed by Dr. Sedat Kacar.

## License

MIT, see `LICENSE.md`.

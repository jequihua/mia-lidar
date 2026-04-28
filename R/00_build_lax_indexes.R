# ============================================================================
# 00_build_lax_indexes.R
#
# One-shot utility: walk a LiDAR archive and build an .lax sidecar for every
# .las / .laz file that doesn't already have one. .lax indexes are tiny files
# (a few MB) that sit next to each LAS file on disk and turn buffered /
# spatially-bounded reads from full-file scans into targeted byte-range reads.
#
# Once built, .lax sidecars accelerate ALL future lidR operations on these
# files — rasterize_terrain, rasterize_canopy, pixel_metrics, clip_*,
# segment_trees, anything that reads less than the whole file. Built once,
# benefits forever.
#
# Cost-free in quality terms: .lax is purely a spatial index. The point cloud
# is unchanged. Lidr's own docs (vignette: "lidR-computation-speed-LAScatalog")
# call out that there is "no reason not to create them" — see notes/001.
#
# Usage:
#   1. Edit `lidar_root` below to point at the archive you want to index.
#   2. Run from the project root:  Rscript R/00_build_lax_indexes.R
#   3. Repeat per reserve (RBMN, RBRL, ...) by editing `lidar_root`.
#
# Behaviour:
#   - Recursive walk under `lidar_root`.
#   - Indexes every .las / .laz file that does not yet have a matching .lax.
#   - Skips files with an existing .lax (idempotent — safe to re-run).
#   - Per-file failures print a [fail] line and do not abort the batch.
# ============================================================================

library(lidR)

# --- config ----------------------------------------------------------------
lidar_root <- "D:/lidar_rbmnn"   # edit per archive: "D:/lidar_rbmn", etc.

# --- discovery -------------------------------------------------------------
las_files <- list.files(
  lidar_root,
  pattern       = "\\.la[sz]$",
  full.names    = TRUE,
  recursive     = TRUE,
  ignore.case   = TRUE
)
cat("Found", length(las_files), "LAS/LAZ file(s) under", lidar_root, "\n\n")

if (length(las_files) == 0L) {
  stop("No .las / .laz files found. Check `lidar_root`.")
}

# --- per-file worker -------------------------------------------------------
build_one <- function(f) {
  lax <- sub("\\.la[sz]$", ".lax", f, ignore.case = TRUE)
  if (file.exists(lax)) {
    cat("[skip]", basename(f), "- .lax exists\n")
    return(invisible(lax))
  }

  cat("[idx ]", basename(f), "... ")
  t0  <- Sys.time()
  ctg <- readLAScatalog(f)
  lidR:::catalog_laxindex(ctg)
  dt  <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("done in %.1fs\n", dt))
  invisible(lax)
}

# --- loop ------------------------------------------------------------------
# tryCatch so a single failed file doesn't abort the batch.
results <- vector("list", length(las_files))
for (i in seq_along(las_files)) {
  results[[i]] <- tryCatch(
    build_one(las_files[[i]]),
    error = function(e) {
      message("[fail] ", basename(las_files[[i]]), ": ", conditionMessage(e))
      NULL
    }
  )
}

ok <- sum(!vapply(results, is.null, logical(1)))
cat(sprintf("\nDone. %d / %d files have a .lax sidecar under %s\n",
            ok, length(las_files), lidar_root))

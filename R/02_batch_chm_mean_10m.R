# ============================================================================
# batch_chm_mean_10m.R
#
# Batch-generate 10 m mean-CHM rasters from every LiDAR triplet under a
# parent directory. Produces one GeoTIFF per area, in the LAS native CRS
# (usually UTM). Used to populate a mean-canopy-height reference layer over
# the whole LiDAR archive without running the full feature pipeline.
#
# Input layout (one folder per area, same naming convention as the smoke
# data — see helpers.R::find_lidar_triplets):
#   D:/lidar_rbmnn/<area_id>/<area_id>_TP.las
#   D:/lidar_rbmnn/<area_id>/<area_id>_GP.las
#   D:/lidar_rbmnn/<area_id>/<area_id>_AP.las
#
# Output:
#   data/rbmn_lidar_chm_mean_10m/<area_id>_chm_mean_10m.tif
#
# Per-area pipeline (minimal — CHM only):
#   AP -> 1 m CHM via p2r(subcircle = 0.2)           (triangulation-free,
#         works around DJI-L1 scale-factor issue, see notes/009)
#   aggregate -> 10 m mean
#   write GeoTIFF in LAS CRS
#
# Files are NOT projected to a common grid — each stays in its own LAS CRS.
# If downstream work needs wall-to-wall alignment (e.g. to the S2 grid), do
# that in a second step.
# ============================================================================

library(lidR)
library(terra)
library(sf)

source("R/helpers.R")

# polyfill
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
}

# --- config ----------------------------------------------------------------
lidar_root <- "D:/lidar_rbrl"
out_dir    <- "data/rbrl_lidar_chm_mean_10m"
overwrite  <- FALSE   # set TRUE to recompute areas that already have output

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- discovery -------------------------------------------------------------
triplets <- find_lidar_triplets(lidar_root)
cat("Found", length(triplets), "area(s) under", lidar_root, "\n\n")

# --- per-area worker -------------------------------------------------------
process_one <- function(trip) {
  out_file <- file.path(out_dir, sprintf("%s_chm_mean_10m.tif", trip$area_id))

  if (file.exists(out_file) && !overwrite) {
    cat("[skip]", trip$area_id, "- already exists\n")
    return(invisible(out_file))
  }

  cat("[proc]", trip$area_id, "... ")
  t0 <- Sys.time()

  # Catalog + filter (kill any sub-ground noise on AP).
  ctg_ap <- readLAScatalog(trip$ap)
  opt_filter(ctg_ap) <- "-drop_z_below 0"

  # Clean 10 m UTM template covering the AP footprint + small buffer.
  # (Same construction as block [6] in R/0_inspect_inputs.R.)
  bb     <- sf::st_bbox(ctg_ap)
  buf    <- 50
  xmin_r <- floor  ((bb[["xmin"]] - buf) / 10) * 10
  ymin_r <- floor  ((bb[["ymin"]] - buf) / 10) * 10
  xmax_r <- ceiling((bb[["xmax"]] + buf) / 10) * 10
  ymax_r <- ceiling((bb[["ymax"]] + buf) / 10) * 10

  ref <- terra::rast(
    xmin = xmin_r, xmax = xmax_r,
    ymin = ymin_r, ymax = ymax_r,
    resolution = 10,
    crs        = sf::st_crs(ctg_ap)$wkt
  )
  fine_ref <- terra::disagg(ref, fact = 10)

  # 1 m CHM via points-to-raster with subcircle, then aggregate to 10 m mean.
  chm_1m <- rasterize_canopy(
    ctg_ap,
    res       = fine_ref,
    algorithm = p2r(subcircle = 0.2)
  )
  chm_mean_10 <- terra::aggregate(chm_1m, fact = 10, fun = mean, na.rm = TRUE)
  names(chm_mean_10) <- "chm_mean"

  terra::writeRaster(
    chm_mean_10, out_file, overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")
  )

  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("done in %.1fs -> %s\n", dt, basename(out_file)))
  invisible(out_file)
}

# --- loop ------------------------------------------------------------------
# tryCatch so a single failed area doesn't abort the batch.
results <- vector("list", length(triplets))
for (i in seq_along(triplets)) {
  results[[i]] <- tryCatch(
    process_one(triplets[[i]]),
    error = function(e) {
      message("[fail] ", triplets[[i]]$area_id, ": ", conditionMessage(e))
      NULL
    }
  )
}

ok <- sum(!vapply(results, is.null, logical(1)))
cat(sprintf("\nDone. %d / %d areas produced a CHM file in %s\n",
            ok, length(triplets), out_dir))

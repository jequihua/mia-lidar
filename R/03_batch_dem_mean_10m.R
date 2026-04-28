# ============================================================================
# 03_batch_dem_mean_10m.R
#
# Batch-generate 10 m mean-DEM rasters from every LiDAR triplet under a
# parent directory. Produces one GeoTIFF per area, in the LAS native CRS
# (usually UTM). Used to populate a mean-ground-level reference layer over
# the whole LiDAR archive without running the full feature pipeline.
#
# Counterpart to R/02_batch_chm_mean_10m.R: same shape, but reads GP (Ground
# Points) instead of AP, and interpolates the ground surface rather than
# rasterising the canopy top.
#
# Caveat (notes/009): on DJI Zenmuse L1 / DJI Terra exports, GP is
# height-normalised — Z is height-above-the-vendor's-ground-reference, not
# elevation-ASL. The output here is therefore a *relative* ground
# micro-topography raster (~1 m of relief in mangrove sites), not an
# absolute DEM. For absolute elevation we still need an external DEM
# (SRTM / national DEM).
#
# Input layout (one folder per area, same naming convention as
# helpers.R::find_lidar_triplets):
#   <lidar_root>/<area_id>/<area_id>_TP.las
#   <lidar_root>/<area_id>/<area_id>_GP.las
#   <lidar_root>/<area_id>/<area_id>_AP.las
#
# Output:
#   out/lidar/lidar/<reserve>/dem_mean_10m/<area_id>_dem_mean_10m.tif
#
# Per-area pipeline (minimal — DEM only):
#   GP -> 1 m DTM via knnidw(k = 10, p = 2)          (triangulation-free,
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
# Edit `reserve` and `lidar_root` together; the output directory derives
# from `reserve` so the layout matches out/lidar/lidar/<reserve>/chm_mean_10m.
reserve    <- "rbmn"            # one of: "rbmn", "rbrl"
lidar_root <- "D:/lidar_rbmn"   # adjust per reserve
overwrite  <- FALSE              # set TRUE to recompute areas that already have output

out_dir <- file.path("out/lidar/lidar", reserve, "dem_mean_10m")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- discovery -------------------------------------------------------------
triplets <- find_lidar_triplets(lidar_root)
cat("Found", length(triplets), "area(s) under", lidar_root, "\n\n")

# --- per-area worker -------------------------------------------------------
process_one <- function(trip) {
  out_file <- file.path(out_dir, sprintf("%s_dem_mean_10m.tif", trip$area_id))

  if (file.exists(out_file) && !overwrite) {
    cat("[skip]", trip$area_id, "- already exists\n")
    return(invisible(out_file))
  }

  cat("[proc]", trip$area_id, "... ")
  t0 <- Sys.time()

  # Catalog. GP is supposed to be ground-only (FUSION-classified), but the
  # vendor leaves ASPRS Classification = 0 on every point (notes/009), so
  # rasterize_terrain's default use_class = c(2L, 9L) would reject the data.
  # Pull the present classes once and pass them explicitly.
  ctg_gp <- readLAScatalog(trip$gp)
  gp_classes <- sort(unique(readLAS(trip$gp, select = "c")@data$Classification))

  # Modest chunk buffer so the interpolator has neighbour points at chunk
  # boundaries (rasterize_terrain is one of the functions that benefits).
  opt_chunk_buffer(ctg_gp) <- 30

  # Clean 10 m UTM template covering the GP footprint + small buffer.
  # (Same construction as block [6] in R/01_inspect_inputs.R and the CHM
  # batch script.)
  bb     <- sf::st_bbox(ctg_gp)
  buf    <- 50
  xmin_r <- floor  ((bb[["xmin"]] - buf) / 10) * 10
  ymin_r <- floor  ((bb[["ymin"]] - buf) / 10) * 10
  xmax_r <- ceiling((bb[["xmax"]] + buf) / 10) * 10
  ymax_r <- ceiling((bb[["ymax"]] + buf) / 10) * 10

  ref <- terra::rast(
    xmin = xmin_r, xmax = xmax_r,
    ymin = ymin_r, ymax = ymax_r,
    resolution = 10,
    crs        = sf::st_crs(ctg_gp)$wkt
  )
  fine_ref <- terra::disagg(ref, fact = 10)

  # 1 m DTM via knnidw, then aggregate to 10 m mean.
  # NB: tin() / kriging() would be alternatives; tin() is blocked by the DJI
  # scale-factor issue (notes/009); kriging() is far slower and offers little
  # for FUSION-cleaned ground points at this density.
  dtm_1m <- rasterize_terrain(
    ctg_gp,
    res       = fine_ref,
    algorithm = knnidw(k = 10, p = 2),
    use_class = gp_classes
  )

  dem_mean_10 <- terra::aggregate(dtm_1m, fact = 10, fun = mean, na.rm = TRUE)
  names(dem_mean_10) <- "dem_mean"

  terra::writeRaster(
    dem_mean_10, out_file, overwrite = TRUE,
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
cat(sprintf("\nDone. %d / %d areas produced a DEM file in %s\n",
            ok, length(triplets), out_dir))

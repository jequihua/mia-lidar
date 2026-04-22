# ============================================================================
# 0_inspect_inputs.R  —  Manual walkthrough of the biomass-from-lidar pipeline
#
# Run block by block (Ctrl+Enter / Cmd+Enter). Each block is marked with a
# banner like "--- [N] SHORT TITLE ---". Outputs print or plot as we go so we
# can sanity-check before factoring into a batch-processing wrapper.
#
# Pipeline (see notes/001, 002, 003 for rationale):
#   [1]  packages + paths + helpers
#   [2]  inspect reference rasters
#   [3]  find + inspect the TP/GP/AP triplet
#   [4]  CRS alignment check
#   [5]  build .lax sidecars (one-time speedup)
#   [6]  define the reference grid (coarse + fine-nested)
#   [7]  DTM from GP + terrain derivatives   -> ref grid
#   [8]  CHM from AP + aggregates            -> ref grid
#   [9]  ABA height metrics from AP          on ref grid
#   [10] pground + return-number stats from TP on ref grid
#   [11] stack and write all features
#
# Run from the project root (so relative paths resolve).
# ============================================================================

# --- [1] packages, paths, helpers ------------------------------------------
library(lidR)
library(terra)
library(sf)

source("R/helpers.R")

# polyfill for R < 4.4
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
}

rs_dir     <- "data/rs"
lidar_root <- "data/lidar"
area_dir   <- file.path(lidar_root, "DJI_202505061638_005_Zenmuse-L1_RBMNN_CC_01")
out_root   <- "out/features"

rs_files <- list.files(rs_dir, pattern = "\\.tif$",
                       full.names = TRUE, ignore.case = TRUE)
stopifnot(length(rs_files) > 0)

trip    <- find_lidar_triplet(area_dir)
out_dir <- file.path(out_root, trip$area_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Area:   ", trip$area_id, "\n")
cat("TP:     ", basename(trip$tp), "\n")
cat("GP:     ", basename(trip$gp), "\n")
cat("AP:     ", basename(trip$ap), "\n")
cat("Output: ", out_dir, "\n\n")

# --- [2] inspect reference rasters -----------------------------------------

inspect_raster <- function(path) {
  cat("-----------------------------------------------------------------\n")
  cat("RASTER:", basename(path), "\n")
  cat("-----------------------------------------------------------------\n")

  r <- terra::rast(path)
  cat("CRS:        ", terra::crs(r, describe = TRUE)$name,
      "(EPSG:", terra::crs(r, describe = TRUE)$code, ")\n")
  cat("Extent:     ", paste(round(as.vector(terra::ext(r)), 2), collapse = ", "), "\n")
  cat("Resolution: ", paste(terra::res(r), collapse = " x "), "\n")
  cat("Dims:       ", nrow(r), "x", ncol(r), "x", terra::nlyr(r), "\n")
  cat("Layer(s):   ", paste(names(r), collapse = ", "), "\n")
  stats <- terra::global(r, fun = c("min", "mean", "max", "isNA"), na.rm = TRUE)
  print(round(stats, 4))
  cat("\n")
  invisible(r)
}

rs_list <- lapply(rs_files, inspect_raster)
names(rs_list) <- basename(rs_files)

# --- [3] inspect the TP/GP/AP triplet --------------------------------------

inspect_lidar_source <- function(file, label, expect_normalized = FALSE) {
  cat("-----------------------------------------------------------------\n")
  cat("LIDAR:", label, "\n")
  cat("-----------------------------------------------------------------\n")

  ctg <- readLAScatalog(file)
  print(ctg); cat("\n")
  cat("las_check(catalog):\n"); las_check(ctg); cat("\n")

  # Direct read for deeper inspection. Fine for smoke files; we'll replace
  # this with header-only checks once we batch.
  las <- readLAS(file)
  cat("Attributes:        ", paste(names(las@data), collapse = ", "), "\n")
  cat("Classification counts (ASPRS 1=unclass, 2=ground, 3-5=veg, 7=noise, 9=water):\n")
  print(table(las@data$Classification, useNA = "ifany"))
  cat("ReturnNumber counts:\n")
  print(table(las@data$ReturnNumber, useNA = "ifany"))

  zr <- range(las@data$Z)
  cat(sprintf("Z range:            [%.2f, %.2f]   spread = %.2f m\n",
              zr[1], zr[2], diff(zr)))

  # AP is supposed to be height-normalized; any Z < 0 is sub-ground noise
  # (multipath etc.) we'll filter downstream via opt_filter("-drop_z_below 0").
  if (expect_normalized) {
    below_zero <- sum(las@data$Z < 0)
    pct        <- 100 * below_zero / nrow(las@data)
    cat(sprintf("Points with Z < 0:  %d  (%.3f%%) — will drop via opt_filter\n",
                below_zero, pct))
  }
  cat(sprintf("Point density:      %.2f pts/m^2 (from header)\n", density(ctg)))
  cat("\n")

  rm(las); invisible(gc(verbose = FALSE))
  invisible(ctg)
}

ctg_tp <- inspect_lidar_source(
  trip$tp,
  "TP — Total Points (ground + vegetation, NOT normalized)")

ctg_gp <- inspect_lidar_source(
  trip$gp,
  "GP — Ground Points (DTM source)")

ctg_ap <- inspect_lidar_source(
  trip$ap,
  "AP — Above Points (height-normalized)",
  expect_normalized = TRUE)

# --- [4] CRS alignment check -----------------------------------------------
# lidR does not reproject points. If the ref rasters and LAS disagree on CRS
# we need to reproject the ref raster into the LAS CRS before using it as a
# grid template (or rewrite the LAS files). Flag + warn here.

cat("-----------------------------------------------------------------\n")
cat("CRS ALIGNMENT\n")
cat("-----------------------------------------------------------------\n")

rs_crs    <- lapply(rs_list, function(r) sf::st_crs(terra::crs(r)))
lidar_crs <- sf::st_crs(ctg_ap)

for (i in seq_along(rs_crs)) {
  cat(sprintf("raster [%s]  EPSG=%s\n",
              names(rs_list)[i], rs_crs[[i]]$epsg %||% NA))
}
cat(sprintf("lidar  (AP)  EPSG=%s\n", lidar_crs$epsg %||% NA))

rs_match <- vapply(rs_crs, function(c) c == lidar_crs, logical(1))
cat("raster vs lidar CRS match: ", paste(rs_match, collapse = ", "), "\n\n")

if (any(!rs_match)) {
  warning(
    "At least one raster CRS differs from the LAS CRS. The downstream blocks ",
    "assume they match. Either reproject the ref raster into the LAS CRS ",
    "before continuing, or update the reference-grid block."
  )
}

# Save the LAS footprint reprojected into the S2 CRS as a GeoJSON. Drop this
# into QGIS / geojson.io / Leaflet to confirm visually that the LAS data
# lies where the S2 raster claims to cover.
ref_path_for_crs <- rs_files[grepl("s2", rs_files, ignore.case = TRUE)][1]
footprint_path <- file.path(out_dir, sprintf("%s_footprint.geojson", trip$area_id))
save_lidar_footprint(
  ctg_ap,
  target_crs = terra::rast(ref_path_for_crs),
  out_path   = footprint_path,
  area_id    = trip$area_id
)
cat("LAS footprint (in S2 CRS) written to:\n  ", footprint_path, "\n\n")

# --- [5] build .lax sidecars -----------------------------------------------
# One-time cost; 2-3x faster buffered reads thereafter (notes/001, speed
# section). `catalog_laxindex` is an internal lidR function — safe to use,
# but noted here so it's not a surprise.

lax_exists <- function(las_path) {
  file.exists(sub("\\.la[sz]$", ".lax", las_path, ignore.case = TRUE))
}

if (!all(vapply(c(trip$tp, trip$gp, trip$ap), lax_exists, logical(1)))) {
  cat("Building .lax indexes...\n")
  lidR:::catalog_laxindex(ctg_tp)
  lidR:::catalog_laxindex(ctg_gp)
  lidR:::catalog_laxindex(ctg_ap)
  cat("done.\n\n")
} else {
  cat(".lax sidecars already present — skipping.\n\n")
}

# --- [6] reference grid ----------------------------------------------------
# Pick one raster as the canonical reference. S2 NDVI is a natural choice:
# native 10 m and typically the grid the rest of the satellite stack is
# already aligned to.
#
# CRS handling: lidR does not reproject points. If the S2 grid is in a CRS
# different from the LAS CRS, we do the lidR work in the LAS CRS and
# reproject the final feature stack onto the S2 grid at write time. That
# keeps LiDAR outputs unambiguously aligned with the point cloud (no resample
# error inside feature computation), while the final product ends up on the
# S2 grid for clean stacking with S1/S2 downstream.
#
# For the fine-then-aggregate pattern (notes/002) we build a 10x-finer
# template perfectly nested inside the coarse one via terra::disagg — same
# extent, cells subdivided by factor 10.

ref_path <- rs_files[grepl("s2", rs_files, ignore.case = TRUE)][1]
stopifnot(!is.na(ref_path))

# ref_out: the final alignment target (original S2 CRS/extent/resolution).
ref_out    <- terra::rast(ref_path)
ref_out_crs <- sf::st_crs(terra::crs(ref_out))

# ref: the lidR-facing template. Two cases:
#   (a) CRSs agree — use the S2 grid directly (rare in our setup).
#   (b) CRSs differ — build a clean 10 m grid in the LAS CRS covering the
#       LAS footprint plus a small buffer. We do NOT project the S2 grid
#       via terra::project: the degree-based S2 yields non-integer
#       resolution (~9.8 m) and a giant extent covering the whole S2 tile
#       (billions of fine-grid cells), neither of which we want.
if (ref_out_crs == lidar_crs) {
  ref <- ref_out
  cat("Reference grid CRS matches LiDAR CRS — using S2 grid directly.\n")
} else {
  cat(sprintf("Reference grid CRS (EPSG:%s) differs from LiDAR CRS (EPSG:%s).\n",
              ref_out_crs$epsg %||% "??", lidar_crs$epsg %||% "??"))
  cat("Building a clean 10 m grid in LAS CRS covering the LAS footprint.\n")

  # LAS footprint + buffer, snapped to a 10 m grid aligned to integer coords.
  ctg_bbox <- sf::st_bbox(ctg_ap)
  buf      <- 50   # metres, small pad to avoid edge loss during rasterize
  xmin_r   <- floor  ((ctg_bbox[["xmin"]] - buf) / 10) * 10
  ymin_r   <- floor  ((ctg_bbox[["ymin"]] - buf) / 10) * 10
  xmax_r   <- ceiling((ctg_bbox[["xmax"]] + buf) / 10) * 10
  ymax_r   <- ceiling((ctg_bbox[["ymax"]] + buf) / 10) * 10

  ref <- terra::rast(
    xmin = xmin_r, xmax = xmax_r,
    ymin = ymin_r, ymax = ymax_r,
    resolution = 10,
    crs        = sf::st_crs(ctg_ap)$wkt
  )
}
fine_ref <- terra::disagg(ref, fact = 10)

cat("Reference grid:", basename(ref_path), "\n")
cat("  ref_out (final alignment, original S2 CRS): ",
    paste(terra::res(ref_out), collapse = " x "),
    " dims:", nrow(ref_out), "x", ncol(ref_out),
    " EPSG:", ref_out_crs$epsg %||% "??", "\n")
cat("  ref     (lidR template, LAS CRS):            ",
    paste(terra::res(ref), collapse = " x "),
    " dims:", nrow(ref),     "x", ncol(ref),
    " EPSG:", sf::st_crs(terra::crs(ref))$epsg %||% "??", "\n")
cat("  fine_ref (10x inside ref):                   ",
    paste(terra::res(fine_ref), collapse = " x "),
    " dims:", nrow(fine_ref), "x", ncol(fine_ref), "\n\n")

# --- [7] DTM from GP + terrain derivatives ---------------------------------
# Build DTM at fine resolution using TIN (fast, accurate for dense ALS).
# Compute slope/aspect/TPI/TRI/roughness at fine resolution, then aggregate
# to the ref grid. Computing these derivatives at native 10 m would smooth
# out the small-scale relief that tends to matter for biomass in
# heterogeneous terrain (notes/002).
#
# GP points may not be classified as class 2 by the vendor's workflow (some
# DJI-derived exports leave class = 1). We read the present classes once and
# pass them via use_class so rasterize_terrain doesn't reject the data.

gp_classes <- sort(unique(readLAS(trip$gp, select = "c")@data$Classification))
cat("[7] GP classifications present:", paste(gp_classes, collapse = ", "), "\n")

cat("    building fine DTM (kNN-IDW)...\n")
# NB: DJI Zenmuse L1 exports use non-standard scale factors (0.0001), which
# breaks lidR's C++ Delaunay path (tin / pitfree / dsmtin all fail with
# "xy coordinates were not converted to integer"). See notes/009.
# knnidw() does not triangulate, so it works with the data as-delivered.
# k = 10 + p = 2 is a sensible default for dense ALS ground points; quality
# difference vs tin() is small at 1 m for a 4 pts/m² ground density.
opt_chunk_buffer(ctg_gp) <- 30
dtm_1m <- rasterize_terrain(ctg_gp, res = fine_ref,
                            algorithm = knnidw(k = 10, p = 2),
                            use_class = gp_classes)

cat("    terrain derivatives at fine resolution...\n")
terr_1m <- terra::terrain(dtm_1m,
                          v = c("slope", "aspect", "TPI", "TRI", "roughness"),
                          unit = "degrees")

cat("    aggregating to ref grid...\n")
agg <- function(r, fun) terra::aggregate(r, fact = 10, fun = fun, na.rm = TRUE)

dtm_10       <- agg(dtm_1m,                 mean)
slope_10     <- agg(terr_1m[["slope"]],     mean)
aspect_10    <- agg(terr_1m[["aspect"]],    mean)
tpi_10       <- agg(terr_1m[["TPI"]],       mean)
tri_10       <- agg(terr_1m[["TRI"]],       mean)
rough_10     <- agg(terr_1m[["roughness"]], mean)
slope_sd_10  <- agg(terr_1m[["slope"]],     sd)   # ruggedness-of-slope

terrain_stack <- c(dtm_10, slope_10, aspect_10, tpi_10, tri_10, rough_10, slope_sd_10)
names(terrain_stack) <- c("dtm", "slope", "aspect", "tpi", "tri", "roughness", "slope_sd")
print(terrain_stack)
cat("\n")

# --- [8] CHM from AP + aggregates ------------------------------------------
# Khosravipour pit-free CHM at fine resolution, aggregated with several stats
# (notes/002). Default thresholds c(0, 2, 5, 10, 15) m and max_edge c(0, 1.5) m
# — revisit if the AP Z quantiles suggest a very tall or very short canopy.

cat("[8] building CHM from AP (p2r with subcircle)...\n")
# Same scale-factor issue as DTM: pitfree / dsmtin trigger the Delaunay
# integer-conversion check and fail on DJI L1 exports (notes/009).
# p2r with a small subcircle replaces each point with a disc of synthetic
# points, which fills sub-pixel pits without triangulation. At 124 pts/m²
# AP density and 1 m resolution, this produces a CHM of comparable quality
# to pitfree for our purposes.
opt_filter(ctg_ap) <- "-drop_z_below 0"    # kill sub-ground noise
# rasterize_canopy does not need a chunk buffer.

chm_1m <- rasterize_canopy(
  ctg_ap, res = fine_ref,
  algorithm = p2r(subcircle = 0.2)
)

cat("    aggregating CHM to ref grid...\n")
chm_max  <- agg(chm_1m,       max)     # dominant height
chm_mean <- agg(chm_1m,       mean)    # mean canopy height
chm_sd   <- agg(chm_1m,       sd)      # canopy rugosity
cov_gt2  <- agg(chm_1m > 2,   mean)    # fractional cover > 2 m
cov_gt5  <- agg(chm_1m > 5,   mean)    # fractional cover > 5 m

chm_stack <- c(chm_max, chm_mean, chm_sd, cov_gt2, cov_gt5)
names(chm_stack) <- c("chm_max", "chm_mean", "chm_sd", "cov_gt2m", "cov_gt5m")
print(chm_stack)
cat("\n")

# --- [9] ABA height metrics from AP ----------------------------------------
# Point-cloud metrics are computed directly at ref resolution (notes/002).
# .stdmetrics_z returns zmean, zsd, zskew, zkurt, zq5..zq95, zpcum1..9,
# pzabove2, pzabovemean.

cat("[9] .stdmetrics_z from AP on ref grid...\n")
# opt_filter("-drop_z_below 0") was set in block [8] and is still active.
aba_z <- pixel_metrics(ctg_ap, .stdmetrics_z, res = ref)
print(aba_z)
cat("\n")

# --- [10] pground + return-number stats from TP ----------------------------
# Only TP preserves the full return profile (notes/003). .stdmetrics_rn
# reports counts and proportions of 1st/2nd/3rd+ returns + ground percent
# (pground). No height normalization needed for these metrics.

cat("[10] .stdmetrics_rn from TP on ref grid...\n")
aba_rn <- pixel_metrics(ctg_tp, .stdmetrics_rn, res = ref)
print(aba_rn)
cat("\n")

# TODO: canonical first-return canopy cover (n_first>h / n_first_all, per
# White et al. 2013) needs height above ground, i.e. TP normalized by the
# DTM. Options: merge_spatial(dtm) inside a custom pixel_metrics, or write a
# normalized TP and reuse the AP pattern. Deferred for its own discussion.

# --- [11] project each group to S2 grid, stack, and write ------------------
# Each lidR operation clips to its own point-cloud hull and inherits whatever
# origin the intermediate 1 m grid happened to land on — so terrain_stack,
# chm_stack, aba_z and aba_rn end up on slightly different 10 m extents
# (same CRS and resolution, but different origins). terra::c() refuses to
# stack them.
#
# Fix: project each group onto a common target. We use `ref_out_sub`, which
# is `ref_out` (the S2 grid) cropped to the LAS footprint — gives a small
# output aligned to the S2 grid, suitable for stacking with S1/S2. All four
# projections snap to the same extent/grid, so the final c() works.
# Bilinear is correct for continuous numeric layers; no categorical layers.

cat("[11] building S2-grid target cropped to LAS footprint...\n")
ctg_bbox_las      <- sf::st_as_sfc(sf::st_bbox(ctg_ap))
ctg_bbox_out_poly <- sf::st_transform(ctg_bbox_las, ref_out_crs)
ctg_bbox_out      <- sf::st_bbox(ctg_bbox_out_poly)

cat(sprintf("    LAS bbox in ref CRS:  xmin=%.4f xmax=%.4f ymin=%.4f ymax=%.4f\n",
            ctg_bbox_out[["xmin"]], ctg_bbox_out[["xmax"]],
            ctg_bbox_out[["ymin"]], ctg_bbox_out[["ymax"]]))
cat(sprintf("    ref_out bbox:         xmin=%.4f xmax=%.4f ymin=%.4f ymax=%.4f\n",
            terra::xmin(ref_out), terra::xmax(ref_out),
            terra::ymin(ref_out), terra::ymax(ref_out)))

# Overlap test.
ref_out_bbox_poly <- sf::st_as_sfc(sf::st_bbox(ref_out))
overlap <- as.logical(
  sf::st_intersects(ctg_bbox_out_poly, ref_out_bbox_poly, sparse = FALSE)
)

if (!overlap) {
  stop(
    "\nLAS footprint does NOT overlap the reference S2 raster.\n",
    "Options:\n",
    "  (a) swap in an S2 tile that actually covers the LAS area,\n",
    "  (b) check the LAS CRS metadata (UTM zone, hemisphere) is correct,\n",
    "  (c) for a smoke test with deliberately non-overlapping data, write\n",
    "      the features in LAS CRS directly (edit this block to skip the\n",
    "      project-to-S2-grid step).\n"
  )
}

crop_ext <- terra::ext(
  ctg_bbox_out[["xmin"]], ctg_bbox_out[["xmax"]],
  ctg_bbox_out[["ymin"]], ctg_bbox_out[["ymax"]]
)
ref_out_sub <- terra::crop(ref_out, crop_ext, snap = "out")
cat("    target dims:", nrow(ref_out_sub), "x", ncol(ref_out_sub), "\n")

cat("    projecting each group...\n")
terrain_out <- terra::project(terrain_stack, ref_out_sub, method = "bilinear")
chm_out     <- terra::project(chm_stack,     ref_out_sub, method = "bilinear")
aba_z_out   <- terra::project(aba_z,         ref_out_sub, method = "bilinear")
aba_rn_out  <- terra::project(aba_rn,        ref_out_sub, method = "bilinear")

cat("    stacking...\n")
features <- c(terrain_out, chm_out, aba_z_out, aba_rn_out)
cat("    ", terra::nlyr(features), "layers:\n    ",
    paste(names(features), collapse = ", "), "\n")

out_file <- file.path(out_dir, sprintf("%s_features.tif", trip$area_id))
cat("    writing ->", out_file, "\n")
terra::writeRaster(
  features, out_file, overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")
)

cat("\nDone.\n")

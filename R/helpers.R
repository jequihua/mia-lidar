# ---------------------------------------------------------------------------
# helpers.R  —  project-wide helpers not tied to the pipeline logic.
# Source with `source("R/helpers.R")`.
# ---------------------------------------------------------------------------

# find_lidar_triplet(area_dir)
#
# Our LiDAR deliveries arrive as three co-located point clouds per area, with
# `_TP.las`, `_GP.las`, `_AP.las` suffixes in the filename (see notes/003).
# This helper locates the three files in one area directory and returns a
# tagged list, failing loudly if the triplet is incomplete.
#
# Returns:
#   list(area_id, dir, tp, gp, ap)   — all paths absolute-normalized.
find_lidar_triplet <- function(area_dir) {
  if (!dir.exists(area_dir)) {
    stop("area_dir does not exist: ", area_dir)
  }
  area_dir <- normalizePath(area_dir, winslash = "/", mustWork = TRUE)

  all_las <- list.files(area_dir, pattern = "\\.la[sz]$",
                        full.names = TRUE, ignore.case = TRUE)
  pick <- function(tag) {
    m <- grep(paste0("_", tag, "\\.la[sz]$"),
              all_las, ignore.case = TRUE, value = TRUE)
    if (length(m) != 1L) {
      stop(sprintf(
        "Expected exactly 1 *_%s.(las|laz) in %s, found %d.",
        tag, area_dir, length(m)))
    }
    m
  }

  list(
    area_id = basename(area_dir),
    dir     = area_dir,
    tp      = pick("TP"),
    gp      = pick("GP"),
    ap      = pick("AP")
  )
}

# find_lidar_triplets(parent_dir)
#
# Batch variant. If `parent_dir` contains subfolders, each is treated as an
# area and `find_lidar_triplet` is applied to it. If `parent_dir` itself
# contains the _TP/_GP/_AP files, it is treated as a single area.
#
# Returns: list of triplet lists, each as from find_lidar_triplet().
find_lidar_triplets <- function(parent_dir) {
  if (!dir.exists(parent_dir)) {
    stop("parent_dir does not exist: ", parent_dir)
  }
  subdirs <- list.dirs(parent_dir, recursive = FALSE)
  if (length(subdirs) == 0L) {
    return(list(find_lidar_triplet(parent_dir)))
  }
  # If the parent itself has a triplet, process it alongside subfolders.
  has_top_level_triplet <- length(list.files(
    parent_dir, pattern = "_(TP|GP|AP)\\.la[sz]$", ignore.case = TRUE)) >= 3L
  candidates <- if (has_top_level_triplet) c(parent_dir, subdirs) else subdirs

  triplets <- list()
  for (d in candidates) {
    triplets[[length(triplets) + 1L]] <- tryCatch(
      find_lidar_triplet(d),
      error = function(e) {
        message("  skipping ", d, ": ", conditionMessage(e))
        NULL
      }
    )
  }
  Filter(Negate(is.null), triplets)
}

# save_lidar_footprint(ctg, target_crs, out_path, area_id = NULL)
#
# Extract the LAS catalog / LAS file footprint as a rectangular polygon,
# reproject it into a target CRS, and write it to disk as a georeferenced
# vector file. Useful for (a) QA of CRS metadata (drop the output in QGIS
# / a web map and confirm the LAS sits where you expect), (b) building an
# index of processed areas over a large project.
#
# Args:
#   ctg        : LAScatalog or LAS object with a valid CRS.
#   target_crs : target CRS. One of: sf::crs object, WKT string, EPSG code
#                (numeric), or a SpatRaster whose CRS will be used.
#   out_path   : output file path. Extension determines format:
#                .geojson, .shp, .gpkg, .fgb, etc. Parent dir is created.
#   area_id    : optional label stored as an attribute column on the feature.
#
# Returns the sf object (invisibly).
save_lidar_footprint <- function(ctg, target_crs, out_path, area_id = NULL) {
  bbox_las <- sf::st_as_sfc(sf::st_bbox(ctg))

  # Accept a SpatRaster directly for convenience.
  if (inherits(target_crs, "SpatRaster")) {
    target_crs <- sf::st_crs(terra::crs(target_crs))
  }

  bbox_target <- sf::st_transform(bbox_las, target_crs)

  src_epsg <- tryCatch(sf::st_crs(ctg)$epsg, error = function(e) NA_integer_)
  attrs <- data.frame(
    area_id  = if (is.null(area_id)) "" else area_id,
    src_epsg = if (is.null(src_epsg)) NA_integer_ else src_epsg,
    stringsAsFactors = FALSE
  )
  sf_obj <- sf::st_sf(attrs, geometry = bbox_target)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(sf_obj, out_path, delete_dsn = TRUE, quiet = TRUE)
  invisible(sf_obj)
}

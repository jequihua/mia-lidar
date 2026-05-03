# ============================================================================
# 04_batch_teacher_grid_v1.R
#
# Build the v1 wall-to-LiDAR consensus dominant-height teacher labels at 10 m
# for every LiDAR triplet under one protected-area root. Output is the dense
# per-pixel label raster the student model is trained against, plus the two
# source-side estimates and full agreement QA. NOT projected to a common S2
# grid -- each parcel stays in its LAS-native CRS (typically UTM).
#
# Spec source: instructions/02_lidar_teacher_products.md and
#              instructions/feature_stack_decision_note.md.
#
# Per pixel u (10 m cell), the active teacher target is:
#   H_dom_consensus(u) = mean(chm_p95_local(u), zq95_local(u))    if qc=1
#                      = NA                                       otherwise
# with the agreement rule (v1.1 -- empirically re-cut, see notes/012):
#   tol(u)       = max(1.5 m, 0.15 * (chm_p95 + zq95) / 2)
#   within  = both present and |chm_p95 - zq95| <= tol
#   outside = both present and |chm_p95 - zq95| >  tol
#
# v1.1 vs v1: the agreement rule's numeric constants were re-cut from
# (0.5 m, 10%) to (1.5 m, 15%) on 2026-05-03 (architect approval, see
# notes/012 and `out/teacher/architect_request_agreement_rule_recut.md`).
# The structure of the rule (absolute floor + relative ceiling) is
# unchanged. The rationale is the empirical statistical offset of ~1 m
# between `Q_{0.95}(per-1 m max)` (chm_p95) and `Q_{0.95}(returns)`
# (zq95) on dense DJI Zenmuse L1 point clouds (~124 pts/m^2). The
# calibration is dataset-specific and should be re-audited for
# materially different LiDAR acquisitions; see notes/011 §3.8.
#
# Per-pixel `agreement_reason` codes (kept separate from `agreement_qc` so we
# do not collapse "actively disagree" into "untestable"):
#   1 = within_tolerance
#   2 = outside_tolerance
#   3 = missing_chm   (chm NA, zq present)
#   4 = missing_zq    (chm present, zq NA)
#   5 = missing_both
#
# `agreement_qc` keeps the binary downstream signal separately:
#   1   = consensus teacher can be formed
#   0   = teacher is withheld pending QA (active disagreement)
#   NA  = teacher cannot be formed (untestable; coverage/missingness)
#
# `chm_p95` minimum-support rule (closes architect review finding P2 #2):
#   Each 10 m cell is supported by 100 underlying 1 m sub-cells. With
#   `aggregate(..., na.rm = TRUE)` alone, a 10 m cell with only a handful of
#   valid 1 m sub-cells would still produce a "valid-looking" chm_p95 from a
#   non-representative sample. To preserve missingness honestly, this script
#   requires a minimum number of valid 1 m sub-cells per 10 m cell
#   (`min_valid_subcells`, default `80` of `100`). Cells failing the rule are
#   set to NA on `chm_p95`, which then propagates through the agreement rule
#   as `agreement_reason = 3 (missing_chm)` and `agreement_qc = NA`. The
#   per-cell support count is also archived as `chm_p95_count` for full
#   coverage transparency.
#
# Archived secondary LiDAR feature library (closes architect review finding
# P2 #1, per `instructions/02_lidar_teacher_products.md` Output #2):
#   The script also produces, on the same 10 m support, the archived secondary
#   feature library required by the frozen contract. They are not part of the
#   active v1 teacher target, but are archived now to support QA, expert
#   review, support-sensitivity analysis, and bounded later ablations:
#     - chm_mean         (mean of 1 m CHM, partial-support OK -- this is the
#                         summary-of-summaries variant; see note below)
#     - chm_iqr          (q75 - q25 of 1 m CHM)
#     - cov_gt2          (CHM occupancy: mean(chm_1m > 2 m) over valid 1 m
#                         sub-cells per 10 m cell -- canonical frozen form)
#     - cov_gt5          (CHM occupancy: mean(chm_1m > 5 m) over valid 1 m
#                         sub-cells per 10 m cell -- canonical frozen form)
#     - zq50             (point-cloud: 50th percentile of Z)
#     - zq95             (point-cloud: 95th percentile of Z -- same raster as
#                         the teacher-side `zq95/` output; not duplicated)
#     - zmax_minus_zq95  (point-cloud: max(Z) - q95(Z))
#   Cover features are CHM-occupancy based (mean of `chm_1m > h` over the 100
#   underlying 1 m sub-cells per 10 m cell). Architect decision: this matches
#   the "cover above threshold" semantics in the project materials and is
#   less sensitive to point-density and acquisition-pattern artifacts than
#   the AP-point occupancy form `mean(Z > h)`. The point-cloud variant is
#   not produced in this slice; if a later ablation needs it, add it as a
#   separately named metric (e.g. `pc_cov_gt2`), never overloading the
#   canonical name.
#
# Note on `chm_mean` partial support: unlike `chm_p95`, `chm_mean` does NOT
# enforce the min-valid-subcell rule. Mean is robust to subsampling in a way
# that a 95th percentile is not, and `chm_mean` is an archived QA/ablation
# layer rather than the active teacher input. The companion `chm_p95_count`
# raster lets a downstream consumer apply any stricter masking if needed.
#
# Per-area output layout (one folder per protected area):
#   out/teacher/<area>/
#     chm_p95/<parcel_id>_chm_p95_10m.tif
#     chm_p95_count/<parcel_id>_chm_p95_count_10m.tif   (Byte, 0..100)
#     zq95/<parcel_id>_zq95_10m.tif
#     agreement/<parcel_id>_agreement_abs_diff_10m.tif
#     agreement_qc/<parcel_id>_agreement_qc_10m.tif
#     agreement_reason/<parcel_id>_agreement_reason_10m.tif
#     teacher_labels/<parcel_id>_H_dom_consensus_10m.tif
#     secondary/<parcel_id>_chm_mean_10m.tif
#     secondary/<parcel_id>_chm_iqr_10m.tif
#     secondary/<parcel_id>_zq50_10m.tif
#     secondary/<parcel_id>_zmax_minus_zq95_10m.tif
#     secondary/<parcel_id>_cov_gt2_10m.tif
#     secondary/<parcel_id>_cov_gt5_10m.tif
#     manifest.csv         (current state: one row per parcel; reruns replace
#                           the parcel's row in place rather than duplicate)
#     manifest_history.csv (append-only history of every successful per-parcel
#                           run, for traceability across reruns)
#
# Per-parcel pipeline:
#   1. Build a clean 10 m UTM template `ref` covering the AP footprint + 50 m.
#   2. Point-cloud-side metrics in one pass:
#        pc = pixel_metrics(ctg_ap, ~pc_metrics(Z), res = ref)
#      Producing zq50, zq95, zmax_minus_zq95, cov_gt2, cov_gt5 in one chunk
#      pass. zq95 is then used both as the teacher's PC-side estimate and as
#      part of the secondary archive.
#   3. CHM-side fine raster:
#        chm_1m = rasterize_canopy(ctg_ap, res = disagg(ref, 10),
#                                  algorithm = p2r(subcircle = 0.2))
#      p2r(subcircle = 0.2) is used (not pitfree) because DJI Terra's
#      non-standard scale factors break lidR's Delaunay path -- see notes/009
#      and notes/010.
#   4. CHM-side aggregates to 10 m:
#        chm_p95_partial   = aggregate(chm_1m, fact = 10, fun = q95)
#        chm_p95_count     = aggregate(!is.na(chm_1m), fact = 10, fun = sum)
#        chm_p95           = chm_p95_partial where chm_p95_count >= threshold,
#                            NA otherwise   (see min-support rule above)
#        chm_mean          = aggregate(chm_1m, fact = 10, fun = mean)
#        chm_iqr           = aggregate(chm_1m, fact = 10, fun = IQR)
#   5. Apply the agreement rule pixel-wise; emit chm_p95, zq95, agreement
#      abs_diff, agreement_qc, agreement_reason, and the teacher-labels
#      raster (mean where qc=1, NA otherwise).
#   6. Write the secondary feature library alongside.
#
# ============================================================================

library(lidR)
library(terra)
library(sf)

source("R/helpers.R")

# polyfill for older R
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
}

# --- config ----------------------------------------------------------------
# Edit these two lines for each protected area you want to process, then
# Source-on-Save in RStudio (or `source("R/04_batch_teacher_grid_v1.R")`).
# The script is intentionally not CLI-driven; same convention as the other
# 02_/03_ batch scripts in this repo.
lidar_root <- "D:/lidar_rbmn"   # parent dir holding one folder per parcel
area_code  <- "rbmn"            # short code; outputs land under out/teacher/<area_code>/
overwrite  <- FALSE             # set TRUE to recompute parcels that already have output

# Minimum number of valid 1 m sub-cells (out of 100) required for a 10 m
# `chm_p95` to be written. Cells below the threshold are set to NA so that
# the agreement rule downstream sees them as `missing_chm` rather than as
# spuriously valid q95s computed from a sparse subset.
#
# Architect decision: 80 is an *implementation default*, not a frozen
# scientific constant. The frozen contract does not pin a numeric threshold;
# 80 is a defensible middle ground (100 is too brittle for edge cells; 50 is
# too permissive for honest CHM-side coverage handling). The chosen value is
# echoed into every manifest row as `min_valid_subcells_threshold` so the
# run state is never hidden.
min_valid_subcells <- 80L

# Agreement rule constants. v1.1 empirical re-cut, architect-approved
# 2026-05-03 (notes/012). v1 values were (0.5, 0.10); v1.1 values are
# (1.5, 0.15). The structure `tol = max(floor, ratio * mean_h)` is
# preserved. Bumping these requires bumping `agreement_rule_version`
# below and updating notes/011 §3.6 / §5.4 and the spec docs in
# `instructions/`.
agreement_floor_m       <- 1.5
agreement_ratio         <- 0.15
agreement_rule_version  <- "v1.1"

if (!nzchar(area_code)) {
  stop("area_code is required (set it in the config block at the top).")
}

# Output layout
out_root           <- file.path("out", "teacher", area_code)
dir_chm_p95        <- file.path(out_root, "chm_p95")
dir_chm_p95_count  <- file.path(out_root, "chm_p95_count")
dir_zq95           <- file.path(out_root, "zq95")
dir_agreement      <- file.path(out_root, "agreement")
dir_agreement_qc   <- file.path(out_root, "agreement_qc")
dir_agreement_rsn  <- file.path(out_root, "agreement_reason")
dir_teacher        <- file.path(out_root, "teacher_labels")
dir_secondary      <- file.path(out_root, "secondary")
manifest_path         <- file.path(out_root, "manifest.csv")
manifest_history_path <- file.path(out_root, "manifest_history.csv")

for (d in c(dir_chm_p95, dir_chm_p95_count, dir_zq95, dir_agreement,
            dir_agreement_qc, dir_agreement_rsn, dir_teacher, dir_secondary)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- discovery -------------------------------------------------------------
triplets <- find_lidar_triplets(lidar_root)
cat("Area:", area_code, "\n")
cat("LiDAR root:", lidar_root, "\n")
cat("Out root:", out_root, "\n")
cat("Found", length(triplets), "parcel(s)\n\n")

# --- pixel_metrics helper --------------------------------------------------
# lidR's formula parser only accepts a single function call on the RHS, and
# dispatches it into a non-global evaluation context. Defining `pc_metrics`
# at the script's top level (rather than inside the per-parcel worker) lets
# the formula `~pc_metrics(Z)` resolve via the global environment on every
# worker. Helpers used only by `terra::aggregate` (q95_fun / iqr_fun) stay
# inside the worker because terra evaluates them in the calling frame.
pc_metrics <- function(z) {
  qs <- stats::quantile(z, c(0.5, 0.95), na.rm = TRUE)
  list(
    zq50            = qs[[1]],
    zq95            = qs[[2]],
    zmax_minus_zq95 = max(z, na.rm = TRUE) - qs[[2]]
  )
}

# --- raster IO helpers -----------------------------------------------------
gdal_float <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")
gdal_byte  <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")

write_float <- function(r, path) {
  terra::writeRaster(
    r, path, overwrite = TRUE,
    datatype = "FLT4S", gdal = gdal_float, NAflag = -9999
  )
}

# Byte writer for the categorical QA layers. NoData = 255 so codes 0..5 are
# all encodable.
write_byte <- function(r, path) {
  terra::writeRaster(
    r, path, overwrite = TRUE,
    datatype = "INT1U", gdal = gdal_byte, NAflag = 255
  )
}

# --- per-parcel worker -----------------------------------------------------
process_one <- function(trip) {
  parcel_id <- trip$area_id
  out_label <- file.path(dir_teacher,
                         sprintf("%s_H_dom_consensus_10m.tif", parcel_id))

  if (file.exists(out_label) && !overwrite) {
    cat("[skip]", parcel_id, "- already exists\n")
    return(invisible(NULL))
  }

  cat("[proc]", parcel_id, "... ")
  t0 <- Sys.time()

  # ----- catalog + filter (kill any sub-ground noise on AP) ---------------
  ctg_ap <- readLAScatalog(trip$ap)
  opt_filter(ctg_ap) <- "-drop_z_below 0"

  # ----- 10 m UTM reference template covering the AP footprint + 50 m ----
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

  # ----- point-cloud-side metrics (one chunk pass) ------------------------
  # Computes the teacher-side zq95 alongside the secondary archive features
  # zq50 and zmax_minus_zq95. Single pass over the catalog to avoid
  # re-reading points. Cover features are CHM-derived (see below) per
  # architect decision and are not computed here.
  # `pc_metrics` is defined at the script top level so lidR's formula
  # dispatch can resolve it via the global environment.
  pc_r <- pixel_metrics(ctg_ap, ~pc_metrics(Z), res = ref)

  zq95_r            <- pc_r[["zq95"]];            names(zq95_r)            <- "zq95"
  zq50_r            <- pc_r[["zq50"]];            names(zq50_r)            <- "zq50"
  zmax_minus_zq95_r <- pc_r[["zmax_minus_zq95"]]; names(zmax_minus_zq95_r) <- "zmax_minus_zq95"

  # ----- 1 m CHM via p2r --------------------------------------------------
  fine_ref <- terra::disagg(ref, fact = 10)
  chm_1m   <- rasterize_canopy(
    ctg_ap,
    res       = fine_ref,
    algorithm = p2r(subcircle = 0.2)
  )

  # ----- CHM-side aggregates to 10 m -------------------------------------
  q95_fun <- function(x, ...) stats::quantile(x, 0.95, na.rm = TRUE)
  iqr_fun <- function(x, ...) {
    qs <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)
    qs[[2]] - qs[[1]]
  }

  # Per-cell count of valid 1 m sub-cells (0..100) -- coverage layer for the
  # chm_p95 min-support rule and an archived QA layer in its own right.
  chm_subcell_count <- terra::aggregate(
    !is.na(chm_1m), fact = 10, fun = sum, na.rm = TRUE
  )
  names(chm_subcell_count) <- "chm_p95_count"

  chm_p95_partial   <- terra::aggregate(chm_1m, fact = 10, fun = q95_fun)
  chm_p95_r         <- terra::ifel(chm_subcell_count >= min_valid_subcells,
                                   chm_p95_partial, NA)
  names(chm_p95_r)  <- "chm_p95"

  chm_mean_r        <- terra::aggregate(chm_1m, fact = 10, fun = mean,
                                        na.rm = TRUE)
  names(chm_mean_r) <- "chm_mean"

  chm_iqr_r         <- terra::aggregate(chm_1m, fact = 10, fun = iqr_fun)
  names(chm_iqr_r)  <- "chm_iqr"

  # CHM-occupancy cover features: fraction of valid 1 m sub-cells whose CHM
  # exceeds the threshold, per 10 m cell. Architect-frozen canonical form
  # for `cov_gt2` / `cov_gt5`. Implemented as `aggregate(mask, fun = mean,
  # na.rm = TRUE)` so the denominator is the count of valid 1 m sub-cells
  # (not the full 100), which matches the "fraction over the supported
  # area" semantics. Cells with zero supported sub-cells become NaN; we
  # convert to NA explicitly for raster cleanliness.
  chm_gt2  <- chm_1m > 2
  chm_gt5  <- chm_1m > 5
  cov_gt2_r <- terra::aggregate(chm_gt2, fact = 10, fun = mean, na.rm = TRUE)
  cov_gt5_r <- terra::aggregate(chm_gt5, fact = 10, fun = mean, na.rm = TRUE)
  cov_gt2_r <- terra::ifel(is.nan(cov_gt2_r), NA, cov_gt2_r)
  cov_gt5_r <- terra::ifel(is.nan(cov_gt5_r), NA, cov_gt5_r)
  names(cov_gt2_r) <- "cov_gt2"
  names(cov_gt5_r) <- "cov_gt5"

  # Defensive: pc_r and the CHM aggregates share `ref` by construction, but
  # `pixel_metrics` can nudge extent/origin during chunking. Snap any drift
  # back to the chm_p95 grid.
  align_to_chm <- function(r) {
    if (!terra::compareGeom(r, chm_p95_r, stopOnError = FALSE)) {
      terra::resample(r, chm_p95_r, method = "near")
    } else {
      r
    }
  }
  zq95_r            <- align_to_chm(zq95_r)
  zq50_r            <- align_to_chm(zq50_r)
  zmax_minus_zq95_r <- align_to_chm(zmax_minus_zq95_r)

  # ----- agreement-rule masks --------------------------------------------
  chm <- chm_p95_r
  zq  <- zq95_r

  miss_chm  <- is.na(chm) & !is.na(zq)
  miss_zq   <- !is.na(chm) & is.na(zq)
  miss_both <- is.na(chm) & is.na(zq)
  both      <- !is.na(chm) & !is.na(zq)

  mean_h  <- (chm + zq) / 2
  abs_dif <- abs(chm - zq)
  tol     <- terra::ifel(agreement_ratio * mean_h > agreement_floor_m,
                         agreement_ratio * mean_h, agreement_floor_m)

  within  <- both & (abs_dif <= tol)
  outside <- both & (abs_dif >  tol)

  # Replace NAs in the disjoint masks with FALSE so the additive build below
  # produces a clean integer code raster with no spurious NAs.
  zero_na <- function(r) terra::ifel(is.na(r), 0, r)
  within_v   <- zero_na(within)
  outside_v  <- zero_na(outside)
  miss_chm_v <- zero_na(miss_chm)
  miss_zq_v  <- zero_na(miss_zq)
  miss_both_v <- zero_na(miss_both)

  # Reason raster: exactly one of the five masks is TRUE per pixel, so the
  # weighted sum is the reason code 1..5.
  reason <- 1 * within_v + 2 * outside_v + 3 * miss_chm_v +
            4 * miss_zq_v + 5 * miss_both_v
  names(reason) <- "agreement_reason"

  # Binary qc: 1 = pass, 0 = withheld (active disagreement), NA = untestable.
  qc <- terra::ifel(within_v == 1, 1L,
        terra::ifel(outside_v == 1, 0L, NA_integer_))
  names(qc) <- "agreement_qc"

  # Teacher labels: mean of the two estimates, NA where qc != 1.
  teacher <- terra::ifel(within_v == 1, mean_h, NA)
  names(teacher) <- "H_dom_consensus"

  # Agreement abs-diff raster (handy QA/visualization layer).
  abs_dif_out <- terra::ifel(both, abs_dif, NA)
  names(abs_dif_out) <- "agreement_abs_diff"

  # ----- write ------------------------------------------------------------
  # Teacher-side primary outputs.
  write_float(chm_p95_r,  file.path(dir_chm_p95,
                                    sprintf("%s_chm_p95_10m.tif", parcel_id)))
  write_byte (chm_subcell_count,
              file.path(dir_chm_p95_count,
                        sprintf("%s_chm_p95_count_10m.tif", parcel_id)))
  write_float(zq95_r,     file.path(dir_zq95,
                                    sprintf("%s_zq95_10m.tif", parcel_id)))
  write_float(abs_dif_out, file.path(dir_agreement,
                                    sprintf("%s_agreement_abs_diff_10m.tif",
                                            parcel_id)))
  write_byte (qc,         file.path(dir_agreement_qc,
                                    sprintf("%s_agreement_qc_10m.tif",
                                            parcel_id)))
  write_byte (reason,     file.path(dir_agreement_rsn,
                                    sprintf("%s_agreement_reason_10m.tif",
                                            parcel_id)))
  write_float(teacher,    out_label)

  # Archived secondary feature library (note: zq95 lives under `zq95/`, not
  # duplicated here -- the secondary archive is the union of `zq95/` and
  # `secondary/*`, plus `chm_p95/` for QA traceability).
  write_float(chm_mean_r,        file.path(dir_secondary,
                                           sprintf("%s_chm_mean_10m.tif",
                                                   parcel_id)))
  write_float(chm_iqr_r,         file.path(dir_secondary,
                                           sprintf("%s_chm_iqr_10m.tif",
                                                   parcel_id)))
  write_float(zq50_r,            file.path(dir_secondary,
                                           sprintf("%s_zq50_10m.tif",
                                                   parcel_id)))
  write_float(zmax_minus_zq95_r, file.path(dir_secondary,
                                           sprintf("%s_zmax_minus_zq95_10m.tif",
                                                   parcel_id)))
  write_float(cov_gt2_r,         file.path(dir_secondary,
                                           sprintf("%s_cov_gt2_10m.tif",
                                                   parcel_id)))
  write_float(cov_gt5_r,         file.path(dir_secondary,
                                           sprintf("%s_cov_gt5_10m.tif",
                                                   parcel_id)))

  # ----- manifest row -----------------------------------------------------
  fr <- terra::freq(reason)
  get_n <- function(v) {
    s <- sum(fr$count[fr$value == v])
    if (length(s) == 0L || is.na(s)) 0L else as.integer(s)
  }
  n_within    <- get_n(1)
  n_outside   <- get_n(2)
  n_miss_chm  <- get_n(3)
  n_miss_zq   <- get_n(4)
  n_miss_both <- get_n(5)
  n_total     <- n_within + n_outside + n_miss_chm + n_miss_zq + n_miss_both

  pct <- function(n) if (n_total == 0L) NA_real_ else 100 * n / n_total

  # Count of 10 m cells where chm_p95 was set to NA by the min-support rule
  # (i.e. cells where at least one 1 m sub-cell had a CHM value but fewer
  # than `min_valid_subcells` did). This is a strict subset of n_missing_chm
  # plus n_missing_both: cells where chm_1m had zero non-NA sub-cells become
  # n_missing_both (or n_missing_chm if zq is present); cells where chm_1m
  # had some but < threshold also land in those buckets, and this counter
  # quantifies the latter portion.
  n_under_support <- as.integer(terra::global(
    (chm_subcell_count > 0L) & (chm_subcell_count < min_valid_subcells),
    fun = "sum", na.rm = TRUE
  )[1, 1])

  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  row <- data.frame(
    run_timestamp                = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    area_code                    = area_code,
    parcel_id                    = parcel_id,
    agreement_rule_version       = agreement_rule_version,
    agreement_floor_m            = agreement_floor_m,
    agreement_ratio              = agreement_ratio,
    min_valid_subcells_threshold = min_valid_subcells,
    n_cells_total                = n_total,
    n_within                     = n_within,
    n_outside                    = n_outside,
    n_missing_chm                = n_miss_chm,
    n_missing_zq                 = n_miss_zq,
    n_missing_both               = n_miss_both,
    n_chm_p95_undersupported     = n_under_support,
    pct_within                   = round(pct(n_within),    2),
    pct_outside                  = round(pct(n_outside),   2),
    pct_missing                  = round(pct(n_miss_chm + n_miss_zq + n_miss_both), 2),
    runtime_s                    = dt,
    chm_p95_path                 = file.path(dir_chm_p95,
                                  sprintf("%s_chm_p95_10m.tif", parcel_id)),
    chm_p95_count_path           = file.path(dir_chm_p95_count,
                                  sprintf("%s_chm_p95_count_10m.tif",
                                          parcel_id)),
    zq95_path                    = file.path(dir_zq95,
                                  sprintf("%s_zq95_10m.tif", parcel_id)),
    teacher_path                 = out_label,
    secondary_dir                = dir_secondary,
    stringsAsFactors             = FALSE
  )
  # Manifest write strategy:
  #   - `manifest.csv` is the *current state*: one row per parcel. A rerun
  #     of the same parcel replaces that parcel's row in place rather than
  #     appending a duplicate.
  #   - `manifest_history.csv` is *append-only*: every successful per-parcel
  #     run is appended for traceability across reruns.
  upsert_manifest <- function(row, path) {
    if (file.exists(path)) {
      cur <- tryCatch(
        utils::read.csv(path, stringsAsFactors = FALSE,
                        check.names = FALSE),
        error = function(e) NULL
      )
      if (!is.null(cur) && nrow(cur) > 0L &&
          all(c("area_code", "parcel_id") %in% names(cur))) {
        # Drop the existing row for this parcel (if any); union columns to
        # tolerate older manifests written before new fields existed.
        keep <- !(cur$area_code == row$area_code &
                  cur$parcel_id == row$parcel_id)
        cur <- cur[keep, , drop = FALSE]
        all_cols <- union(names(cur), names(row))
        for (cn in setdiff(all_cols, names(cur))) cur[[cn]] <- NA
        for (cn in setdiff(all_cols, names(row))) row[[cn]] <- NA
        cur <- cur[, all_cols, drop = FALSE]
        row <- row[, all_cols, drop = FALSE]
        merged <- rbind(cur, row)
      } else {
        merged <- row
      }
    } else {
      merged <- row
    }
    utils::write.table(merged, path, sep = ",", row.names = FALSE,
                       col.names = TRUE, append = FALSE, quote = TRUE)
  }
  upsert_manifest(row, manifest_path)

  append_history <- file.exists(manifest_history_path)
  utils::write.table(row, manifest_history_path, sep = ",", row.names = FALSE,
                     col.names = !append_history, append = append_history,
                     quote = TRUE)

  cat(sprintf("done in %.1fs (within=%d / outside=%d / missing=%d)\n",
              dt, n_within, n_outside,
              n_miss_chm + n_miss_zq + n_miss_both))
  invisible(out_label)
}

# --- loop ------------------------------------------------------------------
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
cat(sprintf("\nDone. %d / %d parcel(s) produced teacher labels in %s\n",
            ok, length(triplets), dir_teacher))

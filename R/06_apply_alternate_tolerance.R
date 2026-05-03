# ============================================================================
# 06_apply_alternate_tolerance.R
#
# Re-derive the rule-dependent teacher outputs from already-archived
# rule-independent rasters, under a chosen `(floor_m, ratio)` tolerance
# pair. Designed for the v1 -> v1.1 migration approved by the architect on
# 2026-05-03 (notes/012), but works for any future re-cut: the durable
# artifacts are `chm_p95/`, `zq95/`, and `chm_p95_count/`; the rule-dependent
# artifacts (`agreement/`, `agreement_qc/`, `agreement_reason/`,
# `teacher_labels/`) are the ones this script regenerates.
#
# This is the cheap path: no LiDAR re-processing, no `pixel_metrics` pass,
# no `rasterize_canopy` call. Reads existing GeoTIFFs and writes new
# GeoTIFFs in seconds per parcel.
#
# RStudio-friendly: edit the config block below and Source.
# ============================================================================

library(terra)

# --- config ----------------------------------------------------------------
area_code              <- "rbmn"
agreement_floor_m      <- 1.5     # was 0.5 in v1; v1.1 architect-approved
agreement_ratio        <- 0.15    # was 0.10 in v1
agreement_rule_version <- "v1.1"
overwrite              <- TRUE    # almost always TRUE for a re-derivation

# Where the producer wrote the rule-independent artifacts and where to
# write the new rule-dependent ones. Same layout as
# `R/04_batch_teacher_grid_v1.R`.
out_root           <- file.path("out", "teacher", area_code)
dir_chm_p95        <- file.path(out_root, "chm_p95")
dir_chm_p95_count  <- file.path(out_root, "chm_p95_count")
dir_zq95           <- file.path(out_root, "zq95")
dir_agreement      <- file.path(out_root, "agreement")
dir_agreement_qc   <- file.path(out_root, "agreement_qc")
dir_agreement_rsn  <- file.path(out_root, "agreement_reason")
dir_teacher        <- file.path(out_root, "teacher_labels")
manifest_path      <- file.path(out_root, "manifest.csv")
manifest_history   <- file.path(out_root, "manifest_history.csv")

if (!dir.exists(out_root)) {
  stop("Area output dir not found: ", out_root,
       "  -- run R/04_batch_teacher_grid_v1.R first to populate the ",
       "rule-independent artifacts.")
}
for (d in c(dir_agreement, dir_agreement_qc, dir_agreement_rsn,
            dir_teacher)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- raster IO helpers (match producer script's conventions) ---------------
gdal_float <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")
gdal_byte  <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")
write_float <- function(r, path) {
  terra::writeRaster(r, path, overwrite = TRUE,
                     datatype = "FLT4S", gdal = gdal_float, NAflag = -9999)
}
write_byte <- function(r, path) {
  terra::writeRaster(r, path, overwrite = TRUE,
                     datatype = "INT1U", gdal = gdal_byte, NAflag = 255)
}

# --- discover parcels ------------------------------------------------------
chm_files  <- list.files(dir_chm_p95, pattern = "_chm_p95_10m\\.tif$",
                         full.names = FALSE)
parcel_ids <- sub("_chm_p95_10m\\.tif$", "", chm_files)
cat(sprintf("Re-deriving %d parcel(s) under %s\n", length(parcel_ids),
            out_root))
cat(sprintf("Rule version: %s   tol = max(%.3f, %.3f * mean_h)\n\n",
            agreement_rule_version, agreement_floor_m, agreement_ratio))

# --- per-parcel re-derivation ---------------------------------------------
process_one <- function(parcel_id) {
  out_label <- file.path(dir_teacher,
                         sprintf("%s_H_dom_consensus_10m.tif", parcel_id))
  if (file.exists(out_label) && !overwrite) {
    cat("[skip]", parcel_id, "- already exists\n")
    return(invisible(NULL))
  }

  cat("[proc]", parcel_id, "... ")
  t0 <- Sys.time()

  chm_path <- file.path(dir_chm_p95,
                        sprintf("%s_chm_p95_10m.tif", parcel_id))
  zq_path  <- file.path(dir_zq95,
                        sprintf("%s_zq95_10m.tif", parcel_id))
  if (!file.exists(chm_path) || !file.exists(zq_path)) {
    cat("MISSING SOURCE\n")
    return(invisible(NULL))
  }
  chm_p95_r <- terra::rast(chm_path)
  zq95_r    <- terra::rast(zq_path)

  # Same per-pixel logic as producer; only the constants differ.
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

  zero_na <- function(r) terra::ifel(is.na(r), 0, r)
  within_v   <- zero_na(within)
  outside_v  <- zero_na(outside)
  miss_chm_v <- zero_na(miss_chm)
  miss_zq_v  <- zero_na(miss_zq)
  miss_both_v <- zero_na(miss_both)

  reason <- 1 * within_v + 2 * outside_v + 3 * miss_chm_v +
            4 * miss_zq_v + 5 * miss_both_v
  names(reason) <- "agreement_reason"

  qc <- terra::ifel(within_v == 1, 1L,
        terra::ifel(outside_v == 1, 0L, NA_integer_))
  names(qc) <- "agreement_qc"

  teacher <- terra::ifel(within_v == 1, mean_h, NA)
  names(teacher) <- "H_dom_consensus"

  abs_dif_out <- terra::ifel(both, abs_dif, NA)
  names(abs_dif_out) <- "agreement_abs_diff"

  write_float(abs_dif_out, file.path(dir_agreement,
                                     sprintf("%s_agreement_abs_diff_10m.tif",
                                             parcel_id)))
  write_byte (qc,          file.path(dir_agreement_qc,
                                     sprintf("%s_agreement_qc_10m.tif",
                                             parcel_id)))
  write_byte (reason,      file.path(dir_agreement_rsn,
                                     sprintf("%s_agreement_reason_10m.tif",
                                             parcel_id)))
  write_float(teacher,     out_label)

  # ----- per-parcel manifest row -----------------------------------------
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

  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  row <- data.frame(
    run_timestamp                = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    area_code                    = area_code,
    parcel_id                    = parcel_id,
    agreement_rule_version       = agreement_rule_version,
    agreement_floor_m            = agreement_floor_m,
    agreement_ratio              = agreement_ratio,
    n_cells_total                = n_total,
    n_within                     = n_within,
    n_outside                    = n_outside,
    n_missing_chm                = n_miss_chm,
    n_missing_zq                 = n_miss_zq,
    n_missing_both               = n_miss_both,
    pct_within                   = round(pct(n_within),  2),
    pct_outside                  = round(pct(n_outside), 2),
    pct_missing                  = round(pct(n_miss_chm + n_miss_zq +
                                              n_miss_both), 2),
    runtime_s                    = dt,
    teacher_path                 = out_label,
    derivation_source            = "R/06_apply_alternate_tolerance.R",
    stringsAsFactors             = FALSE
  )

  # Upsert into manifest.csv (one row per parcel, current state) -- same
  # discipline as the producer script.
  upsert_manifest <- function(row, path) {
    if (file.exists(path)) {
      cur <- tryCatch(
        utils::read.csv(path, stringsAsFactors = FALSE,
                        check.names = FALSE),
        error = function(e) NULL
      )
      if (!is.null(cur) && nrow(cur) > 0L &&
          all(c("area_code", "parcel_id") %in% names(cur))) {
        keep <- !(cur$area_code == row$area_code &
                  cur$parcel_id == row$parcel_id)
        cur <- cur[keep, , drop = FALSE]
        all_cols <- union(names(cur), names(row))
        for (cn in setdiff(all_cols, names(cur))) cur[[cn]] <- NA
        for (cn in setdiff(all_cols, names(row))) row[[cn]] <- NA
        cur <- cur[, all_cols, drop = FALSE]
        row <- row[, all_cols, drop = FALSE]
        merged <- rbind(cur, row)
      } else { merged <- row }
    } else { merged <- row }
    utils::write.table(merged, path, sep = ",", row.names = FALSE,
                       col.names = TRUE, append = FALSE, quote = TRUE)
  }
  upsert_manifest(row, manifest_path)

  append_history <- file.exists(manifest_history)
  utils::write.table(row, manifest_history, sep = ",", row.names = FALSE,
                     col.names = !append_history, append = append_history,
                     quote = TRUE)

  cat(sprintf("done in %.1fs (within=%d / outside=%d / missing=%d)\n",
              dt, n_within, n_outside,
              n_miss_chm + n_miss_zq + n_miss_both))
  invisible(out_label)
}

# --- loop ------------------------------------------------------------------
results <- vector("list", length(parcel_ids))
for (i in seq_along(parcel_ids)) {
  results[[i]] <- tryCatch(
    process_one(parcel_ids[i]),
    error = function(e) {
      message("[fail] ", parcel_ids[i], ": ", conditionMessage(e))
      NULL
    }
  )
}
ok <- sum(!vapply(results, is.null, logical(1)))
cat(sprintf("\nDone. %d / %d parcel(s) re-derived under rule %s in %s\n",
            ok, length(parcel_ids), agreement_rule_version, dir_teacher))
cat(sprintf("Recommended next step: re-run R/05_audit_teacher_agreement.R ",
            "to confirm the v1.1 rule pooled outside-rate is ~31%%.\n"))

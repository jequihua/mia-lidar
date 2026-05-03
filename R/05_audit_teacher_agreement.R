# ============================================================================
# 05_audit_teacher_agreement.R
#
# Diagnostic audit for the consensus dominant-height teacher product produced
# by `R/04_batch_teacher_grid_v1.R`. Quantifies the per-parcel and
# per-height-bin behavior of the agreement rule so the architect can decide
# whether the observed high `outside` rates reflect:
#
#   (a) the intrinsic statistical offset between
#       chm_p95 = Q_{0.95}( per-1 m maxima )    and
#       zq95    = Q_{0.95}( all returns ),
#       which is positive in expectation and grows with canopy height
#       (note 011 §3.1, §3.4);
#   (b) edge / partial-support effects on chm_p95 that survive the
#       80-of-100 sub-cell rule (note 011 §5.1);
#   (c) genuine source-side data quality issues (notes 011 §3.2, §3.3);
#   (d) some combination of the above.
#
# Hypothesis (a) predicts a *systematic positive offset* (chm − zq > 0) that
# is roughly uniform across the parcel interior and grows with canopy
# height. Hypothesis (b) predicts the offset is concentrated in cells with
# low `chm_subcell_count`. This script computes both signals so they can be
# disentangled empirically.
#
# Per-area inputs (one folder per protected area):
#   out/teacher/<area_code>/chm_p95/<parcel_id>_chm_p95_10m.tif
#   out/teacher/<area_code>/zq95/<parcel_id>_zq95_10m.tif
#   out/teacher/<area_code>/agreement_reason/<parcel_id>_agreement_reason_10m.tif
#   out/teacher/<area_code>/chm_p95_count/<parcel_id>_chm_p95_count_10m.tif
#
# Per-area outputs:
#   out/teacher/<area_code>/audit_agreement_per_parcel.csv
#   out/teacher/<area_code>/audit_agreement_by_height_bin.csv
#   out/teacher/<area_code>/audit_summary.txt   (human-readable headline)
#
# RStudio-friendly: edit `area_code` at the top and Source.
# ============================================================================

library(terra)

# --- config ----------------------------------------------------------------
area_code <- "rbmn"   # which area to audit (matches out/teacher/<area_code>/)

# Height-bin breakpoints for the per-bin breakdown. Aligned with the
# validation bins in feature_stack_decision_note.md §"Added validation
# amendments" so any architect-facing comparison uses comparable strata.
height_breaks <- c(-Inf, 3, 6, 10, 15, 20, Inf)
height_labels <- c("<3", "3-6", "6-10", "10-15", "15-20", ">20")

# Min-support threshold used by the producer script (must match for the
# under-support diagnostic to be meaningful).
min_valid_subcells <- 80L

# Active agreement rule. Must match the producer (R/04) and the
# re-derivation script (R/06). v1 was (0.5, 0.10); v1.1 is (1.5, 0.15)
# per architect approval 2026-05-03 (notes/012). All "outside" stats
# below that re-apply the rule on subsets (full-support, short, tall,
# per-bin) use these constants -- not a hardcoded (0.5, 0.10) -- so the
# audit reflects the *current* rule and not the legacy one.
agreement_floor_m      <- 1.5
agreement_ratio        <- 0.15
agreement_rule_version <- "v1.1"

# Tolerance-sensitivity sweep. The agreement rule is
#   tol = max(floor, ratio * mean_h),    outside <=> |chm - zq| > tol
# These two vectors define the candidate grid. The script re-evaluates the
# rule at every (floor, ratio) combination using the pooled both-present
# (chm - zq, mean_h) vectors across all parcels of the area, producing a
# tradeoff curve so the architect can see the empirical cost of each
# tolerance choice without recomputing any LiDAR aggregate.
sweep_floors <- c(0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
sweep_ratios <- c(0.10, 0.125, 0.15, 0.175, 0.20, 0.25)

# --- locations -------------------------------------------------------------
out_root    <- file.path("out", "teacher", area_code)
dir_chm     <- file.path(out_root, "chm_p95")
dir_zq      <- file.path(out_root, "zq95")
dir_rsn     <- file.path(out_root, "agreement_reason")
dir_cnt     <- file.path(out_root, "chm_p95_count")

per_parcel_path <- file.path(out_root, "audit_agreement_per_parcel.csv")
by_bin_path     <- file.path(out_root, "audit_agreement_by_height_bin.csv")
sweep_path      <- file.path(out_root, "audit_tolerance_sensitivity.csv")
sweep_bin_path  <- file.path(out_root, "audit_tolerance_sensitivity_by_bin.csv")
summary_path    <- file.path(out_root, "audit_summary.txt")

if (!dir.exists(out_root)) {
  stop("Area output dir not found: ", out_root,
       "  -- run R/04_batch_teacher_grid_v1.R first.")
}

# --- discover parcels ------------------------------------------------------
chm_files <- list.files(dir_chm, pattern = "_chm_p95_10m\\.tif$",
                        full.names = FALSE)
parcel_ids <- sub("_chm_p95_10m\\.tif$", "", chm_files)
cat("Auditing", length(parcel_ids), "parcel(s) under", out_root, "\n\n")

# --- helpers ---------------------------------------------------------------
read_vec <- function(path) {
  if (!file.exists(path)) return(NULL)
  as.vector(terra::values(terra::rast(path)))
}

q <- function(v, p) {
  v <- v[!is.na(v)]
  if (length(v) == 0L) NA_real_ else as.numeric(stats::quantile(v, p))
}

safe_mean <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0L) NA_real_ else mean(v)
}

# --- per-parcel + per-bin computation --------------------------------------
audit_one <- function(parcel_id) {
  chm <- read_vec(file.path(dir_chm,
                            sprintf("%s_chm_p95_10m.tif", parcel_id)))
  zq  <- read_vec(file.path(dir_zq,
                            sprintf("%s_zq95_10m.tif", parcel_id)))
  rsn <- read_vec(file.path(dir_rsn,
                            sprintf("%s_agreement_reason_10m.tif", parcel_id)))
  cnt <- read_vec(file.path(dir_cnt,
                            sprintf("%s_chm_p95_count_10m.tif", parcel_id)))
  if (is.null(chm) || is.null(zq) || is.null(rsn)) return(NULL)

  both       <- !is.na(chm) & !is.na(zq)
  diff_all   <- chm - zq
  mean_h_all <- (chm + zq) / 2

  diff_b   <- diff_all[both]
  mean_h_b <- mean_h_all[both]
  cnt_b    <- if (is.null(cnt)) rep(NA_real_, sum(both)) else cnt[both]

  # Reapply the agreement rule on the both-present subset so the audit
  # numbers are independent of the producer manifest. Uses the
  # config-block constants (`agreement_floor_m`, `agreement_ratio`) so
  # the audit reflects the *active* rule. The pct_outside_given_both
  # figure (derived from the rsn raster) and the recomputed outside_b
  # should agree to rounding; divergence would be a producer bug.
  tol_b    <- pmax(agreement_floor_m, agreement_ratio * mean_h_b)
  outside_b <- abs(diff_b) > tol_b

  short <- mean_h_b <= 5
  tall  <- mean_h_b >  5

  n_total     <- length(chm)
  n_both      <- sum(both)
  n_within    <- sum(rsn == 1L, na.rm = TRUE)
  n_outside   <- sum(rsn == 2L, na.rm = TRUE)
  n_miss_chm  <- sum(rsn == 3L, na.rm = TRUE)
  n_miss_zq   <- sum(rsn == 4L, na.rm = TRUE)
  n_miss_both <- sum(rsn == 5L, na.rm = TRUE)

  per_parcel <- data.frame(
    area_code              = area_code,
    parcel_id              = parcel_id,
    n_total                = n_total,
    n_both                 = n_both,
    n_within               = n_within,
    n_outside              = n_outside,
    n_missing_chm          = n_miss_chm,
    n_missing_zq           = n_miss_zq,
    n_missing_both         = n_miss_both,
    pct_outside_given_both = round(100 * n_outside / max(n_both, 1L), 2),
    # Canopy-height context (for interpreting the offset)
    mean_h_p05             = round(q(mean_h_b, 0.05), 3),
    mean_h_p50             = round(q(mean_h_b, 0.50), 3),
    mean_h_p95             = round(q(mean_h_b, 0.95), 3),
    # Signed offset diagnostic (key statistic for hypothesis (a))
    mean_offset            = round(safe_mean(diff_b), 3),
    median_offset          = round(q(diff_b, 0.50), 3),
    sd_offset              = round(if (length(diff_b[!is.na(diff_b)]) > 1L)
                                     sd(diff_b, na.rm = TRUE) else NA_real_, 3),
    offset_p05             = round(q(diff_b, 0.05), 3),
    offset_p25             = round(q(diff_b, 0.25), 3),
    offset_p75             = round(q(diff_b, 0.75), 3),
    offset_p95             = round(q(diff_b, 0.95), 3),
    pct_offset_positive    = round(100 * safe_mean(diff_b > 0), 2),
    # Stratified by tolerance regime (floor vs. ratio)
    n_short_h_le_5         = sum(short),
    n_tall_h_gt_5          = sum(tall),
    pct_outside_short      = round(100 * safe_mean(outside_b[short]), 2),
    pct_outside_tall       = round(100 * safe_mean(outside_b[tall]),  2),
    mean_offset_short      = round(safe_mean(diff_b[short]), 3),
    mean_offset_tall       = round(safe_mean(diff_b[tall]),  3),
    # Edge / partial-support diagnostic (hypothesis (b))
    mean_subcell_count     = round(safe_mean(cnt_b), 1),
    pct_full_support       = round(100 * safe_mean(cnt_b == 100), 2),
    pct_undersupported     = round(100 * safe_mean(
                                     !is.na(cnt) & cnt > 0 &
                                     cnt < min_valid_subcells), 2),
    # Outside fraction restricted to fully-supported cells: if (b) dominates,
    # this rate should drop sharply; if (a) dominates, it should be similar
    # to the overall rate.
    pct_outside_full_support = {
      idx_full <- both & !is.na(cnt) & cnt == 100
      if (sum(idx_full) == 0L) NA_real_
      else round(100 *
                 safe_mean(abs(diff_all[idx_full]) >
                           pmax(0.5, 0.10 * mean_h_all[idx_full])), 2)
    },
    stringsAsFactors = FALSE
  )

  # Per-height-bin breakdown for the offset-vs-height plot
  bin <- cut(mean_h_b, breaks = height_breaks, labels = height_labels,
             right = FALSE, include.lowest = TRUE)
  by_bin_rows <- lapply(seq_along(height_labels), function(i) {
    sel <- !is.na(bin) & bin == height_labels[i]
    n   <- sum(sel)
    data.frame(
      area_code     = area_code,
      parcel_id     = parcel_id,
      height_bin    = height_labels[i],
      n_cells       = n,
      n_outside     = sum(outside_b[sel], na.rm = TRUE),
      pct_outside   = if (n == 0L) NA_real_
                      else round(100 * safe_mean(outside_b[sel]), 2),
      mean_offset   = round(safe_mean(diff_b[sel]),    3),
      median_offset = round(q(diff_b[sel], 0.50),       3),
      mean_h        = round(safe_mean(mean_h_b[sel]),  3),
      stringsAsFactors = FALSE
    )
  })
  by_bin <- do.call(rbind, by_bin_rows)

  list(per_parcel = per_parcel, by_bin = by_bin,
       diff_b = diff_b, mean_h_b = mean_h_b)
}

# --- run -------------------------------------------------------------------
results <- vector("list", length(parcel_ids))
for (i in seq_along(parcel_ids)) {
  pid <- parcel_ids[i]
  cat(sprintf("[%3d/%d] %s ...\n", i, length(parcel_ids), pid))
  results[[i]] <- tryCatch(audit_one(pid),
                           error = function(e) {
                             message("  fail: ", conditionMessage(e)); NULL
                           })
}
results <- Filter(Negate(is.null), results)

per_parcel_df <- do.call(rbind, lapply(results, `[[`, "per_parcel"))
by_bin_df     <- do.call(rbind, lapply(results, `[[`, "by_bin"))

# --- pooled tolerance-sensitivity sweep -----------------------------------
# Pool the both-present (chm - zq, mean_h) vectors across all parcels of
# this area, then re-evaluate the agreement rule under every (floor, ratio)
# candidate. Cheap (a few MB of doubles) and answers the architect-facing
# question directly: "what would the global outside-rate look like if we
# moved the threshold here?"
diff_pool   <- unlist(lapply(results, `[[`, "diff_b"))
mean_h_pool <- unlist(lapply(results, `[[`, "mean_h_b"))
abs_pool    <- abs(diff_pool)
n_pool      <- length(diff_pool)

# Pre-compute height-bin labels for the pooled vector so the per-bin sweep
# is one tabulate per candidate.
bin_pool <- cut(mean_h_pool, breaks = height_breaks, labels = height_labels,
                right = FALSE, include.lowest = TRUE)

sweep_rows <- list()
sweep_bin_rows <- list()
for (fl in sweep_floors) {
  for (rt in sweep_ratios) {
    tol_pool  <- pmax(fl, rt * mean_h_pool)
    out_pool  <- abs_pool > tol_pool
    n_out     <- sum(out_pool, na.rm = TRUE)

    sweep_rows[[length(sweep_rows) + 1L]] <- data.frame(
      area_code        = area_code,
      floor_m          = fl,
      ratio            = rt,
      n_both_present   = n_pool,
      n_outside        = n_out,
      pct_outside_pool = round(100 * n_out / max(n_pool, 1L), 2),
      stringsAsFactors = FALSE
    )

    for (lbl in height_labels) {
      sel <- !is.na(bin_pool) & bin_pool == lbl
      n_b <- sum(sel)
      n_o <- sum(out_pool[sel], na.rm = TRUE)
      sweep_bin_rows[[length(sweep_bin_rows) + 1L]] <- data.frame(
        area_code      = area_code,
        floor_m        = fl,
        ratio          = rt,
        height_bin     = lbl,
        n_cells        = n_b,
        n_outside      = n_o,
        pct_outside    = if (n_b == 0L) NA_real_
                         else round(100 * n_o / n_b, 2),
        stringsAsFactors = FALSE
      )
    }
  }
}
sweep_df     <- do.call(rbind, sweep_rows)
sweep_bin_df <- do.call(rbind, sweep_bin_rows)

# Compact wide view for the console (rows = floor, cols = ratio)
sweep_wide <- tapply(sweep_df$pct_outside_pool,
                     list(sweep_df$floor_m, sweep_df$ratio),
                     identity)
dimnames(sweep_wide) <- list(
  floor_m = sprintf("%.2f", sweep_floors),
  ratio   = sprintf("%.3f", sweep_ratios)
)

# --- write -----------------------------------------------------------------
utils::write.csv(per_parcel_df, per_parcel_path, row.names = FALSE)
utils::write.csv(by_bin_df,     by_bin_path,     row.names = FALSE)
utils::write.csv(sweep_df,      sweep_path,      row.names = FALSE)
utils::write.csv(sweep_bin_df,  sweep_bin_path,  row.names = FALSE)

# --- headline summary (printed + saved) ------------------------------------
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
med <- function(x) stats::median(x, na.rm = TRUE)

# Aggregate the per-bin table across parcels for the global view
bin_pooled <- do.call(rbind, lapply(height_labels, function(lbl) {
  rows <- by_bin_df[by_bin_df$height_bin == lbl, , drop = FALSE]
  if (nrow(rows) == 0L) return(NULL)
  total_cells   <- sum(rows$n_cells,   na.rm = TRUE)
  total_outside <- sum(rows$n_outside, na.rm = TRUE)
  data.frame(
    height_bin       = lbl,
    n_cells_pooled   = total_cells,
    pct_outside_pool = if (total_cells == 0L) NA_real_
                       else round(100 * total_outside / total_cells, 2),
    median_offset    = round(med(rows$median_offset), 3),
    median_pct_outside_per_parcel = round(med(rows$pct_outside), 2),
    median_mean_h    = round(med(rows$mean_h), 3),
    stringsAsFactors = FALSE
  )
}))

lines <- c(
  "============================================================",
  sprintf("Teacher-agreement audit  -  area = %s", area_code),
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  sprintf("Parcels audited: %d", nrow(per_parcel_df)),
  sprintf("Active rule: %s   tol = max(%.3f, %.3f * mean_h)",
          agreement_rule_version, agreement_floor_m, agreement_ratio),
  "============================================================",
  "",
  "Per-parcel headline statistics (median across parcels):",
  sprintf("  pct_outside_given_both           : %s %%", fmt(med(per_parcel_df$pct_outside_given_both))),
  sprintf("  pct_outside_full_support         : %s %%", fmt(med(per_parcel_df$pct_outside_full_support))),
  sprintf("  pct_outside_short  (mean_h <=5m) : %s %%", fmt(med(per_parcel_df$pct_outside_short))),
  sprintf("  pct_outside_tall   (mean_h > 5m) : %s %%", fmt(med(per_parcel_df$pct_outside_tall))),
  sprintf("  mean_offset (chm-zq)             : %s m", fmt(med(per_parcel_df$mean_offset), 3)),
  sprintf("  median_offset (chm-zq)           : %s m", fmt(med(per_parcel_df$median_offset), 3)),
  sprintf("  pct_offset_positive              : %s %%", fmt(med(per_parcel_df$pct_offset_positive))),
  sprintf("  mean_offset_short                : %s m", fmt(med(per_parcel_df$mean_offset_short), 3)),
  sprintf("  mean_offset_tall                 : %s m", fmt(med(per_parcel_df$mean_offset_tall),  3)),
  sprintf("  mean_subcell_count (out of 100)  : %s",   fmt(med(per_parcel_df$mean_subcell_count), 1)),
  sprintf("  pct_undersupported               : %s %%", fmt(med(per_parcel_df$pct_undersupported))),
  "",
  "Per-height-bin pooled across all parcels in this area:",
  paste0(capture.output(print(bin_pooled, row.names = FALSE)), collapse = "\n"),
  "",
  "Interpretation guide:",
  "  Hypothesis (a) -- intrinsic Q0.95(maxima) vs Q0.95(returns) offset:",
  "    EXPECTED: median_offset > 0, pct_offset_positive > 80, mean_offset_tall",
  "    notably > mean_offset_short, pct_outside_tall > pct_outside_short.",
  "    pct_outside_full_support remains comparable to overall pct_outside.",
  "  Hypothesis (b) -- edge / partial-support effects:",
  "    EXPECTED: pct_outside_full_support << pct_outside_given_both,",
  "    pct_undersupported is large, pct_full_support is small.",
  "  Hypothesis (c) -- source-side data quality issues:",
  "    EXPECTED: parcel-to-parcel variability dominates, no consistent sign",
  "    or height-dependence of the offset.",
  "",
  "Tolerance-sensitivity sweep -- pooled outside-rate (%) across all parcels",
  "in this area, for each (floor_m, ratio) candidate. Rows = floor (m),",
  "columns = ratio (fraction of mean canopy height).",
  paste0(capture.output(print(round(sweep_wide, 2))), collapse = "\n"),
  "",
  "Tolerance-sensitivity sweep, restricted to the 3-10 m canopy band where",
  "the current rule fails worst (pooled across both 3-6 and 6-10 bins):",
  paste0(capture.output(print({
    band <- sweep_bin_df[sweep_bin_df$height_bin %in% c("3-6", "6-10"), ]
    band_agg <- aggregate(cbind(n_cells, n_outside) ~ floor_m + ratio,
                          data = band, FUN = sum)
    band_agg$pct_outside_3_10m <- round(
      100 * band_agg$n_outside / pmax(band_agg$n_cells, 1L), 2)
    band_wide <- tapply(band_agg$pct_outside_3_10m,
                        list(band_agg$floor_m, band_agg$ratio), identity)
    dimnames(band_wide) <- list(
      floor_m = sprintf("%.2f", sweep_floors),
      ratio   = sprintf("%.3f", sweep_ratios)
    )
    round(band_wide, 2)
  })), collapse = "\n"),
  "",
  "Reading: smaller numbers = more cells pass the agreement rule.",
  sprintf("The active rule is %s: floor_m=%.2f, ratio=%.3f. (v1 was 0.50, 0.100;",
          agreement_rule_version, agreement_floor_m, agreement_ratio),
  "v1.1 architect-approved 2026-05-03, see notes/012). A future re-cut",
  "should keep the structure (absolute floor + relative ceiling) and pick a",
  "cell that meaningfully reduces outside-rate without making the rule so",
  "loose it no longer falsifies anything (>= 20% outside is probably still",
  "meaningful; < 5% means the rule passes nearly everything).",
  "",
  sprintf("Outputs:"),
  sprintf("  per-parcel        : %s", per_parcel_path),
  sprintf("  by-height-bin     : %s", by_bin_path),
  sprintf("  tolerance sweep   : %s", sweep_path),
  sprintf("  sweep by bin      : %s", sweep_bin_path),
  sprintf("  this summary      : %s", summary_path)
)

writeLines(lines, summary_path)
cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")

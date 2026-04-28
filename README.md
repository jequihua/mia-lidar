# R/ — LiDAR processing scripts

R scripts that turn the raw LiDAR deliveries (one TP/GP/AP triplet per area)
into the raster products consumed by the downstream KGML pipeline. Run from
the project root so relative paths resolve.

## Conventions

- **Input layout.** Each LiDAR area is a folder containing three files named
  `<area_id>_TP.las`, `<area_id>_GP.las`, `<area_id>_AP.las`. See
  `notes/003` for the role of each.
- **Reserves.** `rbmn` (Marismas Nacionales) and `rbrl` (Ría Lagartos), each
  with its own archive root on the external drive.
- **Output layout.** Raster products land under
  `out/lidar/lidar/<reserve>/<product>/<area_id>_<product>.tif`,
  in the LAS native CRS (typically UTM). A short product note lives at
  `out/lidar/lidar/<NNN>_<product>_data_product.md`.
- **Vendor quirks.** DJI Zenmuse L1 / DJI Terra files are pre-normalised
  (Z = height above local ground), all classified as ASPRS code 0, and use
  non-standard scale factors that break lidR's Delaunay path. The scripts
  work around all of this. See `notes/009`.

## Scripts

### `helpers.R`
Shared helpers:
- `find_lidar_triplet(area_dir)` / `find_lidar_triplets(parent_dir)` —
  locate the TP/GP/AP files in one or many area folders.
- `save_lidar_footprint(ctg, target_crs, out_path, ...)` — export a LAS
  catalog footprint as a georeferenced vector for QA in QGIS.

Sourced by every other script. Not run on its own.

### `00_build_lax_indexes.R` — one-shot indexer
Walks a LiDAR archive and builds an `.lax` sidecar next to every `.las` /
`.laz` file that doesn't already have one. `.lax` files are tiny spatial
indexes that turn buffered/spatially-bounded reads from full-file scans
into targeted byte-range reads — every later script reads them transparently.

**Usage.** Edit `lidar_root` at the top, then:

```bash
Rscript R/00_build_lax_indexes.R
```

Run once per reserve. Idempotent — safe to re-run.

### `01_inspect_inputs.R` — interactive walkthrough (one area)
Manual block-by-block exploration of the full feature pipeline on a single
area. Used to sanity-check inputs and prototype before factoring into a
batch script. Produces a multi-band feature stack (terrain derivatives,
CHM aggregates, `.stdmetrics_z`, `.stdmetrics_rn`) projected onto the
Sentinel-2 grid.

**Usage.** Open in RStudio and run blocks with Ctrl/Cmd+Enter. Edit
`area_dir` in block [1] to point at the area you want to inspect. Output
goes to `out/features/<area_id>/`.

### `02_batch_chm_mean_10m.R` — batch CHM-mean 10 m
For every triplet under `lidar_root`, builds a 10 m mean canopy-height
raster from AP via `p2r(subcircle = 0.2)` → 1 m → aggregate(mean) → 10 m.
One single-band GeoTIFF (band: `chm_mean`) per area in the LAS native CRS.

**Usage.** Edit `reserve` and `lidar_root` at the top, then:

```bash
Rscript R/02_batch_chm_mean_10m.R
```

Run once per reserve. Per-area `tryCatch` keeps the batch alive on bad
files; existing outputs are skipped unless `overwrite <- TRUE`.

Product note: `out/lidar/lidar/001_chm_mean_10m_data_product.md`.

### `03_batch_dem_mean_10m.R` — batch DEM-mean 10 m
Counterpart to `02`, but reads GP and interpolates the ground surface via
`knnidw(k = 10, p = 2)` → 1 m → aggregate(mean) → 10 m. One single-band
GeoTIFF (band: `dem_mean`) per area in the LAS native CRS.

The output is **relative** ground micro-topography because the vendor
pre-normalises GP. For absolute elevation an external DEM (SRTM /
national) must be brought in separately.

**Usage.** Edit `reserve` and `lidar_root` at the top, then:

```bash
Rscript R/03_batch_dem_mean_10m.R
```

Same conventions as `02`. Much faster after `00_build_lax_indexes.R` has
run, because `rasterize_terrain` is read-heavy.

Product note: `out/lidar/lidar/002_dem_mean_10m_data_product.md`.

## Typical run order

For a fresh archive on a new external drive:

1. `00_build_lax_indexes.R` — once per reserve. Build all `.lax` sidecars
   up-front.
2. `02_batch_chm_mean_10m.R` — once per reserve. Builds `chm_mean` for
   every area.
3. `03_batch_dem_mean_10m.R` — once per reserve. Builds `dem_mean` for
   every area.
4. `01_inspect_inputs.R` — interactive, ad-hoc, when prototyping new
   features or QA-ing a specific area.

## Dependencies

- R ≥ 4.1 (≥ 4.4 preferred; scripts include a `%||%` polyfill for older).
- R packages: `lidR` (≥ 4.3), `terra` (≥ 1.9), `sf`.

## Related notes

- `notes/001` — lidR package overview, processing engine, speed levers.
- `notes/002` — CHM resolution and the fine-then-aggregate pattern.
- `notes/003` — TP / GP / AP division of labour.
- `notes/004` — KGML pipeline design and LiDAR feature selection.
- `notes/009` — DJI Zenmuse L1 vendor LAS quirks.
- `notes/010` — CHM-mean 10 m product (development-time companion to the
  product note in `out/lidar/lidar/`).

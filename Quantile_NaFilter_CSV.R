# =============================================================================
# 01_compute_quantiles_prefilter.R
#
# ALTERNATIVE VERSION — pre-filter NAs before sampling
#
# Strategy: load ALL values into memory first, remove NoData and NAs,
# then sample from the clean values only.
# This guarantees the sample is drawn exclusively from valid pixels,
# which matters for sparse layers (mangroves, biodiversity hotspots)
# where valid pixels cover only a small fraction of the raster extent.
#
# Trade-off: uses more RAM than the original script because the full
# raster is always loaded into memory, even for large rasters.
#
# Use this script alongside 01_compute_quantiles.R and compare the
# Q25/Q50/Q75 values in the two CSVs to see if sparse layers differ.
# =============================================================================

library(terra)

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
input_dir  <- "/Volumes/Andrea_GIS/Hotspot Mapping/R_test/Quantile"
output_csv <- "/Volumes/Andrea_GIS/Hotspot Mapping/R_test/Quantile/Test_prefilter.csv"
n_classes        <- 4        # number of quantile classes
sample_threshold <- 5e6      # rasters larger than this are sampled (after cleaning)
sample_size      <- 500000   # number of CLEAN values to sample
# ──────────────────────────────────────────────────────────────────────────────


tiff_files <- list.files(input_dir,
                         pattern    = "\\.tif{1,2}$",
                         full.names = TRUE,
                         ignore.case = TRUE)

if (length(tiff_files) == 0) stop("No TIFF files found in: ", input_dir)

message("Found ", length(tiff_files), " TIFF file(s) to process.\n")

results <- lapply(tiff_files, function(fp) {
  
  layer_name <- tools::file_path_sans_ext(basename(fp))
  message("  Processing: ", layer_name)
  
  r      <- terra::rast(fp)
  naflag <- terra::NAflag(r)
  
  # ── Load ALL values and clean BEFORE deciding whether to sample ────────────
  # This is the key difference from the original script:
  # we always load everything first so NA pixels are excluded
  # from the pool before any sampling takes place.
  message("    Loading all values for pre-filtering...")
  all_values <- terra::values(r, mat = FALSE)
  
  # Remove GIS-level NoData sentinel (e.g. -9999)
  if (!is.na(naflag)) all_values <- all_values[all_values != naflag]
  
  # Remove R-level NAs
  all_values <- all_values[!is.na(all_values)]
  
  # ── Cap filter: layers whose name contains ASIS or PEy ──────────────
  cap_filter <- grepl("ASIS|PEy", layer_name, ignore.case = FALSE)
  if (cap_filter) {
    all_values <- all_values[all_values <= 100]
    message("    [cap filter] '", layer_name,
            "' — keeping only values <= 100")
  }
  
  n_clean <- length(all_values)
  message("    Clean (non-NA) pixels: ", format(n_clean, big.mark = ","),
          " out of ", format(terra::ncell(r), big.mark = ","), " total")
  
  if (n_clean == 0) {
    warning("  [SKIP] ", layer_name, " — no non-NA values found.")
    return(NULL)
  }
  
  # ── Now sample from CLEAN values only (if the layer is large) ──────────────
  set.seed(1)
  if (n_clean > sample_threshold) {
    message("    (large clean pool: sampling ", format(sample_size, big.mark = ","),
            " values from ", format(n_clean, big.mark = ","), " clean pixels)")
    values <- all_values[sample(n_clean, size = sample_size, replace = FALSE)]
  } else {
    # Small enough — use all clean values directly, no sampling needed
    message("    (small clean pool: using all ", format(n_clean, big.mark = ","),
            " clean values)")
    values <- all_values
  }
  
  unique_vals <- unique(values)
  
  # ── Edge case: only one unique value → assign highest class ──────────────
  if (length(unique_vals) == 1) {
    message("  [NOTE] ", layer_name,
            " has a single unique value (", unique_vals,
            "). Assigning highest class (", n_classes, ").")
    
    return(data.frame(
      layer        = layer_name,
      file_path    = fp,
      n_classes    = n_classes,
      n_clean_pixels = n_clean,
      note         = "single_value",
      Q25          = unique_vals,
      Q50          = unique_vals,
      Q75          = unique_vals
    ))
  }
  
  # ── Normal case: compute the three quantile thresholds ───────────────────
  # Classes after reclassification:
  #   Class 1 (Low)      : value <  Q25
  #   Class 2 (Med-Low)  : Q25 <= value <  Q50
  #   Class 3 (Med-High) : Q50 <= value <  Q75
  #   Class 4 (High)     : value >= Q75
  breaks <- as.numeric(quantile(values, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
  
  data.frame(
    layer          = layer_name,
    file_path      = fp,
    n_classes      = n_classes,
#   n_clean_pixels = n_clean,       # extra column — useful for comparison
    note           = ifelse(cap_filter, "capped_at_100", "ok"),
    Q25            = breaks[1],
    Q50            = breaks[2],
    Q75            = breaks[3]
  )
})

# Drop any NULLs (skipped layers)
results <- Filter(Negate(is.null), results)

if (length(results) == 0) stop("No layers could be processed.")

breaks_df        <- do.call(rbind, results)
rownames(breaks_df) <- NULL

# Always overwrite the CSV on each run
write.csv(breaks_df, output_csv, row.names = FALSE, quote = TRUE)

message("\n✓ CSV written: ", normalizePath(output_csv, mustWork = FALSE))
message("  ", nrow(breaks_df), " layer(s) written to CSV.")
print(breaks_df[, c("layer", "note", "n_clean_pixels", "Q25", "Q50", "Q75")])
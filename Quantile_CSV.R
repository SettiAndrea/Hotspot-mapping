
# =============================================================================
# 01_compute_quantiles.R
# Computes 4-class quantile breaks for every TIFF in a folder and writes
# a CSV that the companion script (02_reclassify.R) will consume.
# =============================================================================

library(terra)

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
input_dir  <- "/Volumes/Andrea_GIS/Hotspot Mapping/R_test/Quantile"   # <-- change this
output_csv <- "/Volumes/Andrea_GIS/Hotspot Mapping/R_test/Quantile/Test.csv"         # written to working directory
n_classes        <- 4                       # number of quantile classes
sample_threshold <- 5e6                     # rasters larger than this are sampled
sample_size      <- 500000                  # number of cells to sample
# ──────────────────────────────────────────────────────────────────────────────

tiff_files <- list.files(input_dir,
                         pattern  = "\\.tif{1,2}$",
                         full.names = TRUE,
                         ignore.case = TRUE)

if (length(tiff_files) == 0) stop("No TIFF files found in: ", input_dir)

# ── Skip layers already in the CSV (incremental update) ───────────────────────
existing_layers <- character(0)
existing_df     <- NULL

if (file.exists(output_csv)) {
  existing_df     <- read.csv(output_csv, stringsAsFactors = FALSE)
  existing_layers <- existing_df$layer
  message("Existing CSV found with ", nrow(existing_df), " layer(s) already processed.")
} else {
  message("No existing CSV found — processing all layers.")
}

# Filter to only TIFFs not yet in the CSV
all_layer_names <- tools::file_path_sans_ext(basename(tiff_files))
new_idx         <- !all_layer_names %in% existing_layers
tiff_files      <- tiff_files[new_idx]

if (length(tiff_files) == 0) {
  message("✓ All layers already processed. Nothing to do.")
  quit(save = "no")
}

message("Found ", length(tiff_files), " new TIFF file(s) to process.\n")

# Probability cut-points that define n_classes equal-frequency classes.
# For 4 classes: 0%, 25%, 50%, 75%, 100%  →  breaks at 25 / 50 / 75 percentiles
probs <- seq(0, 1, length.out = n_classes + 1)

results <- lapply(tiff_files, function(fp) {
  
  layer_name <- tools::file_path_sans_ext(basename(fp))
  message("  Processing: ", layer_name)
  
  r       <- terra::rast(fp)
  n_cells <- terra::ncell(r)
  
  # For large rasters, sample instead of loading all values into RAM.
  # 500k cells is more than enough for stable quantile estimates.
  if (n_cells > sample_threshold) {
    message("    (large raster: ", format(n_cells, big.mark = ","),
            " cells — sampling ", format(sample_size, big.mark = ","), " cells)")
    vals_raw <- terra::spatSample(r, size = sample_size,
                                  method = "regular",
                                  as.df = FALSE)[, 1]
    # Remove both R NAs and any value declared as the raster's NoData flag
    naflag   <- terra::NAflag(r)
    if (!is.na(naflag)) vals_raw <- vals_raw[vals_raw != naflag]
    values   <- vals_raw[!is.na(vals_raw)]
  } else {
    values <- terra::values(r, mat = FALSE)
    naflag <- terra::NAflag(r)
    if (!is.na(naflag)) values <- values[values != naflag]
    values <- values[!is.na(values)]
  }
  
  if (length(values) == 0) {
    warning("  [SKIP] ", layer_name, " — no non-NA values found.")
    return(NULL)
  }
  
  unique_vals <- unique(values)
  
  # ── Edge case: only one unique value → assign highest class ─────────────────
  if (length(unique_vals) == 1) {
    message("  [NOTE] ", layer_name,
            " has a single unique value (", unique_vals,
            "). Assigning highest class (", n_classes, ").")
    
    # breaks 0 … (n_classes-1) = 0, only break_n_classes holds the actual value
    breaks <- c(rep(0, n_classes), unique_vals)
    
    return(data.frame(
      layer      = layer_name,
      file_path  = fp,
      n_classes  = n_classes,
      note       = "single_value",
      # break_0 = lower bound (inclusive min), break_4 = upper bound
      t(setNames(breaks, paste0("break_", 0:n_classes)))
    ))
  }
  
  # ── Normal case: compute quantile breaks ───────────────────────────────────
  breaks <- as.numeric(quantile(values, probs = probs, na.rm = TRUE))
  
  # Guarantee strictly increasing breaks so reclassification is unambiguous:
  # if duplicate break values exist (e.g. highly skewed data), nudge them apart
  # so each class interval is non-empty.
  for (i in seq(2, length(breaks))) {
    if (breaks[i] <= breaks[i - 1]) {
      breaks[i] <- breaks[i - 1] + .Machine$double.eps * abs(breaks[i - 1]) * 1e6
    }
  }
  
  data.frame(
    layer      = layer_name,
    file_path  = fp,
    n_classes  = n_classes,
    note       = "ok",
    t(setNames(breaks, paste0("break_", 0:n_classes)))
  )
})

# Drop any NULLs (skipped layers)
results <- Filter(Negate(is.null), results)

if (length(results) == 0) stop("No new layers could be processed.")

new_df <- do.call(rbind, results)
rownames(new_df) <- NULL

# Combine with existing rows (if any) and write
breaks_df <- if (!is.null(existing_df)) rbind(existing_df, new_df) else new_df

write.csv(breaks_df, output_csv, row.names = FALSE, quote = TRUE)

message("\n✓ CSV updated: ", normalizePath(output_csv, mustWork = FALSE))
message("  ", nrow(new_df), " new layer(s) added | ",
        nrow(breaks_df), " total layer(s) in CSV.")
print(new_df[, c("layer", "note", paste0("break_", 0:n_classes))])
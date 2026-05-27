# =============================================================================
# ALTERNATIVE VERSION — NA-safe quantile estimation with controlled sampling
#
# This script computes Q25, Q50, and Q75 thresholds for raster layers
# after removing NoData and NA values.
#
# For each raster, valid pixel values are first estimated and extracted
# (either fully or via sampling depending on raster size).
#
# If the number of valid pixels is below a defined threshold, all valid
# values are loaded. If the raster is large, a fixed number of valid pixels
# are sampled using systematic (regular) sampling.
#
# Quantiles are then computed from the resulting set of valid pixel values.
#
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
  #naflag <- terra::NAflag(r) #retrieves the internal value used by the raster to represent missing data
  
  # ── Load ALL values and clean BEFORE deciding whether to SAMPLE ────────────
  # (This is the key difference from the original script)
  # we load everything first so NA pixels are excluded
  # from the pool before any sampling takes place.

   n_total <- terra::ncell(r)
  # # Estimate number of non-NA cells without loading everything
  #
   n_clean <- terra::global(!is.na(r), "sum", na.rm=TRUE)[1,1]
  # #creates a TRUE/FALSE raster: #TRUE where pixels are NA  #FALSE where pixels contain data
  # #[1,1] extracts the actual numeric value from the returned table.
  #
     message("    Estimated clean pixels: ",
             format(n_clean, big.mark=","))

    set.seed(1)
    
    if (n_clean > sample_threshold) {
      message("    Sampling valid pixels directly...")
      values <- terra::spatSample(
        r,
        size     = sample_size,
        method   = "regular",
        na.rm    = TRUE,
        values   = TRUE,
        as.points = FALSE
      )[,1]
  
    } else {
  # #
      message("    Small raster: loading valid pixels only")
  # #
      values <- terra::values(r, na.rm = TRUE)
    }
  
  #-----------------
  
  # NO SAMPLING — load all valid pixels
   # values <- terra::values(r, mat = FALSE, 
   #                         na.rm = TRUE)
   # message("    Loaded all valid values.")
  
  #_---------------
  
  # ── Cap filter: layers whose name contains ASIS or PEy ──────────────
  cap_filter <- grepl("ASIS|PEy", layer_name, ignore.case = FALSE)
  
  if (cap_filter) {
    values <- values[values <= 100] #Anything above 100 is removed.
    
    message("    [cap filter] '", layer_name,
            "' — keeping only values <= 100")
  }
  
  n_clean <- length(values) #counts how many usable values remain after:
  
  #removing NA values
  #sampling or loading
  #applying the optional cap filter
  
  message("    Final usable values: ",
          format(n_clean, big.mark=","))
  
  if (n_clean == 0) { #checks whether no usable values remain.
    warning("  [SKIP] ", layer_name,
            " — no valid values found.")
    return(NULL) #skips that raster and moves to the next one.
  }
  
  # Get unique non-NA values
  unique_vals <- unique(as.vector(values)) #extracts all distinct values from the raster.
  
  # ── Edge cases
  # ONLY ONE unique value → assign highest class ──────────────
  if (length(unique_vals) == 1) { #checks whether the raster contains only one unique value
    message("  [NOTE] ", layer_name,
            " has a single unique value (", unique_vals,
            "). Assigning highest class (", n_classes, ").")
    
    return(data.frame(
      layer          = layer_name,
      file_path      = fp,
      n_classes      = n_classes,
      n_clean_pixels = n_clean,
      note           = "single_value",
      Q25            = unique_vals,
      Q50            = unique_vals,
      Q75            = unique_vals
    ))
  } 
  
  # ── LAND CHANGE layers already classified (1–4) ─────────────────────────
  # These layers already contain:
  # 1 = Low
  # 2 = Moderate
  # 3 = High
  # 4 = Very high
  # → skip quantile calculation
  
  landchange_filter <- grepl(
    "cropland_change|pastureland_change|forestland_change",
    layer_name,
    ignore.case = TRUE #uppercase/lowercase differences do not matter.
  )
  
  if (landchange_filter) {
    
    message("  [NOTE] ", layer_name,
            " already classified into 4 categories → skipping quantiles.")
    
    return(data.frame(
      layer          = layer_name,
      file_path      = fp,
      n_classes      = n_classes,
      n_clean_pixels = n_clean,
      note           = "already_classified",
      Q25            = 1,
      Q50            = 2,
      Q75            = 3
    ))
  }
  
  # ── Normal case: compute the three quantile thresholds ───────────────────
  # Classes after reclassification:
  #   Class 1 (Low)      : value <  Q25
  #   Class 2 (Med-Low)  : Q25 <= value <  Q50
  #   Class 3 (Med-High) : Q50 <= value <  Q75
  #   Class 4 (High)     : value >= Q75
  breaks <- as.numeric(quantile(values, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)) #further removes NA values if any somehow remain.
  
  data.frame(
    layer          = layer_name,
    file_path      = fp,
    n_classes      = n_classes,
    n_clean_pixels = n_clean,
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
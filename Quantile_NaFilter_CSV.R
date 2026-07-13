# =============================================================================
# ALTERNATIVE VERSION — NA-safe quantile estimation with controlled sampling

#The script includes filter for already classified tifs

# The script computes Q25, Q50, and Q75 thresholds for raster layers
# after removing NoData and NA values.

# For each raster, valid pixel values are first estimated and extracted
# (either fully or via sampling depending on raster size).

# If the number of valid pixels is below a defined threshold, all valid
# values are loaded. If the raster is large, a fixed number of valid pixels
# are sampled using systematic (regular) sampling.

# Quantiles are then computed from the resulting set of valid pixel values.


#EDITS

#load specific layers - match name - include the specific reclassification****




# =============================================================================

library(terra)

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
input_dir        <- "/home/jovyan/FAO_climate_risks_team/notebooks/Andrea/Hotspot mapping/Repository"
output_csv       <- "/home/jovyan/FAO_climate_risks_team/notebooks/Andrea/Hotspot mapping/Repository/Breaks_summary.csv"
OUT_BASE         <- "/home/jovyan/FAO_climate_risks_team/notebooks/Andrea/Hotspot mapping/Repository"
OUT_RASTER       <- file.path(OUT_BASE, "RASTER")
n_classes        <- 4
sample_threshold <- 5e6
sample_size      <- 500000
# ──────────────────────────────────────────────────────────────────────────────

#For ALL THE TIFF in the repository-------
tiff_files <- list.files(input_dir,
                         pattern     = "\\.tif{1,2}$",
                         full.names  = TRUE,
                         ignore.case = TRUE,
                         recursive   = TRUE)

#For ONLY SPECIFIC tiff-----
#tiff_files <- c(
#  "/home/jovyan/FAO_climate_risks_team/notebooks/Andrea/Hotspot mapping/Repository/Raw_EXP/Livestock_density.tif"
#)

if (length(tiff_files) == 0) stop("No TIFF files found in: ", input_dir)
message("Found ", length(tiff_files), " TIFF file(s) to process.\n")

results <- lapply(tiff_files, function(fp) {
  
  layer_name <- tools::file_path_sans_ext(basename(fp))
  message("  Processing: ", layer_name)
  
  # ── Skip if already processed ───────────────────────────────────────────────
  out_file <- file.path(OUT_RASTER, paste0(layer_name, ".tif"))
  if (file.exists(out_file)) {
    message("  [SKIP] Already processed: ", layer_name)
    return(NULL) #A rerun on an existing output → skipped silently, not written to CSV, not processed again in Reclassification.R
    # return(data.frame(
    #   layer          = layer_name,
    #   file_path      = fp,
    #   n_classes      = n_classes,
    #   n_clean_pixels = NA,
    #   note           = "already_reclassified",
    #   Q25            = NA,
    #   Q50            = NA,
    #   Q75            = NA
    # ))
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # ── Manual reclassification overrides ───────────────────────────────────────
  # Add one entry per layer that needs fixed thresholds instead of quantiles.
  # Copy-paste the block and change the name + Q values as needed.
  
  manual_breaks <- list(
      
      "Low-lying_areas"    = list(Q25 = 20,    Q50 = 10,   Q75 = 5,    inverted = TRUE),
      "Population_density" = list(Q25 = 10,    Q50 = 250,  Q75 = 1000, inverted = FALSE),
      "HDI"                = list(Q25 = 0.8,   Q50 = 0.7,  Q75 = 0.55, inverted = TRUE),
      "Poverty"            = list(Q25 = 0.001, Q50 = 14,   Q75 = 67.8, inverted = FALSE),
      "IDPs"               = list(Q25 = 100,   Q50 = 1000, Q75 = 10000,inverted = FALSE),
      "Food_insecurity"    = list(Q25 = 4.5,   Q50 = 7.2,  Q75 = 27.4, inverted = FALSE),
      "Conflict_intensity" = list(Q25 = 0.01,  Q50 = 4,    Q75 = 5,    inverted = FALSE),
      "Weather_observation"= list(Q25 = 10000, Q50 = 2500, Q75 = 625,  inverted = TRUE),
      "Electricity_access" = list(Q25 = 30,    Q50 = 60,   Q75 = 90,   inverted = FALSE),
      "Irrigated_areas"    = list(Q25 = 2,     Q50 = 10,   Q75 = 50,   inverted = FALSE),
      "Rural_catchment"    = list(Q25 = 29,    Q50 = 15,   Q75 = 10,   inverted = TRUE),
      "Livestock_density"  = list(Q25 = 0,     Q50 = 10,   Q75 = 100,  inverted = FALSE),
      "GDI"                = list(Q25 = 0.966, Q50 = 0.984,Q75 = 1.001,inverted = TRUE)
      
    )
    
  # ────────────────────────────────────────────────────────────────────────────
  
  if (layer_name %in% names(manual_breaks)) {
    mb <- manual_breaks[[layer_name]]
    message("  [MANUAL] Using fixed thresholds for: ", layer_name)
    return(data.frame(
      layer          = layer_name,
      file_path      = fp,
      n_classes      = n_classes,
      n_clean_pixels = NA,
      note           = ifelse(isTRUE(mb$inverted), "manual_inverted", "manual_breaks"),
      Q25            = mb$Q25,
      Q50            = mb$Q50,
      Q75            = mb$Q75
    ))
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # ── Layers already classified at source (use raster values as-is) ────────── so
  #A source-classified layer → always explicitly listed by name, flagged correctly, and handled as-is in Reclassification.R
  already_classified <- c(
    "Cropland_change",
    "Pastureland_change",
    "Forestland_change",
    "Human_pressure"
  )
  
  matched_classified <- already_classified[tolower(already_classified) == tolower(layer_name)]
  
  if (length(matched_classified) == 1) {
    message("  [SOURCE-CLASSIFIED] Using raster values as-is: ", layer_name)
    return(data.frame(
      layer          = layer_name,
      file_path      = fp,
      n_classes      = n_classes,
      n_clean_pixels = NA,
      note           = "already_reclassified",
      Q25            = NA,
      Q50            = NA,
      Q75            = NA
    ))
  }
  
  # ──────────────────────────────────────────────────────────────────────────── 
  # ────────────────────────────────────────────────────────────────────────────
  r      <- terra::rast(fp)
  #naflag <- terra::NAflag(r) #retrieves the internal value used by the raster to represent missing data
  
  # ── Load ALL values and clean BEFORE deciding whether to SAMPLE ────────────
  # (This is the key difference from the original script)
  # we load everything first so NA pixels are excluded
  # from the pool before any sampling takes place.
  
  n_total <- terra::ncell(r)
  
  n_clean <- terra::global(!is.na(r), "sum", na.rm=TRUE)[1,1] #
  #creates a TRUE/FALSE raster: #TRUE where pixels are NA  #FALSE where pixels contain data
  #[1,1] extracts the actual numeric value from the returned table.
  
  message("    Estimated clean pixels: ",
          format(n_clean, big.mark=","))
  
  set.seed(1)
  
  if (n_clean > sample_threshold) {
    message("    Sampling valid pixels...")
    values <- terra::spatSample(
      r,
      size     = sample_size,
      method   = "regular",
      na.rm    = TRUE, #the sampling is done only on valid pixels (NA pixels are skipped).
      values   = TRUE,
      as.points = FALSE
    )[,1]
    
  } else {
    
    # ── Use spatSample instead of terra::values() to avoid loading the full
    #    raster extent into RAM (total cells >> valid cells for sparse layers)
    message("    Loading valid pixels via spatSample...")
    values <- terra::spatSample(
      r,
      size      = max(n_clean, 1L),   # request at most the known valid count
      method    = "regular",
      na.rm     = TRUE,
      values    = TRUE,
      as.points = FALSE
    )[,1]
  }
  
  #-----------------
  
  # NO SAMPLING — load all valid pixels
  # values <- terra::values(r, mat = FALSE,
  #                         na.rm = TRUE)
  # message("    Loaded all valid values.")
  
  #_---------------
  
  # ── Cap filter: layers whose name contains ASIS or PEy ──────────────
  cap_filter <- grepl("ASIS", layer_name, ignore.case = TRUE)
  
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
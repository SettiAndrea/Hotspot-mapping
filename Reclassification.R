# =============================================================================
# 
# Pipeline per layer (loops over every TIFF listed in the quantile CSV - EXCEPT the ones already processed):
#   0. Libraries & settings
#   1. AOI definition  (GAUL0 | GAUL1 | EXTENT | SHAPEFILE)
#   2. Reference layer (sets target resolution)
#   3. For each NEW layer in CSV:
#       3a. Load raw TIFF
#       3b. Reproject to TARGET_CRS
#       3c. NoData cleaning  (metadata flag, fill values, ASIS/PyAEZ filter) - this takes 5 mins for 3billion pixels'layer
#       3d. Crop + mask to AOI
#       3e. Reclassify into 4 bins using CSV quantile breaks (or single value)
#       3f. Resample to reference resolution  (method = "near")
#       3g. Export .tif + .png
# =============================================================================

#Adjustments in QUANTILE_CSV.R script to feed the RECLASSIFICATION.R 
#Different classifications: 

# - Land changes (layers with already 4 values) 

# - Vulnerability layers applying a categorization given by the resource of the data 

# - Option to calculate the quantile at national level - Only if needed in the future

# ── 0. LIBRARIES & SETTINGS ───────────────────────────────────────────────────

library(terra)
library(sf)
library(tidyverse)

# ── Paths ─────────────────────────────────────────────────────────────────────
GAUL0_PATH  <- "/home/jovyan/FAO_climate_risks_team/notebooks/Andrea/Hotspot mapping/GAUL0/g2015_2014_0.shp"
REF_PATH    <- "/home/jovyan/FAO_climate_risks_team/notebooks/Andrea/Hotspot mapping/Reference/Mangroves.tif"
BREAKS_CSV  <- "/home/jovyan/FAO_climate_risks_team/notebooks/Andrea/Hotspot mapping/Repository/Breaks_summary.csv"           

OUT_BASE        <- file.path(getwd(), "OUT")
OUT_EXP         <- file.path(OUT_BASE, "Reclassified_EXP")
OUT_VUL         <- file.path(OUT_BASE, "Reclassified_VUL")
OUT_AC          <- file.path(OUT_BASE, "Reclassified_AC")

dir.create(OUT_EXP, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_VUL, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_AC,  showWarnings = FALSE, recursive = TRUE)

message("✓ Output folders ready under: ", OUT_BASE)

# ── AOI settings ──────────────────────────────────────────────────────────────
aoi_mode     <- "GAUL0"       # "GAUL0" | "GAUL1" | "EXTENT" | "SHAPEFILE"
country_name <- "Uganda"
country_ext  <- ext(-61.1, -60.85, 13.69, 14.12)   # used only if aoi_mode = "EXTENT"
shp_path     <- "Shape/Mauritius.shp"  # SHAPEFILE mode

# ── Classification labels & colours (shared across all layers) ──────────────── To be added ASIS classification? 
bin_labels  <- c("1" = "Low", "2" = "Moderate", "3" = "High", "4" = "Very high")

bin_colours <- c("Low"       = "#d9ef8b",
                 "Moderate"  = "#fee08b",
                 "High"      = "#f46d43",
                 "Very high" = "#a50026")

# =============================================================================
# 1. AOI DEFINITION
# =============================================================================

gaul0_ref  <- st_read(GAUL0_PATH, quiet = TRUE)
TARGET_CRS <- st_crs(gaul0_ref)$wkt

get_aoi <- function(mode, country_name = NULL, shp_path = NULL, country_ext = NULL) {
  
  message("── AOI mode: ", mode, " ─────────────────────────────────────")
  
  if (mode == "GAUL0") {
    
    gaul <- st_read(GAUL0_PATH, quiet = TRUE)
    aoi  <- gaul %>% filter(tolower(ADM0_NAME) == tolower(country_name))
    if (nrow(aoi) == 0) stop("❌ No GAUL0 match for '", country_name, "'")
    if (nrow(aoi) > 1) {
      message("⚠️  ", nrow(aoi), " features matched → using all")
      message("   Matched: ", paste(unique(aoi$ADM0_NAME), collapse = ", "))
    }
    aoi <- st_transform(aoi, TARGET_CRS)
    return(vect(aoi))
    
  } else if (mode == "GAUL1") {
    
    gaul1 <- st_read("/home/jovyan/FAO_climate_risks_team/notebooks/Andrea/Hotspot mapping/GAUL1/gaul_2025_l1.shp", quiet = TRUE)
    aoi   <- gaul1 %>% filter(tolower(gaul1_name) == tolower(country_name))
    if (nrow(aoi) == 0) stop("❌ No GAUL1 match for '", country_name, "'")
    if (nrow(aoi) > 1) message("⚠️  ", nrow(aoi), " features matched → using all")
    aoi <- st_transform(aoi, TARGET_CRS)
    return(vect(aoi))
    
  } else if (mode == "EXTENT") {
    
    if (is.null(country_ext)) stop("❌ country_ext not provided")
    return(as.polygons(country_ext, crs = TARGET_CRS))
    
  } else if (mode == "SHAPEFILE") {
    
    if (is.null(shp_path)) stop("❌ shp_path not provided")
    shp <- st_read(shp_path, quiet = TRUE)
    shp <- st_make_valid(shp)
    if (is.na(st_crs(shp))) {
      st_crs(shp) <- TARGET_CRS
      message("⚠️  CRS missing → TARGET_CRS assigned")
    } else {
      shp <- st_transform(shp, TARGET_CRS)
    }
    return(vect(shp))
    
  } else {
    stop("❌ Invalid aoi_mode. Choose: GAUL0, GAUL1, EXTENT, SHAPEFILE")
  }
}

mask_vect <- get_aoi(aoi_mode, country_name, shp_path, country_ext)
bb        <- ext(mask_vect)
plot_xlim <- c(bb$xmin, bb$xmax)
plot_ylim <- c(bb$ymin, bb$ymax)
message(sprintf("  X: %.4f → %.4f | Y: %.4f → %.4f",
                plot_xlim[1], plot_xlim[2], plot_ylim[1], plot_ylim[2]))
message("✓ AOI ready\n")


# =============================================================================
# 2. REFERENCE LAYER  (resolution target)
# =============================================================================

message("── Loading reference layer ──────────────────────────────")
r_ref <- rast(REF_PATH)
r_ref <- crop(r_ref, mask_vect)
r_ref <- mask(r_ref, mask_vect)
message("✓ Reference layer ready\n")


# =============================================================================
# 3. READ QUANTILE BREAKS CSV
# =============================================================================

breaks_df  <- read.csv(BREAKS_CSV,
                       stringsAsFactors = FALSE)
break_cols <- c("Q25", "Q50", "Q75") #every raster will be reclassified into four categories using the three quantile thresholds.
n_classes  <- 4   # always 4: =Q75

message("Loaded ", nrow(breaks_df),
        " layer(s) from CSV: ", BREAKS_CSV)
message("Using fixed 4-class quantile scheme (Q25/Q50/Q75)\n")


# =============================================================================
# 4. MAIN LOOP — one iteration per layer in the CSV
# =============================================================================

for (i in seq_len(nrow(breaks_df))) {
  row        <- breaks_df[i, ]
  layer_name <- row$layer
  fp         <- row$file_path
  note       <- row$note
  
  # Process ONLY SPECIFIC LAYERS
  if (!(layer_name %in% c("Livestock_density"))) {
    next
  }
  
  message("Processing only: ", layer_name)
  
  message(rep("=", 70))
  message("  Layer ", i, "/", nrow(breaks_df), ": ", layer_name)
  message(rep("=", 70))
  
  if (!file.exists(fp)) { #check if it exist, ! menas NOT
    warning("  [SKIP] File not found: ", fp)
    next #skip this iteration and move to next layer
  }
  
  # ── Determine output folder from source path (AC / EXP / VUL) ─────────────
  if (grepl("Raw_AC", fp, fixed = TRUE)) {
    out_dir <- OUT_AC
  } else if (grepl("Raw_EXP", fp, fixed = TRUE)) {
    out_dir <- OUT_EXP
  } else if (grepl("Raw_VUL", fp, fixed = TRUE)) {
    out_dir <- OUT_VUL
  } else {
    warning("  [WARN] Could not determine category (AC/EXP/VUL) for: ", layer_name,
            " — saving to OUT_BASE instead")
    out_dir <- OUT_BASE
  }
  
  # ── 3a. Load ───────────────────────────────────────────────────────────────
  r <- rast(fp)
  message(sprintf("  Native CRS : %s", crs(r, describe = TRUE)$code))
  message(sprintf("  Native res : %.6f x %.6f", res(r)[1], res(r)[2]))
  message(sprintf("  Native NoData flag: %s",
                  ifelse(is.na(NAflag(r)), "NA", NAflag(r))))
  
  # ── 3b. Reproject to TARGET_CRS ────────────────────────────────────────────
  if (!same.crs(r, TARGET_CRS)) {
    r <- project(r, TARGET_CRS, method = "bilinear")
    message("  ✓ Reprojected to TARGET_CRS")
  }
  
  # ── 3c. NoData cleaning ────────────────────────────────────────────────────
  # Respect declared metadata NoData flag
  target_nodata <- NAflag(r)
  if (!is.na(target_nodata) && !is.nan(target_nodata)) {
    r[r == target_nodata] <- NA
    message(sprintf("  ✓ Metadata NoData (%s) → NA", target_nodata))
  }
  
  # Common fill values
  r[r %in% c(-9999, -32768)] <- NA
  
  # Extreme floating-point fill values (CORDEX / CMIP / NetCDF)
  r[r >  1e20] <- NA
  r[r < -1e20] <- NA
  
  # ASIS / land chenages filter (values > 100 are invalid)
  if (str_detect(layer_name,
                 regex("asis|pey|asy|cropland|pastureland|forestland",
                       ignore_case = TRUE))) {
    
    r[r > 100] <- NA
  } 
  message("  ✓ NoData cleaning complete")
  
  # ── 3d. Crop + mask to AOI ─────────────────────────────────────────────────
  # Check extents overlap before cropping — national layers may not cover the AOI
  layer_ext <- ext(r)
  aoi_ext   <- ext(mask_vect)
  
  x_overlap <- layer_ext$xmin < aoi_ext$xmax && layer_ext$xmax > aoi_ext$xmin
  y_overlap <- layer_ext$ymin < aoi_ext$ymax && layer_ext$ymax > aoi_ext$ymin
  
  if (!x_overlap || !y_overlap) {
    warning(sprintf(
      "  [SKIP] '%s' does not overlap with the AOI. Layer: X[%.2f,%.2f] Y[%.2f,%.2f] | AOI: X[%.2f,%.2f] Y[%.2f,%.2f]",
      layer_name,
      layer_ext$xmin, layer_ext$xmax, layer_ext$ymin, layer_ext$ymax,
      aoi_ext$xmin,   aoi_ext$xmax,   aoi_ext$ymin,   aoi_ext$ymax))
    next
  }
  
  r <- crop(r, mask_vect)
  r <- mask(r, mask_vect)
  message("  ✓ Cropped and masked to AOI")
  
  # ── 3e. Reclassify using CSV quantile breaks ───────────────────────────────
 
  if (note == "already_reclassified") { 
    r_classified <- r
  } else if (note == "single_value") { #Was this raster marked as a single-value layer in the previous script?
    #single_val   <- breaks[1]   # all Q25/Q50/Q75 equal, the script just takes one of them.
    r_classified <- ifel(!is.na(r),
                         n_classes, NA) #the script simply assigns every valid pixel to the highest class:
    message("  ✓ Single-value layer → class ",
            n_classes)
    
  } else if (note == "manual_inverted") {
    q25 <- breaks[1] #reads the thersholds
    q50 <- breaks[2]
    q75 <- breaks[3]
    
    # terra::classify() requires from < to in every row,
    # so we sort breaks ascending but flip the class numbers
    b_lo <- min(q75, q50, q25)   # find and store the lowest threshold e.g. 10
    b_mid <- median(c(q75, q50, q25))  # find and store the mid threshold e.g.15
    b_hi <- max(q75, q50, q25)   # find and store the highest threshold e.g.29
    
    r_classified <- classify(r,
                             matrix(c(-Inf,  b_lo,  4,   # < 10 (LOWEST (Q25) ASSIGNED TO 4 → Very High
                                      b_lo,  b_mid, 3,   # 10–15 → High
                                      b_mid, b_hi,  2,   # 15–29 → Moderate
                                      b_hi,  Inf,   1),  # > 29  → Low
                                    ncol = 3, byrow = TRUE),
                             include.lowest = TRUE, others = NA)
    message("  ✓ Reclassified (INVERTED scale) into 4 bins")
    
  } else { #otherwise RECLASSIFY using the quantiles
    # 3 thresholds → 4 classes:
    #  1: value <  Q25
    #  2: Q25 <= value <  Q50
    #  3: Q50 <= value <  Q75
    #  4: value >= Q75
    q25 <- breaks[1]
    q50 <- breaks[2]
    q75 <- breaks[3]
    r_classified <- classify(r,
                             matrix(c(-Inf, q25, 1, #INF 
                                      q25, q50, 2,
                                      q50, q75, 3,
                                      q75,  Inf, 4), #INF
                                    ncol = 3, byrow = TRUE),
                             include.lowest = TRUE, others = NA)
    message("  ✓ Reclassified into 4 bins",
            ifelse(note == "capped_at_100",
                   " [capped at 100]", ""))
  }
  
  # Frequency table
  freq_tbl <- freq(r_classified)
  for (j in seq_len(nrow(freq_tbl))) {
    code  <- as.character(freq_tbl$value[j])
    label <- ifelse(code %in% names(bin_labels), bin_labels[[code]], paste("code", code))
    message(sprintf("    [%s] %-12s  %s pixels",
                    code, label, format(freq_tbl$count[j], big.mark = ",")))
  }
  
  # ── 3f. Resample to reference resolution ───────────────────────────────────
  r_final <- resample(r_classified, r_ref, method = "near")
  
  # Mask again after resampling
  r_final <- crop(r_final, mask_vect)
  r_final <- mask(r_final, mask_vect)
  
  message(sprintf("  ✓ Resampled to %.6f° x %.6f°", res(r_final)[1], res(r_final)[2]))
  message("  ✓ Cropped and masked again after resampling")
  
  # ── 3g. Export TIFF ────────────────────────────────────────────────────────
  out_tif <- file.path(out_dir, paste0(layer_name, ".tif"))
  writeRaster(r_final, out_tif,
              overwrite = TRUE,
              datatype  = "INT1U",
              NAflag    = 255L)
  message("  ✓ TIFF exported → ", out_tif)
  
  # ── 3g. Export PNG ─────────────────────────────────────────────────────────
  df <- as.data.frame(r_final, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "bin_code"
  
  df <- df %>%
    mutate(
      bin_label = bin_labels[as.character(bin_code)],
      bin_label = factor(bin_label, levels = unname(bin_labels))
    )
  
  mask_sf <- st_as_sf(mask_vect)
  
  p <- ggplot(df, aes(x = x, y = y, fill = bin_label)) +
    geom_tile() +
    geom_sf(data = mask_sf, fill = NA, color = "black",
            linewidth = 0.2, inherit.aes = FALSE) +
    coord_sf(xlim = plot_xlim, ylim = plot_ylim, expand = TRUE) +
    scale_fill_manual(
      values   = bin_colours,
      na.value = "transparent",
      guide    = guide_legend(direction = "horizontal", title.position = "top")
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    labs(
      title = paste(layer_name),
      fill  = "Class",
      x     = "Longitude",
      y     = "Latitude"
    ) +
    theme_bw() +
    theme(panel.grid      = element_blank(),
          legend.position = "bottom",
          plot.title      = element_text(face = "bold"))
  
  out_png <- file.path(out_dir, paste0(layer_name, ".png"))
  ggsave(out_png, plot = p, dpi = 300)
  message("  ✓ PNG  exported → ", out_png, "\n")

  }   # <-- closes the for loop

message(rep("=", 70))
message("✓ All layers processed.")
message("  Reclassified layers → ", OUT_BASE)
message(rep("=", 70))

#-------
  
# ── STEP 5: BUILD CUMULATIVE EXP / VUL / AC LAYERS ─────────────────────────────
# Stack every reclassified TIFF within each category folder and average across
# layers. Since every layer was already resampled onto r_ref inside the main
# loop above, all rasters in a folder already share the exact same grid — no
# additional alignment is needed here, just stack + average.

build_cumulative_layer <- function(folder, name, out_dir) {
  tif_files <- list.files(folder, pattern = "\\.tif{1,2}$", full.names = TRUE, ignore.case = TRUE)
  if (length(tif_files) == 0) {
    cat("No rasters found in", folder, "- skipping", name, "\n")
    return(NULL)
  }
  cat("Stacking", length(tif_files), "layer(s) for", name, ":\n  ",
      paste(basename(tif_files), collapse = ", "), "\n")
  
  stacked    <- terra::rast(tif_files)
  cumulative <- terra::app(stacked, fun = mean, na.rm = TRUE)
  names(cumulative) <- name
  
  out_tif <- file.path(out_dir, paste0(name, ".tif"))
  terra::writeRaster(cumulative, out_tif, overwrite = TRUE, datatype = "FLT4S")
  cat("✓ TIFF saved:", out_tif, "\n")
  
  # ── PNG export ──────────────────────────────────────────────────────────────
  cum_df <- terra::as.data.frame(cumulative, xy = TRUE, na.rm = TRUE)
  names(cum_df)[3] <- "value"
  
  mask_sf <- sf::st_as_sf(mask_vect)
  
  p <- ggplot(cum_df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    geom_sf(data = mask_sf, fill = NA, color = "black",
            linewidth = 0.2, inherit.aes = FALSE) +
    coord_sf(xlim = plot_xlim, ylim = plot_ylim, expand = TRUE) +
    scale_fill_gradientn(colors = c("#1a9850", "#fee08b", "#d73027"),
                         name = "Mean class\n(1=Low, 4=Very high)",
                         limits = c(1, 4)) +
    labs(title = name, x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(panel.grid      = element_blank(),
          legend.position = "bottom",
          plot.title      = element_text(face = "bold"))
  
  out_png <- file.path(out_dir, paste0(name, ".png"))
  ggsave(out_png, plot = p, dpi = 300)
  cat("✓ PNG  saved:", out_png, "\n\n")
  
  cumulative
}

OUT_CUMULATIVE <- file.path(OUT_BASE, "Cumulative")
dir.create(OUT_CUMULATIVE, showWarnings = FALSE, recursive = TRUE)

Cumulative_EXP <- build_cumulative_layer(OUT_EXP, "Cumulative_EXP", OUT_CUMULATIVE)
Cumulative_VUL <- build_cumulative_layer(OUT_VUL, "Cumulative_VUL", OUT_CUMULATIVE)
Cumulative_AC  <- build_cumulative_layer(OUT_AC,  "Cumulative_AC",  OUT_CUMULATIVE)

# ── Sanity check: confirm all three cumulative layers share the exact same grid ──
if (!is.null(Cumulative_EXP) && !is.null(Cumulative_VUL)) {
  cat("EXP vs VUL geometry match:\n")
  print(terra::compareGeom(Cumulative_EXP, Cumulative_VUL, stopOnError = FALSE))
}
if (!is.null(Cumulative_EXP) && !is.null(Cumulative_AC)) {
  cat("EXP vs AC geometry match:\n")
  print(terra::compareGeom(Cumulative_EXP, Cumulative_AC, stopOnError = FALSE))
}

message(rep("=", 70))
message("✓ Cumulative EXP / VUL / AC layers built.")
message("  Output → ", OUT_CUMULATIVE)
message(rep("=", 70))

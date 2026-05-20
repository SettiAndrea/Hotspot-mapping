# =============================================================================
# 
# Pipeline per layer (loops over every TIFF listed in the quantile CSV - EXCEPT the ones already processed):
#   0. Libraries & settings
#   1. AOI definition  (GAUL0 | GAUL1 | EXTENT | SHAPEFILE)
#   2. Reference layer (sets target resolution)
#   3. For each NEW layer in CSV:
#       3a. Load raw TIFF
#       3b. Reproject to TARGET_CRS
#       3c. NoData cleaning  (metadata flag, fill values, ASIS/PyAEZ filter) - this takes 5 mins for 3billion pixels
#       3d. Crop + mask to AOI
#       3e. Reclassify into 4 bins using CSV quantile breaks
#       3f. Resample to reference resolution  (method = "near")
#       3g. Export .tif + .png
# =============================================================================


# ── 0. LIBRARIES & SETTINGS ───────────────────────────────────────────────────

library(terra)
library(sf)
library(tidyverse)

# ── Paths ─────────────────────────────────────────────────────────────────────
GAUL0_PATH  <- "/Volumes/Andrea_GIS/g2015_2014_0/g2015_2014_0/g2015_2014_0.shp"
REF_PATH    <- "/Volumes/Andrea_GIS/Hotspot Mapping/R_test/Test_Buodaries_Standardization/Mangroves.tif"
BREAKS_CSV  <- "/Volumes/Andrea_GIS/Hotspot Mapping/R_test/Quantile/Test.csv"           

OUT_BASE    <- "/Volumes/Andrea_GIS/Hotspot Mapping/R_test/OUTPUT"
OUT_RASTER  <- file.path(OUT_BASE, "RASTER")
OUT_FIGURE  <- file.path(OUT_BASE, "FIGURE")

dir.create(OUT_RASTER, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIGURE, showWarnings = FALSE, recursive = TRUE)

# ── AOI settings ──────────────────────────────────────────────────────────────
aoi_mode     <- "GAUL0"       # "GAUL0" | "GAUL1" | "EXTENT" | "SHAPEFILE"
country_name <- "Uganda"
country_ext  <- ext(-61.1, -60.85, 13.69, 14.12)   # used only if aoi_mode = "EXTENT"
shp_path     <- "/Volumes/Andrea_GIS/Hotspot Mapping/R_test/Boundaries/CostaRica.shp"  # SHAPEFILE mode

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
    
    gaul1 <- st_read("/Volumes/Andrea_GIS/GAUL_2024_L1/GAUL_2024_L1.shp", quiet = TRUE)
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

breaks_df  <- read.csv(BREAKS_CSV, stringsAsFactors = FALSE) #This creates the dataframe in R (table in R with rows and columns) STRING FACTOR FALSE  says "Keep text as normal text."
break_cols <- sort(grep("^break_\\d+$", names(breaks_df), value = TRUE)) #ThiS automatically detects columns named BREAK_NUMBER - 
#grep finds matching text - ^ start of text , \\d+ means one or more digits, $end of text, TRUE give back not the positions but NAMES -SO IF I UPDATE TO 6 CLASSESS IT AUTOAMTICALLY DETECTS
n_classes  <- length(break_cols) - 1 #I NEED 5 breaks to create 4 classes (Classes define spaces between edges)

message("Loaded ", nrow(breaks_df), " layer(s) from CSV: ", BREAKS_CSV) #nrow reads the n of rows (each row = one raster layer)
message("Detected ", n_classes, " classes\n")


# =============================================================================
# 4. MAIN LOOP — one iteration per layer in the CSV
# =============================================================================

for (i in seq_len(nrow(breaks_df))) {
 #The above means "repeat the following code once for every row in the dataframe" 
  row        <- breaks_df[i, ] #This extracts ONE row from the dataframe (i iteration) \ follows the extraction of variables 
  layer_name <- row$layer
  fp         <- row$file_path
  note       <- row$note
  breaks     <- as.numeric(row[, break_cols]) #You select ONLY the break columns from the row.as numeric takes the values not text e.g "10"
  
  message(rep("=", 70))
  message("  Layer ", i, "/", nrow(breaks_df), ": ", layer_name)
  message(rep("=", 70))
  
  if (!file.exists(fp)) { #! menas NOT
    warning("  [SKIP] File not found: ", fp)
    next #skip this iteration and move to next layer
  }
  
  # Skip if output TIFF already exists (avoids reprocessing completed layers)
  out_tif_check <- file.path(OUT_RASTER, paste0(layer_name, ".tif"))
  if (file.exists(out_tif_check)) {
    message("  [SKIP] Output already exists → ", out_tif_check)
    next
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
  
  # ASIS / PyAEZ filter (values > 100 are invalid)
  if (str_detect(tolower(layer_name), "asis|pyaez|pey")) {
    r[r > 100] <- NA
    message("  ✓ ASIS/PyAEZ filter applied (> 100 → NA)")
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
  if (note == "single_value") {
    
    single_val   <- breaks[length(breaks)]
    r_classified <- ifel(r == single_val, n_classes, NA)
    message("  ✓ Single-value layer → all pixels assigned class ", n_classes)
    
  } else {
    
    # Build [from, to, class] matrix
    # First interval extended left so the global minimum is captured
    rcl <- matrix(NA_real_, nrow = n_classes, ncol = 3)
    for (cls in seq_len(n_classes)) {
      lo <- breaks[cls]
      hi <- breaks[cls + 1]
      if (cls == 1) lo <- lo - abs(lo) * 1e-9 - .Machine$double.eps
      rcl[cls, ] <- c(lo, hi, cls)
    }
    
    r_classified <- classify(r, rcl, include.lowest = TRUE, others = NA)
    message("  ✓ Reclassified into ", n_classes, " bins")
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
  message(sprintf("  ✓ Resampled to %.6f° x %.6f°", res(r_final)[1], res(r_final)[2]))
  
  # ── 3g. Export TIFF ────────────────────────────────────────────────────────
  out_tif <- file.path(OUT_RASTER, paste0(layer_name, ".tif"))
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
      title = paste(layer_name, "— standardised"),
      fill  = "Class",
      x     = "Longitude",
      y     = "Latitude"
    ) +
    theme_bw() +
    theme(panel.grid      = element_blank(),
          legend.position = "bottom",
          plot.title      = element_text(face = "bold"))
  
  out_png <- file.path(OUT_FIGURE, paste0(layer_name, ".png"))
  ggsave(out_png, plot = p, dpi = 300)
  message("  ✓ PNG  exported → ", out_png, "\n")
}

message(rep("=", 70))
message("✓ All layers processed.")
message("  TIFFs → ", OUT_RASTER)
message("  PNGs  → ", OUT_FIGURE)
message(rep("=", 70))

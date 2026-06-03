
# =============================================================================
# 01_compute_quantiles.R
# Computes 4-class quantile breaks for every TIFF in a folder
# For layers with only one value e.g.mangroves.tif it assign directly 4 (VH)
# Writes a CSV that the script Reclassification.R will use.
# It won't run tiffs already processed
# =============================================================================

library(terra)

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
input_dir  <- "/Volumes/Andrea_GIS/Hotspot Mapping/R_test/Quantile"   
output_csv <- "/Volumes/Andrea_GIS/Hotspot Mapping/R_test/Quantile/Test.csv"         
n_classes        <- 4                       # number of quantile classes
sample_threshold <- 5e6                     # rasters larger than this are sampled
sample_size      <- 500000                  # number of cells to sample and calculate the quantile
# ──────────────────────────────────────────────────────────────────────────────


tiff_files <- list.files(input_dir,
                         pattern  = "\\.tif{1,2}$",
                         full.names = TRUE,
                         ignore.case = TRUE)

if (length(tiff_files) == 0) stop("No TIFF files found in: ", input_dir) #length give the number of elements

# ── Always produce a new CSV based on th files in the folder ───────────────────────
message("Found ", length(tiff_files),
        " TIFF file(s) to process.\n")

#-----Sampling and NA data--------

results <- lapply(tiff_files, function(fp) {
  
  layer_name <- tools::file_path_sans_ext(basename(fp))
  message("  Processing: ", layer_name)
  
  r       <- terra::rast(fp) #Loads raster.
  naflag <- terra::NAflag(r)   # Reads the NoData flag stored in the raster file's metadata.
  #                              We fetch it once here, before any branching, so we don't repeat the call twice below.
  
  n_cells <- terra::ncell(r)
  set.seed(1)                #Fixes the random number generator so the sample is identical every time 
  if (n_cells > sample_threshold) {
    vals_raw <- terra::spatSample(
      r, size = sample_size,
      method = "regular",    #means a systematic grid, not random scatter 
      as.df  = FALSE)[, 1] #as.df = FALSE returns a plain matrix instead of a data frame — faster. 
                                    #The [, 1] at the end extracts just the first (and usually only) column as a vector. 
    
    # Remove NAs — right after sampling
    if (!is.na(naflag)) vals_raw <-
        vals_raw[vals_raw != naflag]    #If the raster actually has a NoData flag (i.e. naflag is not itself NA), 
    #                                   remove any sample values that equal that sentinel number. 
    values <- vals_raw[!is.na(vals_raw)] #Now remove any remaining R-level NAs (pixels that terra itself marked as missing). 
    
      #                                    The result, values, is a clean numeric vector ready for quantile calculation.
  } else { #FOR SMALL RASTERS — safe to load every pixel into memory
    values <- terra::values(r, mat = FALSE) #FALSE skips the matrix wrapper
    if (!is.na(naflag)) values <- 
        values[values != naflag] ##Same NA value removal as in the sampling branch — naflag <- terra::NAflag(r)   # 
  
    values <- values[!is.na(values)] #Remove R-level NAs. After this line, both branches arrive at the same place: a clean values vector with no missing values of any kind.
  }
    
    # ── Cap filter: layers whose name contains ASIS or PEy ─────────────────
    cap_filter <- grepl("ASIS|PEy", layer_name, ignore.case = FALSE)
    if (cap_filter) {
      values <- values[values <= 100]
      message("    [cap filter] '", layer_name,
              "' — keeping only values <= 100 (",
              format(length(values), big.mark = ","), " values remaining)")
    }
  
  if (length(values) == 0) {
    warning("  [SKIP] ", layer_name, " — no non-NA values found.")
    return(NULL)
  }
  
  unique_vals <- unique(values)
  
  # ── Edge case: only one unique value → assign highest class ─────────────────
  if (length(unique_vals) == 1) {
    message("  [NOTE] ", layer_name,
            " has a single unique value (",
            unique_vals,
            "). Assigning highest class (",
            n_classes, ").")
    
    return(data.frame(
      layer      = layer_name,
      file_path  = fp,
      n_classes  = n_classes,
      note       = "single_value",
      Q25        = unique_vals,
      Q50        = unique_vals,
      Q75        = unique_vals
    ))
  }
  
  # ── Normal case: compute quantile breaks ───────────────────────────────────
  breaks <- as.numeric(
    quantile(values,
             probs = c(0.25, 0.50, 0.75),
             na.rm = TRUE))
  
  data.frame(
    layer     = layer_name,
    file_path = fp,
    n_classes = n_classes,
    note      = ifelse(cap_filter, "capped_at_100", "ok"),
    Q25       = breaks[1],
    Q50       = breaks[2],
    Q75       = breaks[3]
  )
})


if (length(results) == 0) stop("No new layers could be processed.")

new_df <- do.call(rbind, results)
rownames(new_df) <- NULL

breaks_df <- new_df

write.csv(breaks_df, output_csv,
          row.names = FALSE, quote = TRUE)

message("\n✓ CSV written: ",
  normalizePath(output_csv, mustWork = FALSE))
message("  ", nrow(breaks_df), " layer(s) written to CSV.")
print(breaks_df[, c("layer", "note",
  "Q25", "Q50", "Q75")])

library(terra)
library(raster)
library(tidyverse)


#===== COPERNICUS MARINE 

### This one for Copernicus Marine vertical current data to complement 2024 data to Bluelink BRAN 2020 data which will not be updated anymore until new product released

## Note: The NetCDF files need to be downloaded from Copernicus Mairne sepearately in Python using the Copernicus Marine Toolbox

# Everything else in this function is in accordance with the Bluelink function 


extractWz_CM <- function(df, X, Y, datetime, bathy_path, folder_name, export_path, max_depth = -200, 
                         fill_gaps = TRUE, buffer = 10000, verbose = TRUE, export_step = TRUE, keep_nc_files = FALSE) {
  
  if ("Wz" %in% names(df)) {
    message("Previous Wz calculations found, skipping depth extraction")
  } else {
    if ("Depth" %in% names(df)) {
      message("Depth data found, skipping Depth extraction :D")
      df <- df %>%
        dplyr::mutate(Depth_Wz = ifelse(Depth > -5, -6, ifelse(Depth < max_depth, max_depth, Depth)),
                      Depth_Wz = round(Depth_Wz, digits = 1))
    } else {
      message("Loading NOAA bathymetry data")
      aux.bathy <- list.files(bathy_path, pattern = ".tiff")
      for (i in 1:length(aux.bathy)) {
        rast.bathy <- terra::rast(paste(bathy_path, aux.bathy[i], sep = "/"))
        if (i == 1) {
          bathy.tot <- rast.bathy
        }
        if (i > 1) {
          bathy.tot <- terra::merge(bathy.tot, rast.bathy)
        }
      }
      aux.depths <- round(
        terra::extract(
          x = bathy.tot, 
          y = df[,c(X, Y)]), 
        1)
      df$Depth <- aux.depths[,2]
      
      df <- df %>%
        dplyr::mutate(Depth_Wz = ifelse(Depth > -5, -6, ifelse(Depth < max_depth, max_depth, Depth)))
    }
  }
  
  if (dir.exists(folder_name) == FALSE) 
    dir.create(folder_name)
  df$aux.date <- substr(df[,which(names(df) == datetime)], 1, 7)
  
  if ("Wz" %in% names(df)) {
    message("Previous Wz calculations found. Continuing...")
    dates <- unique(df$aux.date[is.na(df$Wz)])  
  } else {
    dates <- unique(df$aux.date)
    df$Wz <- NA
  }
  
  nc_files <- list.files(folder_name, pattern = "\\.nc$", full.names = TRUE)
  
  if (verbose) {
    message("Processing Wz data")
    pb <-  txtProgressBar(min = 0, max = length(dates), initial = 0, style = 3, width = 60)
  }
  
  for (i in 1:length(dates)) {
    aux.date <- dates[i]
    
    for (nc_file in nc_files) {
      # message("Processing file: ", nc_file)
      nc.bran <- terra::rast(nc_file) 
      
      aux.names <- stringr::str_remove(names(nc.bran), pattern = "wo_depth=")
      aux.names <- stringr::str_split(aux.names, pattern = "_")          
      aux.depth <- NULL
      aux.time <- NULL
      for (size.depth in 1:length(aux.names)) {
        aux.depth <- c(aux.depth, aux.names[[size.depth]][1])
        aux.time <- c(aux.time, aux.names[[size.depth]][2])
      }
      aux.names <- data.frame(Depth = as.numeric(aux.depth), Time = as.numeric(aux.time))
      aux.names$Time <- as.Date("2024-01-01", tz = "UTC") + aux.names$Time - 1
      # message("aux.names: ", head(aux.names))
      
      index.day <- which(df$aux.date == dates[i])
      for (ii in 1:length(index.day)) {
        aux.day <- aux.names[which(as.character(aux.names$Time) == as.character(df[index.day[ii], which(names(df) == datetime)])), ]
        # message("Processing index: ", index.day[ii], " Depth_Wz: ", df$Depth_Wz[index.day[ii]])
        # message("aux.day before filtering: ", aux.day)
        
        # # Print depth ranges for debugging
        # message("Range of aux.day$Depth: ", range(aux.day$Depth))
        # message("Depth_Wz: ", df$Depth_Wz[index.day[ii]])
        
        # Correcting the comparison logic
        index.depth <- which(aux.day$Depth <= abs(df$Depth_Wz[index.day[ii]]))
        # message("index.depth: ", index.depth)
        aux.day <- aux.day[index.depth,]
        # message("aux.day after filtering: ", aux.day)
        
        if (nrow(aux.day) > 0) {
          aux.day$W <- NA
          
          for (iii in 1:nrow(aux.day)) {
            index.layer <- which(aux.names$Depth == aux.day$Depth[iii] & 
                                   aux.names$Time == aux.day$Time[iii])
            aux.bran <- nc.bran[[index.layer]]
            aux.val <- terra::extract(
              x = aux.bran, 
              y = df[index.day[ii], c(X, Y)])[1,2]
            # message("Extracted value for layer ", index.layer, ": ", aux.val)
            
            if (is.na(aux.val) & fill_gaps) {
              pos_sf <- 
                df[index.day[ii], c(X, Y)] %>%
                sf::st_as_sf(coords = c(1,2), crs = 4326, remove = FALSE)
              pos_sf <- sf::st_buffer(pos_sf, buffer)
              
              aux.val <- terra::extract(x = aux.bran, y = pos_sf)
              names(aux.val)[2] <- "Buffer"
              aux.val <- aux.val %>%
                dplyr::group_by(ID) %>%
                dplyr::summarise(Buffer_mean = mean(Buffer, na.rm = TRUE))
              aux.val <- aux.val$Buffer_mean
              aux.day$W[iii] <- aux.val
              # message("Filled NA with buffer mean: ", aux.val)
            } else {
              aux.day$W[iii] <- aux.val
            }
          }
          
          if (nrow(aux.day) > 0) {
            aux.day$Height <- NA
            for (iii in 1:nrow(aux.day)) {
              if (iii == 1) {
                aux.day$Height[iii] <- 0 - aux.day$Depth[iii]
              } else {
                aux.day$Height[iii] <- aux.day$Depth[iii - 1] - aux.day$Depth[iii]
              }
              # message("Layer ", iii, " Height: ", aux.day$Height[iii])
            }
          }
          if (nrow(aux.day) == 0) {
            aux.wz <- 0
          } else {
            aux.day$Each <- aux.day$W * aux.day$Height
            aux.wz <- sum(aux.day$Each, na.rm = TRUE) / sum(aux.day$Height, na.rm = TRUE)
          }
          if (is.na(aux.wz))
            aux.wz <- 0
          df$Wz[index.day[ii]] <- aux.wz
          # message("Calculated Wz for index ", index.day[ii], ": ", aux.wz)
          gc()
        } else {
          # message("No valid depth layers found for index ", index.day[ii])
        }
      }
    }
    
    if (export_step)
      write.csv(df[,-which(names(df) == "aux.date")], paste0(export_path, ".csv"), row.names = FALSE)
    
    if (!keep_nc_files) {
      nc.names <- list.files(folder_name, pattern = ".nc")
      nc.names <- paste(folder_name, nc.names, sep = "/")
      invisible(file.remove(nc.names))
    }
    if (verbose)
      setTxtProgressBar(pb, i)     
  }
  if (verbose)
    close(pb) 
  
  if (!keep_nc_files)
    invisible(file.remove(folder_name))
  return(df)
}

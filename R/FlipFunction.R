#' inversionEjeY(muestra1)
#' 
#' This function flips the coordinates of a sample along the Y-axis (newrow = 77 - row).
#' It also flips the image data along the Y-axis.
#' 
#' @param muestra1 A sample object containing the image and coordinates.
#' @return A new sample object with updated coordinates and flipped image data.
#' @export 
inversionEjeY <- function(muestra1) {
  # Save the dataframe of coordinates from the sample
  coord_df <- muestra1@images$slice1@coordinates
  
  # Sort the dataframe by row and column
  sorted_coord_df <- coord_df[order(coord_df$row, coord_df$col),]
  
  # Create lists to store image row/column mappings for each row/col
  rowMapping <- list()
  colMapping <- list()
  
  # Loop through each barcode to store row and column mappings
  for (barcode in rownames(sorted_coord_df)) {
    rowValue <- as.character(sorted_coord_df[barcode, "row"])
    imageRowValue <- sorted_coord_df[barcode, "imagerow"]
    
    colValue <- as.character(sorted_coord_df[barcode, "col"])
    imageColValue <- sorted_coord_df[barcode, "imagecol"]
    
    if (!(imageRowValue %in% rowMapping)) {
      rowMapping[[rowValue]] <- imageRowValue
    }
    
    if (!(imageColValue %in% colMapping)) {
      colMapping[[colValue]] <- imageColValue
    }
  }
  
  # Create a copy of the dataframe
  updated_coord_df <- coord_df
  
  # Update the row and image row values
  for (barcode in rownames(coord_df)) {
    rowValue <- coord_df[barcode, "row"]
    colValue <- coord_df[barcode, "col"]
    
    # Flip the row value along the Y-axis
    newRowValue <- 77 - rowValue
    
    # Update the new row and image row values in the dataframe
    updated_coord_df[barcode, "row"] <- newRowValue
    updated_coord_df[barcode, "imagerow"] <- rowMapping[[as.character(newRowValue)]]
  }
  
  # Create a new sample object with the updated coordinates
  updated_sample <- muestra1
  new_names <- c()
  
  for (i in 1:length(muestra1@meta.data$name)) {
    new_names <- c(new_names, paste0(muestra1@meta.data$name[i], "_Xinvert"))
  }
  
  updated_sample@meta.data$name <- new_names
  updated_sample@images$slice1@coordinates <- updated_coord_df
  
  # Flip the image data along the Y-axis
  flipped_image <- array(NA, dim = dim(muestra1@images$slice1@image))
  
  for (i in 1:dim(muestra1@images$slice1@image)[3]) {
    flipped_image[,,i] <- muestra1@images$slice1@image[nrow(muestra1@images$slice1@image):1, , i]
  }
  
  updated_sample@images$slice1@image <- flipped_image
  
  return(updated_sample)
}

#' inversionEjeX(muestra1)
#' 
#' This function flips the coordinates of a sample along the X-axis (newcol = 127 - col).
#' It also flips the image data along the X-axis.
#' 
#' @param muestra1 A sample object containing the image and coordinates.
#' @return A new sample object with updated coordinates and flipped image data.
#' @export
inversionEjeX <- function(muestra1) {
  # Save the dataframe of coordinates from the sample
  coord_df <- muestra1@images$slice1@coordinates
  
  # Sort the dataframe by row and column
  sorted_coord_df <- coord_df[order(coord_df$row, coord_df$col),]
  
  # Create lists to store image row/column mappings for each row/col
  rowMapping <- list()
  colMapping <- list()
  
  # Loop through each barcode to store row and column mappings
  for (barcode in rownames(sorted_coord_df)) {
    rowValue <- as.character(sorted_coord_df[barcode, "row"])
    imageRowValue <- sorted_coord_df[barcode, "imagerow"]
    
    colValue <- as.character(sorted_coord_df[barcode, "col"])
    imageColValue <- sorted_coord_df[barcode, "imagecol"]
    
    if (!(imageRowValue %in% rowMapping)) {
      rowMapping[[rowValue]] <- imageRowValue
    }
    
    if (!(imageColValue %in% colMapping)) {
      colMapping[[colValue]] <- imageColValue
    }
  }
  
  # Create a copy of the dataframe
  updated_coord_df <- coord_df
  
  # Update the column and image column values
  for (barcode in rownames(coord_df)) {
    rowValue <- coord_df[barcode, "row"]
    colValue <- coord_df[barcode, "col"]
    
    # Flip the column value along the X-axis
    newColValue <- 127 - colValue
    
    # Update the new column and image column values in the dataframe
    updated_coord_df[barcode, "col"] <- newColValue
    updated_coord_df[barcode, "imagecol"] <- colMapping[[as.character(newColValue)]]
  }
  
  # Create a new sample object with the updated coordinates
  updated_sample <- muestra1
  new_names <- c()
  
  for (i in 1:length(muestra1@meta.data$name)) {
    new_names <- c(new_names, paste0(muestra1@meta.data$name[i], "_Yinvert"))
  }
  
  updated_sample@meta.data$name <- new_names
  updated_sample@images$slice1@coordinates <- updated_coord_df
  
  # Flip the image data along the X-axis
  flipped_image <- array(NA, dim = dim(muestra1@images$slice1@image))
  
  for (i in 1:dim(muestra1@images$slice1@image)[3]) {
    flipped_image[,,i] <- muestra1@images$slice1@image[, ncol(muestra1@images$slice1@image):1, i]
  }
  
  updated_sample@images$slice1@image <- flipped_image
  
  return(updated_sample)
}

#' inversionEjeXY(muestra1)
#' 
#' This function flips the coordinates of a sample along both the X- and Y-axes 
#' (newrow = 77 - row, newcol = 127 - col). It also flips the image data along both axes.
#' 
#' @param muestra1 A sample object containing the image and coordinates.
#' @return A new sample object with updated coordinates and flipped image data.
#' @export
inversionEjeXY <- function(muestra1) {
  # Save the dataframe of coordinates from the sample
  coord_df <- muestra1@images$slice1@coordinates
  
  # Sort the dataframe by row and column
  sorted_coord_df <- coord_df[order(coord_df$row, coord_df$col),]
  
  # Create lists to store image row/column mappings for each row/col
  rowMapping <- list()
  colMapping <- list()
  
  # Loop through each barcode to store row and column mappings
  for (barcode in rownames(sorted_coord_df)) {
    rowValue       <- as.character(sorted_coord_df[barcode, "row"])
    imageRowValue  <- sorted_coord_df[barcode, "imagerow"]
    
    colValue       <- as.character(sorted_coord_df[barcode, "col"])
    imageColValue  <- sorted_coord_df[barcode, "imagecol"]
    
    if (!(imageRowValue %in% rowMapping)) {
      rowMapping[[rowValue]] <- imageRowValue
    }
    
    if (!(imageColValue %in% colMapping)) {
      colMapping[[colValue]] <- imageColValue
    }
  }
  
  # Create a copy of the dataframe
  updated_coord_df <- coord_df
  
  # Update the row, column and image row/col values
  for (barcode in rownames(coord_df)) {
    rowValue <- coord_df[barcode, "row"]
    colValue <- coord_df[barcode, "col"]
    
    # Flip the row and column values
    newRowValue <- 77   - rowValue
    newColValue <- 127  - colValue
    
    # Update the new row/col and image row/col in the dataframe
    updated_coord_df[barcode, "row"]      <- newRowValue
    updated_coord_df[barcode, "col"]      <- newColValue
    updated_coord_df[barcode, "imagerow"] <- rowMapping[[as.character(newRowValue)]]
    updated_coord_df[barcode, "imagecol"] <- colMapping[[as.character(newColValue)]]
  }
  
  # Create a new sample object with the updated coordinates
  updated_sample <- muestra1
  new_names <- c()
  for (i in 1:length(muestra1@meta.data$name)) {
    new_names <- c(new_names, paste0(muestra1@meta.data$name[i], "_XYinvert"))
  }
  updated_sample@meta.data$name <- new_names
  updated_sample@images$slice1@coordinates <- updated_coord_df
  
  # Flip the image data along both axes
  flipped_image <- array(NA, dim = dim(muestra1@images$slice1@image))
  for (i in 1:dim(muestra1@images$slice1@image)[3]) {
    # first invert Y (rows), then invert X (cols)
    flipped_image[,,i] <- muestra1@images$slice1@image[
      nrow(muestra1@images$slice1@image):1,
      ncol(muestra1@images$slice1@image):1,
      i
    ]
  }
  updated_sample@images$slice1@image <- flipped_image
  
  return(updated_sample)
}
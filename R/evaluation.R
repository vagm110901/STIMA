#' hex2rgb_table(objetHex)
#' 
#' Convert hexadecimal color data into an RGB array with dimensions (x, y, 4).
#' The function organizes color values (R, G, B, Alpha) into a 3D array.
#' 
#' @param objetHex A matrix or vector of hexadecimal color strings.  
#' @return rgba_array A 3D array of integers with dimensions (rows, columns, 4), containing RGBA color channels.
#' @import SpatialPack
#' @import imager
hex2rgb_table <- function(objetHex) {
  rgb_values <- col2rgb(objetHex, alpha = TRUE)  # Convert hex to RGB and alpha
  
  # Create an empty array for storing RGBA values and populate with color channels.
  x <- ncol(objetHex)
  y <- nrow(objetHex)
  rgba_array <- array(0, dim = c(y, x, 4))
  
  # Assign RGB and alpha values to the array for each color channel.
  rgba_array[,,1] <- matrix(rgb_values[1, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array[,,2] <- matrix(rgb_values[2, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array[,,3] <- matrix(rgb_values[3, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array[,,4] <- matrix(rgb_values[4, ], nrow = y, ncol = x, byrow = TRUE)
  
  return(rgba_array)
}

#' hex2rgb_table_norm(objetHex)
#' 
#' Convert hexadecimal color data into an RGB array with dimensions (x, y, 4), normalizes RGB values to a 0-1 range.
#' The function organizes color values (R, G, B, Alpha) into a 3D array.
#' 
#' @param objetHex A matrix or vector of hexadecimal color strings.  
#' @return rgba_array A 3D array of integers with dimensions (rows, columns, 4), containing RGBA color channels.
#' @import SpatialPack
#' @import imager
hex2rgb_table_norm <- function(objetHex) {
  rgb_values <- col2rgb(objetHex, alpha = TRUE) / 255
  
  # Prepare array to store normalized RGBA values.
  x <- ncol(objetHex)
  y <- nrow(objetHex)
  rgba_array <- array(0, dim = c(y, x, 4))
  
  # Populate array with normalized color channel data.
  rgba_array_norm[,,1] <- matrix(rgb_values[1, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array_norm[,,2] <- matrix(rgb_values[2, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array_norm[,,3] <- matrix(rgb_values[3, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array_norm[,,4] <- matrix(rgb_values[4, ], nrow = y, ncol = x, byrow = TRUE)
  
  return(rgba_array_norm)
}


#' rgb_table_norm2rgb_table(objetRGB)
#' 
#' Generates the 0-255 RGB values for an image with RGB normalized values
#' 
#' @param objetRGB A numeric array of RGB values normalized between 0 and 1.
#' @return An integer array with RGB values scaled to [0, 255].
#' @import SpatialPack
rgb_table_norm2rgb_table <- function(objetRGB) {
  objetRGB <- round(objetRGB * 255)
  storage.mode(objetRGB) <- "integer"
  
  return(objetRGB)
}


#' mse(image1, image2)
#' 
#' Calculate Mean Squared Error (MSE) between two images by comparing pixel values.
#' 
#' @param image1 A 3D array representing the first image (height x width x channels).
#' @param image2 A 3D array representing the second image (same dimensions as image1).
#' @return Numeric value representing the average MSE across RGB channels.
#' @import SpatialPack
#' @import imager
mse <- function(image1, image2) {
  sum_error <- 0  # Initialize total error
  
  # Loop through each color channel to calculate MSE.
  for (channel in seq_along(image1[1,1,]) - 1) {
    error <- sum((image1[,,channel] - image2[,,channel])^2)
    error <- error / (nrow(image1) * ncol(image2))
    sum_error <- sum_error + error
  }
  
  # Average the total error across color channels.
  total_error <- sum_error / 3
  return(total_error)
}

#' mseGS(image1, image2)
#' 
#' Compute the grayscale MSE between two images by first converting them to grayscale.
#' 
#' @param image1 A 3D array representing the first RGB image.
#' @param image2 A 3D array representing the second RGB image.
#' @return Numeric value representing the grayscale MSE between the two images.
#' @import SpatialPack
#' @import imager
mseGS <- function(image1, image2) {
  im1 <- RGB2gray(image1)  # Convert image 1 to grayscale
  im2 <- RGB2gray(image2)  # Convert image 2 to grayscale
  
  # Calculate grayscale MSE by summing squared differences.
  error <- sum((im1[,] - im2[,])^2) / (nrow(im1) * ncol(im2))
  return(error)
}

#' EuclDist(listaCoordenadas, nIm)
#' 
#' Calculate Euclidean distance between two sets of coordinates.
#' 
#' @param listaCoordenadas A list where each element contains a list with elements $x and $y representing coordinates.
#' @param nIm A list or vector of length 2 with indices of the images to compare.
#' @return Numeric value representing the Euclidean distance between the two coordinate sets.
#' @import SpatialPack
EuclDist <- function(listaCoordenadas, nIm) {
  # Retrieve coordinates for images i and j
  i <- nIm[[1]]
  j <- nIm[[2]]
  xi <- listaCoordenadas[[i]]$x
  yi <- listaCoordenadas[[i]]$y
  xj <- listaCoordenadas[[j]]$x
  yj <- listaCoordenadas[[j]]$y
  
  # Calculate the average Euclidean distance across all points.
  distance <- sqrt((xj - xi)^2 + (yj - yi)^2)
  return(distance)
}

#' evalAlign(image1, image2, listaCoordenadas, nIm)
#' 
#' Evaluate alignment between two images using various metrics (MSE, SSIM, etc.).
#' 
#' @param image1 An image array (RGB or hex) to be compared.
#' @param image2 An image array (RGB or hex) to be compared.
#' @param listaCoordenadas A list containing coordinates for Euclidean distance calculation.
#' @param nIm A list with indices of images in `listaCoordenadas` to compare.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{mse_value}: Mean Squared Error between the two images.
#'   \item \code{mse_gray_value}: Grayscale Mean Squared Error.
#'   \item \code{rmse_value}: Root Mean Squared Error.
#'   \item \code{ssim_value}: Structural Similarity Index.
#'   \item \code{Eucl_value}: Euclidean distance between the coordinate sets.
#' }
#'
#' @import SpatialPack
#' @import imager
evalAlign <- function(image1, image2, listaCoordenadas, nIm) {
  # Ensure images have the same dimensions, converting to RGB if necessary.
  if (!all(dim(image1) == dim(image2))) stop("Images must have the same dimensions.")
  if (is.na(dim(image1)[3])) image1 <- hex2rgb_table(image1) 
  if (is.na(dim(image2)[3])) image2 <- hex2rgb_table(image2)
  if (min(image1) >= 0 && max(image1) <= 1) image1 <- rgb_table_norm2rgb_table(image1)
  if (min(image2) >= 0 && max(image2) <= 1) image2 <- rgb_table_norm2rgb_table(image2)
  
  # Calculate alignment metrics: MSE, grayscale MSE, RMSE, SSIM, and Euclidean distance.
  # Calculate the MSE and RMSE
  mse_value <- mse(image1, image2)
  mse_gray_value <- mseGS(image1, image2)
  rmse_value <- sqrt(mse_value)
  # Calculate the SSIM
  im1 <- RGB2gray(image1)
  im2 <- RGB2gray(image2)
  ssim_value <- SpatialPack::SSIM(im1, im2)[["SSIM"]]
  # Calculate the Euclidean distance
  Eucl_value <- EuclDist(listaCoordenadas, nIm)
  
  # Compile metrics into a solution list.
  solution <- list(mse_value = mse_value, mse_gray_value = mse_gray_value,
                   rmse_value = rmse_value, ssim_value = ssim_value,
                   Eucl_value = Eucl_value)
  return(solution)
}

#' controlAlign(image1)
#' 
#' Evaluate alignment metrics for a reference image with various controls.
#' 
#' @param image1 An image array (RGB or hex) used as reference.
#'
#' @return A list with the following control results:
#' \itemize{
#'   \item \code{Positive control solution}: Metrics when image is compared to itself.
#'   \item \code{Negative control solution}: Metrics when image is compared to a random negative control image.
#'   \item \code{Movement control solution}: Metrics when image is compared to a shifted version.
#' }
#'
#' @import SpatialPack
#' @import imager
#' @import jpeg
controlAlign <- function(image1) { 
  # If needed, convert the image to RGB format.
  if (is.na(dim(image1)[3])) image1 <- hex2rgb_table(image1)
  if (min(image1) >= 0 && max(image1) <= 1) image1 <- rgb_table_norm2rgb_table(image1)
  
  # Define the central region for analysis.
  x1 <- round(dim(image1)[1] / 2)
  y1 <- round(dim(image1)[2] / 2)
  image1rec <- image1[(x1-200):(x1+200), (y1-200):(y1+200), ]
  
  # Calculate alignment metrics for positive control (image vs itself).
  mse_value <- mse(image1rec, image1rec)
  mse_gray_value <- mseGS(image1rec, image1rec)
  rmse_value <- sqrt(mse_value)
  im1 <- RGB2gray(image1rec)
  ssim_value <- SSIM(im1, im1)[["SSIM"]]
  positive_solution <- list(mse_value = mse_value, mse_gray_value = mse_gray_value,
                            rmse_value = rmse_value, ssim_value = ssim_value)
  
  # Calculate alignment metrics for a negative control (image vs random image).
  image2 <- readJPEG(system.file("extdata", "image_negative.JPG", package = "STIMA"))
  if (is.na(dim(image2)[3])) {image2 <- hex2rgb_table(image2)}
  
  image2 <- image2[1:dim(image1rec)[[1]], 1:dim(image1rec)[[2]], ]
  mse_value <- mse(image1rec, image2)
  mse_gray_value <- mseGS(image1rec, image2)
  rmse_value <- sqrt(mse_value)
  im2 <- RGB2gray(image2)
  ssim_value <- SSIM(im1, im2)[["SSIM"]]
  negative_solution <- list(mse_value = mse_value, mse_gray_value = mse_gray_value,
                            rmse_value = rmse_value, ssim_value = ssim_value)
  
  # Calculate alignment metrics for movement control (image vs shifted version 10 positions to the rigth).
  image3 <- image1[(x1-190):(x1+210), (y1-200):(y1+200), ]
  mse_value <- mse(image1rec, image3)
  mse_gray_value <- mseGS(image1rec, image3)
  rmse_value <- sqrt(mse_value)
  im3 <- RGB2gray(image3)
  ssim_value <- SSIM(im1, im3)[["SSIM"]]
  movement_solution <- list(mse_value = mse_value, mse_gray_value = mse_gray_value,
                            rmse_value = rmse_value, ssim_value = ssim_value)
  
  # Compile all control solutions into one list.
  solution <- list("Positive control solution" = positive_solution,
                   "Negative control solution" = negative_solution,
                   "Movement control solution" = movement_solution)
  return(solution)
}

#' evaluationComplete(listaCoordenadasNEW, listaCoordenadas, 
#'                    listaRawImages, listaTransImages, 
#'                    patientType = c('unique','multiple'))
#' 
#' Evaluate the alignment of images using various metrics (MSE, SSIM, etc.) for both raw and transformed images.
#'
#' The evaluation covers:
#' - Control alignment for transformed images
#' - Raw images comparisons using largest common square submatrix, same region & same point, and common region & different point
#' - Transformed images comparisons under the same scenarios as raw images
#'  
#' @param listaCoordenadasNEW List of coordinate data for transformed images.
#' @param listaCoordenadas List of coordinate data for raw images.
#' @param listaRawImages List of raw image arrays (typically 3D arrays).
#' @param listaTransImages List of transformed image arrays.
#' @param patientType Character specifying the patient type, affecting evaluation region size. 
#'   Choices are \code{"unique"} (region ~200x200 pixels) or \code{"multiple"} (region ~310x200 pixels).
#'   Default is \code{c("unique", "multiple")}, with \code{"unique"} selected by default.
#'
#' @return A nested list with alignment evaluation metrics for both raw and transformed images, organized as:
#' \describe{
#'   \item{control}{List of control alignment results for transformed images.}
#'   \item{original}{List of raw image comparisons including:
#'     \itemize{
#'       \item \code{max_pixels}: Evaluations on largest common square submatrix.
#'       \item \code{sameRegion_samePoint}: Evaluations on same region and same reference point.
#'       \item \code{commonRegion_differentPoint}: Evaluations on common region but different reference points.
#'     }
#'   }
#'   \item{transformed}{Same structure as \code{original} but for transformed images.}
#' }
#'
#' @import SpatialPack
#' @import imager
evaluationComplete <- function(listaCoordenadasNEW, listaCoordenadas,
                       listaRawImages, listaTransImages, patientType = c('unique','multiple')) {

  patientType <- match.arg(patientType)
  
  # Initialize an evaluation list
  Evaluation <- list()
  # Define the central region for analysis.
  x1 <- round(dim(listaRawImages[[1]])[1] / 2)
  y1 <- round(dim(listaRawImages[[1]])[2] / 2)
  
  # Create a list to store control alignment parameters
  control <- list()
  
  # Loop through the new coordinates and evaluate the alignment
  for (i in 1:length(listaCoordenadasNEW)) {
    control[[i]] <- controlAlign(listaTransImages[[i]])
    print(control)
  }
  
  # Store the control evaluations in the Evaluation list
  Evaluation$control <- control
  
  # Evaluate raw images (shared pixels)
  for (i in 1:length(listaCoordenadas)) {
    if (i != length(listaCoordenadas)) {
      for (j in i:length(listaCoordenadas))  {
        if (i != j) {
          dim1 <- dim(listaRawImages[[i]])
          dim2 <- dim(listaRawImages[[j]])
          
          common_rows <- min(dim1[1], dim2[1])
          common_cols <- min(dim1[2], dim2[2])
          
          img1_crop <- listaRawImages[[i]][1:common_rows, 1:common_cols, ]
          img2_crop <- listaRawImages[[j]][1:common_rows, 1:common_cols, ]
          
          mask1 <- apply(!is.na(img1_crop) & img1_crop != 0 & img1_crop != 1, c(1,2), any)
          mask2 <- apply(!is.na(img2_crop) & img2_crop != 0 & img2_crop != 1, c(1,2), any)
          mask <- mask1 & mask2
          
          # Algorith for largest square submatrix in `mask`
          nr <- nrow(mask); nc <- ncol(mask)
          sq <- matrix(0, nr, nc)
          maxSize <- 0; max_i <- 0; max_j <- 0
          for (r in 1:nr) {
            for (c in 1:nc) {
              if (!mask[r, c]) next
              if (r == 1 || c == 1) {
                sq[r, c] <- 1
              } else {
                sq[r, c] <- min(sq[r-1, c], sq[r, c-1], sq[r-1, c-1]) + 1
              }
              if (sq[r, c] > maxSize) {
                maxSize <- sq[r, c]
                max_i <- r
                max_j <- c
              }
            }
          }
          
          # Coordinates from the square
          rows_sq <- (max_i - maxSize + 1):max_i
          cols_sq <- (max_j - maxSize + 1):max_j
          
          # Evaluate the alignment and store parameters
          parameters <- tryCatch({
            evalAlign(
              listaRawImages[[i]][rows_sq, cols_sq, ],
              listaRawImages[[j]][rows_sq, cols_sq, ],
              listaCoordenadas, c(i,j)) 
          }, error = function(e) {
            message(sprintf("Error when comparing %d-%d: %s", i, j, e$message))
            return(NULL)
          })
          print(parameters)
          imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
          Evaluation$original$max_pixels[[imagescompare]] <- parameters
          
        }}}}
  
  
  
  # Evaluate raw images (same region, same point)
  for (i in 1:length(listaCoordenadas)) {
    if (i != length(listaCoordenadas)) {
      for (j in i:length(listaCoordenadas))  {
        if (i != j) {
          x1 <- round(dim(listaRawImages[[i]])[1] / 2)
          y1 <- round(dim(listaRawImages[[i]])[2] / 2)
          dim1 <- dim(listaRawImages[[i]])
          dim2 <- dim(listaRawImages[[j]])
          
          if (patientType == 'unique') { 
            x1_min <- max(0, x1 - 200)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1 + 200)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            
          } else if (patientType == 'multiple') {
            x1_min <- max(0, x1 - 310)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
          }
          print(paste0("Evaluation of the alignment between images ", 
                as.character(i), " and ", 
                as.character(j), " comparing the same region without selecting a different point."))
          
          # Evaluate the alignment and store parameters
          parameters <- 
            evalAlign(
              listaRawImages[[i]][coordRange[[1]],coordRange[[2]],], 
              listaRawImages[[j]][coordRange[[1]],coordRange[[2]],], 
              listaCoordenadas, c(i,j))
          print(parameters)
          imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
          Evaluation$original$sameRegion_samePoint[[imagescompare]] <- parameters
          
        }}}}
  
  # Evaluate raw images (common region, different point)
  for ( i in 1:length(listaCoordenadas)) {
    if (i != length(listaCoordenadas)) {
      for ( j in i:length(listaCoordenadas))  {
        if (i != j) {
          x1 <- round(dim(listaRawImages[[i]])[1] / 2)
          y1 <- round(dim(listaRawImages[[i]])[2] / 2)
          x2 <- round(dim(listaRawImages[[j]])[1] / 2)
          y2 <- round(dim(listaRawImages[[j]])[2] / 2)
          dim1 <- dim(listaRawImages[[i]])
          dim2 <- dim(listaRawImages[[j]])
          
          if (patientType == 'unique') { 
            x1_min <- max(0, x1 - 200)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1 + 200)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            x2_min <- max(0, x2 - 200)
            y2_min <- max(0, y2 - 200)
            x2_max <- min(dim1[1], dim2[1], x2 + 200)
            y2_max <- min(dim1[2], dim2[2], y2 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            coordRange2 <- list(x2_min:x2_max, y2_min:y2_max)
            
          } else if (patientType == 'multiple') {
            x1_min <- max(0, x1 - 310)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            x2_min <- max(0, x2 - 310)
            y2_min <- max(0, y2 - 200)
            x2_max <- min(dim1[1], dim2[1], x2)
            y2_max <- min(dim1[2], dim2[2], y2 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            coordRange2 <- list(x2_min:x2_max, y2_min:y2_max)
          }
          print(paste0("Evaluation of the alignment between images ", 
                as.character(i), " and ", 
                as.character(j), " selecting an area from a reference point in each image."))
          
          # Evaluate the alignment and store parameters
          parameters <- tryCatch({
            evalAlign(
              listaRawImages[[i]][coordRange[[1]],coordRange[[2]],],
              listaRawImages[[j]][coordRange2[[1]],coordRange2[[2]],],
              listaCoordenadas, c(i,j))
          }, error = function(e) {
            message(sprintf("Error when comparing %d-%d: %s", i, j, e$message))
            return(NULL)
            })
          print(parameters)
          imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
          Evaluation$original$commonRegion_differentPoint[[imagescompare]] <- parameters
        }}}}
  
  
  # Evaluate transformed images (shared pixels)
  for (i in 1:length(listaCoordenadasNEW)) {
    if (i != length(listaCoordenadasNEW)) {
      for (j in i:length(listaCoordenadasNEW))  {
        if (i != j) {
          dim1 <- dim(listaTransImages[[i]])
          dim2 <- dim(listaTransImages[[j]])
          
          common_rows <- min(dim1[1], dim2[1])
          common_cols <- min(dim1[2], dim2[2])
          
          img1_crop <- listaTransImages[[i]][1:common_rows, 1:common_cols, ]
          img2_crop <- listaTransImages[[j]][1:common_rows, 1:common_cols, ]
          
          mask1 <- apply(!is.na(img1_crop) & img1_crop != 0 & img1_crop != 1, c(1,2), any)
          mask2 <- apply(!is.na(img2_crop) & img2_crop != 0 & img2_crop != 1, c(1,2), any)
          mask <- mask1 & mask2
          
          # Algorith for largest square submatrix in `mask`
          nr <- nrow(mask); nc <- ncol(mask)
          sq <- matrix(0, nr, nc)
          maxSize <- 0; max_i <- 0; max_j <- 0
          for (r in 1:nr) {
            for (c in 1:nc) {
              if (!mask[r, c]) next
              if (r == 1 || c == 1) {
                sq[r, c] <- 1
              } else {
                sq[r, c] <- min(sq[r-1, c], sq[r, c-1], sq[r-1, c-1]) + 1
              }
              if (sq[r, c] > maxSize) {
                maxSize <- sq[r, c]
                max_i <- r
                max_j <- c
              }
            }
          }
          
          # Coordinates from the square
          rows_sq <- (max_i - maxSize + 1):max_i
          cols_sq <- (max_j - maxSize + 1):max_j
          
          # Evaluate the alignment and store parameters
          parameters <- tryCatch({
            evalAlign(
              listaTransImages[[i]][rows_sq, cols_sq, ],
              listaTransImages[[j]][rows_sq, cols_sq, ],
              listaCoordenadasNEW, c(i,j)) 
          }, error = function(e) {
            message(sprintf("Error when comapring %d-%d: %s", i, j, e$message))
            return(NULL)
          })
          print(parameters)
          imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
          Evaluation$transformed$max_pixels[[imagescompare]] <- parameters
          
        }}}}
  
  
  # Evaluate transformed images (same region, same point)
  for ( i in 1:length(listaCoordenadasNEW)) {
    if (i != length(listaCoordenadasNEW)) {
      for ( j in i:length(listaCoordenadasNEW))  {
        if (i != j) {
          x1 <- round(dim(listaTransImages[[i]])[1] / 2)
          y1 <- round(dim(listaTransImages[[i]])[2] / 2)
          x1 <- listaCoordenadasNEW[[i]]$x[[5]]
          y1 <- listaCoordenadasNEW[[i]]$y[[5]]
          dim1 <- dim(listaTransImages[[i]])
          dim2 <- dim(listaTransImages[[j]])
          
          if (patientType == 'unique') { 
            x1_min <- max(0, x1 - 200)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1 + 200)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            
          } else if (patientType == 'multiple') {
            x1_min <- max(0, x1 - 310)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
          }
          print(paste0("Evaluation of the alignment between images ", 
                as.character(i), " and ", 
                as.character(j), " comparing the same region without selecting a different point."))
          
          # Evaluate the alignment and store parameters
          parameters <- 
            evalAlign(
              listaTransImages[[i]][coordRange[[1]],coordRange[[2]],],
              listaTransImages[[j]][coordRange[[1]],coordRange[[2]],],
              listaCoordenadasNEW, c(i,j))
          print(parameters)
          imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
          Evaluation$transformed$sameRegion_samePoint[[imagescompare]] <- parameters
          
        }}}}
  
  # Evaluate transformed images (common region, different point)
  for ( i in 1:length(listaCoordenadasNEW)) {
    if (i != length(listaCoordenadasNEW)) {
      for ( j in i:length(listaCoordenadasNEW))  {
        if (i != j) {
          x1 <- round(dim(listaTransImages[[i]])[1] / 2)
          y1 <- round(dim(listaTransImages[[i]])[2] / 2)
          x2 <- round(dim(listaTransImages[[j]])[1] / 2)
          y2 <- round(dim(listaTransImages[[j]])[2] / 2)
          dim1 <- dim(listaTransImages[[i]])
          dim2 <- dim(listaTransImages[[j]])
          
          if (patientType == 'unique') { 
            x1_min <- max(0, x1 - 200)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1 + 200)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            x2_min <- max(0, x2 - 200)
            y2_min <- max(0, y2 - 200)
            x2_max <- min(dim1[1], dim2[1], x2 + 200)
            y2_max <- min(dim1[2], dim2[2], y2 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            coordRange2 <- list(x2_min:x2_max, y2_min:y2_max)
            
          } else if (patientType == 'multiple') {
            x1_min <- max(0, x1 - 310)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            x2_min <- max(0, x2 - 310)
            y2_min <- max(0, y2 - 200)
            x2_max <- min(dim1[1], dim2[1], x2)
            y2_max <- min(dim1[2], dim2[2], y2 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            coordRange2 <- list(x2_min:x2_max, y2_min:y2_max)
          }
          print(paste0("Evaluation of the alignment between images ", 
                as.character(i), " and ", 
                as.character(j), " selecting an area from a reference point in each image."))
          
          # Evaluate the alignment and store parameters
          parameters <- tryCatch({
            evalAlign(
              listaTransImages[[i]][coordRange[[1]],coordRange[[2]],],
              listaTransImages[[j]][coordRange2[[1]],coordRange2[[2]],],
              listaCoordenadasNEW, c(i,j)) 
          }, error = function(e) {
            message(sprintf("Error when comparing %d-%d: %s", i, j, e$message))
            return(NULL)
          })
          print(parameters)
          imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
          Evaluation$transformed$commonRegion_differentPoint[[imagescompare]] <- parameters
        }}}}
  
  return(Evaluation)
}


#' calculateEvaluation(objeto.seurat, mode = c("GTEM", "procrustes", "RVSSimageJ"), 
#'                     listaCoordenadasNEW, listaCoordenadas, 
#'                     patientType = c('unique','multiple'))
#' 
#' Evaluate the alignment of images using various metrics (MSE, SSIM, etc.) for both raw and transformed images.
#' 
#' @param objeto.seurat Seurat object containing images and coordinates.
#' @param mode Evaluation mode: one of "GTEM", "procrustes", or "RVSSimageJ".
#' @param listaCoordenadasNEW List of new coordinates (optional, used mainly if mode != "RVSSimageJ").
#' @param listaCoordenadas List of original coordinates (optional, used mainly if mode != "RVSSimageJ").
#' @param patientType Patient type, affecting region size ("unique" or "multiple").
#' @return Evaluation results as a list.
#' @import SpatialPack
#' @import imager
#' @export 
calculateEvaluation <- function(objeto.seurat, mode = c("GTEM", "procrustes", "RVSSimageJ"), 
                                listaCoordenadasNEW = NULL, listaCoordenadas = NULL, 
                                patientType = c('unique','multiple')) {

  patientType <- match.arg(patientType)
  mode <- match.arg(mode)

  if (!dir.exists("./results/")) {
    dir.create("./results/", recursive = TRUE)
  }
  saveDir <- "./results/"
  
  if (mode == "RVSSimageJ") {
    listaCoordenadas <- list()
    listaCoordenadasNEW <- list()
    for (i in 1:length(objeto.seurat@tools$Staffli@rasterlists$raw)) {
      # Select the Original Coordinates and the Aligned Coordinates
      # Select and round coordinates for the current image
      coordenadas <- selectCoord(objeto.seurat@tools$Staffli@rasterlists$raw[[i]])
      NEWcoordenadas <- selectCoord(objeto.seurat@tools$Staffli@rasterlists$transformed[[i]])
      for (j in seq_along(coordenadas)) { 
        coordenadas[[j]] <- round(coordenadas[[j]])
        NEWcoordenadas[[j]] <- round(NEWcoordenadas[[j]])
      }
      listaCoordenadas[[i]] <- coordenadas
      listaCoordenadasNEW[[i]] <- NEWcoordenadas
    }
  }

  # Load the images
  listaRawImages <- list()
  listaTransImages <- list()
  for (i in 1:length(listaCoordenadasNEW)) {
    listaRawImages[[i]] <- objeto.seurat@images[[paste0("slice1.",i)]]@image
    listaTransImages[[i]] <- objeto.seurat@images[[paste0("slice",i)]]@image
  }
  
  # Evaluate the alignment
  Evaluation <- evaluationComplete(listaCoordenadasNEW, listaCoordenadas,
                                   listaRawImages, listaTransImages, patientType)
  
  saveRDS(Evaluation, paste0(saveDir, "evaluation_merge_",mode,".rds"))

  return(Evaluation)
}

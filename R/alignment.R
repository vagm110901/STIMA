#' solveCoord(coord1, coord2, xmax, ymax)
#'
#' This function calculates the raw transformation parameters (`cos(angle)`, `sin(angle)`, `dx`, and `dy`)
#' required to align two sets of coordinates. It does this by solving an overdetermined system of equations,
#' optimizing the transformation to minimize the difference between the transformed coordinates.
#'
#' @param coord1 A list of coordinates representing the target points. 
#'   Should be a list with components `x` and `y`, each containing numeric vectors of the same length.
#' @param coord2 A list of coordinates representing the source points.
#'   Should be a list with components `x` and `y`, each containing numeric vectors of the same length.
#' @param xmax Numeric. The maximum value for the x-dimension, used for centering the coordinates.
#' @param ymax Numeric. The maximum value for the y-dimension, used for centering the coordinates.
#' 
#' @return A named numeric vector with the estimated transformation parameters:
#' \describe{
#'   \item{coseno}{The cosine of the rotation angle.}
#'   \item{seno}{The sine of the rotation angle.}
#'   \item{dx}{The x-translation component.}
#'   \item{dy}{The y-translation component.}
#' }
#'
#' @export
solveCoord <- function(coord1, coord2, xmax, ymax) {
  cx <- xmax / 2
  cy <- ymax / 2
  
  # Adjust coordinates to be centered based on xmax and ymax
  coord2MOD <- coord2
  for (j in seq_along(coord2$x)) {
    coord2MOD$x[[j]] <- coord2$x[[j]] - cx
    coord2MOD$y[[j]] <- coord2$y[[j]] - cy
  }
  
  coord1MOD <- coord1
  for (j in seq_along(coord2$x)) {
    coord1MOD$x[[j]] <- coord1$x[[j]] - cx
    coord1MOD$y[[j]] <- coord1$y[[j]] - cy
  }
  
  # Optimization of an overdetermined equation system
  fn <- function(x) {
    y <- numeric(0)
    for (i in seq_along(coord1$x)) {
      y[length(y) + 1] <- coord2MOD$x[[i]]*x[1] + coord2MOD$y[[i]]*x[2] + 1*x[3] + 0*x[4] - coord1MOD$x[[i]]
      y[length(y) + 1] <- coord2MOD$y[[i]]*x[1] - coord2MOD$x[[i]]*x[2] + 0*x[3] + 1*x[4] - coord1MOD$y[[i]]
    }
    y[length(y) + 1] <- (x[1])^2 + (x[2])^2 - 1  
    return(sum(y^2))
  }
  
  xstart <- c(cos(pi/4), sin(pi/4), 2, -2)
  solution <- optim(xstart, fn)
  solution <- solution$par
  names(solution) <- c('coseno', 'seno', 'dx', 'dy')
  
  return(solution)
}

#' calcParameters(solution, xmax, ymax)
#'
#' This function takes a set of raw transformation parameters and the dimensions of the image
#' to calculate the complete set of transformation parameters needed for functions in the **semla** package.
#'
#' @param solution A named numeric vector containing the raw transformation parameters, typically
#'   the result from the `solveCoord()` function. Must include at least:
#'   \describe{
#'     \item{coseno}{The cosine of the rotation angle.}
#'     \item{seno}{The sine of the rotation angle.}
#'     \item{dx}{The x-translation component.}
#'     \item{dy}{The y-translation component.}
#'     \item{e}{The scaling factor (optional, default is 1).}
#'     \item{mirrorx}{Mirror flag for the x-axis (optional, default is 0).}
#'     \item{mirrory}{Mirror flag for the y-axis (optional, default is 0).}
#'   }
#' @param xmax Numeric. The maximum value for the x-dimension, used to scale the translation.
#' @param ymax Numeric. The maximum value for the y-dimension, used to scale the translation.
#' @param mode Character. Alignment mode to use. One of "GTEM" or "procrustes".
#'   - "GTEM": Geometric transformations without scaling.
#'   - "procrustes": Uses Procrustes analysis with optional scaling.
#' 
#' @return A named numeric vector containing the calculated transformation parameters:
#' \describe{
#'   \item{angulo}{The rotation angle in degrees.}
#'   \item{trx}{The x-translation component, normalized by `xmax`.}
#'   \item{try}{The y-translation component, normalized by `ymax` (note the sign inversion).}
#'   \item{e}{The scaling factor.}
#'   \item{mirrorx}{Mirror flag for the x-axis.}
#'   \item{mirrory}{Mirror flag for the y-axis.}
#' }
#'
#' @export
calcParameters <- function(solution, xmax, ymax, mode = c("GTEM", "procrustes") ) {
  mode <- match.arg(mode)
  # Calculate angles and restrict cosine and sine values to the range [-1, 1]
  coseno <- solution[["coseno"]]
  seno <- solution[["seno"]]
  if (coseno > 1) {coseno <- 1} else {if (coseno < -1) {coseno <- -1}}
  if (seno > 1) {seno <- 1} else {if (seno < -1) {seno <- -1}}
  
  tangente <- seno / coseno
  angle <- atan2(seno, coseno) 
  angulo <- angle * (180/pi)
  
  # Calculate x and y translation parameters
  trx <-  solution[["dx"]] / xmax
  try <- - solution[["dy"]] / ymax
  
  trx <- max(min(trx, 1), -1)
  try <- max(min(try, 1), -1)
  
  e <- solution[["e"]]  # Scaling factor

  mirrorx <- solution[["mirrorx"]]
  mirrory <- solution[["mirrory"]]
  
  valores <- c(angulo, trx, try, e, mirrorx, mirrory)
  names(valores) <- c('angulo', 'trx', 'try', 'e', 'mirrorx', 'mirrory')
  
  # Replace any NA values with zero
  for (i in seq_along(valores)) {
    if (is.na(valores[[i]])) {valores[[i]] <- 0}
  }
  
  return(valores)
}


#' calcNewCoord(coord, sol, xmax, ymax)
#' 
#' This function applies a series of geometric transformations (mirroring, centering, rotation, and translation)
#' to a set of original coordinates based on precomputed transformation parameters and image dimensions.
#' 
#' @param coord A list containing the original coordinates to be transformed. 
#'   Typically includes two elements:
#'   \describe{
#'     \item{x}{A numeric vector of x-coordinates.}
#'     \item{y}{A numeric vector of y-coordinates.}
#'   }
#' @param sol A named numeric vector containing the transformation parameters, usually generated
#'   by the `calcParameters()` function. Must include:
#'   \describe{
#'     \item{coseno}{The cosine of the rotation angle.}
#'     \item{seno}{The sine of the rotation angle.}
#'     \item{dx}{The x-translation component.}
#'     \item{dy}{The y-translation component.}
#'     \item{mirrorx}{Mirror flag for the x-axis (1 for mirrored, 0 otherwise).}
#'     \item{mirrory}{Mirror flag for the y-axis (1 for mirrored, 0 otherwise).}
#'   }
#' @param xmax Numeric. The maximum x-dimension of the image, used for centering.
#' @param ymax Numeric. The maximum y-dimension of the image, used for centering.
#' 
#' @return A list containing the transformed coordinates, with the same structure as the input `coord`:
#' \describe{
#'   \item{x}{A numeric vector of transformed x-coordinates.}
#'   \item{y}{A numeric vector of transformed y-coordinates.}
#' }
#' 
#' @export
calcNewCoord <- function(coord, sol, xmax, ymax) {
  cx <- xmax / 2
  cy <- ymax / 2
  
  coordMOD <- coord
  
  # Apply mirror transformations based on solution parameters
  if (sol[["mirrorx"]] == 1) {
    for (j in seq_along(coordMOD$x)) {
      coordMOD$x[[j]] <- xmax - coordMOD$x[[j]]
    }
  }
  
  if (sol[["mirrory"]] == 1) {
    for (j in seq_along(coordMOD$y)) {
      coordMOD$y[[j]] <- ymax - coordMOD$y[[j]]
    }
  }
  
  # Center coordinates
  for (j in seq_along(coordMOD$x)) {
    coordMOD$x[[j]] <- coordMOD$x[[j]] - cx
    coordMOD$y[[j]] <- coordMOD$y[[j]] - cy
  }
  
  # Apply rotation and translation transformations
  coordNEW <- coordMOD
  for (j in seq_along(coordMOD$x)) {
    coordNEW$x[[j]] <- round(coordMOD$x[[j]]*sol[["coseno"]] + coordMOD$y[[j]]*sol[["seno"]] + sol[["dx"]] + cx)
    coordNEW$y[[j]] <- round(coordMOD$y[[j]]*sol[["coseno"]] - coordMOD$x[[j]]*sol[["seno"]] + sol[["dy"]] + cy)
  }
  
  return(coordNEW)
}

#' calcScal(coord1, coord2)
#' 
#' Calculate the scale factor between two sets of coordinates.
#' 
#' This function computes the optimal scale factor that minimizes the difference
#' between two sets of coordinates (coord1 and coord2) using a Brent optimization method.
#'
#' @param coord1 A list with named components 'x' and 'y', representing the first set of coordinates.
#' @param coord2 A list with named components 'x' and 'y', representing the second set of coordinates.
#' @return A numeric value representing the optimal scale factor.
#' @export
calcScal <- function(coord1, coord2) {
  fn <- function(x){
    y <- numeric(0)
    for (i in seq_along(coord1$x)) {
      y[length(y) + 1] <- coord2$x[[i]]*x[1] - coord1$x[[i]]
      y[length(y) + 1] <- coord2$y[[i]]*x[1] - coord1$y[[i]]
    }
    return(sum(y^2))
  }
  
  xstart <- 1
  solution <- optim(xstart, fn, method = "Brent", lower = 0, upper = 3)
  solution <- solution$par
  names(solution) <- 'e'
  
  return(solution)
}

#' resultProcrustes(proc, mirrorx, mirrory, scale)
#' 
#' Extracts the key transformation parameters from a Procrustes analysis result,
#' including rotation, translation, reflection (mirroring), and scaling.
#' 
#' @param proc A list containing the results of a Procrustes analysis. Expected to have elements:
#'   - `R`: a rotation matrix,
#'   - `t`: a translation vector,
#'   - `d`: a scaling factor (distance).
#' @param mirrorx 0 or anothe value; whether to apply reflection along the x-axis.
#' @param mirrory 0 or anothe value; whether to apply reflection along the y-axis.
#' @param scale Logical; whether to apply scaling.
#' @return A named numeric vector with the following elements:
#'   - `coseno`: cosine of the rotation angle,
#'   - `seno`: sine of the rotation angle,
#'   - `dx`: translation along the x-axis,
#'   - `dy`: translation along the y-axis,
#'   - `mirrorx`: reflection flag on the x-axis,
#'   - `mirrory`: reflection flag on the y-axis,
#'   - `e`: scaling factor (1 if `scale` is FALSE).
#'
#' @export
resultProcrustes <- function(proc, mirrorx, mirrory, scale) {
  coseno <- proc$R[1,1]
  seno <- proc$R[2,1]
  if (scale == TRUE) {e <- proc$d} else if (scale == FALSE) {e <- 1} 
  
  dx <- proc$t[1]
  dy <- proc$t[2]
  
  solucion <- c(coseno, seno, dx, dy, mirrorx, mirrory, e)
  names(solucion) <- c('coseno', 'seno', 'dx', 'dy', 'mirrorx', 'mirrory', 'e')

  return(solucion)
}
  
#' STIMA(object, mode, scale)
#' 
#' @param object A Seurat object with merged slices/samples containing spatial and image data.
#' @param mode Character. Alignment mode to use. One of "GTEM", "procrustes", or "RVSSimageJ".
#'   - "GTEM": Geometric transformations without scaling.
#'   - "procrustes": Uses Procrustes analysis with optional scaling.
#'   - "RVSSimageJ": Uses ImageJ Register Virtual Stack Slices plugin for manual alignment.
#' @param scale Logical. Whether to compute scaling during alignment (only applies to "procrustes" and "RVSSimageJ" modes).
#'
#' @return A Seurat object updated with alignment transformations.
#' 
#' @return object.seurat aligned
#' @import semla
#' @import Seurat
#' @import tibble
#' @import reticulate
#' @export
STIMA <- function(object, mode = c("GTEM", "procrustes", "RVSSimageJ"), scale = FALSE) {
  mode <- match.arg(mode) 

  if (!dir.exists("./results/")) {
    dir.create("./results/", recursive = TRUE)
  }
  saveDir <- "./results/"

  x_dims <- sapply(object@images, function(img) { dim(img)[[1]] })

  min_image_height <- as.numeric(min(x_dims))

  object.semla <- UpdateSeuratForSemla(object)
  object.semla <- LoadImages(object.semla, image_height = min_image_height)

  #ImagePlot(object.semla)

  if (mode == "RVSSimageJ") {
    # Set up lists for alignment process and coordinate retrieval
    aligned <- list(object.semla@tools$Staffli@rasterlists$raw[[1]])
    original <- list(object.semla@tools$Staffli@rasterlists$raw[[1]])
    listaSoluciones <- list(NA)
    listaTransforms <- list(NA)
    
    # Extract column and row coordinates for sample ID 1 from Seurat object
    CoordinatesSemlaCol <- list(GetCoordinates(object.semla)[which(GetCoordinates(object.semla)["sampleID"] == 1), 2])
    CoordinatesSemlaRow <- list(GetCoordinates(object.semla)[which(GetCoordinates(object.semla)["sampleID"] == 1), 3])
    
    "
    This is the point where you would save the images from the RDS file into a source_dir folder, 
    create a target_dir folder, and then use ImageJ with the Register Virtual Stack Slices plugin for alignment.
    "
    
    # Define source and target directories for ImageJ alignment
    source_dir <- paste0("./original")
    target_dir <- paste0("./aligned")
    #ref_name <- "original"
    ref_name <- "tissue_lowres_image_1"
    
    # Instructions for folder setup before alignment
    cat(paste("\n\nYou should have created in this directory a folder called: ", source_dir,
              "\nAnd a folder called: ", target_dir,
              "\nIn the first one, you should have the images from de rds document.\n"))
    
    # Provide instructions for the manual alignment process in ImageJ
    cat(paste("\nIt is time to open ImageJ and do the alignment of the images.",
              "\nFirst, go to File>Open and select the different images to align.",
              "The images should have only one color channel.",
              "To do this you should go to Image>Color>Split Channels and then Image>Color>Merge Channels unselecting create composite.",
              "\nTo do the align go to Plugins>Registration>Register Virtual Stack Slices.", 
              "Remember you should save the transforms.",
              "\nThe transformation parameters or Transform files should be stored in the same folder as the result images.",
              "\nAnd finally you have to select the reference image.\n"))
    
    # Wait for user confirmation before proceeding
    answer <- "no"
    while (answer != 'yes') {
      cat(paste0("\nWhen you have all this completed, write: 'yes': "))
      answer <- readline()
    }
    if (answer == 'yes') {
      cat(paste0("\nThe alignment is done and the transformation parameters are saved."))
    }
    
    # Extract translation and rotation parameters from XML files in the target directory
    filesTarget <- list.files(target_dir)
    files_xml <- filesTarget[grep("\\.xml$", filesTarget)]
    
    Nimage <- 1
    for (xml in files_xml) {
      # Skip the reference image file and process each alignment XML file
      if (!startsWith(xml, ref_name)) {
        Nimage <- Nimage + 1
        lines <- readLines(paste0(target_dir, "/", xml))
        lines <- lines[startsWith(lines, "\t<iict_transform")]
        lines <- sapply(lines, function(line) substr(line, start = nchar("\\t<iict_transform"), stop = nchar(line) - 3))
        lines <- sapply(lines, function(line) strsplit(line, "\""))
        
        # Extract rotation and translation parameters based on transform class type
        for (line in lines) {
          classN <- which(startsWith(line, " class="))
          if (endsWith(line[[classN + 1]], "transform.RigidModel2D")) { 
            dataN <- which(startsWith(line, " data="))
            transformsR <- line[[dataN + 1]]
            transformsR <- strsplit(transformsR, " ")
          }
          if (endsWith(line[[classN + 1]], "transform.TranslationModel2D")) { 
            dataN <- which(startsWith(line, " data="))
            transformsT <- line[[dataN + 1]]
            transformsT <- strsplit(transformsT, " ")
          }
        }
        if (scale == FALSE) {
          # Save transformation parameters: rotation angle, dx, and dy
          transformsParams <- list()
          transformsParams["angle"] <- as.numeric(transformsR[[1]][1])  # radians
          transformsParams["dx"] <- as.numeric(transformsR[[1]][2]) + as.numeric(transformsT[[1]][1])  # without normalization (-1,1)
          transformsParams["dy"] <- as.numeric(transformsR[[1]][3]) + as.numeric(transformsT[[1]][2])  # without normalization (-1,1)
          listaSoluciones[[Nimage]] <- transformsParams
        } else if (scale == TRUE) {
          # Store the transformation parameters (angle, scale, and translation values)
          transformsParams <- list()
          transformsParams["s*cos"] <- as.numeric(transformsS[[1]][1])
          transformsParams["s*sin"] <- as.numeric(transformsS[[1]][2])  # radianes
          transformsParams["dx"] <- as.numeric(transformsS[[1]][3]) + as.numeric(transformsT[[1]][1])     # sin relativizar (-1,1)
          transformsParams["dy"] <- as.numeric(transformsS[[1]][4]) + as.numeric(transformsT[[1]][2])     # sin relativizar (-1,1)
          listaSoluciones[[Nimage]] <- transformsParams
        } 
      }
    }
  } else if (mode == "GTEM" | mode == "procrustes") {
    
    # Select and round the coordinates from the reference image (first sample)
    coordenadas1 <- selectCoord(object.semla@tools$Staffli@rasterlists$raw[[1]])
    for (j in seq_along(coordenadas1)) {
      coordenadas1[[j]] <- round(coordenadas1[[j]])
    }
    print(coordenadas1)  # Display rounded coordinates for the reference image
    matrixCoord1 <- matrix(data = unlist(coordenadas1), ncol = 2)

    # Initialize lists for storing image data, coordinates, transformations, and other parameters
    aligned <- list(object.semla@tools$Staffli@rasterlists$raw[[1]])
    original <- list(object.semla@tools$Staffli@rasterlists$raw[[1]])
    listaCoordenadas <- list(coordenadas1)
    listaSoluciones <- list(NA)
    listaOpciones <- list(NA)
    listaOpcionesCalc <- list(NA)
    listaValores <- list(NA)
    listaTransforms <- list(NA)
    listaCoordenadasNEW <- list(coordenadas1)

    # Extract column and row coordinates for the first sample
    CoordinatesSemlaCol <- list(GetCoordinates(object.semla)[which(GetCoordinates(object.semla)["sampleID"] == 1), 2])
    CoordinatesSemlaRow <- list(GetCoordinates(object.semla)[which(GetCoordinates(object.semla)["sampleID"] == 1), 3])
  
    # Iterate through each image starting from the second sample
    for (i in 2:length(object.semla@tools$Staffli@rasterlists$raw)) {
      # Select and round coordinates for the current image
      coordenadas2 <- selectCoord(object.semla@tools$Staffli@rasterlists$raw[[i]])
      for (j in seq_along(coordenadas2)) {
        coordenadas2[[j]] <- round(coordenadas2[[j]])
      }
      print(coordenadas2)  # Display rounded coordinates for the current image
      listaCoordenadas[[i]] <- coordenadas2
    
      # Define dimensions for the current image
      xmax2 <- nrow(object.semla@tools$Staffli@rasterlists$raw[[i]])
      ymax2 <- ncol(object.semla@tools$Staffli@rasterlists$raw[[i]])

      coordCalc <- list()

      if (mode == "GTEM") {

        # Solve transformations without mirroring (original orientation)
        solucionOrig <- solveCoord(coordenadas1, coordenadas2, xmax2, ymax2)
        solucionOrig[["mirrorx"]] <- 0
        solucionOrig[["mirrory"]] <- 0
        solucionOrig[["e"]] <- 1
        coordCalc[["solucionOrig"]] <- calcNewCoord(coordenadas2, solucionOrig, xmax2, ymax2)
        
    
        # Solve transformations with mirroring on the x-axis
        coordenadas2X <- coordenadas2
        for (j in seq_along(coordenadas2X$x)) {
          coordenadas2X$x[[j]] <- xmax2 - coordenadas2X$x[[j]]
        }
        solucionMirrorX <- solveCoord(coordenadas1, coordenadas2X, xmax2, ymax2)
        solucionMirrorX[["mirrorx"]] <- 10
        solucionMirrorX[["mirrory"]] <- 0
        solucionMirrorX[["e"]] <- 1
        coordCalc[["solucionMirrorX"]] <- calcNewCoord(coordenadas2, solucionMirrorX, xmax2, ymax2)
    
        # Solve transformations with mirroring on the y-axis
        coordenadas2Y <- coordenadas2
        for (j in seq_along(coordenadas2Y$y)) {
          coordenadas2Y$y[[j]] <- ymax2 - coordenadas2Y$y[[j]]
        }
        solucionMirrorY <- solveCoord(coordenadas1, coordenadas2Y, xmax2, ymax2)
        solucionMirrorY[["mirrorx"]] <- 0
        solucionMirrorY[["mirrory"]] <- 10
        solucionMirrorY[["e"]] <- 1
        coordCalc[["solucionMirrorY"]] <- calcNewCoord(coordenadas2, solucionMirrorY, xmax2, ymax2)
    
        # Solve transformations with mirroring on both x and y axes
        coordenadas2XY <- coordenadas2
        for (j in seq_along(coordenadas2XY$x)) {
          coordenadas2XY$x[[j]] <- xmax2 - coordenadas2XY$x[[j]]
        }
        for (j in seq_along(coordenadas2XY$y)) {
          coordenadas2XY$y[[j]] <- ymax2 - coordenadas2XY$y[[j]]
        }
        solucionMirrorXY <- solveCoord(coordenadas1, coordenadas2XY, xmax2, ymax2)
        solucionMirrorXY[["mirrorx"]] <- 10
        solucionMirrorXY[["mirrory"]] <- 1
        solucionMirrorXY[["e"]] <- 1
        coordCalc[["solucionMirrorXY"]] <- calcNewCoord(coordenadas2, solucionMirrorXY, xmax2, ymax2)

      } else if (mode == "procrustes") {
        
        # Solve for the original orientation (without mirroring)
        matProb <- matrix(data = unlist(coordenadas2), ncol = 2)
        proc <- IMIFA::Procrustes(X = matProb, Xstar = matrixCoord1, translate = TRUE, dilate = scale, sumsq = FALSE)
        solucionOrig <- resultProcrustes(proc, 0, 0, scale)
        coordCalc[["solucionOrig"]] <- proc$X.new
    
        # Solve with mirror on x-axis
        coordenadas2X <- coordenadas2
        for (j in seq_along(coordenadas2X$x)) {
          coordenadas2X$x[[j]] <- xmax2 - coordenadas2X$x[[j]]
        }
        matProb <- matrix(data = unlist(coordenadas2X), ncol = 2)
        proc <- IMIFA::Procrustes(X = matProb, Xstar = matrixCoord1, translate = TRUE, dilate = scale, sumsq = FALSE) 
        solucionMirrorX <- resultProcrustes(proc, 10, 0, scale)
        coordCalc[["solucionMirrorX"]] <- proc$X.new
    
        # Solve with mirror on y-axis
        coordenadas2Y <- coordenadas2
        for (j in seq_along(coordenadas2Y$y)) {
          coordenadas2Y$y[[j]] <- ymax2 - coordenadas2Y$y[[j]]
        }
        matProb <- matrix(data = unlist(coordenadas2Y), ncol = 2)
        proc <- IMIFA::Procrustes(X = matProb, Xstar = matrixCoord1, translate = TRUE, dilate = scale, sumsq = FALSE)
        solucionMirrorY <- resultProcrustes(proc, 0, 10, scale)
        coordCalc[["solucionMirrorY"]] <- proc$X.new
        
        # Solve with mirror on both x and y axes
        coordenadas2XY <- coordenadas2
        for (j in seq_along(coordenadas2XY$x)) {
          coordenadas2XY$x[[j]] <- xmax2 - coordenadas2XY$x[[j]]
        }
        for (j in seq_along(coordenadas2XY$y)) {
          coordenadas2XY$y[[j]] <- ymax2 - coordenadas2XY$y[[j]]
        }
        matProb <- matrix(data = unlist(coordenadas2XY), ncol = 2)
        proc <- IMIFA::Procrustes(X = matProb, Xstar = matrixCoord1, translate = TRUE, dilate = scale, sumsq = FALSE)
        solucionMirrorXY <- resultProcrustes(proc, 1, 1, scale)
        coordCalc[["solucionMirrorXY"]] <- proc$X.new
      }
    
      # Store each transformation option in a list for the current image
      option <- list(
        solucionOrig     = solucionOrig,
        solucionMirrorX  = solucionMirrorX,
        solucionMirrorY  = solucionMirrorY,
        solucionMirrorXY = solucionMirrorXY
      )
      listaOpciones[[i]] <- option

      # Select the transformation option with the lowest sum of squares
        
      # Calculate the minimum dimensions between the reference image and the current image
      xmax <- min(nrow(object.semla@tools$Staffli@rasterlists$raw[[1]]),
                  nrow(object.semla@tools$Staffli@rasterlists$raw[[i]]))
      ymax <- min(ncol(object.semla@tools$Staffli@rasterlists$raw[[1]]),
                  ncol(object.semla@tools$Staffli@rasterlists$raw[[i]]))
      
      # Calculate transformation parameters for each option
      todosvalores <- lapply(option, function(solucion) calcParameters(solucion, xmax, ymax, mode))
      # Store calculated transformation values for each option in a separate list
      listaOpcionesCalc[[i]] <- todosvalores
      
      # Calculate the sum of squares for each transformation option to find the optimal alignment
      suma_de_cuadrados <- sapply(todosvalores, function(valores) {
        if (abs(valores[["trx"]]) >= 1 || abs(valores[["try"]]) >= 1) {
          return(Inf)
        } else {return(sum(valores^2))}
      })
      
      indice_fila_minima <- names(which.min(suma_de_cuadrados)) # Name of minimum sum of squares
      print(indice_fila_minima)
      
      # Determine mirroring settings based on the index of the optimal transformation
      if (indice_fila_minima == "solucionOrig")     { mirrorx <- FALSE;  mirrory <- FALSE }
      if (indice_fila_minima == "solucionMirrorX")  { mirrorx <- TRUE ;  mirrory <- FALSE }
      if (indice_fila_minima == "solucionMirrorY")  { mirrorx <- FALSE;  mirrory <- TRUE  }
      if (indice_fila_minima == "solucionMirrorXY") { mirrorx <- TRUE ;  mirrory <- TRUE  }
      
      # Update the solution with the optimal transformation option and mirror settings
      solucion <- listaOpciones[[i]][[indice_fila_minima]]
      solucion[["mirrorx"]] <- mirrorx
      solucion[["mirrory"]] <- mirrory
      listaSoluciones[[i]] <- solucion
      
      # Store the selected transformation parameters for later use
      listaValores[[i]] <- listaOpcionesCalc[[i]][[indice_fila_minima]]
      print(solucion)
      print(listaValores[[i]])
        
      # Update coordinates based on selected transformation
      if (mode == "GTEM") { 
        listaCoordenadasNEW[[i]] <- coordCalc[[indice_fila_minima]]
      } else if (mode == "procrustes") {
          coordNEWmatrix <- coordCalc[[indice_fila_minima]]
          coordNEW <- list()
          coordNEW$x <- round(coordNEWmatrix[,1])
          coordNEW$y <- round(coordNEWmatrix[,2])
          listaCoordenadasNEW[[i]] <- coordNEW 
      }
    }
  }

  # Iterate over each image in the dataset (starting from the second image)
  for (i in 2:length(object.semla@tools$Staffli@rasterlists$raw)) {
    
    if (mode == "GTEM" ) {
      if (scale == TRUE) {
        # Up to this point, TR has been applied correctly, 
        # now it will be calculated if any scaling modifications would be necessary.
        valueScal <- round(calcScal(listaCoordenadas[[1]], listaCoordenadasNEW[[i]]),2) 
        if (valueScal != 1) {
          listaCoordenadasNEW[[i]][['x']] <- listaCoordenadasNEW[[i]][['x']] * valueScal
          listaCoordenadasNEW[[i]][['y']] <- listaCoordenadasNEW[[i]][['y']] * valueScal
        }
        listaValores[[i]][["e"]] <- valueScal
      } else if (scale == FALSE) {listaValores[[i]][["e"]] <- 1}
      # Display the selected transformation parameters
      valores <- listaValores[[i]]

    } else if (mode == "procrustes") {
      # Display the selected transformation parameters
      valores <- listaValores[[i]]

    } else if (mode == "RVSSimageJ") {
      if (scale == FALSE) {
        # Retrieve and convert transformation parameters
        solucion <- listaSoluciones[[i]]
        valores <- solucion
        valores$angulo <- solucion$angle * (180/pi)  # Convert angle to degrees
        
        # Adjust dx and dy based on image center for transformation normalization
        valores$trx <-  (solucion$dx - (xmax/2 - xmax/2 * cos(solucion$angle) + ymax/2 * sin(solucion$angle))) / xmax
        valores$try <-  (solucion$dy - (ymax/2 - ymax/2 * cos(solucion$angle) - xmax/2 * sin(solucion$angle))) / ymax

	      valores$e <- 1        

      } else if (scale == TRUE) {
        # Retrieve and convert transformation parameters
        solucion <- listaSoluciones[[i]]
        valores <- list()
        
        # Calculate the angle based on the cosine and sine values
        angle_rad <- atan2(solucion[["s*sin"]], solucion[["s*cos"]])
        valores$angulo <- angle_rad * (180/pi)  # Convert angle to degrees
        
        # Calculate the scale factor extracting from the estimated values
        valores$e <- sqrt(solucion[["s*sin"]]^2 + solucion[["s*cos"]]^2)
        if (valores$e > 3) valores$e <- 3
        
        # Adjust dx and dy based on image center for transformation normalization
        valores$trx <-  (solucion$dx - (xmax/2 - xmax/2 * cos(angle_rad) + ymax/2 * sin(angle_rad))) / xmax
        valores$try <-  (solucion$dy - (ymax/2 - ymax/2 * cos(angle_rad) - xmax/2 * sin(angle_rad))) / ymax
        
        valores$dx <- ifelse(valores$dx > 1, 1, valores$dx)
        valores$dy <- ifelse(valores$dy > 1, 1, valores$dy)
      }
    }

    print(valores)
  
    # Initialize a list to store all transformations applied to the image
    alltransforms <- list()
    original[[i]] <- object.semla@tools$Staffli@rasterlists$raw[[i]]
  
    # Apply mirror transformation if the mirror values are non-zero
    if (mode != "RVSSimageJ") { if (valores[['mirrorx']] != 0 || valores[['mirrory']] != 0) {
      transforms_mirror <- generate_rigid_transform(sampleID = i, 
                                                    mirror_x = as.logical(valores[['mirrorx']]),
                                                    mirror_y = as.logical(valores[['mirrory']]))
      object.semla <- RigidTransformImages(object.semla, transforms = transforms_mirror)
      alltransforms[[1]] <- transforms_mirror
      object.semla@tools$Staffli@rasterlists$raw[[i]] <- object.semla@tools$Staffli@rasterlists$transformed[[i]]
    
      # Update meta_data with new transformations
      object.semla@tools[["Staffli"]]@meta_data[
        which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- object.semla@tools[["Staffli"]]@meta_data[
          which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      object.semla@tools[["Staffli"]]@meta_data[
        which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- object.semla@tools[["Staffli"]]@meta_data[
          which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }}
  
    # Apply rotation if the angle is non-zero
    angulo <- valores[['angulo']]
    if ( angulo != 0 && angulo != 360 ){
      transforms_angle <- generate_rigid_transform(sampleID = i, 
                                                  angle = angulo)
      object.semla <- RigidTransformImages(object.semla, transforms = transforms_angle)
      alltransforms[[2]] <- transforms_angle
      object.semla@tools$Staffli@rasterlists$raw[[i]] <- object.semla@tools$Staffli@rasterlists$transformed[[i]]
    
      # Update meta_data with new transformations
      object.semla@tools[["Staffli"]]@meta_data[
        which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- object.semla@tools[["Staffli"]]@meta_data[
          which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      object.semla@tools[["Staffli"]]@meta_data[
        which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- object.semla@tools[["Staffli"]]@meta_data[
          which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
  
    # Apply translation if dx or dy is non-zero
    if ( valores[['trx']] != 0 || valores[['try']] != 0 ) {
      transforms_trans <- generate_rigid_transform(sampleID = i, 
                                                  tr_x = valores[['trx']], #round(valores[[2]], digits = 2), 
                                                  tr_y = valores[['try']]) #round(valores[[3]], digits = 2),
      object.semla <- RigidTransformImages(object.semla, transforms = transforms_trans)
      alltransforms[[3]] <- transforms_trans
      object.semla@tools$Staffli@rasterlists$raw[[i]] <- object.semla@tools$Staffli@rasterlists$transformed[[i]]
    
      # Update meta_data with new transformations
      object.semla@tools[["Staffli"]]@meta_data[
        which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- object.semla@tools[["Staffli"]]@meta_data[
          which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      object.semla@tools[["Staffli"]]@meta_data[
        which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- object.semla@tools[["Staffli"]]@meta_data[
          which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
  
    # Apply scaling if the scale factor is non-zero
    if ( valores[['e']] != 1 ) {
      transforms_scale <- generate_rigid_transform(sampleID = i, 
                                                  scalefactor = valores[['e']])
      object.semla <- RigidTransformImages(object.semla, transforms = transforms_scale)
      alltransforms[[4]] <- transforms_scale
      object.semla@tools$Staffli@rasterlists$raw[[i]] <- object.semla@tools$Staffli@rasterlists$transformed[[i]]
    
      # Update meta_data with new transformations
      object.semla@tools[["Staffli"]]@meta_data[
        which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- object.semla@tools[["Staffli"]]@meta_data[
          which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      object.semla@tools[["Staffli"]]@meta_data[
        which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- object.semla@tools[["Staffli"]]@meta_data[
          which(object.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
  
    # Store aligned images and transformation lists
    aligned[[i]] <- object.semla@tools$Staffli@rasterlists$transformed[[i]]
    object.semla@tools$Staffli@rasterlists$raw[[i]] <- original[[i]]
    listaTransforms[[i]] <- alltransforms
  
    # Save transformed coordinates for each sample
    CoordinatesSemlaCol[[i]] <- GetCoordinates(object.semla)[which(GetCoordinates(object.semla)[6] == i), 4]
    CoordinatesSemlaRow[[i]] <- GetCoordinates(object.semla)[which(GetCoordinates(object.semla)[6] == i), 5]
  }

  # Update transformed image list in Seurat object
  object.semla@tools$Staffli@rasterlists[["transformed"]] <- aligned

  "Since I do the transformations sequentially, every time you do one 
  new, only the modification of the last one you make is saved. That's why I'm 
  saving at the end of the transformations of each image the new values of 
  the positions, once done on all the images, to save the 
  definitive versions of each of the images in the Seurat object in Staffli by semla.
  "

  # Aggregate transformed column and row coordinates
  CoordCol <- integer()
  CoordRow <- integer()
  for (i in seq_along(CoordinatesSemlaCol)) {
    CoordCol <- c(CoordCol, CoordinatesSemlaCol[[i]][[1]])
    CoordRow <- c(CoordRow, CoordinatesSemlaRow[[i]][[1]])
  }
  object.semla@tools$Staffli@meta_data$pxl_col_in_fullres_transformed <- CoordCol
  object.semla@tools$Staffli@meta_data$pxl_row_in_fullres_transformed <- CoordRow

  png(paste0(saveDir, "imagesOriginal",mode,".png"))
  ImagePlot(object.semla)
  dev.off()

  png(paste0(saveDir, "imagesTransformed",mode,".png"))
  ImagePlot(object.semla, image_use = "transformed")
  dev.off()

  object.seurat <- UpdateSeuratFromSemla(object.semla, image_use = "transformed")
  saveRDS(object.seurat, paste0(saveDir, "objectAligned_merge_",mode,".rds"))

  if (mode == "RVSSimageJ") {
    saveObj <- list(alignedObj = object.seurat)
  } else if (mode == "GTEM" | mode == "procrustes") {
    saveObj <- list(alignedObj = object.seurat,
                    listCoord = listaCoordenadas,
                    listCoordNew = listaCoordenadasNEW)
  }
  saveRDS(saveObj, paste0(saveDir, "objectAligned_merge_",mode,".rds"))
  return(saveObj)
                
}    


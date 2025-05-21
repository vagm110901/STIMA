#' utils 
#' selectCoordPython(image)
#'
#' Allows the user to select multiple points on the provided image by clicking on it.
#' It uses the image of an object created by python (STalign or PASTE2).
#'
#' @param image Image where landmarks are selected 
#' @return list containing the coordinates of the selected points with two elements: x and y.
#' @import imager
#' @export
selectCoordPython <- function(image) {
  par(mfrow=c(1,1))
  plot(c(0,0), xlim=c(0,ncol(image)), ylim=c(0,nrow(image)),
       xlab=NA, ylab=NA, axes = FALSE, asp=1)
  rasterImage(image, xleft = 0, xright = ncol(image),
              ytop = 0, ybottom = nrow(image), interpolate = FALSE)
  coordinates <- locator(type = "p")
  return(coordinates)
}

#' utils
#' escale(object.seurat, new_coords, patientType = c('unique','multiple'), N) 
#' 
#' This function rescales the coordinates of the images to match the original image dimensions after PASTE2 alignment.
#' 
#' @param object.seurat Seurat object containing the images.
#' @param new_coords Data frame containing the new coordinates to be rescaled.
#' @param patientType Type of patient data. Options are 'unique' or 'multiple'.
#' @param N Slice number to be rescaled.
#' @return Data frame with rescaled coordinates.
escale <- function(object.seurat, new_coords, patientType = c('unique','multiple'), N) {
  if (patient == "multiple") {
    # Select the max values changed (row = col & col = row)
    newCOLmax <- max(new_coords$imagerow)
    newROWmax <- max(new_coords$imagecol)
    
    origCOLmax <- max(objSeurat@images[[paste0("slice1.", N)]]@coordinates$imagecol)
    origROWmax <- max(objSeurat@images[[paste0("slice1.", N)]]@coordinates$imagerow)
    
    new_coords$imagerow <- new_coords$imagerow / newCOLmax * origCOLmax
    new_coords$imagecol <- new_coords$imagecol / newROWmax * origROWmax
  }
  return(new_coords)
}

#' seurat2pythonCommunication
#' saveSeurat_forAnnData(seuratObj, patientType = c('unique','multiple'), samples = c("19","43","47","45"))
#' 
#' This function saves the Seurat object in a format compatible with AnnData.
#' It saves the expression matrix, cell metadata, image coordinates, scale factors,
#' and gene metadata for each slice of the Seurat object.
#' 
#' @param object Seurat object to be saved.
#' @param patientType Type of patient data to be saved. Options are 'unique' or 'multiple'.
#' @param samples Vector of sample identifiers. Default is c("19","43","47","45"). It is needed for multiple patients.
#' @return None
#' @export
saveSeurat_forAnnData <- function(object, patientType = c('unique','multiple'), samples = c("19","43","47","45")) {
  patientType <- match.arg(patientType)
  
  if (!dir.exists("./results/AnnData/")) {
    dir.create("./results/AnnData/", recursive = TRUE)
  }
  saveDir <- "./results/AnnData/"
  
  for (i in seq_along(object@images)) {

    # Expression data = counts
    counts_layer <- object@assays$Spatial@layers[[paste0("counts.", i)]]
    counts_layer_T <- t(counts_layer)
    write.csv(as.matrix(counts_layer_T), paste0(saveDir,"expression_matrix_slice", i, ".csv"))

    # Cell metadata 
    if (patientType == "multiple") {
      sample <- samples[i]
      cell_metadata <- object@meta.data[object@meta.data$name == paste0("slice",sample,"_1"), ]
    } else if (patient == "unique") {
      sample <- samples[1]
      cell_metadata <- object@meta.data[object@meta.data$name == paste0("slice",sample,"_", i), ]
    }
    write.csv(cell_metadata, paste0(saveDir, "cell_metadata_slice", i, ".csv"))
  
    # Spatial coordinates (row, col, imagerow, imagecol)
    image_coords <- object@images[[paste0("slice1.", i)]]@coordinates
    image_coords_df <- as.data.frame(image_coords)  # Convertir en un data.frame
    write.csv(image_coords_df, paste0(saveDir, "image_coordinates_slice", i, ".csv"))
    
    #  scale.factors
    scale_factors <- list()
    scale_factors$tissue_lowres_scalef <- object@images[[paste0("slice1.", i)]]@scale.factors$lowres
    scale_factors$tissue_hires_scalef <- object@images[[paste0("slice1.", i)]]@scale.factors$hires
    scale_factors$spot_diameter_fullres <- object@images[[paste0("slice1.", i)]]@scale.factors$spot
    scale_factors$fiducial_diameter_fullres <- object@images[[paste0("slice1.", i)]]@scale.factors$fiducial
    write.csv(scale_factors, paste0(saveDir, "scale_factors_slice", i, ".csv"))
    
    # Gene metadata (shared for all the slices)
    write.csv(rownames(object@assays$Spatial@features[which(object@assays$Spatial@features[,i] == TRUE),]), 
                paste0(saveDir, "gene_metadata", i, ".csv"))
  }

  cat("Images for each slice have been manually copied from the original directory in the directory: ", saveDir, "\n")
}

#' seurat2pythonCommunication
#' PASTE2toSeurat(object, patientType = c('unique','multiple'))
#' 
#' This function loads the PASTE2 data and creates a Seurat object.
#' 
#' @param object Original object to be modified.
#' @param patientType Type of patient data. Options are 'unique' or 'multiple'.
#' @export 
PASTE2toSeurat <- function(object, patientType = c('unique','multiple')) {
  patientType <- match.arg(patientType)
  
  if (!dir.exists("./results/PASTE2/")) {
    dir.create("./results/PASTE2/", recursive = TRUE)
  }
  saveDir <- "./results/PASTE2/"
  
  Ns <- as.character(seq(2, length(object@images)))

  for (N in Ns) {
    # Load the new coordinates saved and generated as csv files by PASTE2.  
    new_coordsProb <- read.csv(file = paste0(saveDir, "PASTE2_", N, "_align1", N, "_h_coord.csv"))
  
    new_coordsProb <- escale(new_coordsProb, patient, N, object)
  
    object@images[[paste0("slice",N)]] <- object@images[[paste0("slice1.", N)]]
    object@images[[paste0("slice", N)]]@coordinates$imagerow <- unlist(new_coordsProb$imagecol)
    object@images[[paste0("slice", N)]]@coordinates$imagecol <- unlist(new_coordsProb$imagerow)
  }  

  saveRDS(object, file = paste0(saveDir,"objectAligned_merge_PASTE2.rds"))
}

#' seurat2pythonCommunication
#' STaligntoSeurat(object.STalign, object, patientType = c('unique','multiple'))
#' 
#' This function loads the STalign data and creates a Seurat object.
#' 
#' @param object.STalign List object resulting from the STalign pipeline.
#' @param object Original object to be modified.
#' @param patientType Type of patient data. Options are 'unique' or 'multiple'.
#' @export 
STaligntoSeurat <- function(object.STalign, object, patientType = c('unique','multiple')) {
  patientType <- match.arg(patientType)
  
  if (!dir.exists("./results/STalign/")) {
    dir.create("./results/STalign/", recursive = TRUE)
  }
  saveDir <- "./results/STalign/"

  object@images[["slice1"]] <- object@images[["slice1.1"]]

  for (i in c(2:4)) {
    object@images[[paste0("slice",i)]] <- object@images[[paste0("slice1.",i)]]
    object@images[[paste0("slice",i)]]@image <- object.STalign[[i]][["affine+diffeo"]][["image"]]
    object@images[[paste0("slice",i)]]@coordinates[,4:5] <- object.STalign[[i]][["affine+diffeo"]][["coords"]]
  }

  saveRDS(object, file = paste0(saveDir,"objectAligned_merge_STalign.rds"))
}

#' evaluation 
#' calculateEvaluationPython(objeto.seurat, mode = c("STalign", "PASTE2"), 
#'                           listaCoordenadasNEW, listaCoordenadas, 
#'                           patientType = c('unique','multiple'))
#' 
#' Evaluate the alignment of images using various metrics (MSE, SSIM, etc.) for both raw and transformed images.
#' It uses the image of an object created by python (STalign or PASTE2).
#' 
#' @param objeto.seurat Seurat object containing images and coordinates.
#' @param mode Evaluation mode: one of "GTEM", "procrustes", or "RVSSimageJ".
#' @param listaCoordenadasNEW List of new coordinates (optional, used mainly if mode != "RVSSimageJ").
#' @param listaCoordenadas List of original coordinates (optional, used mainly if mode != "RVSSimageJ").
#' @param patientType Patient type, affecting region size ("unique" or "multiple").
#' @return Evaluation
#' @import SpatialPack
#' @import imager
#' @export 
calculateEvaluationPython <- function(objeto.seurat, mode = c("STalign", "PASTE2"), 
                                listaCoordenadasNEW, listaCoordenadas, 
                                patientType = c('unique','multiple')) {

  patientType <- match.arg(patientType)
  mode <- match.arg(mode)

  if (!dir.exists(paste0("./results/",mode,"/"))) {
    dir.create(paste0("./results/",mode,"/"), recursive = TRUE)
  }
  saveDir <- paste0("./results/",mode,"/")
  
  if (modeAbr == "PASTE2") {
    for (i in seq_along(1:4)) {
      coordenadas <- selectCoordPython(objeto.seurat@images[[paste0("slice1.",i)]]@image)
      for (j in seq_along(coordenadas)) {
        coordenadas[[j]] <- round(coordenadas[[j]])
      }
      listaCoordenadas[[i]] <- coordenadas
    }
  listaCoordenadasNEW <- listaCoordenadas

  } else if (modeAbr == "STalign") {
    listaCoordenadas <- list()
    for (i in seq_along(1:4)) {
      coordenadas <- selectCoordPython(objeto.seurat@images[[paste0("slice1.",i)]]@image)
      for (j in seq_along(coordenadas)) {
        coordenadas[[j]] <- round(coordenadas[[j]])
      }
        listaCoordenadas[[i]] <- coordenadas
    }

    listaCoordenadasNEW <- list()
    for (i in seq_along(1:4)) {
      coordenadas <- selectCoordPython(objeto.seurat@images[[paste0("slice",i)]]@image)
      for (j in seq_along(coordenadas)) {
        coordenadas[[j]] <- round(coordenadas[[j]])
      }
      listaCoordenadasNEW[[i]] <- coordenadas
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

#' deconvolution
#' createDeconvolutionListsPython(object.seurat, split_by = "name", mode = c("STalign", "PASTE2"))
#'
#' Creates a list of Seurat objects for deconvolution analysis.
#' It uses an objetc created by python (STalign or PASTE2).
#'
#' @param object.seurat A Seurat object after the alignment. 
#'                      If mode == STalign, it should be the Seurat object after the STalign alignment. {STaligntoSeurat()}
#'                      If mode == PASTE2, it should be the Seurat object after the PASTE2 alignment. {PASTE2toSeurat()}
#' @param object A Seurat object with the original data.
#' @param split_by A character string indicating the metadata field to split the Seurat object by. Default is "name".
#' @param mode A character string indicating the deconvolution method to be used. Options are "GTEM", "procrustes", or "RVSSimageJ". Default is "GTEM". 
#' @return A list of Seurat objects, each containing a reference and two problematic objects for deconvolution analysis.
#' @import Seurat
#' @export 
createDeconvolutionListsPython <- function(object.seurat, object,
                                           split_by = "name",  
                                           mode = c("STalign", "PASTE2")) {

  mode <- match.arg(mode)

  if (!dir.exists(paste0("./results/",mode,"/"))) {
    dir.create(paste0("./results/",mode,"/"), recursive = TRUE)
  }
  saveDir <- paste0("./results/",mode,"/")
  
  # Obtain the individual object of the slide
  # Splits the Seurat object into a list of Seurat objects based on the 'name' metadata
  object.seurat.split.orig <- SplitObject(object, split.by = split_by)
  
  # Create a copy of the split objects for transformed images
  object.seurat.split.trans <- object.seurat.split.orig

  if (mode == "STalign") {
    # Assign transformed images and coordinates to the corresponding objects in the new split list
    for (i in 2:length(object.seurat.split.orig)) {
      object.seurat.split.trans[[i]]@images[[paste0("slice1.",i)]]@image <- 
        object.seurat@images[[paste0("slice",i)]]@image

      object.seurat.split.trans[[i]]@images[[paste0("slice1.",i)]]@coordinates[,4:5] <- 
        object.seurat@images[[paste0("slice",i)]]@coordinates[,4:5]
    }
  } else if (mode == "PASTE2") {
    # Assign transformed images to the corresponding objects in the new split list
    for (i in 2:length(object.seurat.split.orig)) {
      object.seurat.split.trans[[i]]@images[[paste0("slice1.",i)]] <- 
        object.seurat@images[[paste0("slice",i)]]
    }
  }

  # Loop through each of the split original objects starting from the second one
  for (i in 2:length(object.seurat.split.orig)) {
    listaObjDeconv <- list()
    
    # Create a list to hold reference and problematic objects for deconvolution
    listaObjDeconv[["reference"]] <- object.seurat.split.orig[[1]]
    listaObjDeconv[["problemNOTalign"]] <- object.seurat.split.orig[[i]]
    listaObjDeconv[["problemYESalign"]] <- object.seurat.split.trans[[i]]
    
    # Update the name in the metadata for each object to reflect its status
    listaObjDeconv[["reference"]]@meta.data$name <- gsub(
      listaObjDeconv[["reference"]]@meta.data$name[[1]], 
      "reference", 
      listaObjDeconv[["reference"]]@meta.data$name
    )
    listaObjDeconv[["problemNOTalign"]]@meta.data$name <- gsub(
      listaObjDeconv[["problemNOTalign"]]@meta.data$name[[1]], 
      "problemNOTalign", 
      listaObjDeconv[["problemNOTalign"]]@meta.data$name
    )
    listaObjDeconv[["problemYESalign"]]@meta.data$name <- gsub(
      listaObjDeconv[["problemYESalign"]]@meta.data$name[[1]], 
      "problemYESalign", 
      listaObjDeconv[["problemYESalign"]]@meta.data$name
    )
    
    # Guardar la list
    saveRDS(listaObjDeconv, paste0(saveDir, "objectAligned_merge_",mode,"_list_im",i,".rds"))
  }
}


#´ alignment
#' alignmentSTalign(object)
#' 
#' @param object first object with the slices/samples merged together
#' patientType Type of patient data. Options are 'unique' or 'multiple'.
#' @import reticulate
#' @export 
alignmentSTalign <- function(object, patientType = c('unique','multiple')) {
  patientType <- match.arg(patientType)
  pkg_path <- system.file("python", package = "STIMA")
  reticulate::py_run_file(file.path(pkg_path, "publishedMethodsComparison.py"))
  funct_py <- reticulate::import_from_path("publishedMethodsComparison", path = pkg_path)

  funct_py$pythonEnviroment("STalign")

  if (!dir.exists("./results/STalign/")) {
    dir.create("./results/STalign/", recursive = TRUE)
  }
  saveDir <- "./results/STalign/"

  listSTalignResults <- list()
  listSTalignResults[[1]] <- list()

  ref_image_coords <- object@images[["slice1.1"]]@coordinates
  ref_pos <- as.data.frame(ref_image_coords)  # Convert to a data frame
  ref_barcodes <- rownames(ref_pos)
  ref_scalefactor <- object@images[["slice1.1"]]@scale.factors$lowres
  ref_image <- object@images[["slice1.1"]]@image

  posCalc_ref <- ref_pos[ref_barcodes,5:4] * ref_scalefactor
  colnames(posCalc_ref) <- c("x", "y")
  max_yref <- max(posCalc_ref$y, na.rm = TRUE)
  min_yref <- min(posCalc_ref$y, na.rm = TRUE)

  listSTalignResults[[1]][["image"]] <- ref_image
  listSTalignResults[[1]][["coords"]] <- ref_pos

  # source (reference)
  im_ref_norm <- funct_py$normalize_images(ref_image)

  # Reference landmarks selection
  #
  # To assist with the spatial alignment, I will place a few landmark points.
  # This will help mitigate the likelihood of us falling into a local minimum 
  # in the gradient descent and arrive at a suboptimal solution. I will just
  # manually create them by eyeballing the image. They can be very approximate
  # as STalign will integrate these landmarks with other imaging features in 
  # its optimization.
  
  coordenadas_im_ref <- selectCoordPython(py$im_ref_norm)
  for (j in seq_along(coordenadas_im_ref)) {
    coordenadas_im_ref[[j]] <- round(coordenadas_im_ref[[j]])
  }

  points_im_ref <- t(data.frame(a = c(coordenadas_im_ref$x[[1]], coordenadas_im_ref$y[[1]]), 
                                b = c(coordenadas_im_ref$x[[2]], coordenadas_im_ref$y[[2]]), 
                                c = c(coordenadas_im_ref$x[[3]], coordenadas_im_ref$y[[3]]),
                                d = c(coordenadas_im_ref$x[[4]], coordenadas_im_ref$y[[4]]),
                                e = c(coordenadas_im_ref$x[[5]], coordenadas_im_ref$y[[5]])))
  colnames(points_im_ref) <- c('x', 'y')

  points_im_ref_rc <- points_im_ref[,2:1]
  x_im_ref <- posCalc_ref[,1]
  y_im_ref <- posCalc_ref[,2]


  for (i in seq_along(2:length(object@images))) {
    listSTalignResults[[im]] <- list()
  
    prob_image_coords <- object@images[[paste0("slice1.",im)]]@coordinates
    prob_pos <- as.data.frame(prob_image_coords)  # Convertir en un data.frame
    prob_barcodes <- rownames(prob_pos)
    prob_scalefactor <- object@images[[paste0("slice1.",im)]]@scale.factors$lowres
    prob_image <- object@images[[paste0("slice1.",im)]]@image
      
    posCalc_prob <- prob_pos[prob_barcodes,5:4] * prob_scalefactor
    colnames(posCalc_prob) <- c("x", "y")
    max_yprob <- max(posCalc_prob$y, na.rm = TRUE)
    min_yprob <- min(posCalc_prob$y, na.rm = TRUE)

    # target (problem)
    im_prob_norm <- funct_py$normalize_images(prob_image) 

    coordenadas_im_prob <- selectCoordPython(py$im_prob_norm)
    for (j in seq_along(coordenadas_im_prob)) {
      coordenadas_im_prob[[j]] <- round(coordenadas_im_prob[[j]])
    }

    ## target (problem)
    points_im_prob <- t(data.frame(a = c(coordenadas_im_prob$x[[1]], coordenadas_im_prob$y[[1]]), 
                                  b = c(coordenadas_im_prob$x[[2]], coordenadas_im_prob$y[[2]]), 
                                  c = c(coordenadas_im_prob$x[[3]], coordenadas_im_prob$y[[3]]),
                                  d = c(coordenadas_im_prob$x[[4]], coordenadas_im_prob$y[[4]]),
                                  e = c(coordenadas_im_prob$x[[5]], coordenadas_im_prob$y[[5]])))
    colnames(points_im_prob) <- c('x', 'y')

    points_im_prob_rc <- points_im_prob[,2:1]
    x_im_prob <- posCalc_prob[,1]
    y_im_prob <- posCalc_prob[,2]

    result <- funct_py$STalign_transformation(points_im_prob_rc, points_im_ref_rc, 
                                              im_prob_norm, im_ref_norm,
                                              y_im_prob, x_im_prob)

    align_image_im_prob <- result$align_image_im_prob
    align_points_im_prob <- result$align_points_im_prob

    transf <- "affine+diffeo"

    listSTalignResults[[im]][[transf]][["image"]] <- align_image_im_prob
  
    posAligned <- align_points_im_prob ## currently in row-col order
    posAligned <- data.frame(posAligned[,2:1]) ## put into x-y order
    rownames(posAligned) <- prob_barcodes
    colnames(posAligned) <- c('x', 'y')
    
    posAligned_adapt <- posAligned[prob_barcodes,2:1] / prob_scalefactor
    colnames(posAligned_adapt) <- c("imagerow", "imagecol")
    
    listSTalignResults[[im]][[transf]][["coords"]] <- posAligned_adapt
  }

  # Save the list of aligned images and coordinates
  saveRDS(listSTalignResults, paste0(saveDir,"objectAligned_list_STalign.rds"))

  STaligntoSeurat(listSTalignResults, object, patientType)
}

#´ alignment
#' alignmentPASTE2(object)
#' 
#' @param object first object with the slices/samples merged together
#' @param patientType Type of patient data. Options are 'unique' or 'multiple'.
#' @param samples Vector of sample identifiers. Default is c("19","43","47","45"). It is needed for multiple patients. 
#' @import reticulate
#' @export 
alignmentPASTE2 <- function(object, patientType = c('unique','multiple'), samples = c("19","43","47","45")) {
  patientType <- match.arg(patientType)
  pkg_path <- system.file("python", package = "STIMA")
  reticulate::py_run_file(file.path(pkg_path, "publishedMethodsComparison.py"))
  funct_py <- reticulate::import_from_path("publishedMethodsComparison", path = pkg_path)

  funct_py$pythonEnviroment("PASTE2")

  if (!dir.exists("./results/AnnData/")) {
    dir.create("./results/AnnData/", recursive = TRUE)
  }
  if (!dir.exists("./results/PASTE2/")) {
    dir.create("./results/PASTE2/", recursive = TRUE)
  }
  carpetaData <- "./results/AnnData/" 
  saveDir <- "./results/PASTE2/"

  saveSeurat_forAnnData(object, patientType = patientType, samples = samples)
  funct_py$h5ad_create()
  funct_py$RGBvalues()
  funct_py$PASTE2_align()
  funct_py$saveAnnData_forSeurat()
  PASTE2toSeurat(object, patientType)
}

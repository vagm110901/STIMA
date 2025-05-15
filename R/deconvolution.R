#' createDeconvolutionLists(object.seurat, split_by = "name", mode = c("GTEM", "procrustes", "RVSSimageJ"))
#'
#' Creates a list of Seurat objects for deconvolution analysis.
#'
#' @param object.seurat A Seurat object after the alignment containing the data to be split. 
#' @param split_by A character string indicating the metadata field to split the Seurat object by. Default is "name".
#' @param mode A character string indicating the deconvolution method to be used. Options are "GTEM", "procrustes", or "RVSSimageJ". Default is "GTEM". 
#' @return A list of Seurat objects, each containing a reference and two problematic objects for deconvolution analysis.
#' @import Seurat
#' @export 
createDeconvolutionLists <- function(object.seurat, 
                                     split_by = "name",  
                                     mode = c("GTEM", "procrustes", "RVSSimageJ")) {

  mode <- match.arg(mode)

  if (!dir.exists("./results/")) {
    dir.create("./results/", recursive = TRUE)
  }
  saveDir <- "./results/"
  
  # Obtain the individual object of the slide
  # Splits the Seurat object into a list of Seurat objects based on the 'name' metadata
  object.seurat.split.orig <- SplitObject(object.seurat, split.by = split_by)
  
  # Create a copy of the split objects for transformed images
  object.seurat.split.trans <- object.seurat.split.orig
  
  # Assign transformed images to the corresponding objects in the new split list
  for (i in 2:length(object.seurat.split.orig)) {
    object.seurat.split.trans[[i]]@images[[paste0("slice1.",i)]] <- object.seurat@images[[paste0("slice",i)]]
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


#' as_AssayObject(object)
#' 
#' Converts an RCTD object to a Seurat Assay object.
#' 
#' @param object An RCTD object.
#' @return A Seurat Assay object containing the deconvolution results.
#' @import Seurat
#' @import Matrix
#' @import data.table
#' @import spacexr
as_AssayObject <- function(object) {
  if (is(object, "RCTD")) {
    if (!requireNamespace("spacexr", quietly = TRUE)) {
      stop("Install spacexr.")
    }
  } else {
    stop("Only RCTD objects supported")
  }
  r <- object@results
  if (length(r) > 1) {
    if (!is.null(r[[1]]$sub_weights)) {
      sw <- data.table::rbindlist(lapply(seq_along(r), function(i)
        data.table::data.table(
          barcode = colnames(object@spatialRNA@counts)[i],
          cell_type = r[[i]]$cell_type_list,
          weight = r[[i]]$sub_weights
        )), fill = TRUE)
      sw$cell_type[is.na(sw$cell_type)] <- "unassigned"
      swd <- data.table::dcast(sw, barcode ~ cell_type, value.var = "weight", fill = 0)
      swm <- as.matrix(swd[, -1])
      rownames(swm) <- swd$barcode
      swm <- t(swm)
      #swm <- spacexr::normalize_weights(swm)
      #swm <- rbind(swm, max = apply(swm[!rownames(swm) %in% "unassigned", ], 2, max))
      swm <- as(swm, "sparseMatrix")
      return(CreateAssayObject(data = swm))
    }
  } else if (length(r) == 1) {
    m <- t(spacexr::normalize_weights(as.matrix(r$weights)))
    m <- rbind(m, max = apply(m, 2, max))
    return(CreateAssayObject(data = m))
  }
}

#' as_AssayObject_complete(object)
#' 
#' Converts an RCTD object to a Seurat Assay object with complete deconvolution results.
#' 
#' @param object An RCTD object.
#' @return A Seurat Assay object containing the deconvolution results.
#' @import Seurat
#' @import Matrix
#' @import data.table
#' @import spacexr
as_AssayObject_complete <- function(object) {
  if (is(object, "RCTD")) {
    if (!requireNamespace("spacexr", quietly = TRUE)) {
      stop("Install spacexr.")
    }
  } else {
    stop("Only RCTD objects supported")
  }
  r <- object@results
  if (length(r) > 1) {
    if (!is.null(r[[1]]$all_weights)) {
      sw <- data.table::rbindlist(lapply(seq_along(r), function(i)
        data.table::data.table(
          barcode = colnames(object@spatialRNA@counts)[i],
          cell_type = names(r[[i]]$all_weights),
          weight = r[[i]]$all_weights
        )), fill = TRUE)
      sw$cell_type[is.na(sw$cell_type)] <- "unassigned"
      swd <- data.table::dcast(sw, barcode ~ cell_type, value.var = "weight", fill = 0)
      swm <- as.matrix(swd[, -1])
      rownames(swm) <- swd$barcode
      swm <- t(swm)
      #swm <- spacexr::normalize_weights(swm)
      #swm <- rbind(swm, max = apply(swm[!rownames(swm) %in% "unassigned", ], 2, max))
      swm <- as(swm, "sparseMatrix")
      return(CreateAssayObject(data = swm))
    }
  } else if (length(r) == 1) {
    m <- t(spacexr::normalize_weights(as.matrix(r$weights)))
    m <- rbind(m, max = apply(m, 2, max))
    return(CreateAssayObject(data = m))
  }
}

#' deconvolutionRCTD()
#' 
#' Performs deconvolution using the RCTD package for a list of reference-NOTalign-YESalign.
#' 
#' @param object.list A list of Seurat objects containing the data to be deconvoluted.
#' @param im A numeric value indicating the image number.
#' @param mode A character string indicating the deconvolution method to be used. Options are "GTEM", "procrustes", or "RVSSimageJ". Default is "GTEM".
#' @return A list of Seurat objects with deconvolution results added as assays.
#' @import Seurat
#' @import spacexr
#' @import Matrix
#' @import data.table
#' @import spacexr
#' @export
deconvolutionRCTD <- function(object.list, im,
                              mode = c("GTEM", "procrustes", "RVSSimageJ")) {
  mode <- match.arg(mode)

  if (!dir.exists("./results/")) {
    dir.create("./results/", recursive = TRUE)
  }
  saveDir <- "./results/"

  for (i in 1:length(object.list)) {
    objectST <- object.list[[i]]
    # we normalize the objects
    objectST <- SCTransform(objectST, assay = "Spatial", verbose = FALSE) %>%
                RunPCA(verbose = FALSE)
        
    objname <- names(object.list[i])
    object.list[[objname]] <- objectST
  }

  listaObjST_RCTD <- list()
      
  for (i in 1:length(object.list)) {
    # 1. coords: A numeric data.frame (or matrix) representing the spatial pixel 
    # locations. rownames are barcodes/pixel names, and there should be two columns 
    # for ‘x’ and for ‘y’.
    coords <- object.list[[i]]@images[[1]]@coordinates[,c("imagecol", "imagerow")]
    colnames(coords) <- c("x", "y") 
    
    # 2. counts: A matrix (or dgCmatrix) representing Digital Gene Expression (DGE). 
    # Rownames should be genes and colnames represent barcodes/pixel names. 
    # Counts should be untransformed count-level data.
    counts <- object.list[[i]]@assays[["SCT"]]@counts
    
    # 3. nUMI: Optional, a named (by pixel barcode) list of total counts or UMI’s 
    # appearing at each pixel. If not provided, nUMI will be assumed to be the total
    # counts appearing on each pixel.
    nUMI <- colSums(counts)
    
    STsample <- SpatialRNA(coords, counts, nUMI)
    
    objname <- names(object.list[i])
    listaObjST_RCTD[[objname]] <- STsample
    
    # RCTD
    # Create RCTD object
    myRCTD <- create.RCTD(listaObjST_RCTD[[i]], reference, max_cores = 8, UMI_min = 10)
    # Run RCTD
    myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')
    
    # results
    assay_myRCTD <- as_AssayObject(myRCTD)
    
    object.list[[i]]@assays[["deconvolution.RCTD"]] <- assay_myRCTD

    assay_myRCTD_complete <- as_AssayObject_complete(myRCTD)

    object.list[[i]]@assays[["deconvolution.RCTD.complete"]] <- assay_myRCTD_complete
  }
  
  # Save the object
  saveRDS(object.list, file.path(saveDir,paste0("objectAligned_RCTD_SCAnnot_merge_", 
                                                mode,"_list_im", im, ".rds")))
                                              
}

#' deconvolutionRCTD_mergeFiles()
#' 
#' Merges the deconvolution results from multiple files into a single list of Seurat objects.
#' 
#' @param modes A character vector indicating the deconvolution methods to be used. Default is c("GTEM", "procrustes", "RVSSimageJ").
#' @param ims A character vector indicating the image numbers to be processed. Default is c("2","3","4").
#' @return A list of Seurat objects with deconvolution results added as assays.
#' @import Seurat
#' @import spacexr
#' @import Matrix
#' @import data.table
#' @import spacexr
#' @export
deconvolutionRCTD_mergeFiles <- function(modes = c("GTEM", "procrustes", "RVSSimageJ"),
                                         ims = c("2","3","4")) {

  if (!dir.exists("./results/")) {
    dir.create("./results/", recursive = TRUE)
  }
  saveDir <- "./results/"

  listaObjDeconv <- list()
  patient <- "1"
  listaObjDeconv[[patient]] <- list()
  for (mode in modes) {
    print(mode)
    listaObjDeconv[[patient]][[mode]] <- list()
    for  (i in ims)  {
      print(i)
      listaObjAnnot <- readRDS(paste0(saveDir, 'objectAligned_RCTD_SCAnnot_merge_', 
                                      mode, '_list_im', i, '.rds'))
      for (j in 1:length(listaObjAnnot)) {
        DefaultAssay(listaObjAnnot[[j]]) <- "deconvolution.RCTD.complete"
        RCTDmode <- "deconvolution.RCTD.complete"
      }
      
      if (i == ims[[1]]) { 
        listaObjDeconv[[patient]][[mode]][["1"]] <- listaObjAnnot$reference 
      }
      listaObjDeconv[[patient]][[mode]][[i]] <- listaObjAnnot$problemYESalign
      
    } # end_for_Nimage 
  } # end_for_mode

  return(listaObjDeconv)
}

#' matrixComparison(listaObjAnnot)
#' 
#' Generates a matrix comparison of cell types across different images.
#' 
#' @param listaObjAnnot list of objects after deconvolution
#' @return rv_long A data frame containing the RV coefficients and p-values for the comparison of cell types across different images.
#' @import grDevices
#' @import data.table
#' @import stats
#' @import FactoMineR
#' @import reshape2
#' @import ggplot2
#' @import ggpubr
#' @import dplyr
#' @import tidyr
#' @export
matrixComparison <- function(listaObjAnnot) {
  library(data.table)
  # Build empty matrixes
  imageMatrixRef <- data.table::data.table(matrix(numeric(0), ncol = 11))
  imageMatrixNOT <- data.table::data.table(matrix(numeric(0), ncol = 11))
  imageMatrixYES <- data.table::data.table(matrix(numeric(0), ncol = 11))
  setnames(imageMatrixRef, cell.types)  
  setnames(imageMatrixNOT, cell.types)  
  setnames(imageMatrixYES, cell.types)  
  
  # Coordinates from the reference image
  objectST1 <- listaObjAnnot[[1]]
  imagerow1 <- objectST1@images[[1]]@coordinates$imagerow
  imagecol1 <- objectST1@images[[1]]@coordinates$imagecol
  
  max_imagerow1 <- max(imagerow1)
  min_imagerow1 <- min(imagerow1)
  max_imagecol1 <- max(imagecol1)
  min_imagecol1 <- min(imagecol1)
  
  # Useful variables
  group <- 0
  listRef <- list()
  listNOT <- list()
  listYES <- list()
  
  colsRef <- listaObjAnnot[[1]]@images[[1]]@coordinates$col
  rowsRef <- listaObjAnnot[[1]]@images[[1]]@coordinates$row
  min_colsRef <- min(colsRef)
  max_colsRef <- max(colsRef)
  min_rowsRef <- min(rowsRef)
  max_rowsRef <- max(rowsRef)
  
  colmin <- min_colsRef
  colmax <- colmin+5
  
  while (colmin <= max_colsRef) {
    #print(paste0("columna ", colmin))
    rowmin <- min_rowsRef
    rowmax <- rowmin+3
    
    while (rowmin <= max_rowsRef) {
      group <- group + 1
      #print(paste0("grupo ", group, " y fila ", rowmin))
      
      ################################ Selection of row and col values ##########
      imagerowmin <- min(imagerow1[objectST1@images[[1]]@coordinates$row >= rowmin & 
                                     objectST1@images[[1]]@coordinates$row < rowmax])
      imagerowmax <- max(imagerow1[objectST1@images[[1]]@coordinates$row >= rowmin & 
                                     objectST1@images[[1]]@coordinates$row < rowmax])
      imagecolmin <- min(imagecol1[objectST1@images[[1]]@coordinates$col >= colmin & 
                                     objectST1@images[[1]]@coordinates$col < colmax])
      imagecolmax <- max(imagecol1[objectST1@images[[1]]@coordinates$col >= colmin & 
                                     objectST1@images[[1]]@coordinates$col < colmax])
      ########################################################################### 
      
      rel_imagerowmin <- (imagerowmin - min_imagerow1) / (max_imagerow1 - min_imagerow1)
      rel_imagerowmax <- (imagerowmax - min_imagerow1) / (max_imagerow1 - min_imagerow1) 
      rel_imagecolmin <- (imagecolmin - min_imagecol1) / (max_imagecol1 - min_imagecol1) 
      rel_imagecolmax <- (imagecolmax - min_imagecol1) / (max_imagecol1 - min_imagecol1) 
      
      for (i in seq_along(listaObjAnnot)) {
        objectST <- listaObjAnnot[[i]]
        imagerow <- objectST@images[[1]]@coordinates$imagerow
        imagecol <- objectST@images[[1]]@coordinates$imagecol
        
        # If the image has been aligned, the coordinates to make a subset must be 
        # selected from the image before the alignment
        if (i == 3) {
          objectST_orig <- listaObjAnnot[[2]]
          imagerow_orig <- objectST_orig@images[[1]]@coordinates$imagerow
          imagecol_orig <- objectST_orig@images[[1]]@coordinates$imagecol
          
          min_imagerow_orig <- min(imagerow_orig)
          max_imagerow_orig <- max(imagerow_orig)
          min_imagecol_orig <- min(imagecol_orig)
          max_imagecol_orig <- max(imagecol_orig)
          
          NEWimagerowmin <- rel_imagerowmin * (max_imagerow_orig - min_imagerow_orig) + min_imagerow_orig
          NEWimagerowmax <- rel_imagerowmax * (max_imagerow_orig - min_imagerow_orig) + min_imagerow_orig
          NEWimagecolmin <- rel_imagecolmin * (max_imagecol_orig - min_imagecol_orig) + min_imagecol_orig
          NEWimagecolmax <- rel_imagecolmax * (max_imagecol_orig - min_imagecol_orig) + min_imagecol_orig
          
        } else {
          min_imagerow <- min(imagerow)
          max_imagerow <- max(imagerow)
          min_imagecol <- min(imagecol)
          max_imagecol <- max(imagecol)
          
          NEWimagerowmin <- rel_imagerowmin * (max_imagerow - min_imagerow) + min_imagerow
          NEWimagerowmax <- rel_imagerowmax * (max_imagerow - min_imagerow) + min_imagerow
          NEWimagecolmin <- rel_imagecolmin * (max_imagecol - min_imagecol) + min_imagecol
          NEWimagecolmax <- rel_imagecolmax * (max_imagecol - min_imagecol) + min_imagecol
        } # end if-else
        
        coordValues <- which(imagerow >= NEWimagerowmin & imagerow <= NEWimagerowmax &
                               imagecol >= NEWimagecolmin & imagecol <= NEWimagecolmax)
        
        if (length(coordValues) != 0) {
          grupo <- rownames(objectST@images[[1]]@coordinates[coordValues,])
          objectST.group <- subset(objectST, cells = grupo, invert = FALSE) 
          
          # We take the cell types weigths from the complete deconvolution data (full mode)
          list_rows <- as.list(sapply(cell.types, function(cell_type) {
            mean(objectST.group@assays[["deconvolution.RCTD.complete"]]@data[cell_type, ])
          }))
        } else {
          grupo <- NA
          list_rows <- as.list(sapply(cell.types, function(cell_type) { NA }))
        } # end if-else
        
        if (i == 1) {
          listRef[[group]] <- list_rows
        } else if (i == 2) {
          listNOT[[group]] <- list_rows
        } else if (i == 3) {
          listYES[[group]] <- list_rows
        }
        
      } # end for listaObjAnnot
      rowmin <- rowmax
      rowmax <- rowmin+3
    } # end for while row
    colmin <- colmax
    colmax <- colmin+5
  } # end for while col
  
  imageMatrixRef <- rbindlist(lapply(listRef, as.list), fill = TRUE)
  imageMatrixNOT <- rbindlist(lapply(listNOT, as.list), fill = TRUE)
  imageMatrixYES <- rbindlist(lapply(listYES, as.list), fill = TRUE)
  
  #table(complete.cases(imageMatrixRef))
  #table(complete.cases(imageMatrixNOT))
  #table(complete.cases(imageMatrixYES))
  
  #sum(imageMatrixRef[complete.cases(imageMatrixRef) & complete.cases(imageMatrixYES)])
  #sum(imageMatrixNOT[complete.cases(imageMatrixNOT) & complete.cases(imageMatrixYES)])
  shared_groups <- complete.cases(imageMatrixRef) & complete.cases(imageMatrixNOT) & complete.cases(imageMatrixYES)
  #table(shared_groups)
  
  coeffEV_Ref_NOT <- FactoMineR::coeffRV(imageMatrixRef[shared_groups],
                                         imageMatrixNOT[shared_groups])
  
  coeffEV_Ref_YES <- FactoMineR::coeffRV(imageMatrixRef[shared_groups],
                                         imageMatrixYES[shared_groups])
  
  coeffEV_NOT_YES <- FactoMineR::coeffRV(imageMatrixNOT[shared_groups],
                                         imageMatrixYES[shared_groups])
  
  
  coeffEV_Ref <- FactoMineR::coeffRV(imageMatrixRef[shared_groups],
                                     imageMatrixRef[shared_groups])
  
  coeffEV_NOT <- FactoMineR::coeffRV(imageMatrixNOT[shared_groups],
                                     imageMatrixNOT[shared_groups])
  
  coeffEV_YES <- FactoMineR::coeffRV(imageMatrixNOT[shared_groups],
                                     imageMatrixNOT[shared_groups])
  
  rv_matrix <- matrix(c(coeffEV_Ref$rv,     coeffEV_Ref_NOT$rv, coeffEV_Ref_YES$rv,
                        coeffEV_Ref_NOT$rv, coeffEV_NOT$rv,     coeffEV_NOT_YES$rv,
                        coeffEV_Ref_YES$rv, coeffEV_NOT_YES$rv, coeffEV_YES$rv),    
                      nrow = 3, byrow = TRUE)
  rownames(rv_matrix) <- c("Reference", "Original", "Transformed")
  colnames(rv_matrix) <- c("Reference", "Original", "Transformed")
  
  pval_matrix <- matrix(c(coeffEV_Ref$p.value,   coeffEV_Ref_NOT$p.value, coeffEV_Ref_YES$p.value,
                          coeffEV_Ref_NOT$p.value, coeffEV_NOT$p.value,     coeffEV_NOT_YES$p.value,
                          coeffEV_Ref_YES$p.value, coeffEV_NOT_YES$p.value, coeffEV_YES$p.value),    
                        nrow = 3, byrow = TRUE)
  rownames(pval_matrix) <- c("Reference", "Original", "Transformed")
  colnames(pval_matrix) <- c("Reference", "Original", "Transformed")
  
  #rvstd_matrix <- matrix(c(coeffEV_Ref$rvstd,     coeffEV_Ref_NOT$rvstd, coeffEV_Ref_YES$rvstd,
  #                         coeffEV_Ref_NOT$rvstd, coeffEV_NOT$rvstd,     coeffEV_NOT_YES$rvstd,
  #                         coeffEV_Ref_YES$rvstd, coeffEV_NOT_YES$rvstd, coeffEV_YES$rvstd),    
  #                       nrow = 3, byrow = TRUE)
  #rownames(rvstd_matrix) <- c("Reference", "Original", "Transformed")
  #colnames(rvstd_matrix) <- c("Reference", "Original", "Transformed")
  
  rv_long <- melt(rv_matrix, varnames = c("Matrix_A", "Matrix_B"), value.name = "RV")
  pval_long <- melt(pval_matrix, varnames = c("Matrix_A", "Matrix_B"), value.name = "p.value")
  #rvstd_long <- melt(rvstd_matrix, varnames = c("Matrix_A", "Matrix_B"), value.name = "RVstd")
  
  rv_long <- merge(rv_long, pval_long, by = c("Matrix_A", "Matrix_B"))
  
  return(rv_long)
}

#' calculateRVcoeff()
#' 
#' Calculates the RV coefficients for different deconvolution methods and images.
#' 
#' @param modes A character vector indicating the deconvolution methods to be used. Default is c("GTEM", "procrustes", "RVSSimageJ").
#' @param ims A character vector indicating the image numbers to be processed. Default is c("2","3","4").
#' @return A list of RV coefficients for each deconvolution method and image.
#' @import Seurat
#' @import data.table
#' @import reshape2
#' @import foreach
#' @import doParallel
#' @import stats
#' @export
calculateRVcoeff <- function(modes = c("GTEM", "procrustes", "RVSSimageJ"),
                             ims = c("2","3","4")) {

  if (!dir.exists("./results/")) {
    dir.create("./results/", recursive = TRUE)
  }
  saveDir <- "./results/"

  cl <- makeCluster(detectCores() - 4)
  registerDoParallel(cl)
  
  patient <- "1"
  print(patient)
  results_list <- foreach(mode = modes, 
                          .packages = c("Seurat", "data.table","reshape2")) %:%
                  foreach(i = ims, 
                          .packages = c("Seurat", "data.table","reshape2")) %dopar% {
      cat(mode)
      cat(i)
      listaObjAnnot <- readRDS(paste0(saveDir,"objectAligned_merge_",mode,"_list_im",i,".rds"))
      
      for (j in 1:length(listaObjAnnot)) {DefaultAssay(listaObjAnnot[[j]]) <- "deconvolution.RCTD.complete"; RCTDmode <- "deconvolution.RCTD.complete"}
      
      # calcular
      rvMatrix <- matrixComparison(listaObjAnnot = listaObjAnnot)
      list(mode = mode, image = i, rv = rvMatrix)
    }
  
  RV_table <- list()
  for (result in unlist(results_list, recursive = FALSE)) {
    if (!is.null(result)) {
      RV_table[[result$mode]][[result$image]] <- result$rv
    }
  }
  saveRDS(RV_table, paste0(saveDir, "objectAligned_RV_merge.rds"))

  stopCluster(cl)

  return(RV_table)
}

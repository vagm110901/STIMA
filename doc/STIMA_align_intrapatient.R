## -----------------------------------------------------------------------------
#remove.packages("STIMA")
#remotes::install_github("vagm110901/STIMA")
library(STIMA)

## ----eval=FALSE---------------------------------------------------------------
# # Define the folder where raw data is located
# inputDir <- system.file("extdata", "exampleData", package = "STIMA")
# 
# paciente_merge <- readRDS(paste0(inputDir,"/Paciente19_merge.rds"))

## ----eval=FALSE---------------------------------------------------------------
# listAlignedGTEM <- STIMA::STIMA(paciente_merge, mode = "GTEM", scale = FALSE)
# listAlignedprocrustes <- STIMA::STIMA(paciente_merge, mode = "procrustes", scale = FALSE)
# listAlignedRVSSimageJ <- STIMA::STIMA(paciente_merge, mode = "RVSSimageJ", scale = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# listAlignedGTEM <- readRDS(paste0(inputDir,"/results/objectAligned_merge_GTEM.rds"))

## ----eval=FALSE---------------------------------------------------------------
# paciente_merge_aligned <- listAlignedGTEM$alignedObj
# listCoordinates <- listAlignedGTEM$listCoord
# listCoordinatesNew <- listAlignedGTEM$listCoordNew

## ----eval=FALSE---------------------------------------------------------------
# # Do it for each method
# EvaluationGTEM <- STIMA::calculateEvaluation(paciente_merge_aligned, mode = "GTEM",
#                                              listCoordinatesNew, listCoordinates,
#                                              patientType = "unique")

## ----eval=FALSE---------------------------------------------------------------
# # Do it for each method
# STIMA::createDeconvolutionLists(paciente_merge_aligned, mode = "GTEM")

## ----eval=FALSE---------------------------------------------------------------
# # Do it for each method and image
# reference <- readRDS(paste0(inputDir,"/snRNA_QC_Yadav.rds"))
# 
# paciente_merge_aligned_list_im2 <- readRDS("./results/objectAligned_merge_GTEM_list_im2.rds")
# STIMA::deconvolutionRCTD(paciente_merge_aligned_list_im2, reference,
#                          im = "2", mode = "GTEM")
# 

## ----eval=FALSE---------------------------------------------------------------
# listDeconv <- STIMA::deconvolutionRCTD_mergeFiles(modes = c("GTEM","procrustes"), ims = c("2","3","4"))

## ----eval=FALSE---------------------------------------------------------------
# setwd("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab")
# RVcoeff <- STIMA::calculateRVcoeff(modes = c("GTEM","procrustes"), ims = c("2","3","4"))


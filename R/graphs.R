#' graph_cellTypesDeconvolution(listaObjDeconv, cell.types, cell_colors)
#' 
#' This function generates and saves spatial feature plots for cell type deconvolution results.
#' 
#' @param listaObjDeconv A list of Seurat objects containing deconvolution results.
#' @param cell.types A vector of cell type names to be plotted.
#' @param cell_colors A vector of colors corresponding to the cell types.
#' @return A list of plots for each cell type and image.
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#' @import scales
#' @import grDevices
#' @import cowplot
#' @export
graph_cellTypesDeconvolution <- function(listaObjDeconv,
                                         cell.types = c("Astrcytes", "Microglia", "OPC", "Endothelial", 
                                                         "Neurons", "Pericytes", "Schwann", "Lymphocytes", 
                                                         "Oligodendrocytes", "Ependymal Cells", "Meninges"),
                                         cell_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
                                                          "#984EA3", "#FFFF33", "#A65628", "#F781BF",
                                                          "#999999", "#66C2A5", "#FC8D62")
                                        ) {

  if (!dir.exists("./results/graphs/")) {
    dir.create("./results/graphs/", recursive = TRUE)
  }
  saveDir <- "./results/graphs/"
                                         
  listaGraphs <- list() 
  listaPoints <- list()
  patient <- "1"
  modos <- names(listaObjDeconv[[patient]])
  for (i in 1:length(cell.types)) {
    for (mode in modos) {
      for (j in names(listaObjDeconv[[patient]][[mode]])) {
        if ((j == "1" & mode == modos[[1]]) | j != "1") {
          listaGraphs[[j]][[cell.types[[i]]]][[mode]] <- SpatialFeaturePlot(listaObjDeconv[[patient]][[mode]][[j]], 
                                                                            features = cell.types[[i]], 
                                                                            pt.size.factor = 1.2, 
                                                                            crop = FALSE, 
                                                                            alpha = c(0,1)) + 
            scale_fill_gradientn(colors = c("transparent", cell_colors[[i]]), 
                                 limits = c(0, 1), 
                                 oob = scales::squish)
          
          listaPoints[[j]][[cell.types[[i]]]][[mode]] <- SpatialFeaturePlot(listaObjDeconv[[patient]][[mode]][[j]], 
                                                                            features = cell.types[[i]], 
                                                                            pt.size.factor = 1.2, 
                                                                            crop = FALSE, 
                                                                            alpha = c(0,1),
                                                                            image.alpha = 0) + 
            scale_fill_gradientn(colors = c("transparent", cell_colors[[i]]), 
                                 limits = c(0, 1), 
                                 oob = scales::squish)
        }
      }
    }
  }


  listaFinalPlots <- list()

  for (cell_type in cell.types) {
    filas <- list()  # Cada elemento será una fila con 2 gráficos (sin y con histología)
    for (im in names(listaObjDeconv[[patient]][[modos[[1]]]])) {
      for (modo in modos) {
        if (modo == names(listaObjDeconv[[patient]])[[1]]) {
          # Image 1 as the reference 
          p1 <- listaPoints[["1"]][[cell_type]][[modo]] + ggtitle("") + theme_void() + theme(legend.position = "none")
          p2 <- listaGraphs[["1"]][[cell_type]][[modo]] + ggtitle("") + theme_void() + theme(legend.position = "none")
          fila <- p1 + p2
          filas[["ref"]] <- fila
        }
        # Image 3 as the problem
        p2 <- listaPoints[[im]][[cell_type]][[modo]] + ggtitle("") + theme_void() + theme(legend.position = "none")
        p21 <- listaGraphs[[im]][[cell_type]][[modo]] + ggtitle("") + theme_void() + theme(legend.position = "none")
       
        # Una fila con las dos columnas (sin y con histología)
        fila <- p2 + p21
        filas[[modo]] <- fila
      }
    
      # Apilar las filas verticalmente
      panel_completo <- wrap_plots(filas, ncol = 1)
      listaFinalPlots[[cell_type]][[im]] <- panel_completo
    }
  }

  for (cell_type in cell.types) {
    for (im in names(listaObjDeconv[[patient]][[modos[[1]]]])) {
      plot <- listaFinalPlots[[cell_type]][[im]]
      ggsave(paste0(saveDir, "objectAligned_RCTD_SCAnnot_",cell_type,"_",im,".png"), plot = plot, width = 10, height = 30, limitsize = FALSE)
    
    }
  }

  # Save the legend
  legendlist <- list()
  for (i in seq_along(cell.types)) {
    legendlist[[cell.types[[i]]]] <- SpatialFeaturePlot(listaObjAnnot$reference, 
                                                  features = cell.types[[i]], 
                                                  pt.size.factor = 0, 
                                                  crop = FALSE, 
                                                  alpha = c(0,1), image.alpha = 0) + 
      scale_fill_gradientn(colors = c("transparent", cell_colors[[i]]), 
                           limits = c(0, 1), 
                           oob = scales::squish)
  }
  leyenda_combined <- plot_grid(plotlist = legendlist, ncol = 2, align = "v")
  ggsave(paste0(saveDir, "/legend.png"), plot = leyenda_combined)

}

#' graph_evalMetrics(modes = c("GTEM", "procrustes", "RVSSimageJ","PASTE2","STalign"), patientType = c('unique','multiple'))
#' 
#' @param modes A character vector indicating the deconvolution methods to be used. Default is c("GTEM","procrustes","RVSSimageJ","PASTE2","STalign").
#' @param patientType Type of patient data. Options are 'unique' or 'multiple'.
#' @import ggplot2
#' @import ggpubr
#' @import tidyr
#' @import gridExtra
#' @import grid
#' @import lme4
#' @import lmerTest
#' @export
graph_evalMetrics <- function(modes = c("GTEM", "procrustes", "RVSSimageJ","PASTE2","STalign"), patientType = c('unique','multiple'),
                              colors = c("green","red","lightblue","orange","pink","violet","purple","mediumpurple","purple4")) {
  patientType <- match.arg(patientType)
  if (!dir.exists("./results/graphs/")) {
    dir.create("./results/graphs/", recursive = TRUE)
  }
  saveDir <- "./results/graphs/"

  dataDir <- "./results/"
  typeComparison <- "max_pixels"
  listaEvalFiles <- setNames(vector("list", length(modes)), modes)
  for (mode in modes) { 
    listaEvalFiles[["mode"]] <- readRDS(paste0(dataDir, "evaluation_merge_",mode,".rds"))
  }
  RV_table <- readRDS(paste0(dataDir, "objectAligned_RV_merge.rds"))

  metrics <- c("1 MSE", "2 MSE Gray", "3 SSIM", "4 Euclidean", "5 RV")
  
  listaDataFrames <- list()
  listaStats <- list()
  listaControl <- data.frame()
  datosEucl <- data.frame()
  datosRV <- data.frame()

  for (elem in seq_along(listaEvalFiles)) {
    Evaluation <- listaEvalFiles[[elem]]
    RVvalue <- RV_table[[elem]]
    group <- names(listaEvalFiles)[[elem]]
    
    # Initialize vectors to store MSE, SSIM, and Euclidean values for original and transformed images
    MSE_controlPos <- c()
    MSE_controlNeg <- c()
    MSE_controlMov <- c()
    MSE_original <- c()
    MSE_transformado <- c()
    MSE_gray_controlPos <- c()
    MSE_gray_controlNeg <- c()
    MSE_gray_controlMov <- c()
    MSE_gray_original <- c()
    MSE_gray_transformado <- c()
    SSIM_controlPos <- c()
    SSIM_controlNeg <- c()
    SSIM_controlMov <- c()
    SSIM_original <- c()
    SSIM_transformado <- c()
    Eucl_original_mean <- c()
    Eucl_transformado_mean <- c()
    Eucl_original <- c()
    Eucl_transformado <- c()
    RV_original <- c()
    RV_transformado <- c()
    Image <- c()
    
    # Loop through the first three images to extract MSE, SSIM, Euclidean and RV values
    for (i in 1:(length(Evaluation$control) - 1)) {
      MSE_original <- c(MSE_original, Evaluation$original[[typeComparison]][[i]]$mse_value)
      MSE_transformado <- c(MSE_transformado, Evaluation$transformed[[typeComparison]][[i]]$mse_value)
      MSE_gray_original <- c(MSE_gray_original, Evaluation$original[[typeComparison]][[i]]$mse_gray_value)
      MSE_gray_transformado <- c(MSE_gray_transformado, Evaluation$transformed[[typeComparison]][[i]]$mse_gray_value)
      SSIM_original <- c(SSIM_original, Evaluation$original[[typeComparison]][[i]]$ssim_value)
      SSIM_transformado <- c(SSIM_transformado, Evaluation$transformed[[typeComparison]][[i]]$ssim_value)
      Eucl_original_mean <- c(Eucl_original_mean, mean(Evaluation$original[[typeComparison]][[i]]$Eucl_value))
      Eucl_transformado_mean <- c(Eucl_transformado_mean, mean(Evaluation$transformed[[typeComparison]][[i]]$Eucl_value))
      Eucl_original <- c(Eucl_original, Evaluation$original[[typeComparison]][[i]]$Eucl_value)
      Eucl_transformado <- c(Eucl_transformado, Evaluation$transformed[[typeComparison]][[i]]$Eucl_value)
      RV_original <- c(RV_original, RVvalue[[i]][2,"RV"])
      RV_transformado <- c(RV_transformado, RVvalue[[i]][8,"RV"])
    }
    
    # Loop through the control data to extract MSE and SSIM values
    for (i in 2:length(Evaluation$control)) {
      MSE_controlPos <- c(MSE_controlPos, Evaluation$control[[i]]$`Positive control solution`$mse_value)
      MSE_controlNeg <- c(MSE_controlNeg, Evaluation$control[[i]]$`Negative control solution`$mse_value)
      MSE_controlMov <- c(MSE_controlMov, Evaluation$control[[i]]$`Movement control solution`$mse_value)
      MSE_gray_controlPos <- c(MSE_gray_controlPos, Evaluation$control[[i]]$`Positive control solution`$mse_gray_value)
      MSE_gray_controlNeg <- c(MSE_gray_controlNeg, Evaluation$control[[i]]$`Negative control solution`$mse_gray_value)
      MSE_gray_controlMov <- c(MSE_gray_controlMov, Evaluation$control[[i]]$`Movement control solution`$mse_gray_value)
      SSIM_controlPos <- c(SSIM_controlPos, Evaluation$control[[i]]$`Positive control solution`$ssim_value)
      SSIM_controlNeg <- c(SSIM_controlNeg, Evaluation$control[[i]]$`Negative control solution`$ssim_value)
      SSIM_controlMov <- c(SSIM_controlMov, Evaluation$control[[i]]$`Movement control solution`$ssim_value)

      Image <- c(Image, paste0("Image",i))
    }
    
    data_control <- data.frame(
      Image = Image,
      Group = group,
      MSE_controlPos = MSE_controlPos,
      MSE_controlNeg = MSE_controlNeg,
      MSE_controlMov = MSE_controlMov,
      MSE_gray_controlPos = MSE_gray_controlPos,
      MSE_gray_controlNeg = MSE_gray_controlNeg,
      MSE_gray_controlMov = MSE_gray_controlMov,
      SSIM_controlPos = SSIM_controlPos,
      SSIM_controlNeg = SSIM_controlNeg,
      SSIM_controlMov = SSIM_controlMov,
      Eucl_controlPos = 0,
      Eucl_controlNeg = NaN,
      Eucl_controlMov = NaN,
      RV_controlPos = 1,
      RV_controlNeg = 0,
      RV_controlMov = NaN
    )
    
    # Create a data frame to hold the values from the same region and same point
    data_comparison <- data.frame(
      Image = Image,
      MSE_original     = MSE_original,
      MSE_transformado = MSE_transformado,
      MSE_gray_original     = MSE_gray_original,
      MSE_gray_transformado = MSE_gray_transformado,
      SSIM_original     = SSIM_original,
      SSIM_transformado = SSIM_transformado,
      Eucl_original     = Eucl_original_mean,
      Eucl_transformado = Eucl_transformado_mean,
      RV_original     = RV_original,
      RV_transformado = RV_transformado
    )
    
    # Create a data frame to group the Euclidean distance for each mode, image and landmark
    landmark_id <- paste0("lmk",seq_along(Evaluation$original[[typeComparison]][[i]]$Eucl_value))
    df_temp <- data.frame(
      image_id = rep(Image, each = length(landmark_id)),
      method = group,
      landmark_id = landmark_id,
      control = Eucl_original,
      transformed = Eucl_transformado
    )
    datosEucl <- rbind(datosEucl, df_temp)
    
    
    df_temp2 <- data.frame(
      image_id = Image,
      method = group,
      control = RV_original,
      transformed = RV_transformado
    )
    datosRV <- rbind(datosRV, df_temp2)
    
    # Set row names for the data frame
    rownames(data_control) <- Image
    rownames(data_comparison) <- Image
    
    listaDataFrames[[names(listaEvalFiles[elem])]] <- data_comparison
    
    
    # Calculate the mean and standard deviation
    mean_values <- colMeans(data_comparison[, -1])  # Ignorar la columna "Image"
    sd_values <- apply(data_comparison[, -1], 2, sd)  # Ignorar la columna "Image"
    error_values <- sapply(sd_values, function(x) {x/sqrt(3)})
    
    # Create a new dataframe for statistics
    stats_df <- data.frame(
      Parameter = c("Original","Transformed"),
      Mean = mean_values,
      SD = sd_values,
      error = error_values,
      group = c("0 original", paste0(elem, " ", group)),
      metric = rep(metrics, each = 2)
    )
    
    # Perform t-tests
    mse_test <- t.test(data_comparison$MSE_original,
                      data_comparison$MSE_transformado,
                      paired = TRUE)
    mse_p_value <- mse_test$p.value
    
    mse_gray_test <- t.test(data_comparison$MSE_gray_original,
                            data_comparison$MSE_gray_transformado,
                            paired = TRUE)
    mse_gray_p_value <- mse_gray_test$p.value
    
    ssim_test <- t.test(data_comparison$SSIM_original,
                        data_comparison$SSIM_transformado,
                        paired = TRUE)
    ssim_p_value <- ssim_test$p.value
    
    eucl_p_value <- NaN
    
    RV_p_value <- NaN
    
    # Groupt the t.test results
    p_values_df <- data.frame(
      Metric = metrics,
      P_value = c(mse_p_value, mse_gray_p_value, ssim_p_value, eucl_p_value, RV_p_value)
    )
    
    stats_df$p_value <- NaN
    stats_df$p_value[c(2, 4, 6, 8, 10)] <- c(mse_p_value, mse_gray_p_value, ssim_p_value, eucl_p_value, RV_p_value)

    listaStats[[names(listaEvalFiles[elem])]] <- stats_df
    listaControl <- rbind(listaControl,data_control)
  }

  # Add the control parameters

  # Calculate the mean and standard deviation
  mean_values <- colMeans(listaControl[,c(-1,-2)])  # Ignorar la columna "Image"
  sd_values <- apply(listaControl[,c(-1,-2)], 2, sd)  # Ignorar la columna "Image"
  error_values <- sapply(sd_values, function(x) {x/sqrt(3)})

  # "Positive Control","Negative Control","Movement Control",
  # Create a new dataframe for statistics
  stats_df <- data.frame(
    Parameter = c("Transformed"),
    Mean = mean_values,
    SD = sd_values,
    error = error_values,
    group = c("0 1 positive","0 2 negative","0 3 movement"),
    metric = rep(metrics, each = 3),
    p_value = NaN
  )

  listaStats[["control"]] <- stats_df

  datosEucl_long <- pivot_longer(
    datosEucl,
    cols = c(control, transformed),
    names_to  = "condition",
    values_to = "distance"
  )
  # ‘condition’ and ‘method’ as factors
  datosEucl_long$condition <- factor(datosEucl_long$condition, levels = c("control","transformed"))
  datosEucl_long$method    <- factor(datosEucl_long$method)

  # Linear Mixed Model
  # distance ~ condition + (1 | image_id) 
  #      (1|image) variability between images
  #      (1 | image_id:landmark_id) variability by landmarks in each image

  resultsEucl <- lapply(levels(datosEucl_long$method), function(method_i) {
    dd <- subset(datosEucl_long, method == method_i)
    mod <- lmer(
      distance ~ condition + (1 | image_id),
      data = dd)
    # p‑value for ‘condition’ with Anova 
    pcond <- anova(mod)["condition", "Pr(>F)"]
    coef  <- summary(mod)$coefficients["conditiontransformed","Estimate"]
    data.frame(
      method      = method_i,
      estimate    = coef,
      p_value     = pcond
    )
  })
  resultsEucl <- do.call(rbind, resultsEucl)

  for (metric in metrics) {
    listaStats[[metric]][["Eucl_transformado","p_value"]] <- resultsEucl[resultsEucl$method == metric,"p_value"]
  }


  datosRV_long <- pivot_longer(
    datosRV,
    cols = c(control, transformed),
    names_to  = "condition",
    values_to = "RV"
  )
  # ‘condition’ and ‘method’ as factors
  datosRV_long$condition <- factor(datosRV_long$condition, levels = c("control","transformed"))
  datosRV_long$method    <- factor(datosRV_long$method)

  # Linear Mixed Model
  # distance ~ condition + (1 | image_id) 
  #      (1|image) variability between images
  #      (1 | image_id:landmark_id) variability by landmarks in each image

  resultsRV <- lapply(levels(datosRV_long$method), function(method_i) {
    dd <- subset(datosRV_long, method == method_i)
    mod <- lmer(
      RV ~ condition + (1 | image_id),
      data = dd)
    # p‑value for ‘condition’ with Anova 
    pcond <- anova(mod)["condition", "Pr(>F)"]
    coef  <- summary(mod)$coefficients["conditiontransformed","Estimate"]
    data.frame(
      method      = method_i,
      estimate    = coef,
      p_value     = pcond
    )
  })
  resultsRV <- do.call(rbind, resultsRV)

  for (metric in metrics) {
    listaStats[[metric]][["RV_transformado","p_value"]] <- resultsRV[resultsRV$method == metric,"p_value"]
  }


  # Raw data for boxplots
  metrics_list <- c("MSE", "MSEgray", "SSIM", "Eucl", "RV")

  data_list <- lapply(metrics_list, function(metric) {
    # 'listaControl'
    control_data <- unlist(listaControl[, grep(metric, colnames(listaControl))], use.names = FALSE)
    # 'listaDataFrames'
    model_data <- do.call(c, lapply(listaDataFrames, function(df) df[, grep(metric, colnames(df))]))

    c(control_data, model_data)
  })

  # From list to dataframe
  data <- data.frame(matrix(unlist(data_list), ncol = length(metrics_list), byrow = FALSE))
  colnames(data) <- metrics_list

  # Add 'image' and 'metric'
  images_control <- rep(rownames(listaControl), each = 3)
  metrics_control <- rep(c("0 1 positive", "0 2 negative", "0 3 movement"), length.out = nrow(listaControl))

  images_models <- rep(names(listaDataFrames), sapply(listaDataFrames, nrow))
  metrics_models <- metrics

  data$image <- c(images_control, images_models)
  data$metric <- as.factor(c(metrics_control, metrics_models))
  data$image <- as.factor(data$image)

  # Boxplots graphs
  listGraphs <- list()
  label_list <- c("Positive Control","Negative Control","Original",modes)
  for (metric in metrics) {

    p_values <- sapply(modes, function(mode) {
      val <- listaStats[[mode]]$p_value[[metrics[metric]]]
      if (is.nan(val)) 1 else val
    })

    max_mean <- max(sapply(modes, function(mode) {
      listaStats[[mode]]$Mean[[metrics[metric] - 1]]
    }))
    max_error <- max(sapply(modes, function(mode) {
      listaStats[[mode]]$error[[metrics[metric] - 1]]
    }))

    # Define plot limits and y positions
    if (patientType == "unique") {
      maxY <- ifelse(metric == "SSIM" || metric == "RV", 1, 2500)
      y_positions <- seq(0.05 * maxY, 0.45 * maxY, length.out = length(p_values))
    } else if (patientType == "multiple") {
      maxY <- ifelse(metric == "SSIM" || metric == "RV", 1, 35000)
      y_positions <- seq(0.45 * maxY, 0.85 * maxY, length.out = length(p_values))
    }

    # Generate Plot
    plot <- ggplot(data, aes_string(x = "metric", y = metric)) +
      geom_boxplot(aes(fill = metric), outlier.shape = NA) +
      geom_jitter(aes(fill = image), width = 0.1, shape = 21, size = 3, stroke = 0.2) +
      scale_fill_manual(values = colors, labels = label_list) +
      labs(y = "Mean value", x = "", title = paste(metric, "- Pixel similarity")) +
      coord_cartesian(ylim = c(0, maxY)) +
      theme_minimal() +
      theme(plot.title = element_text(size = 42, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 0, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 34),
            axis.title.y = element_text(size = 34),
            legend.position = "none") +
      mapply(function(comp, pval, ypos) {
        geom_signif(comparisons = list(comp), annotations = sprintf("%.4f", pval), y_position = ypos, 
                    tip_length = 0.02 / maxY, textsize = 10, color = ifelse(pval <= 0.05, "red", "black"))
      }, comparisons, p_values, y_positions)

    # Save Plot
    ggsave(paste0(saveDir, "STIMA_merge_", metric, "_boxplot.png"), plot = plot)

    listaGraphs[[metric]] <- plot
  }

  combined_plot <- grid.arrange(
    unlist(listaGraphs),
    ncol = 1, 
    heights = c(1, 1, 1, 1, 1.6)
  )

  # Guarda la figura combinada
  ggsave(paste0(saveDir, "STIMA_combined_boxplot_5.png"), plot = combined_plot, width = 11, height = 40)

}
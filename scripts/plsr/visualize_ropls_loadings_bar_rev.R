
#' Barplot for loadings
#'
#' @param model ropls object
#' @param ind_LV latent variable index
#'
#' @return
#' @export
#'
visualize_ropls_loadings_bar_rev <- function(model, options = list()) {
  
  # ----------------- BEGIN OPTIONS ----------------- #
  if (!("LV_ind" %in% names(options))) {
    options$LV_ind <- c(1)
  }
  if ("y" %in% names(options)) {
    y <- options$y
    if (is.factor(y)) {
      n_groups <- nlevels(y)
    } else {
      n_groups <- NA
    }
  } else {
    n_groups <- NA
  }
  if (!("mark_enrichment" %in% names(options)) | is.na(n_groups)) {
    options$mark_enrichment <- FALSE
  }
  if (options$mark_enrichment & (is.na(n_groups) | !("X" %in% names(options))))  {
    stop("Enrichment only works for classification and when X and y are provided")
  }
  
  # color for the scores and name of the grouping
  if (!("y_name" %in% names(options))) {
    y_name <- "y"
  } else {
    y_name <- options$y_name
  }
  # color for the scores and name of the grouping
  if (!("colors" %in% names(options)) | length(grep(y_name, names(options$colors))) == 0) {
    if (is.factor(y)) {
      tmp <- rep(NA, length = nlevels(y))
      names(tmp) <- levels(y)
      for (ind in 1:nlevels(y)) {
        tmp[ind] <- RColorBrewer::brewer.pal(n = max(3, nlevels(y)), name = 'Dark2')[ind]
      }
      options$colors <- list()
      options$colors[[y_name]] <- tmp
    } else {
      # For regression, a color palette needs to be provided
      options$colors$y <- list(low = "#C7E4F9", high = "#004D7F")
    }
  }
  
  if (ropls::getSummaryDF(model)$pre +
      ropls::getSummaryDF(model)$ort < options$LV_ind) {
    stop("required LV exceed existing LVs")
  }
  
  if (!("df_features" %in% names(options))) {
    options$df_features <- data.frame(name = rownames(model@loadingMN),
                                      label = rownames(model@loadingMN))
  }
  
  # ----------------- END OPTIONS ----------------- #
  
  # check first whether its a orthogonal PLS or a regular PLS
  if (ropls::getSummaryDF(model)$ort > 0) {
    stop("orthogonal PLS-DA not supported yet")
    # if (options$LV_ind[1] == 1) {
    #   df_loadings <- data.frame(LV1 = ropls::getLoadingMN(model),
    #                             LV2 = ropls::getLoadingMN(model, orthoL = TRUE)[,options$LV_ind[2] - 1])
    # } else {
    #   df_loadings <- data.frame(LV1 = ropls::getLoadingMN(model, orthoL = TRUE)[, options$LV_ind[1] - 1],
    #                             LV2 = ropls::getLoadingMN(model, orthoL = TRUE)[, options$LV_ind[2] - 1])
    # }
  } else {
    df_loadings <- data.frame(LV = ropls::getLoadingMN(model)[,options$LV_ind[1]],
                              vip_scores = ropls::getVipVn(model))
    df_loadings$features <- rownames(df_loadings)
    df_loadings$labels <- options$df_features$label[match(rownames(df_loadings), options$df_features$name)]
  }
  
  
  
  # TODO: catch if its an orthogonal
  
  if (options$mark_enrichment & !is.na(n_groups)) {
    df_loadings$mark <- NA
    X <- options$X
    
    for (ind_feat in 1:nrow(df_loadings)) {
      tmp_mean <- rep(NA, length = nlevels(y))
      for (ind_class in 1:nlevels(y)) {
        tmp_mean[ind_class] <- mean(X[which(y == levels(y)[ind_class]),
                                      which(colnames(X) == df_loadings$features[ind_feat])])
      }
      # change color for loadings
      #1. using median
      df_loadings$mark[ind_feat] <- levels(y)[which.max(tmp_mean)]
      #2. using loading direction
      # if (df_loadings$LV[ind_feat]<0){
      #   df_loadings$mark[ind_feat] <-"p"
      # }else{
      #   df_loadings$mark[ind_feat] <-"b"
      # }
      
    }
    df_loadings$mark  <- factor(df_loadings$mark, levels = levels(y))
  }
  
  # for plsr
  if (options$mark_enrichment & is.na(n_groups)) {
    df_loadings$mark <- NA
    X <- options$X
    
    for (ind_feat in 1:nrow(df_loadings)) {
      
      if(df_loadings$LV[ind_feat]>0|df_loadings$LV[ind_feat]==0){
        df_loadings$mark[ind_feat] <- "high"
        
      }else{
        df_loadings$mark[ind_feat] <- "low"
        
      }

    }
    df_loadings$mark  <- factor(df_loadings$mark, levels = c("low", "high"))
  }
  
  
  df_loadings <- df_loadings[order(df_loadings$vip_scores), ]
  df_loadings$features <- factor(df_loadings$features, levels = unique(df_loadings$features))
  
  # df_loadings <- df_loadings[df_loadings$vip_scores>1, ]
  
  
  
  # plot loadings sorted according to the VIP score and color coding it
  # according to enrichent in classes
  if (options$mark_enrichment) {
    plt_bar <- ggplot2::ggplot(data = df_loadings, ggplot2::aes(x = features, y = LV, fill = mark)) +
      ggplot2::scale_fill_manual(values = options$colors[[y_name]])
  } else {
    plt_bar <- ggplot2::ggplot(data = df_loadings, ggplot2::aes(x = features, y = LV))
  }
  plt_bar <- plt_bar +
    ggplot2::geom_bar(stat = "identity", color = "black") +
    ggplot2::scale_fill_manual(values = c("red3", "royalblue3"), labels = c("naive", "SIVmac239∆GY"))+
    ggplot2::coord_flip() +
    ggplot2::xlab("") +
    ggplot2::ylab(paste("LV", options$LV_ind, " loadings", sep = "")) +
    ggplot2::labs(fill = "enriched in") +
    ggplot2::scale_x_discrete(labels = df_loadings$labels) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=16,color = "black"), 
                   axis.text.y = ggplot2::element_text(size=16,color = "black"),
                   axis.title.x = ggplot2::element_text(size=16,color = "black"),
                   axis.title.y = ggplot2::element_text(size=16,color = "black"),
                   legend.key.size = unit(1, 'cm'),
                   legend.title = ggplot2::element_text(size=16,color = "black"),
                   legend.text = ggplot2::element_text(size=16,color = "black")
                   )#,
  #axis.text.y = element_text(colour = as.character(feature_annot$useColor[match(dfBar$features[order(dfBar$vipScores)],
  # rownames(feature_annot))])))
  
}
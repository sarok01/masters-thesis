# data: dataframe 
# colresponse: (int) column for response variable (e.g. viral load)
# colstart: (int) column from which system serology are recorded
# thresh: (0~1) threshold for lasso feature selection (1 means features will be selected if they were selected from 100/100 rounds of lasso, 0.9 means 90/100 rounds)
# alphal: 0: ridge, 1: lasso, 0~1: elastic net 
# nfold: (int) number of folds during lasso selection (usually 5)
# cvtry: (int) number of cross-validation (usually 10~100)
# foldername: (string) name of the folder where output figures will be 
# title: (string) name of figure files that will be saved

plsr_novalidation<-function(data, colresponse, colstart, thresh, alphal, nfold, cvtry, foldername, title){

  source("select_lasso_rev.R")
  source("select_repeat_rev.R")
  source("visualize_ropls_loadings_bar_rev.R")
  source("visualize_validate_largefont.R")
  source("correlation_network_rev.R")
  source("visualize_ropls_revs.R")
  dir.create(foldername)
  setwd(foldername)  
  
  
  responsename<-colnames(data)[colresponse]
  dataname <- comment(data) 
  
  X <- as.matrix(as.data.frame(lapply(data[, colstart:(ncol(data))], as.numeric)))
  
  # zscore data
  X <- scale(X, center = TRUE, scale = TRUE)
  y <- data[,colresponse]
  print(responsename)
  
  if (length(which(is.na(y))>0)){

    X<-X[-which(is.na(y)),]
    y<-y[-which(is.na(y))]
  }
  
  # select features by lasso
  opts_sel <- list(n_trials = 100, threshold =thresh , alpha=alphal, subfolds=nfold, return_count = FALSE)
  sel_features <- select_repeat_rev(X, y, selector = select_lasso_rev, options = opts_sel)
  print(sel_features)
  
  if(length(sel_features)>1){     
  X_sel <- X[,sel_features]

  # build plsr model
  model <- train_ropls(X_sel, y, options  = list(n_LV = 2)) # MODEL MODEL MODEL MODEL MODEL
  print(model)
  # opts_plot for score
  df_features <- data.frame(name = colnames(X))
  df_features$label <- df_features$name
  opts_plot <- list(df_features = df_features,
                    loading_alpha = 1,# transparency for the loadings
                    score_alpha = 1,# transparency for the scores
                    LV_ind = c(1,2), # which LVs to plot
                    # colors = my_colors,
                    size =4,
                    y_name = responsename) 
  
  # plot score plot
  plt_scores <- visualize_ropls_scores_revs(model, y, options = opts_plot)
  print(plt_scores)
  ggsave(filename = paste0("score_", title, ".pdf"), plot = plt_scores, height=5, width=5, device="pdf")
  
  # opts plot for model
  opts_plot$X <- X
  opts_plot$y <- y
  opts_plot$LV_ind <- 1
  opts_plot$mark_enrichment <- T
  plt_loadings_bar <- visualize_ropls_loadings_bar_rev(model, options = opts_plot)
  print(plt_loadings_bar)
  ggsave(filename = paste0("loading_", title, ".pdf"), plot = plt_loadings_bar, height=9, width=15, device="pdf")
  
  opts_sel <- list(n_trials = 100, threshold = thresh, alpha=alphal, subfolds=nfold,  return_count = FALSE)
  select <- function(X, y) { return(select_repeat_rev(X, y, selector = select_lasso_rev, options = opts_sel)) }
  method = list(select=select,
                train = train_ropls,
                predict = predict_ropls,
                score = score_accuracy) #score_r2
  # change rf_trials = 10 or 100 and pt_trials = 10 or 100 for validation
  opts = list(n_folds = 5, rf_trials = 5, pt_trials = 5, save = TRUE)
  
  # change cvtry = 100 or 10 for validation
  return_cv <- validate_repeat(X, y, method, opts, n_trials=cvtry)
  #saveRDS(return_cv, file = "return_cv.RDS")
  cv_all=unlist(return_cv)
  
  print(mean(cv_all))
  
  # plot accuracy
  plt_val<-visualize_validate_largefont(return_cv)
  print(plt_val)
  ggsave(filename = paste0("validation_", title, ".pdf"), plot = plt_val, height=9, width=7, device="pdf")

  # plot co-correlate network
  #opts<-list(FDR=0.005)
  #co<-correlation_network_rev(X, sel_features, opts)
  #ggsave(filename = paste0("correlation_", title, ".pdf"), plot = co, height=10, width=10, dpi =300, device="pdf")

  }else{
    
    print("only one feature selected by regularization")
  }

  
}
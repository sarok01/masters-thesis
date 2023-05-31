library(circlize)
library(fields)
library(RColorBrewer)
library(pheatmap)
library(dendsort)
# library(bootcorci)
# library(parallel)
# library(psych)

# correlation differences -------

correlation_matrix_compare <- function(
  X1, X2, FDR=0.05, idx1 = NULL, idx2 = NULL, p.adjust.method = "BH",
  filter_vals = FALSE, star = "*", bootstrap = F,
  bootstrap.method = c("lapply", "for", "mclapply"), mc.cores = NA, fig=i, ...) {

  # X1 and X2 are two different groups with the same measurements taken in each group.
  # This function computes and plots the Fisher r-to-z transformation for the
  # different in correlation between the two groups (for all pairs of features).
  # Statistically significant differences (after multiple-hypothesis correction)
  # are marked with a '*', at the specified FDR level. By default, significance
  # is determined by bootstrapping methods. The function also includes
  # options for doing computations in parallel across multiple cores.

  if (ncol(X1) != ncol(X2)) {
    stop("X1 and X2 must have same number of variables")
  }

  if (any(colnames(X1) != colnames(X2))) {
    stop("X1 and X2 must have same variable names")
  }

  # variables to correlate with each other
  if (is.null(idx1)) {
    idx1 <- 1:ncol(X1)
  }
  if (is.null(idx2)) {
    idx2 <- 1:ncol(X1)
  }

  # browser()
  # compute correlation matrices and Fisher transformation
  r1 <- as.matrix(X1 %>% correlate(method = "spearman", diagonal=1) %>% column_to_rownames("term"))
  z1 <- pmin(pmax(0.5 * log((1+r1) / (1-r1)), -4), 4)
  r2 <- as.matrix(X2 %>% correlate(method = "spearman", diagonal=1) %>% column_to_rownames("term"))
  z2 <- pmin(pmax(0.5 * log((1+r2) / (1-r2)), -4), 4)

  # browser()
  # compute fisher Z test statistic
  z <- (z1 - z2) / sqrt(1/(nrow(X1) - 3) + 1/(nrow(X2) - 3))
  z <- pmin(pmax(z, -4), 4)
  p_fisher_all <- 2*pmin(pnorm(z), pnorm(-z))
  if (length(which(is.na(p_fisher_all))) > 0) {
    stop("NAs are not expected here")
  }
  # select pairs of variables of interest
  p_fisher <- p_fisher_all
  p_fisher[,] <- NA
  p_fisher[idx1,idx2] <- p_fisher_all[idx1,idx2]
  p_fisher[idx2,idx1] <- p_fisher_all[idx2,idx1]
  diag(p_fisher) <- NA
  # multiple-hypothesis correction
  flat_fisher <- as.numeric(p_fisher)
  print(sum(!is.na(flat_fisher)))
  flat_fisher <- p.adjust(flat_fisher, method = p.adjust.method)
  p_fisher.adj <- matrix(flat_fisher, ncol = ncol(p_fisher),
                     dimnames = dimnames(p_fisher))
  
  # z <-z[which(colnames(z)=="FcgR2A_HPV16"|colnames(z)=="FcgR3A_HPV16"|colnames(z)=="FcgR2B_HPV16"|colnames(z)=="FcgR3B_HPV16"|colnames(z)=="totalIgG_HPV18"|colnames(z)=="totalIgG_HPV16"|colnames(z)=="FcgR3A_HPV18"|colnames(z)=="FcgR3B_HPV18"|colnames(z)=="FcgR2A_HPV18"|colnames(z)=="FcgR2B_HPV18"),which(colnames(z)=="FcgR2A_HPV16"|colnames(z)=="FcgR3A_HPV16"|colnames(z)=="FcgR2B_HPV16"|colnames(z)=="FcgR3B_HPV16"|colnames(z)=="totalIgG_HPV18"|colnames(z)=="totalIgG_HPV16"|colnames(z)=="FcgR3A_HPV18"|colnames(z)=="FcgR3B_HPV18"|colnames(z)=="FcgR2A_HPV18"|colnames(z)=="FcgR2B_HPV18")]
  # 
  # p_fisher.adj <-p_fisher.adj[which(colnames(p_fisher.adj)=="FcgR2A_HPV16"|colnames(p_fisher.adj)=="FcgR3A_HPV16"|colnames(p_fisher.adj)=="FcgR2B_HPV16"|colnames(p_fisher.adj)=="FcgR3B_HPV16"|colnames(p_fisher.adj)=="totalIgG_HPV18"|colnames(p_fisher.adj)=="totalIgG_HPV16"|colnames(p_fisher.adj)=="IgA2_HPV16"|colnames(p_fisher.adj)=="IgA2_HPV18"|colnames(p_fisher.adj)=="FcgR3A_HPV18"|colnames(p_fisher.adj)=="FcgR3B_HPV18"|colnames(p_fisher.adj)=="FcgR2A_HPV18"|colnames(p_fisher.adj)=="FcgR2B_HPV18"),which(colnames(p_fisher.adj)=="FcgR2A_HPV16"|colnames(p_fisher.adj)=="FcgR3A_HPV16"|colnames(p_fisher.adj)=="FcgR2B_HPV16"|colnames(p_fisher.adj)=="FcgR3B_HPV16"|colnames(p_fisher.adj)=="totalIgG_HPV18"|colnames(p_fisher.adj)=="totalIgG_HPV16"|colnames(p_fisher.adj)=="IgA2_HPV16"|colnames(p_fisher.adj)=="IgA2_HPV18"|colnames(p_fisher.adj)=="FcgR3A_HPV18"|colnames(p_fisher.adj)=="FcgR3B_HPV18"|colnames(p_fisher.adj)=="FcgR2A_HPV18"|colnames(p_fisher.adj)=="FcgR2B_HPV18")]

  # bootstrap testing
  bootstrapper <- function(idx) {
    i <- idx[1]
    j <- idx[2]
    if ((i > j) && (((i %in% idx1) && (j %in% idx2)) || ((j %in% idx1) && (i %in% idx2)))) {
      print(c(i,j))
      tmp <- twocorci(X1[,i], X1[,j], X2[,i], X2[,j], method = "spearman")
      return(tmp$p.value)
    } else {
      return(NA)
    }
  }
  p_bootstrap <- NULL
  p_bootstrap.adj <- NULL

  # browser()
  if (bootstrap) {
    if (bootstrap.method == "for") {
      p_bootstrap <- matrix(NA, nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
      for (i in 1:nrow(p_bootstrap)) {
        for (j in 1:ncol(p_bootstrap)) {
          print(c(i, j))
          tmp <- twocorci(X1[,i], X1[,j], X2[,i], X2[,j], method = "spearman")
          p_bootstrap[i, j] <- tmp$p.value
        }
      }
    } else if (bootstrap.method == "lapply") {
      idxs <- split(as.matrix(unname(expand.grid(1:nrow(z), 1:ncol(z)))),
                    seq(length(z)))
      # results in out2 don't match out...
      p_bootstrap <- matrix(as.numeric(lapply(idxs, bootstrapper)),
                            nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
    } else if (bootstrap.method == "mclapply") {
      idxs <- split(as.matrix(unname(expand.grid(1:nrow(z), 1:ncol(z)))),
                    seq(length(z)))
      # results in out2 don't match out...
      if (is.na(mc.cores)) {
        mc.cores <- detectCores() - 1
        print(mc.cores)
      }
      p_bootstrap <- matrix(as.numeric(mclapply(idxs, bootstrapper, mc.cores = mc.cores)),
                            nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
    }
    p_bootstrap[is.na(p_bootstrap) & !is.na(t(p_bootstrap))] <- t(p_bootstrap)[is.na(p_bootstrap) & !is.na(t(p_bootstrap))]
    flat_bootstrap <- as.numeric(p_bootstrap)
    print(sum(!is.na(flat_bootstrap)))
    flat_bootstrap <- p.adjust(flat_bootstrap, method = p.adjust.method)
    p_bootstrap.adj <- matrix(flat_bootstrap, ncol = ncol(p_bootstrap),
                              dimnames = dimnames(p_bootstrap))
    display_text <- ifelse(!is.na(p_bootstrap.adj) & p_bootstrap.adj < FDR, star, "")
  } else {
    display_text <- ifelse(!is.na(p_fisher.adj) & p_fisher.adj < FDR, star, "")
  }

  min_idx <- round((min(z) + max(abs(z))) / (2*max(abs(z))) * 99)+1
  max_idx <- round((max(z) + max(abs(z))) / (2*max(abs(z))) * 99)+1
  color <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)[min_idx:max_idx]

  callback = function(hc, mat) {
    # View(mat)
    dend = reorder(as.dendrogram(hc), -rowMeans(mat))
    as.hclust(dend)
  }

  # out <- pheatmap(z, color = color, display_numbers = display_text,
  #          clustering_callback = callback, ...)

  # show_rownames = F, show_colnames = F, legend = FALSE, 
  # out<-pheatmap(z, color = colorRampPalette(rev(brewer.pal(n = 7, name = "PuOr")))(100), breaks=seq(-4, 4, length.out=100),display_numbers = display_text, show_rownames = T, show_colnames = T, legend = F, fontsize_number = 70, fontsize_row=40, fontsize_col=40,cluster_rows =FALSE, cluster_cols=FALSE)

  z[] 
  # figurename = paste("C:/Users/WJ/Documents/research/A_immunology/A_melanoma/f", toString(fig), ".tif")
  # ggsave(plot = out, height = 20, width = 20, dpi = 100,
  #        filename = figurename, device='tiff')
  
  # corrdiff<-list(r1 = r1, r2 = r2, z1 = z1, z2 = z2, z = z,
  #      p_fisher = p_fisher, p_fisher.adj = p_fisher.adj,
  #      p_bootstrap = p_bootstrap, p_bootstrap.adj = p_bootstrap.adj,
  #      pheatmap.res = out)
  
  corrdiff<-list(r1 = r1, r2 = r2, z1 = z1, z2 = z2, z = z,
                 p_fisher = p_fisher, p_fisher.adj = p_fisher.adj,
                 p_bootstrap = p_bootstrap, p_bootstrap.adj = p_bootstrap.adj)
  
  return(corrdiff)
}


correlation_matrix_compare_circle <- function(
  X1, X2, FDR=0.05, idx1 = NULL, idx2 = NULL, p.adjust.method = "BH",
  filter_vals = FALSE, star = "*", bootstrap = T,
  bootstrap.method = c("lapply", "for", "mclapply"), mc.cores = NA, fig=i, ...) {
  
  # X1 and X2 are two different groups with the same measurements taken in each group.
  # This function computes and plots the Fisher r-to-z transformation for the
  # different in correlation between the two groups (for all pairs of features).
  # Statistically significant differences (after multiple-hypothesis correction)
  # are marked with a '*', at the specified FDR level. By default, significance
  # is determined by bootstrapping methods. The function also includes
  # options for doing computations in parallel across multiple cores.
  
  if (ncol(X1) != ncol(X2)) {
    stop("X1 and X2 must have same number of variables")
  }
  
  if (any(colnames(X1) != colnames(X2))) {
    stop("X1 and X2 must have same variable names")
  }
  
  # variables to correlate with each other
  if (is.null(idx1)) {
    idx1 <- 1:ncol(X1)
  }
  if (is.null(idx2)) {
    idx2 <- 1:ncol(X1)
  }
  
  # browser()
  # compute correlation matrices and Fisher transformation
  r1 <- as.matrix(X1 %>% correlate(method = "spearman", diagonal=1) %>% column_to_rownames("term"))
  z1 <- pmin(pmax(0.5 * log((1+r1) / (1-r1)), -4), 4)
  r2 <- as.matrix(X2 %>% correlate(method = "spearman", diagonal=1) %>% column_to_rownames("term"))
  z2 <- pmin(pmax(0.5 * log((1+r2) / (1-r2)), -4), 4)
  
  # browser()
  # compute fisher Z test statistic
  z <- (z1 - z2) / sqrt(1/(nrow(X1) - 3) + 1/(nrow(X2) - 3))
  z <- pmin(pmax(z, -4), 4)
  p_fisher_all <- 2*pmin(pnorm(z), pnorm(-z))
  if (length(which(is.na(p_fisher_all))) > 0) {
    stop("NAs are not expected here")
  }
  # select pairs of variables of interest
  p_fisher <- p_fisher_all
  p_fisher[,] <- NA
  p_fisher[idx1,idx2] <- p_fisher_all[idx1,idx2]
  p_fisher[idx2,idx1] <- p_fisher_all[idx2,idx1]
  diag(p_fisher) <- NA
  # multiple-hypothesis correction
  flat_fisher <- as.numeric(p_fisher)
  print(sum(!is.na(flat_fisher)))
  flat_fisher <- p.adjust(flat_fisher, method = p.adjust.method)
  p_fisher.adj <- matrix(flat_fisher, ncol = ncol(p_fisher),
                         dimnames = dimnames(p_fisher))
  
  # bootstrap testing
  bootstrapper <- function(idx) {
    i <- idx[1]
    j <- idx[2]
    if ((i > j) && (((i %in% idx1) && (j %in% idx2)) || ((j %in% idx1) && (i %in% idx2)))) {
      print(c(i,j))
      tmp <- twocorci(X1[,i], X1[,j], X2[,i], X2[,j], method = "spearman")
      return(tmp$p.value)
    } else {
      return(NA)
    }
  }
  p_bootstrap <- NULL
  p_bootstrap.adj <- NULL
  
  # browser()
  if (bootstrap) {
    if (bootstrap.method == "for") {
      p_bootstrap <- matrix(NA, nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
      for (i in 1:nrow(p_bootstrap)) {
        for (j in 1:ncol(p_bootstrap)) {
          print(c(i, j))
          tmp <- twocorci(X1[,i], X1[,j], X2[,i], X2[,j], method = "spearman")
          p_bootstrap[i, j] <- tmp$p.value
        }
      }
    } else if (bootstrap.method == "lapply") {
      idxs <- split(as.matrix(unname(expand.grid(1:nrow(z), 1:ncol(z)))),
                    seq(length(z)))
      # results in out2 don't match out...
      p_bootstrap <- matrix(as.numeric(lapply(idxs, bootstrapper)),
                            nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
    } else if (bootstrap.method == "mclapply") {
      idxs <- split(as.matrix(unname(expand.grid(1:nrow(z), 1:ncol(z)))),
                    seq(length(z)))
      # results in out2 don't match out...
      if (is.na(mc.cores)) {
        mc.cores <- detectCores() - 1
        print(mc.cores)
      }
      p_bootstrap <- matrix(as.numeric(mclapply(idxs, bootstrapper, mc.cores = mc.cores)),
                            nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
    }
    p_bootstrap[is.na(p_bootstrap) & !is.na(t(p_bootstrap))] <- t(p_bootstrap)[is.na(p_bootstrap) & !is.na(t(p_bootstrap))]
    flat_bootstrap <- as.numeric(p_bootstrap)
    print(sum(!is.na(flat_bootstrap)))
    flat_bootstrap <- p.adjust(flat_bootstrap, method = p.adjust.method)
    p_bootstrap.adj <- matrix(flat_bootstrap, ncol = ncol(p_bootstrap),
                              dimnames = dimnames(p_bootstrap))
    display_text <- ifelse(!is.na(p_bootstrap.adj) & p_bootstrap.adj < FDR, star, "")
  } else {
    display_text <- ifelse(!is.na(p_fisher.adj) & p_fisher.adj < FDR, star, "")
  }
  
  min_idx <- round((min(z) + max(abs(z))) / (2*max(abs(z))) * 99)+1
  max_idx <- round((max(z) + max(abs(z))) / (2*max(abs(z))) * 99)+1
  color <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)[min_idx:max_idx]
  
  callback = function(hc, mat) {
    # View(mat)
    dend = reorder(as.dendrogram(hc), -rowMeans(mat))
    as.hclust(dend)
  }
  
  # out <- pheatmap(z, color = color, display_numbers = display_text,
  #          clustering_callback = callback, ...)
  
  out<-pheatmap(z, color = colorRampPalette(rev(brewer.pal(n = 7, name = "PuOr")))(100), display_numbers = display_text, fontsize_number = 44, fontsize_row=34, fontsize_col=34,cluster_rows =FALSE, cluster_cols=FALSE)
  
  z[] 
  figurename = paste("C:/Users/WJ/Documents/research/A_immunology/B-hpv/figure/z", toString(fig), ".tif")
  ggsave(plot = out, height = 15, width = 15, dpi = 100,
         filename = figurename, device='tiff')
  
  list(r1 = r1, r2 = r2, z1 = z1, z2 = z2, z = z,
       p_fisher = p_fisher, p_fisher.adj = p_fisher.adj,
       p_bootstrap = p_bootstrap, p_bootstrap.adj = p_bootstrap.adj,
       pheatmap.res = out)
}

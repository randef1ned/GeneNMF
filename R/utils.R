#Calculate Jaccard Index
jaccardIndex <- function(a, b) {
  i = length(intersect(a, b))
  u = length(a) + length(b) - i
  return (i/u)
}

jaccardSimilarity <- function(gene.vectors) {
  nprogs <- length(gene.vectors)
  J <- matrix(data=0, ncol=nprogs, nrow = nprogs)
  colnames(J) <- names(gene.vectors)
  rownames(J) <- names(gene.vectors)
  for (i in 1:nprogs) {
    for (j in 1:nprogs) {
      J[i,j] <- jaccardIndex(names(gene.vectors[[i]]), names(gene.vectors[[j]]))
    }  
  }
  return(J)
}

#Calculate cosine similarity
cosineSimilarity <- function(gene.table) {
  cosine.matrix <- lsa::cosine(as.matrix(gene.table))
  return(cosine.matrix)
}

#Downsample vector, while keeping at least one element per class
downsampleMin <- function(vector, size=500) {
  if (length(vector) <= size) {
    return(vector)
  }
  vecout <- names(vector)[!duplicated(vector)] #first match
  vector <- names(vector)
  vector <- setdiff(vector, vecout)
  if (length(vecout) >= size) {
    return(vecout)
  }
  ss <- sample(vector, size = size - length(vecout))
  return(c(vecout, ss))
}

#From list of genes to complete table, with imputed zeros
geneList2table <- function(gene.vectors) {
  names <- names(gene.vectors)
  gene.vectors <- lapply(names, function(n) {
    g <- gene.vectors[[n]]
    colnames(g) <- paste(n,seq(1,ncol(g)),sep=".")
    g
  })
  gene.table <- Reduce(f=cbind, x=gene.vectors)
  return(gene.table)
}

#Calculate entropy
getEntropy <- function(tab, pseudo=0.1) {
  p <- (tab+pseudo) / sum(tab+pseudo)
  entropy <- -sum(p * log2(p), na.rm = TRUE)
  return(entropy)
}

#Find highly variable genes in a list of Seurat objects
findHVG <- function(obj.list, nfeatures=2000,
                    min.exp=0.01, max.exp=3.0, hvg.blocklist=NULL) {
  obj.list <- lapply(obj.list, function(x){
    
    ncalc <- min(5*nfeatures, nrow(x))
    x <- findVariableFeatures_wfilters(x, nfeatures=ncalc, min.exp=min.exp, max.exp=max.exp,
                                       genesBlockList=hvg.blocklist)
    x
  })
  hvg <- Seurat::SelectIntegrationFeatures(obj.list,nfeatures = nfeatures,
                                           verbose = FALSE)
  return(hvg)
}

#Calculate metrics for meta-programs
get_metaprogram_consensus <- function(nmf.wgt,
                                      nMP=10,
                                      min.confidence=0.5,
                                      weight.explained=0.5,
                                      max.genes=200,
                                      cl_members=NULL) {
  
  #calculate genes that explain 80% of weight in individual samples
  #this is used to calculate gene confidence 
  nmf.genes.single <- getNMFgenes(nmf.res=nmf.wgt,
                                     specificity.weight=NULL,
                                     weight.explained=0.9,
                                     max.genes=1000) 
  
  markers.consensus <- lapply(seq(1, nMP), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    gene.table <- geneList2table(nmf.wgt)[,which.samples]
    
    genes.avg <- apply(as.matrix(gene.table), 1, function(x){
      mean <- mean(x)
      if (length(x) >=3) { #remove outliers (SD only with 3 or more points)
        sd <- sd(x)
        x.out <- x[x>mean-3*sd & x<mean+3*sd]  
      } else {
        x.out <- mean
      }
      
      mean(x.out)
    })
    #first criterion: explain x% of total weight
    genes.avg <- sort(genes.avg, decreasing = T)
    genes.pass <- weightCumul(genes.avg, weight.explained=weight.explained)
    
    #second criterion: consistently detect a given gene across runs
    this <- nmf.genes.single[which.samples]
    genes.only <- lapply(this, names)
    genes.sum <- sort(table(unlist(genes.only)), decreasing=T)
    genes.confidence <- genes.sum/length(this)
    genes.confidence <- genes.confidence[genes.confidence > min.confidence]
    
    genes.pass <- genes.pass[names(genes.pass) %in% names(genes.confidence)]
    
    head(genes.pass, min(length(genes.pass), max.genes))
  })
  
  names(markers.consensus) <- paste0("MetaProgram",seq(1,nMP))
  return(markers.consensus)
}

#Calculate metrics for meta-programs
get_metaprogram_metrics <- function(J=NULL, Jdist=NULL,
                                   markers.consensus=NULL,
                                   cl_members=NULL) {
  nMP <- length(markers.consensus)
  all.samples <- unique(gsub("\\.k\\d+\\.\\d+","",colnames(J)))
  sample.coverage <- lapply(seq(1, nMP), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    ss <- gsub("\\.k\\d+\\.\\d+","",which.samples)
    ss <- factor(ss, levels=all.samples)
    ss.tab <- table(ss)
    #Percent samples represented
    sum(ss.tab>0)/length(ss.tab)
  })
  names(sample.coverage) <- paste0("MetaProgram",seq(1,nMP))
  
  #calculate MP silhouettes
  sil <- cluster::silhouette(cl_members, dist=Jdist)
  sil.widths <- summary(sil)$clus.avg.widths
  names(sil.widths) <- paste0("MetaProgram",seq(1,nMP))
  
  #calculate MP internal average similarity
  clusterSim <- rep(NA,nMP)
  for(i in seq_len(nMP)){
    selectMP <- which(cl_members==i)
    if (length(selectMP) > 1) { #needs at least two values
      selectJ <- J[selectMP,selectMP]
      value <- round(mean(selectJ[upper.tri(selectJ)]),3)
    } else {
      value <- 0
    }
    clusterSim[i] <- value
  }
  #number of genes in each meta-program
  metaprograms.length <- unlist(lapply(markers.consensus,length))
  
  #number of programs in meta-program
  metaprograms.size <- as.character(table(cl_members))
  
  metaprograms.metrics <- data.frame(
    sampleCoverage=unlist(sample.coverage),
    silhouette=sil.widths,
    meanSimilarity=clusterSim,
    numberGenes=metaprograms.length,
    numberPrograms=metaprograms.size)
  
  rownames(metaprograms.metrics) <- paste0("MetaProgram",seq(1,nMP))
  
  return(metaprograms.metrics)
}

#In which samples was a MP detected?
get_metaprogram_composition <- function(J=NULL,
                                    markers.consensus=NULL,
                                    cl_members=NULL) {
  nMP <- length(markers.consensus)
  all.samples <- unique(gsub("\\.k\\d+\\.\\d+","",colnames(J)))
  sample.comp <- lapply(seq(1, nMP), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    ss <- gsub("\\.k\\d+\\.\\d+","",which.samples)
    ss <- factor(ss, levels=all.samples)
    table(ss)
  })
  names(sample.comp) <- paste0("MetaProgram",seq(1,nMP))
  composition <- do.call(rbind, sample.comp)
  return(composition)
}

#Split positive and negative components of PCA, and reorder by variance
nonNegativePCA <- function(pca, k) {
  
  pcaP <- pca$rotation
  pcaN <- pca$rotation
  colnames(pcaP) <- paste0(colnames(pcaP),"p")
  colnames(pcaN) <- paste0(colnames(pcaN),"n")
  
  pcaP[pcaP<0] <- 0
  pcaN[pcaN>0] <- 0
  pcaN <- abs(pcaN)
  
  sumP <- apply(pcaP, 2, sum)
  sumN <- apply(pcaN, 2, sum)
  
  wP <- pca$sdev * sumP/(sumP+sumN)
  wN <- pca$sdev * sumN/(sumP+sumN)
  wSort <- sort(c(wP,wN), decreasing = T)
  
  #Collate, re-rank components, and rescale coefficients
  pca_abs <- cbind(pcaP, pcaN)
  pca_abs <- pca_abs[,names(wSort)[1:k]]
#  pca_abs <- apply(pca_abs, 2, function(x){x/sum(x)})
  return(pca_abs)
} 

#Weighting factor matrix by feature specificity
weightedLoadings <- function(nmf.res,
                        specificity.weight=5) {
  
  nmf.wgt <- lapply(nmf.res, function(model) {
    wgtLoad(model$w, w=specificity.weight)
  })
  return(nmf.wgt)
}

wgtLoad <- function(matrix, w) {
  rownorm <- apply(matrix, 1, normVector)
  spec <- apply(rownorm, 2, max)
  spec.w <- spec^w
  matrix <- matrix * spec.w
  #renormalize
  apply(matrix, 2, normVector)
}

normVector <- function(vector) {
  s <- sum(vector)
  if (s>0) {
    vector <- vector/s
  }
  return(vector)
}

weightCumul <- function(vector, weight.explained=0.5) {
    x.sorted <- sort(vector, decreasing = T)
    norm.x <- normVector(x.sorted)
    cs <- cumsum(norm.x)
    norm.x[cs<weight.explained]
}

check_cpp_version <- function(model) {
  cppversion <- packageVersion("RcppML")
  if(cppversion >= "0.5.6"){
    model<-list(w=model@w, 
                d=model@d, 
                h=model@h, 
                tol=model@misc$tol,
                iter=model@misc$iter)
  }
  return(model)
}

#' @importFrom reticulate py_to_r use_condaenv import
#' @importFrom dplyr `%>%`
cuda_nmf <- function(mat, k, tol = 1e-4, maxit = 100L, L1 = c(0, 0), L2 = c(0, 0), 
                     seed = seed,
                     cuda.device = 'cuda:0', cuda.env = 'torchnmf', 
                     cuda.conda_binary = '/opt/mambaforge/bin/conda', 
                     cuda.verbose = TRUE) {
  tol <- as.integer(tol)
  maxit <- as.integer(maxit)
  k <- as.integer(k)
  seed <- as.integer(seed)
  
  use_condaenv(cuda.env, conda = cuda.conda_binary)
  np       <- import('numpy', convert = FALSE)
  torch    <- import('torch', convert = FALSE)
  torchnmf <- import('torchnmf', convert = FALSE)
  torch$manual_seed(seed)
  
  mat_mat <- as.matrix(mat) %>% 
    t() %>%
    pmax(0) %>%
    {np$array(.)} %>%
    {torch$from_numpy(.)$float()$to(cuda.device)}
  
  nmf_model <- torchnmf$nmf$NMF(mat_mat$shape, rank = k)$to(cuda.device)
  nmf_model$fit(mat_mat, max_iter = maxit, verbose = cuda.verbose, tol = tol)
  
  eps <- 1e-12
  W_t <- nmf_model$W$detach()   # likely (n_rows, k)
  H_t <- nmf_model$H$detach()   # likely (n_cols, k)
  
  if (sum(abs(L1)) + sum(abs(L2)) > 0) {
    res <- nmf_fit_pg(mat_mat, W_t, H_t, torch = torch,
                      maxit = maxit, tol = tol,
                      L1 = L1, L2 = L2,
                      inner_iters = maxit, verbose = cuda.verbose)
    # write back into the model (in-place on GPU)
    nmf_model$W$copy_(res$W)
    nmf_model$H$copy_(res$H)
  }
  
  # --- scaleH analogue (normalize "component" axis) ---
  # Since H is (n_cols, k), normalize columns (dim=0)
  d_scaleh <- H_t$sum(dim = 0L)$add(eps)             # (k,)
  H_scaleh <- H_t$div(d_scaleh$unsqueeze(0L))        # column-normalize
  W_scaleh <- W_t$mul(d_scaleh$unsqueeze(0L))        # push scaling into W (broadcast over rows)
  
  # --- scaleW analogue ---
  d_scalew <- W_scaleh$sum(dim = 0L)$add(eps)        # (k,)  <-- this is the Rcpp-style d
  W_scalew <- W_scaleh$div(d_scalew$unsqueeze(0L))   # columns sum to 1
  
  # export
  model <- list(
    w = py_to_r(W_scalew$cpu()$numpy()) %>%
      `rownames<-`(rownames(mat)) %>%
      `colnames<-`(paste0('nmf', seq_len(ncol(.)))) %>%
      as.matrix(),
    d = py_to_r(d_scalew$cpu()$numpy()) %>%
      as.integer(),
    h = py_to_r(H_scaleh$cpu()$numpy()) %>%
      `rownames<-`(paste0('Col_', seq_len(nrow(.)))) %>%
      `colnames<-`(paste0('nmf', seq_len(ncol(.)))) %>%
      as.matrix() %>%
      t(),
    tol = tol,
    iter = maxit,
    mse = 0
  )
  rm(list = c('nmf_model', 'mat_mat', 'd_scalew', 'W_scalew', 'd_scaleh', 'H_scaleh', 'W_scaleh', 'H_t', 'W_t'))
  torch$cuda$empty_cache()
  
  return(model)
}

# ---- Projected GD updates mimicking RcppML L1/L2 inside predict() ----
update_H_pg <- function(A, W, H, torch, L1 = 0, L2 = 0, iters = 50L, eps = 1e-12) {
  # A: (n_rows, n_cols)  CUDA
  # W: (n_rows, k)       CUDA
  # H: (n_cols, k)       CUDA  where Ahat = W @ H^T
  
  with(torch$no_grad(), {
    G <- W$transpose(0L, 1L)$matmul(W)      # (k,k)   = W^T W
    C <- A$transpose(0L, 1L)$matmul(W)      # (n_cols,k) = A^T W
    
    # Lipschitz constant for gradient: ||G||_2 + L2
    # k is small, so SVD is fine.
    svals <- torch$linalg$svdvals(G)
    Llip  <- svals[[0L]]$add(L2)$add(eps)
    step  <- torch$reciprocal(Llip)
    
    for (t in seq_len(iters)) {
      grad <- H$matmul(G)$sub(C)            # H G - C
      if (L2 != 0) grad <- grad$add(H$mul(L2))
      if (L1 != 0) grad <- grad$add(L1)     # + L1 (since H>=0, ||H||1 = sum(H))
      
      H$sub_(grad$mul(step))
      H$clamp_(min = 0)
    }
    H
  })
}

update_W_pg <- function(A, W, H, torch, L1 = 0, L2 = 0, iters = 50L, eps = 1e-12) {
  # Update W given H, same convention Ahat = W @ H^T
  
  with(torch$no_grad(), {
    G <- H$transpose(0L, 1L)$matmul(H)      # (k,k) = H^T H
    C <- A$matmul(H)                        # (n_rows,k) = A H
    
    svals <- torch$linalg$svdvals(G)
    Llip  <- svals[[0L]]$add(L2)$add(eps)
    step  <- torch$reciprocal(Llip)
    
    for (t in seq_len(iters)) {
      grad <- W$matmul(G)$sub(C)            # W G - C
      if (L2 != 0) grad <- grad$add(W$mul(L2))
      if (L1 != 0) grad <- grad$add(L1)
      
      W$sub_(grad$mul(step))
      W$clamp_(min = 0)
    }
    W
  })
}

# ---- One ALS loop with L1/L2 like RcppML (L1[0]=W, L1[1]=H) ----
nmf_fit_pg <- function(A, W, H, torch,
                       maxit = 100L, tol = 1e-4,
                       L1 = c(0, 0), L2 = c(0, 0),
                       inner_iters = 50L, verbose = TRUE) {
  # Correlation-like stopping (Rcpp uses cor(w, w_prev))
  # Here we use cosine similarity on flattened W as a proxy.
  cos_sim <- function(x, y, eps = 1e-12) {
    num <- (x * y)$sum()
    den <- torch$sqrt((x * x)$sum() * (y * y)$sum())$add(eps)
    num$div(den)
  }
  
  with(torch$no_grad(), {
    for (it in seq_len(maxit)) {
      W_prev <- W$clone()
      
      # update H then W (regularized NNLS via projected GD)
      H <- update_H_pg(A, W, H, L1 = L1[[2]], L2 = L2[[2]], iters = inner_iters)
      W <- update_W_pg(A, W, H, L1 = L1[[1]], L2 = L2[[1]], iters = inner_iters)
      
      # stopping criterion
      sim <- cos_sim(W$reshape(-1L), W_prev$reshape(-1L))
      tol_now <- (1 - sim)$item()
      
      if (verbose) cat(sprintf("%4d | tol %.3e\n", it, tol_now))
      if (tol_now < tol) break
    }
    list(W = W, H = H)
  })
}

utils::globalVariables(".")

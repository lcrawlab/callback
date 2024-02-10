

#' @title Returns the genes selected by the knockoff filter
#'
#' @description Given two Seurat objects, returns the  the genes selected by the knockoff filter and their W statistics.
#' @details 
#'
#' @param seurat_obj A Seurat object 
#' @param all_markers The output from Seurat::FindAllMarkers.
#' @param cluster1 The Idents of the cluster of interest in seurat_obj1
#' @param cluster2 The Idents of the cluster of interest in seurat_obj2
#' @param q The desired rate to control the FDR at
#' @param return_all Determines if the returned object will contain all genes or just the selected genes.
#' @param threshold One of "fdr", "kfwer", or "heuristic".
#' @param num_cores The number of cores for computing marker genes in parallel.
#' @returns The 
#' @name compute_knockoff_filter_precomputed_markers
compute_knockoff_filter_precomputed_markers <- function(seurat_obj, all_markers, cluster1, cluster2, q, return_all=FALSE, threshold="fdr", num_cores=1) {
  #library(future)
  options(future.globals.maxSize = 8000 * 1024^2) # todo note what this is for, figure this out as a parameter or programmatically
  future::plan("multicore", workers = as.numeric(num_cores)) # todo log number of cores being used
  markers <- subset(all_markers, cluster %in% c(cluster1, cluster2))
  

  # FindMarkers orders by p-value, so we can't rely on position to know which genes are which
  knockoff.indices <- grepl("^knockoff", rownames(markers))
  original.indices <- !knockoff.indices
  
  # subset the markers data.frame into originals and knockoffs
  knockoff.markers <- markers[knockoff.indices, ]
  original.markers <- markers[original.indices, ]
  
  all.genes <- rownames(seurat_obj)
  
  # get indices of knockoffs and originals from seurat_obj, should be [FALSE, ..., FALSE, TRUE, ..., TRUE]
  knockoff.indices.sorted <- grepl("^knockoff", all.genes)
  original.indices.sorted <- !knockoff.indices.sorted
  
  knockoff_names_sorted <- all.genes[knockoff.indices.sorted]
  original_names_sorted <- all.genes[original.indices.sorted]

  # sort markers data.frames by their original orderings
  knockoff.markers.sorted <- knockoff.markers[knockoff_names_sorted, ]
  original.markers.sorted <- original.markers[original_names_sorted, ]
  
  original_p_values <- original.markers.sorted$p_val
  knockoff_p_values <- knockoff.markers.sorted$p_val
  
  log_original_p_values <- -log10(original_p_values)
  log_knockoff_p_values <- -log10(knockoff_p_values)
  
  W <- log_original_p_values - log_knockoff_p_values

  if (threshold == "fdr") {
    thres = knockoff::knockoff.threshold(W, fdr=q, offset=1)
  }

  if (threshold == "kfwer") {
    k <- 10
    alpha <- q
    thres <- knockoff.kfwer.threshold(W, k, alpha)
  }

  if (threshold == "heuristic") {
    k <- 10
    alpha <- q
    thres <- knockoff.heuristic.threshold(W, q, k, alpha)
  }

  message(paste("threshold:", thres))

  if (return_all) {
    all_features <- as.data.frame(list("gene" = original_names_sorted, "W" = W))

    ret <-  list("all_features"=all_features, "threshold"=thres)

    return(ret)
  }
  selected_indices = which(W >= thres) # todo check if this should be > (case where threshold is Inf, but there are still some Inf -log p)
  #selected_indices = which(W > thres) # todo check if this should be > (case where threshold is Inf, but there are still some Inf -log p)
  

  message("Num selected indices:")
  message(length(selected_indices))
  selected_genes <- original_names_sorted[selected_indices]
  selected_Ws <- W[selected_indices]
  
  selected_features <- as.data.frame(list("selected_gene" = selected_genes, "W" = selected_Ws))
  
  selected_features <- selected_features[order(selected_features$W, decreasing = TRUE),]

  ret <-  list("selected_features"=selected_features, "threshold"=thres)
  
  return(ret)
}






#' @title Runs a typical Seurat workflow on a Seurat object (up to dimensionality reduction and clustering).
#'
#' @description Given a Seurat object, returns a new Seurat that has been normalized, had variable features identified,
#' scaled, had principal components computed, had clusters identified, and had tSNE and UMAP embeddings determined.
#'
#' @details 
#'
#' @param seurat_obj The Seurat object that will be analyzed.
#' @param resolution_start The starting resolution to be used for the clustering algorithm (Louvain and Leiden algorithms).
#' @param num_clusters_start The starting number of clusters to be used for the clustering algorithm (K-means and Hierarchical clustering algorithms).
#' @param reduction_percentage The amount that the starting parameter will be reduced by after each iteration (between 0 and 1).
#' @param dims The dimensions to use as input features (i.e. 1:10).
#' @param algorithm The clustering algorithm to be used.
#' @param cores The number of cores to compute marker genes in parallel.
#' @returns Returns a Seurat object where the idents have been updated with the clusters determined via the KCC algorithm. 
#' Latest clustering results will be stored in the object metadata under 'kcc_clusters'. 
#' Note that 'kcc_clusters' will be overwritten every time FindClustersKC is run.
#' @name FindClustersKCCFast
#' @export
FindClustersKCCFast <- function(seurat_obj,
                            resolution_start=0.8,
                            reduction_percentage=0.2,
                            num_clusters_start=20,
                            dims=1:10,
                            algorithm="louvain", # todo implement all algorithms in one function
                            cores=1) {

  # todo check function arguments for validity
  

  knockoff_seurat_obj <- get_seurat_obj_with_knockoffs(seurat_obj)

  num_variable_features <- 2*length(Seurat::VariableFeatures(seurat_obj))


  # Pre-process data

  #library(future)
  options(future.globals.maxSize = 8000 * 1024^2)
  future::plan("multicore", workers = as.numeric(cores)) # todo log number of cores being used
  #options(future.globals.maxSize = 8000 * 1024^2)
  #plan("multicore", workers = as.numeric(cores)) # todo log number of cores being used

  knockoff_seurat_obj <- Seurat::NormalizeData(knockoff_seurat_obj)
   
  knockoff_seurat_obj <- Seurat::FindVariableFeatures(knockoff_seurat_obj, selection.method = "vst", nfeatures = num_variable_features)
    
  all.genes <- rownames(knockoff_seurat_obj)
    
  knockoff_seurat_obj <- Seurat::ScaleData(knockoff_seurat_obj)
    
  knockoff_seurat_obj <- Seurat::RunPCA(knockoff_seurat_obj, features = Seurat::VariableFeatures(object = knockoff_seurat_obj))
    
  knockoff_seurat_obj <- Seurat::FindNeighbors(knockoff_seurat_obj, dims = dims) # todo check if i should use all dims for knockoffs

  resolution_param <- resolution_start

  while(TRUE) {




    message("####################################################################")
    message("Resolution param:")
    message(resolution_param)
    message("####################################################################")


    message("Finding clusters")

    if (algorithm == "louvain") {
      message("Louvain")
      knockoff_seurat_obj <- Seurat::FindClusters(knockoff_seurat_obj, resolution = resolution_param)
    }

    if (algorithm == "leiden") {
      message("Leiden")
      #plan("sequential") # todo log number of cores being used # this is a weird one because leiden has a forked job hanging
      knockoff_seurat_obj <- Seurat::FindClusters(knockoff_seurat_obj, resolution = resolution_param, algorithm=4, method = "igraph")
    }
    message("Found clusters")

    # Reduce resolution for next iteration of the loop
    resolution_param <- (1 -reduction_percentage) * resolution_param



    k <- length(levels(Seurat::Idents(knockoff_seurat_obj)))
    #knock_idents <- 0:(k-1)

    message("Num clusters:")
    message(k)

    knock_idents <- levels(Seurat::Idents(knockoff_seurat_obj))


    num_selected_matrix <- matrix(nrow=k, ncol=k)

    found_no_sign_diff <- FALSE

    all_markers <- Seurat::FindAllMarkers(knockoff_seurat_obj,
                                          features = Seurat::VariableFeatures(knockoff_seurat_obj),
                                          logfc.threshold=0.0,
                                          min.pct=0.0,
                                          return.thresh = 1)


    
    m <- 0
    for (i in 1:length(knock_idents)) {
      for (j in 1:length(knock_idents)) {
        if (j >= i) {
          next
        }
        
        m <- m + 1
        
        message("Pair:")
        message(paste(i,j))
        
        message("Knockoff Pair:")
        message(paste(knock_idents[i], knock_idents[j]))
        
        markers.selected <- compute_knockoff_filter_precomputed_markers(knockoff_seurat_obj, knock_idents[i], knock_idents[j], 0.05, num_cores=cores)
        
        num.selected <- nrow(markers.selected$selected_features)
        
        if (num.selected == 0) {
          found_no_sign_diff <- TRUE
          break
        }

        num_selected_matrix[i,j] <- num.selected
        num_selected_matrix[j,i] <- num.selected

      }
      if (found_no_sign_diff) {
        message("Found clusters with no significant differences. Progressing to next clustering iteration.")
        break
      }
    }
    
    if (found_no_sign_diff) {
      next
    }
    break
  }


  seurat_obj@meta.data$kcc_clusters <- Seurat::Idents(knockoff_seurat_obj)

  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data$kcc_clusters

  return(seurat_obj)
}


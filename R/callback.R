

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
#' @param assay The assay to generate knockoffs from.
#' @param cores The number of cores to compute marker genes in parallel.
#' @returns Returns a Seurat object where the idents have been updated with the clusters determined via the callback algorithm. 
#' @param verbose Whether or not to show all logging.
#' Latest clustering results will be stored in the object metadata under 'callback_clusters'. 
#' Note that 'callback_clusters' will be overwritten every time FindClustersKC is run.
#' @name FindClustersCallback
#' @export
FindClustersCallback <- function(seurat_obj,
                            resolution_start=0.8,
                            reduction_percentage=0.2,
                            num_clusters_start=20,
                            dims=1:10,
                            algorithm="louvain", # todo implement all algorithms in one function
                            assay="RNA",
                            cores=1,
                            verbose=TRUE) {

  # todo check function arguments for validity
  

  knockoff_seurat_obj <- get_seurat_obj_with_knockoffs(seurat_obj, assay=assay)

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




    if (verbose) {
      message("####################################################################")
      message("Resolution param:")
      message(resolution_param)
      message("####################################################################")


      message("Finding clusters")
    }

    if (algorithm == "louvain") {
      if (verbose) {
        message("Louvain")
      }
      knockoff_seurat_obj <- Seurat::FindClusters(knockoff_seurat_obj, resolution = resolution_param)
    }

    if (algorithm == "leiden") {
      if (verbose) {
        message("Leiden")
      }
      #plan("sequential") # todo log number of cores being used # this is a weird one because leiden has a forked job hanging
      knockoff_seurat_obj <- Seurat::FindClusters(knockoff_seurat_obj, resolution = resolution_param, algorithm=4, method = "igraph")
    }
    message("Found clusters")

    # Reduce resolution for next iteration of the loop
    resolution_param <- (1 -reduction_percentage) * resolution_param



    k <- length(levels(Seurat::Idents(knockoff_seurat_obj)))
    #knock_idents <- 0:(k-1)

    if (verbose) {
      message("Num clusters:")
      message(k)
    }

    knock_idents <- levels(Seurat::Idents(knockoff_seurat_obj))


    num_selected_matrix <- matrix(nrow=k, ncol=k)

    found_no_sign_diff <- FALSE
    
    #pb <- progress::progress_bar$new(clear = FALSE)
    cli::cli_progress_bar("Cleaning data", total = 100)

    m <- 0
    for (i in 1:length(knock_idents)) {
      for (j in 1:length(knock_idents)) {
        if (j >= i) {
          next
        }
        
        m <- m + 1
        
        #pb$tick(100 * 1 / (length(knock_idents) * (length(knock_idents) - 1) / 2))
        cli::cli_progress_update(100 * 1 / (length(knock_idents) * (length(knock_idents) - 1) / 2))

        if (verbose) {
          #message("Pair:")
          #message(paste(i,j))
          
          #message("Knockoff Pair:")
          #message(paste(knock_idents[i], knock_idents[j]))
        }
        
        markers.selected <- compute_knockoff_filter(knockoff_seurat_obj, knock_idents[i], knock_idents[j], 0.05, num_cores=cores)
        
        num.selected <- nrow(markers.selected$selected_features)
        
        if (num.selected == 0) {
          found_no_sign_diff <- TRUE
          break
        }

        num_selected_matrix[i,j] <- num.selected
        num_selected_matrix[j,i] <- num.selected

      }
      if (found_no_sign_diff) {
        if (verbose) {
          cli::cli_progress_done()
          message("Found clusters with no significant differences. Progressing to next clustering iteration.")
        }
        break
      }
    }
    
    if (found_no_sign_diff) {
      next
    }
    break
  }


  seurat_obj@meta.data$callback_clusters <- Seurat::Idents(knockoff_seurat_obj)

  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data$callback_clusters

  return(seurat_obj)
}



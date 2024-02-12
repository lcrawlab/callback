#' @title Runs a typical Seurat workflow on a Seurat object (up to
#' dimensionality reduction and clustering).
#'
#' @description Given a Seurat object, returns a new Seurat that has been
#' normalized, had variable features identified, scaled, had principal
#' components computed, hadclusters identified, and had tSNE and UMAP
#' embeddings determined.
#'
#' @details
#'
#' @param seurat_obj The Seurat object that will be analyzed.
#' @param resolution_start The starting resolution to be used for the
#' clustering algorithm (Louvain and Leiden algorithms).
#' @param num_clusters_start The starting number of clusters to be used for the
#' clustering algorithm (K-means and Hierarchical clustering algorithms).
#' @param reduction_percentage The amount that the starting parameter will be
#' reduced by after each iteration (between 0 and 1).
#' @param dims The dimensions to use as input features (i.e. 1:10).
#' @param algorithm The clustering algorithm to be used.
#' @param assay The assay to generate knockoffs from.
#' @param cores The number of cores to compute marker genes in parallel.
#' @returns Returns a Seurat object where the idents have been updated with the
#' clusters determined via the callback algorithm.
#' @param verbose Whether or not to show all logging.
#' Latest clustering results will be stored in the object metadata under
#' callback_clusters'. Note that 'callback_clusters' will be overwritten ever
#' time FindClustersKC is run.
#' @name FindClustersCallback
#' @export
FindClustersCallback <- function(seurat_obj,
                                 resolution_start = 0.8,
                                 reduction_percentage = 0.2,
                                 num_clusters_start = 20,
                                 dims = 1:10,
                                 algorithm = "louvain", # todo implement all algos
                                 assay = "RNA",
                                 cores = 1,
                                 verbose = TRUE) {

  # todo check function arguments for validity

  knockoff_seurat_obj <- get_seurat_obj_with_knockoffs(seurat_obj, assay = assay)

  num_variable_features <- 2 * length(Seurat::VariableFeatures(seurat_obj))


  # Pre-process data

  #library(future)
  options(future.globals.maxSize = 8000 * 1024^2)
  # todo log number of cores being used
  future::plan("multicore", workers = as.numeric(cores))
  #options(future.globals.maxSize = 8000 * 1024^2)
  # todo log number of cores being used
  #plan("multicore", workers = as.numeric(cores))

  knockoff_seurat_obj <- Seurat::NormalizeData(knockoff_seurat_obj,
                                               verbose = FALSE)
   
  knockoff_seurat_obj <- Seurat::FindVariableFeatures(knockoff_seurat_obj,
                                                      selection.method = "vst",
                                                      nfeatures = num_variable_features,
                                                      verbose = FALSE)
    
  all.genes <- rownames(knockoff_seurat_obj)

  knockoff_seurat_obj <- Seurat::ScaleData(knockoff_seurat_obj, verbose = FALSE)
  knockoff_seurat_obj <- Seurat::RunPCA(knockoff_seurat_obj,
                                        features = Seurat::VariableFeatures(object = knockoff_seurat_obj),
                                        verbose = FALSE)
  # todo check if i should use all dims for knockoffs
  knockoff_seurat_obj <- Seurat::FindNeighbors(knockoff_seurat_obj,
                                               dims = dims,
                                               verbose = FALSE)

  resolution_param <- resolution_start

  while (TRUE) {
    if (verbose) {
      message("####################################################################")
      message(paste("Finding clusters with", stringr::str_to_title(algorithm), "algorithm"))
      message(paste("Resolution param:", resolution_param))
    }

    if (algorithm == "louvain") {
      knockoff_seurat_obj <- Seurat::FindClusters(knockoff_seurat_obj,
                                                  resolution = resolution_param,
                                                  verbose = FALSE)
    }

    if (algorithm == "leiden") {
      #plan("sequential") # todo log number of cores being used # this is a weird one because leiden has a forked job hanging
      knockoff_seurat_obj <- Seurat::FindClusters(knockoff_seurat_obj,
                                                  resolution = resolution_param,
                                                  algorithm = 4,
                                                  method = "igraph",
                                                  verbose = FALSE)
    }

    # Reduce resolution for next iteration of the loop
    resolution_param <- (1 - reduction_percentage) * resolution_param

    k <- length(levels(Seurat::Idents(knockoff_seurat_obj)))
    #knock_idents <- 0:(k-1)

    if (verbose) {
      message("Num clusters:")
      message(k)
    }

    knock_idents <- levels(Seurat::Idents(knockoff_seurat_obj))

    num_selected_matrix <- matrix(nrow = k, ncol = k)

    found_no_sign_diff <- FALSE

    num_clusters <- length(knock_idents)


    if (verbose) {
      progress_bar_length <- num_clusters * (num_clusters - 1) / 2
      cli::cli_progress_bar("Processing cluster pairs:",
                            total = progress_bar_length,
                            clear = FALSE)
    }

    m <- 0
    for (i in 1:num_clusters) {
      for (j in 1:num_clusters) {
        if (j >= i) {
          next
        }

        m <- m + 1

        if (verbose) {
          cli::cli_progress_update()
        }

        markers_selected <- compute_knockoff_filter(knockoff_seurat_obj,
                                                    knock_idents[i],
                                                    knock_idents[j],
                                                    0.05, # todo add var names?
                                                    num_cores = cores)

        num_selected <- nrow(markers_selected$selected_features)

        if (num_selected == 0) {
          found_no_sign_diff <- TRUE
          break
        }

        num_selected_matrix[i, j] <- num_selected
        num_selected_matrix[j, i] <- num_selected

      }
      if (found_no_sign_diff) {
        if (verbose) {
          cli::cli_progress_done()
          message("Found clusters with no significant differences.")
          message("Progressing to next clustering iteration.")
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






#' @title Returns a Seurat object that contains additional (fake) RNA
#' expression counts in the form of knockoffs.
#'
#' @description Given a Seurat object, returns a new Seurat object whose RNA
#' expression counts includes the
#' variable features from the original object and an equal number of knockoff
#' features.
#'
#' @details
#'
#' @param seurat_obj A Seurat object containing RNA expression counts.
#' @param assay The assay to generate knockoffs from.
#' @returns A Seurat object that contains the original variable features and an
#' equal number of knockoff features.
#' @name get_seurat_obj_with_knockoffs
get_seurat_obj_with_knockoffs <- function(seurat_obj, assay = "RNA") {
  var_features <- Seurat::VariableFeatures(seurat_obj)
  #seurat_obj_data <- as.data.frame(t(as.matrix(seurat_obj@assays$RNA@counts)))

  #seurat_obj_data <- seurat_obj_data[var_features]

  message("Pulling data from Seurat object")
  #seurat_obj_data <- as.data.frame(t(as.matrix(seurat_obj@assays$RNA@counts[Seurat::VariableFeatures(seurat_obj),])))
  seurat_obj_data <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts")[Seurat::VariableFeatures(seurat_obj), ])))

  message("Computing MLE for zero-inflated poisson")
  ml_estimates <- lapply(seurat_obj_data, estimate_zi_poisson)

  message("Computing knockoffs")
  ko <- as.data.frame(lapply(ml_estimates,
                             function(x) {
                                           rzipoisson(nrow(seurat_obj_data),
                                           x$lambda.hat,
                                           x$pi.hat)
                             }))


  num_variable_features <- length(var_features)
  colnames(ko) <- paste0(rep("knockoff", num_variable_features), 1:num_variable_features)
  combined_data <- cbind(seurat_obj_data, ko)

  new_project_name <- paste0(seurat_obj@project.name, "_with_knockoffs")
  new_seurat_obj <- Seurat::CreateSeuratObject(counts = t(combined_data), project = new_project_name)

  return(new_seurat_obj)
}

#' Functions for Adjacency Matrix
#' From https://github.com/lfranzen/scwat-st/blob/e3352598412b7a8a7feff59e5cd0b999548b2179/scripts/visium_nbs_analysis.between-clusters.R

library("GGally")
library("ggnet")
library("ggplot2")
library("ggraph")
library("igraph")
library("magrittr")
library("network")
library("pheatmap")
library("scico")
library("sna")
library("STutility")

#' Create Adjacency Matrix
#' 
#' @param nbs.df Data frame with output from STUtility::GetSpatNet() data for all clusters of interest. Rows should correspond to spot ID and columns should include cluster IDs as well as colums starting with "nbs_".
#' @param column.clusters.id Column name in nbs.df corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param cluster.include Vector of cluster IDs to include in your analysis.
#' @return Adjacency matrix with number of neighbours present between each cluster pair
#' @export
CreateAdjMatrix <- function(
    nbs.df,
    column.clusters.id = "seurat_clusters",
    cluster.include
) {
  nbs_adjmat <- matrix(
    0L, nrow = length(cluster.include), ncol = length(cluster.include)
  )
  for (i in seq_along(cluster.include)) {
    c <- cluster.include[i]
    c_nbs_df <- nbs.df[nbs.df[,column.clusters.id] == c, ]
    for (j in seq_along(cluster.include)) {
      nbs <- paste0("nbs_", cluster.include[j])
      if (j == i) {
        n_nbs <- sum(!is.na(c_nbs_df[, nbs] == c))
      } else {
        n_nbs <- sum(!is.na(c_nbs_df[, nbs] == nbs))
      }
      nbs_adjmat[i,j] <- n_nbs
      if (nbs_adjmat[j, i] > 0) {
        nbs_adjmat[i, j] <- nbs_adjmat[j, i] <- max(c(nbs_adjmat[i, j], nbs_adjmat[j, i]))
      }
    }
    nbs_adjmat[i, i] <- sum(nbs_adjmat[i, ][-i])
  }
  return(nbs_adjmat)
}

### Adjacency Matrix per Subject

#' Create Adjacency Matrix per Subject
#' 
#' @param nbs.df Data frame with output from STUtility::GetSpatNet() data for all clusters of interest. Rows should correspond to spot ID and columns should include cluster IDs as well as colums starting with "nbs_".
#' @param se.metadata Metadata dataframe from seurat object containing Sample IDs for each spot.
#' @param column.clusters.id Column name in nbs.df corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param cluster.include Vector of cluster IDs to include in your analysis.
#' @param column.subject.id Column name in se.metadata corresponding to Sample ID of the spots.
#' @return List of adjacency matrices with number of neighbours present between each cluster pair. Each matrix in the list corresponds to data from one sample.
#' @export
CreateAdjMatrixPerSubject <- function(
    nbs.df,
    se.metadata,
    column.clusters.id = "seurat_clusters",
    cluster.include,
    column.subject.id
){
  nbs_adjmat_list <- list()
  subjects_include <- unique(as.character(se.metadata[, column.subject.id]))
  
  for (subject_id in subjects_include) {
    rows_subject <- rownames(
      se.metadata[se.metadata[,column.subject.id] %in% subject_id, ]
    )
    nbs.df_subject <- nbs.df[rows_subject, ]
    
    nbs_adjmat <- CreateAdjMatrix(
      nbs.df = nbs.df_subject,
      cluster.include = cluster.include,
      column.clusters.id = column.clusters.id
    )
    
    nbs_adjmat_list[[subject_id]] <- nbs_adjmat
  }
  return(nbs_adjmat_list)
}

### Randomise Cluster IDs within Subject data

#' Randomise Cluster IDs within Subject data
#' 
#' @param se.object Seurat (STUtility) object containing cluster and sample identities for each spot in the metadata.
#' @param column.clusters.id Column name in metadata corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param column.subject.id Column name in metadata corresponding to Sample ID of the spots.
#' @return New Seurat object with shuffled cluster identities per sample
#' @export
RandomiseClusteringIDs <- function(
  se.object,
  column.cluster.id,
  column.sample.id,
  random.seed = NA
) {
  if (!is.na(random.seed)) {
    message(paste("Setting random seed to", random.seed))
    set.seed(random.seed)
  }
  #' Shuffle cluster ids for each sample
  se_metadata <- se.object@meta.data[, c(column.cluster.id, column.sample.id)]
  se_metadata$clusters_original <- se_metadata[, column.cluster.id]
  se_metadata$sample_id <- se_metadata[, column.sample.id]

  se_metadata_perm <- se_metadata %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(
      clusters_perm = clusters_original[sample(dplyr::row_number())]
    )
  #' Add shuffled clusters to se object metadata
  se.object <- AddMetaData(
    se.object,
    as.character(se_metadata_perm$clusters_perm),
    col.name = "clusters_perm"
  )
  return(se.object)
}

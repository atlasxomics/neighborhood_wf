library("dplyr")
library("GGally")
library("ggnet")
library("ggplot2")
library("ggraph")
library("igraph")
library("Matrix")
library("magrittr")
library("network")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")
library("pheatmap")
library("purrr")
library("scico")
library("Seurat")
library("sna")
library("STutility")

source("wf/colors.R")

dir.create("/root/neighborhood", showWarnings = FALSE)
setwd("/root/neighborhood")

find_func <- function(tempdir, pattern) {

  list.files(
    path = tempdir, # replace with the directory you want
    pattern = pattern, # has "test", ".csv", and then nothing else ($)
    full.names = TRUE, # include the directory in the result
    recursive = TRUE
  )
}

args <- commandArgs(trailingOnly = TRUE)

project_name <- args[1]
genome <- args[2]

runs <- strsplit(args[3:length(args)], ",")
runs

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

all <-  list()
for (i in seq_along(runs)) {
  all[[i]] <- readRDS(runs[[i]][2])
  all[[i]] <- RenameCells(
    all[[i]],
    new.names = paste0(
      unique(all[[i]]@meta.data$Sample),
      "#",
      colnames(all[[i]]), "-1"
    )
  )
}
all

saveRDS(all, "/root/neighborhood/all.rds")

#===========================#===========================

main_func <- function(seurat_lst) {
  find_samples_name <- function(seurat_lst) {
    sapply(
      seq_along(seurat_lst),
      function(i) unique(seurat_lst[[i]]@meta.data$Sample))
  }

  samples <- find_samples_name(seurat_lst)

  D00_fun <- function(seurat_lst) {
    toRemove <- lapply(
      seurat_lst,
      function(x) {
        names(which(colSums(is.na(x@assays[[1]]@counts)) > 0))
      }
    )
    mapply(function(x, y) x[, !colnames(x) %in% y], seurat_lst, toRemove)
  }

  D00 <- D00_fun(seurat_lst)
  Spatial_D00_fun <- function(D00) {

    Spatial_D00 <- lapply(
      D00,
      function(x) as.data.frame(x@images[[1]]@coordinates[, c(5, 4)])
    )
    Spatial_D00 <- lapply(
      Spatial_D00,
      function(x) {
        colnames(x) <- paste0("Spatial_", 1:2)
        x
      }
    )
    lapply(
      Spatial_D00,
      function(x) {
        x$Spatial_2 <- -(x$Spatial_2)
        x
      }
    )
  }
  Spatial_D00 <- Spatial_D00_fun(D00)

  Spatial_D00_all_fun <- function(Spatial_D00) {

    tmp <- lapply(
      seq_along(Spatial_D00),
      function(i) {
        bind_rows(Spatial_D00[-i])
      }
    )

    tmp <- lapply(
      tmp,
      function(x) {
        x$Spatial_1 <- 0
        x
      }
    )
    tmp <- lapply(
      tmp,
      function(x) {
        x$Spatial_2 <- 0
        x
      }
    )
    tmp <- lapply(
      seq_along(Spatial_D00),
      function(i) {
        as.matrix(rbind(Spatial_D00[[i]], tmp[[i]]))
      }
    )
  }

  Spatial_D00_all <- Spatial_D00_all_fun(Spatial_D00)

  temp_fun <- function(D00) {
    temp <- lapply(D00, function(x) as.data.frame(x@assays[[1]]@counts))
    temp <- lapply(
      temp,
      function(x) {
        x$region <- rownames(x)
        x
      }
    )
    lapply(
      temp,
      function(x) {
        rownames(x) <- NULL
        x
      }
    )
  }

  temp <- temp_fun(D00)

  # merge seurat objects
  combined_mat <- reduce(temp, full_join, by = "region")

  rownames(combined_mat) <- combined_mat$region
  combined_mat$region <- NULL

  # removd extra cells
  extra_cells <- setdiff(colnames(combined_mat), rownames(Spatial_D00_all[[1]]))
  combined_mat <- combined_mat[, which(
    !colnames(combined_mat) %in% extra_cells
  )
  ]
  combined_mat <- as.matrix(combined_mat)

  # clean columns of metadata per sample attached sample's name before rbind
  l <- D00
  l <- lapply(
    l,
    function(x) {
      colnames(x@meta.data) <- gsub(
        paste0("_", Assays(x)), "", colnames(x@meta.data)
      )
      x
    }
  )
  D00 <- l

  # first get the list of meta data
  list_of_metadata <- lapply(D00, function(x) x@meta.data)
  # rbind meta data per samples
  meta.data <- do.call("rbind", list_of_metadata)
  write.csv(meta.data, "req_meta_data.csv", row.names = TRUE)

  combined <- CreateSeuratObject(
    counts = combined_mat,
    assay = "scATAC",
    meta.data = meta.data
  )

  combined@meta.data$Clusters <- factor(
    combined@meta.data$Clusters,
    levels = c(
      paste0("C", seq_along(unique(combined@meta.data$Clusters)))
    )
  )

  Spatial_D00 <- list()
  for (i in seq_along(samples)) {
    Spatial_D00[[i]] <- Spatial_D00_all[[i]][colnames(combined), ]
    combined[[paste0(samples[i], "Spatial")]] <- CreateDimReducObject(
      embeddings = Spatial_D00[[i]],
      key = paste0(samples[i], "Spatial_"),
      assay = DefaultAssay(combined)
    )
  }

  # we need to run Variable Features
  combined <- NormalizeData(
    combined,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  combined <- FindVariableFeatures(
    combined,
    selection.method = "vst",
    nfeatures = 2000
  )
  return(combined)
}

#===========================================

combined <- main_func(all)
combined
saveRDS(combined, "/root/neighborhood/combined.rds")

#================================

for (run in runs) {

  run_id <- run[1]
  spatial_path <- run[5]
  print(run_id)

  cidr <- "/root"

  # Make dir /root/neighborhood/data/[run_id] for each run.
  run_dir <- paste0(cidr, "/neighborhood/data/", run_id)
  dir.create(run_dir, recursive = TRUE)

  # Copy spatial dir to run dir.
  file.copy(spatial_path, run_dir, recursive = TRUE)

  # Rename spatial dir from flyte id to 'spatial/'.
  spatial_dir <- list.files(run_dir)
  file.rename(
    paste0(run_dir, "/", spatial_dir),
    paste0(run_dir, "/spatial")
  )

  # Make dir /root/neighborhood/data/[run_id]/ctmat
  ctmat <- paste0(run_dir, "/ctmat")
  dir.create(ctmat, recursive = TRUE)
}

# make dir /root/neighborhood/figures
figs_dir <- "neighborhood/figures"
dir.create(file.path(cidr, figs_dir), recursive = TRUE)

if (genome == "hg38") {
  species <- org.Hs.eg.db
} else if (genome == "mm10") {
  species <- org.Mm.eg.db
} else if (genome == "rnor6") {
  species <- org.Rn.eg.db 
}
species

se_base <- combined
strsplits <- strsplit(rownames(se_base@meta.data), "\\s|#")

fun1 <- function(lst, nthelement) {
  sapply(lst, `[`, 1)
}

fun2 <- function(lst,nthelement) {
  sapply(lst, `[`, 2)
}

se_base@meta.data$sample <- fun1(strsplits)
se_base@meta.data$barcode <- fun2(strsplits)

samples <- sort(unique(se_base@meta.data$sample))

my_vec <- c()

for (sple in sort(samples)) {
  # Head of for-loop
  se_base_D_ordered <- se_base[, which(
    se_base@meta.data$Sample == sple
  )]@meta.data

  req_order <- intersect(
    read.csv(
      paste0(
        "/root/neighborhood/data/",
        sple,
        "/spatial/tissue_positions_list.csv"
      ),
      header = FALSE)$V1,
      se_base_D_ordered$barcode) # Creating new value

  se_base_D_ordered$barcode <- factor(
    se_base_D_ordered$barcode,
    levels = req_order
  )

  se_base_D_ordered <- se_base_D_ordered[order(se_base_D_ordered$barcode), ]

  new_value <- rownames(se_base_D_ordered)

  print(sple)
  my_vec <- c(my_vec, new_value)    # Appending new value to vector
}

saveRDS(se_base, "/root/neighborhood/se_base.rds")

ct <- as.matrix(se_base@assays$scATAC@layers$counts)
colnames(ct)<- colnames(se_base)
rownames(ct)<- rownames(se_base)
ct <- ct[,match(my_vec,colnames(ct))]
meta.data <- se_base@meta.data[match(my_vec,rownames(se_base@meta.data)), ]

se_base_reorder <- CreateSeuratObject(
  counts = as.data.frame(ct),
  assay = "scATAC",
  meta.data = meta.data
)

for (i in sort(names(se_base@reductions))) {
  embeddings <- Embeddings(se_base, reduction = i)

  # All cells in reductions must be in the same order as in the Seurat object
  embeddings <- embeddings[
    match(colnames(se_base_reorder), rownames(embeddings)),
  ]

  # this is for seurat >= v5.0
  se_base_reorder@reductions[[i]] <- CreateDimReducObject(
    embeddings = embeddings,
    key = paste0(i, "_"),
    assay = DefaultAssay(se_base_reorder)
  )
}

se_base_reorder <- NormalizeData(
  se_base_reorder,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
se_base_reorder <- FindVariableFeatures(
  se_base_reorder,
  selection.method = "vst",
  nfeatures = nrow(se_base_reorder)
)
se_base_reorder@meta.data$barcodes1 <- rownames(se_base_reorder@meta.data)

saveRDS(se_base_reorder, "/root/neighborhood/se_base_reorder.rds")

for (sple in samples) {

  D_ <- se_base[, which(se_base@meta.data$sample == sple)]

  write.table(
    data.frame(
      unname(
        mapIds(
          species,
          rownames(D_),
          keytype="SYMBOL",
          column="ENSEMBL",
          multiVals = "first"
        )
      ),
      rownames(D_),
      rep(
        "Gene Expression",
        length(rownames(D_))
      )
    ),
    gzfile(
      paste0(
        "/root/neighborhood/data/", "/", sple, "/ctmat/features.tsv.gz"
      )
    ),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote= FALSE
  )

  # We compress it deleting the .csv
  write.table(
    D_@meta.data$barcode,
    gzfile(
      paste0(
        "/root/neighborhood/data/", "/", sple, "/ctmat/barcodes.tsv.gz"
      )
    ),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )

  writeMM(
    obj = D_@assays$scATAC@layers$counts,
    file=(
      paste0("/root/neighborhood/data/", "/", sple, "/ctmat/matrix.mtx.gz")
    )
  )
}

# also add -1 at the name of barcodes in tissue_position_list.csv file
for (sple in samples) {

  tiss_pos <- read.csv(
    paste0(
      "/root/neighborhood/data/",
      "/",
      sple,
      "/spatial/tissue_positions_list.csv"
    ),
    header = FALSE
  )

  tiss_pos$V1 <- paste0(tiss_pos$V1, "-1")

  names(tiss_pos) <- NULL

  write.csv(
    tiss_pos,
    paste0(
      "/root/neighborhood/data/",
      "/",
      sple,
      "/spatial/tissue_positions_list.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )
}

infoTable = data.frame(
  "samples" = paste0(
    getwd(),
    "/",
    list.files(pattern = "ctmat", recursive = TRUE, include.dirs =TRUE)
  ),
  "spotfiles" = paste0(
    getwd(),
    "/",
    list.files(
      pattern = "tissue_positions_list.csv",
      recursive = TRUE,
      include.dirs = TRUE
    )
  ),
  "imgs" = paste0(
    getwd(),
    "/",
    list.files(
      pattern = "tissue_hires_image.png", recursive = TRUE, include.dirs = TRUE
    )
  ),
  "json" = paste0(
    getwd(),
    "/",
    list.files(
      pattern = "scalefactors_json.json", recursive = TRUE, include.dirs = TRUE
    )
  )
)

se <- InputFromTable(
  infotable = infoTable,
  transpose = FALSE,
  platform =  "Visium"
)
saveRDS(se, "/root/neighborhood/se.rds")

new.names <- colnames(se)
renamed.assay <- RenameCells(se_base_reorder, new.names = new.names)

se_base_reorder <- renamed.assay

head(se@tools$Staffli@imgs)

# Add Staffli 
se_base_reorder@tools$Staffli <- se@tools$Staffli

se_base_reorder <- LoadImages(
  se_base_reorder,
  time.resolve = TRUE,
  verbose = TRUE,
  xdim = 100
)

se_base_reorder@meta.data$seurat_clusters <- se_base_reorder@meta.data$Clusters
se_base_reorder@meta.data$seurat_clusters <- gsub(
  "C",
  "",
  se_base_reorder@meta.data$seurat_clusters
)

c_include <- seq(
  1, max(as.numeric(se_base_reorder@meta.data$seurat_clusters)), 1
)

se_base_reorder@meta.data$sample_id <- se_base_reorder@meta.data$Sample

samples <- sort(unique(se_base_reorder@meta.data$Sample))
samples

for (i in seq_along(samples)) {
  se_base_reorder@meta.data$sample_id <- gsub(
    samples[i],
    paste0("S", i),
    se_base_reorder@meta.data$sample_id
  )
}

saveRDS(se_base_reorder, "se_base_reorder2.rds")

n_clusters <- length(unique(se_base_reorder$Clusters))
cols <- colorRampPalette(
  c(
    "blue", "red", "dark green", "orange", "purple", "brown", "cyan", "yellow"
  )
)(n_clusters)

names(cols) <- paste0("C", seq_len(n_clusters))

######### Define functions #########

### Create Adjacency Matrix

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

minmax_norm <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# RUN NBS ANALYSIS
print("RUN NBS ANALYSIS")
se_base <- se_base_reorder

# ======================================================================
#' RUN NBS ANALYSIS
#'
#' "EXPECTED" VALUES: RegionNeighbours() with permuted cluster IDs for each
#' sample
n_perm <- 50
perm_adj_mat_list <- list()
for (i in 1:n_perm) {

  se_base <- RandomiseClusteringIDs(
    se.object = se_base,
    column.cluster.id = "seurat_clusters",
    column.sample.id = "sample_id",
    random.seed = i
  )
  se_base <- SetIdent(se_base, value = "clusters_perm")
  for (
    column_rm in grep(pattern = "nbs_", colnames(se_base[[]]), value = TRUE)
  ) {
    se_base[[column_rm]] <- NULL
  }
  for (c in c_include) {
    se_base <- RegionNeighbours(
      se_base, id = c, keep.within.id = TRUE, verbose = TRUE
    )
  }
  perm_nbs_df <- se_base[[]][,
  c(
    "clusters_perm",
    grep(pattern = "nbs_", colnames(se_base[[]]), value = TRUE)
  )
]
  perm_adj_mat <- CreateAdjMatrix(
    nbs.df = perm_nbs_df,
    cluster.include = c_include,
    column.clusters.id = "clusters_perm"
  )
  perm_adj_mat_list[[i]] <- perm_adj_mat
}

n_cols <- dim(perm_adj_mat_list[[1]])[1]
perm_adj_mat_avg <- perm_adj_mat_sd <- matrix(0L, nrow = n_cols, ncol = n_cols)
for (i in 1:n_cols) {
  for (j in 1:n_cols) {
    list_ij <- c()
    for (list_n in 1:length(perm_adj_mat_list)) {
      list_ij <- c(list_ij, perm_adj_mat_list[[list_n]][i, j])
    }
    perm_adj_mat_avg[i, j] <- mean(list_ij)
    perm_adj_mat_sd[i, j] <- sd(list_ij)
  }
}

#' Check if random data is normally dist.
perm_adj_mat_df <- data.frame(
  matrix(unlist(perm_adj_mat_list),
  nrow=length(perm_adj_mat_list),
  byrow=TRUE)
)
perm_adj_mat_df$perm_n <- rownames(perm_adj_mat_df)
perm_adj_mat_df_long <- reshape2::melt(perm_adj_mat_df, id.vars = c("perm_n"))

#####
#' OBSERVED VALUES: Run RegionNeighbours()
se_base <- SetIdent(se_base, value = "seurat_clusters")
for (column_rm in grep(pattern = "nbs", colnames(se_base[[]]), value = TRUE)) {
  se_base[[column_rm]] <- NULL
}
for (c in c_include) {
  se_base <- RegionNeighbours(
    se_base,
    id = c,
    keep.within.id = TRUE,
    verbose = TRUE
  )
}

# Create df with nbs info
nbs_df <- se_base[[]][, c(
  "seurat_clusters",
  grep(pattern = "nbs",
  colnames(se_base[[]]),
  value = TRUE
))]

c_count <- nbs_df %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::count()
colnames(c_count) <- c("cluster", "n")
c_count$id <- paste0("C", c_count$cluster)

#' Adjacency matrix
nbs_adjmat <- CreateAdjMatrix(
  nbs.df = nbs_df,
  cluster.include = c_include,
  column.clusters.id = "seurat_clusters"
)

#' Adjacency matrix
nbs_adjmat <- CreateAdjMatrix(
  nbs.df = nbs_df,
  cluster.include = c_include,
  column.clusters.id = "seurat_clusters"
)

#' Normalize method 1
nbs_adjmat_norm <- nbs_adjmat
for (i in seq(1:dim(nbs_adjmat_norm)[1])) {
  for (j in seq(1:dim(nbs_adjmat_norm)[1])) {
    e_sum <- sum(nbs_adjmat_norm[i,i], nbs_adjmat_norm[j,j])
    if (i != j) {
      nbs_adjmat_norm[i, j] <- round(nbs_adjmat_norm[i, j] / e_sum, 4) * 100
    }
  }
}

nbs_adjmat_permscore <- round(
  ((nbs_adjmat - perm_adj_mat_avg) / perm_adj_mat_sd), digits = 3
)

# replace NaN value with zero 
nbs_adjmat_permscore[is.nan(nbs_adjmat_permscore)] <- 0

max_diagonal_val <- max(nbs_adjmat_permscore[nbs_adjmat_permscore > 0])

test <- rep(max_diagonal_val, dim(nbs_adjmat)[1])
names(test) <- names(diag(nbs_adjmat))
test

#' "Normalize"/Standardise method 2: Using permuted values
nbs_adjmat_permscore <- round(
  ((nbs_adjmat - perm_adj_mat_avg) / perm_adj_mat_sd), digits = 3
)
nbs_adjmat_permscore[is.nan(nbs_adjmat_permscore)] <- 0
diag(nbs_adjmat_permscore) <- test

#' Set row/column names
rownames(nbs_adjmat) <- rownames(nbs_adjmat_norm) <- rownames(nbs_adjmat_permscore) <- paste0("C", c_include)
colnames(nbs_adjmat) <- colnames(nbs_adjmat_norm) <- colnames(nbs_adjmat_permscore) <- paste0("C", c_include)

rownames(perm_adj_mat_avg) <- rownames(perm_adj_mat_sd) <- colnames(perm_adj_mat_avg) <- colnames(perm_adj_mat_sd) <- paste0("C", c_include)

#####
#' PER SUBJECT, OBSERVED VALUES

#' Adjacency matrix
nbs_adjmat_list_subject <- CreateAdjMatrixPerSubject(
  nbs.df = nbs_df,
  se.metadata = se_base@meta.data,
  column.subject.id = "sample_id",
  cluster.include = c_include,
  column.clusters.id = "seurat_clusters"
)

#' PER SUBJECT, EXPECTED VALUES
n_perm <- 50
perm_adj_mat_list_subject <- list()
for (i in 1:n_perm) {

  se_base <- RandomiseClusteringIDs(
    se.object = se_base,
    column.cluster.id = "seurat_clusters",
    column.sample.id = "sample_id",
    random.seed = i
  )
  se_base <- SetIdent(se_base, value = "clusters_perm")
  for (
    column_rm in grep(pattern = "nbs_", colnames(se_base[[]]), value = TRUE)
  ) {
    se_base[[column_rm]] <- NULL
  }
  for (c in c_include) {
    se_base <- RegionNeighbours(
      se_base,
      id = c,
      keep.within.id = TRUE,
      verbose = TRUE
    )
  }
  perm_nbs_df <- se_base[[]][,
    c(
      "clusters_perm",
      grep(pattern = "nbs_", colnames(se_base[[]]), value = TRUE)
    )
  ]
  perm_adj_mat_list_subject_id <- CreateAdjMatrixPerSubject(
    nbs.df = perm_nbs_df, 
    se.metadata = se_base@meta.data, 
    cluster.include = c_include, 
    column.clusters.id = "clusters_perm",
    column.subject.id = "sample_id"
  )
  perm_adj_mat_list_subject[[i]] <- perm_adj_mat_list_subject_id
}

n_cols <- dim(perm_adj_mat_list_subject[[1]][[1]])[1]
perm_adj_mat_avg_list_subject <- perm_adj_mat_sd_list_subject <- list()
perm_adj_mat_avg_subject <- perm_adj_mat_sd_subject <- matrix(0L, nrow = n_cols, ncol = n_cols)

for (subject_id in names(perm_adj_mat_list_subject[[1]])) {
  perm_adj_mat_list_subject_id <- list()
  for (i in 1:length(perm_adj_mat_list_subject)) {
    perm_adj_mat_list_subject_id[[i]] <- perm_adj_mat_list_subject[[i]][[subject_id]]
  }
  for (i in 1:n_cols) {
    for (j in 1:n_cols) {
      list_ij <- c()
      for (list_n in 1:length(perm_adj_mat_list_subject_id)) {
        list_ij <- c(list_ij, perm_adj_mat_list_subject_id[[list_n]][i,j])
      }
      perm_adj_mat_avg_subject[i,j] <- mean(list_ij)
      perm_adj_mat_sd_subject[i,j] <- sd(list_ij)
    }
    perm_adj_mat_avg_list_subject[[subject_id]] <- perm_adj_mat_avg_subject
    perm_adj_mat_sd_list_subject[[subject_id]] <- perm_adj_mat_sd_subject
  }
}

#' PER SUBJECT, OBS - EXP_AVG
nbs_adjmat_permscore_subject <- list()
for (subject_id in names(nbs_adjmat_list_subject)) {
  nbs_adjmat_permscore_subject[[subject_id]] <- round(
    (nbs_adjmat_list_subject[[subject_id]] -
      perm_adj_mat_avg_list_subject[[subject_id]]) /
        perm_adj_mat_sd_list_subject[[subject_id]],
    digits = 3
  )
  diag(nbs_adjmat_permscore_subject[[subject_id]]) <- diag(
    nbs_adjmat_list_subject[[subject_id]]
  )
  
  #' Set row/column names
  names_vector <- paste0("C", c_include)
  colnames(nbs_adjmat_permscore_subject[[subject_id]]) <- names_vector
  rownames(nbs_adjmat_permscore_subject[[subject_id]]) <- names_vector
}
# ===================================
#' PLOTS
color_low2 <- "#62376e"

#' PLOT 1
#' Make graph
nbs_adjmat_df <- nbs_adjmat_norm
g <- graph.adjacency(
  nbs_adjmat_df, mode = "undirected", weighted = TRUE, diag = FALSE
)

e_sum <- diag(nbs_adjmat_df)
df <- data.frame(
  name = names(e_sum),
  size = e_sum,
  norm_size = e_sum / max(e_sum) * 40,
  stringsAsFactors = FALSE
)
links <- data.frame(
  source = as_edgelist(g, names = TRUE)[, 1],
  target = as_edgelist(g, names = TRUE)[, 2])

g_df_norm <- data.frame(
  links,
  weight = E(g)$weight,
  weight_minmax = minmax_norm(abs(E(g)$weight))
)

g2 <- graph_from_data_frame(d = links, vertices = df, directed = FALSE)

g2 <- set_edge_attr(g2, "weight", value = E(g)$weight)
E(g2)$width <- E(g2)$weight * 1.5


#' PLOT 3: Permuted score values
#' Make graph
nbs_adjmat_df <- nbs_adjmat_permscore
g <- graph.adjacency(
  nbs_adjmat_df, mode = "undirected", weighted = TRUE, diag = FALSE
)

e_sum <- diag(nbs_adjmat_df)
df <- data.frame(
  name = names(e_sum),
  size = e_sum,
  norm_size = e_sum / max(e_sum) * 40,
  stringsAsFactors = FALSE
)
links <- data.frame(
  source = as_edgelist(g, names = TRUE)[, 1],
  target = as_edgelist(g, names = TRUE)[, 2]
)

g_df_permscore <- data.frame(
  links,
  weight = E(g)$weight,
  weight_minmax = minmax_norm(abs(E(g)$weight))
)

g2 <- graph_from_data_frame(d = links, vertices = df, directed = FALSE)

g2 <- set_edge_attr(g2, "weight", value = minmax_norm(abs(E(g)$weight)))
E(g2)$width <- (E(g2)$weight + 0.1) * 14

fname <- paste0("nbs_analysis.permscore.", project_name, ".pdf")
pdf(
  file = file.path("/root", figs_dir, fname),
  width = 5.5,
  height = 5.5,
  useDingbats = FALSE
)
for (cluster_plot in names(e_sum)) {
  e_pairs <- c()
  for (c in names(e_sum)[!names(e_sum) == cluster_plot]) {
    e_pairs <- c(e_pairs, cluster_plot, c)
  }
  e_weights <- E(g)$weight
  ecol <- rep("grey90", ecount(g2))
  ecol[get.edge.ids(g2, vp = e_pairs)] <- "grey20"
  ecol[e_weights > 0 & ecol == "grey20"] <- color_high
  ecol[e_weights < 0 & ecol == "grey20"] <- color_low2
  ecol[ecol == "grey90"] <- NA
  node_ids <- names(V(g2))
  node_e_widths <- E(g2)$width[get.edge.ids(g2, vp = e_pairs)]
  node_sizes <- data.frame(
    node = node_ids,
    size = V(g2)$norm_size,
    stringsAsFactors = FALSE
  )
  LS <- layout_as_star(
    g2,
    center = cluster_plot,
    order = order(node_sizes$size, decreasing = TRUE)
  )
  plot(
    g2,
    layout = LS,
    vertex.label.family = "Helvetica",
    vertex.label.color = "black",
    vertex.label.cex = 1.5,
    vertex.size = V(g2)$norm_size,
    vertex.color = adjustcolor(V(g2)$color, alpha.f = 1),
    vertex.frame.color = "white",
    edge.color = adjustcolor(ecol, alpha.f = .8),
    edge.curved = 0
  )
}
dev.off()

#########################

#' PLOT HEATMAP
hm_df <- nbs_adjmat_permscore
hm_df <- hm_df[, rev(colnames(hm_df))]

breaksList <- seq(-10, 10, by = .5)

p1 <- pheatmap(
  hm_df,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  breaks = breaksList,
  scale = "none",
  color = colorRampPalette(
    c(color_low3, "white", color_high2))(length(breaksList) - 1
  ),
  border_color = "white",
  na_col = "white"
)

fname <- paste0("nbs_analysis.permscore_hm.", project_name, ".pdf")
pdf(
  file = file.path("/root", figs_dir, fname),
  width = 5,
  height = 4.8,
  useDingbats = FALSE
)
p1

dev.off()

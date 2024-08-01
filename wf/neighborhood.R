library("BPCells")
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
source("wf/scwat_funcs.R")
source("wf/utils.R")


# Set globals -----------------------------------------------------------------
dir.create("/root/neighborhood", showWarnings = FALSE)
setwd("/root/neighborhood")

args <- commandArgs(trailingOnly = TRUE)

project_name <- args[1]
genome <- args[2]

runs <- strsplit(args[3: length(args)], ",")
runs

annotations <- list(
  "mm10" = org.Mm.eg.db, "hg38" = org.Hs.eg.db, "rnor6" = org.Rn.eg.db
)
species <- annotations[[genome]]

# Create list of SeuratObjs, rename cells with run_id to be unique -----
seuratobj_paths <- sapply(runs, `[`, 2)
seurat_list <-  lapply(seuratobj_paths, readRDS)
all <- rename_cells(seurat_list) # from utils.R
saveRDS(all, "/root/neighborhood/all.rds")

# Combine SeuratObject for each run into 'combined' SeuratObject --------------
samples <- find_sample_names(all) # from utils.R

# Extract image coordinates as -(imagecols) | imagerow -----
spatial <- lapply(all, function(x) {
  df <- as.data.frame(x@images[[1]]@coordinates[, c(5, 4)])
  colnames(df) <- paste0("Spatial_", 1:2)
  df$Spatial_2 <- -df$Spatial_2
  df
})

print("Creating combined SeuratObject...")
combined <- combine_objs(all, samples, spatial, project_name)  # from utils.R
saveRDS(combined, "/root/neighborhood/combined.rds")

# Main script -----------------------------------------------------------------

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

se_base <- combined
strsplits <- strsplit(rownames(se_base@meta.data), "\\s|#")

se_base@meta.data$sample <- extract_nth_ele(strsplits, 1) # from utils.R
se_base@meta.data$barcode <- extract_nth_ele(strsplits, 2) # from utils.R

samples <- sort(unique(se_base@meta.data$sample))

my_vec <- c()

for (sample in samples) {

  se_base_D_ordered <- se_base[, which(
    se_base@meta.data$Sample == sample
  )]@meta.data

  # Get intersect of barcodes bt obj metadata and tissue_positions
  req_order <- intersect(
    read.csv(
      paste0(
        "/root/neighborhood/data/", sample, "/spatial/tissue_positions_list.csv"
      ),
      header = FALSE
    )$V1,
    se_base_D_ordered$barcode
  )

  se_base_D_ordered$barcode <- factor(
    se_base_D_ordered$barcode,
    levels = req_order
  )

  se_base_D_ordered <- se_base_D_ordered[order(se_base_D_ordered$barcode), ]

  new_value <- rownames(se_base_D_ordered)

  print(sample)
  my_vec <- c(my_vec, new_value) # Appending new value to vector
}

saveRDS(se_base, "/root/neighborhood/se_base.rds")

ct <- as.matrix(se_base@assays$scATAC@layers$counts)
colnames(ct) <- colnames(se_base)
rownames(ct) <- rownames(se_base)
ct <- ct[, match(my_vec, colnames(ct))]
meta.data <- se_base@meta.data[match(my_vec, rownames(se_base@meta.data)), ]

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

for (sample in samples) {

  D_ <- se_base[, which(se_base@meta.data$sample == sample)]

  write.table(
    data.frame(
      unname(
        mapIds(
          species,
          rownames(D_),
          keytype = "SYMBOL",
          column = "ENSEMBL",
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
        "/root/neighborhood/data/", sample, "/ctmat/features.tsv.gz"
      )
    ),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )

  # We compress it deleting the .csv
  write.table(
    D_@meta.data$barcode,
    gzfile(
      paste0(
        "/root/neighborhood/data/", sample, "/ctmat/barcodes.tsv.gz"
      )
    ),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )

  writeMM(
    obj = D_@assays$scATAC@layers$counts,
    file = (
      paste0("/root/neighborhood/data/", sample, "/ctmat/matrix.mtx.gz")
    )
  )
}

# also add -1 at the name of barcodes in tissue_position_list.csv file
for (sample in samples) {

  tiss_pos <- read.csv(
    paste0(
      "/root/neighborhood/data/", sample, "/spatial/tissue_positions_list.csv"
    ),
    header = FALSE
  )

  tiss_pos$V1 <- paste0(tiss_pos$V1, "-1")

  names(tiss_pos) <- NULL

  write.csv(
    tiss_pos,
    paste0(
      "/root/neighborhood/data/",
      sample,
      "/spatial/tissue_positions_list.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )
}

infoTable <- data.frame(
  "samples" = paste0(
    getwd(),
    "/",
    list.files(pattern = "ctmat", recursive = TRUE, include.dirs = TRUE)
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

# RUN NBS ANALYSIS
print("RUN NBS ANALYSIS")
se_base <- se_base_reorder

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
  matrix(
    unlist(perm_adj_mat_list), nrow = length(perm_adj_mat_list), byrow = TRUE
  )
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
nbs_df <- se_base[[]][,
  c(
    "seurat_clusters",
    grep(pattern = "nbs", colnames(se_base[[]]), value = TRUE)
  )
]

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
    e_sum <- sum(nbs_adjmat_norm[i, i], nbs_adjmat_norm[j, j])
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
        list_ij <- c(list_ij, perm_adj_mat_list_subject_id[[list_n]][i, j])
      }
      perm_adj_mat_avg_subject[i, j] <- mean(list_ij)
      perm_adj_mat_sd_subject[i, j] <- sd(list_ij)
    }
    perm_adj_mat_avg_list_subject[[subject_id]] <- perm_adj_mat_avg_subject
    perm_adj_mat_sd_list_subject[[subject_id]] <- perm_adj_mat_sd_subject
  }
}

#' PER SUBJECT, OBS - EXP_AVG
nbs_adjmat_permscore_subject <- list()
for (subject_id in names(nbs_adjmat_list_subject)) {
  nbs_adjmat_permscore_subject[[subject_id]] <- round(
    (
      nbs_adjmat_list_subject[[subject_id]] -
        perm_adj_mat_avg_list_subject[[subject_id]]
    ) /
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
  target = as_edgelist(g, names = TRUE)[, 2]
)

g_df_norm <- data.frame(
  links,
  weight = E(g)$weight,
  weight_minmax = minmax_norm(abs(E(g)$weight)) # from utils.R
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
  weight_minmax = minmax_norm(abs(E(g)$weight)) # from utils.R
)

g2 <- graph_from_data_frame(d = links, vertices = df, directed = FALSE)

g2 <- set_edge_attr(  # from utils.R
  g2, "weight", value = minmax_norm(abs(E(g)$weight))
)
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

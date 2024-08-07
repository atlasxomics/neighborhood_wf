library("BPCells")
library("dplyr")
library("purrr")
library("Seurat")

chunk_seurat_object <- function(seurat_obj, max_elements = 2^31 - 1) {
  #' Function to chunk SeuratObject into minimum number chunks, such that no
  #' chunck exceeds the integer limit for sparse matrices; written mostly by
  #' GPT4.

  # Calculate the number of cells per chunk
  num_cells <- ncol(seurat_obj)
  num_features <- nrow(seurat_obj)
  max_cells_per_chunk <- floor(max_elements / num_features)

  # Split the cells into chunks
  cell_indices <- split(
    1:num_cells,
    ceiling(seq_along(1:num_cells) / max_cells_per_chunk)
  )

  # Create a list to store the chunks
  seurat_chunks <- list()

  # Loop through each chunk of cells
  for (i in seq_along(cell_indices)) {
    cells <- cell_indices[[i]]
    chunk <- subset(seurat_obj, cells = colnames(seurat_obj)[cells])
    seurat_chunks[[i]] <- chunk
  }

  return(seurat_chunks)
}

extract_nth_ele <- function(lst, n = 1) {
  #' Extract nth element for each item in an array; return as a vector.
  sapply(lst, `[`, n)
}

find_func <- function(tempdir, pattern) {
  list.files(
    path = tempdir,     #' Replace tempdir with the directory you want, followed
    pattern = pattern,  #' by 0 or more characters, then ".csv", and then
    full.names = TRUE,  #' nothing else ($) include the directory in the result.
    recursive = TRUE
  )
}

find_sample_names <- function(seurat_lst) {
  #' Extract list of sample names from list of SeuratObjs.
  sapply(
    seq_along(seurat_lst), function(i) {
      unique(seurat_lst[[i]]@meta.data$Sample)
    }
  )
}

combine_objs <- function(seurat_lst, samples, spatial, project_name) {

  # ----------- Filter SeuratObjs  -----------
  # Remove cell with NA counts from each SeuratObj
  to_remove <- lapply(seurat_lst, function(x) {
    names(which(colSums(is.na(x@assays[[1]]@counts)) > 0))
  })
  filtered <- mapply(
    function(x, y) x[, !colnames(x) %in% y], seurat_lst, to_remove
  )

  # ----------- Prepare coordinates -----------
  # Make combinations of spatial coordinates, return as matrix; if n=1, simply
  # return df as matix.
  spatial_all <- lapply(seq_along(spatial), function(i) {
    tmp <- dplyr::bind_rows(spatial[-i])
    tmp$Spatial_1 <- 0
    tmp$Spatial_2 <- 0
    as.matrix(rbind(spatial[[i]], tmp))
  })

  # ----------- Prepare Metadata -----------
  # remove assay-specific suffixes in obj meta.data (nCount_Spatial -> nCount)
  filtered <- lapply(filtered, function(x) {
    colnames(x@meta.data) <- gsub(
      paste0("_", Seurat::Assays(x)),
      "",
      colnames(x@meta.data)
    )
    x
  })

  # get the list of metadata from seurat objects
  list_of_metadata <- lapply(filtered, function(x) {
    x@meta.data
  })

  # combine meta data
  meta.data <- do.call("rbind", list_of_metadata)
  write.csv(meta.data, "req_meta_data.csv", row.names = TRUE)

  # ----------- Create Combined SeuratObj -----------

  ncells  <- as.double(sum(sapply(filtered, ncol)))
  nfeats <- as.double(nrow(seurat_lst[[1]]))
  matrix_size <- ncells * nfeats

  if (matrix_size < 2^31 - 1) {

    print("Feature matrix less than 2^31 -1...")

    # Convert list of SeuratObjs to list of feat counts as data frames.
    filtered_dfs <- lapply(filtered, function(x) {
      df <- as.data.frame(x@assays[[1]]@counts)
      colnames(df) <- Seurat::Cells(x)
      df$region <- rownames(df)
      df
    })

    # Merge counts dfs, drop region
    combined_mat <- purrr::reduce(filtered_dfs, full_join, by = "region")
    rownames(combined_mat) <- combined_mat$region
    combined_mat$region <- NULL

    combined <- Seurat::CreateSeuratObject(
      counts = combined_mat,
      assay = "scATAC",
      meta.data = meta.data
    )

  } else {

    print("Feature matrix greater than 2^31 -1; using BPCells...")

    # Use BPCells to handle large matrices
    # Convert SeuratObjs to BPCells files
    counts_mat <- c()
    for (i in seq_along(filtered)) {

      path <- paste0(project_name, "/", samples[[i]], "_BP")

      BPCells::write_matrix_dir(
        mat = filtered[[i]]@assays[["Spatial"]]@counts,
        dir = path,
        overwrite = TRUE
      )
      counts_mat[[i]] <- BPCells::open_matrix_dir(dir = path)
    }

    # Create combiend SeuratObject
    combined <- Seurat::CreateSeuratObject(
      counts = counts_mat,
      assay = "scATAC",
      meta.data = meta.data
    )

    combined <- SeuratObject::JoinLayers(combined)
  }

  # Normalize and calculate variable features
  combined <- Seurat::NormalizeData(
    combined, normalization.method = "LogNormalize", scale.factor = 10000
  )

  combined <- Seurat::FindVariableFeatures(
    combined, selection.method = "vst", nfeatures = 2000, slot = "counts"
  )

  # Make clusters factors
  combined@meta.data$Clusters <- factor(
    combined@meta.data$Clusters,
    levels = c(paste0("C", seq_along(unique(combined@meta.data$Clusters))))
  )

  # Add spatial coordinates as embeddings
  # Remove cells that remain in coordinates but not in combined (had NA counts)
  include <- intersect(Seurat::Cells(combined), rownames(spatial_all[[1]]))
  combined <- subset(combined, cells = include)

  # Add embeddings to combined
  for (i in seq_along(samples)) {
    spatial_all[[i]] <- spatial_all[[i]][colnames(combined), ]
    combined[[samples[i]]] <- Seurat::CreateDimReducObject(
      embeddings = spatial_all[[i]],
      key = samples[i],
      assay = Seurat::DefaultAssay(combined)
    )
  }

  return(combined)
}

minmax_norm <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

rename_cells <- function(seurat_list) {
  #' For each SeuratObj is a array of SeuratObjs, rename the cells from barcode
  #' to run_id#barcode-1; ensures cell names are unique when combining
  #' SeuratObjs

  all <-  list()
  for (i in seq_along(seurat_list)) {
    all[[i]] <- seurat_list[[i]]
    all[[i]] <- Seurat::RenameCells(
      all[[i]],
      new.names = paste0(
        unique(all[[i]]@meta.data$Sample),
        "#",
        colnames(all[[i]]),
        "-1"
      )
    )
  }
  return(all)
}
packages <- c(
  "Seurat", "SeuratDisk", "dplyr", "R.utils", "ggplot2", 
  "ggExtra", "RColorBrewer", "openxlsx", "scales", 
  "HGNChelper", "dittoSeq", "harmony", "batchelor", 
  "zellkonverter", "SingleCellExperiment", "data.table"
)
lapply(packages, library, character.only = T)

setwd("/media/senthilkumar/New/scRNA-Seq_Workshop/Bonus_data")

inputDir <- "/media/senthilkumar/New/scRNA-Seq_Workshop/Bonus_data/input/"
outputDir <- "/media/senthilkumar/New/scRNA-Seq_Workshop/Bonus_data/output/"
plotDir <- "/media/senthilkumar/New/scRNA-Seq_Workshop/Bonus_data/plots/"
study_id <- "GSE183276"

gc()
# QC-Count matrix ---------------------------------------------------------

counts <- readRDS("/media/senthilkumar/New/scRNA-Seq_Workshop/Bonus_data/GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Counts_03282022.RDS")
seu_obj <- CreateSeuratObject(counts = counts, assay = "RNA", project = "Kidney_Atlas", min.cells = 3, min.features = 500)

metadata <- fread("/media/senthilkumar/New/scRNA-Seq_Workshop/Bonus_data/GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt.gz")
colnames(metadata)[1] <- "barcodes"
seu_obj <- AddMetaData(seu_obj, metadata = metadata)
seu_obj[["percent.rb"]] <- PercentageFeatureSet(seu_obj, pattern = "^RP[SL]")
remove(counts)
gc()

violin_plot <- function(seurat_object, plotDir, study_id, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")) {
  
  if(!dir.exists(plotDir)){
    dir.create(plotDir, recursive = TRUE)
  }
 
  meta <- as.data.frame(seurat_object@meta.data)
  sample_id <- as.factor(meta$orig.ident)
  
  for (feat in  features){
    
    if (!feat %in% colnames(meta)){
       warning(feat, " not available in meta data ! skipping..")
       next
    }
      cat("plotting feature ", feat, "\n")
    
   
    
      p = ggplot(data = meta, mapping = aes(x = orig.ident, y = .data[[feat]]))+
        geom_violin(trim = TRUE, na.rm = TRUE, alpha = 0.7, fill = "steelblue") +
        geom_boxplot(width = 0.1, na.rm = TRUE, alpha = 0.6, outlier.shape = NA) + 
        labs(x = "Sample", y = feat, title = feat)+ 
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              axis.text.y = element_text(vjust = 0.5))

      plot_name = paste0(feat, "_PreQC_violin_plot.png")
      
      cat("Saving feature ", plot_name," in " ,plotDir,  "\n")
      
      ggsave(filename = plot_name, plot = p, device = png, path = plotDir, dpi = 300, width = 15, height = 9)        
  }
  
}

violin_plot(seu_obj, plotDir = plotDir, study_id = study_id)

density_plot <- function(seurat_object, plotDir, file_prefix) {
  density_df <- data.frame(RNA_Count = log10(seurat_object$nCount_RNA), feature_Count = log10(seurat_object$nFeature_RNA), condition = seurat_object$condition.l1)
  
  p = ggplot(data = density_df, aes(x = RNA_Count, y = feature_Count, colour = condition)) + 
    geom_point(alpha = 0.5, size = 0.5) +
    theme_minimal() +
    theme(legend.position.inside = c(0.05, 0.95), legend.justification = c("left", "top"), legend.key.size = unit(0.5, 'cm')) +
    guides(color = guide_legend(override.aes = list(size = 5))) + # 'size' here controls the symbol size
    labs(x = "log1p(nCount_RNA)", y = "log1p(nFeature_RNA)")
  
  filename = paste0(plotDir,file_prefix,".png")
  
  cat(filename)
  
  
  p <- ggMarginal(p, type = "histogram", fill = "skyblue", bins = 40)
  ggsave(
    filename,
    plot = p,
    width = 8,
    height = 10,
    dpi = 400,
    bg = "white"
  )
}

density_plot(seu_obj, file_prefix = "preQC_Density_plot", plotDir = plotDir)

# Grid search optimal feature ---------------------------------------------

feature_filter <- function(seurat_object, min_feature, max_feature){
  
  subset(seurat_object, subset = nFeature_RNA > min_feature & nFeature_RNA < max_feature )
}

mt_filter <- function(seurat_object, mt_threshold) { 
  subset(seurat_object, subset = percent.mt > mt_threshold)
  
}

mahalanobis_filter <- function(seurat_object, mahalanobis_distance){
  
  mahalanobis_threshold = qchisq(mahalanobis_distance, df = 1)
  cat("Mahalanobis distance for threshold", mahalanobis_distance, "is", round(mahalanobis_threshold,2), "\n")
  
  counts_df <- data.frame(nCount_RNA = seurat_object@meta.data$nCount_RNA)
  
  mahalanobis_values <- mahalanobis(x = counts_df, center = mean(counts_df$nCount_RNA), cov = var(counts_df$nCount_RNA))
  
  seurat_object$mahalanobis_values <- mahalanobis_values
  names_to_keep <- rownames(seurat_object@meta.data)[seurat_object$mahalanobis_values  < mahalanobis_threshold]
  names_length <- (length(names_to_keep))
  cat(names_length, "\n")
  subset(x = seurat_object, cells = names_to_keep)
  
  
}
seurat_filtered <- mahalanobis_filter(seu_obj, mahalanobis_distance = 0.98)
ncol(seurat_filtered) 
min_features = 200
max_features = 7000
mt_threshold = c(20,30)
mahalanobis_threshold = c(0.95, 0.98)

grid_search_features <- function(seurat_object, min_features, max_features, mt_threshold, mahalanobis_threshold) { 
  
    pre_qc_count <- as.data.frame(table(seurat_object$orig.ident))
    colnames(pre_qc_count) <- c("orig.ident", "Counts_pre_QC")
    print(colnames(pre_qc_count))
    cat("moving to QC")
    results_list <- list()
    
    for (mahal_val in mahalanobis_threshold) {
      for (mt in mt_threshold) {
        cat("running grid search with mahalanobis_threshold", mahal_val , "and mt_threshold", mt, "\n")
        
        tmp_obj <- feature_filter(seurat_object, min_features, max_features)
        cat("Ran feature filter \n")
        tmp_obj <- mahalanobis_filter(tmp_obj, mahal_val)
        cat("Ran mahal filter \n")
        tmp_obj <- mt_filter(tmp_obj, mt)
        cat("Ran mt filter \n")
        
        post_qc_count <- as.data.frame(table(tmp_obj$orig.ident))
        post_qc_count
        colnames(post_qc_count) <- c("orig.ident", "Counts_post_QC")
        print(colnames(post_qc_count))
        cat("moving to post QC \n")
        post_qc_count$Counts_post_QC[is.na(post_qc_count$Counts_post_QC)] <- 0
        
        QC_table <- merge(pre_qc_count, post_qc_count,by = "orig.ident", all.x = TRUE )
        QC_table$Count_cells_removed <- (QC_table$Counts_pre_QC -  QC_table$Counts_post_QC)
        QC_table$percent_cells_removed <- round(((QC_table$Counts_post_QC/QC_table$Counts_pre_QC)*100),3)
        
        QC_table$min_feature <- min_features
        QC_table$max_feature <- max_features
        QC_table$mahal_value <- mahal_val
        QC_table$mt_threshold <- mt
        
        Total_cells_pre_qc  <- sum(pre_qc_count$Counts_pre_QC, na.rm = TRUE)
        Total_cells_post_qc <- sum(post_qc_count$Counts_post_QC, na.rm =TRUE)
        Total_count_removed <- sum(QC_table$Count_cells_removed)
        Total_percent_cells_removed <- round(x = mean(QC_table$percent_cells_removed),digits = 3)
        
        

        combo_name <- paste0("mahal_", mahal_val, "_mt_", mt, "_run")
        
        total_row <- data.frame(orig.ident = "Total",
                                Count_cells_removed = Total_count_removed,
                                percent_cells_removed = Total_percent_cells_removed,
                                min_feature = min_features,
                                max_feature = max_features,
                                mahal_value = mahal_val,
                                mt_threshold = mt,
                                Counts_pre_QC = Total_cells_pre_qc,
                                Counts_post_QC = Total_cells_post_qc)
        
        
        QC_table <- rbind(QC_table, total_row)
        
        results_list[[combo_name]] <- QC_table 
        
      
      }
    }
    
    grid_search_results <- do.call(rbind, results_list)
    #return(grid_search_results)
    file_name <- paste0(plotDir,study_id,"_grid_search_results.csv")
    
    write.csv(x = grid_search_results, file = file_name)
    gc()
}


grid_search_features(seurat_object = seu_obj, min_features = 200, max_features = 7000, mt_threshold = c(5, 10, 20), mahalanobis_threshold = c(0.95, 0.98))

seurat_filtered <- feature_filter(seurat_object = seu_obj, min_feature = 200, max_feature = 7000)
seurat_filtered <- mahalanobis_filter(seurat_object = seurat_filtered, mahalanobis_distance = 0.98)
seurat_filtered <- mt_filter(seurat_object = seurat_filtered, mt_threshold = 5)
gc()

violin_plot_post <- function(seurat_object, plotDir, study_id, features = c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
  
  if(!dir.exists(plotDir)){
    dir.create(plotDir, recursive = TRUE)
  }
  
  meta <- as.data.frame(seurat_object@meta.data)
  sample_id <- as.factor(meta$orig.ident)
  
  for (feat in  features){
    
    if (!feat %in% colnames(meta)){
      warning(feat, " not available in meta data ! skipping..")
      next
    }
    cat("plotting feature ", feat, "\n")
    
    
    
    p = ggplot(data = meta, mapping = aes(x = orig.ident, y = .data[[feat]]))+
      geom_violin(trim = TRUE, na.rm = TRUE, alpha = 0.7, fill = "steelblue") +
      geom_boxplot(width = 0.1, na.rm = TRUE, alpha = 0.6, outlier.shape = NA) + 
      labs(x = "Sample", y = feat, title = feat)+ 
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            axis.text.y = element_text(vjust = 0.5))
    
    plot_name = paste0(feat, "_PostQC_violin_plot.png")
    
    cat("Saving feature ", plot_name," in " ,plotDir,  "\n")
    
    ggsave(filename = plot_name, plot = p, device = png, path = plotDir, dpi = 300, width = 15, height = 9)        
  }
  
}

violin_plot_post(seurat_object = seurat_filtered, plotDir = plotDir, study_id = study_id)
density_plot(seurat_object = seurat_filtered, plotDir = plotDir, file_prefix = "PostQC_Density_plot")

seurat_processed <- NormalizeData(seurat_filtered, normalization.method = "LogNormalize")
seurat_processed <- FindVariableFeatures(seurat_processed, selection.method = "vst", nfeatures = 2500)
seurat_processed <- ScaleData(seurat_processed)
seurat_processed <- RunPCA(seurat_processed, npcs = 100, dims=1:100)


stdev_seu <- seurat_processed[["pca"]]@stdev
variance_seu <- (stdev_seu^2)/sum(stdev_seu^2)
percentile_95 <- min(which(cumsum(variance_seu) >= 0.95))
cat("The first", percentile_95, "PC's explain 95 % of total variantion")

ElbowPlot(seurat_processed, ndims = 50, reduction = "pca") + labs(title = "Elbow plot of first 50 PC's from QC processed scRNAseq data")



cumulative_variance_pc40 <- cumsum((variance_seu)* 100)[40]
cat("Elbow plot pleatues at PC 40 and variance explained by first 40 PC's was", round(cumulative_variance_pc40, 3), "%" )
DimHeatmap(seurat_processed, dims = 35:45)

seurat_processed <- RunUMAP(seurat_processed, dims = 1:40)
seurat_processed <- FindNeighbors(seurat_processed, dims = 1:40, reduction = "pca", graph.name = "pca_nn")
seurat_processed <- FindClusters(seurat_processed, graph.name = "pca_nn", cluster.name = "kidney_UMAP", resolution = 0.8)

nclusters_umap <- length(unique(seurat_processed$kidney_UMAP))

orig_ident <- unique(seurat_processed@meta.data$orig.ident)

Hormony_processed <- RunHarmony(seurat_processed, c("orig.ident"), plot_convergence = TRUE)
Hormony_processed <- RunUMAP(Hormony_processed, reduction = "harmony", dims = 1:50)
Hormony_processed <- FindNeighbors(Hormony_processed, reduction = "harmony",  dims = 1:50, graph.name = "kidney_harmony")
Hormony_processed <- FindClusters(Hormony_processed, graph.name = "kidney_harmony", resolution = 0.8, group.name = "kidney_clusters")
gc()


compute_knn_batch_mixing <- function(seurat_object, batch_var="orig.ident", reduction = "pca", dims=1:50, k=20)
  {
  
  if(!batch_var %in% colnames(seurat_object@meta.data)){
    cat("Column ", batch_var, " not found in seurat object")
  }
  
  if(!reduction %in% Reductions(seurat_object)){
    seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
    seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst",  nfeatures = 3000)
    seurat_object <- ScaleData(seurat_object, verbose = FALSE)
    seurat_object <- RunPCA(seurat_object, npcs = max(dims), verbose = FALSE)
  }
  
  labs <- seurat_object@meta.data[[batch_var]]
  
  emb <- Embeddings(seurat_object, reduction = reduction)[,dims, drop = FALSE]
  nn <- FNN::get.knn(emb, k = k)$nn.index
  
  samples <- sapply(seq_len(nrow(nn)), function(i){
    mean(labs[nn[i,]]==labs[i], na.rm = TRUE)
    })
  
  df <- data.frame(batch = labs, nearest_neighbors = samples)
  aggregate(nearest_neighbors ~ batch, df, mean)
  
}

pre <- compute_knn_batch_mixing(seurat_object = seurat_processed, reduction = "pca", batch_var = "orig.ident", dims = 1:50, k = 20)
post <- compute_knn_batch_mixing(seurat_object = Hormony_processed, reduction = "harmony")


data_compared <- rbind(transform(pre, status = "pre"), transform(post, status = "post"))
data_compared
transform(pre, status = "post")

batch_plot <- ggplot(data = data_compared, aes(x = batch, y = nearest_neighbors, fill = factor(status, levels = c("pre", "post")))) + 
  geom_col(position = position_dodge(width = 0.9), width = 0.7) +
         labs(
           title = "Batch Mixing",
           x = "Samples",
           y = "Mean same-batch fraction (higher = stronger batch effect)"
         ) +
         scale_fill_manual(
           values = c(
             "pre" = "#00CED1",     # turquoise
             "post" = "#9370DB"     # new purple shade (example 3rd color)
           )
         ) +
         theme_minimal(base_size = 15) +
         theme(
           axis.text.x = element_text(angle = 45, hjust = 1),
           legend.title = element_blank(),
           plot.title = element_text(face = "bold", hjust = 0.5)
         ) 

+  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                    axis.text.y = element_text(vjust = 0.5))

ggsave(filename = paste0(plotDir,"batch_mixing_comparision.png"), plot = batch_plot, dpi = 300, height = 9, width = 12)

sce <- as.SingleCellExperiment(seurat_processed)
writeH5AD(sce, paste0(outputDir,study_id,"seurat_processed.h5ad"))   
sce <- as.SingleCellExperiment(Hormony_processed)
writeH5AD(sce, paste0(outputDir,study_id,"harmony_processed.h5ad"))   

harmony_umap_coordinates <- Embeddings(Hormony_processed, reduction = "umap")
write.csv(harmony_umap_coordinates, file = paste0(outputDir, study_id, "harmony_umap_coordinates.csv"))
pca_umap_coordinates <- Embeddings(seurat_processed, reduction = "umap")   
write.csv(pca_umap_coordinates, file = paste0(outputDir, study_id, "pca_umap_coordinates.csv"))

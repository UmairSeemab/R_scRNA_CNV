library(Seurat)
setwd("bc/")

load("seurat_obj_2.Robj", verbose = TRUE)

count_matrix <- as.matrix(GetAssayData(seurat_obj[['RNA']], slot = "counts"))

ann <- cbind(colnames(seurat_obj), as.character(seurat_obj@meta.data$seurat_clusters), as.character(seurat_obj@meta.data$stim))
write.table(ann, file = "culster_annotation_for_merged_data.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")
library(infercnv)
infercnv_obj <-CreateInfercnvObject(raw_counts_matrix = count_matrix ,
                                    gene_order_file = "gene_position_file.txt",
                                    annotations_file = "culster_annotation_for_merged_data.txt",
                                    delim = "\t",
                                    ref_group_names = c("3","5","9"))
#?CreateInfercnvObject

options(bitmapType='cairo') 

out_dir <- "results_bc_HMMi6_3_5_9_CbG-F/"

cnv_default <- run(infercnv_obj, cutoff = 0.1,
                   min_cells_per_gene = 3,
                   out_dir = out_dir,
                   cluster_by_groups = FALSE,
                   plot_steps = TRUE,
                   denoise = TRUE,
                   HMM = TRUE,
                   #HMM_type = c("i3"),
                   #HMM_i3_pval = 0.05,
                   num_threads = 16,
                   png_res = 300)

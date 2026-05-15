library(Banksy)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)

library(ggplot2)
library(gridExtra)
library(pals)
library(tidyverse)

seurat_obj <- LoadXenium(
  data.dir = "~/Downloads/XeniumOutputs2",
  assay = "Xenium",
  segmentations = "cell"
)

seurat_obj@active.assay <- "Xenium"

FeatureScatter(
  object = seurat_obj,
  feature1 = "nCount_Xenium",
  feature2 = "nFeature_Xenium",
  raster = F
)

VlnPlot(seurat_obj,
        features = "nFeature_Xenium",
        ncol = 1,
        raster = F,
        pt.size = 0
) +
  xlab("")

VlnPlot(seurat_obj,
        features = "nCount_Xenium",
        ncol = 1, raster = F,
        pt.size = 0
)

seurat_obj@meta.data %>%
  ggplot(aes(x = nFeature_Xenium)) +
  geom_histogram(binwidth = 1, col = "black", fill = "white") +
  theme_bw() +
  ylab("") +
  ggtitle("Distribution of Genes per Cell")

seurat_obj@meta.data %>% ggplot(aes(x = nCount_Xenium)) +
  geom_histogram(binwidth = 1, col = "black", fill = "white") +
  theme_bw() +
  ylab("") +
  ggtitle("Distribution of Transcripts per Cell")

gene_quantiles <- quantile(seurat_obj$nFeature_Xenium, probs = c(0.05, 0.95))

count_quantiles <- quantile(seurat_obj$nCount_Xenium, probs = c(0.05, 0.95))


seurat_obj <- subset(seurat_obj,
                     subset = nFeature_Xenium > gene_quantiles[1] &
                       nFeature_Xenium < gene_quantiles[2] &
                       nCount_Xenium > count_quantiles[1] &
                       nCount_Xenium < count_quantiles[2]
)

seurat_obj@meta.data %>% ggplot(aes(x = nFeature_Xenium)) +
  geom_histogram(binwidth = 1, col = "black", fill = "white") +
  theme_bw() +
  ylab("") +
  ggtitle("Distribution of Genes per Cell")

seurat_obj@meta.data %>% ggplot(aes(x = nCount_Xenium)) +
  geom_histogram(binwidth = 1, col = "black", fill = "white") +
  theme_bw() +
  ylab("") +
  ggtitle("Distribution of Transcripts per Cell")

selected_genes <- c("ABCC11", "ALAS2", "ASCL3", "CA4", "CD1A", "CD70", "CLCA2", "CTSG", 
                    "DIRAS3", "ERBB2", "FOXA1", "GNLY", "HIGD1B", "KCNMA1", "LILRB2", 
                    "MET", "MYC", "PDGFRB", "PRF1", "S100A1", "SLC18A2", "SPIB", "THAP2", "TREM2")

seurat_obj@assays$Xenium$counts %>%
  rowSums() %>%
  broom::tidy() %>%
  rename(Gene = names,Total = x) %>% 
  filter(Gene %in% selected_genes) %>% 
  arrange(-Total) %>% 
  ggplot(aes(reorder(Gene,Total),Total))+
  geom_col(aes(fill=Gene),col="black")+
  coord_flip()+
  theme_minimal()+
  xlab("")+
  ylab("Expression counts")+
  ggtitle("Expression Counts for Genes of interest")+
  scale_fill_viridis_d(alpha = 0.6)

#############################RUN ONLY ONCE, FIXING SPATIAL COORDS##########################################
#seurat_obj %>% Seurat::GetTissueCoordinates() %>% select(x,y) -> spatial_coords
#seurat_obj@meta.data <- seurat_obj@meta.data %>% bind_cols(spatial_coords)

seurat_obj@meta.data %>%
  ggplot(aes(x = x, y = y, col = nCount_Xenium)) +
  geom_point(alpha = 0.5, size = 0.05) +
  labs(
    title = "Spatial Map of Total Transcripts per Cell",
    subtitle = "Colored by Total Transcript Count",
    x = "Spatial X (µm)",
    y = "Spatial Y (µm)",
    color = "Total Transcripts"
  ) +
  theme_minimal() +
  scale_y_reverse() +
  coord_equal() +
  scale_color_viridis_c()

seurat_obj@meta.data %>% 
  ggplot(aes(x = x, y = y, col = nFeature_Xenium)) +
  geom_point(alpha = 0.5, size = 0.05) +
  labs(
    title = "Spatial Map of Gene Diversity per Cell",
    subtitle = "Colored by Number of Genes Detected",
    x = "Spatial X (µm)",
    y = "Spatial Y (µm)",
    color = "Detected Genes"
  ) +
  theme_minimal() +
  coord_equal() +
  scale_y_reverse() +
  scale_color_viridis_c()


seurat_obj <- NormalizeData(seurat_obj, assay = "Xenium")
invisible(gc())

seurat_obj <- FindVariableFeatures(seurat_obj)
invisible(gc())

seurat_obj <- RunBanksy(seurat_obj,
                        lambda = 0.2,
                        dimx = "x",
                        dimy = "y",
                        assay = "Xenium",
                        slot = "data",
                        features = "all",
                        k_geom = 1500,
                        use_agf = T,
                        verbose = T,
                        assay_name = "BANKSY",
                        spatial_mode = "rNN_gauss"
)
invisible(gc())

DefaultAssay(seurat_obj) <- "BANKSY"

seurat_obj <- RunPCA(seurat_obj,
                     assay = "BANKSY",
                     features = rownames(seurat_obj),
                     npcs = 20,
                     verbose = TRUE
)
invisible(gc())

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, verbose = T)
invisible(gc())

seurat_obj <- FindClusters(seurat_obj, resolution = 1.2, group.singletons = T, verbose = T, algorithm = 4)
invisible(gc())

library(FNN) # For get.knn()

compute_ncc_score <- function(coords, cluster_labels, k = 10) {
  
  knn_result <- get.knn(coords, k = k)
  
  matches <- sapply(1:nrow(coords), function(i) {
    neighbor_ids <- knn_result$nn.index[i, ]
    same_cluster <- sum(cluster_labels[neighbor_ids] == cluster_labels[i])
    same_cluster / k
  })
  
  ncc_score <- mean(matches)
  return(ncc_score)
}

compute_ncc_score(
  coords = tibble(
    x = seurat_obj@meta.data$x,
    y = seurat_obj@meta.data$y
  ),
  cluster_labels = seurat_obj@meta.data$BANKSY_snn_res.1.2,
  k = 1500
)

All_Markers <- FindAllMarkers(object = seurat_obj, assay = "BANKSY", test.use = "wilcox")
invisible(gc())

All_Markers %>%
  arrange(-avg_log2FC) %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  write_csv("~/Documents/top10.csv")

driver_genes <- All_Markers %>% 
  filter(gene %in% c("EGFR", "KRAS", "ALK", "MET")) %>% 
  select(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2,p_val)

driver_genes %>% write_csv("~/Documents/EGFRMET.csv")

cluster_annotations <- tibble(
  cluster = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24),
  cluster_annotation = c(
    "Proliferative epithelial cells",
    "Lung epithelial cells",
    "T cells",
    "Club cells or lung secretory cells",
    "Fibroblasts or stromal cells",
    "Lung epithelial cells",
    "Endothelial cells",
    "Macrophages",
    "Tumor epithelial cells",
    "Proliferating cells",
    "Macrophages or monocytes",
    "Club cells or epithelial cells",
    "Dendritic cells",
    "B cells or plasma cells",
    "Lymphoid immune cells",
    "B cells",
    "Macrophages",
    "Stromal cells",
    "Tumor epithelial cells",
    "Mast cells",
    "Smooth muscle cells",
    "Alveolar type I cells",
    "Epithelial cells (uncertain subtype)",
    "NK cells or cytotoxic T cells"
  )
) %>% mutate(cluster = as.factor(cluster))

# RUN THIS CODE ONLY ONCE
#seurat_obj@meta.data <-  seurat_obj@meta.data %>% left_join(cluster_annotations,by = c("BANKSY_snn_res.1.2" = "cluster"))

my_pal <- c(
  "#D81B60", 
  "#F06292", 
  "#1E88E5", 
  "#C2185B", 
  "#AD1457", 
  "#F48FB1",
  "#EC407A", 
  "#43A047", 
  "#8D6E63", 
  "#B0BEC5", 
  "#4CAF50", 
  "#880E4F", 
  "#0288D1", 
  "#039BE5", 
  "#4FC3F7", 
  "#FF8A80",
  "#A1887F", 
  "#0288D1",
  "#81D4FA", 
  "#6D4C41",
  "#66BB6A", 
  "#90A4AE",
  "#FF80AB",
  "#5D4037"
)

seurat_obj@meta.data %>%
  ggplot(aes(x, y, col = BANKSY_snn_res.1.2)) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  theme_minimal() +
  scale_color_manual(values = my_pal) +
  theme(legend.position = "bottom") +
  guides(
    color = guide_legend(
      ncol = 4, nrow = 6, byrow = TRUE,
      override.aes = list(shape = 15, size = 4)
    )
  ) +
  xlab("X axis coordinates") +
  ylab("Y axis coordinates") +
  ggtitle("Cell typing") +
  labs(color = "Cell Type")

seurat_obj@meta.data %>%
  ggplot(aes(x, y, col = cluster_annotation)) +
  geom_point(size = 0.01, show.legend = F) +
  scale_y_reverse() +
  theme_minimal() +
  scale_color_manual(values = my_pal) +
  theme(legend.position = "bottom") +
  guides(
    color = guide_legend(
      ncol = 4, nrow = 6, byrow = TRUE,
      override.aes = list(shape = 15, size = 4)
    )
  ) +
  xlab("X axis coordinates") +
  ylab("Y axis coordinates") +
  ggtitle("Cell typing") +
  facet_wrap(~cluster_annotation, scales = "free") +
  labs(color = "Cell Type")

seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, n.components = 3, verbose = T)
invisible(gc())

seurat_obj@reductions$umap@cell.embeddings %>%
  as_tibble() %>%
  mutate(
    cluster = seurat_obj@meta.data$cluster_annotation,
    cluster_n = seurat_obj@meta.data$BANKSY_snn_res.1.2
  ) %>%
  ggplot(aes(umap_1, umap_2, col = cluster), alpha = 0.1) +
  geom_point() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(
    color = guide_legend(
      ncol = 4, nrow = 6, byrow = TRUE,
      override.aes = list(shape = 15, size = 4)
    )
  ) +
  scale_color_manual(values = my_pal) +
  ggtitle("Uniform Manifold Approximation and Projection ( UMAP ) of our Data") +
  labs(color = "Cell Type")

seurat_obj@reductions$umap@cell.embeddings %>% as_tibble() %>% mutate(cluster = seurat_obj@meta.data$cluster_annotation) -> umap_data

library(plotly)
plot_ly(
  data = umap_data %>% sample_n(10000),
  x = ~umap_1, y = ~umap_2, z = ~umap_3,
  color = ~cluster, 
  #colors = colorRamp(c("blue", "red")),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 2)
)

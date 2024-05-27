if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("BiocManager", quietly = TRUE))
  install.packages('dendextend') # stable CRAN version
  install.packages("ape")

BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(ape)

######### PCA
df <- read.csv("C:/Users/46705/Documents/SpiderSilk/data/post_filtering/filtered_all_fam_pca.csv")
full_df <- df
heatmap_df <- df[, (ncol(df) - 9):ncol(df)]



########## all expression
# df <- read.csv("C:/Users/46705/Documents/SpiderSilk/data/post_filtering/filtered_all_fam.csv")


# full_df <- df[5:nrow(df), ]
# heatmap_df <- df[5:nrow(df), 8:20]
# heatmap_df <- df[5:nrow(df), 8:(ncol(df)-28)]


########## Handpicked Orthogroups
# list_orthogroups <- list('OG0000282','OG0008243','OG0008607', 'OG0008312', 'OG0000264','OG0004717', 'OG0007747','OG0007860', 'OG0004409', 'OG0000189')
# full_df <- df[5:nrow(df), ]
# heatmap_df <- df[5:nrow(df), 8:20]

# heatmap_df <- full_df[, names(full_df) %in% list_orthogroups]

full_df$Family <- factor(full_df$Family)
family_freq <- table(full_df$Family)
top_family_values <- names(sort(family_freq, decreasing = TRUE)[1:6])
top_family_colors <- colorRampPalette(c("blue", "red"))(length(top_family_values))


numeric_df <- apply(heatmap_df, 2, as.numeric)
matrix_data <- as.matrix(numeric_df)

# color scheme
specific_colors <- c("Araneidae" = "#1f77b4", "Tetragnathidae" = "#2ca02c", "Theridiidae"= "#ff7f0e", "Salticidae"= "#d62728", "Thomisidae"= "#8c564b", "Pisauridae" = "#9467bd")
gray_color <- "gray"

all_colors <- sapply(levels(full_df$Family), function(x) {
  if (x %in% names(specific_colors)) {
    specific_colors[[x]]
  } else {
    gray_color
  }
}, USE.NAMES = TRUE) 


clust_object <-hclust(dist(matrix_data))
clusters <- cutree(clust_object, k = 5)
colors = c("red", "blue", "green", "black", "yellow")

# plot(as.phylo(clust_object), type = "fan", tip.color = colors[clusters],
#   label.offset = 1, cex = 0.7)
####

 
cluster_to_heatmap <- 1

# Identify the samples falling within the specified cluster
samples_in_cluster <- which(clusters == cluster_to_heatmap)

# Subset the matrix_data to include only the samples in the specified cluster
subset_matrix_data <- matrix_data[samples_in_cluster, ]
df_sub <- full_df[samples_in_cluster, ]
clust_object <- hclust(dist(subset_matrix_data))

clusters_2 <- cutree(clust_object, k = 5)

row_dend = as.dendrogram(clust_object)
row_dend = color_branches(row_dend, k = 6)

# Create the heatmap using the subset of samples
heatmap_object <- Heatmap(subset_matrix_data,
                          name = "Expression",
                          cluster_columns = FALSE, 
                          row_names_side = "left", 
                          row_title = "Spider individuals",  
                          show_column_names = TRUE,
                          show_row_names = TRUE,
                          cluster_rows = row_dend
)


text_list <-  as.character(df_sub$species_x)


row_anno <- rowAnnotation(Family = df_sub$Family, col= list(Family= all_colors))
#  foo = anno_text(df_sub$species_x)
# , Genus= full_df$outliers_genus, Spicies= full_df$outliers_species)
# clu2 <- draw(heatmap_object + row_anno)

clu1 <- draw(heatmap_object + row_anno)
plot(clu1)

# Define the cluster number you want to focus on
cluster_to_heatmap <- 2

# Identify the samples falling within the specified cluster
samples_in_cluster <- which(clusters_2 == cluster_to_heatmap)

# Subset the matrix_data to include only the samples in the specified cluster
subset_matrix_data_2 <- subset_matrix_data[samples_in_cluster, ]
df_sub_2 <- df_sub[samples_in_cluster, ]
clust_object <- hclust(dist(subset_matrix_data_2))

######

row_dend = as.dendrogram(clust_object)
row_dend = color_branches(row_dend, k = 1)

# Create the heatmap using the subset of samples
heatmap_object <- Heatmap(subset_matrix_data_2,
                          name = "Expression",
                          cluster_columns = FALSE, 
                          row_names_side = "left", 
                          row_title = "Spider individuals",  
                          show_column_names = TRUE,
                          cluster_rows = row_dend
)


row_anno <- rowAnnotation(Family = df_sub_2$Family, col= list(Family= all_colors), Genus = anno_text(df_sub_2$Genus), Spicies = anno_text(df_sub_2$species_x), ID = anno_text(df_sub_2$ID))
# , Genus= full_df$outliers_genus, Spicies= full_df$outliers_species)
clus1_2 <- draw(heatmap_object + row_anno)

plot(clus1_2)


cluster_to_heatmap <- 5
samples_in_cluster <- which(clusters_2 == cluster_to_heatmap)

subset_matrix_data_2 <- subset_matrix_data[samples_in_cluster, ]
df_sub_2 <- df_sub[samples_in_cluster, ]
clust_object <- hclust(dist(subset_matrix_data_2))

######

row_dend = as.dendrogram(clust_object)
row_dend = color_branches(row_dend, k = 1)

heatmap_object <- Heatmap(subset_matrix_data_2,
                          name = "Expression",
                          cluster_columns = FALSE, 
                          row_names_side = "left", 
                          row_title = "Spider individuals",  
                          show_column_names = TRUE,
                          cluster_rows = row_dend
)

row_anno <- rowAnnotation(Family = df_sub_2$Family, col= list(Family= all_colors),Genus = anno_text(df_sub_2$Genus), Spicies = anno_text(df_sub_2$species_x))
# , Genus= full_df$outliers_genus, Spicies= full_df$outliers_species)
clus1_5 <- draw(heatmap_object + row_anno)


#pdf("C:/Users/46705/Documents/SpiderSilk/output/Heatmaps_R/Ara_outliers.pdf")
# Your plotting code here, e.g.,

plot(big)
plot(clu2)
plot(clu3)
dev.off()

pdf("C:/Users/46705/Documents/SpiderSilk/output/Heatmaps_R/AraGroup_outliers.pdf")
# Your plotting code here, e.g.,

plot(big)
plot(clu1)
plot(clus1_2)
plot(clus1_5)

dev.off()

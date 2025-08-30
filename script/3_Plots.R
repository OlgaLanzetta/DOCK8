
library("pheatmap")  
library("gplots")    
library("viridis")   
library("SCpubr")
library("ComplexHeatmap")  
library(Seurat)
library(pheatmap)
library(RColorBrewer)

# Figure 5

# PANEL C
# load the integrated object

unique_colors <- c(
  "0"  =  "#619CFF",  # Bright Blue
  "1"  = "#33a02c",  # Bright Green
  "2"  = "#984ea3",  # Bright Purple
  "3"  = "#ff7f00",  # Bright Red
  "4"  = "#fb9a99",  # Bright Orange
  "5"  = "#a65628",  # Chestnut Brown
  "6"  = "#E76BF3",  # Hot Pink
  "7"  = "#999999",  # Medium Gray
  "8"  = "#dede00",  # Neon Yellow
  "9"  = "darkslateblue",  # Bright Cyan
  "10" = "gold",  # Blue (alt)
  "11" = "#00BCD8",  # Green (alt)
  "12" = "#e41a1c",  # Coral Pink
  "13" = "#cab2d6",  # Lilac
  "14" = "#00FF00",  # Deep Purple
  "15" = "#c49c94",  # Bright Pink
  "16" = "#0000FF",  # Cocoa
  "17" = "#00FFFF",  # Pale Aqua
  "18" = "#911EB4",  # Light Yellow
  "19" = "#fdb462"   # Bright Peach
)
 DimPlot(integrated, reduction = "umap_harmony", group.by = "harmony_0.6", cols = unique_colors, label = T, pt.size = 0.7, label.size = 10, split.by = "orig.ident")

# PANEL D
#load the Cd4+ Foxp3+ subset

# Define gene list
genes <- c("Rorc", "Ikzf3", "Maf", "Ahr", "Batf", "Ctla4", "Gata3", "Il1rl1", "Ikzf2", "Setd2", "Il23r","Il10", "Irf4", "Tnfrsf18", "Gzmb", "Prdm1", "Itgae", "Icos", "Tgfbr1", "Samd3", "Ccr6", "Ccr9", "Ccr4", "Gpr15", "Bcl6")

# Check if genes exist in Seurat object
genes_available <- genes[genes %in% rownames(subset)]
genes_missing <- genes[!genes %in% rownames(subset)]

# Print missing genes (if any)
if (length(genes_missing) > 0) {
  message("The following genes were not found in the dataset: ", paste(genes_missing, collapse = ", "))
}

# Extract average expression for available genes
avg_expr <- AverageExpression(subset, features = genes_available, group.by = "orig.ident")$RNA

# Normalize by Z-score (Standardization)
scaled_expr <- t(scale(t(avg_expr)))  
# **Limit Z-score range to [-2, 2]**
scaled_expr[scaled_expr > 2] <- 2   
scaled_expr[scaled_expr < -2] <- -2  

# Define color palette (Red-White-Blue)
heatmap_colors <- colorRampPalette(c("navy", "ivory", "firebrick"))(100)
pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)
p1 <- pheatmap(scaled_expr, 
               color = pal, 
               cluster_rows = TRUE,   
               cluster_cols = FALSE,  
               scale = "none",        
               show_rownames = TRUE, 
               show_colnames = TRUE, 
               border_color = NA,
               main = "Gene Expression Heatmap (Z-Score Normalized)")
p1

# PANEL E

# Load the subset

# Define custom colors
color_map <- c("DOCK8_WT" = "gray",
               "DOCK8_H962R" = "blue",               
               "DOCK8_N1987Y" = "purple",
               "DOCK8_stat3flox" = "orange")

# Reorder factor levels to place WT first
subset$orig.ident <- factor(subset$orig.ident,
                                   levels = c("DOCK8_WT", "DOCK8_H962R", "DOCK8_N1987Y", "DOCK8_stat3flox"))

# Define genes to plot
genes <- c("Rorc","Ikzf3","Maf","Ahr","Batf","Ctla4","Gata3","Il1rl1", "Ikzf2")

# Generate violin plots for each gene separately
p_list <- VlnPlot(subset, features = genes, group.by = "orig.ident", pt.size = 0.5, combine = FALSE)

# Apply custom colors to each plot
p_list <- lapply(p_list, function(plot) {
  plot + scale_fill_manual(values = color_map)
})

# Combine all plots using patchwork
final_plot <- patchwork::wrap_plots(p_list) 

# Display the final plot
print(final_plot)

# Extended data FIGURE 1
# Panel A

# load the integrated object

unique_colors <- c(
  "0"  =  "#619CFF",  # Bright Blue
  "1"  = "#33a02c",  # Bright Green
  "2"  = "#984ea3",  # Bright Purple
  "3"  = "#ff7f00",  # Bright Red
  "4"  = "#fb9a99",  # Bright Orange
  "5"  = "#a65628",  # Chestnut Brown
  "6"  = "#E76BF3",  # Hot Pink
  "7"  = "#999999",  # Medium Gray
  "8"  = "#dede00",  # Neon Yellow
  "9"  = "darkslateblue",  # Bright Cyan
  "10" = "gold",  # Blue (alt)
  "11" = "#00BCD8",  # Green (alt)
  "12" = "#e41a1c",  # Coral Pink
  "13" = "#cab2d6",  # Lilac
  "14" = "#00FF00",  # Deep Purple
  "15" = "#c49c94",  # Bright Pink
  "16" = "#0000FF",  # Cocoa
  "17" = "#00FFFF",  # Pale Aqua
  "18" = "#911EB4",  # Light Yellow
  "19" = "#fdb462"   # Bright Peach
)
 DimPlot(integrated, reduction = "umap_harmony", group.by = "harmony_0.6", cols = unique_colors, label = T, pt.size = 0.7, label.size = 10)

# Panel B
# load the integrated object
# Define custom colors for each cluster
unique_colors <- c(
  "0"  = "#619CFF", "1"  = "#33a02c", "2"  = "#984ea3", "3"  = "#ff7f00", "4"  = "#fb9a99",
  "5"  = "#a65628", "6"  = "#E76BF3", "7"  = "#999999", "8"  = "#dede00", "9"  = "darkslateblue",
  "10" = "gold",     "11" = "#00BCD8", "12" = "#e41a1c", "13" = "#cab2d6", "14" = "#00FF00",
  "15" = "#c49c94", "16" = "#0000FF", "17" = "#00FFFF", "18" = "#911EB4", "19" = "#fdb462"
)

# Set default assay
DefaultAssay(integrated) <- "RNA"

# Marker genes to visualize
phenotype_TF_marker <- c(
  "Cd3e", "Cd3d", "Cd3g", "Cd8a", "Cd4", "Trgv1", "Cd19", "Ms4a1", "Ncr1", "H2-Ab1",
  "H2-Eb2", "Klrb1", "Itgam", "Itgax", "Cd14", "Lyz2", "Il7r", "Itgae", "Cd69", "Pdcd1",
  "Foxp3", "Eomes", "Tbx21", "Gata3", "Rorc", "Rora", "Il23r"
)


# Extract expression and metadata
df <- FetchData(integrated, vars = c("harmony_0.6", phenotype_TF_marker)) %>%
  mutate(cell = rownames(.)) %>%
  pivot_longer(cols = all_of(phenotype_TF_marker), names_to = "gene", values_to = "expression")

# Ensure harmony_0.6 is a factor with proper levels

cluster_levels <- as.character(sort(as.numeric(unique(df$harmony_0.6))))

df$harmony_0.6 <- factor(df$harmony_0.6)
df$gene <- factor(df$gene, levels = phenotype_TF_marker)

# Match color to clusters
stopifnot(length(cluster_levels) <= length(unique_colors))
cluster_colors <- setNames(unique_colors[seq_along(cluster_levels)], cluster_levels)

# Save as high-resolution TIFF
tiff("phenotype_TF_marker.tiff", width = 8, height = 12, units = "in", res = 600)

# Plot
p <- ggplot(df, aes(x = harmony_0.6, y = expression, fill = harmony_0.6)) +
  geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 1) +
  geom_jitter(size = 0.01, width = 0.25, alpha = 0.2) +
  facet_wrap(~gene, nrow = length(phenotype_TF_marker), strip.position = "right", scales = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 0),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    strip.text.y.right = element_text(angle = 0, size = 10, face = "italic"),
    strip.placement = "outside",
    legend.position = "none",
    panel.grid = element_blank()
  )

print(p)
dev.off()

#######################

cytokines <- c("Ifng", "Il4", "Il5", "Il17a", "Il22", "Tnf", "Il10", "Il21",
               "Gzma", "Gzmb", "Fasl", "Prf1", "Tgfb1", "Il1b", "Ccl3", "Ccl4")

# Extract expression and metadata
df <- FetchData(integrated, vars = c("harmony_0.6", cytokines)) %>%
  mutate(cell = rownames(.)) %>%
  pivot_longer(cols = all_of(cytokines), names_to = "gene", values_to = "expression")

# Ensure harmony_0.6 is a factor with proper levels

cluster_levels <- as.character(sort(as.numeric(unique(df$harmony_0.6))))

df$harmony_0.6 <- factor(df$harmony_0.6)
df$gene <- factor(df$gene, levels = cytokines)

# Match color to clusters
stopifnot(length(cluster_levels) <= length(unique_colors))
cluster_colors <- setNames(unique_colors[seq_along(cluster_levels)], cluster_levels)

# Save as high-resolution TIFF
tiff("cytokines.tiff", width = 8, height = 12, units = "in", res = 600)

# Plot
p <- ggplot(df, aes(x = harmony_0.6, y = expression, fill = harmony_0.6)) +
  geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 1) +
  geom_jitter(size = 0.01, width = 0.25, alpha = 0.2) +
  facet_wrap(~gene, nrow = length(cytokines), strip.position = "right", scales = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 0),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    strip.text.y.right = element_text(angle = 0, size = 10, face = "italic"),
    strip.placement = "outside",
    legend.position = "none",
    panel.grid = element_blank()
  )

print(p)
dev.off()


# Panel C
# load the integrated object

# Reorder factor levels to place WT first
integrated$orig.ident <- factor(integrated$orig.ident,
                                   levels = c("DOCK8_WT", "DOCK8_H962R", "DOCK8_N1987Y", "DOCK8_stat3flox"))

# Define custom colors
color_map <- c("DOCK8_WT" = "gray",
               "DOCK8_H962R" = "blue",
               
               "DOCK8_N1987Y" = "purple",
               "DOCK8_stat3flox" = "orange")

genes <- c("Ifng", "Tnf", "Il1b", "Gzmb", "Prf1", "Ccl4", "Il17a", "Il5")
# Generate violin plots for each gene separately
p_list <- VlnPlot(integrated, features = genes, group.by = "orig.ident", pt.size = 0.5, combine = FALSE)

# Apply custom colors to each plot
p_list <- lapply(p_list, function(plot) {
  plot + scale_fill_manual(values = color_map)
})

# Combine all plots using patchwork
final_plot <- patchwork::wrap_plots(p_list) 

# Display the final plot
print(final_plot)

# Extended data FIGURE 2 
# Panel C
# load the integrated object


# Define custom colors
color_map <- c("DOCK8_WT" = "gray",
               "DOCK8_H962R" = "blue",
               
               "DOCK8_N1987Y" = "purple",
               "DOCK8_stat3flox" = "orange")
# Reorder factor levels to place WT first
integrated$orig.ident <- factor(integrated$orig.ident,
                                   levels = c("DOCK8_WT", "DOCK8_H962R", "DOCK8_N1987Y", "DOCK8_stat3flox"))

genes <- c("Il2ra", "Il2rb", "Tigit", "Stat5a", "Stat5b", "Pik3cd", "Mitor", "Socs1", "Socs3", "Cish", "Vim", "Tnfrsf18")

# Generate violin plots for each gene separately
p_list <- VlnPlot(integrated, features = genes, group.by = "orig.ident", pt.size = 0.5, combine = FALSE)

# Apply custom colors to each plot
p_list <- lapply(p_list, function(plot) {
  plot + scale_fill_manual(values = color_map)
})

# Combine all plots using patchwork
final_plot <- patchwork::wrap_plots(p_list) 

# Display the final plot
print(final_plot)



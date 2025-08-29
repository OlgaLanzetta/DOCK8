library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(SingleR)
library(DoubletFinder)
library(glue)
library(tibble)
atlas <- celldex::ImmGenData()

so_list= readRDS(file = "/path/rds/so_list.rds")

library(SeuratWrappers)

# set seed numbers
seeds <- c(000, 111, 222, 333, 444, 555, 666, 777, 888, 999)

cd4_counts_list <- vector("list", length(seeds))
names(cd4_counts_list) <- seeds

# loop over seeds
results_list <- lapply(seeds, function(seed) {
  set.seed(seed)
  message(paste("Running with seed:", seed))
  
  # Step 1: ALRA for each object in the list
  so_list_seeded <- lapply(so_list, function(x) {
    x <- RunALRA(x)
    features.plot <- "Foxp3"
    invisible(gc())
    
    # Plot RNA assay
    DefaultAssay(x) <- "RNA"
    plot1 <- FeaturePlot(x, features.plot, ncol = 2)
    
    # Plot ALRA assay
    DefaultAssay(x) <- "alra"
    plot2 <- FeaturePlot(x, features.plot, ncol = 2, cols = c("lightgrey", "red"))
    print(CombinePlots(list(plot1, plot2), ncol = 1))
    
    return(x)
  })
  
  # Count total Cd4+ cells across objects
  total_cd4 <- sum(sapply(so_list_seeded, function(x) {
    cd4_expression <- GetAssayData(object = x, assay = "alra", slot = "data")["Cd4", ]
    length(which(cd4_expression > 0))
  }))
    
  cd4_counts_list[[as.character(seed)]] <<- total_cd4
  
  # Step 2: Subset Foxp3+ within Cd4+ and preprocess
  so_list_seeded <- lapply(so_list_seeded, function(x) {
    DefaultAssay(x) <- "alra"
    
    # Filter Cd4+ cells
    cd4_expression <- GetAssayData(object = x, assay = "alra", slot = "data")["Cd4", ]
    pos_ids_cd4 <- names(which(cd4_expression > 0))
    x <- subset(x, cells = pos_ids_cd4)

    # Within Cd4+, filter Foxp3+ cells
    foxp3_expression <- GetAssayData(object = x, assay = "alra", slot = "data")["Foxp3", ]
    pos_ids_foxp3 <- names(which(foxp3_expression > 0))
    x <- subset(x, cells = pos_ids_foxp3)
  
    # Standard preprocessing
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    HVG <- VariableFeatures(x)
    all.genes <- rownames(x)
    x <- ScaleData(x, features = all.genes)
    x <- RunPCA(x, features = HVG, verbose = FALSE)
    
    print(ElbowPlot(x, ndims = 30))
    print(DimPlot(x, reduction = "pca"))
    
    return(x)
  })
  
  # Step 3: Integration
  combined <- merge(x = so_list_seeded[[1]], y = so_list_seeded[2:length(so_list_seeded)])
  combined[["alra5"]] <- as(object = combined[["alra"]], Class = "Assay5")
  DefaultAssay(combined) <- "RNA"
  
  combined <- NormalizeData(combined)
  combined <- FindVariableFeatures(combined, nfeatures = 3000)
  combined <- ScaleData(combined)
  combined <- RunPCA(combined)
  
  obj <- IntegrateLayers(
    object = combined,
    method = HarmonyIntegration,
    new.reduction = "harmony",
    verbose = FALSE
  )
  
  # Save integrated object
  saveRDS(obj, file = paste0("/path/rds/seurat_integrated_alra_harmony_cd4_foxp3_10seed", seed, ".rds"))
  return(obj)
})

# PART 2 DE ANALYSIS

library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(readr)
library(ggrepel)
library(glue)
library(VennDiagram)
library(gridExtra)
library(UpSetR)
library(grid)

# === Initial settings ===
wd <- "/path/rds/"
setwd(wd)

rds_files <- list.files(path = wd, pattern = "seurat_integrated_alra_harmony_cd4_foxp3_10seed.*\\.rds$", full.names = TRUE)

cell_counts <- sapply(rds_files, function(file) {
  obj <- readRDS(file)       # Read Seurat object
  ncol(obj)                  # Number of cells
})

cell_counts_df <- data.frame(
  file = basename(rds_files),
  n_cells = cell_counts
)

print(cell_counts_df)

padj.cutoff <- 0.05
lfc.cutoff <- log2(1.2)

wd <- "/path/"
setwd(wd)
dir.create("Figure_alra_DE_seed_foxp3_cd4_seed10", showWarnings = FALSE)
dir.create("Figure_alra_DE_seed_foxp3_cd4_seed10/Venn", showWarnings = FALSE)

# === Lists for saving DEGs and Venn ===
up_genes_list <- list()
down_genes_list <- list()
up_genes_table <- list()
down_genes_table <- list()
all_list <- list()
all_genes_list <- list()

# === Comparisons ===
comparisons <- list(
  H962R_vs_WT = c("DOCK8_H962R", "DOCK8_WT"),
  N1987Y_vs_WT = c("DOCK8_N1987Y", "DOCK8_WT"),
  STAT3Flox_vs_WT = c("DOCK8_stat3flox", "DOCK8_WT")
)

# === Main function ===
analyze_seurat_DE <- function(file_path) {
  message(glue("Processing: {basename(file_path)}"))
  
  seurat.integrated <- readRDS(file_path)
  Idents(seurat.integrated) <- seurat.integrated$orig.ident
  DefaultAssay(seurat.integrated) <- "alra"
  
  # Select expressed genes (1st quartile > 1)
  assay_data <- as.matrix(seurat.integrated@assays$alra@data)
  assay_1stq <- apply(assay_data, 1, quantile)[2,]
  selected_genes <- names(assay_1stq[assay_1stq > 1])
  
  seed_name <- tools::file_path_sans_ext(basename(file_path))
  
  for (comp_name in names(comparisons)) {
    ident.1 <- comparisons[[comp_name]][1]
    ident.2 <- comparisons[[comp_name]][2]
    
    message(glue("  DE: {ident.1} vs {ident.2}"))
    
    table_res <- tryCatch({
      FindMarkers(
        object = seurat.integrated,
        ident.1 = ident.1,
        ident.2 = ident.2,
        assay = "alra",
        slot = "data",
        features = selected_genes,
        test.use = "wilcox",
        min.pct = -Inf,
        min.cells.group = 0,
        logfc.threshold = -Inf
      )
    }, error = function(e) {
      warning(glue("  Skipping {comp_name}: {e$message}"))
      return(NULL)
    })
    
    if (!is.null(table_res)) {
      results <- table_res %>% rownames_to_column("genes")
      
      results$gene_type <- case_when(
        results$avg_log2FC > lfc.cutoff & results$p_val_adj < padj.cutoff ~ "up",
        results$avg_log2FC < -lfc.cutoff & results$p_val_adj < padj.cutoff ~ "down",
        TRUE ~ "ns"
      )
      
      # Save up/down tables
      up_table <- results %>% filter(gene_type == "up")
      down_table <- results %>% filter(gene_type == "down") 
      
      up_genes_table[[comp_name]][[seed_name]] <<- up_table
      down_genes_table[[comp_name]][[seed_name]] <<- down_table
      
      # Save up/down for Venn
      up_genes <- results %>% filter(gene_type == "up") %>% pull(genes)
      down_genes <- results %>% filter(gene_type == "down") %>% pull(genes)
      
      up_genes_list[[comp_name]][[seed_name]] <<- up_genes
      down_genes_list[[comp_name]][[seed_name]] <<- down_genes
 
      all_list[[comp_name]][[seed_name]] <<- results
      all_genes_list[[comp_name]][[seed_name]] <<- results %>% pull(genes)
      
      # Volcano plot
      volcano_plot <- ggplot(results, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = gene_type)) +
        geom_point() +
        geom_hline(yintercept = -log10(padj.cutoff), linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(-lfc.cutoff, lfc.cutoff), linetype = "dashed", color = "red") +
        scale_color_manual(values = c("up" = "salmon", "down" = "turquoise", "ns" = "grey")) +
        labs(
          title = glue("{basename(file_path)}: {ident.1} vs {ident.2}"),
          x = "log2(Fold Change)",
          y = "-log10(Adjusted P-value)",
          colour = "Expression Change"
        ) +
        theme_classic()
      
      pdf_name <- glue("Figure_alra_DE_seed_foxp3_cd4_seed10/{seed_name}_{comp_name}.pdf")
      ggsave(pdf_name, plot = volcano_plot, width = 8, height = 6)
    }
  }
}

# === Analyze all files ===
for (rds_file in rds_files) {
  analyze_seurat_DE(rds_file)
}

plot_upset_for_type <- function(gene_list, label) {
  for (comp in names(gene_list)) {
    comp_sets <- gene_list[[comp]]
    
    if (length(comp_sets) >= 2) {
      all_genes <- unique(unlist(comp_sets))
      
      upset_data <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
      for (set_name in names(comp_sets)) {
        upset_data[[set_name]] <- as.integer(all_genes %in% comp_sets[[set_name]])
      }
      
      rownames(upset_data) <- upset_data$gene
      upset_data$gene <- NULL
      
      pdf(glue("Figure_alra_DE_seed_foxp3_cd4_seed10/Venn/{comp}_{label}_upset.pdf"), 
          width = 10, height = 6, onefile = FALSE)
      print(
        upset(
          upset_data,
          sets = names(comp_sets),
          order.by = "freq",
          mainbar.y.label = "Intersection Size",
          sets.x.label = "Set Size"
        )
      )
      grid.text(glue("{comp} ({label})"), 
                x = 0.5, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"))
      dev.off()
    }
  }
}

plot_upset_for_type(up_genes_list, "up")
plot_upset_for_type(down_genes_list, "down")

# === STAT3 vs WT ===
common_up_genes_STAT3 <- Reduce(intersect, up_genes_list$STAT3Flox_vs_WT)
combined_up_stat3 <- bind_rows(up_genes_table$STAT3Flox_vs_WT, .id = "seed") %>% dplyr::filter(genes %in% common_up_genes_STAT3)
summary_up_STAT3 <- combined_up_stat3 %>%
  group_by(genes) %>%
  summarise(
    max_padj = max(p_val_adj),
    min_log2FC = min(avg_log2FC),
    .groups = "drop"
  ) %>%
  arrange(max_padj)

common_dw_genes_STAT3 <- Reduce(intersect, down_genes_list$STAT3Flox_vs_WT)
combined_dw_stat3 <- bind_rows(down_genes_table$STAT3Flox_vs_WT, .id = "seed") %>% dplyr::filter(genes %in% common_dw_genes_STAT3)
summary_dw_STAT3 <- combined_dw_stat3 %>%
  group_by(genes) %>%
  summarise(
    max_padj = max(p_val_adj),
    min_log2FC =  avg_log2FC[which.min(abs(avg_log2FC))],
    .groups = "drop"
  ) %>%
  arrange(max_padj)

common_genes_STAT3 <- Reduce(intersect, all_genes_list$STAT3Flox_vs_WT)
combined_stat3 <- bind_rows(all_list$STAT3Flox_vs_WT, .id = "seed") %>% dplyr::filter(genes %in% common_genes_STAT3)
summary_STAT3 <- combined_stat3 %>%
  group_by(genes) %>%
  summarise(
    max_padj = max(p_val_adj),
    min_log2FC =  avg_log2FC[which.min(abs(avg_log2FC))],
    .groups = "drop"
  ) %>%
  arrange(max_padj)

summary_STAT3$gene_type <- case_when(
  summary_STAT3$min_log2FC > lfc.cutoff & summary_STAT3$max_padj < padj.cutoff ~ "up",
  summary_STAT3$min_log2FC < -lfc.cutoff & summary_STAT3$max_padj < padj.cutoff ~ "down",
  TRUE ~ "ns"
)

volcano_plot_stat3 <- ggplot(summary_STAT3, aes(x = min_log2FC, y = -log10(max_padj), colour = gene_type)) +
  geom_point() +
  geom_hline(yintercept = -log10(padj.cutoff), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-lfc.cutoff, lfc.cutoff), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("up" = "salmon", "down" = "turquoise", "ns" = "grey")) +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-value)",
    colour = "Expression Change"
  ) + 
  theme_classic()

# === H962R vs WT ===
common_genes_h962R <- Reduce(intersect, all_genes_list$H962R_vs_WT)
combined_h962r <- bind_rows(all_list$H962R_vs_WT, .id = "seed") %>% dplyr::filter(genes %in% common_genes_h962R)
summary_h962r <- combined_h962r %>%
  group_by(genes) %>%
  summarise(
    max_padj = max(p_val_adj),
    min_log2FC = avg_log2FC[which.min(abs(avg_log2FC))],
    .groups = "drop"
  ) %>%
  arrange(max_padj)

summary_h962r$gene_type <- case_when(
  summary_h962r$min_log2FC > lfc.cutoff & summary_h962r$max_padj < padj.cutoff ~ "up",
  summary_h962r$min_log2FC < -lfc.cutoff & summary_h962r$max_padj < padj.cutoff ~ "down",
  TRUE ~ "ns"
)

volcano_plot_h962r <- ggplot(summary_h962r, aes(x = min_log2FC, y = -log10(max_padj), colour = gene_type)) +
  geom_point() +
  geom_hline(yintercept = -log10(padj.cutoff), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-lfc.cutoff, lfc.cutoff), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("up" = "salmon", "down" = "turquoise", "ns" = "grey")) +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-value)",
    colour = "Expression Change"
  ) +
  theme_classic()

# === N1987Y vs WT ===
common_genes_N1987Y <- Reduce(intersect, all_genes_list$N1987Y_vs_WT)
combined_N1987Y <- bind_rows(all_list$N1987Y_vs_WT, .id = "seed") %>% dplyr::filter(genes %in% common_genes_N1987Y)
summary_N1987Y <- combined_N1987Y %>%
  group_by(genes) %>%
  summarise(
    max_padj = max(p_val_adj),
    min_log2FC =  avg_log2FC[which.min(abs(avg_log2FC))],
    .groups = "drop"
  ) %>%
  arrange(max_padj)

summary_N1987Y$gene_type <- case_when(
  summary_N1987Y$min_log2FC > lfc.cutoff & summary_N1987Y$max_padj < padj.cutoff ~ "up",
  summary_N1987Y$min_log2FC < -lfc.cutoff & summary_N1987Y$max_padj < padj.cutoff ~ "down",
  TRUE ~ "ns"
)

volcano_plot_N1987Y <- ggplot(summary_N1987Y, aes(x = min_log2FC, y = -log10(max_padj), colour = gene_type)) +
  geom_point() +
  geom_hline(yintercept = -log10(padj.cutoff), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-lfc.cutoff, lfc.cutoff), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("up" = "salmon", "down" = "turquoise", "ns" = "grey")) +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-value)",
    colour = "Expression Change"
  ) +
  theme_classic()

# === Save results ===
dir.create("/path/DE_genes_alra_Foxp3_cd4/10_seed_result/")
write.csv(summary_STAT3,
          file= ("/path/DE_genes_alra_Foxp3_cd4/10_seed_result/stat3_DE_cd4_foxp3_alra_10seeds.csv"))
write.csv(summary_h962r,
          file= ("/path/DE_genes_alra_Foxp3_cd4/10_seed_result/h962R_DE_cd4_foxp3_alra_10seeds.csv"))
write.csv(summary_N1987Y,
          file= ("/path/DE_genes_alra_Foxp3_cd4/10_seed_result/n1987y_DE_cd4_foxp3_alra_10seeds.csv"))

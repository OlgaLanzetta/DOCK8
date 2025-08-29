# Figure 5

# PANEL C
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


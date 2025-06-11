```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")  
install.packages(c("Seurat", "tidyverse", "dplyr", "ggplot2", "patchwork", "ggrepel", "ggforce", "RColorBrewer", "circlize", "GGally", "data.table")) 

library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(ggforce)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(GGally)
library(data.table)
library(stats)
library(utils)
library(grDevices)
library(grid)
devtools::install_github("002377/TfuncViz2")
library(TfuncViz)
```

# example
```r
seurat_obj <- ...  # yours Seurat object
results <- analyze_tcell_functionality(seurat_obj)
plot <- plot_tfunc_bubble(results)
print(plot)
```
#pipeline
```r
## Using default gene sets
results <- analyze_tcell_functionality(seurat_obj)

# Using custom gene sets

my_gene_sets <- list(
  Cytotoxic = c("GZMB", "PRF1", "GNLY"),
  Helper = c("IL2", "IL4", "IL21"),
  Regulatory = c("FOXP3", "IL10")
)
results <- analyze_tcell_functionality(seurat_obj, gene_sets = my_gene_sets)

# View specific plots
print(results$plots$bubble)
print(results$plots$vector)

# Save all plots to directory
results <- analyze_tcell_functionality(seurat_obj, output_dir = "plots/")

# Extract the updated Seurat object
seurat_obj <- results$seurat_obj
```

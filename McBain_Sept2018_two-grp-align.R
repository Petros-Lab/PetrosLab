#Remove packages and objects
pkg = names(sessionInfo()$otherPkgs)
pkgs = paste('package:', pkg, sep = "")
lapply(pkgs, detach, character.only = TRUE, unload = TRUE)
rm(list=ls(all=TRUE))

require(tidyverse)
#require(plotly)
library(Seurat)
library(dplyr)

#Sys.setenv("plotly_username"="USERNAME")
#Sys.setenv("plotly_api_key"="KEY")

#set.seed(123)

#make sure everything groovy
search() #find attached datasets, objects, etc
######

getwd()
list.files()

# change sample names as needed
sampleDir = c("WTGFPP", "MUTGFPP")

baseDir = "."
matDirs = file.path(baseDir, sampleDir, "outs/filtered_gene_bc_matrices/mm10_GFP")
matDirs

# Load the PBMC dataset
mats = list()
conds = list()
for (i in seq_along(sampleDir)) {
  mats[[i]] <- Read10X(data.dir = matDirs[i])
  conds[[i]] <- CreateSeuratObject(raw.data = mats[[i]], min.cells = 3, min.genes = 200, project = sampleDir[i])
}
names(conds) = tolower(sampleDir)
list2env(conds ,.GlobalEnv)
# rm(mats)


# prep functions
qc_fast = function(x){
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = x@data), value = TRUE)
  percent.mito <- Matrix::colSums(x@raw.data[mito.genes, ]) / Matrix::colSums(x@raw.data)

  # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
  x <- AddMetaData(object = x, metadata = percent.mito, col.name = "percent.mito")
  x
}

plot_fast = function(x) {
  VlnPlot(object = x, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  par(mfrow = c(1, 2))
  GenePlot(object = x, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = x, gene1 = "nUMI", gene2 = "nGene")
}

prep_fast = function(x,y){
  x@meta.data$stim <- y
  x <- FilterCells(x, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
  x <- NormalizeData(x)
  x <- ScaleData(x, display.progress = F)
}


#prep control/WT data
wtgfpp = qc_fast(wtgfpp)
VlnPlot(object = wtgfpp, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
plot_fast(wtgfpp)
wtgfpp = prep_fast(wtgfpp, "CTRL")


#prep stim/mut data
mutgfpp = qc_fast(mutgfpp)
VlnPlot(object = mutgfpp, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
plot_fast(mutgfpp)
mutgfpp = prep_fast(mutgfpp, "STIM")


# Gene selection for input to CCA
wtgfpp <- FindVariableGenes(wtgfpp, do.plot = F)
mutgfpp <- FindVariableGenes(mutgfpp, do.plot = F)
g.1 <- head(rownames(wtgfpp@hvg.info), 1000)
g.2 <- head(rownames(mutgfpp@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(wtgfpp@scale.data))
genes.use <- intersect(genes.use, rownames(mutgfpp@scale.data))

neuro.combined <- RunCCA(wtgfpp, mutgfpp, genes.use = genes.use, num.cc = 30, add.cell.id1 = "wt", add.cell.id2 = "mut")

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = neuro.combined, reduction.use = "cca", group.by = "stim", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = neuro.combined, features.plot = "CC1", group.by = "stim", 
              do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = neuro.combined, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)

p3 <- MetageneBicorPlot(neuro.combined, grouping.var = "stim", dims.eval = 1:30, 
                        display.progress = FALSE)

DimHeatmap(object = neuro.combined, reduction.type = "cca", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)


#Align subspaces
neuro.combined <- AlignSubspace(neuro.combined, reduction.type = "cca", grouping.var = "stim", 
                                 dims.align = 1:20)

p1 <- VlnPlot(object = neuro.combined, features.plot = "ACC1", group.by = "stim", 
              do.return = TRUE)
p2 <- VlnPlot(object = neuro.combined, features.plot = "ACC2", group.by = "stim", 
              do.return = TRUE)
plot_grid(p1, p2)

# t-SNE and Clustering
neuro.combined <- RunTSNE(neuro.combined, reduction.use = "cca.aligned", dims.use = 1:20, 
                           do.fast = T)
neuro.combined <- FindClusters(neuro.combined, reduction.type = "cca.aligned", 
                                resolution = 0.6, dims.use = 1:20)

p1 <- TSNEPlot(neuro.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(neuro.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)

g0.markers <- FindConservedMarkers(neuro.combined, ident.1 = 0, grouping.var = "stim", 
                                   print.bar = FALSE)
head(g0.markers)

wt.subset <- SubsetData(object = neuro.combined, 
                        cells.use = rownames(neuro.combined@meta.data)[neuro.combined@meta.data$orig.ident %in% c("WTGFPP")])

mut.subset <- SubsetData(object = neuro.combined, 
                         cells.use = rownames(neuro.combined@meta.data)[neuro.combined@meta.data$orig.ident %in% c("MUTGFPP")])


genes_reporter = c("EGFP")
genes_nsc = c("Gfap", "Nes", "Vim", "Pax6", "Sox2")
genes_astro = c("Slc1a3", "Slc1a2", "Glul", "S100b", "Aldh1l1")
genes_oligo = c("Olig1", "Olig2", "Cldn11", "Mbp", "Sox10")
genes_microglia = c("Itgam", "Aif1")
genes_npc = c("Ascl1", "Dcx")
genes_immature = c("Tubb3", "Tbr1")
genes_mature = c("Rbfox3", "Calb1", "Dlg4")
genes_glut = c("Slc17a7", "Slc17a6", "Grin1", "Grin2b", "Gls")
genes_inter = c("Slc6a1", "Gabrb1", "Gabrb2", "Gad1", "Gad2")
genes_dopa = c("Th", "Pitx3", "Nr4a2", "Slc6a3", "Kcnj6")

VlnPlot(object = neuro.combined, features.plot = genes_reporter, size.x.use = 8, x.lab.rot = TRUE)
VlnPlot(object = neuro.combined, features.plot = genes_nsc, size.x.use = 8, x.lab.rot = TRUE)
VlnPlot(object = neuro.combined, features.plot = genes_astro, size.x.use = 8, x.lab.rot = TRUE)
VlnPlot(object = neuro.combined, features.plot = genes_oligo, size.x.use = 8, x.lab.rot = TRUE)
VlnPlot(object = neuro.combined, features.plot = genes_microglia, size.x.use = 8, x.lab.rot = TRUE)
VlnPlot(object = neuro.combined, features.plot = genes_npc, size.x.use = 8, x.lab.rot = TRUE)
VlnPlot(object = neuro.combined, features.plot = genes_immature, size.x.use = 8, x.lab.rot = TRUE)
VlnPlot(object = neuro.combined, features.plot = genes_mature, size.x.use = 8, x.lab.rot = TRUE)
VlnPlot(object = neuro.combined, features.plot = genes_glut, size.x.use = 8, x.lab.rot = TRUE)
VlnPlot(object = neuro.combined, features.plot = genes_inter, size.x.use = 8, x.lab.rot = TRUE)
VlnPlot(object = neuro.combined, features.plot = genes_dopa, size.x.use = 8, x.lab.rot = TRUE)


genes_subtypes = c("Htr3a", "Vip", "Calb2", "Reln", "Cck")
v1 = VlnPlot(object = wt.subset, features.plot = genes_subtypes)
v2 = VlnPlot(object = mut.subset, features.plot = genes_subtypes)
plot_grid(v1, v2)

f1 = FeaturePlot(object = wt.subset, features.plot = genes_subtypes, cols.use = c("grey", "blue"), reduction.use = "tsne")
f2 = FeaturePlot(object = mut.subset, features.plot = genes_subtypes, cols.use = c("grey", "red"), reduction.use = "tsne")

Groups
# 1-11 = interneurons, 12 = immature/glutamatergic neurons, 13 = astrocytes, 14 oligos, 15 ???? --> Reln and CCK appear to
# increase in muts, no dopa clusters

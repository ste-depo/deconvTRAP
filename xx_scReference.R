library(Seurat)

integrated_only_brain <- readRDS(
  '../../../external_datasets/scRNA_Brain/fromCluster/integrated_only_brain.rds')
integrated_only_brain <- integrated_only_brain[
  ,integrated_only_brain$cell_names != 'unknown']

gexp_brain <- AverageExpression(integrated_only_brain, 
                                group.by = 'cell_names')$RNA
gspec_brain <- gexp_brain/rowSums(gexp_brain)
saveRDS(gspec_brain, file = '_rdata/xx_gspec_brain.rds')
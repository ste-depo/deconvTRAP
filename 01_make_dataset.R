library(DESeq2)
library(pheatmap)

source('../common/00_generic_functions.R')

# load spinal dataset
load(#rowdata_sc, coldata_sc, cts_sc, cg, 
  file='../common/_rdata/cts_coldata_spinalcord.RData')

# Bonanomi dataset
idx <- coldata_sc$condition %in% c('EC_TRAP_90_WT', 'spinal_cord_whole_tissue_90_WT')
spinal_dds <- DESeqDataSetFromMatrix(countData = cts_sc[,idx],
                                     colData = coldata_sc[idx,],
                                     design = ~ selection)

# Cleuren dataset
cleuren_counts <- loadCounts('../_raw_counts/TRAP_Cleuren_counts.gz')
cleuren_rowdata <- getRowdata('../_raw_counts/TRAP_Cleuren_counts.gz')
cleuren_coldata <- read.table('../_raw_counts/TRAP_Cleuren_metadata.txt', 
                              sep=',', header = TRUE)
rownames(cleuren_coldata) <- cleuren_coldata$Run
cleuren_coldata <- cleuren_coldata[colnames(cleuren_counts),]
#--------- remove un-needed samples
idx <- cleuren_coldata$tissue != 'Blood' & 
  !cleuren_coldata$Run %in% paste0('SRR1025088', 0:5)
cleuren_coldata <- cleuren_coldata[idx,]
cleuren_counts <- cleuren_counts[,idx]
cleuren_dds <- DESeqDataSetFromMatrix(cleuren_counts, 
                                      cleuren_coldata, ~ source_name)


## create DSeq2 dataset
all_counts <- cbind(counts(spinal_dds), counts(cleuren_dds))
all_coldata <- rbind(
  data.frame(tissue = 'SpinalCord',
             selection = colData(spinal_dds)$selection,
             dataset='Bonanomi',
             row.names = colnames(spinal_dds)),
  data.frame(tissue = colData(cleuren_dds)$tissue,
             selection = ifelse(colData(cleuren_dds)$molecule_subtype == "total mRNA", 'whole_tissue', 'TRAP'),
             dataset='Cleuren',
             row.names = colnames(cleuren_dds)))
all_coldata$libsize <- colSums(all_counts)
all_dds <- DESeqDataSetFromMatrix(all_counts, all_coldata, ~ tissue + selection)
stopifnot(identical(rownames(rowdata_sc), rownames(all_dds)))
rowData(all_dds) <- rowdata_sc
all_dds <- estimateSizeFactors(all_dds)
all_vsd <- vst(all_dds)

saveRDS(all_dds, file = '_rdata/01_CleurenBonanomi_raw.rds')

# make palette
tissue_palette <- ggplotColours(6)
names(tissue_palette) <- levels(colData(all_dds)$tissue)
tissue_palette['SpinalCord'] <- 'dodgerblue3'

## PCA
pca_raw_dsets <- plotPCA.san(all_vsd, intgroup = c('dataset'), 
                             fix.coord = FALSE, labels = F) + 
  theme_bw() 
pca_raw <- plotPCA.san(all_vsd, intgroup = 'tissue', 
                       pchgroup = 'selection', fix.coord = FALSE, labels = F) + 
  theme_bw() + 
  scale_shape_manual(values = c('TRAP'=16,'whole_tissue'=2)) + 
  scale_color_manual(values = tissue_palette)

pdf('_panel_figures/01_rawPCA_byCondition.pdf', 
    width = 5, height = 4)
pca_raw
dev.off()

pdf('_panel_figures/01_rawPCA_byDataset.pdf', 
    width = 5, height = 4)
pca_raw_dsets
dev.off()

## plot 500 main features before deconvolution
ntop <- 500
rv <- rowVars(assay(all_vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
phOut <- pheatmap(assay(all_vsd)[select, ], show_rownames = F, 
                  show_colnames=F, treeheight_row = 0,
                  annotation_col = data.frame(colData(all_vsd))[,c('tissue','selection')],
                  annotation_colors = list(tissue=tissue_palette, 
                                           selection=c(TRAP='grey40', whole_tissue='grey80')),
                  file='_panel_figures/01_rawHeatmap.pdf', height=3.5, width=6)

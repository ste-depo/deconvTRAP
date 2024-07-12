library(DESeq2)
library(pheatmap)
library(enrichR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(patchwork)
library(egg)
library(tidyverse)
library(reshape2)

source('../common/00_generic_functions.R')

# make palette
dtrap_dds <- readRDS(file = '_rdata/02_CleurenBonanomi_deconv.rds')
tissue_palette <- ggplotColours(6)
names(tissue_palette) <- levels(colData(dtrap_dds)$tissue)
tissue_palette['SpinalCord'] <- 'dodgerblue3'

# get degs using specificity data

mine_degs <- function(res, p_adj=1e-3, l2fc=2, spec=.5) {
  list(
    down = res %>% filter(padj<p_adj, log2FoldChange<(-l2fc), 
                          Macrophages_VS_Endothelium < (-spec)) %>% 
      rownames(),
    up = res %>% filter(padj<p_adj, log2FoldChange>l2fc, 
                        Fibroblasts_VS_Endothelium < (-spec)) %>% 
      rownames()
  )
}

mine_degs_df <- function(res, p_adj=1e-3, l2fc=2, spec=.5) {
  list(
    down = res %>% filter(padj<p_adj, log2FoldChange<(-l2fc), 
                          Macrophages_VS_Endothelium < (-spec)) %>% 
      mutate(log10_padj = log10(padj)) %>% 
      dplyr::select(starts_with('fpkm'), log10_padj, 
                    log2FoldChange, Macrophages_VS_Endothelium),
    up = res %>% filter(padj<p_adj, log2FoldChange>l2fc, 
                        Fibroblasts_VS_Endothelium < (-spec)) %>% 
      mutate(log10_padj = log10(padj)) %>% 
      dplyr::select(starts_with('fpkm'), log10_padj, 
                    log2FoldChange, Fibroblasts_VS_Endothelium)
  )
}

###################################################################
########### call markers on the deconvoluted dataset ##############
###################################################################

dtrap_dds <- readRDS(file = '_rdata/02_CleurenBonanomi_deconv.rds')
dtrap_cleuren_markers_dds <- dtrap_dds[,dtrap_dds$dataset != 'Bonanomi' &
                                         dtrap_dds$selection == 'TRAP']
dtrap_cleuren_markers_dds <- DESeqDataSetFromMatrix(counts(dtrap_cleuren_markers_dds), 
                                                    colData(dtrap_cleuren_markers_dds), 
                                                    ~ tissue)
rowData(dtrap_cleuren_markers_dds) <- rowData(dtrap_dds)
#-------------------------
coldata <- colData(dtrap_cleuren_markers_dds)
#-------------------------
for(tissue in unique(coldata$tissue)) {
  coldata[[tissue]] <- 1*(coldata$tissue == tissue)
}
coldata$libsize <- colSums(counts(dtrap_cleuren_markers_dds))
dtrap_cleuren_markers_dds <- DESeqDataSetFromMatrix(counts(dtrap_cleuren_markers_dds), coldata, design = ~ tissue)
rowData(dtrap_cleuren_markers_dds) <- rowData(dtrap_dds)[rownames(dtrap_cleuren_markers_dds),]
sizeFactors(dtrap_cleuren_markers_dds) <- dtrap_cleuren_markers_dds$libsize/median(dtrap_cleuren_markers_dds$libsize)
dtrap_cleuren_markers_dds <- DESeq(dtrap_cleuren_markers_dds)
rowData(dtrap_cleuren_markers_dds)$baseLog10Mean <- log10(rowData(dtrap_cleuren_markers_dds)$baseMean)
rowData(dtrap_cleuren_markers_dds)$baseLog10Mean[!is.finite(rowData(dtrap_cleuren_markers_dds)$baseLog10Mean)] <- NA
dtrap_markers_vsd <- vst(dtrap_cleuren_markers_dds)
#-------------------------
design(dtrap_cleuren_markers_dds) <- ~ Brain
Brain_degs_deconv  <- makeStrictComparison('Brain_1_vs_0', dtrap_cleuren_markers_dds, counts_filter = 20, fpkms_filter = 3, f = all)
design(dtrap_cleuren_markers_dds) <- ~ Heart
Heart_degs_deconv  <- makeStrictComparison('Heart_1_vs_0', dtrap_cleuren_markers_dds, counts_filter = 20, fpkms_filter = 3, f = all)
design(dtrap_cleuren_markers_dds) <- ~ Kidney
Kidney_degs_deconv  <- makeStrictComparison('Kidney_1_vs_0', dtrap_cleuren_markers_dds, counts_filter = 20, fpkms_filter = 3, f = all)
design(dtrap_cleuren_markers_dds) <- ~ Liver
Liver_degs_deconv  <- makeStrictComparison('Liver_1_vs_0', dtrap_cleuren_markers_dds, counts_filter = 20, fpkms_filter = 3, f = all)
design(dtrap_cleuren_markers_dds) <- ~ Lung
Lung_degs_deconv  <- makeStrictComparison('Lung_1_vs_0', dtrap_cleuren_markers_dds, counts_filter = 20, fpkms_filter = 3, f = all)

res_deconv <- list(
  Brain=Brain_degs_deconv,
  Heart=Heart_degs_deconv,
  Kidney=Kidney_degs_deconv,
  Liver=Liver_degs_deconv,
  Lung=Lung_degs_deconv
)

saveRDS(res_deconv, file = '_rdata/04_res_deconv.rds')
res_deconv <- readRDS(file = '_rdata/04_res_deconv.rds')

cleuren_heatmap <- 
  mypheatcts(dtrap_cleuren_markers_dds, # row_ann = 'baseLog10Mean',
            genotypic_res = res_deconv, annont_clustering = 10, cluster_cols = F, 
            annot_colors = list(tissue=tissue_palette),
            only.pos = TRUE, log2fc_thr_min = 2, padj_thr_min = 1e-3,
            treeheight_col = 0, treeheight_row = 0,
            border_color = NA, cutree_rows = 5,
            col_ann = c('tissue'), show_rownames = F, show_colnames = F, 
#            filename = '_panel_figures/04_heatmap_deconv.pdf', 
            height = 4.4, width = 4.4,
            annotation_legend = TRUE)


marker_genes_deconv <- cleuren_heatmap$clusters_gene_list
names(marker_genes_deconv) <- 
  unique(apply(cbind(ifelse(cleuren_heatmap$report$Lung==1, 'Lung', ''),
                     ifelse(cleuren_heatmap$report$Liver==1, 'Liver', ''),
                     ifelse(cleuren_heatmap$report$Kidney==1, 'Kidney', ''),
                     ifelse(cleuren_heatmap$report$Heart==1, 'Heart', ''),
                     ifelse(cleuren_heatmap$report$Brain==1, 'Brain', '')), 
               1, paste, collapse=''))
saveRDS(marker_genes_deconv, file = '_rdata/04_marker_genes_deconv.rds')

# setEnrichrSite("Enrichr")
# dbs <- 'ARCHS4_Tissues'
# deconv_mark_ARCHS4_Tissues <- lapply(1:length(marker_genes_deconv), function(i)
#   dbsres <- enrichr(marker_genes_deconv[[i]], dbs)[[dbs]]
# )

setEnrichrSite("Enrichr")
dbs <- 'GO_Biological_Process_2023'
deconv_mark_GO_Biological_Process_2023 <- lapply(1:length(marker_genes_deconv), function(i)
  dbsres <- enrichr(marker_genes_deconv[[i]], dbs)[[dbs]]
)
names(deconv_mark_GO_Biological_Process_2023) <- 
  names(marker_genes_deconv)
deconv_mark_GO_Biological_Process_2023 <- 
  lapply(deconv_mark_GO_Biological_Process_2023, function(x) {
    x$Term <- sapply(strsplit(x$Term, ' \\('), '[', 1)
    x$P.value <- log10(x$P.value)
    return(x)
  })

for(tissue in names(marker_genes_deconv)) {
  gg <- plotEnrich(deconv_mark_GO_Biological_Process_2023[[tissue]], 
             showTerms = 10)
  ggsave(paste0('_panel_figures/04_BP_',tissue,'_deconv.pdf'), gg,
         height = 5, width = 6)
}

##########################################################
########### call markers on the RAW dataset ##############
##########################################################

all_dds <- readRDS(file = '_rdata/01_CleurenBonanomi_raw.rds')
trap_cleuren_markers_dds <- all_dds[,all_dds$dataset != 'Bonanomi' &
                                      all_dds$selection == 'TRAP']
trap_cleuren_markers_dds <- DESeqDataSetFromMatrix(counts(trap_cleuren_markers_dds), 
                                                   colData(trap_cleuren_markers_dds), 
                                                   ~ tissue)
rowData(trap_cleuren_markers_dds) <- rowData(all_dds)
#-------------------------
coldata <- colData(trap_cleuren_markers_dds)
#-------------------------
for(tissue in unique(coldata$tissue)) {
  coldata[[tissue]] <- 1*(coldata$tissue == tissue)
}
coldata$libsize <- colSums(counts(trap_cleuren_markers_dds))
trap_cleuren_markers_dds <- DESeqDataSetFromMatrix(counts(trap_cleuren_markers_dds), coldata, design = ~ tissue)
rowData(trap_cleuren_markers_dds) <- rowData(all_dds)[rownames(trap_cleuren_markers_dds),]
sizeFactors(trap_cleuren_markers_dds) <- trap_cleuren_markers_dds$libsize/median(trap_cleuren_markers_dds$libsize)
trap_cleuren_markers_dds <- DESeq(trap_cleuren_markers_dds)
rowData(trap_cleuren_markers_dds)$baseLog10Mean <- log10(rowData(trap_cleuren_markers_dds)$baseMean)
rowData(trap_cleuren_markers_dds)$baseLog10Mean[!is.finite(rowData(trap_cleuren_markers_dds)$baseLog10Mean)] <- NA
dtrap_markers_vsd <- vst(trap_cleuren_markers_dds)
#-------------------------
design(trap_cleuren_markers_dds) <- ~ Brain
Brain_degs_raw  <- makeStrictComparison('Brain_1_vs_0', trap_cleuren_markers_dds, counts_filter = 20, fpkms_filter = 3, f = all)
design(trap_cleuren_markers_dds) <- ~ Heart
Heart_degs_raw  <- makeStrictComparison('Heart_1_vs_0', trap_cleuren_markers_dds, counts_filter = 20, fpkms_filter = 3, f = all)
design(trap_cleuren_markers_dds) <- ~ Kidney
Kidney_degs_raw  <- makeStrictComparison('Kidney_1_vs_0', trap_cleuren_markers_dds, counts_filter = 20, fpkms_filter = 3, f = all)
design(trap_cleuren_markers_dds) <- ~ Liver
Liver_degs_raw  <- makeStrictComparison('Liver_1_vs_0', trap_cleuren_markers_dds, counts_filter = 20, fpkms_filter = 3, f = all)
design(trap_cleuren_markers_dds) <- ~ Lung
Lung_degs_raw  <- makeStrictComparison('Lung_1_vs_0', trap_cleuren_markers_dds, counts_filter = 20, fpkms_filter = 3, f = all)

res_raw <- list(
  Brain=Brain_degs_raw,
  Heart=Heart_degs_raw,
  Kidney=Kidney_degs_raw,
  Liver=Liver_degs_raw,
  Lung=Lung_degs_raw
)

saveRDS(res_raw, file = '_rdata/04_res_raw.rds')
res_raw <- readRDS(file = '_rdata/04_res_raw.rds')

cleuren_heatmap <- 
  mypheatcts(trap_cleuren_markers_dds, # row_ann = 'baseLog10Mean',
            genotypic_res = res_raw, annont_clustering = 10, cluster_cols = F, 
            annot_colors = list(tissue=tissue_palette),
            only.pos = TRUE, log2fc_thr_min = 2, padj_thr_min = 1e-3,
            treeheight_col = 0, treeheight_row = 0,
            border_color = NA, cutree_rows = 5,
            col_ann = c('tissue'), show_rownames = F, show_colnames = F, 
#            filename = '_panel_figures/04_heatmap_raw.pdf', 
            height = 4.4, width = 4.4,
            annotation_legend = TRUE)


marker_genes_raw <- cleuren_heatmap$clusters_gene_list
names(marker_genes_raw) <- 
  unique(apply(cbind(ifelse(cleuren_heatmap$report$Lung==1, 'Lung', ''),
                     ifelse(cleuren_heatmap$report$Liver==1, 'Liver', ''),
                     ifelse(cleuren_heatmap$report$Kidney==1, 'Kidney', ''),
                     ifelse(cleuren_heatmap$report$Heart==1, 'Heart', ''),
                     ifelse(cleuren_heatmap$report$Brain==1, 'Brain', '')), 
               1, paste, collapse=''))
saveRDS(marker_genes_raw, file = '_rdata/04_marker_genes_raw.rds')

# setEnrichrSite("Enrichr")
# dbs <- 'ARCHS4_Tissues'
# deconv_mark_ARCHS4_Tissues <- lapply(1:length(marker_genes_raw), function(i)
#   dbsres <- enrichr(marker_genes_raw[[i]], dbs)[[dbs]]
# )

setEnrichrSite("Enrichr")
dbs <- 'GO_Biological_Process_2023'
raw_mark_GO_Biological_Process_2023 <- lapply(1:length(marker_genes_raw), function(i)
  dbsres <- enrichr(marker_genes_raw[[i]], dbs)[[dbs]]
)

names(raw_mark_GO_Biological_Process_2023) <- 
  names(marker_genes_raw)
raw_mark_GO_Biological_Process_2023 <- 
  lapply(raw_mark_GO_Biological_Process_2023, function(x) {
    x$Term <- sapply(strsplit(x$Term, ' \\('), '[', 1)
    x$P.value <- log10(x$P.value)
    return(x)
  })

for(tissue in names(marker_genes_deconv)) {
  gg <- plotEnrich(raw_mark_GO_Biological_Process_2023[[tissue]], 
                   showTerms = 10)
  ggsave(paste0('_panel_figures/04_BP_',tissue,'_raw.pdf'), gg,
         height = 5, width = 6)
}

###################################################################################
###### print out markers (both raw and deconvoluted) ##############################
###################################################################################

marker_genes_raw <- readRDS(file = '_rdata/04_marker_genes_raw.rds')
marker_genes_deconv <- readRDS(file = '_rdata/04_marker_genes_deconv.rds')

marker_genes_list <- 
  lapply(names(marker_genes_raw), 
         function(k) {
           marker_list <- list(raw=marker_genes_raw[[k]], deconv=marker_genes_deconv[[k]])
           marker_intersections <- attr(venn::venn(marker_list), 'intersections')
           do.call('rbind', 
                   lapply(names(marker_intersections), 
                          function(j) {
                            data.frame(
                              gene=marker_intersections[[j]], 
                              type=j
                            )
                          }))         
         })
names(marker_genes_list) <- names(marker_genes_raw)
openxlsx::write.xlsx(marker_genes_list, 
                     file = '_tables/marker_genes_list.xlsx')


###################################################################################
##### visualize marker genes on the whole_tissue vs TRAP expression plot ##########
######## + make BP gene ontology with compareCluster ##############################
###################################################################################

all_vsd <- vst(all_dds)
condition <- apply(data.frame(colData(all_vsd))[,c('tissue','selection')], 1, 
                   function(x) paste(as.character(x), collapse = '__'))
all_vsd_means <- data.frame(sapply(unique(condition), function(cond)
  rowMeans(assay(all_vsd)[,condition == cond])))

for(tissue in names(marker_genes_deconv)) {
  
  print(tissue)

  pdf(paste0('_panel_figures/04_venn_',tissue,'.pdf'), height = 2, width = 2)
  vennOut <-venn::venn(
    list(raw=marker_genes_raw[[tissue]], 
         deconv=marker_genes_deconv[[tissue]]), box = F
    )
  dev.off()
  intersections <- attr(vennOut, 'intersections')
  tissue_vsd_means <- all_vsd_means[,grep(tissue, colnames(all_vsd_means))]
  tissue_vsd_means <- data.frame(tissue_vsd_means)
  colnames(tissue_vsd_means) <- sub(tissue, '', colnames(tissue_vsd_means))
  colnames(tissue_vsd_means) <- sub('__', '', colnames(tissue_vsd_means))
  # add marker info
  tissue_vsd_means$class <- 'no_marker'
  for(cl in names(intersections))
    tissue_vsd_means[intersections[[cl]], 'class'] <- cl
  tissue_vsd_means$class <- factor(tissue_vsd_means$class)
  tissue_vsd_means$class <- relevel(tissue_vsd_means$class, 'no_marker')
  tissue_vsd_means <- tissue_vsd_means[order(tissue_vsd_means$class),]
  # dot plot
  gg <- ggplot(tissue_vsd_means, 
         aes(whole_tissue, TRAP, color=class, size=class)) + 
    geom_point() + theme_bw() + 
    scale_color_manual(values = c(no_marker='grey50', deconv='blue', 
                                  raw='red', 'raw:deconv'='purple')) + 
    scale_size_manual(values = c(no_marker=.5, deconv=2, raw=2, 'raw:deconv'=2))
  ggsave(paste0('_panel_figures/04_dot_',tissue,'.pdf'), gg, height = 5, width = 6)
  
  # compareCluster BP
  intersections_BP <- 
    compareCluster(intersections, fun='enrichGO', 
                   OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', ont = 'BP')
  intersections_BP@compareClusterResult$log10.p.adjust <- 
    log10(intersections_BP@compareClusterResult$p.adjust)
  intersections_BP@compareClusterResult$Cluster <- 
    factor(intersections_BP@compareClusterResult$Cluster, 
           levels = c('deconv','raw:deconv','raw'))
  bp_gg <- dotplot(intersections_BP, color = 'log10.p.adjust')
  ggsave(paste0('_panel_figures/04_BP_',tissue,'_compare.pdf'), bp_gg, 
         height = 5, width = 6.5)
  

}

tissue <- 'SpinalCord'

tissue_vsd_means <- all_vsd_means[,grep(tissue, colnames(all_vsd_means))]
tissue_vsd_means <- data.frame(tissue_vsd_means)
colnames(tissue_vsd_means) <- sub(tissue, '', colnames(tissue_vsd_means))
colnames(tissue_vsd_means) <- sub('__', '', colnames(tissue_vsd_means))
# # add marker info
# tissue_vsd_means$class <- 'no_marker'
# for(cl in names(intersections))
#   tissue_vsd_means[intersections[[cl]], 'class'] <- cl
# tissue_vsd_means$class <- factor(tissue_vsd_means$class)
# tissue_vsd_means$class <- relevel(tissue_vsd_means$class, 'no_marker')
# tissue_vsd_means <- tissue_vsd_means[order(tissue_vsd_means$class),]
# dot plot

orange_genes <- c("Slc16a9","Plxna2","Pard6g","Plekhg1","Akr1b8","Mc5r","Slc24a5","Rgs16")
tissue_vsd_means$label <- NA
tissue_vsd_means[orange_genes, 'label'] <- orange_genes
tissue_vsd_means$color <- 'no'
tissue_vsd_means[orange_genes, 'color'] <- 'orange'
tissue_vsd_means <- tissue_vsd_means %>% arrange(color)
tissue_vsd_means$size <- .1
tissue_vsd_means[orange_genes,'size'] <- 1

ggplot(tissue_vsd_means, 
             aes(whole_tissue, TRAP, label=label, size = size, color = color)) + #, color=class, size=class)) + 
  geom_point() + theme_bw() + geom_label_repel() + 
  geom_abline(intercept = 0, slope = 1)#+ 
#  scale_color_manual(values = c(no_marker='grey50', deconv='blue', 
#                                raw='red', 'raw:deconv'='purple')) + 
#  scale_size_manual(values = c(no_marker=.5, deconv=2, raw=2, 'raw:deconv'=2))

yellow_genes <- c("Dchs1","Bmp6","Lsr","Col4a1","Vwa1","Sema7a","Shank3","Grrp1","Dpp4","Npr1","Notch4","Mmp28","Gipc3","Efnb2","Adamts9","Bmp2","Hoxb7","Slco2a1","Selp")
tissue_vsd_means$label <- NA
tissue_vsd_means[yellow_genes, 'label'] <- yellow_genes
tissue_vsd_means$color <- 'no'
tissue_vsd_means[yellow_genes, 'color'] <- 'yellow'
tissue_vsd_means <- tissue_vsd_means %>% arrange(color)
tissue_vsd_means$size <- 'no'
tissue_vsd_means[yellow_genes,'size'] <- 'yellow'
  
  ggplot(tissue_vsd_means, 
         aes(whole_tissue, TRAP, label=label, size = size, color = color)) + #, color=class, size=class)) + 
    geom_point() + theme_bw() + geom_label_repel() + 
    geom_abline(intercept = 0, slope = 1) + 
    scale_color_manual(values = c(no='grey50', yellow='blue')) + 
    scale_size_manual(values = c(no=.1, yellow=3))
  

###################################################################################
######## make a global comparison between all Cleuren samples (Tek+) and ##########
#################### Bonanomi samples (Cdh5+) #####################################
###################################################################################

# -------- retrieve cell annotations from TabulaMuris

tabula_path <- '~/OneDrive - Ospedale San Raffaele/TabulaMuris/'

annotations_facs <- readRDS(paste0(tabula_path, 'FACS-Sorted-sc/annotations_facs.rds'))
message('spinal cord cells:')
str(spinal_cord_cells <- c("astrocyte","endothelial cell",
                           "fibroblast","oligodendrocyte","macrophage",
                           "neuron","oligodendrocyte precursor cell",
                           "brain pericyte","microglial cell"))

message('cleuren tissues cells:')
str(cleuren_tissues_cells <- setdiff(unique(subset(annotations_facs, tissue %in% c('Brain_Myeloid','Brain_Non-Myeloid', 'Lung', 'Kidney', 'Heart', 'Liver'))$cell_ontology_class), ''))

# -------- retrieve (maxAvg) specificity tables from TabulaMuris

maxavg_expr  <- readRDS(paste0(tabula_path, '/alltissues_cells_maxavgExpr.rds'))
cleuren_expr <- maxavg_expr[,cleuren_tissues_cells]
# calculate specificity
cleuren_spec <- cleuren_expr/rowSums(cleuren_expr)

# hist(log10(apply(cleuren_expr, 1, max)), 50)
# abline(v=log10(5e-2))

# filter for genes expressed over a certain threshold
expressed_genes <- apply(cleuren_expr, 1, max)>5e-2
table(expressed_genes)
cleuren_spec <- cleuren_spec[expressed_genes, ]

cleuren_spec <- cleuren_spec[intersect(rownames(cleuren_spec), rownames(dtrap_dds)),]

# modify cell type names
colnames(cleuren_spec) <- stringr::str_to_title(paste0(colnames(cleuren_spec), 's'))
colnames(cleuren_spec) <- gsub('-| ', '_', colnames(cleuren_spec))
colnames(cleuren_spec)[colnames(cleuren_spec)=='Endothelial_Cells'] <- 'Endothelium' 

# --------- compare SpCord vs ALL

dtrap_markers_dds <- dtrap_dds[,dtrap_dds$selection == 'TRAP']
dtrap_dataset_markers_dds <- DESeqDataSetFromMatrix(counts(dtrap_markers_dds), 
                                                    colData(dtrap_markers_dds), 
                                                    ~ dataset)
dtrap_dataset_markers_vsd <- vst(dtrap_dataset_markers_dds)
rowData(dtrap_dataset_markers_dds) <- rowData(dtrap_markers_dds)

# MFSD2A case

g_d <- plotCounts(dtrap_dataset_markers_dds, gene = 'Mfsd2a', intgroup = 'dataset', returnData = T) %>% 
  mutate(dataset=ifelse(dataset=='Bonanomi', 'SpinalCord','Other')) %>% 
  ggplot(aes(dataset, log10(count))) + geom_boxplot() + theme_bw()

g_t <- plotCounts(dtrap_dataset_markers_dds, gene = 'Mfsd2a', intgroup = 'tissue', returnData = T) %>% 
  ggplot(aes(tissue, log10(count))) + geom_boxplot() + theme_bw()

gg <- ggarrange(g_d, g_t, widths = c(1.5, 4), ncol = 2)
ggsave('_panel_figures/04_Mfsd2a_counts.pdf', gg, height = 3, width = 6)

degs_dataset <- makeStrictComparison('dataset_Bonanomi_vs_Cleuren', 
                                     dtrap_dataset_markers_dds,
                                     fpkms_filter = 3, 
                                     counts_filter = 20, f = any)

degs_dataset <- na.omit(degs_dataset)
degs_dataset$saturated <- degs_dataset$padj < 1e-50
degs_dataset$padj_saturated <- degs_dataset$padj
degs_dataset$padj_saturated[degs_dataset$saturated] <- 1e-50

# in this case any might be required (instead of all) because the Cleuren condition is too heterogeneous to 
# take only samples with ALL samples above a threshold

# add cell type specificities
degs_dataset_specs <- cbind(data.frame(degs_dataset), data.frame(cleuren_spec)[rownames(degs_dataset),])
degs_dataset_specs$max_spec <- apply(data.frame(cleuren_spec)[rownames(degs_dataset),],1,max)
degs_dataset_specs$Fibroblasts_VS_Macrophages <- 
  degs_dataset_specs$Fibroblasts - degs_dataset_specs$Macrophages
degs_dataset_specs$Fibroblasts_VS_Endothelium <- 
  degs_dataset_specs$Fibroblasts - degs_dataset_specs$Endothelium
degs_dataset_specs$Macrophages_VS_Endothelium <- 
  degs_dataset_specs$Macrophages - degs_dataset_specs$Endothelium

gg <- ggplot(degs_dataset_specs[order(abs(degs_dataset_specs$Fibroblasts_VS_Macrophages), na.last = F),], 
             aes_string('log2FoldChange', '-log10(padj_saturated)', 
                        color='Fibroblasts_VS_Macrophages', 
                        shape='saturated')) + 
  geom_point() + 
  # geom_point(data=orderAndSubset(degs_dataset_specs, 'Fibroblasts_VS_Macrophages', abs_val = T)) +
  scale_colour_gradient2(low='darkblue', high = 'darkgreen', limits = c(-1,1),
                         na.value = 'grey90') + 
  theme_bw() + ggtitle('Cdh5+ vs Tek+ EC')

gg

pdf('_panel_figures/04_Tek_vs_Cdh5_specificity.pdf', height = 4.5)
gg
dev.off()

# ---------- view fibro vs endo specificities

cmp <- degs_dataset_specs

gene_labels <- mine_degs(cmp, spec = .5, l2fc = 2) %>% unlist()
cmp$labels <- NA
cmp[gene_labels, 'labels'] <- gene_labels

gg_right <- cmp %>% 
  filter(!is.na(Fibroblasts_VS_Endothelium)) %>% 
  arrange(abs(Fibroblasts_VS_Endothelium)) %>% 
  ggplot(aes_string('log2FoldChange', '-log10(padj_saturated)', 
                    color='Fibroblasts_VS_Endothelium', 
                    shape='saturated', label='labels')) +
  geom_point() + xlim(0,12.5) + ylab('') + 
  scale_colour_gradient2(low='darkred', high = 'darkgreen', limits = c(-1,1), na.value = 'grey90') +
  geom_text_repel() + 
  theme_bw() + theme(legend.position = 'none', 
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank())

gg_left <- cmp %>% 
  filter(!is.na(Macrophages_VS_Endothelium)) %>% 
  arrange(abs(Macrophages_VS_Endothelium)) %>% 
  ggplot(aes_string('log2FoldChange', '-log10(padj_saturated)', 
                    color='Macrophages_VS_Endothelium', 
                    shape='saturated', label='labels')) +
  geom_point() + xlim(-12.5,0) +
  scale_colour_gradient2(low='darkred', high = 'darkblue', limits = c(-1,1), na.value = 'grey90') +
  geom_text_repel() + 
  theme_bw() + theme(legend.position = 'none')

gg <- gg_left + gg_right + plot_layout()

ggsave('_panel_figures/04_SpinalCord_VS_ALL_volcano_wLabels.pdf', gg, height = 5, width = 6)

degs_list_vsALL <- mine_degs(degs_dataset_specs, spec = .2, l2fc = 1)

# ---------- view celltype specificities once at a time

dir.create('_panel_figures/04_Tek_vs_Cdh5_specificities')

N <- nrow(degs_dataset_specs)
n_labels <- 100
for(j in colnames(cleuren_spec)) {
  degs_dataset_ct <- degs_dataset_specs[order(degs_dataset_specs[,j], na.last = F),]
  degs_dataset_ct$label <- NA
  degs_dataset_ct[(N-n_labels+1):N,'label'] <- rownames(degs_dataset_ct)[(N-n_labels+1):N]
  non_significant <- abs(degs_dataset_ct$log2FoldChange)<2 | degs_dataset_ct$padj > 1e-3
  degs_dataset_ct[non_significant,'label'] <- NA
  pdf(paste0('_panel_figures/04_Tek_vs_Cdh5_specificities/',j,'.pdf'), height = 4, width = 5)
  print(ggplot(degs_dataset_ct, 
               aes_string('log2FoldChange', '-log10(padj_saturated)', 
                          color=j, label = 'label',
                          shape='saturated')) +
          geom_point() + 
          geom_text_repel(max.overlaps = Inf) +
          ggtitle(j) +
          scale_colour_gradient2(low='blue', high = 'red',
                                 na.value = 'grey90', limits = c(0,1)) +
          theme_bw())
  dev.off()
}

# ---------- specificity barplot

cdh5_vs_tek_celltypes <- sapply(colnames(cleuren_spec), function(j) {
  score <- degs_dataset_specs[,j] * -log10(degs_dataset_specs$padj) * sign(degs_dataset_specs$log2FoldChange)
  sum(score, na.rm = T)
})

cdh5_vs_tek_df <- data.frame(cell_type=names(cdh5_vs_tek_celltypes), score=cdh5_vs_tek_celltypes)
cdh5_vs_tek_df$cell_type <- factor(cdh5_vs_tek_df$cell_type, levels = names(sort(cdh5_vs_tek_celltypes)))

gg <- ggplot(cdh5_vs_tek_df, aes(score, cell_type, fill=score)) + geom_bar(stat='identity') + 
    scale_fill_gradient2(low='darkblue', high = 'darkgreen') + theme_bw()

pdf('_panel_figures/04_Tek_vs_Cdh5_specificity_barplot.pdf')
gg
dev.off()

# ----------- display GO by Ilaria

tissue_names <- c('brain','heart','kidney','liver','lung')
sat <- c(-5,-10,-3,-5,-3)
names(sat) <- tissue_names

for(tn in tissue_names)
{

  go_list <- list(
    common = 
      cbind(openxlsx::read.xlsx(
        paste0('_GO_ECmarkers/GO_',tn,' common.xlsx'), sheet = 'Summary'), 
            method='common'),
    raw = 
      cbind(openxlsx::read.xlsx(
        paste0('_GO_ECmarkers/GO_',tn,' raw.xlsx'), sheet = 'Summary'),
            method='raw'
            ))
  
  go_df <- do.call('rbind', go_list)
  colnames(go_df) <- make.names(colnames(go_df))
  go_df$log10.p.adjusted <- log10(go_df$Adjusted.P.value)
  go_df$log10.p.adjusted <- scales::oob_squish(go_df$log10.p.adjusted, 
                                               range = c(sat[tn],0))
  go_df$method <- factor(go_df$method, levels = c('raw','common'))
  # go_df$Term <- stringr::str_wrap(go_df$Term, width = 15)
  
  term_order <- reshape2::acast(go_df, Term ~ method, value.var = 'log10.p.adjusted')
  term_levels <- names(sort(term_order[,'common'] - term_order[,'raw']))
  go_df$Term <- factor(go_df$Term, levels = rev(term_levels))
  
  ggplot(go_df, aes(method, Term, 
                    color=log10.p.adjusted, 
                    size=Odds.Ratio)) +
    geom_point() + scale_color_gradient(high='black', low = 'red') + 
    xlab('') + ylab('') +
    scale_y_discrete(labels = scales::label_wrap(30)) +
    theme_bw() -> gg
  
  gg_fixed <- set_panel_size(gg, width  = unit(4.5, "cm"),
                             height = unit(7, "cm"))
  
  ggsave(paste0('_panel_figures/GO_',tn,'_markers.pdf'), gg_fixed,
         height = 5, width = 5)
  
}

# ----------- 

library(viridis)

GO_tissues <- openxlsx::read.xlsx('_GO_ECmarkers/GO tissues_summary.xlsx')
colnames(GO_tissues) <- make.names(colnames(GO_tissues))
GO_tissues <- GO_tissues[!apply(is.na(GO_tissues[,1:4]),1,all),]
GO_tissues$N <- strsplit(GO_tissues$Overlap, '/') %>% sapply('[',1) %>% as.numeric
GO_tissues$tissue <- paste('vs.', GO_tissues$tissue)
GO_tissues <- GO_tissues %>% mutate(
  log10p=ifelse(class=='up', log10(Adjusted.P.value), -log10(Adjusted.P.value)))
GO_tissues$class <- factor(GO_tissues$class, levels = c('up','down'))
GO_tissues$Description <- factor(GO_tissues$Description, 
                                 levels = GO_tissues$Description[order(-as.numeric(GO_tissues$class), GO_tissues$Description)] %>% unique)

# volendo clusterizzare i tessuti e processi
dcast(GO_tissues %>% filter(!str_detect(tissue, 'All')), 
      tissue ~ Description, value.var = 'log10p') %>% 
  column_to_rownames('tissue') %>% pheatmap -> tissueClust
GO_tissues$tissue <- factor(GO_tissues$tissue, 
                            c('vs. All', tissueClust$tree_row$labels[tissueClust$tree_row$order]))
GO_tissues$Description <- factor(GO_tissues$Description, 
                                 tissueClust$tree_col$labels[tissueClust$tree_col$order])
GO_tissues$log10p[GO_tissues$log10p < (-12)] <- -12

GO_tissues %>% filter(class=='up') %>% 
  ggplot(aes(Description, tissue, size=N, color=log10p)) + 
    geom_point() + theme_bw() + 
    # scale_color_gradient2(low='red', mid = 'black', high = 'blue') + 
    scale_colour_gradientn(limits = c(-12, 0), 
                           colors = RColorBrewer::brewer.pal(10,'Spectral')) + 
                           #colours = rev(c('white','pink1','red','darkred'))) + 
    theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1)) +
    facet_wrap(~class, scales = 'free_x') + xlab('') + ylab('') -> gg

gg_fixed <- set_panel_size(gg, width  = unit(4.5, "cm"), height = unit(7, "cm"))

ggsave('_panel_figures/04_Tissue_markers_up_spectral.pdf', gg_fixed, height = 7, width = 5)

GO_tissues %>% filter(class=='down') %>% 
  ggplot(aes(Description, tissue, size=N, color=log10p)) + 
  geom_point() + theme_bw() + scale_color_gradient2(limits = c(0,2.5), low='red', mid = 'white', high = 'blue') + 
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1)) +
  facet_wrap(~class, scales = 'free_x') + xlab('') + ylab('')  -> gg

gg_fixed <- set_panel_size(gg, width  = unit(2.5, "cm"), height = unit(7, "cm"))

ggsave('_panel_figures/04_Tissue_markers_down.pdf', gg_fixed, height = 7, width = 5)

# ------------ compare SpinalCord against each tissue 

cleuren_tissues <- setdiff(levels(dtrap_dataset_markers_dds$tissue), 'SpinalCord')
names(cleuren_tissues) <- cleuren_tissues

spCord_vs <- lapply(cleuren_tissues, function(tn) {
  degs_dataset <- makeStrictComparison(paste0('tissue_SpinalCord_vs_',tn), 
                                       dtrap_dataset_markers_dds,
                                       fpkms_filter = 3, 
                                       counts_filter = 20, f = any)
  
  degs_dataset <- na.omit(degs_dataset)
  degs_dataset$saturated <- degs_dataset$padj < 1e-100
  degs_dataset$padj_saturated <- degs_dataset$padj
  degs_dataset$padj_saturated[degs_dataset$saturated] <- 1e-100
  
  degs_dataset_specs <- cbind(data.frame(degs_dataset), data.frame(cleuren_spec)[rownames(degs_dataset),])
  degs_dataset_specs$max_spec <- apply(data.frame(cleuren_spec)[rownames(degs_dataset),],1,max)
  degs_dataset_specs$Macrophages_VS_Endothelium <- 
    degs_dataset_specs$Macrophages - degs_dataset_specs$Endothelium
  degs_dataset_specs$Fibroblasts_VS_Endothelium <- 
    degs_dataset_specs$Fibroblasts - degs_dataset_specs$Endothelium  
  return(degs_dataset_specs)
})

Bioplanet_pathways <- read.csv('~/OneDrive - Ospedale San Raffaele/COSR/annot/Bioplanet_pathways.csv')
TERM2GENE <- Bioplanet_pathways %>% mutate(GENE_SYMBOL=stringr::str_to_title(GENE_SYMBOL)) %>% dplyr::select(PATHWAY_NAME, GENE_SYMBOL)

for(k in cleuren_tissues) {
  mine_degs(spCord_vs[[k]], spec = .2, l2fc = 1) %>% #str()
    compareCluster(fun = 'enricher', TERM2GENE = TERM2GENE) %>% dotplot() -> gg
  gg <- gg + ggtitle(paste('SpinalCord VS', k))
  ggsave(paste0('_panel_figures/04_SpinalCord_VS_',k,'.pdf'), gg, height = 5, width = 8)
}

# --- write DEGs on excel

openxlsx::write.xlsx(mine_degs_df(spCord_vs$Brain, spec = .2, l2fc = 1), 
                     rowNames = T, overwrite = T,
                     file = '_tables/04_SpinalCord_VS_Brain.xlsx')

openxlsx::write.xlsx(mine_degs_df(degs_dataset_specs, spec = .2, l2fc = 1), 
                     rowNames = T, overwrite = T,
                     file = '_tables/04_SpinalCord_VS_ALL.xlsx')

# -------- old DEGs on excel

degs_list_vsTissue <- lapply(cleuren_tissues, function(k) {
  mine_degs(spCord_vs[[k]], spec = .2, l2fc = 1)
})
names(degs_list_vsTissue) <- paste0('vs_', names(degs_list_vsTissue))
endo_degs_list <- c(vs_ALL=list(degs_list_vsALL), degs_list_vsTissue)

lapply(endo_degs_list, function(degs_list) {
  N <- sapply(degs_list, length) %>% max()
  lapply(degs_list, function(x) {
    if(length(x)<N) x <- c(x, rep('', N-length(x)))
    return(x)
  }) %>% data.frame() -> degs_df
}) -> endo_degs
#   
openxlsx::write.xlsx(endo_degs, file = '_tables/04_SpinalCord_VS_each_geneList.xlsx')

# ----- common DEGs

# lapply(endo_degs, '[[', 'up') %>% venn::venn()
# lapply(endo_degs, '[[', 'down') %>% venn::venn()

venn_up <- lapply(endo_degs, '[[', 'up') %>% venn::venn() %>% attr('intersections')
venn_down <- lapply(endo_degs, '[[', 'down') %>% venn::venn() %>% attr('intersections')

list(
  always_down = venn_down$`Brain:Heart:Kidney:Liver:Lung`,
  always_up = venn_up$`Brain:Heart:Kidney:Liver:Lung`) -> recurrent_degs
  
recurrent_enrich <- compareCluster(fun = 'enricher', TERM2GENE = TERM2GENE)

dotplot(recurrent_enrich) 

ggsave('_panel_figures/SpianlCord_recurrent_degs.pdf', height = 5, width = 7)

for(k in cleuren_tissues) {
  
  cmp <- spCord_vs[[k]]
  
  gene_labels <- mine_degs(spCord_vs[[k]], spec = .5, l2fc = 4) %>% unlist()
  cmp$labels <- NA
  cmp[gene_labels, 'labels'] <- gene_labels
  
  gg_right <- cmp %>% 
    filter(!is.na(Fibroblasts_VS_Endothelium)) %>% 
    arrange(abs(Fibroblasts_VS_Endothelium)) %>% 
    ggplot(aes_string('log2FoldChange', '-log10(padj_saturated)', 
                      color='Fibroblasts_VS_Endothelium', 
                      shape='saturated', label='labels')) +
    geom_point() + xlim(0,20) + ylab('') + 
    scale_colour_gradient2(low='darkred', high = 'darkgreen', limits = c(-1,1), na.value = 'grey90') +
    geom_text_repel() + 
    theme_bw() + theme(legend.position = 'none', 
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank())
  
  gg_left <- cmp %>% 
    filter(!is.na(Macrophages_VS_Endothelium)) %>% 
    arrange(abs(Macrophages_VS_Endothelium)) %>% 
    ggplot(aes_string('log2FoldChange', '-log10(padj_saturated)', 
                      color='Macrophages_VS_Endothelium', 
                      shape='saturated', label='labels')) +
    geom_point() + xlim(-20,0) +
    scale_colour_gradient2(low='darkred', high = 'darkblue', limits = c(-1,1), na.value = 'grey90') +
    geom_text_repel() + 
    theme_bw() + theme(legend.position = 'none')
  
  gg <- gg_left + gg_right + plot_layout()
  
  ggsave(paste0('_panel_figures/04_SpinalCord_VS_',k,'_volcano.pdf'), 
         gg, height = 5, width = 6)
  
}


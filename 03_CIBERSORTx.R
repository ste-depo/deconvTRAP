library(ggplot2)
library(patchwork)
library(reshape2)
library(DESeq2)
library(dplyr)

source('../common/00_generic_functions.R')

## -------------------------------------------------------------------------
## ----- create input files for CIBERSORTx ---------------------------------
## -------------------------------------------------------------------------

all_dds <- readRDS(file = '_rdata/01_CleurenBonanomi_raw.rds')
dtrap_dds <- readRDS(file = '_rdata/02_CleurenBonanomi_deconv.rds')

# raw dataset

cibersortx_in <- counts(all_dds)
cibersortx_in <- as.matrix(cibersortx_in)
cibersortx_in <- as.data.frame(cibersortx_in)
cibersortx_in <- cbind(GeneSymbol=rownames(cibersortx_in), cibersortx_in)
write.table(cibersortx_in, file='_cibersortx_in/03_EC_TRAP_raw_counts.txt', quote=F, sep='\t', row.names = F)

# deconvoluted dataset

cibersortx_in <- counts(dtrap_dds)
cibersortx_in <- as.matrix(cibersortx_in)
cibersortx_in <- as.data.frame(cibersortx_in)
cibersortx_in <- cbind(GeneSymbol=rownames(cibersortx_in), cibersortx_in)
write.table(cibersortx_in, file='_cibersortx_in/03_EC_TRAP_deconv_counts.txt', quote=F, sep='\t', row.names = F)

# per tissue raw + deconv

for(k in unique(dtrap_dds$tissue)) {
  tissue_dds_raw <- counts(all_dds[,all_dds$tissue == k])
  colnames(tissue_dds_raw) <- paste0(colnames(tissue_dds_raw),'_r')
  tissue_dds_deconv <- counts(dtrap_dds[,dtrap_dds$tissue == k])
  colnames(tissue_dds_deconv) <- paste0(colnames(tissue_dds_deconv),'_d')
  # -------
  cibersortx_in <- cbind(tissue_dds_raw, tissue_dds_deconv)
  cibersortx_in <- as.matrix(cibersortx_in)
  cibersortx_in <- as.data.frame(cibersortx_in)
  cibersortx_in <- cbind(GeneSymbol=rownames(cibersortx_in), cibersortx_in)
  write.table(cibersortx_in, file=paste0('_cibersortx_in/03_',k,'_counts.txt'),
              quote=F, sep='\t', row.names = F)
}

## -------------------------------------------------------------------------
## ----- create plots from the output of CIBERSORTx ------------------------
## -------------------------------------------------------------------------

all_dds <- readRDS(file = '_rdata/01_CleurenBonanomi_raw.rds')

samples_ann <- data.frame(colData(all_dds))
samples_ann$Mixture <- rownames(samples_ann)

# file_in <- list.files('_cibersortx_out', full.names = T)[1]

cibx_plot <- function(file_in, samples_ann, what, average, col=NA,
                         filename = NA, width = 7, height = 5) {
  
  cib_in <- read.csv(file_in, row.names = 1)
  cib_res <- cib_in[,1:(ncol(cib_in)-3)]
  cib_res$Mixture <- sapply(strsplit(rownames(cib_res), '_'), 
                            function(x) 
                              paste(x[1:(length(x)-1)], collapse = '_'))
  cib_res$Method <- sapply(strsplit(rownames(cib_res), '_'), 
                           function(x) 
                             x[length(x)])
  cib_res <- merge(cib_res, samples_ann[,c('selection','Mixture')], by='Mixture')
  
  if(average) {
    
    cib_res <- cib_res[cib_res$selection == what,]
    cib_samples_ann <- cib_res[,(ncol(cib_res)-1):ncol(cib_res)]
    cib_res_mean <- apply(cib_res[,2:(ncol(cib_res)-2)], 2, 
                          function(x) 
                            tapply(x, cib_samples_ann$Method, mean))
    cib_res_mean <- data.frame(cib_res_mean)
    cib_res_mean$Method <- factor(rownames(cib_res_mean), levels = c('r', 'd'))
    levels(cib_res_mean$Method) <- c('raw', 'deconv')
    cib_res_long <- melt(cib_res_mean, id.vars = c('Method'))
    gg <- ggplot(cib_res_long, aes(Method, value, fill=variable)) + 
      geom_bar(stat = 'identity') + 
      theme_bw()
    if(!is.na(col[1])) gg <- gg + scale_fill_manual(values = col)
    
  } else {
    
    cib_ann <- cib_in[,(ncol(cib_in)-2):ncol(cib_in)]
    cib_res_long <- melt(cib_res, id.vars = c('Mixture','Method','selection'))
    gr <- ggplot(subset(cib_res_long, Method=='r' & selection == what), 
                 aes(Mixture, value, fill=variable)) + geom_bar(stat = 'identity') +
      theme_bw()
    if(!is.na(col[1])) gr <- gr + scale_fill_manual(values = col)
    gd <- ggplot(subset(cib_res_long, Method=='d' & selection == what), 
                 aes(Mixture, value, fill=variable)) + geom_bar(stat = 'identity') +
      theme_bw()
    if(!is.na(col[1])) gd <- gd + scale_fill_manual(values = col)
    gg <- gr + gd + plot_layout(guides = 'collect')  
    
  }
  
  if(is.na(filename)) print(gg) else {
    ggsave(filename, plot = gg, width=width, height = height)
  }

}

# set.seed(0)
# cns_cols <- sample(depo_palette, 14)
# names(cns_cols) <- c("Endothelium","OPC","Early.Oligodendrocytes",
#                      "OPC.Proliferating","Neurons","Oligodendrocytes","Fibroblasts",
#                      "Interneurons","pyramidal.SS","Astrocytes","pyramidal.CA1",
#                      "Microglia","SMC","Pericytes")

cns_cols <- c(
"Endothelium"=  '#FF0010', #Red  
"OPC"=  '#191919', #Ebony
"Early.Oligodendrocytes"=    '#FFA8BB', #Pink
"OPC.Proliferating"='#4C005C', #Damson
"Neurons"=  '#FFFF80', #Xanthin
"Oligodendrocytes"='#C20088', #Mallow
"Fibroblasts"=  '#740AFF', #Violet
"Interneurons"=  '#808080', #Iron
"pyramidal.SS"=  '#E0FF66', #Uranium
"Astrocytes"='#FF5005', #Zinnia
"pyramidal.CA1"=  '#FFFF00', #yellow
"Microglia"=    '#FFCC99', #Honeydew
"SMC"=  '#5EF1F2', #Sky
"Pericytes"=    '#F0A3FF' #Amethyst
)

for(tissue in c('Brain','Heart','Kidney','Lung','Liver','SpinalCord')) {
  file_in <- paste0('_cibersortx_out/CIBERSORTx_', tissue, '_Results.csv')
  cibx_plot(file_in, samples_ann,
            what = 'TRAP', average = FALSE, 
            col = if(tissue %in% c('Brain', 'SpinalCord')) cns_cols else NA,
            filename = paste0('_panel_figures/03_CIBERSORTx_',tissue,'_TRAP.pdf'))
  cibx_plot(file_in, samples_ann,
            what = 'TRAP', average = TRUE, 
            col = if(tissue %in% c('Brain', 'SpinalCord')) cns_cols else NA,
            filename = paste0('_panel_figures/03_CIBERSORTx_',tissue,'_TRAP_avg.pdf'))
  cibx_plot(file_in, samples_ann,
            what = 'whole_tissue', average = TRUE, width = 5,
            col = if(tissue %in% c('Brain', 'SpinalCord')) cns_cols else NA,
            filename = paste0('_panel_figures/03_CIBERSORTx_',tissue,'_wt.pdf'))
}



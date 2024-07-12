library(DESeq2)
library(pheatmap)
library(reshape2)
library(tidyverse)

source('../common/00_generic_functions.R')
source('../common/00_trap_modeling_functions.R')
all_dds <- readRDS(file = '_rdata/01_CleurenBonanomi_raw.rds')

# deconvolve TRAP

(dtypes <- unique(data.frame(colData(all_dds))[,c('tissue','dataset')]))

dcounts <- matrix(0, nrow = nrow(all_dds), ncol = ncol(all_dds), dimnames = dimnames(all_dds))
alphas <- betas <- matrix(NA, nrow=nrow(dtypes), 3, dimnames = list(as.character(dtypes$tissue)))
n_left_peaks <- matrix(c(
  2, 1, 1,
  1, 1, 1,
  1, 1, 1,
  1, 1, 1,
  1, 1, 1,
  2, 1, 1
), byrow = T, ncol = 3)
n_right_peaks <- matrix(c(
  2, 1, 2,
  1, 1, 1,
  1, 2, 1,
  1, 1, 1,
  1, 1, 1,
  1, 1, 1
), byrow = T, ncol = 3)


## -------------------------------------------------------- ##
## -------------- perform the deconvolution --------------- ##
## -------------------------------------------------------- ##

res <- 150
png('_panel_figures/02_deconvMethod.png', height = 60*res, width = 10*res, res=res)

par(mfrow=c(nrow(dtypes)*3,3), mar=c(2,2,4,1))
for(i in 1:nrow(dtypes)) {
  idxTRAP <- which(all_dds$tissue == dtypes[i,'tissue'] &
                     all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'TRAP')
  idxPare <- which(all_dds$tissue == dtypes[i,'tissue'] &
                     all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'whole_tissue')
  deconv <- deconvolveTRAPddsNew(all_dds, idxTRAP, idxPare, minVal = 25, 
                                 n_left_peak = n_left_peaks[i,], n_right_peak = n_right_peaks[i,],
                                 method = 'peak', plot = TRUE, 
                                 min_obs = if(i %in% c(3,6)) 25 else 100, 
                                 alphas = if(i==6) c(1.35,1,1) else NULL,
                                 paired = T, use.cv = FALSE
  )
  dcounts[,idxTRAP] <- deconv$TRAP
  dcounts[,idxPare] <- deconv$lysate
  alphas[i,] <- deconv$alphas
  betas[i,]  <- deconv$betas
}
dev.off()

## store the deconvolution results
save(dcounts, alphas, betas, file='_rdata/02_deconv_results.rda')

## store the DESeq2 dataset for the deconvoluted datasets
dcounts[!is.finite(dcounts)] <- 0
dtrap_dds <- DESeqDataSetFromMatrix(dcounts, 
                                    colData(all_dds)[colnames(dcounts),], ~ tissue + selection)
dtrap_dds <- estimateSizeFactors(dtrap_dds)
rowData(dtrap_dds) <- rowData(all_dds)
dtrap_vsd <- vst(dtrap_dds)
saveRDS(dtrap_dds, file = '_rdata/02_CleurenBonanomi_deconv.rds')

# ------------- reload

load(# dcounts, alphas, betas, 
  file='_rdata/02_deconv_results.rda')
dtrap_dds <- readRDS(file = '_rdata/02_CleurenBonanomi_deconv.rds')
dtrap_vsd <- vst(dtrap_dds)

## ----------------------------------------------------------------- ##
## ------------------- visualize the deconvolution results --------- ##
## ---------------------------- on raw counts ---------------------- ##
## ----------------------------------------------------------------- ##

res <- 150
png('_panel_figures/02_deconvMethod_1.png', height = 20*res, width = 10*res, res=res)
cnt <- 0
dotpl <- list()
par(mfrow=c(nrow(dtypes),3), mar=c(2,2,2,1))
for(i in 1:nrow(dtypes)) {
  idxTRAP <- which(all_dds$tissue == dtypes[i,'tissue'] &
                     all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'TRAP')
  idxPare <- which(all_dds$tissue == dtypes[i,'tissue'] &
                     all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'whole_tissue')
  for(j in 1:3) {
    cnt <- cnt + 1
    idx <- counts(all_dds)[,idxPare[j]] > 5 & counts(all_dds)[,idxTRAP[j]] > 5
    x <- log2((1+counts(all_dds, normalized=T)[idx,idxPare[j]]))
    y <- log2((1+counts(all_dds, normalized=T)[idx,idxTRAP[j]]))
    cols <- densCols(x,y,colramp = colorRampPalette(viridis::inferno(5)))
    specific_genes_1 <- (y-x) > log2(exp(alphas[i,j]))
    specific_genes_2 <- (x-y) > log2(exp(betas[i,j]))
    cols[specific_genes_1|specific_genes_2] <- '#09954F'
      plot(x, y, col=cols,
           cex=.1+.8*grepl('^mt-', rownames(all_dds)[idx]),
           pch=19)
      title(paste0(colnames(all_dds)[idxTRAP[j]], ' vs ', 
                   colnames(all_dds)[idxPare[j]], ' (', 
                   dtypes[i,'tissue'] ,')'))
      if(!is.na(alphas[i,j])) abline(log2(exp(alphas[i,j])), 1, lty=2, lwd=2)
      if(!is.na(betas[i,j])) abline(-log2(exp(betas[i,j])), 1, lty=2, lwd=2)
  }
}
dev.off()


## ----------------------------------------------------------------- ##
## ------------------- visualize the deconvolution results --------- ##
## ---------------------------- on deconv counts ------------------- ##
## ----------------------------------------------------------------- ##

res <- 150

png('_panel_figures/02_deconvMethod_2.png', height = 20*res, width = 10*res, res=res)
cnt <- 0
dotpl <- list()
par(mfrow=c(nrow(dtypes),3), mar=c(2,2,2,1))
for(i in 1:nrow(dtypes)) {
  idxTRAP <- which(dtrap_dds$tissue == dtypes[i,'tissue'] &
                     dtrap_dds$dataset == dtypes[i,'dataset'] & dtrap_dds$selection == 'TRAP')
  idxPare <- which(dtrap_dds$tissue == dtypes[i,'tissue'] &
                     dtrap_dds$dataset == dtypes[i,'dataset'] & dtrap_dds$selection == 'whole_tissue')
  for(j in 1:3) {
    cnt <- cnt + 1
    idx <- counts(dtrap_dds)[,idxPare[j]] > 5 | counts(dtrap_dds)[,idxTRAP[j]] > 5
    x <- log2((1+counts(dtrap_dds, normalized=T)[idx,idxPare[j]]))
    y <- log2((1+counts(dtrap_dds, normalized=T)[idx,idxTRAP[j]]))
    cols <- densCols(x,y,colramp = colorRampPalette(viridis::inferno(5)))
    x_raw <- log2((1+counts(all_dds, normalized=T)[idx,idxPare[j]]))
    y_raw <- log2((1+counts(all_dds, normalized=T)[idx,idxTRAP[j]]))
    specific_genes_1 <- (y_raw-x_raw) > log2(exp(alphas[i,j]))
    specific_genes_2 <- (x_raw-y_raw) > log2(exp(betas[i,j]))
    cols[specific_genes_1|specific_genes_2] <- '#09954F'
      plot(x, y, col=cols,
           cex=.1+.8*grepl('^mt-', rownames(dtrap_dds)[idx]),
           pch=19)
      title(paste0(colnames(dtrap_dds)[idxTRAP[j]], ' vs ', 
                   colnames(dtrap_dds)[idxPare[j]], ' (', 
                   dtypes[i,'tissue'] ,')'))
      if(!is.na(alphas[i,j])) abline(log2(exp(alphas[i,j])), 1, lty=2, lwd=2)
      if(!is.na(betas[i,j])) abline(-log2(exp(betas[i,j])), 1, lty=2, lwd=2)
  }
}
dev.off()


## -----------------------------------------------------------
## ------------------ focus on BRAIN -------------------------
## -----------------------------------------------------------

i <- 2; j <- 1

# -------------------- raw ------------------------------- #

pdf('_panel_figures/02_deconvMethod_BrainFOCUS_1.pdf', height = 5, width=5)

idxTRAP <- which(all_dds$tissue == dtypes[i,'tissue'] &
                   all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'TRAP')
idxPare <- which(all_dds$tissue == dtypes[i,'tissue'] &
                   all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'whole_tissue')

filtered_brain <- counts(all_dds)[,idxPare[j]] > 5 & counts(all_dds)[,idxTRAP[j]] > 5
x <- log2((1+counts(all_dds, normalized=T)[filtered_brain,idxPare[j]]))
y <- log2((1+counts(all_dds, normalized=T)[filtered_brain,idxTRAP[j]]))
cols <- densCols(x,y,colramp = colorRampPalette(viridis::inferno(5)))
specific_genes_trap <- (y-x) > log2(exp(alphas[i,j]))
specific_genes_brain <- (x-y) > log2(exp(betas[i,j]))
cols[specific_genes_trap|specific_genes_brain] <- '#09954F'
plot(x, y, col=cols,
     cex=.1+.8*grepl('^mt-', rownames(all_dds)[filtered_brain]),
     xlab = paste0('Parenchyma (', dtypes[i,'tissue'], ')'),
     ylab = paste0('EC-TRAP (', dtypes[i,'tissue'], ')'),
     pch=19)
if(!is.na(alphas[i,j])) abline(log2(exp(alphas[i,j])), 1, lty=2, lwd=2)
if(!is.na(betas[i,j])) abline(-log2(exp(betas[i,j])), 1, lty=2, lwd=2)

legend('topleft',legend = c('genomically encoded','mt encoded'),
       pch=c(20), col=c('black'), pt.cex = c(.2,.8), bty = 'n')


dev.off()

pdf('_panel_figures/02_deconvMethod_BrainFOCUS_1_Snap25.pdf', height = 5, width=5)

# top10_contaminants <- 
#   head(sort(y[specific_genes_brain & !grepl('^mt-', names(x))], decreasing = T), 10) %>% names()

top_contaminants <- c("Snap25", "Plp1", "Mbp", "Map1b")
idx_gene <- names(x) %in% top_contaminants
cols[idx_gene] <- 'darkgoldenrod3'
plot(x, y, col=cols,
     cex=.1+.8*(idx_gene | grepl('^mt-', names(x))),
     xlab = paste0('Parenchyma (', dtypes[i,'tissue'], ')'),
     ylab = paste0('EC-TRAP (', dtypes[i,'tissue'], ')'),
     pch=19) #c(19,17)[1+1*grepl('^mt-', names(x))])
if(!is.na(alphas[i,j])) abline(log2(exp(alphas[i,j])), 1, lty=2, lwd=2)
if(!is.na(betas[i,j])) abline(-log2(exp(betas[i,j])), 1, lty=2, lwd=2)
# for(idx in which(idx_gene)) {
#   segments(x[idx], 0, x[idx], y[idx], lty=3)
#   segments(0, y[idx], x[idx], y[idx], lty=3)  
# }
text(x[idx_gene], y[idx_gene], labels = top_contaminants, col = 'darkgoldenrod3', 
     pos = c(3,3,3,4), cex = .7)
legend('topleft',legend = c('genomically encoded','top contaminants','chrM encoded'),
       pch=c(20), col=c('black','darkgoldenrod3','black'), pt.cex = c(.2,.8,.8), bty = 'n')

dev.off()

# ----------------- deconvoluted ------------------------- #  

pdf('_panel_figures/02_deconvMethod_BrainFOCUS_2.pdf', height = 5, width=5)

idxTRAP <- which(all_dds$tissue == dtypes[i,'tissue'] &
                   all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'TRAP')
idxPare <- which(all_dds$tissue == dtypes[i,'tissue'] &
                   all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'whole_tissue')


x <- assay(dtrap_vsd)[filtered_brain,idxPare[j]]
y <- assay(dtrap_vsd)[filtered_brain,idxTRAP[j]]
cols <- densCols(x,y,colramp = colorRampPalette(viridis::inferno(5)))
cols[specific_genes_trap|specific_genes_brain] <- '#09954F'
plot(x, y, col=cols,
     cex=.1+.8*grepl('^mt-', rownames(all_dds)[filtered_brain]),
     xlab = paste0('Parenchyma (', dtypes[i,'tissue'], ')'),
     ylab = paste0('EC-TRAP (', dtypes[i,'tissue'], ')'),
     pch=19)
if(!is.na(alphas[i,j])) abline(log2(exp(alphas[i,j])), 1, lty=2, lwd=2)
if(!is.na(betas[i,j])) abline(-log2(exp(betas[i,j])), 1, lty=2, lwd=2)

legend('topleft',legend = c('genomically encoded','mt encoded'),
       pch=c(20), col=c('black'), pt.cex = c(.2,.8), bty = 'n')


dev.off()

pdf('_panel_figures/02_deconvMethod_BrainFOCUS_2_Snap25.pdf', height = 5, width=5)

cols[idx_gene] <- 'darkgoldenrod3'
plot(x, y, col=cols,
     cex=.1+.8*(idx_gene | grepl('^mt-', names(x))),
     xlab = paste0('Parenchyma (', dtypes[i,'tissue'], ')'),
     ylab = paste0('EC-TRAP (', dtypes[i,'tissue'], ')'),
     pch=19)
points(x[idx_gene], y[idx_gene], cex=.8, col='darkgoldenrod3', pch=19)
text(x[idx_gene], y[idx_gene], labels = top_contaminants, col = 'darkgoldenrod3', 
     srt=45, adj = c(0,-.5), cex = .7)

dev.off()


## -----------------------------------------------------------
## ------------------ focus on SpCord ------------------------
## -----------------------------------------------------------

i <- 1; j <- 1

# -------------------- raw ------------------------------- #

pdf('_panel_figures/02_deconvMethod_SpCordFOCUS_1.pdf', height = 5, width=5)

idxTRAP <- which(all_dds$tissue == dtypes[i,'tissue'] &
                   all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'TRAP')
idxPare <- which(all_dds$tissue == dtypes[i,'tissue'] &
                   all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'whole_tissue')

filtered_brain <- counts(all_dds)[,idxPare[j]] > 5 & counts(all_dds)[,idxTRAP[j]] > 5
x <- log2((1+counts(all_dds, normalized=T)[filtered_brain,idxPare[j]]))
y <- log2((1+counts(all_dds, normalized=T)[filtered_brain,idxTRAP[j]]))
cols <- densCols(x,y,colramp = colorRampPalette(viridis::inferno(5)))
specific_genes_trap <- (y-x) > log2(exp(alphas[i,j]))
specific_genes_brain <- (x-y) > log2(exp(betas[i,j]))
cols[specific_genes_trap|specific_genes_brain] <- '#09954F'
plot(x, y, col=cols,
     cex=.1+.8*grepl('^mt-', rownames(all_dds)[filtered_brain]),
     xlab = paste0('Parenchyma (', dtypes[i,'tissue'], ')'),
     ylab = paste0('EC-TRAP (', dtypes[i,'tissue'], ')'),
     pch=19)
if(!is.na(alphas[i,j])) abline(log2(exp(alphas[i,j])), 1, lty=2, lwd=2)
if(!is.na(betas[i,j])) abline(-log2(exp(betas[i,j])), 1, lty=2, lwd=2)

legend('topleft',legend = c('genomically encoded','mt encoded'),
       pch=c(20), col=c('black'), pt.cex = c(.2,.8), bty = 'n')
  
  
dev.off()
  
pdf('_panel_figures/02_deconvMethod_SpCordFOCUS_1_Snap25.pdf', height = 5, width=5)

top10_contaminants <-
  head(sort(y[specific_genes_brain & !grepl('^mt-', names(x))], decreasing = T), 10) %>% names()

top_contaminants <- c("Snap25", "Plp1", "Mbp", "Map1a")
idx_gene <- names(x) %in% top_contaminants
cols[idx_gene] <- 'darkgoldenrod3'
plot(x, y, col=cols,
     cex=.1+.8*(idx_gene | grepl('^mt-', names(x))),
     xlab = paste0('Parenchyma (', dtypes[i,'tissue'], ')'),
     ylab = paste0('EC-TRAP (', dtypes[i,'tissue'], ')'),
     pch=19) #c(19,17)[1+1*grepl('^mt-', names(x))])
if(!is.na(alphas[i,j])) abline(log2(exp(alphas[i,j])), 1, lty=2, lwd=2)
if(!is.na(betas[i,j])) abline(-log2(exp(betas[i,j])), 1, lty=2, lwd=2)
# for(idx in which(idx_gene)) {
#   segments(x[idx], 0, x[idx], y[idx], lty=3)
#   segments(0, y[idx], x[idx], y[idx], lty=3)  
# }
text(x[idx_gene], y[idx_gene], labels = top_contaminants, col = 'darkgoldenrod3', 
     pos = c(3,3,3,4), cex = .7)
legend('topleft',legend = c('genomically encoded','top contaminants','chrM encoded'),
       pch=c(20), col=c('black','darkgoldenrod3','black'), pt.cex = c(.2,.8,.8), bty = 'n')

dev.off()
  
# ----------------- deconvoluted ------------------------- #  

pdf('_panel_figures/02_deconvMethod_SpCordFOCUS_2.pdf', height = 5, width=5)

idxTRAP <- which(all_dds$tissue == dtypes[i,'tissue'] &
                   all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'TRAP')
idxPare <- which(all_dds$tissue == dtypes[i,'tissue'] &
                   all_dds$dataset == dtypes[i,'dataset'] & all_dds$selection == 'whole_tissue')


x <- assay(dtrap_vsd)[filtered_brain,idxPare[j]]
y <- assay(dtrap_vsd)[filtered_brain,idxTRAP[j]]
cols <- densCols(x,y,colramp = colorRampPalette(viridis::inferno(5)))
cols[specific_genes_trap|specific_genes_brain] <- '#09954F'
plot(x, y, col=cols,
     cex=.1+.8*grepl('^mt-', rownames(all_dds)[filtered_brain]),
     xlab = paste0('Parenchyma (', dtypes[i,'tissue'], ')'),
     ylab = paste0('EC-TRAP (', dtypes[i,'tissue'], ')'),
     pch=19)
if(!is.na(alphas[i,j])) abline(log2(exp(alphas[i,j])), 1, lty=2, lwd=2)
if(!is.na(betas[i,j])) abline(-log2(exp(betas[i,j])), 1, lty=2, lwd=2)

legend('topleft',legend = c('genomically encoded','mt encoded'),
       pch=c(20), col=c('black'), pt.cex = c(.2,.8), bty = 'n')


dev.off()

pdf('_panel_figures/02_deconvMethod_SpCordFOCUS_2_Snap25.pdf', height = 5, width=5)

cols[idx_gene] <- 'darkgoldenrod3'
  plot(x, y, col=cols,
       cex=.1+.8*(idx_gene | grepl('^mt-', names(x))),
       xlab = paste0('Parenchyma (', dtypes[i,'tissue'], ')'),
       ylab = paste0('EC-TRAP (', dtypes[i,'tissue'], ')'),
       pch=19)
points(x[idx_gene], y[idx_gene], cex=.8, col='darkgoldenrod3', pch=19)
text(x[idx_gene], y[idx_gene], labels = top_contaminants, col = 'darkgoldenrod3', 
     srt=45, adj = c(0,-.5), cex = .7)

dev.off()

## ------------------------------------------------------------------------
## ------------- show estimated fractions of cross-contamination ----------
## ------------------------------------------------------------------------


pdf('_panel_figures/02_deconvMethod_ECFraction_1.pdf', height = 5, width = 5)
ggplot(melt(alphas), aes(x=Var1, y=exp(-value))) + 
  geom_jitter() + geom_boxplot() + ggtitle('EC fraction in parenchyma') + 
  theme_bw() + ylab('percentage')
dev.off()

pdf('_panel_figures/02_deconvMethod_parenchymaFraction.pdf', height = 5, width = 5)
ggplot(melt(betas), aes(x=Var1, y=exp(-value))) + 
  geom_jitter() + geom_boxplot(alpha=.1) + ggtitle('parenchyma fraction in EC') + 
  theme_bw() + ylab('percentage')
dev.off()

pdf('_panel_figures/02_deconvMethod_ECFraction_2.pdf',
    height = 3, width = 5)
alpha_means <- apply(exp(-alphas[-1,]),1,mean,na.rm=T)*100
alpha_sd <- apply(exp(-alphas[-1,]),1,sd,na.rm=T)*100
beta_means <- apply(exp(-betas[-1,]),1,mean,na.rm=T)*100
beta_sd <- apply(exp(-betas[-1,]),1,sd,na.rm=T)*100
tek_pos_mean <- c(Brain=5.2, Heart=38.5, Kidney=16, Liver=13.2, Lung=48)
tek_pos_sd <- c(Brain=1, Heart=3.2, Kidney=2.2, Liver=1.5, Lung=6)
par(mar=c(2,4,.1,1))
plot(seq(tek_pos_mean)+.1, tek_pos_mean, xaxt='n',ylab='% Tek positive cells',
     xlab='', lty=2, log='', lwd=2, col='grey',
     xlim=c(.9,length(tek_pos_mean)+.1),
     ylim=range(c(tek_pos_mean+tek_pos_sd, 
                  tek_pos_mean-tek_pos_sd,
                  na.omit(alpha_means-alpha_sd),
                  na.omit(alpha_means+alpha_sd),
                  na.omit(beta_means-beta_sd),
                  na.omit(beta_means+beta_sd)
     )))
axis(1, at=seq(tek_pos_mean), labels = names(tek_pos_mean))
points(seq(tek_pos_mean)-.1, alpha_means, col=4, type='p', lty=2, lwd=2)
arrows(seq(tek_pos_mean)-.1,
       alpha_means-alpha_sd, 
       seq(tek_pos_mean)-.1,
       alpha_means+alpha_sd, 
       angle=90, code=3, length=0.1, col=4, lwd=2)
arrows(seq(tek_pos_mean)+.1,
       tek_pos_mean-tek_pos_sd, 
       seq(tek_pos_mean)+.1,
       tek_pos_mean+tek_pos_sd, 
       angle=90, code=3, length=0.1, col='grey', lwd=2)
dev.off()

## PCA

# make palette
tissue_palette <- ggplotColours(6)
names(tissue_palette) <- levels(colData(all_dds)$tissue)
tissue_palette['SpinalCord'] <- 'dodgerblue3'

dtrap_vsd <- vst(dtrap_dds)
(pca_deconv <- plotPCA.san(dtrap_vsd, intgroup = 'tissue', pchgroup = 'selection', 
                           labels = F, #paste(dtrap_vsd$tissue, dtrap_vsd$selection, sep=':'), 
                           fix.coord = F) + theme_bw() + 
    scale_color_manual(values = tissue_palette) + 
    scale_shape_manual(values = c('TRAP'=16,'whole_tissue'=2)))

pdf('_panel_figures/02_deconvPCA_byCondition.pdf', height = 4, width = 5)
pca_deconv
dev.off()

## Top500 heatmap

ntop <- 500
rv <- rowVars(assay(dtrap_vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
phOut <- pheatmap(assay(dtrap_vsd)[select, ], show_rownames = F, 
                  show_colnames=F, treeheight_row = 0,
                  annotation_col = data.frame(colData(dtrap_vsd))[,c('tissue','selection')],
                  annotation_colors = list(tissue=tissue_palette, 
                                           selection=c(TRAP='grey40', whole_tissue='grey80')),
                  file='_panel_figures/02_deconvHheatmap.pdf', height=3.5, width=6)

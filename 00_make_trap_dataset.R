source('00_generic_functions.R')

library(DESeq2)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

cts_Oct2018 <- loadCounts('../_raw_counts/Oct2018_all.counts.gz')
cts_Apr2019 <- loadCounts('../_raw_counts/Apr2019_all.counts.gz')
cts_Jan2021 <- loadCounts('../_raw_counts/Jan2021_all.counts.gz')
cts_Jun2021 <- loadCounts('../_raw_counts/Jun2021_all.counts.gz')
cts_Feb2022 <- loadCounts('../_raw_counts/Feb2022_all.counts.gz')
cts_all <- cbind(cts_Oct2018, cts_Apr2019, cts_Jan2021, cts_Jun2021, cts_Feb2022)
# rowdata
rowdata_all <- getRowdata('../_raw_counts/Oct2018_all.counts.gz')
# metadata
coldata_all <- read.csv('../_raw_counts/TRAP_metadata.csv', stringsAsFactors = FALSE)
coldata_all$genotype <- factor(coldata_all$genotype, levels = c('WT','SOD','TDP43'))
coldata_all$age <- as.character(coldata_all$age)
# filter for spinal cord
coldata_sc <- coldata_all[coldata_all$tissue == 'spinal_cord',]
rowdata_sc <- rowdata_all
cts_sc <- cts_all[,coldata_sc$sample]
coldata_sc$celltype[coldata_sc$celltype=='whole_tissue'] <- 'spinal_cord'
coldata_sc$tissue <- NULL
coldata_sc$condition <- with(coldata_sc, paste(celltype, selection, age, genotype, sep='_'))
coldata_sc$libsize <- colSums(cts_sc)

# save
save(rowdata_sc, coldata_sc, cts_sc,
     file='_rdata/cts_coldata_spinalcord.RData')

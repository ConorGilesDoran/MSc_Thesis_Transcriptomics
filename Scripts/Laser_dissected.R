# MSc Bioinformatics & Computational Biology 
# University College Cork
# Name: Conor Giles-Doran 
# Submitted: October 2022

# Analysis of Laser Captured Microdissection (LCM) Dataset - GSE20141


# Libs/Dirs/Functions ####

libs<-c('affy', 'GEOquery', 'dplyr', 'pheatmap', 'RankProd', 'tidyverse', 'limma',
        'ggplot2', 'stringr', 'limma', 'sva', 'RColorBrewer','genefilter','biomaRt',
        'clusterProfiler', 'org.Hs.eg.db', 'enrichplot', 'ggnewscale','WGCNA')

lapply(libs, library, character.only = TRUE)

# script directory
SCRIPTDIR <- "G:/My Drive/Masters/Thesis/MSc_Thesis/Scripts"

# data directory
DATADIR <- "G:/My Drive/Masters/Thesis/MSc_Thesis/Data/Laser_dissected"

# output directory
OUTDIR <- "G:/My Drive/Masters/Thesis/MSc_Thesis/Results/Laser_dissected"

if(file.exists(OUTDIR)==F){
  dir.create(OUTDIR)
}

# source necessary functions
source(paste0(SCRIPTDIR, '/Helper_functions.R'))


# Process files and merge ####
processed_list <- process_affy(DATADIR)

# merge data - in this case it is just one dataset but function will still work the same
exprs_data <- processed_list[["indiv_processed"]][["indiv_processed"]]
pheno_data <- pData(exprs_data[[1]])
feat_annotation <- processed_list[["indiv_processed"]][["feat_annotation"]]

laser_merged <- merge_process(pheno_data, exprs_data, feat_annotation)

saveRDS(laser_merged,file=paste0(OUTDIR,'/laser_merged.rds'))

# Detecting Outliers ####

lcm_array <- as.matrix(laser_merged[["expression"]])

detect_outliers(laser_merged, output_dir = OUTDIR, filename = 'lasers')

# PCA #

# analysing batch effects
pc <- prcomp(t(lcm_array), scale = TRUE)

sampleInfo <- clean_merged[['metadata']][,c('geo_accession','disease state:ch1')]

## Join the PCs to the sample information
pdf(file=paste0(OUTDIR, '/PCA_laser.pdf'), width = 10)
cbind(sampleInfo, pc$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=`disease state:ch1`, label=geo_accession)) + geom_point() + geom_text_repel()  
dev.off()


# Correlation Heatmap 
samples <- sampleInfo[order(pheno_data$`disease state:ch1`), ]
gene_matrix <- lcm_array[,rownames(samples)]
samples <- subset(samples, select = `disease state:ch1`)
colnames(samples) <- 'Disease_state'

corMatrix <- cor(gene_matrix, use = 'c')

pdf(file=paste0(OUTDIR, '/corMatrix_laser.pdf'), width = 15)
pheatmap(corMatrix,
         annotation_col = samples, fontsize = 8)#,cluster_rows = F, cluster_cols = F)
dev.off()

# Sample clustering 
sampleTree = hclust(dist(t(lcm_array)), method = "average")

pdf(file=paste0(OUTDIR, '/Sample_clust_laser.pdf'))
par(mfrow=c(1,1))
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = paste("Sample clustering on all genes"),
     xlab="", sub="", cex = 0.7)
dev.off()


# There appears to be no apparent outliers, proceeded with limma DE analysis


# Limma DE analysis ####

cleaned_annotation <- laser_merged[["annotation"]]
pheno_data <- laser_merged[["metadata"]]
disease_state <- pheno_data[,36]
colnames(pheno_data)[36] <- 'disease_state'
disease_state <- factor(pheno_data$disease_state)

pheno_data$sample_id <- rownames(pheno_data)
design <- model.matrix(~0+disease_state)

colnames(design) <- c("Control","PD")

contrasts <- makeContrasts(PD - Control, levels=design)

laser_limma <- run_limma(lcm_array, pheno_data, cleaned_annotation, design=design, contrasts=contrasts, 
                        pval=1, lfc=0, filter_dups = T, OUTDIR = OUTDIR, file_name = 'Laser_limma')

laser_limma_sig <- laser_limma[which(laser_limma$P.Value < 0.05),]


# RET & Dopaminergic MArkers ####

genes_of_int <- c('RET','TH', 
                  'ALDH1A1', 'SLC6A3', 
                  'SLC18A2', 'DDC',
                  'KCNJ6','NR4A2')


# sample_id col 
pheno_data$sample_id <- rownames(pheno_data)

# control samples
controls <- pheno_data[which(pheno_data[,'disease_state'] %in% c('control','Control','CTRL','ctrl')),]

# disease samples
disease <- pheno_data[which(!pheno_data[,'disease_state'] %in% c('control','Control','CTRL','ctrl')),]

# get DE data and order by logFC
of_interest <- laser_limma[which(laser_limma$Gene_symbol %in% genes_of_int),] 
of_interest <- of_interest[order(of_interest$logFC),]

# get array expression data
of_interest_exp <- lcm_array[of_interest$ID,]

# plot expression heatmap and save as pdf
plot_heatmap(of_interest, of_interest_exp, metadata = pheno_data, 
             topN = length(genes_of_int), status_col = 'disease_state', file_tag = 'laser_dopaminergic', OUTDIR = OUTDIR)

write.csv(of_interest, file=paste0(OUTDIR, '/dopaminergic_comparison.csv'))


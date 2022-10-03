# MSc Bioinformatics & Computational Biology 
# University College Cork
# Name: Conor Giles-Doran 
# Submitted: October 2022

# WORKFLOW FOR:
# Weighted Gene Co-expression Network Analysis
# Co-expression module profiling

# LIBs/DIRs/FUNs ####

libs<-c('dplyr', 'pheatmap', 'tidyverse', 'limma', 'ggplot2', 
        'stringr', 'limma', 'sva', 'RColorBrewer','genefilter',
        'clusterProfiler', 'org.Hs.eg.db', 'enrichplot', 'ggnewscale',
        'ggrepel','WGCNA','openxlsx','cowplot')

lapply(libs, library, character.only = TRUE)

# script directory
SCRIPTDIR <- "G:/My Drive/Masters/Thesis/MSc_Thesis/Scripts"

# data directory
DATADIR <- "G:/My Drive/Masters/Thesis/MSc_Thesis/Results"

# output directory
OUTDIR <- "G:/My Drive/Masters/Thesis/MSc_Thesis/Results/WGCNA/Combined_array"

if (file.exists(OUTDIR)==F) {
  dir.create(OUTDIR)
}

# source necessary functions
source(paste0(SCRIPTDIR, '/Helper_functions.R'))

# Read in merged_data object generated from 'DE_GO_GSEA_workflow.R' script
merged_data <- readRDS(paste0(DATADIR, '/Processed/clean_merged.rds'))


# BATCH EFFECTS ####

# PCA to analyse batch effects
pc <- prcomp(t(merged_data[['expression']]), scale = TRUE)

sampleInfo <- merged_data[['metadata']][,c('platform_id','sample_id')]

## Join the PCs to the sample information
pdf(file=paste0(OUTDIR, '/PCA_platform.pdf'), width = 10)
cbind(sampleInfo, pc$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=platform_id, label=sample_id)) + geom_point() + geom_text_repel()   
dev.off()

# remove batch effects using ComBat
mod = model.matrix(~merged_data[['metadata']]$disease_state)

batch = merged_data$metadata$platform_id

batch = as.numeric(as.factor(batch))

cleaned_array <- ComBat(merged_data$expression,batch,mod)
saveRDS(cleaned_array, file=paste0(OUTDIR, '/cleaned_array.rds'))

merged_data[['expression']] <- cleaned_array

# used to plot PCA, corMatrix, clustering
detect_outliers(merged_data, output_dir = OUTDIR, filename = 'WGCNA_batch_analysis')


# PROCESSING OUTLIERS ####

# If other outliers present, can remove them here.
# Sample clustering dendrogram is presented and cutoff can be selected.

exprs_data <- t(merged_data[['expression']])
anno <- merged_data[['annotation']]

# any further outliers to be removed?
# in this case - NO - so select a value larger than the dendrogram tree height
exprs_data <- filter_outliers(exprs_data, output_dir = paste0(OUTDIR,'/Plots'),nametag = 'samples_dendro')

# metadata samples must match array data
metadata <- merged_data[['metadata']][rownames(exprs_data),]

# Assess the disease status of clustered samples
# make binary vector
status <- binarizeCategoricalVariable(metadata$disease_state)
colnames(status) <- 'PD vs CTRL'
rownames(status) <- rownames(metadata)

# Dendro with trait (status) info
pdf(file=paste0(OUTDIR, '/Sample_dendro_with_trait.pdf'), height = 10, width = 15)
plotClusterTreeSamples(datExpr=exprs_data, y=status)
dev.off()


# SOFT THRESHOLD ####

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))

sft <- choose_sft(exprs_data, powers, network_type = "signed", 
                  output_dir=OUTDIR, file_name='sft_choice')


# selected power threshold
power = sft$powerEstimate


# NETWORK CONSTRUCTION ####

picked_power <- power
temp_cor <- cor       
cor <- WGCNA::cor         # Fixes a namespace conflict issue

# check scale-free topology 
k=softConnectivity(datE=exprs_data,power=picked_power)

# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
pdf(paste0(OUTDIR,'/scale_free_top.pdf'), width=15)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
dev.off()

# network construction code inspired by https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
# automatic network construction and module detection 
# Note: long run time but < 1 hour
full_net <- blockwiseModules(exprs_data,               
                          
                          # Adjacency Function 
                          power = picked_power,                
                          networkType = "signed",
                          corType = "bicor",
                          
                          # Tree and Block Options 
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          minModuleSize = 30,
                          maxBlockSize = 20000, #20k to do all in one block
                          
                          # Module Adjustments
                          reassignThreshold = 0,
                          mergeCutHeight = 0.2,
                        
                          # TOM Archive the run results in TOM file 
                          saveTOMs = T,
                          saveTOMFileBase = paste0(OUTDIR,"/combined_array_modules"),
                          
                          # Output Options
                          numericLabels = T,
                          verbose = 3)

# Save network
saveRDS(full_net, file=paste0(OUTDIR,'/full_net.rds'))

# Convert labels to colors for plotting
mergedColors = labels2colors(full_net$colors)

# Plot the dendrogram and the module colors underneath
pdf(file=paste0(OUTDIR, '/combined_array_modules.pdf'), width = 12)
plotDendroAndColors(
  full_net$dendrograms[[1]],
  mergedColors,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
dev.off()


# MODULE INFO ####

# extract MEs, Module Membership(Signed Eigengene Connectivity), Hub genes
module_details <- module_details(exprs_data,full_net)

# Save RDS object
saveRDS(module_details, paste0(OUTDIR,'/module_details.rds'))

# Write each object to csv file also
for(x in names(module_details)){
  openxlsx::write.xlsx(module_details[[x]], file = paste0(OUTDIR,'/',x,'.xlsx')) 
}

# example of module eigengene plot
MEs <- module_details[['MEs']]

ME_magenta <- MEs[,'MEmagenta']

pdf(paste0(OUTDIR,'/Plots/Magenta_ME.pdf'), width=10)
barplot(ME_magenta, col='magenta', main="", cex.main=2,
        ylab="Eigengene expression",xlab="Array sample")
dev.off()


# COR WITH DISEASE STATUS #### 

# Define numbers of genes and samples
nSamples = nrow(exprs_data)

# Module colours
moduleColors <- labels2colors(full_net$colors)

# Binarize disease status column
status <- binarizeCategoricalVariable(metadata$disease_state)
colnames(status) <- 'PD vs CTRL'
rownames(status) <- rownames(metadata)

# Corelation with MEs
moduleTraitCor = WGCNA::cor(MEs, status, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Plotting results

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
pdf(file=paste0(OUTDIR,'/Plots/Module_vs_status.pdf'), width = 6)
par(mar = c(6, 12, 4, 4))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = 'Disease Status: PD vs Control',
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               xLabelsAdj = 0.5,
               xLabelsAngle = 0,
               main = paste("Module-Disease status relationships"))
dev.off()


# VISUALIZATION ####

MET = orderMEs(cbind(MEs, status))

# Plot the relationships among the eigengenes and the trait

# Plot the dendrogram
pdf(file=paste0(OUTDIR, '/Eigengene_trait_relationships.pdf'), height=10,width=12)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(10,10,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

dev.off()


# GO ANALYSIS ####

cleaned_annotation <- merged_data[["annotation"]]

module_genes <- module_details[['module_genes']]

# Find modules with absolute correlation > 0.25 with PD Disease Status.
mods_interest <- gsub(pattern = 'ME', '', names(moduleTraitCor[which(abs(moduleTraitCor) > 0.25),]))

# Perform Go over-rep analysis (CC, BP, MF) for each module of interest.
coExp_GO <- lapply(mods_interest, function(x) GO_overep(module_genes[[x]][,'ENTREZ_GENE_ID'],ref_list = unique(anno[,'ENTREZ_GENE_ID']),FDR = 0.2))
names(coExp_GO) <- mods_interest
saveRDS(coExp_GO,paste0(OUTDIR,'/coExpGO.rds'))

# Write to file
CoExp_BP <- lapply(coExp_GO, function(x) x$BP@result)
CoExp_CC <- lapply(coExp_GO, function(x) x$CC@result)
CoExp_MF <- lapply(coExp_GO, function(x) x$MF@result)

openxlsx::write.xlsx(CoExp_BP, file = paste0(OUTDIR,'/CoExpBP_results.xlsx'))
openxlsx::write.xlsx(CoExp_CC, file = paste0(OUTDIR,'/CoExpCC_results.xlsx'))
openxlsx::write.xlsx(CoExp_MF, file = paste0(OUTDIR,'/CoExpMF_results.xlsx'))


# KEGG ANALYSIS ####

# Same over-representation analysis but with KEGG

coExp_KEGG <- lapply(mods_interest, function(x) enrichKEGG(gene = module_genes[[x]][,'ENTREZ_GENE_ID'],
                         universe=unique(anno[,'ENTREZ_GENE_ID']),
                         organism = 'hsa',
                         pvalueCutoff = 0.05))


names(coExp_KEGG) <- mods_interest
saveRDS(coExp_KEGG, paste0(OUTDIR,'/coExpKEGG.rds'))

KEGG_coExp <- lapply(coExp_KEGG, function(x) x@result)
openxlsx::write.xlsx(KEGG_coExp, file = paste0(OUTDIR,'/CoExpKEGG_results.xlsx'))


# PLOTTING PURPLE MODULE GO ####

# PURPLE #
BP_pairs_purple <- pairwise_termsim(coExp_GO[["purple"]][["BP"]])

pdf(paste0(OUTDIR,'/Plots/Purple_emapplot.pdf'), height = 10, width = 12)
emapplot(BP_pairs_purple, node_label='none', cex_line = 0.2, showCategory = 50,
         cex_category=0.8, group_category = T, nWords = 5, group_legend = T,
         cex_label_category = 0.5)
dev.off()






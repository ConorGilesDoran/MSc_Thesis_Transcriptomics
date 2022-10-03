# MSc Bioinformatics & Computational Biology 
# University College Cork
# Name: Conor Giles-Doran 
# Submitted: October 2022

# WORKFLOW FOR:
# Data pre-processing
# Data merging
# Differential gene expression analysis 
# GO Over-representation analysis
# Gene set erichment analysis


# Libraries/Directories/Functions ####

libs<-c('affy', 'GEOquery', 'dplyr', 'pheatmap', 'tidyverse', 'limma',
        'ggplot2', 'stringr', 'limma', 'sva', 'RColorBrewer','genefilter','biomaRt',
        'clusterProfiler', 'org.Hs.eg.db', 'enrichplot', 'ggupset', 'ggnewscale','ggrepel',
        'WGCNA','openxlsx','cowplot')

lapply(libs, library, character.only = TRUE)


# script directory
SCRIPTDIR <- "G:/My Drive/Masters/Thesis/MSc_Thesis/Scripts"

# data directory
DATADIR <- "G:/My Drive/Masters/Thesis/MSc_Thesis/Data/Microarray"

# output directory
OUTDIR <- "G:/My Drive/Masters/Thesis/MSc_Thesis/Results"

# source necessary functions
source(paste0(SCRIPTDIR, '/Helper_functions.R'))


# Pre-processing ####

# read in .CEL files, merge same platform, perform RMA
processed_list <- process_affy(file_dir = DATADIR)
saveRDS(processed_list, file = paste0(OUTDIR, '/Processed/full_processed_list.rds'))

# source 'process_meta' script
# This script manually extracts and combines all relevant metadata for downstream analysis.
# 144 samples remaining - All SN tissue, all either control or PD.
# Should return and save 'cleaned_pheno.rds' object
source(paste0(SCRIPTDIR, '/process_meta.R'))

# Merging expression data and filtering of problematic probes
array_data <- processed_list[['rma_platform']]

feat_annotation <- processed_list[["indiv_processed"]][["feat_annotation"]]

merged <- merge_process(cleaned_pheno=cleaned_pheno, exprs_data=array_data, feat_annotation=feat_annotation)

saveRDS(merged, file=paste0(OUTDIR, '/Processed/merged_data.rds'))


# Detecting Outliers ####

detect_outliers(merged, output_dir = paste0(OUTDIR,'/Batch_analysis'),filename = 'outlier')

# Re-run pre-processing but with outliers removed
outliers <- c('GSM506020','GSM508732')

clean_processed_list <- process_affy(file_dir = DATADIR,outliers = outliers)
saveRDS(clean_processed_list, file=paste0(OUTDIR, '/Processed/clean_processed_list.rds'))

array_data <- clean_processed_list[['rma_platform']]
feat_annotation <- clean_processed_list[["indiv_processed"]][["feat_annotation"]]

# Remove outliers from metadata
cleaned_pheno <- cleaned_pheno[!rownames(cleaned_pheno) %in% c('GSM506020','GSM508732'),]

# Re-run merging procedure
clean_merged <- merge_process(cleaned_pheno=cleaned_pheno, exprs_data=array_data, feat_annotation=feat_annotation)
saveRDS(clean_merged, file=paste0(OUTDIR, '/Processed/clean_merged.rds'))

# repeat outlier detection just to see difference
detect_outliers(clean_merged, output_dir = paste0(OUTDIR,'/Batch_analysis'),filename = 'outlier_rem')


# Limma DE analysis ####

# Extract combined array data
arrays_combined <- as.matrix(clean_merged[["expression"]])
saveRDS(arrays_combined, file = paste0(OUTDIR, '/Processed/arrays_combined.rds'))

# Create new directory for limma results
limma_dir <- paste0(OUTDIR, '/limma_results')

if (file.exists(limma_dir)==F) {
  dir.create(limma_dir)
}

# Extract probe annotation and metadata 
cleaned_annotation <- clean_merged[["annotation"]]
cleaned_pheno <- clean_merged[['metadata']]

# Factorize variables to be included in limma design
platform <- factor(cleaned_pheno$platform_id)
disease_state <- factor(cleaned_pheno$disease_state)

# Establish limma design matrix
design <- model.matrix(~0+disease_state+platform)

colnames(design) <- c("Control","PD", 'GPL96')

# PD vs control is the comparison being investigated
contrasts <- makeContrasts(PD - Control, levels=design)

# Run the limma analysis

# full_limma will consist of the entire set of results (significant and non-significant)
full_limma <- run_limma(arrays_combined, cleaned_pheno, cleaned_annotation, design=design, contrasts=contrasts, 
                           pval=1, lfc=0, filter_dups = T, OUTDIR=limma_dir, file_name='full_limma' )

# signif_limma will consist of just significant results (FDR=0.05, lfc=0.4)
signif_limma <- run_limma(arrays_combined, cleaned_pheno, cleaned_annotation, design=design, contrasts=contrasts, 
                        pval=0.05, lfc=0.4, filter_dups = T, OUTDIR=limma_dir, file_name='signif_limma')



# PLOT TOP 30 GENES AND SAVE #

# ComBat Removal of batch effect for plotting expression of significant genes
mod = model.matrix(~clean_merged[['metadata']]$disease_state)

batch = clean_merged$metadata$platform_id

batch = as.numeric(as.factor(batch))

cleaned_array <- ComBat(clean_merged$expression,batch,mod)

# plotting
plot_heatmap(signif_limma, cleaned_array, cleaned_pheno, 
             topN = 30, file_tag = 'combined', status_col = 'disease_state',OUTDIR = limma_dir)



# RET Correlation Analysis ####

# perform correlation analysis and linear regression between RET and dopaminergic markers
main_gene <- 'RET'
genes_of_interest <- c('RET','TH','ALDH1A1',
                       'SLC6A3','SLC18A2', 'DDC',
                       'KCNJ6','NR4A2', 'SNCA')

output_dir <- paste0(OUTDIR, '/Dopaminergic_markers/')
pheno_data <- clean_merged$metadata
file_tag <- 'RET_dopaminergic'

RET_lm <- run_correlation(signif_limma, genes_of_interest, main_gene,
                          cleaned_array, pheno_data, file_tag, output_dir, status_col = 'disease_state')



# GO Over-representation ####

# Split significant limma results into up and downregulated DEGs
upreg <- signif_limma[which(signif_limma$logFC>0),]

downreg <- signif_limma[which(signif_limma$logFC<0),]

# reference list is all unique Entrez gene IDs from the cleaned array
ref_list <- unique(cleaned_annotation[,'ENTREZ_GENE_ID'])

# Upregulated GO over-representation
GO_up <- GO_overep(upreg$ENTREZ_GENE_ID, ref_list=ref_list, keytype='ENTREZID', organism_db=org.Hs.eg.db,
          FDR=0.05)

# Downregulated GO over-representation
GO_down <- GO_overep(downreg$ENTREZ_GENE_ID, ref_list=ref_list, keytype='ENTREZID', organism_db=org.Hs.eg.db,
                     FDR=0.05)

# GO results list
GO_results <- list(up = GO_up,
                   down = GO_down)


# Save results
# define sheet names for each data frame
dataset_names <- unlist(GO_results)
dataset_names <- lapply(dataset_names, function(x) x@result)

# export each data frame to separate sheets in same Excel file
GO_dir <- paste0(OUTDIR, '/GO_over_rep')

if (file.exists(GO_dir)==F) {
  dir.create(GO_dir)
}

openxlsx::write.xlsx(dataset_names, file = paste0(GO_dir,'/GO_results.xlsx')) 
saveRDS(GO_results, file=paste0(GO_dir, '/over_rep.rds'))


# PLOTTING GO BP RESULTS

# DOWNREG DEGs
down_logFC <- downreg$logFC
names(down_logFC) <- downreg$Gene_symbol

pdf(file=paste0(GO_dir,'/Downreg_BP_plots.pdf'), height = 10, width=14)
# calculate pairwise distances
to_plot <- pairwise_termsim(GO_results[["down"]][["BP"]])

# barplot
barplot(GO_results[["down"]][["BP"]],showCategory=20, 
        title='Top 20 over-represented biological processes for downregulated DEGs')

# enrichment map
emapplot(to_plot, cex_label_category = 1,
         cex_category = 1,
         cex_line = 1, layout = 'kk')

# enriched GO directed acyclic graph
goplot(to_plot, showCategory = 5)

# gene-concept network plot
cnetplot(to_plot, categorySize="pvalue", foldChange=down_logFC, showCategory = 5,
         cex_gene=0.2)

# grouping tree diagram
treeplot(to_plot, showCategory = 50, cex_category = 0.7)

dev.off()


# UPREG DEGs

pdf(file=paste0(GO_dir,'/Upreg_BP_plots.pdf'), height = 10, width=14)
up_logFC <- upreg$logFC
names(up_logFC) <- upreg$Gene_symbol

# calculate pairwise distances
to_plot <- pairwise_termsim(GO_results[["up"]][["BP"]])

# barplot
barplot(GO_results[["up"]][["BP"]],showCategory=20, 
        title='Top 20 over-represented biological processes for upregulated DEGs')

# enrichment map
emapplot(to_plot, cex_label_category = 1,
         cex_category = 1,
         cex_line = 1, layout = 'kk')

# enriched GO directed acyclic graph
goplot(to_plot, showCategory = 5)

# gene-concept network plot
cnetplot(to_plot, categorySize="pvalue", foldChange=up_logFC, showCategory = 5,
         cex_gene=0.2)

# grouping tree diagram
treeplot(to_plot, showCategory = 50, cex_category = 0.7, nCluster=6)

dev.off()


# GSEA ####

GSEA_results <- list()

# full categories

full_categories <- c('H','C5')

for(cat in full_categories){

  gsea <- run_gsea(full_limma, entrez_col='ENTREZ_GENE_ID',rank_col='t',species='Homo sapiens', cat=cat)
  GSEA_results[[cat]] <- gsea
}


# sub categories
specific_cat <- c('C3','C2')
sub_cat <- c('TFT:GTRD','CP:KEGG')

for(cat in 1:length(specific_cat)){
  
  gsea <- run_gsea(full_limma, entrez_col='ENTREZ_GENE_ID',rank_col='t',species='Homo sapiens', 
                   cat=specific_cat[cat], subcategory=sub_cat[cat])
  
  name <- paste0(specific_cat[cat],'_',strsplit(sub_cat[cat],':')[[1]][[2]])
  GSEA_results[[name]] <- gsea
}


# saving 
GSEA_dir <- paste0(OUTDIR, '/GSEA')

if (file.exists(GSEA_dir)==F) {
  dir.create(GSEA_dir)
}

# convert Entrez ID number to gene symbol
GSEA_readable <- lapply(GSEA_results, function(x) setReadable(x, 'org.Hs.eg.db','ENTREZID'))

# save
GSEA_files <- lapply(GSEA_readable, function(x) x@result)
saveRDS(GSEA_results, file=paste0(GSEA_dir,'/GSEA_results.rds'))
openxlsx::write.xlsx(GSEA_files, file = paste0(GSEA_dir,'/GSEA_results.xlsx')) 


# Hallmark GSEA plots ####

# GSEA PLOT
of_interest <- c(1,3,4,5)

pdf(file=paste0(GSEA_dir,'/Hall_GSEA_plot.pdf'), height=15, width=18)

gseaplot2(GSEA_results[["H"]], geneSetID = of_interest,
          base_size = 20)

pp <- lapply(c(1,3,4,10), function(i) {
  anno <- GSEA_results[["H"]][i, c("NES", "pvalue", "p.adjust")]
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  
  gsearank(GSEA_results[["H"]], i, GSEA_results[["H"]][i, 2]) + xlab(NULL) +ylab(NULL) +
    annotate("text", 10000, GSEA_results[["H"]][i, "enrichmentScore"] * .75, label = lab, hjust=0.5, vjust=0.5)
})
plot_grid(plotlist=pp, ncol=1)

dev.off()


# Ridgeline plot
pdf(file=paste0(GSEA_dir,'/Hall_ridgeplot.pdf'), height=10, width=15)
ridgeplot(GSEA_results[["H"]], label_format = 100)
dev.off()


# GTRD and KEGG ####

# GTRD
pdf(file=paste0(GSEA_dir,'/GTRD_plots.pdf'), height=10, width=20)
p1 = dotplot(GSEA_results[["C3_GTRD"]], showCategory = 17, font.size = 20, orderBy= 'NES')

p2 <- ggplot(GSEA_results[["C3_GTRD"]]@result, ggplot2::aes(reorder(ID, NES), NES)) +
  geom_col(aes(fill = p.adjust)) +
  scale_fill_gradient(low='red', high='blue') +
  coord_flip() +
  labs(x = '',
       y = "Normalized Enrichment Score") +
  theme_minimal() +
  theme(text=element_text(size=20,colour = 'black'))


cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])
dev.off()

# KEGG
pdf(file=paste0(GSEA_dir,'/KEGG_plots.pdf'), height=10, width=15)

dotplot(GSEA_results[["C2_KEGG"]], showCategory = 33, label_format = 50, font.size= 20, orderBy= 'NES')

ggplot(GSEA_results[["C2_KEGG"]]@result, ggplot2::aes(reorder(ID, NES), NES)) +
  geom_col(aes(fill = p.adjust)) +
  scale_fill_gradient(low='red', high='blue') +
  coord_flip() +
  labs(x = '',
       y = "Normalized Enrichment Score") +
  theme_minimal() +
  theme(text=element_text(size=20,colour = 'black'))

ridgeplot(GSEA_results[["C2_KEGG"]], label_format = 100)

dev.off()


# GO Plots ####

# GO has 1079 terms, break down into sections

BP_enrich <- GSEA_results[["C5"]]@result[grep('GOBP_',GSEA_results[["C5"]]@result$ID),]
CC_enrich <- GSEA_results[["C5"]]@result[grep('GOCC_',GSEA_results[["C5"]]@result$ID),]
MF_enrich <- GSEA_results[["C5"]]@result[grep('GOMF_',GSEA_results[["C5"]]@result$ID),]
HP_enrich <- GSEA_results[["C5"]]@result[grep('HP_',GSEA_results[["C5"]]@result$ID),]

GO_breakdown <- list(BP_enrich,CC_enrich,MF_enrich,HP_enrich)
names(GO_breakdown) <- c("BP_enrich","CC_enrich","MF_enrich","HP_enrich")

openxlsx::write.xlsx(GO_breakdown, file = paste0(GSEA_dir,'/GSEA_GO.xlsx')) 

# Plots for BP enrichment
BP_signif <- BP_enrich[BP_enrich$p.adjust<0.01,]
BP_signif <- BP_signif[order(BP_signif$NES, decreasing = TRUE),]

# Count number of genes enriched for each term and add to table for plotting
core_enrichments <- BP_signif$core_enrichment
enrich_count <- c()

for(e in 1:length(core_enrichments)){
  
  count <- length(unlist(strsplit(core_enrichments[e],'/')))
  
  enrich_count <- c(enrich_count, count)
}

BP_signif$EnrichedGenes <- enrich_count

# lollipop plot
# overall view of GO results (<0.01 FDR)
# colour of points based on adjusted p value
# size of points based on number of enriched genes

pdf(paste0(GSEA_dir,'/GSEA_GOBP_plot.pdf'), height = 8, width = 15)
ggplot(BP_signif, aes(x=ID, y=NES)) +
  geom_segment(aes(x=ID, xend=ID, y=0, yend=NES),colour='grey',
               size=0.6, alpha=0.9) +
  geom_point(aes(colour=p.adjust, size=EnrichedGenes)) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x=element_blank(),
    text = element_text(size = 15)
  ) +
  geom_hline(yintercept=0, colour='red') +
  xlab("GO Biological Process") +
  ylab("NES")
dev.off()


# GO Text Pattern Search ####
BP_pos_NES <- BP_signif[BP_signif$NES>0,]
BP_neg_NES <- BP_signif[BP_signif$NES<0,]

prompt <- c('IMMUN', 'INFLAM', 'INTERLEUKIN', 'LYMPHOCYTE','LEUKOCYTE', '_T_CELL', '_B_CELL',
            'DIFFERENTIATION', 'METAB','DNA','ORGANIZATION', 'NEGATIVE_REGULATION',
            'GLIAL', 'INTERFERON','OSSIFICATION', 'CYTOKINE', 'POSITIVE_REGULATION', 'REMODELLING',
            'RESORPTION', 'DEVELOP', 'SKELET', 'DOPA', 'SYNAP','MITOCHOND', 'BIOSYNTH', 
            'TRANSPORT', 'VESIC','BEHAVIOR', 'PROTEIN', 'GROWTH', 
            'TOLL_LIKE', 'RESPONSE', 'EXOCYTOSIS', 'LYSIS', 'ACTIVATION')

# occurrence of text pattern in results with NES > 0
BP_pos_NES_terms <- list()
for(p in prompt){
  
  count <- length(grep(p, BP_pos_NES$ID))
  
  BP_pos_NES_terms[[p]] <- count
  
}

# occurrence of text pattern in results with NES < 0
BP_neg_NES_terms <- list()
for(p in prompt){
  
  count <- length(grep(p, BP_neg_NES$ID))
  
  BP_neg_NES_terms[p] <- count
  
}

# tally the results for each and plot
BP_neg_NES_terms <- unlist(BP_neg_NES_terms)
BP_pos_NES_terms <- unlist(BP_pos_NES_terms)

counts <- c(BP_neg_NES_terms, BP_pos_NES_terms)

enrichment_tally <- data.frame(BP_term = rep(names(BP_neg_NES_terms),2),
                               NES = c(rep('negative', length(BP_pos_NES_terms)), rep('positive', length(BP_pos_NES_terms))),
                               tally = counts)



# STACKED BARPLOT

my_color <- c("#B4464B", "#4682B4")

enrichment_tally$NES <- as.factor(enrichment_tally$NES)

names(my_color) <- levels(enrichment_tally$NES)

pdf(paste0(GSEA_dir,'/GSEA_GOBP_tally.pdf'), height = 8, width = 10)
ggplot(enrichment_tally, aes(fill=NES, y=tally, x=BP_term)) + 
  geom_bar(position='stack', stat='identity') +
  scale_fill_manual(values = my_color,
                    name = "NES Sign")+ 
  labs(y = "\nNo. of enriched BP terms (FDR<0.01)", x = "Text Pattern Searched\n") +
  theme_minimal() +
  theme(text = element_text(size=20))+
        #axis.text.x = element_text(size = 12),
        #axis.text.y = element_text(size=12),
        #axis.title = element_text(size = 13)) 
  coord_flip()
dev.off()


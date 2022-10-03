# MSc Bioinformatics & Computational Biology 
# University College Cork
# Name: Conor Giles-Doran 
# Submitted: October 2022

# Functions for analyses performed in 'WGCNA_workflow.R', 'DE_GO_GSEA_workflow.R' and 'Laser-dissected' scripts.

# PRE-PROCESSING ####
# Untars raw .tar files, reads in .CEL files and creates Affybatch objects.
# Identifies platforms used and merges raw data from same platform.
# Option to filter known outliers before RMA.
# Performs RMA of merged data objects. 
# Performs RMA of individual raw data objects. 
# Annotates each individual RMA data object with GEO metadata.
# Extracts one comprehensive probe feature annotation object. 
# Returns a list containing platform specific RMA processed data, individual
# RMA processed data with GEO sample metadata included and one probe feature annotation dataset.

# file_dir = directory where raw .tar files are located.
# outliers option = if outliers are known, a vector of outlier sample IDs can be included here.

process_affy <- function(file_dir, outliers =  NULL){
  
  # untar raw files and create directory for each set of .CEL files
  print('Processing .CEL files...')
  dir_files <- list.files(file_dir)
  
  raw_data <- dir_files[grep('_RAW.tar', dir_files)]
  
  new_dirs <- gsub('_RAW.tar', '', raw_data)
  
  for(i in 1:length(raw_data)){
    
    untar(paste0(file_dir, '/', raw_data[i]), exdir = paste0(file_dir, '/', new_dirs[i]))
    
  }
  
  # read in .CEL files and create AffyBatch object for each
  raw_micro <- list()
  for(i in new_dirs){
    raw_micro[[i]] <- ReadAffy(celfile.path = paste0(file_dir, '/', i))
  }
  
  
  # remove outliers before RMA if they are known
  if(is.null(outliers)==F){
    
    print(paste0('Removing outlier: ', outliers))
    out <- paste0(outliers, '.CEL.gz')
    
    samples <- lapply(raw_micro, colnames)
    
    for(i in 1:length(out)){

      dataset <- grep(out[i], samples)
      
      remove_ind <- which(colnames(raw_micro[[dataset]]) == out[i])
      
      raw_micro[[dataset]] <- raw_micro[[dataset]][,-remove_ind]
      
    }
  }
  
  
  # what platforms were used across the datasets
  platforms <- unlist(lapply(raw_micro, function(x) x@cdfName))
  print('Platforms:')
  print(platforms)
  
  # get a list of each platform, within which is another list of the associated AffyBatch objects 
  same_platform <- list()
  for(i in unique(platforms)){
    
    same_platform[[i]] <- list()
    
    for(p in 1:length(platforms)){
      
      same_platform[[platforms[[p]]]][[names(platforms)[p]]] <- raw_micro[[names(platforms)[p]]]
      
    }
  }
  
  
  print('Merging data of same platforms...')
  # merge each set of platform specific AffyBatch objects together
  merged_platform <- lapply(same_platform, function(p) Reduce(function(x,y) merge(x = x, y = y, 
                                                                      annotation = paste(annotation(x), annotation(y))),p))
  
  processed_list <- list()
  if(length(merged_platform) != 1){
    # rma for matching platforms
    print('Robust multi-array averaging (RMA) of matching platforms...')
    rma_platform <- lapply(merged_platform, rma)
    processed_list[['rma_platform']] <- rma_platform
  }
  
  print('Robust multi-array averaging (RMA) of individual datasets...')
  # rma with each indiviudal dataset
  indiv_rma <- lapply(raw_micro, rma)
  
  # add pheno and feature data 
  print('Adding GEO metadata...')
  
  if(is.null(outliers)==T){
    indiv_processed <- add_GEO_meta(new_dirs, indiv_rma)
  }else{
    indiv_processed <- add_GEO_meta(new_dirs, indiv_rma, outliers = outliers)
  }
  processed_list[['indiv_processed']] <- indiv_processed
                      
  return(processed_list)
}


# GEO METADATA ####
# Extracts sample metadata/probe feature data.
# Used within the process_affy function.
# Certain GEO data series contain more than one microarray dataset, need to 
# extract the correct information.
# Option to remove outliers
# Returns individually processed dataset(s) with GEO metadata included and
# a probe feature annotation dataset. 

# GEO_ids = Vector containing accession numbers of the GEO data series 
# indiv_rma = list of RMA processed raw data from individual GEO data series
# outliers = if outliers are known, a vector of outlier sample IDs can be included here.

add_GEO_meta <- function(GEO_ids, indiv_rma, outliers=NULL){
  
  GEO_info <- list()
  
  # format GDS and GSE differently
  for(id in GEO_ids){
    
    if(grepl("GSE", id, fixed=TRUE)==T){
      
      GEO_info[[id]] <- getGEO(id)
      
    }else if(grepl("GDS", id, fixed=TRUE)==T){
      
      gds <- getGEO(id)
      
      GEO_info[[id]] <- GDS2eSet(gds,do.log2=F)
      
    }
    
  }
  
  
  # extract all datasets into a list
  datasets <- unlist(GEO_info, use.names = T)
  
  # remove outliers if they are known
  if(is.null(outliers)==F){
    
    samples <- lapply(datasets, colnames)
    
    for(i in 1:length(outliers)){
      
      dataset <- grep(outliers[i], samples)
      
      remove_ind <- which(colnames(datasets[[dataset]]) == outliers[i])
      
      datasets[[dataset]] <- datasets[[dataset]][,-remove_ind]
      
    }
  }
  
  
  # edit column names so all are the same format
  for(i in 1:length(indiv_rma)){
    
    colnames(indiv_rma[[i]]) <- gsub('.CEL.gz|.cel.gz', "", colnames(indiv_rma[[i]]))
    colnames(indiv_rma[[i]]) <- gsub('_.*', "", colnames(indiv_rma[[i]]))
  }
  
  # extract column names and identify matching datasets
  GEO_cols <- lapply(datasets, function(x) sort(colnames(exprs(x))))
  data_cols <- lapply(indiv_rma, function(x) sort(colnames(exprs(x))))
  
  print(paste0('Identified ', length(GEO_cols),' GEO datasets.'))
  print('Finding the matching datasets...')
 
  # find matching samples and add associated metadata to match data list 
  match_data <- list()
  print('Adding phenotype and feature metadata...')
  for(i in 1:length(data_cols)){
    
    matching_set <- which(lapply(GEO_cols, function(x) identical(x, data_cols[[i]])) == TRUE)
    
    match_data[[names(data_cols[i])]] <- datasets[[matching_set]]
    
    pData(indiv_rma[[names(match_data[i])]]) <- pData(match_data[[i]])
  }
  
  probe_set <- unlist(lapply(match_data, nrow))
  
  feat_annotation <- fData(match_data[[which.max(probe_set)]])
  fData(indiv_rma[[which.max(probe_set)]]) <- feat_annotation
  
  indiv_processed <- list(indiv_processed = indiv_rma,
                          feat_annotation = feat_annotation)
  
  return(indiv_processed)
}


# MERGING ARRAYS ####
# Performs merging procedure of platform specific RMA processed data
# Removes probes mapping to multiple gene symbols.
# Attempts to annotate probes mapping to no gene symbol or Entrez ID.
# Removes probes still mapping to no gene symbol or Entrez ID.
# Returns list containing cleaned metadata, cleaned combined array data and 
# cleaned probe feature annotation data.

# cleaned_pheno = pre-processed sample metadata 
# exprs_data = list containing platform specific RMA processed data
# feat_annotation = probe feature annotation data

merge_process <- function(cleaned_pheno, exprs_data, feat_annotation){
  
  print('Merging expression data...')
  # extract expression data
  array_data <- lapply(exprs_data, exprs)
  
  # what rows must be kept for the combined dataset
  keep_samples <- rownames(cleaned_pheno)
  
  # create id column for merging purposes
  array_data <- lapply(array_data, function(x){x <- as.data.frame(x)
                                               x$id <- rownames(x)
                                               return(x)})
  
  arrays_combined <- Reduce(function(x,y) merge(x = x, y = y, by = "id"), 
                            array_data)
  
  # set id as rownames
  rownames(arrays_combined) <- arrays_combined[,"id"]
  arrays_combined <- dplyr::select(.data = arrays_combined, !"id")
  
  # format sample names
  colnames(arrays_combined) <- gsub('.CEL.gz|.cel.gz', '', colnames(arrays_combined))
  colnames(arrays_combined) <- gsub("_.*", "", colnames(arrays_combined))
  
  # retain the samples from cleaned pheno
  arrays_combined <- arrays_combined[,keep_samples]
  
  # Feature annotation data 
  print('Processing feature annotation data...')
  anno <- feat_annotation[which(feat_annotation$ID %in% rownames(arrays_combined)),]
  
  # reformat Gene Symbol column name
  gene_name_col <- grep('Gene symbol|gene symbol|Gene Symbol|Gene_symbol|gene_symbol', colnames(anno))
  
  colnames(anno)[gene_name_col] <- 'Gene_symbol'
  
  # Try to annotate missing probe information using ensemb biomaRt
  missing <- anno[which(anno$Gene_symbol == '' & anno$ENTREZ_GENE_ID == ''),]$ID
  
  print(paste0('Identified ', length(missing), ' probes with no gene symbol AND no Entrez ID.'))
  print('Attempting to annotate using ensembl biomaRt...')
  
  mart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
  
  entrez_IDs <- select(mart, keys=missing, columns=c('affy_hg_u133a','hgnc_symbol','entrezgene_id'),
                       keytype='affy_hg_u133a')
  
  # remove NAs
  entrez_IDs <- na.omit(entrez_IDs)
  
  # remove the ids associated with duplicate ids as these will indicate multiple mappings
  dup_IDs <- entrez_IDs[duplicated(entrez_IDs$affy_hg_u133a),]$affy_hg_u133a
  
  if(length(dup_IDs)!=0){
    
    entrez_IDs <- entrez_IDs[-which(entrez_IDs$affy_hg_u133a %in% dup_IDs),]
  }
  
  print(paste0(nrow(entrez_IDs), ' complete annotations, ', length(missing)-nrow(entrez_IDs),' remaining blank entries will be removed...'))
  
  # append to anno 
  rownames(anno) <- anno$ID
  rownames(entrez_IDs) <- entrez_IDs$affy_hg_u133a
  anno[rownames(entrez_IDs),]$ENTREZ_GENE_ID <- entrez_IDs[rownames(entrez_IDs),]$entrezgene_id
  anno[rownames(entrez_IDs),]$Gene_symbol <- entrez_IDs[rownames(entrez_IDs),]$hgnc_symbol
  
  # remove those that are still blank and those that map to multiple probe IDs
  still_blank <- anno[which(anno$Gene_symbol == '' & anno$ENTREZ_GENE_ID == ''),]$ID
  
  multi_maps <- anno[grep('/', anno$Gene_symbol),]$ID
  
  print(paste0('Removing ', length(multi_maps), ' probes that map to multiple genes...'))
 
  anno_filt <- anno[-which(anno$ID %in% c(still_blank, multi_maps)),]
  
  entrez_blank <- anno_filt[which(anno_filt$ENTREZ_GENE_ID == ''),]$ID  
  
  print(paste0('Removing ', length(entrez_blank), ' probes that still have no Entrez ID...'))
  
  anno_filt <- anno_filt[-which(anno_filt$ID %in% entrez_blank),]
  
  arrays_combined <- arrays_combined[which(rownames(arrays_combined) %in% anno_filt$ID),]
  
  clean_merged <- list()
  clean_merged[['metadata']] <- cleaned_pheno
  clean_merged[['expression']] <- arrays_combined
  clean_merged[['annotation']] <- anno_filt
  
  return(clean_merged)
}


# DETECT OUTLIERS ####

# NOTE: FOR MERGED DATA ONLY #
# Function to aid in the visual identification of sample outliers.
# Performs PCA and plots the result.
# Plots a sample correlation matrix.
# Plots a sample dendrogram.

# clean merged = list containing cleaned (pre-processed) metadata, array data and feature data
# outpur_dir = directory to which plots should be saved.
# filename = name tag associated with each pdf plot file. 

# CODE AIDED BY: https://sbc.shef.ac.uk/geo_tutorial/tutorial.nb.html

detect_outliers <- function(clean_merged, output_dir, filename){
  
  arrays_combined <- as.matrix(clean_merged[["expression"]])
  
  # analysing batch effects
  pc <- prcomp(t(arrays_combined), scale = TRUE)
  
  sampleInfo <- clean_merged[['metadata']][,c('platform_id','sample_id','dataset')]
  
  pdf(file=paste0(output_dir, '/',filename, '_clustering.pdf'), height=10, width = 15)
  # Correlation Heatmap - obvious batch effect between platforms
  samples <- sampleInfo[order(sampleInfo$platform_id), ]
  gene_matrix <- arrays_combined[,rownames(samples)]
  samples <- subset(samples, select = platform_id)
  colnames(samples) <- 'Platform'
  
  corMatrix <- cor(gene_matrix, use = 'c')
  pheatmap(corMatrix,
           annotation_col = samples, fontsize = 7)
  
  
  # Sample clustering #
  # make binary vector for platform
  platform <- binarizeCategoricalVariable(sampleInfo$platform_id)
  colnames(platform) <- 'GPL570 (Black) vs\n\nGPL96 (Red)'
  rownames(platform) <- rownames(sampleInfo)
  plotClusterTreeSamples(datExpr=t(arrays_combined), y=platform, 
                         main = 'Sample Dendrogram with Platform Indicator', cex.traitLabels = 0.9)
  dev.off()
  
  
  ## Join the PCs to the sample information
  cbind(sampleInfo, pc$x) %>% 
    ggplot(aes(x = PC1, y=PC2, col=platform_id, label=sample_id, shape=dataset)) + 
    geom_point(size=3) + scale_color_manual(values = c('#00BFC4','#F8766D')) + geom_text_repel()
  
  ggsave(filename = paste0(output_dir,'/',filename,'_PCA.pdf'),height = 8, width=10)
  
}


# LIMMA DE ANALYSIS ####
# Function to perform limma differential expression analysis.
# Plots a density plot of median values.
# Prompts user for median cutoff point.
# Probes with a median expression value < cutoff in > x samples removed here. 
# X samples = number of samples in smallest experimental group.
# Annotates results with 'Gene_symbol' and 'ENTREZ_GENE_ID' columns.
# Returns annotated limma results table and saves to file. 

# array_data = pre-processed gene expression data (in this case the full merged array).
# metadata = pre-processed sample metadata.
# feature_data = pre-processed probe feature annotation data.
# design = limma design matrix.
# contrasts = object returned from makeContrasts function, detailing which variables are being compared. 
# pval = adjusted p-value cutoff.
# lfc = log fold change cut off.
# filter_dups = if set to True, probes mapping to the same gene will be filtered, with the most.
# significant probe being retained.
# OUTDIR = output directory.
# file_name = name tag for the output file.

run_limma <- function(array_data, metadata, feature_data, design, contrasts, 
                      pval = 1, lfc = 0, filter_dups = T, OUTDIR, file_name){
  
  # get cutoff from user
  repeat{
    
    probe_medians <- rowMedians(array_data)
    
    hist(probe_medians, 100, col = 'cornsilk1', freq = F, 
         main = 'Histogram of median intensities', 
         border = 'antiquewhite4',
         xlab = 'Median intensities')
    lines(density(probe_medians), col = 'black')
    
    cutoff = readline(prompt = 'Please select a cut-off value: ');
    
    abline(v = as.numeric(cutoff), lwd=2,col='red')
    
    happy <- readline(prompt = 'Are you satisfied with this cut-off? [y/n]')
    
    if(happy == 'y'){
      
      break
    }
  }
  
  print('Performing limma DE analysis...')
  
  # filter by selected expression cutoff
  is_expressed <- array_data > as.numeric(cutoff)
  
  # find minimum number of samples out of all experiment groups
  group_info <- group_by(metadata, disease_state) %>%
    tally()
  
  # select only probes that have this level of expression in at least this number of samples
  sample_cutoff <- min(group_info$n)
  keep <- rowSums(is_expressed) >= sample_cutoff
  expressed <- array_data[keep,]
  
  # limma linear model fit
  fit <- lmFit(expressed, design)
  
  fit2 <- contrasts.fit(fit, contrasts)
  
  eBayes_fit <- eBayes(fit2)
  
  print(table(decideTests(fit2)))
  
  res <- topTable(eBayes_fit, num = Inf, p.value = pval, lfc = lfc)
  
  # if running on one gene only
  if(nrow(res) ==1){
    
    return(res)
  }
  
  # merge annotation data
  res <- tibble::rownames_to_column(res,"ID")
   
  feat_data <- feature_data[,c('ID','Gene_symbol', 'ENTREZ_GENE_ID')]

  annotated_res <- merge(res, feat_data, by = 'ID')

  # order by FDR  
  annotated_res <- annotated_res[order(annotated_res$adj.P.Val), ]
  
  if(filter_dups == T){
    
    annotated_res <- filter_limma(annotated_res)
  }
  
  write.csv(annotated_res, file=paste0(OUTDIR,'/',file_name,'.csv'))
  saveRDS(annotated_res, file=paste0(OUTDIR,'/',file_name,'.rds'))
  
  return(annotated_res)
  
}


# Filter Limma Function
# Filters out duplicate gene symbols and takes the most significant result.

# results_table = annotated limma results table.

filter_limma <- function(results_table){
  
  results_table <- results_table[order(results_table$adj.P.Val),] 
  
  # duplicated genes (WITH NO DUPLICATED BLANKS)
  dups <- results_table[which(duplicated(results_table$ENTREZ_GENE_ID)),]$ID
  
  print(paste0('Identified ', length(dups), ' duplicated Gene ID(s).'))
  
  if(length(dups) !=0){
    
    print('Identifying most significant duplicates...')
    # ordered by p so just remove all duplicates
    results_table <- results_table[-which(results_table$ID %in% dups),]
  }
  
  return(results_table)
  
}


# Plot heatmap Function
# Used for plotting a heatmap of differential expression results.

# limma_results = results table from 'run_limma' function.
# array_data = gene expression array data.
# metadata = sample metadata.
# topN = number of genes to plot.
# status_col = column name or index of variable being compared in limma analysis.

# NOTE: if using directly after 'run_limma' function, note that the expression data 
# should be pre-processed and batch-adjusted prior to implementation (if necessary). 

plot_heatmap <- function(limma_results, array_data, metadata, topN = 25, status_col, file_tag, OUTDIR){
  
  ids_of_interest <- mutate(limma_results, Rank = 1:n()) %>% 
    filter(Rank <= topN) %>% 
    pull(ID)
  
  gene_names <- mutate(limma_results, Rank = 1:n()) %>% 
    filter(Rank <= topN) %>% 
    pull('Gene_symbol')
  
  # process to get relevant info, ordering of groups, etc
  gene_matrix <- array_data[ids_of_interest,]
  samples <- metadata[colnames(gene_matrix),]
  samples <- samples[order(samples[,status_col]),] # divide control/PD used to be sample$disease_state
  gene_matrix <- gene_matrix[,rownames(samples)]
  samples <- subset(samples, select = status_col) # used to be disease_state
  colnames(samples) <- 'Group'
  
  # plot and save as pdf
  pdf(paste0(OUTDIR, "/", file_tag, "_expression.pdf"), height=7, width=15)
  pheatmap(gene_matrix,
           labels_row = gene_names,
           annotation_col = samples,
           fontsize_col = 6,
           scale = 'row',
           main = paste0('Heatmap of expression levels for ', length(gene_names), ' selected DE genes of interest'),
           cluster_rows = F, cluster_cols =  F)
  dev.off()
  
}


# CORRELATION ANALYSIS ####
# Function for running Pearson's correlation analysis and simple linear regression analysis 
# between one specific gene and other genes of interest between control/disease groups.
# Plots expression heatmap, scatter plot and correlation heatmaps.

# DE_results = results table from 'run_limma' function.
# genes_of_interest = vector of all gene symbols, including 'main_gene', that will be used in the analysis.
# main_gene = the specific gene of interest that will be correlated with all other genes. 
# array_data = gene expression array data.
# pheno_data = sample metadata.
# file_tag = name tag for output file.
# output_dir = directory to which output files should be saved.
# status_col = column name or index of variable being compared in limma analysis.

# NOTE: text plotting paramaters were hard coded here and may need to be adjusted.

run_correlation <- function(DE_results, genes_of_interest, main_gene,
                            array_data, pheno_data, file_tag, output_dir, status_col){
  
  if (file.exists(output_dir)==F) {
    dir.create(output_dir)
  }
  
  # sample_id col 
  pheno_data$sample_id <- rownames(pheno_data)
  
  #status_col <- grep('control|Control|CTRL|ctrl', pheno_data)
  
  # control samples
  controls <- pheno_data[which(pheno_data[,status_col] %in% c('control','Control','CTRL','ctrl')),]
  
  # disease samples
  disease <- pheno_data[which(!pheno_data[,status_col] %in% c('control','Control','CTRL','ctrl')),]
  
  # get DE data and order by logFC
  of_interest <- DE_results[which(DE_results$Gene_symbol %in% genes_of_interest),] 
  of_interest <- of_interest[order(of_interest$logFC),]
  
  # get array expression data
  of_interest_exp <- array_data[of_interest$ID,]
  
  # plot expression heatmap and save as pdf
  plot_heatmap(of_interest, of_interest_exp, metadata = pheno_data, 
               topN = length(genes_of_interest), status_col = status_col, file_tag = file_tag, OUTDIR = output_dir)
  
  # SIMPLE LINEAR REGRESSION ANALYSIS
  # modify expression layout
  rownames(of_interest_exp) <- of_interest$Gene_symbol 
  of_interest_exp <- as.data.frame(t(of_interest_exp))
  
  # split to samples of interest
  control_exp <- of_interest_exp[controls$sample_id,]
  disease_exp <- of_interest_exp[disease$sample_id,]
  
  # SAVE AS PDF FOR EACH PLOT
  pdf(paste0(output_dir, "/", file_tag, "_regression.pdf"), height=12, width=14)
  par(mfrow = c(4,4))
  lm_list <- list()
  for(i in 1:ncol(of_interest_exp)){
    
    # Controls
    lm_list[[paste0(main_gene, 'vs', colnames(control_exp)[i],'_ctrl')]] <- lm(control_exp[,main_gene] ~ control_exp[,i])
    names(lm_list[[paste0(main_gene, 'vs', colnames(control_exp)[i],'_ctrl')]]$coefficients)[2] <- colnames(control_exp)[i]
    
    plot(control_exp[,i],control_exp[,main_gene], main=paste0(main_gene,' vs ', colnames(control_exp)[i], ' Expression - Controls'),
         pch=16, xlab = colnames(control_exp)[i], ylab = main_gene, col = 'blue') 
    abline(lm_list[[paste0(main_gene, 'vs', colnames(control_exp)[i],'_ctrl')]])
    
    control_cor <- cor.test(control_exp[,main_gene],control_exp[,i])
    corr <- control_cor$estimate
    pval <- control_cor$p.value
    
    text(paste("Cor:", round(corr,3)), x = 10.5, y = 6.5)
    text(paste("p-value:", round(pval,8)), x = 10.5, y= 6.1)
    
    # PD
    lm_list[[paste0(main_gene, 'vs', colnames(disease_exp)[i],'_disease')]] <- lm(disease_exp[,main_gene] ~ disease_exp[,i])
    names(lm_list[[paste0(main_gene, 'vs', colnames(disease_exp)[i],'_disease')]]$coefficients)[2] <- colnames(disease_exp)[i]
    
    plot(disease_exp[,i], disease_exp[,main_gene], main=paste0(main_gene,' vs ', colnames(disease_exp)[i], ' Expression - PD'),
         pch=16, xlab = colnames(disease_exp)[i], ylab = main_gene, col = 'red') 
    abline(lm_list[[paste0(main_gene, 'vs', colnames(disease_exp)[i],'_disease')]])
    
    disease_cor <- cor.test(disease_exp[,main_gene],disease_exp[,i])
    PD_corr <- disease_cor$estimate
    PD_pval <- disease_cor$p.value
    
    text(paste("Cor:", round(PD_corr,3)), x = 7.5, y = 10)
    text(paste("p-value:", round(PD_pval,8)), x = 7.5, y= 9.6)
    
  }
  
  dev.off()
  
  pdf(paste0(output_dir, "/", file_tag, "_cor_heatmap.pdf"))
  par(mfrow=c(1,1))
  pheatmap(cor(control_exp), main = 'Correlation Matrix Controls',
           cluster_rows = F, cluster_cols =  F)
  
  pheatmap(cor(disease_exp), main = 'Correlation Matrix PD',
           cluster_rows = F, cluster_cols =  F)
  dev.off()
  
  lm_list[['control_exp']] <- control_exp
  lm_list[['disease_exp']] <- disease_exp
  
  saveRDS(lm_list, file=(paste0(output_dir,'/',main_gene,'_lm_list.rds')))
  
  return(lm_list)
}


# GO OVER-REP ####
# Performs GO over-representation analysis for DEGs.
# CC, BP and MF ontologies. 

# genes_of_interest = vector containing DEGs of interest
# ref_list = all unique genes included in combined array
# keytype = format in which genes are identified 
# organism_db = organism annotation database
# FDR = False discovery rate cutoff

GO_overep <- function(genes_of_interest, ref_list, keytype = 'ENTREZID', organism_db = org.Hs.eg.db,
                      FDR = 0.05){
  
  GO_results <- list()
  for(ont in c('CC','BP', 'MF')){
    
    go_over <- enrichGO(gene = genes_of_interest,
                        universe = ref_list,
                        OrgDb = organism_db, 
                        keyType = keytype,
                        readable = T,
                        ont = ont)
    
    go_over@result <- go_over@result[which(go_over@result$p.adjust<FDR),]
    GO_results[[paste0(ont)]] <- go_over
    
  }
  
  return(GO_results)
}



# GSEA ####
# Performs gene set enrichment analysis for entire list of DEGs (significant and non-significant).
# Ranks all DEGs by t-value.
# Extracts gene set data from The Broad Institute of Molecular Signatures Database (MSigDB).

# full_limma_res = full results table from run-limma containing significant and non-significant results. 
# entrez_col = name of column containing Entrez ID numbers. 
# rank_col = name of column being used for ranking.
# species = name of species being analysed e.g. 'Homo Sapiens'.
# cat = MSigDB main gene set category code.
# subcategory = MSigDb gene set subcategory code.
# pval = adjusted p-value cutoff.

run_gsea <- function(full_limma_res, entrez_col, rank_col, species, cat, subcategory=NULL, pval=0.05){

  # ranking method taken from:
  # https://stephenturner.github.io/deseq-to-fgsea/#load_some_results_from_deseq2 
  # https://github.com/brucemoran/RNAseqon/blob/master/R/fgsea.R 
  
  # double exclamation mark selects just that column, single ! is complement of this column
  res_rank <- dplyr::select(.data = full_limma_res, !!entrez_col, !!rank_col) %>%
    na.omit() %>%
    dplyr::distinct() %>%
    dplyr::group_by(!!as.symbol(entrez_col)) %>%
    dplyr::summarize(rank = mean(!!as.symbol(rank_col)))
  
  rank_vec <- tibble::deframe(res_rank)
  rank_vec <- sort(rank_vec, decreasing = T)
  
  if(is.null(subcategory)==T){
    
    t2g <- msigdbr::msigdbr(species = species, category = cat) %>%
      dplyr::select(gs_name, entrez_gene)
    
  }else{
    
    t2g <- msigdbr::msigdbr(species = species, category = cat, subcategory = subcategory) %>% 
      dplyr::select(gs_name, entrez_gene)
  }
  
  enrich_res <- GSEA(rank_vec, TERM2GENE = t2g, eps = 0, pvalueCutoff = pval, seed=T)
  
  return(enrich_res)
  
}


# WGCNA ####

# Filter Outliers 
# Function for helping the user identify and remove outliers in expression data.
# Mainly code from the WGCNA tutorial documentation.
# Returns a filtered expression array.

# exprs_data = gene expression array data.
# output_dir = directory to which files should be saved.
# nametag = file name tag.

filter_outliers <- function(exprs_data,output_dir, nametag){
  
  # trees of clustered samples
  sampleTree = hclust(dist(exprs_data), method = "average")
  
  # get cutoff from user
  par(mfrow=c(1,1))
  par(mar = c(0, 4, 2, 0))
  plot(sampleTree, main = paste("Sample clustering on all genes"),
       xlab="", sub="", cex = 0.7)
  
  cutoff = readline(prompt = 'Please select a cut-off height value: ')
  
  abline(h = as.numeric(cutoff), col = 'red')
  
  remove <- readline(prompt = 'Remove all samples above cutoff? [y/n] ')
  
  if(remove == 'y'){
    
    # Find clusters cut by the line
    labels = cutreeStatic(sampleTree, cutHeight = as.numeric(cutoff))
    
    # Keep the largest one (labeled by the number 1)
    keep = (labels==1)
    exprs_data = exprs_data[keep, ]
    
    pdf(file=paste0(output_dir,'/',nametag,'.pdf'), width=12)
    par(mar = c(0, 4, 2, 0))
    plot(sampleTree, main = paste("Sample clustering on all genes"),
         xlab="", sub="", cex = 0.7)
    abline(h = as.numeric(cutoff), col='red')
    
    dev.off()
  }
  
  return(exprs_data)   
}


# Choose Soft Threshold 
# Function for helping the user choose soft thresholding power.
# Mainly code from the WGCNA tutorial documentation.
# Generates plots to aid in the assessment of scale free topology 

# exprs_data = gene expression array data.
# powers = vector of integers representing the range of powers to investigate.
# network_type = type of weighted network to use.
# output_dir = directory to which files should be saved.
# file_name = desired name of file to be saved


choose_sft <- function(exprs_data, powers, network_type="signed", output_dir, file_name){
  
  # Call the network topology analysis function
  sft = pickSoftThreshold(
    exprs_data,           
    networkType = network_type,
    powerVector = powers,
    verbose = 5
  )
  
  pdf(file=paste0(output_dir, '/', file_name, '.pdf'), width=12)
  par(mfrow = c(1,2))
  cex1 = 0.9
  
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2",
       main = paste("Scale independence")
  )
  
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = cex1, col = "red"
  )
  
  abline(h = 0.9, col = "red")
  plot(sft$fitIndices[, 1],
       sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       type = "n",
       main = paste("Mean connectivity")
  )
  text(sft$fitIndices[, 1],
       sft$fitIndices[, 5],
       labels = powers,
       cex = cex1, col = "red")
  
  dev.off()
  
  return(sft)
}


# Module Details 
# Function for calculating MEs, Module Membership (Signed Eigengene Connectivity), Hub genes
# Returns list containing MEs, module genes, signedKME for all genes with all modules,
# signedKME of each gene for its associated module and hub genes.

# exprs_data = gene expression array data.
# net = network object returned from blockwiseModules function
# metadata = sample metadata


module_details <- function(exprs_data, net, metadata){
  
  # extract gene module colours
  moduleColors <- labels2colors(full_net$colors)
  
  
  # OVERALL MODULE MEMBERSHIP #
  
  # Calculate Module Eigengenes (1st principal component) of modules
  MEs0 = moduleEigengenes(exprs_data, moduleColors)$eigengenes
  
  # Reorder (eigen-)vectors such that similar ones are next to each other.
  MEs = orderMEs(MEs0)
  
  # Full gene set in each module
  module_genes <- list()
  intModules <- unique(moduleColors)
  for (module in intModules){
    # Select module probes
    modGenes <- (moduleColors==module)
    # Get their entrez ID codes
    mod_IDs <- anno[modGenes, c('Gene_symbol','ENTREZ_GENE_ID')]
    
    module_genes[[module]] <- mod_IDs
  }
  
  geneModuleMembership <- signedKME(exprs_data, MEs)
  # ensure correct order so the KMEtable is correct below
  colors <- substring(colnames(MEs),3)
  colnames(geneModuleMembership)=paste("PC",colors,".cor",sep="")
  geneModuleMembership <- geneModuleMembership[,order(colnames(geneModuleMembership))]
  
  # ensure correct order as before
  MMPvalue <- corPvalueStudent(as.matrix(geneModuleMembership),dim(exprs_data)[[2]])
  colnames(MMPvalue)=paste("PC",colors,".pval",sep="")
  MMPvalue <- MMPvalue[,order(colnames(MMPvalue))]
  
  # create kME table
  probe <- colnames(exprs_data)
  Gene <- anno[probe, 'Gene_symbol'] 
  kMEtable = cbind(probe,Gene,moduleColors)
  
  for (i in 1:length(colors)){
    kMEtable = cbind(kMEtable, geneModuleMembership[,i], MMPvalue[,i])
  }
  colnames(kMEtable) <- c("ProbeID","Gene","Module",sort(c(colnames(geneModuleMembership),
                                                           colnames(MMPvalue))))
  kMEtable <- as.data.frame(kMEtable)
  
 
  # INDIVIDUAL MODULE MEMBERSHIP #
  # adapted from https://github.com/Bioinformatics-rookie/WGCNA 
  sign_eig_connect <- list()
  
  for(module in substring(colnames(MEs),3)){
    
    if(module == "grey") next  # grey are unclustered
    
    ME=as.data.frame(MEs[,paste("ME",module,sep="")])
    colnames(ME)=module
    
    datModExpr=exprs_data[,moduleColors==module] # gene exprs for that module gene set
    datKME = signedKME(datModExpr, ME, outputColumnName =paste0('KME',substring(module,1,2)))
    
    datKME=cbind(datKME, color=rep(module,length(datKME)),Gene_symbol=anno[rownames(datKME), 'Gene_symbol'])
    
    sign_eig_connect[[module]] <- datKME
  }
  
  # HUB GENES #
  # extract top 10 signedKME values
  module_hubs <- lapply(sign_eig_connect, function(x) x[order(-x[,1])[1:10],3])
  
  # extract top 10 pvals from overall KME table
  signedKME_pvals <- list()
  for(module in names(sign_eig_connect)){
    
    KME_col <- grep(paste0('PC',module,'.cor'), colnames(kMEtable))
    pval_col <- grep(paste0('PC',module,'.pval'), colnames(kMEtable))
    
    genes <- kMEtable[which(kMEtable$Module == module),]
    
    ordered <- genes[order(genes[,KME_col], decreasing = T)[1:10],pval_col]
    
    signedKME_pvals[[paste0(module,'.pval')]] <- ordered
    
  }
  
  # bind gene symbols and associated p-values side-by-side
  top_hubs <- do.call(cbind, module_hubs)
  hub_pval <- do.call(cbind, signedKME_pvals)
  hub_info <- cbind(top_hubs, hub_pval)
  hub_info <- hub_info[,order(colnames(hub_info))]
  
  module_info <- list(module_genes = module_genes,
                      MEs = MEs,
                      kMEtable = as.data.frame(kMEtable),
                      moduleKME = sign_eig_connect,
                      hub_genes = as.data.frame(hub_info))
  
  return(module_info)
  
}
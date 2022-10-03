# MSc Bioinformatics & Computational Biology 
# University College Cork
# Name: Conor Giles-Doran 
# Submitted: October 2022


# SCRIPT FOR MANUALLY PROCESSING PHENOTYPE DATA #
# Sourced in 'WGCNA_workflow.R' and 'DE_GO_GSEA_workflow.R' scripts #

# METADATA ASSEMBLY MUST BE DONE MANUALLY #

# read in individually processed data 
processed_list <- readRDS(file = paste0(OUTDIR, '/Processed/full_processed_list.rds'))

pheno <- lapply(processed_list[[2]][[1]],pData)

# modified phenotype data - all performed manually as column names/content differ
mod_pheno <- list()

# extract relevant metadata for each dataset, in the correct order of columns
mod_pheno[['GSE20186']] <- pheno[['GSE20186']][,c('platform_id', 'characteristics_ch1.2', 
                                                  'tissue:ch1','age:ch1','gender:ch1')]

mod_pheno[['GSE20292']] <- pheno[['GSE20292']][,c('platform_id', 'disease state:ch1', 
                                                  'brain region:ch1','age:ch1','gender:ch1')]

mod_pheno[['GSE49036']] <- pheno[['GSE49036']][,c('platform_id','disease state:ch1', 'tissue:ch1')]
mod_pheno[['GSE49036']]$age <- NA 
mod_pheno[['GSE49036']]$sex <- NA


mod_pheno[["GSE7621"]] <- pheno[["GSE7621"]][,c('platform_id','characteristics_ch1','source_name_ch1')]
mod_pheno[["GSE7621"]]$age <- NA
mod_pheno[["GSE7621"]]$sex <- pheno[["GSE7621"]]$characteristics_ch1.1


mod_pheno[["GSE8397"]] <- pheno[["GSE8397"]][,c('platform_id','title','source_name_ch1',
                                                'age:ch1')]
mod_pheno[["GSE8397"]][c('age', 'sex')] <- str_split_fixed(mod_pheno[["GSE8397"]]$`age:ch1`, ';', 2)
mod_pheno[["GSE8397"]] <- mod_pheno[["GSE8397"]][,-4]

# add dataset and sample_id columns to each and rename columns
for(i in 1:length(mod_pheno)){
  
  mod_pheno[[i]]$dataset <- names(mod_pheno[i])
  mod_pheno[[i]]$sample_id <- rownames(mod_pheno[[i]])
  colnames(mod_pheno[[i]]) <- c('platform_id', 'disease_state', 'tissue', 'age','sex', 'dataset', 'sample_id')
}

# ensure columns are in a tidy order
mod_pheno <- lapply(mod_pheno, function(x) dplyr::select(.data = x, !!'sample_id', !!'dataset', !!'platform_id', 
                                                         !!'tissue', !!'disease_state', !!'sex', !!'age'))

# combine metadata
sample_names <- unlist(lapply(mod_pheno, rownames))
sample_names <- unlist(sample_names)
pheno_combined <- do.call(rbind, mod_pheno)
rownames(pheno_combined) <- sample_names

# make all same format in combined table

# control
control_ind <- c(grep('control', pheno_combined$disease_state),
                 grep('Control', pheno_combined$disease_state), 
                 grep('normal', pheno_combined$disease_state))

pheno_combined$disease_state[control_ind] <- 'control'

# affected
affected_ind <- c(grep('inson', pheno_combined$disease_state))
pheno_combined$disease_state[affected_ind] <- "Parkinson's Disease"

# get rid of text N/As
n_a <- grep('N/A', pheno_combined$age)

pheno_combined$age[n_a] <- NA

#pheno_tally <- group_by(pheno_combined, dataset,disease_state) %>% tally()

# format sex
F_change <- grep('female|Female', pheno_combined$sex)

pheno_combined$sex[F_change] <- 'F'

M_change <- grep('male|Male', pheno_combined$sex)

pheno_combined$sex[M_change] <- 'M'

pheno_combined$sex <- gsub(".* ", '', pheno_combined$sex)

# Lewy body samples must be removed 
lewy_body <- grep('body', pheno_combined$disease_state) # indices

# Tissues other than SN must be removed
not_SN <- grep('substantia nigra', pheno_combined[,'tissue'], invert = T)

cleaned_pheno <- pheno_combined[-c(lewy_body, not_SN),]

rm(mod_pheno,pheno, pheno_combined)

saveRDS(cleaned_pheno, file = paste0(OUTDIR, '/Processed/cleaned_pheno.rds'))


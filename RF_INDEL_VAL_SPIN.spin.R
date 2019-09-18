#'---
#' title: "Random Forest Prediciton of Somatic Indels from Marek's Disease Lymphomas"
#' author: Alec Steep
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
#' output: 
#'     html_document:
#'         code_folding: hide
#'         toc: true
#'         highlight: zenburn
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)

#' ## Explanation of Samples and Data 
#' ### Explaination of Biological Samples and Sequencing
#' A random forest model is used to predict somatic indels in Marek's Disease Lymphomas seeded in gonadal tissue. Twenty-two gonadal tumors and four biological replicates from four tumors underwent whole genome sequencing analysis. Four indel callers were used in tandem on WGS data. A similar machine learning model (somaticseq) was trained with in-silico indels and predicitons were used as an additional feature titled "SOMATIC". In-silico indels were generated with BAMSurgeon (not biological but created by us artifically in bam files). High-confidence indels underwent targetted deep sequencing in 2 rounds to build a high-confidence trianing set. Validated indels were used to train our random forst model on raw indels calls from the original wgs data (4 callers) to predict true indels. Finally, high-confidence predicted and validated indels were combined to represent a final set of indels in Marek's Disease Lymphomas.
#'

#' ## Setup the Environment

#+ Setup Environment

################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# For information about the 'feather' file format
# https://blog.rstudio.com/2016/03/29/feather/

# Set the working directory
#setwd('/Users/Alec/Documents/Bioinformatics/MDV_Project/random_forest_indels')

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("R.utils")


# Load dependencies
pacs...man <- c("dplyr", "feather", "caret", "Metrics", "randomForest", "mlbench","e1071","tibble","magrittr","stringr", "data.table","ROCR","SomaticSignatures","GenomicRanges","pROC","readr","R.utils")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) 
})

############################################################
##### Functions ############################################
############################################################

# Make the 'not in' operator
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}

# Convert charcter columns with '; delimited numerical values to max numerical value
###########################################################
char2maxnum_df <- function(df, column) {
        for ( i in seq_along(df[[column]]) ) {
                # Get the value as a string
                value = as.character(unlist(df[[column]][i]))
                if (!is.na(value)) {
                        vals = strsplit(value, ";")[[1]]
                        new_val = max(as.numeric(vals))
                        df[i,column] = new_val
                }
        }
        
        df
}
###########################################################

# Capture the Date
date <- format.Date( Sys.Date(), '%Y%m%d' )
auth <- "steep"

################################################################################

#' ## Load Data
#' ##### Data Files to Load:
#' * Sample Annotation File
#' * Raw Variant Calls (Indels called by 2 or more callers)
#' * Germline Indels (Called from WGS data in same 22 tumors)
#'
#+ Load the Data

################################################################################
#####     Load Data for Random Forest      #####################################
################################################################################

# Load in the sample annotation table
sample_df <- read.table(file = './data/sample_annotation.txt', header = TRUE, sep = '\t')

# Load the dataframe of all variant calls
gunzip('./data/bigboy_snv_indel_processed_ml.feather.gz')
feather_file <- "./data/bigboy_snv_indel_processed_ml.feather"
# Copy text file saved in ./data/raw_calls_snv_indel_process_ml.txt
df <- read_feather(feather_file)
gzip('./data/bigboy_snv_indel_processed_ml.feather')
#gunzip('./data/raw_calls_snv_indel_process_ml.txt.gz')
# This will take a fair amount of time, many data points. Also, don't mind the warning: "Unknown or uninitialised column: 'STATUS'."
#df <- read.table(file = './data/raw_calls_snv_indel_process_ml.txt', header = TRUE, sep = '\t')
dim(df)
# Gzip the file again
#gzip('./data/raw_calls_snv_indel_process_ml.txt')
# We will only be applying machine learning to INDELS called by 2 or more callers (all true positives are called by 2 or more callers)
df <- filter(df, NUM_TOOLS >= 2 & VAR_TYPE == 'INDEL')
stats_2plusindels <- dim(df)[1]

# Load germline indels
gunzip('./data/F1_indels_lenient_from_L6_L7.txt.gz')
F1_INDELS <- read.table('./data/F1_indels_lenient_from_L6_L7.txt', sep = '\t', header = TRUE)
F1_INDELS <- as_tibble(F1_INDELS)
F1_INDELS <- F1_INDELS %>% dplyr::select(-FILTER)
stats_germline_indels <- dim(F1_INDELS)[1]
# Gzip the file again
gzip('./data/F1_indels_lenient_from_L6_L7.txt')

# Load the validated calls from val1
dfval1 <- read.table('./data/somatic_snvs_and_indels_validated_agriplex_2014_2017.txt', 
                     header = FALSE, sep ='\t', comment = "#")
names(dfval1) <- c('CHROM','POS','REF','ALT','VAR_ID','MUT','IMPACT','SYMBOL',
                   'GENE_ID','ORTHOLOGUE','TSN_VAR','SAMPLE','VAC','VAF','NUM_TOOLS',
                   'CGC_STATUS','FILTER')

# Load the somatic variants from the 2nd round of validation (file comes from Ipython notebook script)
infile <- "./data/somatic_INDELS_val2_VAFs.txt"
dfval2 <- read.table(infile, sep = '\t', header = TRUE)

# Load power table
power_stats <- read.table('./data/power_table.txt',sep = '\t', header = TRUE)
names(power_stats) <- str_replace(names(power_stats), 'sam', 'SAMPLE')
names(power_stats) <- str_replace(names(power_stats), 'power', 'POWER')

################################################################################

#' ##### Examine Data:
#' Raw Indels called by 2 or more callers
head(df)
#' High confidence germline indels
head(F1_INDELS)
#' Validated somatic Indels from 1nd round of validation
head(dfval1)
#' Validated somatic Indels from 2nd round of validation
head(dfval2)
#' Feature to add: The relative power to detect somatic variants per sample (high, medium, low) ~ tumor purity
print(summary(power_stats))

#' ## Remove germline indels from somatic indel call cohort

################################################################################
####### Remove Germline INDELS from Raw Somatic INDEL Calls ####################
################################################################################

# Collect the unique somatic variants
som_vars <- df %>% dplyr::select(CHROM,POS,REF,ALT) %>% unique
som_vars <- as_tibble(som_vars)

# Combine the germline and somatic variants
germ_som <- rbind(F1_INDELS,som_vars)
germ_som <- as_tibble(germ_som)

# Extract any somatic variants that are actually germline varinats
# Takes a minute or 2
duplicated_som <- germ_som[duplicated(germ_som),]

# Mark germline variants in somatic cohort
duplicated_som$STATUS <- 'GERMLINE'

# Perform a left join, we only add the 'VAL2' column
df <- df %>% left_join(duplicated_som, by=c("CHROM","POS","REF","ALT"))

# Adjust the status name
df$STATUS[is.na(df$STATUS)] <- 'SOMATIC'

# Remove germline samples
df <- filter(df, STATUS != 'GERMLINE')

# Drop the unneccessary val columns
df <- df %>% dplyr::select(-STATUS)
df <- as_tibble(df)
# Stat: Raw Indels (2 or more callers) after germline variant removed
som_post_germ <- dim(df)[1]

# Write the output file
#outfile = './data/somatic_indels_2_callers.txt'
#write.table(df, file = outfile, quote = FALSE, sep = '\t', row.names = FALSE)

#' ##### Stats:
#' Raw Indels called by 2 or more callers
print(stats_2plusindels)
#' High confidence germline indels
print(stats_germline_indels)
#' Raw Indels (2 or more callers) after germline variant removed
print(som_post_germ)

################################################################################

#' ## Clean & Process Validated Somatic Indels (Validated from targetted deep sequencing)
#' 
#' #' Clean & Process Indels from Validation Round 2 
#' (Process Round 2 before Round 1 for scripting convenience)

################################################################################
####### Process Validation Round 2 Calls  ######################################
################################################################################

# Set all VALs to 0
dfval2$VAL <-  0
#unique(dfval2$SAMPLE)

# Adjust all TRUE positives
# IL1RAPL1
dfval2[dfval2$SYMBOL == 'IL1RAPL1' & dfval2$SAMPLE %in% c('017834-2_2', '017834-2'),]$VAF <- 0.2251553
dfval2[dfval2$SYMBOL == 'IL1RAPL1' & dfval2$SAMPLE %in% c('017842-2_2', '017842-2'),]$VAL <- 1
dfval2[dfval2$SYMBOL == 'IL1RAPL1' & dfval2$SAMPLE %in% c('017834-2_2', '017834-2'),]$VAL <- 1
# WNT5A
dfval2[dfval2$SYMBOL == 'WNT5A' & dfval2$SAMPLE %in% c('017911-1_2', '017911-1'),]$VAL <- 1
# DGKB
dfval2[dfval2$SYMBOL == 'DGKB' & dfval2$SAMPLE %in% c('017756-1','017756-2','017756-3','017766-1','017834-2_2',
                                                      '017834-2','017911-1','017911-1_2','017911-2'),]$VAL <- 0
dfval2[dfval2$SYMBOL == 'DGKB' & dfval2$SAMPLE %in% c('017911-1','017911-1_2','017911-2'),]$VAL <- 1
# RARB
dfval2[dfval2$SYMBOL == 'RARB' & dfval2$SAMPLE %in% c('017842-2','017842-2_2'),]$VAL <- 1
# ITGA9
dfval2[dfval2$SYMBOL == 'ITGA9' & dfval2$SAMPLE %in% c('017911-1','017911-1_2','017911-2','017911-3'),]$VAL <- 1
# HECW1
dfval2[dfval2$SYMBOL == 'HECW1' & dfval2$SAMPLE %in% c('017834-2','017834-2_2'),]$VAL <- 1

# Remove variants that were not properly measured
dfval2 <- dfval2[!is.na(dfval2$VAF),]
# IKZF1 INDELS were not properly measured
dfval2 <- dfval2[dfval2$SYMBOL != 'IKZF1',]
#dim(dfval2)
# Label variants that failed validation at 'FP'
dfval2$VAL2 <- 'EMPTY'
dfval2[dfval2$VAL == 0,]$VAL2 <- 'FP'
dfval2[dfval2$VAL == 1,]$VAL2 <- 'TP'

# Add annotation if sample is germline or not
dfval2 <- mutate(dfval2, TUM_NORM = ifelse(str_detect(SAMPLE,"-0"), 'NORMAL', 'TUMOR'))
# Do germline samples exhbit VAF greater than 0.05?
#filter(dfval2, TUM_NORM == 'NORMAL' & VAF >= 0.05) %>% dplyr::select(SYMBOL, VAF, SAMPLE, VAL2)
# Note: Suggests that certain control samples contained contaminating neoplastic cells or that limitations of bioinformatic analysis may contribute to these measurements

# Extract somatic variants 
dfval2 <- filter(dfval2, TUM_NORM != 'NORMAL')
# Examine mutated gene frequencies
dfval2$SAMPLE <- as.character(dfval2$SAMPLE)
# Stats: Genes mutated with validated variants 
stats_genes_val_mut2 <- sort(table(filter(dfval2, VAL2 == 'TP') %>% dplyr::select(SYMBOL)))
# Stats: Samples mutated with validated variants 
stats_sam_val_mut2 <- sort(table(filter(dfval2, VAL2 == 'TP') %>% dplyr::select(SAMPLE)))

# Remove duplicated values
dfval2 <- dplyr::distinct(dfval2, CHROM,POS,REF,ALT,SAMPLE, .keep_all= TRUE)
dfval2 <- dfval2 %>% dplyr::select(CHROM,POS,REF,ALT,SAMPLE,SYMBOL,NUM_TOOLS,VAR_TYPE,CON,IMPACT,SYMBOL.1,GENE,HUGO,HOTSPOT_V2,VAL,VAF,TUM_NORM,VAL2)
# Take only the validated indels
dfval2_out <- filter(dfval2, VAL2 == 'TP')

# Save validated indel calls (both FP and TP) from round 2
#write.table(dfval2, file = './data/validated_somatic_indels_round2.txt', quote = FALSE, sep = '\t', row.names = FALSE)

################################################################################

#' Clean & Process Indels from Validation Round 1

################################################################################
####### Clean & Process Indels from Validation Round 1  ########################
################################################################################

# Annotate and seperate INDELs
nucs <- c('A','T','C','G')
dfval1 <- dfval1 %>% mutate(VAR_TYPE = ifelse(ALT %in% nucs & REF %in% nucs,'SNV','INDEL'))
dfval1 <- filter(dfval1, VAR_TYPE == 'INDEL')

# Splits the column by a delimiter and replicates rows (SPLIT_COL)
dfval1 <- dfval1 %>% mutate(SAMPLE = strsplit(as.character(SAMPLE), ";")) %>% 
        tidyr::unnest(SAMPLE)

# Select columns of interest
dfval1 <- dfval1 %>% dplyr::select('CHROM','POS','REF','ALT','SAMPLE')

# Add Validation annotation
dfval1$VAL2 <- 'TP'

################################################################################

#' Further Process then Combine Validation Data Sets (Both Rounds)

################################################################################
####### Further Process then Combine Validation Data Sets  #####################
################################################################################

# Subset val2 dataset to have matching names
dfval2 <- dfval2[,(names(dfval1)[names(dfval1) %in% names(dfval2)])]
# Adjust the sample names of validated samples from second round of validation
dfval2$SAMPLE <- as.character(dfval2$SAMPLE)
dfval1$SAMPLE <- as.character(dfval1$SAMPLE)
# Rename the samples (dfval2)
for (x in c(1:26)) {
        sample_old <- as.character(sample_df[x,'SAMPLES'])
        sample_new <- as.character(sample_df[x,'ALL'])
        #print(sample_old)
        if (sample_old %in% dfval2$SAMPLE){
                dfval2[dfval2$SAMPLE == sample_old,]$SAMPLE <- sample_new
        }
}

# Adjust the names of a few samples individually
dfval1[dfval1$SAMPLE == "911-1-2",]$SAMPLE <- '911-1_2_S13'
dfval1[dfval1$SAMPLE == "834-2-2",]$SAMPLE <- '834-2_2_S12'
dfval1[dfval1$SAMPLE == "901-2-2",]$SAMPLE <- "901-2_2_S26"

# Rename the samples (dfval1)
for (x in c(1:26)) {
        #print(x)
        sample_old <- as.character(sample_df[x,'SHORT'])
        sample_new <- as.character(sample_df[x,'ALL'])
        #print(sample_old)
        #print(sample_new)
        if (sample_old %in% dfval1$SAMPLE){
                print('Replace')
                dfval1[dfval1$SAMPLE == sample_old,]$SAMPLE <- sample_new
        }
}
# Filter samples that were not in our original 22 tumor cohort, which we will make predictions on
#dfval1$SAMPLE[dfval1$SAMPLE %!in% sample_df$ALL]
#table(dfval2$SAMPLE[dfval2$SAMPLE %!in% sample_df$ALL])
dfval1 <- filter(dfval1, SAMPLE %in% unique(sample_df$ALL))
dfval2 <- filter(dfval2, SAMPLE %in% unique(sample_df$ALL))

# Grab only the necessary columns
dfval2 <- dfval2 %>% dplyr::select('CHROM','POS','REF','ALT','SAMPLE','VAL2')

# Combine validation results
val_df <- rbind(dfval1,dfval2)
#table(duplicated(dfval2))
val_df <- as_tibble(val_df)
val_df$CHROM <- as.character(val_df$CHROM)
#dim(val_df)

# Save the validated INDELs from both rounds
#write.table(val_df, file = './data/validated_somatic_indels_both_rounds.txt', quote = FALSE, sep = '\t', row.names = FALSE)

################################################################################

#' Validated Variants Combined
print(table(val_df$VAL2))
print(val_df %>% dplyr::filter(VAL2 == 'TP'))

################################################################################

#' Combine the Validated Indel Annotation with Raw Calls

################################################################################
####### Combine the Validated Indel Annotation with Raw Calls ##################
################################################################################

# Ensure all samples are in sync between validated calls and raw calls
#unique(val_df$SAMPLE)[unique(val_df$SAMPLE) %!in% unique(df$SAMPLE)]

# Perform a left join, we only add the 'VAL2' column
df <- df %>% left_join(val_df, by=c("CHROM","POS","REF","ALT","SAMPLE"))
#'VAL2' %in% names(df)
# Apply this to the validation column
#table(val_df$VAL2)
df$VALIDATION <- df$VAL2
# Drop the unneccessary val columns
df <- df %>% dplyr::select(-VAL2)
# Check how many samples this applies to
#as.vector(unique(df$SAMPLE))
# Dedup the df dataframe
#table(duplicated(df))
df <- dplyr::distinct(df)
#table(duplicated(df))

################################################################################

#' Validated Variants in Raw Calls
print(table(df$VALIDATION))
print(head( df %>% dplyr::filter(VALIDATION == 'TP') %>% dplyr::select(CHROM,POS,REF,ALT,SAMPLE,VALIDATION,MVJSDU,SYMBOL) ) )

#' Add Features to Dataset
#' Added Features:
#' * Column for each caller
#' * The power to detect a single nucleotide variant per sample based on average coverage
#' 

################################################################################
####### FEATURE ADDITION  ######################################################
################################################################################

# Add a column for each caller
#table(df$MVJSDU)
# Add the Mutect column
m_calls <- c()
v_calls <- c()
j_calls <- c()
s_calls <- c()
for (call in names(table(df$MVJSDU))) {
        # M
        if (substr(call, start = 1, stop = 1) == "1") {
                m_calls <- append(m_calls,call)
        }
        # V
        if (substr(call, start = 3, stop = 3) == "1") {
                v_calls <- append(v_calls,call)
        }
        # J
        if (substr(call, start = 5, stop = 5) == "1") {
                j_calls <- append(j_calls,call)
        }
        # S
        if (substr(call, start = 7, stop = 7) == "1") {
                s_calls <- append(s_calls,call)
        }
}

df <- df %>% mutate(M = ifelse(MVJSDU %in% m_calls, 1, 0))
df <- df %>% mutate(V = ifelse(MVJSDU %in% v_calls, 1, 0))
df <- df %>% mutate(J = ifelse(MVJSDU %in% j_calls, 1, 0))
df <- df %>% mutate(S = ifelse(MVJSDU %in% s_calls, 1, 0))

# Adjust to factors
df$M <- as.factor(df$M)
df$V <- as.factor(df$V)
df$J <- as.factor(df$J)
df$S <- as.factor(df$S)

# Add power of each sample
# Adjust names
power_stats$SAMPLE <- as.character(power_stats$SAMPLE)
for (x in c(1:26)) {
        sample_old <- as.character(sample_df[x,'SAMPLES'])
        sample_new <- as.character(sample_df[x,'ALL'])
        #print(sample_old)
        if (sample_old %in% power_stats$SAMPLE){
                #print(sample_new)
                power_stats[power_stats$SAMPLE == sample_old,]$SAMPLE <- sample_new
        }
}

# Join the power df to main df
df <- left_join(df, power_stats, by="SAMPLE")
# Remove NA's
df$POWER <- as.character(df$POWER)
df[is.na(df$POWER),]$POWER <- 'other'
df$POWER <- as.factor(df$POWER)

################################################################################

################################################################################
#### Exploratory Analysis of Features; Most Influenctial FEATURES (Top 20) #####
################################################################################

# Remove duplicates from raw data
#dim(df)
df <- dplyr::distinct(df, CHROM,POS,REF,ALT,SAMPLE, .keep_all= TRUE)
# Save the original
df_original <- df
#dim(df)

# Round 1
# 1) AC_T. Resembles guassian dist before 15. Set a cutoff for 15.
summary(df$AD_T)
hist(df$AD_T, xlim=c(0,20), breaks=100)
#df[df$AD_T > 15,]$AD_T <- 15
hist(df$AC_T, xlim=c(0,3), breaks=100)
# 2) AD_T. Is exactly the same as AC_T, remove.
#df <- df %>% dplyr::select(-AD_T)
# 3) NUM_TOOLS.
df$NUM_TOOLS <- ordered(df$NUM_TOOLS)
# 4) SAMPLE. Leave as is for now.
# 5) SOMATIC. Adjust to numeric as an experiment
df$SOMATIC <- as.character(df$SOMATIC)
#df[df$SOMATIC == 'SOMATIC',]$SOMATIC <- '1'
#df[df$SOMATIC == 'NONSOMATIC',]$SOMATIC <- '0'

# Round 2
# 6) TNCF_TSN. 49845 NA's, remove. 
df <- df %>% dplyr::select(-TNCF_TSN)
# 7) CON. 22 levels, 320 NA's. Label these values as intergenic
df[is.na(df$CON),]$CON <- 'intergenic_variant'
# 8) AN. AN is the number of alleles in the entire cohort. (320 NAs)
# AN = AC * total samples with alt allele. AC is also NA in all cases. Assume 2.
df[is.na(df$AN),]$AN <- 2 * (df[is.na(df$AN),]$SAM_NUM)

# 9) SAM_NUM. 
# 10) AC. 320 NA's. Assign all of them to 2.
# AC is the number of alleles within a sample. Let's just assume it's 2 for 51 NA values
df[is.na(df$AC),]$AC <- 2
#df[df$AC > 10,]$AC <- 10
hist(df$AC, xlim=c(0,10), breaks=10)
# Reduce numeric to max of 1
hist(df$AC, xlim=c(0,10), breaks=20)
# Round 3
# 11) MVJSDU. Convert all 5 caller factors to one facot to reduce levels to 53 or lower
df$MVJSDU <- as.factor(df$MVJSDU)
# 12) S. Leave as is.
# 13) IMPACT. 320 NA's. These values have been replaced with intergenic--MODIFIER
df[is.na(df$IMPACT),]$IMPACT <- 'MODIFIER'
# 14) ALTERATION. Leave as is.
# 15) CONTEXT. Leave as is.
# 16) CpG. Leave as is.
# 17) POWER. Leave as is.

################################################################################
####### Adjust Features and DataSet to Comply with Random Forest Algorithm #####
################################################################################

## Formating
# Rename TUM_TYPES
names(df)[names(df) == 'CGC_TUMTYPES-S'] <- 'CGC_TUMTYPES_S'
names(df)[names(df) == 'CGC_TUMTYPES-G'] <- 'CGC_TUMTYPES_G'
#summary(is.na(df))
# consider: 'CGC_TUMTYPES_G','CGC_SYNDROME','CGC_TISSUE','CGC_MOLGEN','CGC_ROLE','CGC_MUT','CGC'
# Remove predictors that are deemed irrelevent or redundent
lose <-c('ALLELE')
df <- df[ , !(names(df) %in% lose)]
# Predictors with too many NA's that cannot be compensated for
lose <- c('ID','altMQ_N','altBQ_N','altNM_N','zMQ_N','zBQ_N','SYMBOL','GENE','FEAT_TYPE',
          'FEAT','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_POS','CDS_POS','PRO_POS',
          'AA','CODONS','DIST','STRAND','SYMBOL_SOURCE','HGNC_ID','ENSP','SWISSPROT','TREMBL',
          'UNIPARC','SIFT','DOMAINS','HGVS_OFFSET','HUGO','CGC_TUMTYPES_S','CGC_MUT','TNCF_TSN',
          'TNCF_SAM','TNCF_VAC','TNCF_VAF','TNCF',
          "GT_N","CON","CGC_TUMTYPES_G","CGC_SYNDROME")
#df <- df[ , !(names(df) %in% lose)]

# Compensate for NA values for categorical variables
# CHROM has amny chromosome, reduce levels to major chromosomes
#df$CHROM <- as.character(df$CHROM)
#df[grepl("NT_", df$CHROM),]$CHROM <- 'Micro'
#df$CHROM <- as.factor(df$CHROM)

# Adjust NA values in high throughput manner when possible
easy_adjust <- c('CGC_TUMTYPES_G','CGC_SYNDROME','CGC_TISSUE','CGC_MOLGEN','CGC_ROLE','CGC')
#easy_adjust <- c('CGC_SYNDROME','CGC_TISSUE','CGC_MOLGEN','CGC_ROLE','CGC')

#summary(is.na(df))
for (p in easy_adjust) {
        df[[p]] <- as.character(df[[p]])
        df[is.na(df[[p]]),][[p]] <- 'NO'
        df[[p]] <- as.factor(df[[p]])
}

# Compensate for NA values for numerical variables

# Mapping quality variables (refMQ_N, altMQ_N,altMQ_T,refMQ_T)
# refMQ_N has 222 NA values, substitute with the mean reference allele MQ from tumor
# HC samples: refMQ_N ~ refMQ_T <- SAFE, low NA freq
df[is.na(df$refMQ_N),]$refMQ_N <- as.numeric(summary(df[is.na(df$refMQ_N),]$refMQ_T)[3])
# altMQ_T has 4 NA values. Replace values with altMQ_N.
#df[is.na(df$altMQ_T),]$altMQ_T <- df[is.na(df$altMQ_T),]$altMQ_N
# refMQ_T has 886 NA values. Replace values with refMQ_N
df[is.na(df$refMQ_T),]$refMQ_T <- df[is.na(df$refMQ_T),]$refMQ_N
# Base quality variables
# refBQ_N has 211 NA values. 
# Step 1. Replace refBQ_N NA values with refBQ_T (147 NA values remianing)
df[is.na(df$refBQ_N),]$refBQ_N <- df[is.na(df$refBQ_N),]$refBQ_T
# Step 2. Replace refBQ_N NA values with altBQ_T (maintains the MQ unique to locus)
df[is.na(df$refBQ_N),]$refBQ_N <- df[is.na(df$refBQ_N),]$altBQ_T
# refMQ_T has NA when VAF_T is 0. Replace those 51 values with refMQ_N
#df[is.na(df$refMQ_T),]$refMQ_T <- df[is.na(df$refMQ_T),]$refMQ_N
# refBQ_T has NA when there are no reads supporting the reference allele, replace values with BQ of alt allele
df[is.na(df$refBQ_T),]$refBQ_T <- df[is.na(df$refBQ_T),]$altBQ_T
# refNM_T has NA when there are no reference reads in tumor. Replace those 51 values with 0.
df[is.na(df$refNM_T),]$refNM_T <- 0
# AN is the number of alleles in the entire cohort. AN = AC * total samples with alt allele
#df[is.na(df$AN),]$AN <- (df[is.na(df$AN),]$AC) * (df[is.na(df$AN),]$SAM_NUM)
# refMN_N is missing one value, just put the mean in there
df[is.na(df$refNM_N),]$refNM_N <- as.numeric(summary(df$refNM_N)[[4]])
# zMQ_T replace with 0 (crude)
df[is.na(df$zMQ_T),]$zMQ_T <- 0
# zBQ_T replace with 0 (crude)
df[is.na(df$zBQ_T),]$zBQ_T <- 0

# Adjust characters to facotrs
df$POWER <- as.factor(df$POWER)
#summary(df)
#summary(is.na(df))
# Convert numeric character columns to numeric
df$VAF_T <- as.numeric(df$VAF_T)

# Convert character vectors to factors
cvecs <- as.vector(names(df[, sapply(df, class) == 'character']))
# Iterate through each column and transform
for (char_vec in cvecs) {
        #print(char_vec)
        #print(table(factor(df[[char_vec]])))
        df[[char_vec]] <- as.factor(df[[char_vec]])
}

#relabel certain columns to ordinal
#ordinals <- c('NUM_TOOLS','FILTER')
#for (ord in ordinals) {
#        df[[ord]] <- ordered(as.factor(df[[ord]])) 
#}

################################################################################

################################################################################
############## Choose the Dependent and Response Variables #####################
################################################################################

# Determine which columns should be considered in the model (response and dependents)
response <- 'VALIDATION'
dependents <- names(df)
# Custom column removal
drops <- c('VALIDATION',"VAR_TYPE","POS","CHROM","REF","ALT")
drops <- append(drops, lose)
dependents <- dependents[dependents %!in% drops]

################################################################################

#' ## Build Training Set

################################################################################
################# Build: Training Set ##########################################
################################################################################

# Generate the training set
train_df <- dplyr::filter(df, VALIDATION %in% c('TP', 'FP'))
# Already deduped above
#train_df <- train_df[!duplicated(train_df[,c('CHROM','POS','REF','ALT','SAMPLE')]),]

################################################################################

#' Training set is small and unbalanced
dim(train_df)
table(train_df$VALIDATION)

#' ## Build Model Formulas

################################################################################
################# Build: Formula ###############################################
################################################################################

# Create the appropriate formula: (only important features)
formula <- as.formula(paste(response, paste(dependents, collapse=" + "), sep=" ~ "))
#feats <- pc2[pc2 %!in% c('CHROM','POS','REF','ALT')]
#formula <- as.formula(paste(response, paste(dependents, collapse=" + "), sep=" ~ "))

################################################################################

#' Model Formula
print(formula)

#' ## Build Test Set
#' The testing set was reduced to only indels called by 3 or 4 callers. The trianing data, the model, and my skills are not strong enough to differentiate accurate calls from indels of 2 or more callers. 

################################################################################
################# Build: Test Dataset ##########################################
################################################################################

# Extract the values in the data that have undergone validation
test_df <- dplyr::filter(df, NUM_TOOLS >= 3 & VALIDATION %!in% c('TP', 'FP'))

################################################################################

#' Test Set
print(dim(test_df))

#' Ensure trianing set and test set in good order

################################################################################
##########  Ensure trianing set and test set are in good order  ################
################################################################################

# Ensure the levels in the trianing set are the only levels represented in the test set
# This is not a fool-proof command, make sure to manually inspect or update command as needed
for (f in dependents) {
        # Only apply to columns that are factors
        if (is.factor(test_df[[f]])) {
                if (sort(unique(train_df[[f]])) != sort(unique(test_df[[f]]))) {
                        print(f)
                        print('training')
                        print(unique(train_df[[f]]))
                        print('testing')
                        print(unique(test_df[[f]]))
                }
        }
}


# Double check factor level length is under 53 for all factors
lvl_cols <- c()
# Remove columns with more than 53 categories
for (col in dependents) {
        factor_len <- length(levels(train_df[[col]]))
        if (factor_len >= 53) {
                lvl_cols <- append(lvl_cols, col)
                print(factor_len)
                print(col)
        }
}

# This approach suggests 6 variables is most sutable
# extract the 8 most suitable variables
#priority_cols_all <- priority_cols
#red_deps <- names(priority_cols) 
#formula <- as.formula(paste(response, paste(red_deps, collapse=" + "), sep=" ~ "))
#length(priority_cols)

################################################################################

#' ## Model Development
#' 
################################################################################
############# Model Development ################################################
################################################################################

# Make sure test dataset has no trianing data values
#table(test_df$VALIDATION)

# Perform a random forest analysis on the training data 
# (random forest does not handle missing values well at all)
# (here we will use the entire validation data)
set.seed(123)
summary(train_df)
rf_model <- randomForest(formula = formula,
                         data = train_df,
                         na.action = na.omit,
                         type = classification,
                         ntree = 1000,
                         mtry = 5)

# Examine general stats about the model
print(rf_model)

# Examine Variable Importance
#class(rf_model)
priority_cols <- rev(importance(rf_model)[order(importance(rf_model)),])
# Remove predictors with no affect and update formula
priority_cols <- priority_cols[priority_cols > 0.01]
red_deps <- names(priority_cols) 
formula <- as.formula(paste(response, paste(red_deps, collapse=" + "), sep=" ~ "))
length(priority_cols)
names(priority_cols)
varImpPlot(rf_model)

# Random forest provides a built-in validation set without any extra work
# Grab the OOB error matrix and examine
err <- rf_model$err.rate
head(err)

# Look at final OOB error rate (last row in err matrix)
oob_err <- err[nrow(err), "OOB"]
print(oob_err)

# Plot the model trained in the previous exercise
plot(rf_model)

# Add a legend since it doesn't have one by default
legend(x = "topright", 
       legend = colnames(err),
       fill = 1:ncol(err))

# Construct object needed to tune the model
train_predictor_df <- train_df[ ,dependents]
train_response_vector <- as.factor(train_df$VALIDATION)

# Tune the model
# Train the mtry parameter based on the OOB error
trees <- c(200, 400, 1000, 2000)
for (ntree in trees) {
        set.seed(1)
        res <- tuneRF(x = train_predictor_df,
                      y = train_response_vector,
                      ntreeTry = ntree,
                      stepFactor = .5,
                      na.action = na.omit)
}

# Examine the results
print(res)

# Manual Search (best seems to be 14 and 400 trees or 1000 trees)
control <- trainControl(method="repeatedcv", number=10, repeats=5, search="grid")
mtrys = c(3,4,5,6,7)
trees <- c(1000,2000)
modellist <- list()
for (m in seq_along(mtrys)) {
        mtry = mtrys[m]
        print(mtry)
        tunegrid <- expand.grid(.mtry=mtry)
        for (t in seq_along(trees)) {
                ntree = trees[t]
                print(ntree)
                set.seed(1)
                fit <- caret::train(formula, data=train_df, method="rf", metric="Accuracy", 
                                    tuneGrid=tunegrid, trControl=control, ntree=ntree,
                                    na.action = na.omit)
                key <- paste0('T_',toString(ntree),'MTRY_',toString(mtry))
                modellist[[key]] <- fit
        } 
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results)

# Landis and Koch considers 0-0.20 as slight, 0.21-0.40 as fair, 0.41-0.60 as moderate, 
# 0.61-0.80 as substantial, and 0.81-1 as almost perfect. Fleiss considers kappas > 0.75 as
# excellent, 0.40-0.75 as fair to good, and < 0.40 as poor.

#' ##Final Generation of Model

################################################################################
########## Final Generation of Model ###########################################
################################################################################

set.seed(123)
# Generate the model with the optimal parameters
rf_model <- randomForest(formula = formula,
                         data = train_df,
                         na.action = na.omit,
                         type = classification,
                         ntree = 1000,
                         mtry = 5,
                         keep.forest=TRUE,
                         importance=TRUE)

# Examine general stats about the model
print(rf_model)
# Random forest provides a built-in validation set without any extra work
# Grab the OOB error matrix and examine
err <- rf_model$err.rate
head(err)
# Look at final OOB error rate (last row in err matrix)
oob_err <- err[nrow(err), "OOB"]
print(oob_err)
# Plot the model trained in the previous exercise
plot(rf_model)

# Add a legend since it doesn't have one by default
legend(x = "right", 
       legend = colnames(err),
       fill = 1:ncol(err))

#' ##Generate a ROC CURVE and AUC

################################################################################
########## Generate a ROC CURVE and AUC ########################################
################################################################################

# http://scg.sdsu.edu/rf_r/
# https://blog.revolutionanalytics.com/2016/11/calculating-auc.html
# https://statquest.org/2018/12/17/roc-and-auc-in-r/ (This site is used for this script)

set.seed(123)
# Generate the model with the optimal parameters
rf_model <- randomForest(formula = formula,
                         data = train_df,
                         na.action = na.omit,
                         type = classification,
                         ntree = 1000,
                         mtry = 5,
                         keep.forest=TRUE,
                         importance=TRUE)

## pty sets the aspect ratio of the plot region. Two options:
##                "s" - creates a square plotting region
##                "m" - (the default) creates a maximal plotting region
par(pty = "s") 

# Save the plot
#pdf('./figures/ROC_RF_all_indels_final.pdf', width = 5, height = 5)
roc(train_df$VALIDATION, rf_model$votes[,1], plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)
#dev.off()

################################################################################

#' Predict Somatic Indels from the Prediction Dataset

################################################################
######## Prepare the dataset for prediction ####################
################################################################

# Extract the values in the data that have not undergone validation
test_df2 <- dplyr::filter(df, NUM_TOOLS >= 3, VALIDATION %!in% c('TP', 'FP'))
dim(test_df2)

# Make sure to set aside the identifier columns
restore <- c("CHROM","POS","REF","ALT",'SAMPLE')
results_keys <- test_df[ , (names(test_df) %in% restore)]

# Convert character vectors to factors
#cvecs <- as.vector(names(test_df[, sapply(test_df, class) == 'character']))
# Iterate through each column and transform
#for (char_vec in cvecs) {
#        print(char_vec)
#print(table(factor(df[[char_vec]])))
#        test_df[[char_vec]] <- as.factor(test_df[[char_vec]])
#}
#priority_cols

# Use the remainder of the dataset and the model to predict
set.seed(123)
pred <- predict(object = rf_model,
                newdata = test_df,
                na.action = na.omit)
# Apply the predicted values to the dataset
test_df$RF_PREDICTED <- pred

#' No Values were ignored (Resulting Predictions)
# Detemrine how many values were ignored
table(test_df$RF_PREDICTED)
table(is.na(test_df$RF_PREDICTED))

#Extract only the predicted true positives
test_tp <- filter(test_df, RF_PREDICTED == 'TP')
test_fp <- filter(test_df, RF_PREDICTED == 'FP')
dim(test_tp)
test_tp['VALIDATION'] <- test_tp['RF_PREDICTED']

# Join the predicted true calls back into the original raw df
test_tp <- as_tibble(test_tp)
dim(test_tp)
test_join <- test_tp %>% dplyr::select(CHROM,POS,REF,ALT,SAMPLE)
dim(test_join)
dim(df)
test_orig <- inner_join(df, test_join, by= c("CHROM","POS","REF","ALT","SAMPLE"))
dim(test_orig)

# Test set
spc <- data.frame()
for (s in unique(test_tp$SAMPLE)){
        print(s)
        pow <- as.character(unique(filter(test_tp, SAMPLE == s) %>% dplyr::select(POWER) %>% unlist()))
        c <- nrow(filter(test_tp, SAMPLE == s))
        spc2 <- data.frame(s,pow,c)
        spc <- rbind(spc,spc2)
}

names(spc) <- c('SAMPLE','POWER','FREQ') 
spc
# Examine the mutation freq relationship
ggplot(spc, aes(x=POWER, y=FREQ)) +
        geom_boxplot()

# Comnine the validated and 
test_tp <- test_tp[,names(test_tp) %!in% 'RF_PREDICTED']
train_tp <- filter(train_df, VALIDATION == 'TP')

# Combine the dataframes
results_tp <- rbind(train_tp,test_tp)
dim(results_tp)

# No duplicated results
#head(results_tp[duplicated(results_tp[,c('CHROM','POS','REF','ALT','SAMPLE')]),])

# Extract significantly muated indel genes
results_tp$SYMBOL <- as.character(results_tp$SYMBOL)
sort(table(results_tp$SYMBOL))
indel_genes <- names(table(results_tp$SYMBOL)[table(results_tp$SYMBOL) >= 1])
# Save indel genes to file
write.table(indel_genes, './data/indel_1_or_more_genes.txt', sep ='\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Check for empty gene symbols
#filter(results_tp, SYMBOL == '')

# Save the results to a tab seperated file
out_file <- "./data/indel_rf_tp_pred_final.txt"
write.table(results_tp, out_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

# Save the results to a csv file
out_file <- "./data/snv_indel_rf_tp_pred_final.csv"
write.table(results_tp, out_file, sep = ",", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

# Important details about model
# 52 Features used (43 features were considered significant influencers of the model)
# Most important features were the agreement between callers, samples, altNM_T, zMQ_T, altMQ_T. refNM_N, refNM_T, NUM_TOOLS, VAF_T
# Type of random forest: classification
# No. of variables tried at each split: 5
# Number of trees: 1000
# OOB estimate of  error rate: 7.94%
# AUC: 95.3%
# training n: 63
# Testing n: 4,966
# Predicted Indels: 1,632
# File of predcited and validated indels: "./data/indel_rf_tp_pred_final.txt"

#' #### End

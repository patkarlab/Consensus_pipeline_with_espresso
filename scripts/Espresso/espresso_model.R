## Function for loading data as VRanges object
load_as_VRanges <-
function(sample_name, sample_path, genome, metadata = TRUE) {

  # Load in as dataframe
  varscan_df <- data.table::fread(sample_path) %>%
    dplyr::filter(Reads2 != 0)

  if(metadata == "VAF"){

  # Convert to VRanges but without refDepth or altDepth
  varscan_output <- VariantAnnotation::VRanges(
    if(length(grep("chr",varscan_df$Chrom))==0) {seqnames = paste0("chr",varscan_df$Chrom)} else {seqnames = varscan_df$Chrom} ,
    ranges = IRanges::IRanges(varscan_df$Position, varscan_df$Position),
    ref = varscan_df$Ref,
    alt = varscan_df$VarAllele,
    sampleNames = sample_name)

  # Add metadata
  S4Vectors::mcols(varscan_output) <- varscan_df %>%
    dplyr::mutate(VAF = Reads2 / (Reads1 + Reads2)) %>%
    dplyr::select(VAF)

  } else {

    # Convert to VRanges
    varscan_output <- VariantAnnotation::VRanges(
      if(length(grep("chr",varscan_df$Chrom))==0) {seqnames = paste0("chr",varscan_df$Chrom)} else {seqnames = varscan_df$Chrom} ,
      ranges = IRanges::IRanges(varscan_df$Position, varscan_df$Position),
      ref = varscan_df$Ref,
      alt = varscan_df$VarAllele,
      refDepth = varscan_df$Reads1,
      altDepth = varscan_df$Reads2,
      sampleNames = sample_name)

    if(metadata==TRUE) {

      # Add metadata
      S4Vectors::mcols(varscan_output) <- varscan_df %>%
        dplyr::mutate(VAF = Reads2 / (Reads1 + Reads2)) %>%
        dplyr::select(VAF, Qual1, Qual2, MapQual1, MapQual2, Reads1Plus, Reads1Minus, Reads2Plus, Reads2Minus, "Varscan_Pval" = Pvalue, GT)

      # Convert quality scores to Rle to save memory
      varscan_output$Qual1 <- methods::as(varscan_output$Qual1, "Rle")
      varscan_output$Qual2 <- methods::as(varscan_output$Qual2, "Rle")
      varscan_output$MapQual1 <- methods::as(varscan_output$MapQual1, "Rle")
      varscan_output$MapQual2 <- methods::as(varscan_output$MapQual2, "Rle")

    }
  }

  # specify genome
  GenomeInfoDb::genome(varscan_output) = genome

  # clean up and remove dataframe
  rm(varscan_df)

  return(varscan_output)
}


# Load required packages
library(VariantAnnotation)
library(Espresso)

# Load other packages
library(doParallel)   # processing samples in parallel
library(magrittr)   # piping and readibility
maf_database="MafDb.gnomADex.r2.1.hs37d5" #If you are using a different MafDB annotation package (see "Installation") change it here accordingly
library(maf_database,character.only = TRUE)

# Specify the genome
genome = "hg19"

#args <- commandArgs()
#print(args)
#dir <- args[1]

#dir <- "/home/tuhina/MRD_DUPLEX/test/variant_calling/"
#files=list.files(dir, pattern = "*.final.cns")
files=list.files(pattern = "*.final.cns")
files
sample_paths <- paste0("./", files)
sample_paths
sample_names <- basename(files)
sample_names

#download.file(url = "https://raw.githubusercontent.com/abelson-lab/Espresso_paper/master/heme_COSMIC/COSMIC_heme_freq10.txt", 
# destfile = "COSMIC_heme_freq10.txt")
#download.file(url = "https://raw.githubusercontent.com/abelson-lab/Espresso_paper/master/heme_COSMIC/COSMIC_heme_freq3.txt", 
# destfile = "COSMIC_heme_freq3.txt")

hemeCOSMIC_10 <- load_recurrent_mutations("/home/pipelines/Consensus_pipeline_with_espresso/scripts/Espresso/COSMIC_heme_freq10.txt", genome = "hg19")
hemeCOSMIC_3 <- load_recurrent_mutations("/home/pipelines/Consensus_pipeline_with_espresso/scripts/Espresso/COSMIC_heme_freq3.txt", genome = "hg19")
#download.file(url = "https://raw.githubusercontent.com/abelson-lab/Espresso_paper/master/heme_COSMIC/hemeCOSMIC_hotspot_range.bed", 
#               destfile = "hemeCOSMIC_hotspot_range.bed")
COSMIC_hotspot_range <- load_bed("mrd_bed_ranges.bed", genome)

# Make cluster and start
cl <- makeCluster(4)
registerDoParallel(cl)

# initialize variant calls
variant_calls <- VRangesList()

# Try for all files
for(i in 1:length(sample_paths)){

  # get sample name and path
  samp_path <- sample_paths[i]
  samp_name <- sample_names[i]

  print(samp_name)

  # get sample as VRanges and annotate with sequence context and MAF
  samp <- load_as_VRanges(samp_name, samp_path, genome) %>%
    sequence_context(., genome, context = 3) %>%
    annotate_MAF(., maf_database, genome)
  
# use sample to generate the error models
   samp_models <- samp %>%
        filter_model_input(., MAF_cutoff = 0.001, VAF_cutoff = 0.05, MAPQ_cutoff = 59, recurrent_mutations = hemeCOSMIC_10) %>%    # preprocessing to clean training set (error models)
       generate_all_models()
#samp_models <- generate_all_models(samp)

# call variants using error models, and aggregate together
    variant_calls[[samp_name]] <- samp %>% 
        intersect_VRanges(., COSMIC_hotspot_range) %>% # only keep positions within the hotspot range we want to call in
      filter_MAPQ(MAPQ_cutoff = 59) %>% 
        call_all_variants(., samp_models)

# cleanup
    rm(samp_name, samp, samp_models)
}

# Unregister Cluster
stopCluster(cl)
rm(cl)
registerDoSEQ()

# Output is in vranges list
variant_calls

# unlist and correct pvalue (bonferroni)
variant_calls_unlisted <- variant_calls %>% unlist() %>% correct_pvalues(., method = "bonferroni")
## if we want to correct sample-by-sample
variant_calls_unlisted <- variant_calls %>% correct_pvalues(., method = "bonferroni") %>% unlist()

# show significant calls 
significant_calls <- variant_calls_unlisted[which(variant_calls_unlisted$corrected_pvalue <= 0.05)]
significant_calls
library('tidyverse')
significant_calls_df <- significant_calls %>% 
  tidyr::as_tibble() %>% 
  dplyr::select(-width, -strand, -totalDepth) %>% 
  dplyr::rename("chr" = seqnames) %>% 
  dplyr::mutate(chr = sub(pattern="chr", replacement="",x=chr))
significant_calls_df %>% readr::write_delim("espresso_calls.txt", delim = "\t")

significant_calls_df %>% head() %>% print()

significant_calls %>% 
  VariantAnnotation::asVCF() %>% 
  VariantAnnotation::writeVcf("espresso_calls.vcf")

#library("cellbaseR")
#anno <- annotate_variants(significant_calls, genome='hg19')
#anno %>% head()

#df1 <- select(anno, chromosome, start, reference, alternate, id, displayConsequenceType)
#write.csv(df1, file = "anno_3.txt")

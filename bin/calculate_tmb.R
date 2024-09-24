#!/usr/bin/env Rscript
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("-m", "--maf"), type="character", default=NULL, 
              help="maf file"),
  make_option(c("-o", "--out"), type="character", default="tmb.tsv", 
              help="output file name [default= %default]"),
  make_option(c("-d", "--depth"), type="integer", default=250, 
              help="tumor depth"),
  make_option(c("-c", "--count"), type="integer", default=5, 
              help="tumor alt count"),
  make_option(c("-v", "--vaf"), type="double", default=0.05, 
              help="variant allele frequency")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

calculate_tmb <- function(maf_file, out_tsv, tumor_coverage, alt_count, tumor_vaf) {
  # read input files
  tmb_maf <- read.delim(maf_file, header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  tmb_df <- tmb_maf %>% 
                    select(Hugo_Symbol, Center, Chromosome, Start_Position, End_Position,
                               Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
                               Strand, Variant_Classification, Variant_Type,Tumor_Sample_Barcode,
                               HGVSc, HGVSp, t_depth, t_ref_count, t_alt_count,
                               n_depth, n_ref_count, n_alt_count, Consequence,
                               CLIN_SIG, Existing_variation, IMPACT) %>%
                    filter(Variant_Classification %in% c('Missense_Mutation',
                                         'Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins',
                                         'Nonsense_Mutation','Stop_Gained', 'Synonymous_Variant') 
                          & t_depth >= tumor_coverage & t_alt_count >= alt_count ) %>%
                    mutate(var_id=paste(Chromosome, Start_Position, End_Position,
                            Reference_Allele,Tumor_Seq_Allele1,sep="-"),
                            t_vaf=t_alt_count/t_depth,
                            n_vaf=n_alt_count/n_depth) %>%
                    filter( t_vaf >= tumor_vaf) %>%
                    arrange(Hugo_Symbol, var_id, Center) 
  # remove duplicates
  tmb_df <- tmb_df[!duplicated(tmb_df$var_id),]
  # calcualte TMB
  tmb <- data.frame(unclass(table(tmb_df$Tumor_Sample_Barcode)))
  t_value <- apply(tmb,1, function(x) round(x/round(5387125/10^6,2),2)  )
  write.table(data.frame(t_value), file=out_tsv, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
}
calculate_tmb(opt$maf, opt$out, opt$depth, opt$count, opt$vaf)
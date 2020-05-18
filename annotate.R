# Tempus Bioinformatics Technical Challenge
# For this challenge, you are asked to prototype a variant annotation tool. We will provide you with
# a VCF file, and you will create a small software program to output a table annotating each
# variant in the file. Each variant must be annotated with the following pieces of information:

# 1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple
#    possibilities, annotate with the most deleterious possibility.
# 2. Depth of sequence coverage at the site of variation.
# 3. Number of reads supporting the variant.
# 4. Percentage of reads supporting the variant versus those supporting reference reads.
# 5. Allele frequency of variant from Broad Institute ExAC Project API
#    (API documentation is available here: http://exac.hms.harvard.edu/)
# 6. Additional optional information from ExAC that you feel might be relevant.

# For this project please upload all relevant code (written in whatever language you like) along
# with the annotated VCF file to a Github account and provide the link to the below email address.
# Please note that work will be assessed based on quality of code and documentation more-so
# than the annotation.
# Please direct any queries to Stephen Bush (stephen.bush@tempus.com).


### OPTS/FRESH START
rm(list = ls(all = TRUE))
options(warn = -1)
getwd()



### REQUIRED LIBS & DBS
# devtools::install_github(repo="knausb/vcfR") 
# BiocManager::install("MafDb.ExAC.r1.0.nonTCGA.hs37d5") # MAF data from ExAC release 1.0 (subset of nonTCGA exomes for hs37d5)
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37") # SNP locations for Homo sapiens (dbSNP Build 144)
pacman::p_load(vcfR, # Efficient VCF parser
               MafDb.ExAC.r1.0.nonTCGA.hs37d5, # Database of minor allele frequency data from the Exome Aggregation Consortium
               SNPlocs.Hsapiens.dbSNP144.GRCh37, # Database of SNP locations, needed for lookup in ExAC by SNP ID
               GenomicRanges # Used for efficient range overlap of SNPs
               ) 



### DATA
# VCF is small enought to load in memory, so let's do that
vcf_file <- file.path("data", "Challenge_data.vcf")
vcf <- read.vcfR(vcf_file, verbose = TRUE, cols = seq(1, 11))
vcf_info <- INFO2df(vcf) # parses INFO column into a data.frame 
vcf_info_metadata <- metaINFO2df(vcf) # metadata for columns above
to_return <- list() 



### 1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.)
# Assuming most deleterious type listed first.
to_return$TYPE <- sapply(strsplit(vcf_info$TYPE, ",", fixed=TRUE), function(x) (x[1]))



### 2. Depth of sequence coverage at the site of variation.
# vcf_info_metadata[grep("Total read depth at the locus", vcf_info_metadata$Description), ] # DP
to_return$DEPTH <- vcf_info$DP



### 3. Number of reads supporting the variant.
# vcf_info_metadata[grep("Total number of alternate alleles in called genotypes", vcf_info_metadata$Description), ] # AC
to_return$SUPPORT <- vcf_info$AC



### 4. Percentage of reads supporting the variant versus those supporting reference reads.
# vcf_info_metadata[grep("Estimated allele frequency in the range", vcf_info_metadata$Description), ] # AF
to_return$PCT_VARIANT <- as.numeric(vcf_info$AF) * 100 # convert to pct, NAs included



### 5. Allele frequency of variant from Broad Institute ExAC Project API
mafdb <- MafDb.ExAC.r1.0.nonTCGA.hs37d5
snpdb <- SNPlocs.Hsapiens.dbSNP144.GRCh37 
coords <- data.frame(getFIX(vcf), stringsAsFactors = FALSE)
coords$END <- as.numeric(coords$POS) + 1
coords_gr <- GenomicRanges::makeGRangesFromDataFrame(coords, 
                                                     keep.extra.columns = FALSE, 
                                                     seqnames.field = "CHROM", 
                                                     start.field = "POS", 
                                                     end.field = "END")
annotated_snps <- snpsByOverlaps(snpdb, coords_gr) # overlap SNPs in VCF with SNP database 
annotated_snps_af <- gscores(mafdb, annotated_snps)
mcols(coords_gr)$AF <- annotated_snps_af$AF[match(start(coords_gr), start(annotated_snps_af))]
mcols(coords_gr)$SNP_ID <- annotated_snps_af$RefSNP_id[match(start(coords_gr), start(annotated_snps_af))]
mcols(coords_gr)$ALLELES_AS_AMBIG <- annotated_snps_af$alleles_as_ambig[match(start(coords_gr), start(annotated_snps_af))]
to_return$AF <- coords_gr$AF
to_return$SNP_ID <- coords_gr$SNP_ID
to_return$ALLELES_AS_AMBIG <- coords_gr$ALLELES_AS_AMBIG



### CONSTRUCT FINAL DF
init_info <- coords[, c("CHROM", "POS", "REF", "ALT")]
final_df <- cbind.data.frame(init_info, do.call(cbind.data.frame, to_return))



### EXPORT
out <- file.path("output", "annotated.rds")
saveRDS(final_df, compress = TRUE, out)
write.table(final_df, gsub(".rds", ".txt", out), sep = "\t", quote = FALSE)
save(list = ls(all=TRUE), file = gsub("annotated.rds", "env.Rdata", out))


# Tempus Challenge

The goal of this challenge is to annotate each variant in the [supplied VCF](https://github.com/ShanSabri/tempus/blob/master/data/Challenge_data.vcf) file with the following features: 


1. Type of variation
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API


Challenge details posted [here](https://github.com/ShanSabri/tempus/blob/master/data/Tempus_Bioinformatics_Challenge.pdf). 

Variants are queried by locations and alleles (for Homo sapiens) from [NCBI dbSNP Build 144](http://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP144.GRCh37.html). SNP IDs are then used to query [ExAC release 1.0 subset of nonTCGA exomes](http://bioconductor.org/packages/release/data/annotation/html/MafDb.ExAC.r1.0.nonTCGA.hs37d5.html) for minor allele frequency data. 

## Installation & Usage

```bash
git clone https://github.com/ShanSabri/tempus.git
cd tempus
R CMD BATCH annotate.R
```

## Output
Included in the [`output`](https://github.com/ShanSabri/tempus/tree/master/output) directory is the annotated data.frame object (rds) and tab-seperated text annotation file (txt). 

```
tempus/output
├── [151K]  annotated.rds
└── [390K]  annotated.txt

 540K used in 0 directories, 2 files
```

The columns of these data correspond to: 

1. Chromosome
2. Position 
3. Reference
4. Alternative 
5. Type of variation 
6. Depth of coverage
7. Number of reads supporting variant 
8. Percentage of reads supporting the variant
9. Allele frequency of variant from ExAC; NA if not found in SNPdb or annotated in `MafDb.ExAC.r1.0.nonTCGA.hs37d5_3.10.0`
10. SNP ID; NA if not found in SNPdb
11. Alleles for each SNP represented by an IUPAC nucleotide ambiguity code; NA if not found in SNP db

<details><summary>Session info</summary>
<p>

``` r
> sessionInfo()
R version 3.6.2 (2019-12-12)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Catalina 10.15.4

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] SNPlocs.Hsapiens.dbSNP144.GRCh37_0.99.20 BSgenome_1.54.0                          rtracklayer_1.46.0                      
 [4] Biostrings_2.54.0                        XVector_0.26.0                           MafDb.ExAC.r1.0.nonTCGA.hs37d5_3.10.0   
 [7] GenomicScores_1.10.0                     GenomicRanges_1.38.0                     GenomeInfoDb_1.22.0                     
[10] IRanges_2.20.2                           S4Vectors_0.24.3                         BiocGenerics_0.32.0                     
[13] vcfR_1.11.0                             

loaded via a namespace (and not attached):
 [1] Biobase_2.46.0                httr_1.4.1                    bit64_0.9-7                   viridisLite_0.3.0            
 [5] AnnotationHub_2.18.0          splines_3.6.2                 shiny_1.4.0                   assertthat_0.2.1             
 [9] memuse_4.1-0                  interactiveDisplayBase_1.24.0 BiocManager_1.30.10           BiocFileCache_1.10.2         
[13] blob_1.2.1                    GenomeInfoDbData_1.2.2        Rsamtools_2.2.3               yaml_2.2.1                   
[17] BiocVersion_3.10.1            pillar_1.4.4                  RSQLite_2.2.0                 lattice_0.20-40              
[21] glue_1.4.1                    digest_0.6.25                 promises_1.1.0                htmltools_0.4.0              
[25] httpuv_1.5.2                  Matrix_1.2-18                 XML_3.99-0.3                  pkgconfig_2.0.3              
[29] zlibbioc_1.32.0               purrr_0.3.4                   xtable_1.8-4                  later_1.0.0                  
[33] BiocParallel_1.20.1           tibble_3.0.1                  mgcv_1.8-31                   ellipsis_0.3.1               
[37] pacman_0.5.1                  SummarizedExperiment_1.16.1   magrittr_1.5                  crayon_1.3.4                 
[41] mime_0.9                      memoise_1.1.0                 nlme_3.1-145                  MASS_7.3-51.5                
[45] vegan_2.5-6                   tools_3.6.2                   lifecycle_0.2.0               matrixStats_0.55.0           
[49] cluster_2.1.0                 DelayedArray_0.12.2           AnnotationDbi_1.48.0          compiler_3.6.2               
[53] rlang_0.4.6                   grid_3.6.2                    RCurl_1.98-1.1                rstudioapi_0.11              
[57] rappdirs_0.3.1                bitops_1.0-6                  DBI_1.1.0                     curl_4.3                     
[61] R6_2.4.1                      GenomicAlignments_1.22.1      dplyr_0.8.5                   pinfsc50_1.1.0               
[65] fastmap_1.0.1                 bit_1.1-15.2                  permute_0.9-5                 ape_5.3                      
[69] Rcpp_1.0.4.6                  vctrs_0.3.0                   dbplyr_1.4.2                  tidyselect_1.1.0     
```
</p>
</details>

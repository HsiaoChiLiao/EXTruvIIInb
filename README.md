# EXTruvIIInb
### Author: Hsiao-Chi Liao
### Date: 11 October 2024
<br>

# Introduction
Single-cell multimodal technologies provide an opportunity to study biological mechanisms more comprehensively. The integrated analysis of matched single-cell transcriptomics (mRNA expression) and proteomics (protein abundance) data can help reveal biological insights that would not have been possible from separate analyses of each modality.

## Motivation
Unwanted variation from sources such as shared batches and modality-specific library size effects inevitably exists in the data from both domains. If not properly corrected, the unwanted variation can potentially lead to misleading conclusions being made from the downstream analyses. However, none of the existing methods adequately takes biological factors into account in their models, often leading to the removal of both biological and unwanted variations when they are associated. To address the limitations, we propose the Extended RUV-III-NB. Our model accounts for biology and decomposes unwanted variation into joint and modality-specific components, providing a clearer understanding of the unwanted variation present in the data.

## The RUV-III-NB framework and the Extended RUV-III-NB model
RUV-III-NB (Salim et al., 2022) is a method developed for removing unwanted variation from single-cell mRNA UMI counts, using a subset of annotated cells (as low as 5\%) for model training. RUV-III-NB infers unwanted factors without requiring them to be specified in the model, and provides two types of corrected counts: logPAC (percentile adjusted counts on the natural logarithm scale) and Pearson residuals. logPAC is derived from the probability mass function of the negative binomial distribution, while Pearson residuals are calculated using the first and second moments  (the expected value and the variance) of the negative binomial distribution. Pearson residuals tend to be more robust when count data deviates from the negative binomial distribution, whereas logPAC provides more accurate results when the data approximates negative binomial.

The model of the Extended RUV-III-NB which shares the same general structure with the RUV-III-NB model is as follows.

$$
\log(\boldsymbol{\mu_f})=\zeta_f\boldsymbol{1}+\boldsymbol{M}\boldsymbol{\beta_f}+\boldsymbol{W}\boldsymbol{\alpha_f},
$$

The Extended RUV-III-NB has a different design for the unwanted component $\boldsymbol{W\alpha_f}$, where
$\boldsymbol{W}=
\begin{bmatrix}
\boldsymbol{W_1}|&\boldsymbol{W_2}|&\boldsymbol{W_3}
\end{bmatrix}
_{n \times (k_1+k_2+k_3)}$

$\boldsymbol{\alpha_f}=
\begin{bmatrix}
\boldsymbol{\alpha_{1,f}} \\
\boldsymbol{0} \\
\boldsymbol{\alpha_{3,f}}
\end{bmatrix}
_{(k_1+k_3) \times 1}$ for $f$ from the mRNA domain, and 

$\boldsymbol{\alpha_f}=
\begin{bmatrix}
\boldsymbol{0} \\
\boldsymbol{\alpha_{2,f}} \\
\boldsymbol{\alpha_{3,f}}
\end{bmatrix}
_{(k_2+k_3) \times 1}$ for $f$ from the protein domain.
<br>


# Analysing Data with the EXTruvIIInb Package
## Installation
We firstly install the Extended RUV-III-NB package with the following.
``` r
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
devtools::install_github(""HsiaoChiLiao/EXTruvIIInb", build_vignettes=FALSE)

library(EXTruvIIInb)
```

Below packages are required for the analysis in this Vignette and are loaded below.
``` r
library(ruvIIInb)
library(DelayedArray)
library(SingleCellExperiment)
library(ggplot2)
library(GGally)
library(RColorBrewer)
library(uwot)
library(lisi)
library(gridExtra)
library(parallel)
library(cowplot)
library(metapod)
```

# Anslysing the Triana dataset with EXTruvIIInb
## The Triana dataset
This human bone marrow dataset was published by Triana et al. (2021, Nature Immunology). After quality control (QC), this dataset contains 422 targeted mRNA features, and 97 proteins measured by oligo-tagged antibodies (Antibody-Derived Tags, ADT). A total of 9,663 cells were randomly selected from 49,057 cells from 3 healthy young adults and 3 healthy aged adults with each sample from a different batch, are used for the below demonstration.

The build-in data `triana_data` can be called by the function `data()` and it contains (1) the raw counts of mRNA and ADT features (`rna_raw` and `adt_raw`), (2) the metadata of the cells (cell annotation and batch information), (3) the negative control sets for inferring ${\boldsymbol W_1}$, ${\boldsymbol W_2}$, and ${\boldsymbol W_3}$

``` r
#load the build-in data
data("triana_data")
str(triana_data)

#the components in the "triana_data" list
rna_raw <- triana_data$rna_raw
adt_raw <- triana_data$adt_raw
meta <- triana_data$meta
nc1 <- triana_data$nc1
nc2 <- triana_data$nc2 
nc3 <- triana_data$nc3

#the raw counts
rna_raw[1:5,1:3]
adt_raw[1:5,1:3]

#the metadata for each cell
head(triana_data$meta)

# check the details
help(triana_data, package = "EXTruvIIInb")
```

## Running the vanilla EXTruvIIInb (the version without constraints)
Based on our experience, a majority of the library size effect of each modality can usually be captured with one unwanted factor. Thus, we set $k_1=k_2=1$ for the inference of the modality-specific library size (and other unknown unwanted variation), and $k_3=2$ for capturing the joint unwanted effects.

``` r
# Using known cell types to define pseudo-replicates
anno <- meta$cellType_broad

# Construct the replicate matrix M using the known cell-types
unique_ctype<- unique(anno)[!is.na(unique(anno))]
M <- matrix(0, nrow(meta), length(unique_ctype))
dim(M) 
for(i in 1:length(unique_ctype)){
    M[anno == unique_ctype[i], i] <- 1
}

# Running with the vanilla Extended RUV-III-NB        
a <- Sys.time()
extruv3nb_broad <- extFastruvIIInb_vanilla(Y1=DelayedArray(t(rna_raw)),
                                           Y2=DelayedArray(t(adt_raw)),
                                           M=M,
                                           ctl1=nc1,
                                           ctl2=nc2,
                                           ctl3=nc3,
                                           k1=1,
                                           k2=1,
                                           k3=2,
                                           ncores = 2)
b <- Sys.time()
print("finished!")
print(b-a)

# check the details
help(extFastruvIIInb_vanilla, package = "EXTruvIIInb")
```

The `extruv3nb_broad` is an R object that contains the corrected data and the parameters estimated by the model.

# THANK YOU

Thanks for carrying out analyses with our `EXTruvIIInb` package. Please feel free to raise any issues and/or questions on our GitHub page and/or send them to Hsiao-Chi [hsiaochi.liao@student.unimelb.edu.au](mailto:hsiaochi.liao@student.unimelb.edu.au).

## Information about the R session when this vignette was built

``` r
> sessionInfo()
R version 4.3.3 (2024-02-29)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.6.3

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Australia/Melbourne
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] devtools_2.4.5 usethis_3.0.0 

loaded via a namespace (and not attached):
  [1] splines_4.3.3               later_1.3.2                 bitops_1.0-9               
  [4] tibble_3.2.1                R.oo_1.26.0                 rpart_4.1.23               
  [7] lifecycle_1.0.4             rstatix_0.7.2               rprojroot_2.0.4            
 [10] edgeR_4.0.16                doParallel_1.0.17           mirt_1.42                  
 [13] globals_0.16.3              lattice_0.22-5              MASS_7.3-60.0.1            
 [16] backports_1.5.0             magrittr_2.0.3              limma_3.58.1               
 [19] Hmisc_5.1-3                 rmarkdown_2.28              yaml_2.3.10                
 [22] remotes_2.5.0               metapod_1.10.1              httpuv_1.6.15              
 [25] MultBiplotR_23.11.0         sessioninfo_1.2.2           pkgbuild_1.4.4             
 [28] reticulate_1.39.0           cowplot_1.1.3               pbapply_1.7-2              
 [31] keras_2.15.0                abind_1.4-8                 pkgload_1.4.0              
 [34] zlibbioc_1.48.2             audio_0.1-11                GenomicRanges_1.54.1       
 [37] purrr_1.0.2                 R.utils_2.12.3              BiocGenerics_0.48.1        
 [40] RCurl_1.98-1.16             nnet_7.3-19                 rgl_1.3.1                  
 [43] GenomeInfoDbData_1.2.11     IRanges_2.36.0              S4Vectors_0.40.2           
 [46] irlba_2.3.5.1               listenv_0.9.1               testthat_3.2.1.1           
 [49] RPushbullet_0.3.4           vegan_2.6-8                 parallelly_1.38.0          
 [52] permute_0.9-7               DelayedMatrixStats_1.24.0   dae_3.2.28                 
 [55] codetools_0.2-19            DelayedArray_0.28.0         xml2_1.3.6                 
 [58] tidyselect_1.2.1            futile.logger_1.4.3         matrixStats_1.4.1          
 [61] stats4_4.3.3                base64enc_0.1-3             roxygen2_7.3.2             
 [64] jsonlite_1.8.9              polycor_0.8-1               e1071_1.7-16               
 [67] ellipsis_0.3.2              progressr_0.14.0            Formula_1.2-5              
 [70] iterators_1.0.14            foreach_1.5.2               tools_4.3.3                
 [73] snow_0.4-4                  Rcpp_1.0.13                 glue_1.8.0                 
 [76] mnormt_2.1.1                gridExtra_2.3               Rttf2pt1_1.3.12            
 [79] tfruns_1.5.3                SparseArray_1.2.4           xfun_0.48                  
 [82] mgcv_1.9-1                  admisc_0.36                 MatrixGenerics_1.14.0      
 [85] GenomeInfoDb_1.38.8         dplyr_1.1.4                 withr_3.0.1                
 [88] formatR_1.14                beepr_2.0                   fastmap_1.2.0              
 [91] fansi_1.0.6                 caTools_1.18.3              digest_0.6.37              
 [94] R6_2.5.1                    mime_0.12                   colorspace_2.1-1           
 [97] gtools_3.9.5                R.methodsS3_1.8.2           tidyr_1.3.1                
[100] utf8_1.2.4                  generics_0.1.3              data.table_1.16.0          
[103] class_7.3-22                tryCatchLog_1.3.1           SimDesign_2.17.1           
[106] htmlwidgets_1.6.4           S4Arrays_1.2.1              whisker_0.4.1              
[109] pkgconfig_2.0.3             gtable_0.3.5                tensorflow_2.16.0          
[112] SingleCellExperiment_1.24.0 XVector_0.42.0              brio_1.1.5                 
[115] htmltools_0.5.8.1           carData_3.0-5               profvis_0.4.0              
[118] scales_1.3.0                Biobase_2.62.0              png_0.1-8                  
[121] geometry_0.5.0              lambda.r_1.2.4              knitr_1.48                 
[124] rstudioapi_0.16.0           checkmate_2.3.2             nlme_3.1-164               
[127] magic_1.6-1                 curl_5.2.3                  proxy_0.4-27               
[130] cachem_1.1.0                stringr_1.5.1               KernSmooth_2.23-22         
[133] parallel_4.3.3              miniUI_0.1.1.1              RcppZiggurat_0.1.6         
[136] foreign_0.8-86              extrafont_0.19              desc_1.4.3                 
[139] pillar_1.9.0                grid_4.3.3                  ruvIIInb_0.8.2.0           
[142] vctrs_0.6.5                 gplots_3.2.0                ggpubr_0.6.0               
[145] urlchecker_1.0.1            promises_1.3.0              car_3.1-3                  
[148] dunn.test_1.3.6             xtable_1.8-4                Deriv_4.1.6                
[151] cluster_2.1.6               dcurver_0.9.2               htmlTable_2.4.3            
[154] extrafontdb_1.0             GPArotation_2024.3-1        evaluate_1.0.0             
[157] zeallot_0.1.0               futile.options_1.0.1        mvtnorm_1.3-1              
[160] cli_3.6.3                   locfit_1.5-9.10             compiler_4.3.3             
[163] rlang_1.1.4                 crayon_1.5.3                ThreeWay_1.1.3             
[166] ggsignif_0.6.4              future.apply_1.11.2         matlib_1.0.0               
[169] ZIM_1.1.0                   plyr_1.8.9                  fs_1.6.4                   
[172] psych_2.4.6.26              stringi_1.8.4               deldir_2.0-4               
[175] BiocParallel_1.36.0         munsell_0.5.1               EXTruvIIInb_0.0.0.9000     
[178] Matrix_1.6-5                sparseMatrixStats_1.14.0    future_1.34.0              
[181] ggplot2_3.5.1               statmod_1.5.0               shiny_1.9.1                
[184] SummarizedExperiment_1.32.0 Rfast_2.1.0                 broom_1.0.7                
[187] memoise_2.0.1               RcppParallel_5.1.9          xgboost_1.7.8.1     
```

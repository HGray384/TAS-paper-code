# TAS-paper-code

## About

The code reposity for Gray, H., Leday, G.G., Vallejos, C.A. and Richardson, S., 2018. Shrinkage estimation of large covariance matrices using multiple shrinkage targets. arXiv preprint arXiv:1809.08024. The code generates all figures from the paper and was a collaborative effort between Gray, H., Leday, G.G., Vallejos, C.A. The corresponding author may be contacted at gwenael.leday@mrc-bsu.cam.ac.uk (GGRL). Any bug reports should be submitted as issues to https://github.com/HGray384/TAS-paper-code/issues. 

## How to use

Pay careful attention to install the necessary R packages detailed at the beginning of each script and setwd() to the base directory of the folder. No intermediary results are stored since each script should take no more than an hour to run with basic computing resources.

Each file name corresponds to the figures in the paper that are generated in its output. The scripts are best run section-by-section to generate the stated output of a given section, however they can be run fully provided a base directory is specified in the basedir variable near the start of the script.

For Figs-5-6.r the data must first be downloaded from http://tcpaportal.org/tcpa/download.html then selecting the /Downloads/4.2/TCGA/Pan-Can 32/Level 4/TCGA-PANCAN32-L4.zip path option and downloading this dataset. The directory of this dataset on your local machine must then be set in the datadir variable in the script.

## Session info

The sessionInfo() of the last successful run:

R version 3.5.0 (2018-04-23)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.6

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] corrplot_0.84        BDgraph_2.51         data.table_1.11.4   
 [4] ggplot2_3.0.0        ShrinkCovMat_1.2.0   corpcor_1.6.9       
 [7] MVN_5.5              org.Hs.eg.db_3.6.0   AnnotationDbi_1.42.1
[10] IRanges_2.14.10      S4Vectors_0.18.3     Biobase_2.40.0      
[13] BiocGenerics_0.26.0  cgdsr_1.2.10         Matrix_1.2-14       
[16] TAS_1.0              MCMCpack_1.4-3       MASS_7.3-50         
[19] coda_0.19-1          abind_1.4-5         

loaded via a namespace (and not attached):
  [1] colorspace_1.3-2      mvoutlier_2.0.9       class_7.3-14         
  [4] modeltools_0.2-22     rio_0.5.10            mclust_5.4.1         
  [7] rprojroot_1.3-2       pls_2.7-0             rstudioapi_0.7       
 [10] MatrixModels_0.4-1    cvTools_0.3.2         flexmix_2.3-14       
 [13] bit64_0.9-7           mvtnorm_1.0-8         xml2_1.2.0           
 [16] splines_3.5.0         sROC_0.1-2            R.methodsS3_1.7.1    
 [19] mnormt_1.5-5          robustbase_0.93-2     knitr_1.20           
 [22] robCompositions_2.0.8 mcmc_0.9-5            kernlab_0.9-27       
 [25] cluster_2.0.7-1       R.oo_1.22.0           rrcov_1.4-4          
 [28] readr_1.1.1           compiler_3.5.0        httr_1.3.1           
 [31] backports_1.1.2       assertthat_0.2.0      lazyeval_0.2.1       
 [34] htmltools_0.3.6       quantreg_5.36         tools_3.5.0          
 [37] igraph_1.2.2          bindrcpp_0.2.2        gtable_0.2.0         
 [40] glue_1.3.0            dplyr_0.7.6           Rcpp_0.12.19         
 [43] carData_3.0-1         cellranger_1.1.0      trimcluster_0.1-2.1  
 [46] zCompositions_1.1.1   sgeostat_1.0-27       nlme_3.1-137         
 [49] fpc_2.1-11.1          psych_1.8.4           lmtest_0.9-36        
 [52] laeken_0.4.6          stringr_1.3.1         openxlsx_4.1.0       
 [55] rvest_0.3.2           DEoptimR_1.0-8        zoo_1.8-3            
 [58] scales_1.0.0          VIM_4.7.0             hms_0.4.2            
 [61] SparseM_1.77          RColorBrewer_1.1-2    yaml_2.2.0           
 [64] curl_3.2              NADA_1.6-1            memoise_1.1.0        
 [67] reshape_0.8.7         stringi_1.2.4         RSQLite_2.1.1        
 [70] nortest_1.0-4         pcaPP_1.9-73          e1071_1.7-0          
 [73] energy_1.7-5          boot_1.3-20           zip_1.0.0            
 [76] truncnorm_1.0-8       rlang_0.2.2           pkgconfig_2.0.2      
 [79] moments_0.14          prabclus_2.2-6        matrixStats_0.54.0   
 [82] evaluate_0.11         lattice_0.20-35       purrr_0.2.5          
 [85] bindr_0.1.1           bit_1.1-14            tidyselect_0.2.4     
 [88] GGally_1.4.0          plyr_1.8.4            magrittr_1.5         
 [91] R6_2.2.2              DBI_1.0.0             withr_2.1.2          
 [94] pillar_1.3.0          haven_1.1.2           foreign_0.8-71       
 [97] survival_2.42-6       sp_1.3-1              nnet_7.3-12          
[100] tibble_1.4.2          crayon_1.3.4          car_3.0-2            
[103] rmarkdown_1.10        grid_3.5.0            readxl_1.1.0         
[106] blob_1.1.1            forcats_0.3.0         vcd_1.4-4            
[109] digest_0.6.15         diptest_0.75-7        munsell_0.5.0        
[112] viridisLite_0.3.0     kableExtra_0.9.0  

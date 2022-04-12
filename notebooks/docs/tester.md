```R
print(getwd())
library(data.table)
print(renv::status())
```

    [1] "/Users/quangnguyen/research/microbe_set_trait/notebooks"
    * The project is already synchronized with the lockfile.
    $library
    $R
    $R$Version
    [1] "4.1.3"
    
    $R$Repositories
    $R$Repositories$CRAN
    [1] "https://cran.rstudio.com"
    
    
    
    $Bioconductor
    $Bioconductor$Version
    [1] "3.14"
    
    
    $Python
    $Python$Version
    [1] "3.8.9"
    
    $Python$Type
    [1] "virtualenv"
    
    $Python$Name
    [1] "./renv/python/virtualenvs/renv-python-3.8"
    
    
    $Packages
    $Packages$HMP16SData
    $Packages$HMP16SData$Package
    [1] "HMP16SData"
    
    $Packages$HMP16SData$Version
    [1] "1.14.0"
    
    $Packages$HMP16SData$Source
    [1] "Bioconductor"
    
    $Packages$HMP16SData$Depends
    [1] "R (>= 4.1.0)"         "SummarizedExperiment"
    
    $Packages$HMP16SData$Imports
     [1] "AnnotationHub" "ExperimentHub" "S4Vectors"     "assertthat"   
     [5] "dplyr"         "kableExtra"    "knitr"         "magrittr"     
     [9] "methods"       "readr"         "stringr"       "tibble"       
    [13] "utils"        
    
    $Packages$HMP16SData$git_url
    [1] "https://git.bioconductor.org/packages/HMP16SData"
    
    $Packages$HMP16SData$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$HMP16SData$git_last_commit
    [1] "557973d"
    
    $Packages$HMP16SData$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$HMP16SData$Hash
    [1] "45981b9d777b5fead36ab81ae34a1312"
    
    
    $Packages$tinytex
    $Packages$tinytex$Package
    [1] "tinytex"
    
    $Packages$tinytex$Version
    [1] "0.37"
    
    $Packages$tinytex$Source
    [1] "Repository"
    
    $Packages$tinytex$Imports
    [1] "xfun (>= 0.29)"
    
    $Packages$tinytex$Repository
    [1] "CRAN"
    
    $Packages$tinytex$Hash
    [1] "a80abeb527a977e4bef21873d29222dd"
    
    
    $Packages$clisymbols
    $Packages$clisymbols$Package
    [1] "clisymbols"
    
    $Packages$clisymbols$Version
    [1] "1.2.0"
    
    $Packages$clisymbols$Source
    [1] "Repository"
    
    $Packages$clisymbols$Repository
    [1] "CRAN"
    
    $Packages$clisymbols$Hash
    [1] "96c01552bfd5661b9bbdefbc762f4bcd"
    
    
    $Packages$diffobj
    $Packages$diffobj$Package
    [1] "diffobj"
    
    $Packages$diffobj$Version
    [1] "0.3.5"
    
    $Packages$diffobj$Source
    [1] "Repository"
    
    $Packages$diffobj$Depends
    [1] "R (>= 3.1.0)"
    
    $Packages$diffobj$Imports
    [1] "crayon (>= 1.3.2)" "tools"             "methods"          
    [4] "utils"             "stats"            
    
    $Packages$diffobj$Repository
    [1] "CRAN"
    
    $Packages$diffobj$Hash
    [1] "bcaa8b95f8d7d01a5dedfd959ce88ab8"
    
    
    $Packages$lattice
    $Packages$lattice$Package
    [1] "lattice"
    
    $Packages$lattice$Version
    [1] "0.20-45"
    
    $Packages$lattice$Source
    [1] "Repository"
    
    $Packages$lattice$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$lattice$Imports
    [1] "grid"      "grDevices" "graphics"  "stats"     "utils"    
    
    $Packages$lattice$Repository
    [1] "CRAN"
    
    $Packages$lattice$Hash
    [1] "b64cdbb2b340437c4ee047a1f4c4377b"
    
    
    $Packages$GSVA
    $Packages$GSVA$Package
    [1] "GSVA"
    
    $Packages$GSVA$Version
    [1] "1.42.0"
    
    $Packages$GSVA$Source
    [1] "Bioconductor"
    
    $Packages$GSVA$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$GSVA$Imports
     [1] "methods"              "stats"                "utils"               
     [4] "graphics"             "S4Vectors"            "IRanges"             
     [7] "Biobase"              "SummarizedExperiment" "GSEABase"            
    [10] "Matrix"               "parallel"             "BiocParallel"        
    [13] "SingleCellExperiment" "sparseMatrixStats"    "DelayedArray"        
    [16] "DelayedMatrixStats"   "HDF5Array"            "BiocSingular"        
    
    $Packages$GSVA$git_url
    [1] "https://git.bioconductor.org/packages/GSVA"
    
    $Packages$GSVA$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$GSVA$git_last_commit
    [1] "c99b10b"
    
    $Packages$GSVA$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$GSVA$Hash
    [1] "95684adcece76d9d48c78a933ae01f87"
    
    
    $Packages$haven
    $Packages$haven$Package
    [1] "haven"
    
    $Packages$haven$Version
    [1] "2.4.3"
    
    $Packages$haven$Source
    [1] "Repository"
    
    $Packages$haven$Depends
    [1] "R (>= 3.2)"
    
    $Packages$haven$Imports
    [1] "forcats (>= 0.2.0)" "hms"                "methods"           
    [4] "readr (>= 0.1.0)"   "rlang (>= 0.4.0)"   "tibble"            
    [7] "tidyselect"         "vctrs (>= 0.3.0)"  
    
    $Packages$haven$LinkingTo
    [1] "cpp11"
    
    $Packages$haven$Repository
    [1] "CRAN"
    
    $Packages$haven$Hash
    [1] "10bec8a8264f3eb59531e8c4c0303f96"
    
    
    $Packages$vctrs
    $Packages$vctrs$Package
    [1] "vctrs"
    
    $Packages$vctrs$Version
    [1] "0.3.8"
    
    $Packages$vctrs$Source
    [1] "Repository"
    
    $Packages$vctrs$Depends
    [1] "R (>= 3.3)"
    
    $Packages$vctrs$Imports
    [1] "ellipsis (>= 0.2.0)" "glue"                "rlang (>= 0.4.10)"  
    
    $Packages$vctrs$Repository
    [1] "CRAN"
    
    $Packages$vctrs$Hash
    [1] "ecf749a1b39ea72bd9b51b76292261f1"
    
    
    $Packages$castor
    $Packages$castor$Package
    [1] "castor"
    
    $Packages$castor$Version
    [1] "1.7.2"
    
    $Packages$castor$Source
    [1] "Repository"
    
    $Packages$castor$Depends
    [1] "Rcpp (>= 0.12.10)"
    
    $Packages$castor$Imports
    [1] "parallel"    "naturalsort" "stats"       "Matrix"      "RSpectra"   
    
    $Packages$castor$LinkingTo
    [1] "Rcpp"
    
    $Packages$castor$Repository
    [1] "CRAN"
    
    $Packages$castor$Hash
    [1] "0280b6b1793e9b6c98b727e967b4dd37"
    
    
    $Packages$mgcv
    $Packages$mgcv$Package
    [1] "mgcv"
    
    $Packages$mgcv$Version
    [1] "1.8-39"
    
    $Packages$mgcv$Source
    [1] "Repository"
    
    $Packages$mgcv$Depends
    [1] "R (>= 3.6.0)"     "nlme (>= 3.1-64)"
    
    $Packages$mgcv$Imports
    [1] "methods"  "stats"    "graphics" "Matrix"   "splines"  "utils"   
    
    $Packages$mgcv$Repository
    [1] "CRAN"
    
    $Packages$mgcv$Hash
    [1] "055265005c238024e306fe0b600c89ff"
    
    
    $Packages$blob
    $Packages$blob$Package
    [1] "blob"
    
    $Packages$blob$Version
    [1] "1.2.2"
    
    $Packages$blob$Source
    [1] "Repository"
    
    $Packages$blob$Imports
    [1] "methods"          "rlang"            "vctrs (>= 0.2.1)"
    
    $Packages$blob$Repository
    [1] "CRAN"
    
    $Packages$blob$Hash
    [1] "dc5f7a6598bb025d20d66bb758f12879"
    
    
    $Packages$survival
    $Packages$survival$Package
    [1] "survival"
    
    $Packages$survival$Version
    [1] "3.3-1"
    
    $Packages$survival$Source
    [1] "Repository"
    
    $Packages$survival$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$survival$Imports
    [1] "graphics" "Matrix"   "methods"  "splines"  "stats"    "utils"   
    
    $Packages$survival$Repository
    [1] "CRAN"
    
    $Packages$survival$Hash
    [1] "f6189c70451d3d68e0d571235576e833"
    
    
    $Packages$prodlim
    $Packages$prodlim$Package
    [1] "prodlim"
    
    $Packages$prodlim$Version
    [1] "2019.11.13"
    
    $Packages$prodlim$Source
    [1] "Repository"
    
    $Packages$prodlim$Depends
    [1] "R (>= 2.9.0)"
    
    $Packages$prodlim$Imports
    [1] "Rcpp (>= 0.11.5)" "stats"            "grDevices"        "graphics"        
    [5] "survival"         "KernSmooth"       "lava"            
    
    $Packages$prodlim$LinkingTo
    [1] "Rcpp"
    
    $Packages$prodlim$Repository
    [1] "CRAN"
    
    $Packages$prodlim$Hash
    [1] "c243bf70db3a6631a0c8783152fb7db9"
    
    
    $Packages$nloptr
    $Packages$nloptr$Package
    [1] "nloptr"
    
    $Packages$nloptr$Version
    [1] "2.0.0"
    
    $Packages$nloptr$Source
    [1] "Repository"
    
    $Packages$nloptr$LinkingTo
    [1] "testthat"
    
    $Packages$nloptr$Repository
    [1] "CRAN"
    
    $Packages$nloptr$Hash
    [1] "6e45c045fea34a9d0d1ceaa6fb7c4e91"
    
    
    $Packages$later
    $Packages$later$Package
    [1] "later"
    
    $Packages$later$Version
    [1] "1.3.0"
    
    $Packages$later$Source
    [1] "Repository"
    
    $Packages$later$Imports
    [1] "Rcpp (>= 0.12.9)" "rlang"           
    
    $Packages$later$LinkingTo
    [1] "Rcpp"
    
    $Packages$later$Repository
    [1] "CRAN"
    
    $Packages$later$Hash
    [1] "7e7b457d7766bc47f2a5f21cc2984f8e"
    
    
    $Packages$DBI
    $Packages$DBI$Package
    [1] "DBI"
    
    $Packages$DBI$Version
    [1] "1.1.2"
    
    $Packages$DBI$Source
    [1] "Repository"
    
    $Packages$DBI$Depends
    [1] "methods"      "R (>= 3.0.0)"
    
    $Packages$DBI$Repository
    [1] "CRAN"
    
    $Packages$DBI$Hash
    [1] "dcd1743af4336156873e3ce3c950b8b9"
    
    
    $Packages$R.utils
    $Packages$R.utils$Package
    [1] "R.utils"
    
    $Packages$R.utils$Version
    [1] "2.11.0"
    
    $Packages$R.utils$Source
    [1] "Repository"
    
    $Packages$R.utils$Depends
    [1] "R (>= 2.14.0)"    "R.oo (>= 1.24.0)"
    
    $Packages$R.utils$Imports
    [1] "methods"                "utils"                  "tools"                 
    [4] "R.methodsS3 (>= 1.8.1)"
    
    $Packages$R.utils$Repository
    [1] "CRAN"
    
    $Packages$R.utils$Hash
    [1] "a7ecb8e60815c7a18648e84cd121b23a"
    
    
    $Packages$tarchetypes
    $Packages$tarchetypes$Package
    [1] "tarchetypes"
    
    $Packages$tarchetypes$Version
    [1] "0.4.1"
    
    $Packages$tarchetypes$Source
    [1] "Repository"
    
    $Packages$tarchetypes$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$tarchetypes$Imports
     [1] "digest (>= 0.6.25)"    "dplyr (>= 1.0.0)"      "fs (>= 1.4.2)"        
     [4] "rlang (>= 0.4.7)"      "targets (>= 0.6.0)"    "tibble (>= 3.0.1)"    
     [7] "tidyselect (>= 1.1.0)" "utils"                 "vctrs (>= 0.3.4)"     
    [10] "withr (>= 2.1.2)"     
    
    $Packages$tarchetypes$Repository
    [1] "CRAN"
    
    $Packages$tarchetypes$Hash
    [1] "24c9eb58a70cef665befbcc9891572de"
    
    
    $Packages$SingleCellExperiment
    $Packages$SingleCellExperiment$Package
    [1] "SingleCellExperiment"
    
    $Packages$SingleCellExperiment$Version
    [1] "1.16.0"
    
    $Packages$SingleCellExperiment$Source
    [1] "Bioconductor"
    
    $Packages$SingleCellExperiment$Depends
    [1] "SummarizedExperiment"
    
    $Packages$SingleCellExperiment$Imports
    [1] "methods"       "utils"         "stats"         "S4Vectors"    
    [5] "BiocGenerics"  "GenomicRanges" "DelayedArray" 
    
    $Packages$SingleCellExperiment$git_url
    [1] "https://git.bioconductor.org/packages/SingleCellExperiment"
    
    $Packages$SingleCellExperiment$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$SingleCellExperiment$git_last_commit
    [1] "bb27609"
    
    $Packages$SingleCellExperiment$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$SingleCellExperiment$Hash
    [1] "27d907d4865fa44f6f17f757352ae7a3"
    
    
    $Packages$rappdirs
    $Packages$rappdirs$Package
    [1] "rappdirs"
    
    $Packages$rappdirs$Version
    [1] "0.3.3"
    
    $Packages$rappdirs$Source
    [1] "Repository"
    
    $Packages$rappdirs$Depends
    [1] "R (>= 3.2)"
    
    $Packages$rappdirs$Repository
    [1] "CRAN"
    
    $Packages$rappdirs$Hash
    [1] "5e3c5dc0b071b21fa128676560dbe94d"
    
    
    $Packages$selectr
    $Packages$selectr$Package
    [1] "selectr"
    
    $Packages$selectr$Version
    [1] "0.4-2"
    
    $Packages$selectr$Source
    [1] "Repository"
    
    $Packages$selectr$Depends
    [1] "R (>= 3.0)"
    
    $Packages$selectr$Imports
    [1] "methods" "stringr" "R6"     
    
    $Packages$selectr$Repository
    [1] "CRAN"
    
    $Packages$selectr$Hash
    [1] "3838071b66e0c566d55cc26bd6e27bf4"
    
    
    $Packages$jpeg
    $Packages$jpeg$Package
    [1] "jpeg"
    
    $Packages$jpeg$Version
    [1] "0.1-9"
    
    $Packages$jpeg$Source
    [1] "Repository"
    
    $Packages$jpeg$Depends
    [1] "R (>= 2.9.0)"
    
    $Packages$jpeg$Repository
    [1] "CRAN"
    
    $Packages$jpeg$Hash
    [1] "441ee36360a57b363f4fa3df0c364630"
    
    
    $Packages$zlibbioc
    $Packages$zlibbioc$Package
    [1] "zlibbioc"
    
    $Packages$zlibbioc$Version
    [1] "1.40.0"
    
    $Packages$zlibbioc$Source
    [1] "Bioconductor"
    
    $Packages$zlibbioc$git_url
    [1] "https://git.bioconductor.org/packages/zlibbioc"
    
    $Packages$zlibbioc$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$zlibbioc$git_last_commit
    [1] "3f116b3"
    
    $Packages$zlibbioc$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$zlibbioc$Hash
    [1] "3598bf6c766b0d2c15d1eefefa26d9ef"
    
    
    $Packages$fontawesome
    $Packages$fontawesome$Package
    [1] "fontawesome"
    
    $Packages$fontawesome$Version
    [1] "0.2.2"
    
    $Packages$fontawesome$Source
    [1] "Repository"
    
    $Packages$fontawesome$Depends
    [1] "R (>= 3.3.0)"
    
    $Packages$fontawesome$Imports
    [1] "rlang (>= 0.4.10)"      "htmltools (>= 0.5.1.1)"
    
    $Packages$fontawesome$Repository
    [1] "CRAN"
    
    $Packages$fontawesome$Hash
    [1] "55624ed409e46c5f358b2c060be87f67"
    
    
    $Packages$htmlwidgets
    $Packages$htmlwidgets$Package
    [1] "htmlwidgets"
    
    $Packages$htmlwidgets$Version
    [1] "1.5.4"
    
    $Packages$htmlwidgets$Source
    [1] "Repository"
    
    $Packages$htmlwidgets$Imports
    [1] "grDevices"            "htmltools (>= 0.3)"   "jsonlite (>= 0.9.16)"
    [4] "yaml"                
    
    $Packages$htmlwidgets$Repository
    [1] "CRAN"
    
    $Packages$htmlwidgets$Hash
    [1] "76147821cd3fcd8c4b04e1ef0498e7fb"
    
    
    $Packages$future
    $Packages$future$Package
    [1] "future"
    
    $Packages$future$Version
    [1] "1.24.0"
    
    $Packages$future$Source
    [1] "Repository"
    
    $Packages$future$Imports
    [1] "digest"                 "globals (>= 0.14.0)"    "listenv (>= 0.8.0)"    
    [4] "parallel"               "parallelly (>= 1.30.0)" "tools"                 
    [7] "utils"                 
    
    $Packages$future$Repository
    [1] "CRAN"
    
    $Packages$future$Hash
    [1] "5cc7addaa73372fbee0a7d06c880068e"
    
    
    $Packages$ANCOMBC
    $Packages$ANCOMBC$Package
    [1] "ANCOMBC"
    
    $Packages$ANCOMBC$Version
    [1] "1.4.0"
    
    $Packages$ANCOMBC$Source
    [1] "Bioconductor"
    
    $Packages$ANCOMBC$Imports
    [1] "stats"      "MASS"       "nloptr"     "Rdpack"     "phyloseq"  
    [6] "microbiome"
    
    $Packages$ANCOMBC$git_url
    [1] "https://git.bioconductor.org/packages/ANCOMBC"
    
    $Packages$ANCOMBC$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$ANCOMBC$git_last_commit
    [1] "b9a7fb1"
    
    $Packages$ANCOMBC$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$ANCOMBC$Hash
    [1] "01d19305af91bf6cf081f8a9dfdb18d1"
    
    
    $Packages$scater
    $Packages$scater$Package
    [1] "scater"
    
    $Packages$scater$Version
    [1] "1.22.0"
    
    $Packages$scater$Source
    [1] "Bioconductor"
    
    $Packages$scater$Depends
    [1] "SingleCellExperiment" "scuttle"              "ggplot2"             
    
    $Packages$scater$Imports
     [1] "stats"                "utils"                "methods"             
     [4] "grid"                 "gridExtra"            "Matrix"              
     [7] "BiocGenerics"         "S4Vectors"            "SummarizedExperiment"
    [10] "DelayedArray"         "DelayedMatrixStats"   "beachmat"            
    [13] "BiocNeighbors"        "BiocSingular"         "BiocParallel"        
    [16] "rlang"                "ggbeeswarm"           "viridis"             
    [19] "Rtsne"                "RColorBrewer"         "ggrepel"             
    
    $Packages$scater$git_url
    [1] "https://git.bioconductor.org/packages/scater"
    
    $Packages$scater$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$scater$git_last_commit
    [1] "ea2c95c"
    
    $Packages$scater$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$scater$Hash
    [1] "8e229c5126b386c658047881fe2e5299"
    
    
    $Packages$irlba
    $Packages$irlba$Package
    [1] "irlba"
    
    $Packages$irlba$Version
    [1] "2.3.5"
    
    $Packages$irlba$Source
    [1] "Repository"
    
    $Packages$irlba$Depends
    [1] "R (>= 3.6.2)" "Matrix"      
    
    $Packages$irlba$Imports
    [1] "stats"   "methods"
    
    $Packages$irlba$LinkingTo
    [1] "Matrix"
    
    $Packages$irlba$Repository
    [1] "CRAN"
    
    $Packages$irlba$Hash
    [1] "066c11bb9bc75b343f3de1ecaf3b7ba2"
    
    
    $Packages$BiocSet
    $Packages$BiocSet$Package
    [1] "BiocSet"
    
    $Packages$BiocSet$Version
    [1] "1.8.1"
    
    $Packages$BiocSet$Source
    [1] "Bioconductor"
    
    $Packages$BiocSet$Depends
    [1] "R (>= 3.6)" "dplyr"     
    
    $Packages$BiocSet$Imports
     [1] "methods"       "tibble"        "utils"         "rlang"        
     [5] "plyr"          "S4Vectors"     "BiocIO"        "AnnotationDbi"
     [9] "KEGGREST"      "ontologyIndex" "tidyr"        
    
    $Packages$BiocSet$git_url
    [1] "https://git.bioconductor.org/packages/BiocSet"
    
    $Packages$BiocSet$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocSet$git_last_commit
    [1] "1f7d340"
    
    $Packages$BiocSet$git_last_commit_date
    [1] "2021-11-02"
    
    $Packages$BiocSet$Hash
    [1] "1f006a9fb0988ff10b02c4378a67b78b"
    
    
    $Packages$tidymodels
    $Packages$tidymodels$Package
    [1] "tidymodels"
    
    $Packages$tidymodels$Version
    [1] "0.1.4"
    
    $Packages$tidymodels$Source
    [1] "Repository"
    
    $Packages$tidymodels$Depends
    [1] "R (>= 3.1)"
    
    $Packages$tidymodels$Imports
     [1] "broom (>= 0.7.9)"        "cli (>= 3.0.1)"         
     [3] "conflicted (>= 1.0.4)"   "dials (>= 0.0.10)"      
     [5] "dplyr (>= 1.0.7)"        "ggplot2 (>= 3.3.5)"     
     [7] "hardhat (>= 0.1.6)"      "infer (>= 1.0.0)"       
     [9] "modeldata (>= 0.1.1)"    "parsnip (>= 0.1.7)"     
    [11] "purrr (>= 0.3.4)"        "recipes (>= 0.1.17)"    
    [13] "rlang (>= 0.4.11)"       "rsample (>= 0.1.0)"     
    [15] "rstudioapi (>= 0.13)"    "tibble (>= 3.1.5)"      
    [17] "tidyr (>= 1.1.4)"        "tune (>= 0.1.6)"        
    [19] "workflows (>= 0.2.3)"    "workflowsets (>= 0.1.0)"
    [21] "yardstick (>= 0.0.8)"   
    
    $Packages$tidymodels$Repository
    [1] "CRAN"
    
    $Packages$tidymodels$Hash
    [1] "74cfd3e1de860e43453b7b36675f6480"
    
    
    $Packages$DEoptimR
    $Packages$DEoptimR$Package
    [1] "DEoptimR"
    
    $Packages$DEoptimR$Version
    [1] "1.0-10"
    
    $Packages$DEoptimR$Source
    [1] "Repository"
    
    $Packages$DEoptimR$Imports
    [1] "stats"
    
    $Packages$DEoptimR$Repository
    [1] "CRAN"
    
    $Packages$DEoptimR$Hash
    [1] "455af0b32e7e45047cf7f495d462e456"
    
    
    $Packages$Rcpp
    $Packages$Rcpp$Package
    [1] "Rcpp"
    
    $Packages$Rcpp$Version
    [1] "1.0.8"
    
    $Packages$Rcpp$Source
    [1] "Repository"
    
    $Packages$Rcpp$Imports
    [1] "methods" "utils"  
    
    $Packages$Rcpp$Repository
    [1] "CRAN"
    
    $Packages$Rcpp$Hash
    [1] "22b546dd7e337f6c0c58a39983a496bc"
    
    
    $Packages$readr
    $Packages$readr$Package
    [1] "readr"
    
    $Packages$readr$Version
    [1] "2.1.2"
    
    $Packages$readr$Source
    [1] "Repository"
    
    $Packages$readr$Depends
    [1] "R (>= 3.1)"
    
    $Packages$readr$Imports
     [1] "cli (>= 3.0.0)"       "clipr"                "crayon"              
     [4] "hms (>= 0.4.1)"       "lifecycle (>= 0.2.0)" "methods"             
     [7] "R6"                   "rlang"                "tibble"              
    [10] "utils"                "vroom (>= 1.5.6)"    
    
    $Packages$readr$LinkingTo
    [1] "cpp11"           "tzdb (>= 0.1.1)"
    
    $Packages$readr$Repository
    [1] "CRAN"
    
    $Packages$readr$Hash
    [1] "9c59de1357dc209868b5feb5c9f0fe2f"
    
    
    $Packages$KernSmooth
    $Packages$KernSmooth$Package
    [1] "KernSmooth"
    
    $Packages$KernSmooth$Version
    [1] "2.23-20"
    
    $Packages$KernSmooth$Source
    [1] "Repository"
    
    $Packages$KernSmooth$Depends
    [1] "R (>= 2.5.0)" "stats"       
    
    $Packages$KernSmooth$Repository
    [1] "CRAN"
    
    $Packages$KernSmooth$Hash
    [1] "8dcfa99b14c296bc9f1fd64d52fd3ce7"
    
    
    $Packages$DT
    $Packages$DT$Package
    [1] "DT"
    
    $Packages$DT$Version
    [1] "0.21"
    
    $Packages$DT$Source
    [1] "Repository"
    
    $Packages$DT$Imports
    [1] "htmltools (>= 0.3.6)" "htmlwidgets (>= 1.3)" "jsonlite (>= 0.9.16)"
    [4] "magrittr"             "crosstalk"            "jquerylib"           
    [7] "promises"            
    
    $Packages$DT$Repository
    [1] "CRAN"
    
    $Packages$DT$Hash
    [1] "45fa28dbf288cd606e13ca35d3d72437"
    
    
    $Packages$promises
    $Packages$promises$Package
    [1] "promises"
    
    $Packages$promises$Version
    [1] "1.2.0.1"
    
    $Packages$promises$Source
    [1] "Repository"
    
    $Packages$promises$Imports
    [1] "R6"       "Rcpp"     "later"    "rlang"    "stats"    "magrittr"
    
    $Packages$promises$LinkingTo
    [1] "later" "Rcpp" 
    
    $Packages$promises$Repository
    [1] "CRAN"
    
    $Packages$promises$Hash
    [1] "4ab2c43adb4d4699cf3690acd378d75d"
    
    
    $Packages$DelayedArray
    $Packages$DelayedArray$Package
    [1] "DelayedArray"
    
    $Packages$DelayedArray$Version
    [1] "0.20.0"
    
    $Packages$DelayedArray$Source
    [1] "Bioconductor"
    
    $Packages$DelayedArray$Depends
    [1] "R (>= 4.0.0)"              "methods"                  
    [3] "stats4"                    "Matrix"                   
    [5] "BiocGenerics (>= 0.37.0)"  "MatrixGenerics (>= 1.1.3)"
    [7] "S4Vectors (>= 0.27.2)"     "IRanges (>= 2.17.3)"      
    
    $Packages$DelayedArray$Imports
    [1] "stats"
    
    $Packages$DelayedArray$LinkingTo
    [1] "S4Vectors"
    
    $Packages$DelayedArray$git_url
    [1] "https://git.bioconductor.org/packages/DelayedArray"
    
    $Packages$DelayedArray$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$DelayedArray$git_last_commit
    [1] "829b529"
    
    $Packages$DelayedArray$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$DelayedArray$Hash
    [1] "8400bc4ac5cf78a44eefa2395a2ed782"
    
    
    $Packages$pkgload
    $Packages$pkgload$Package
    [1] "pkgload"
    
    $Packages$pkgload$Version
    [1] "1.2.4"
    
    $Packages$pkgload$Source
    [1] "Repository"
    
    $Packages$pkgload$Imports
    [1] "cli"        "crayon"     "desc"       "methods"    "rlang"     
    [6] "rprojroot"  "rstudioapi" "utils"      "withr"     
    
    $Packages$pkgload$Repository
    [1] "CRAN"
    
    $Packages$pkgload$Hash
    [1] "7533cd805940821bf23eaf3c8d4c1735"
    
    
    $Packages$vegan
    $Packages$vegan$Package
    [1] "vegan"
    
    $Packages$vegan$Version
    [1] "2.5-7"
    
    $Packages$vegan$Source
    [1] "Repository"
    
    $Packages$vegan$Depends
    [1] "permute (>= 0.9-0)" "lattice"            "R (>= 3.4.0)"      
    
    $Packages$vegan$Imports
    [1] "MASS"    "cluster" "mgcv"   
    
    $Packages$vegan$Repository
    [1] "CRAN"
    
    $Packages$vegan$Hash
    [1] "01771d8de354fa30c87c68e599e2429b"
    
    
    $Packages$graph
    $Packages$graph$Package
    [1] "graph"
    
    $Packages$graph$Version
    [1] "1.72.0"
    
    $Packages$graph$Source
    [1] "Bioconductor"
    
    $Packages$graph$Depends
    [1] "R (>= 2.10)"               "methods"                  
    [3] "BiocGenerics (>= 0.13.11)"
    
    $Packages$graph$Imports
    [1] "stats"  "stats4" "utils" 
    
    $Packages$graph$git_url
    [1] "https://git.bioconductor.org/packages/graph"
    
    $Packages$graph$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$graph$git_last_commit
    [1] "7afbd26"
    
    $Packages$graph$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$graph$Hash
    [1] "008757628d8e6bacfe8f122954d8f024"
    
    
    $Packages$RcppParallel
    $Packages$RcppParallel$Package
    [1] "RcppParallel"
    
    $Packages$RcppParallel$Version
    [1] "5.1.5"
    
    $Packages$RcppParallel$Source
    [1] "Repository"
    
    $Packages$RcppParallel$Depends
    [1] "R (>= 3.0.2)"
    
    $Packages$RcppParallel$Repository
    [1] "CRAN"
    
    $Packages$RcppParallel$Hash
    [1] "f3e94e34ff656a7c8336ce01207bc2b8"
    
    
    $Packages$ShortRead
    $Packages$ShortRead$Package
    [1] "ShortRead"
    
    $Packages$ShortRead$Version
    [1] "1.52.0"
    
    $Packages$ShortRead$Source
    [1] "Bioconductor"
    
    $Packages$ShortRead$Depends
    [1] "BiocGenerics (>= 0.23.3)"      "BiocParallel"                 
    [3] "Biostrings (>= 2.47.6)"        "Rsamtools (>= 1.31.2)"        
    [5] "GenomicAlignments (>= 1.15.6)"
    
    $Packages$ShortRead$Imports
     [1] "Biobase"                   "S4Vectors (>= 0.17.25)"   
     [3] "IRanges (>= 2.13.12)"      "GenomeInfoDb (>= 1.15.2)" 
     [5] "GenomicRanges (>= 1.31.8)" "hwriter"                  
     [7] "methods"                   "zlibbioc"                 
     [9] "lattice"                   "latticeExtra"             
    
    $Packages$ShortRead$LinkingTo
    [1] "S4Vectors"  "IRanges"    "XVector"    "Biostrings" "Rhtslib"   
    [6] "zlibbioc"  
    
    $Packages$ShortRead$git_url
    [1] "https://git.bioconductor.org/packages/ShortRead"
    
    $Packages$ShortRead$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$ShortRead$git_last_commit
    [1] "4d7304d"
    
    $Packages$ShortRead$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$ShortRead$Hash
    [1] "cc412b33be7046db49d38794133ddfc4"
    
    
    $Packages$RSpectra
    $Packages$RSpectra$Package
    [1] "RSpectra"
    
    $Packages$RSpectra$Version
    [1] "0.16-0"
    
    $Packages$RSpectra$Source
    [1] "Repository"
    
    $Packages$RSpectra$Depends
    [1] "R (>= 3.0.2)"
    
    $Packages$RSpectra$Imports
    [1] "Matrix (>= 1.1-0)" "Rcpp (>= 0.11.5)" 
    
    $Packages$RSpectra$LinkingTo
    [1] "Rcpp"                     "RcppEigen (>= 0.3.3.3.0)"
    
    $Packages$RSpectra$Repository
    [1] "CRAN"
    
    $Packages$RSpectra$Hash
    [1] "a41329d24d5a98eaed2bd0159adb1b5f"
    
    
    $Packages$fs
    $Packages$fs$Package
    [1] "fs"
    
    $Packages$fs$Version
    [1] "1.5.2"
    
    $Packages$fs$Source
    [1] "Repository"
    
    $Packages$fs$Depends
    [1] "R (>= 3.1)"
    
    $Packages$fs$Imports
    [1] "methods"
    
    $Packages$fs$Repository
    [1] "CRAN"
    
    $Packages$fs$Hash
    [1] "7c89603d81793f0d5486d91ab1fc6f1d"
    
    
    $Packages$brio
    $Packages$brio$Package
    [1] "brio"
    
    $Packages$brio$Version
    [1] "1.1.3"
    
    $Packages$brio$Source
    [1] "Repository"
    
    $Packages$brio$Repository
    [1] "CRAN"
    
    $Packages$brio$Hash
    [1] "976cf154dfb043c012d87cddd8bca363"
    
    
    $Packages$fastmatch
    $Packages$fastmatch$Package
    [1] "fastmatch"
    
    $Packages$fastmatch$Version
    [1] "1.1-3"
    
    $Packages$fastmatch$Source
    [1] "Repository"
    
    $Packages$fastmatch$Depends
    [1] "R (>= 2.3.0)"
    
    $Packages$fastmatch$Repository
    [1] "CRAN"
    
    $Packages$fastmatch$Hash
    [1] "dabc225759a2c2b241e60e42bf0e8e54"
    
    
    $Packages$conflicted
    $Packages$conflicted$Package
    [1] "conflicted"
    
    $Packages$conflicted$Version
    [1] "1.1.0"
    
    $Packages$conflicted$Source
    [1] "Repository"
    
    $Packages$conflicted$Depends
    [1] "R (>= 3.2)"
    
    $Packages$conflicted$Imports
    [1] "memoise"          "rlang (>= 0.3.4)"
    
    $Packages$conflicted$Repository
    [1] "CRAN"
    
    $Packages$conflicted$Hash
    [1] "c6bb5e1ef58f2f1c84f238f55bd2e56a"
    
    
    $Packages$digest
    $Packages$digest$Package
    [1] "digest"
    
    $Packages$digest$Version
    [1] "0.6.29"
    
    $Packages$digest$Source
    [1] "Repository"
    
    $Packages$digest$Depends
    [1] "R (>= 3.3.0)"
    
    $Packages$digest$Imports
    [1] "utils"
    
    $Packages$digest$Repository
    [1] "CRAN"
    
    $Packages$digest$Hash
    [1] "cf6b206a045a684728c3267ef7596190"
    
    
    $Packages$png
    $Packages$png$Package
    [1] "png"
    
    $Packages$png$Version
    [1] "0.1-7"
    
    $Packages$png$Source
    [1] "Repository"
    
    $Packages$png$Depends
    [1] "R (>= 2.9.0)"
    
    $Packages$png$Repository
    [1] "CRAN"
    
    $Packages$png$Hash
    [1] "03b7076c234cb3331288919983326c55"
    
    
    $Packages$mia
    $Packages$mia$Package
    [1] "mia"
    
    $Packages$mia$Version
    [1] "1.2.7"
    
    $Packages$mia$Source
    [1] "Bioconductor"
    
    $Packages$mia$Depends
    [1] "R (>= 4.0)"                          
    [2] "SummarizedExperiment"                
    [3] "SingleCellExperiment"                
    [4] "TreeSummarizedExperiment (>= 1.99.3)"
    [5] "MultiAssayExperiment"                
    
    $Packages$mia$Imports
     [1] "methods"              "stats"                "utils"               
     [4] "MASS"                 "ape"                  "decontam"            
     [7] "vegan"                "BiocGenerics"         "S4Vectors"           
    [10] "IRanges"              "Biostrings"           "DECIPHER"            
    [13] "BiocParallel"         "DelayedArray"         "DelayedMatrixStats"  
    [16] "scuttle"              "scater"               "DirichletMultinomial"
    [19] "rlang"                "dplyr"                "tibble"              
    [22] "tidyr"               
    
    $Packages$mia$git_url
    [1] "https://git.bioconductor.org/packages/mia"
    
    $Packages$mia$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$mia$git_last_commit
    [1] "217dc43"
    
    $Packages$mia$git_last_commit_date
    [1] "2022-02-06"
    
    $Packages$mia$Hash
    [1] "824cd6adeffd1bad73f2b6ed19f3967c"
    
    
    $Packages$gh
    $Packages$gh$Package
    [1] "gh"
    
    $Packages$gh$Version
    [1] "1.3.0"
    
    $Packages$gh$Source
    [1] "Repository"
    
    $Packages$gh$Imports
    [1] "cli (>= 2.0.1)" "gitcreds"       "httr (>= 1.2)"  "ini"           
    [5] "jsonlite"      
    
    $Packages$gh$Repository
    [1] "CRAN"
    
    $Packages$gh$Hash
    [1] "38c2580abbda249bd6afeec00d14f531"
    
    
    $Packages$Rhtslib
    $Packages$Rhtslib$Package
    [1] "Rhtslib"
    
    $Packages$Rhtslib$Version
    [1] "1.26.0"
    
    $Packages$Rhtslib$Source
    [1] "Bioconductor"
    
    $Packages$Rhtslib$Imports
    [1] "zlibbioc"
    
    $Packages$Rhtslib$LinkingTo
    [1] "zlibbioc"
    
    $Packages$Rhtslib$git_url
    [1] "https://git.bioconductor.org/packages/Rhtslib"
    
    $Packages$Rhtslib$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$Rhtslib$git_last_commit
    [1] "f5b20e9"
    
    $Packages$Rhtslib$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$Rhtslib$Hash
    [1] "b88f760ac9993107b94388100f136f43"
    
    
    $Packages$here
    $Packages$here$Package
    [1] "here"
    
    $Packages$here$Version
    [1] "1.0.1"
    
    $Packages$here$Source
    [1] "Repository"
    
    $Packages$here$Imports
    [1] "rprojroot (>= 2.0.2)"
    
    $Packages$here$Repository
    [1] "CRAN"
    
    $Packages$here$Hash
    [1] "24b224366f9c2e7534d2344d10d59211"
    
    
    $Packages$janeaustenr
    $Packages$janeaustenr$Package
    [1] "janeaustenr"
    
    $Packages$janeaustenr$Version
    [1] "0.1.5"
    
    $Packages$janeaustenr$Source
    [1] "Repository"
    
    $Packages$janeaustenr$Depends
    [1] "R (>= 3.1.2)"
    
    $Packages$janeaustenr$Repository
    [1] "CRAN"
    
    $Packages$janeaustenr$Hash
    [1] "8b07a4b9d0a0d97d9fe12de8af6d219e"
    
    
    $Packages$pkgconfig
    $Packages$pkgconfig$Package
    [1] "pkgconfig"
    
    $Packages$pkgconfig$Version
    [1] "2.0.3"
    
    $Packages$pkgconfig$Source
    [1] "Repository"
    
    $Packages$pkgconfig$Imports
    [1] "utils"
    
    $Packages$pkgconfig$Repository
    [1] "CRAN"
    
    $Packages$pkgconfig$Hash
    [1] "01f28d4278f15c76cddbea05899c5d6f"
    
    
    $Packages$modeldata
    $Packages$modeldata$Package
    [1] "modeldata"
    
    $Packages$modeldata$Version
    [1] "0.1.1"
    
    $Packages$modeldata$Source
    [1] "Repository"
    
    $Packages$modeldata$Depends
    [1] "R (>= 2.10)"
    
    $Packages$modeldata$Repository
    [1] "CRAN"
    
    $Packages$modeldata$Hash
    [1] "30951e1c2022c6a8e1809d4224c0a5ab"
    
    
    $Packages$DelayedMatrixStats
    $Packages$DelayedMatrixStats$Package
    [1] "DelayedMatrixStats"
    
    $Packages$DelayedMatrixStats$Version
    [1] "1.16.0"
    
    $Packages$DelayedMatrixStats$Source
    [1] "Bioconductor"
    
    $Packages$DelayedMatrixStats$Depends
    [1] "MatrixGenerics (>= 1.5.3)" "DelayedArray (>= 0.17.6)" 
    
    $Packages$DelayedMatrixStats$Imports
    [1] "methods"                 "matrixStats (>= 0.60.0)"
    [3] "sparseMatrixStats"       "Matrix"                 
    [5] "S4Vectors (>= 0.17.5)"   "IRanges (>= 2.25.10)"   
    
    $Packages$DelayedMatrixStats$git_url
    [1] "https://git.bioconductor.org/packages/DelayedMatrixStats"
    
    $Packages$DelayedMatrixStats$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$DelayedMatrixStats$git_last_commit
    [1] "d44a3d7"
    
    $Packages$DelayedMatrixStats$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$DelayedMatrixStats$Hash
    [1] "9ae660d2827845adeaa2fbdc3e6ad468"
    
    
    $Packages$warp
    $Packages$warp$Package
    [1] "warp"
    
    $Packages$warp$Version
    [1] "0.2.0"
    
    $Packages$warp$Source
    [1] "Repository"
    
    $Packages$warp$Depends
    [1] "R (>= 3.2)"
    
    $Packages$warp$Repository
    [1] "CRAN"
    
    $Packages$warp$Hash
    [1] "2982481615756e24e79fee95bdc95daa"
    
    
    $Packages$gower
    $Packages$gower$Package
    [1] "gower"
    
    $Packages$gower$Version
    [1] "1.0.0"
    
    $Packages$gower$Source
    [1] "Repository"
    
    $Packages$gower$Repository
    [1] "CRAN"
    
    $Packages$gower$Hash
    [1] "e30deda901954a80e035ad2e972ba7fd"
    
    
    $Packages$ggbeeswarm
    $Packages$ggbeeswarm$Package
    [1] "ggbeeswarm"
    
    $Packages$ggbeeswarm$Version
    [1] "0.6.0"
    
    $Packages$ggbeeswarm$Source
    [1] "Repository"
    
    $Packages$ggbeeswarm$Depends
    [1] "R (>= 3.0.0)"     "ggplot2 (>= 2.0)"
    
    $Packages$ggbeeswarm$Imports
    [1] "beeswarm" "vipor"   
    
    $Packages$ggbeeswarm$Repository
    [1] "CRAN"
    
    $Packages$ggbeeswarm$Hash
    [1] "dd68b9b215b2d3119603549a794003c3"
    
    
    $Packages$gt
    $Packages$gt$Package
    [1] "gt"
    
    $Packages$gt$Version
    [1] "0.4.0"
    
    $Packages$gt$Source
    [1] "Repository"
    
    $Packages$gt$Depends
    [1] "R (>= 3.2.0)"
    
    $Packages$gt$Imports
     [1] "base64enc (>= 0.1-3)"  "bitops (>= 1.0.6)"     "checkmate (>= 2.0.0)" 
     [4] "commonmark (>= 1.7)"   "dplyr (>= 1.0.8)"      "fs (>= 1.5.2)"        
     [7] "ggplot2 (>= 3.3.5)"    "glue (>= 1.6.1)"       "htmltools (>= 0.5.2)" 
    [10] "magrittr (>= 2.0.2)"   "rlang (>= 1.0.1)"      "sass (>= 0.4.0)"      
    [13] "scales (>= 1.1.1)"     "stringr (>= 1.4.0)"    "tibble (>= 3.1.6)"    
    [16] "tidyselect (>= 1.1.1)"
    
    $Packages$gt$Repository
    [1] "CRAN"
    
    $Packages$gt$Hash
    [1] "4905765870343c66704c9e495e4e3bc0"
    
    
    $Packages$yardstick
    $Packages$yardstick$Package
    [1] "yardstick"
    
    $Packages$yardstick$Version
    [1] "0.0.9"
    
    $Packages$yardstick$Source
    [1] "Repository"
    
    $Packages$yardstick$Depends
    [1] "R (>= 2.10)"
    
    $Packages$yardstick$Imports
    [1] "dplyr (>= 1.0.0)" "generics"         "pROC (>= 1.15.0)" "rlang (>= 0.4.2)"
    [5] "tidyselect"       "utils"            "vctrs (>= 0.3.6)"
    
    $Packages$yardstick$Repository
    [1] "CRAN"
    
    $Packages$yardstick$Hash
    [1] "fd1588dbcbb85aacd0d9b9009351c79d"
    
    
    $Packages$GPfit
    $Packages$GPfit$Package
    [1] "GPfit"
    
    $Packages$GPfit$Version
    [1] "1.0-8"
    
    $Packages$GPfit$Source
    [1] "Repository"
    
    $Packages$GPfit$Imports
    [1] "lhs (>= 0.5)"        "lattice (>= 0.18-8)"
    
    $Packages$GPfit$Repository
    [1] "CRAN"
    
    $Packages$GPfit$Hash
    [1] "29a7dccade1fd037c8262c2a239775eb"
    
    
    $Packages$iterators
    $Packages$iterators$Package
    [1] "iterators"
    
    $Packages$iterators$Version
    [1] "1.0.14"
    
    $Packages$iterators$Source
    [1] "Repository"
    
    $Packages$iterators$Depends
    [1] "R (>= 2.5.0)" "utils"       
    
    $Packages$iterators$Repository
    [1] "CRAN"
    
    $Packages$iterators$Hash
    [1] "8954069286b4b2b0d023d1b288dce978"
    
    
    $Packages$pixmap
    $Packages$pixmap$Package
    [1] "pixmap"
    
    $Packages$pixmap$Version
    [1] "0.4-12"
    
    $Packages$pixmap$Source
    [1] "Repository"
    
    $Packages$pixmap$Imports
    [1] "methods"   "graphics"  "grDevices"
    
    $Packages$pixmap$Repository
    [1] "CRAN"
    
    $Packages$pixmap$Hash
    [1] "ff0e0b97265ced0db4d13d0005334d78"
    
    
    $Packages$reticulate
    $Packages$reticulate$Package
    [1] "reticulate"
    
    $Packages$reticulate$Version
    [1] "1.24"
    
    $Packages$reticulate$Source
    [1] "Repository"
    
    $Packages$reticulate$Depends
    [1] "R (>= 3.0)"
    
    $Packages$reticulate$Imports
     [1] "Matrix"           "Rcpp (>= 0.12.7)" "RcppTOML"         "graphics"        
     [5] "here"             "jsonlite"         "methods"          "png"             
     [9] "rappdirs"         "utils"            "withr"           
    
    $Packages$reticulate$LinkingTo
    [1] "Rcpp"
    
    $Packages$reticulate$Repository
    [1] "CRAN"
    
    $Packages$reticulate$Hash
    [1] "ffdf27627a3c1537478073c43b6e7980"
    
    
    $Packages$SummarizedExperiment
    $Packages$SummarizedExperiment$Package
    [1] "SummarizedExperiment"
    
    $Packages$SummarizedExperiment$Version
    [1] "1.24.0"
    
    $Packages$SummarizedExperiment$Source
    [1] "Bioconductor"
    
    $Packages$SummarizedExperiment$Depends
    [1] "R (>= 4.0.0)"              "methods"                  
    [3] "MatrixGenerics (>= 1.1.3)" "GenomicRanges (>= 1.41.5)"
    [5] "Biobase"                  
    
    $Packages$SummarizedExperiment$Imports
    [1] "utils"                     "stats"                    
    [3] "tools"                     "Matrix"                   
    [5] "BiocGenerics (>= 0.37.0)"  "S4Vectors (>= 0.27.12)"   
    [7] "IRanges (>= 2.23.9)"       "GenomeInfoDb (>= 1.13.1)" 
    [9] "DelayedArray (>= 0.15.10)"
    
    $Packages$SummarizedExperiment$git_url
    [1] "https://git.bioconductor.org/packages/SummarizedExperiment"
    
    $Packages$SummarizedExperiment$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$SummarizedExperiment$git_last_commit
    [1] "d37f193"
    
    $Packages$SummarizedExperiment$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$SummarizedExperiment$Hash
    [1] "619feb3635b27a198477316a033b1ea8"
    
    
    $Packages$beeswarm
    $Packages$beeswarm$Package
    [1] "beeswarm"
    
    $Packages$beeswarm$Version
    [1] "0.4.0"
    
    $Packages$beeswarm$Source
    [1] "Repository"
    
    $Packages$beeswarm$Imports
    [1] "stats"     "graphics"  "grDevices" "utils"    
    
    $Packages$beeswarm$Repository
    [1] "CRAN"
    
    $Packages$beeswarm$Hash
    [1] "0f4e9d8caa6feaa7e409ae6c30f2ca66"
    
    
    $Packages$xfun
    $Packages$xfun$Package
    [1] "xfun"
    
    $Packages$xfun$Version
    [1] "0.30"
    
    $Packages$xfun$Source
    [1] "Repository"
    
    $Packages$xfun$Imports
    [1] "stats" "tools"
    
    $Packages$xfun$Repository
    [1] "CRAN"
    
    $Packages$xfun$Hash
    [1] "e83f48136b041845e50a6658feffb197"
    
    
    $Packages$bslib
    $Packages$bslib$Package
    [1] "bslib"
    
    $Packages$bslib$Version
    [1] "0.3.1"
    
    $Packages$bslib$Source
    [1] "Repository"
    
    $Packages$bslib$Depends
    [1] "R (>= 2.10)"
    
    $Packages$bslib$Imports
    [1] "grDevices"            "htmltools (>= 0.5.2)" "jsonlite"            
    [4] "sass (>= 0.4.0)"      "jquerylib (>= 0.1.3)" "rlang"               
    
    $Packages$bslib$Repository
    [1] "CRAN"
    
    $Packages$bslib$Hash
    [1] "56ae7e1987b340186a8a5a157c2ec358"
    
    
    $Packages$tidyselect
    $Packages$tidyselect$Package
    [1] "tidyselect"
    
    $Packages$tidyselect$Version
    [1] "1.1.2"
    
    $Packages$tidyselect$Source
    [1] "Repository"
    
    $Packages$tidyselect$Depends
    [1] "R (>= 3.2)"
    
    $Packages$tidyselect$Imports
    [1] "ellipsis"         "glue (>= 1.3.0)"  "purrr (>= 0.3.2)" "rlang (>= 1.0.1)"
    [5] "vctrs (>= 0.3.0)"
    
    $Packages$tidyselect$Repository
    [1] "CRAN"
    
    $Packages$tidyselect$Hash
    [1] "17f6da8cfd7002760a859915ce7eef8f"
    
    
    $Packages$reshape2
    $Packages$reshape2$Package
    [1] "reshape2"
    
    $Packages$reshape2$Version
    [1] "1.4.4"
    
    $Packages$reshape2$Source
    [1] "Repository"
    
    $Packages$reshape2$Depends
    [1] "R (>= 3.1)"
    
    $Packages$reshape2$Imports
    [1] "plyr (>= 1.8.1)" "Rcpp"            "stringr"        
    
    $Packages$reshape2$LinkingTo
    [1] "Rcpp"
    
    $Packages$reshape2$Repository
    [1] "CRAN"
    
    $Packages$reshape2$Hash
    [1] "bb5996d0bd962d214a11140d77589917"
    
    
    $Packages$purrr
    $Packages$purrr$Package
    [1] "purrr"
    
    $Packages$purrr$Version
    [1] "0.3.4"
    
    $Packages$purrr$Source
    [1] "Repository"
    
    $Packages$purrr$Depends
    [1] "R (>= 3.2)"
    
    $Packages$purrr$Imports
    [1] "magrittr (>= 1.5)" "rlang (>= 0.3.1)" 
    
    $Packages$purrr$Repository
    [1] "CRAN"
    
    $Packages$purrr$Hash
    [1] "97def703420c8ab10d8f0e6c72101e02"
    
    
    $Packages$sourcetools
    $Packages$sourcetools$Package
    [1] "sourcetools"
    
    $Packages$sourcetools$Version
    [1] "0.1.7"
    
    $Packages$sourcetools$Source
    [1] "Repository"
    
    $Packages$sourcetools$Depends
    [1] "R (>= 3.0.2)"
    
    $Packages$sourcetools$Repository
    [1] "CRAN"
    
    $Packages$sourcetools$Hash
    [1] "947e4e02a79effa5d512473e10f41797"
    
    
    $Packages$viridisLite
    $Packages$viridisLite$Package
    [1] "viridisLite"
    
    $Packages$viridisLite$Version
    [1] "0.4.0"
    
    $Packages$viridisLite$Source
    [1] "Repository"
    
    $Packages$viridisLite$Depends
    [1] "R (>= 2.10)"
    
    $Packages$viridisLite$Repository
    [1] "CRAN"
    
    $Packages$viridisLite$Hash
    [1] "55e157e2aa88161bdb0754218470d204"
    
    
    $Packages$snow
    $Packages$snow$Package
    [1] "snow"
    
    $Packages$snow$Version
    [1] "0.4-4"
    
    $Packages$snow$Source
    [1] "Repository"
    
    $Packages$snow$Depends
    [1] "R (>= 2.13.1)" "utils"        
    
    $Packages$snow$Repository
    [1] "CRAN"
    
    $Packages$snow$Hash
    [1] "40b74690debd20c57d93d8c246b305d4"
    
    
    $Packages$rlang
    $Packages$rlang$Package
    [1] "rlang"
    
    $Packages$rlang$Version
    [1] "1.0.2"
    
    $Packages$rlang$Source
    [1] "Repository"
    
    $Packages$rlang$Depends
    [1] "R (>= 3.4.0)"
    
    $Packages$rlang$Imports
    [1] "utils"
    
    $Packages$rlang$Repository
    [1] "CRAN"
    
    $Packages$rlang$Hash
    [1] "04884d9a75d778aca22c7154b8333ec9"
    
    
    $Packages$jquerylib
    $Packages$jquerylib$Package
    [1] "jquerylib"
    
    $Packages$jquerylib$Version
    [1] "0.1.4"
    
    $Packages$jquerylib$Source
    [1] "Repository"
    
    $Packages$jquerylib$Imports
    [1] "htmltools"
    
    $Packages$jquerylib$Repository
    [1] "CRAN"
    
    $Packages$jquerylib$Hash
    [1] "5aab57a3bd297eee1c1d862735972182"
    
    
    $Packages$isoband
    $Packages$isoband$Package
    [1] "isoband"
    
    $Packages$isoband$Version
    [1] "0.2.5"
    
    $Packages$isoband$Source
    [1] "Repository"
    
    $Packages$isoband$Imports
    [1] "grid"  "utils"
    
    $Packages$isoband$Repository
    [1] "CRAN"
    
    $Packages$isoband$Hash
    [1] "7ab57a6de7f48a8dc84910d1eca42883"
    
    
    $Packages$glue
    $Packages$glue$Package
    [1] "glue"
    
    $Packages$glue$Version
    [1] "1.6.2"
    
    $Packages$glue$Source
    [1] "Repository"
    
    $Packages$glue$Depends
    [1] "R (>= 3.4)"
    
    $Packages$glue$Imports
    [1] "methods"
    
    $Packages$glue$Repository
    [1] "CRAN"
    
    $Packages$glue$Hash
    [1] "4f2596dfb05dac67b9dc558e5c6fba2e"
    
    
    $Packages$waldo
    $Packages$waldo$Package
    [1] "waldo"
    
    $Packages$waldo$Version
    [1] "0.3.1"
    
    $Packages$waldo$Source
    [1] "Repository"
    
    $Packages$waldo$Imports
    [1] "cli"                "diffobj (>= 0.3.4)" "fansi"             
    [4] "glue"               "methods"            "rematch2"          
    [7] "rlang (>= 0.4.10)"  "tibble"            
    
    $Packages$waldo$Repository
    [1] "CRAN"
    
    $Packages$waldo$Hash
    [1] "ad8cfff5694ac5b3c354f8f2044bd976"
    
    
    $Packages$RColorBrewer
    $Packages$RColorBrewer$Package
    [1] "RColorBrewer"
    
    $Packages$RColorBrewer$Version
    [1] "1.1-2"
    
    $Packages$RColorBrewer$Source
    [1] "Repository"
    
    $Packages$RColorBrewer$Depends
    [1] "R (>= 2.0.0)"
    
    $Packages$RColorBrewer$Repository
    [1] "CRAN"
    
    $Packages$RColorBrewer$Hash
    [1] "e031418365a7f7a766181ab5a41a5716"
    
    
    $Packages$hunspell
    $Packages$hunspell$Package
    [1] "hunspell"
    
    $Packages$hunspell$Version
    [1] "3.0.1"
    
    $Packages$hunspell$Source
    [1] "Repository"
    
    $Packages$hunspell$Depends
    [1] "R (>= 3.0.2)"
    
    $Packages$hunspell$Imports
    [1] "Rcpp"   "digest"
    
    $Packages$hunspell$LinkingTo
    [1] "Rcpp (>= 0.12.12)"
    
    $Packages$hunspell$Repository
    [1] "CRAN"
    
    $Packages$hunspell$Hash
    [1] "3987784c19192ad0f2261c456d936df1"
    
    
    $Packages$lhs
    $Packages$lhs$Package
    [1] "lhs"
    
    $Packages$lhs$Version
    [1] "1.1.4"
    
    $Packages$lhs$Source
    [1] "Repository"
    
    $Packages$lhs$Depends
    [1] "R (>= 3.4.0)"
    
    $Packages$lhs$Imports
    [1] "Rcpp"
    
    $Packages$lhs$LinkingTo
    [1] "Rcpp"
    
    $Packages$lhs$Repository
    [1] "CRAN"
    
    $Packages$lhs$Hash
    [1] "452453e7fc05361e7c8cf98e3972d22c"
    
    
    $Packages$slider
    $Packages$slider$Package
    [1] "slider"
    
    $Packages$slider$Version
    [1] "0.2.2"
    
    $Packages$slider$Source
    [1] "Repository"
    
    $Packages$slider$Depends
    [1] "R (>= 3.3.0)"
    
    $Packages$slider$Imports
    [1] "ellipsis (>= 0.3.1)" "glue"                "rlang (>= 0.4.5)"   
    [4] "vctrs (>= 0.3.6)"    "warp"               
    
    $Packages$slider$LinkingTo
    [1] "vctrs (>= 0.3.6)"
    
    $Packages$slider$Repository
    [1] "CRAN"
    
    $Packages$slider$Hash
    [1] "5237bd176dc0c4dd7eb8dcdafe514de3"
    
    
    $Packages$modelr
    $Packages$modelr$Package
    [1] "modelr"
    
    $Packages$modelr$Version
    [1] "0.1.8"
    
    $Packages$modelr$Source
    [1] "Repository"
    
    $Packages$modelr$Depends
    [1] "R (>= 3.2)"
    
    $Packages$modelr$Imports
    [1] "broom"            "magrittr"         "purrr (>= 0.2.2)" "rlang (>= 0.2.0)"
    [5] "tibble"           "tidyr (>= 0.8.0)" "tidyselect"       "vctrs"           
    
    $Packages$modelr$Repository
    [1] "CRAN"
    
    $Packages$modelr$Hash
    [1] "9fd59716311ee82cba83dc2826fc5577"
    
    
    $Packages$lambda.r
    $Packages$lambda.r$Package
    [1] "lambda.r"
    
    $Packages$lambda.r$Version
    [1] "1.2.4"
    
    $Packages$lambda.r$Source
    [1] "Repository"
    
    $Packages$lambda.r$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$lambda.r$Imports
    [1] "formatR"
    
    $Packages$lambda.r$Repository
    [1] "CRAN"
    
    $Packages$lambda.r$Hash
    [1] "b1e925c4b9ffeb901bacf812cbe9a6ad"
    
    
    $Packages$matrixStats
    $Packages$matrixStats$Package
    [1] "matrixStats"
    
    $Packages$matrixStats$Version
    [1] "0.61.0"
    
    $Packages$matrixStats$Source
    [1] "Repository"
    
    $Packages$matrixStats$Depends
    [1] "R (>= 2.12.0)"
    
    $Packages$matrixStats$Repository
    [1] "CRAN"
    
    $Packages$matrixStats$Hash
    [1] "b8e6221fc11247b12ab1b055a6f66c27"
    
    
    $Packages$MatrixGenerics
    $Packages$MatrixGenerics$Package
    [1] "MatrixGenerics"
    
    $Packages$MatrixGenerics$Version
    [1] "1.6.0"
    
    $Packages$MatrixGenerics$Source
    [1] "Bioconductor"
    
    $Packages$MatrixGenerics$Depends
    [1] "matrixStats (>= 0.60.1)"
    
    $Packages$MatrixGenerics$Imports
    [1] "methods"
    
    $Packages$MatrixGenerics$git_url
    [1] "https://git.bioconductor.org/packages/MatrixGenerics"
    
    $Packages$MatrixGenerics$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$MatrixGenerics$git_last_commit
    [1] "4588a60"
    
    $Packages$MatrixGenerics$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$MatrixGenerics$Hash
    [1] "506b92cb263d9e014a391b41939badcf"
    
    
    $Packages$stringr
    $Packages$stringr$Package
    [1] "stringr"
    
    $Packages$stringr$Version
    [1] "1.4.0"
    
    $Packages$stringr$Source
    [1] "Repository"
    
    $Packages$stringr$Depends
    [1] "R (>= 3.1)"
    
    $Packages$stringr$Imports
    [1] "glue (>= 1.2.0)"    "magrittr"           "stringi (>= 1.1.7)"
    
    $Packages$stringr$Repository
    [1] "CRAN"
    
    $Packages$stringr$Hash
    [1] "0759e6b6c0957edb1311028a49a35e76"
    
    
    $Packages$lava
    $Packages$lava$Package
    [1] "lava"
    
    $Packages$lava$Version
    [1] "1.6.10"
    
    $Packages$lava$Source
    [1] "Repository"
    
    $Packages$lava$Depends
    [1] "R (>= 3.0)"
    
    $Packages$lava$Imports
     [1] "future.apply" "progressr"    "grDevices"    "graphics"     "methods"     
     [6] "numDeriv"     "stats"        "survival"     "SQUAREM"      "utils"       
    
    $Packages$lava$Repository
    [1] "CRAN"
    
    $Packages$lava$Hash
    [1] "4c31a28528978d8689145f5274ce9058"
    
    
    $Packages$workflows
    $Packages$workflows$Package
    [1] "workflows"
    
    $Packages$workflows$Version
    [1] "0.2.4"
    
    $Packages$workflows$Source
    [1] "Repository"
    
    $Packages$workflows$Depends
    [1] "R (>= 3.2)"
    
    $Packages$workflows$Imports
     [1] "cli (>= 2.0.0)"        "ellipsis (>= 0.2.0)"   "generics (>= 0.1.0)"  
     [4] "glue"                  "hardhat (>= 0.1.6)"    "lifecycle (>= 1.0.0)" 
     [7] "parsnip (>= 0.1.5)"    "rlang (>= 0.4.1)"      "tidyselect (>= 1.1.0)"
    [10] "vctrs (>= 0.3.6)"     
    
    $Packages$workflows$Repository
    [1] "CRAN"
    
    $Packages$workflows$Hash
    [1] "a087bd9bcc22f572730a8949a2ab4a47"
    
    
    $Packages$recipes
    $Packages$recipes$Package
    [1] "recipes"
    
    $Packages$recipes$Version
    [1] "0.2.0"
    
    $Packages$recipes$Source
    [1] "Repository"
    
    $Packages$recipes$Depends
    [1] "dplyr"      "R (>= 3.1)"
    
    $Packages$recipes$Imports
     [1] "ellipsis"                 "generics (>= 0.1.0.9000)"
     [3] "glue"                     "gower"                   
     [5] "hardhat (>= 0.1.6.9001)"  "ipred (>= 0.9-12)"       
     [7] "lifecycle"                "lubridate"               
     [9] "magrittr"                 "Matrix"                  
    [11] "purrr (>= 0.2.3)"         "rlang (>= 0.4.0)"        
    [13] "stats"                    "tibble"                  
    [15] "tidyr (>= 1.0.0)"         "tidyselect (>= 1.1.0)"   
    [17] "timeDate"                 "utils"                   
    [19] "vctrs"                    "withr"                   
    
    $Packages$recipes$Repository
    [1] "CRAN"
    
    $Packages$recipes$Hash
    [1] "2a0f64b148f3064f980404287b511f98"
    
    
    $Packages$kableExtra
    $Packages$kableExtra$Package
    [1] "kableExtra"
    
    $Packages$kableExtra$Version
    [1] "1.3.4"
    
    $Packages$kableExtra$Source
    [1] "Repository"
    
    $Packages$kableExtra$Depends
    [1] "R (>= 3.1.0)"
    
    $Packages$kableExtra$Imports
     [1] "knitr (>= 1.16)"      "magrittr"             "stringr (>= 1.0)"    
     [4] "xml2 (>= 1.1.1)"      "rvest"                "rmarkdown (>= 1.6.0)"
     [7] "scales"               "viridisLite"          "stats"               
    [10] "grDevices"            "htmltools"            "rstudioapi"          
    [13] "glue"                 "tools"                "webshot"             
    [16] "digest"               "graphics"             "svglite"             
    
    $Packages$kableExtra$Repository
    [1] "CRAN"
    
    $Packages$kableExtra$Hash
    [1] "49b625e6aabe4c5f091f5850aba8ff78"
    
    
    $Packages$phyloseq
    $Packages$phyloseq$Package
    [1] "phyloseq"
    
    $Packages$phyloseq$Version
    [1] "1.38.0"
    
    $Packages$phyloseq$Source
    [1] "Bioconductor"
    
    $Packages$phyloseq$Depends
    [1] "R (>= 3.3.0)"
    
    $Packages$phyloseq$Imports
     [1] "ade4 (>= 1.7.4)"          "ape (>= 5.0)"            
     [3] "Biobase (>= 2.36.2)"      "BiocGenerics (>= 0.22.0)"
     [5] "biomformat (>= 1.0.0)"    "Biostrings (>= 2.40.0)"  
     [7] "cluster (>= 2.0.4)"       "data.table (>= 1.10.4)"  
     [9] "foreach (>= 1.4.3)"       "ggplot2 (>= 2.1.0)"      
    [11] "igraph (>= 1.0.1)"        "methods (>= 3.3.0)"      
    [13] "multtest (>= 2.28.0)"     "plyr (>= 1.8.3)"         
    [15] "reshape2 (>= 1.4.1)"      "scales (>= 0.4.0)"       
    [17] "vegan (>= 2.5)"          
    
    $Packages$phyloseq$git_url
    [1] "https://git.bioconductor.org/packages/phyloseq"
    
    $Packages$phyloseq$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$phyloseq$git_last_commit
    [1] "1e2409a"
    
    $Packages$phyloseq$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$phyloseq$Hash
    [1] "cf198d6103ff6270c347f475ce325c14"
    
    
    $Packages$labeling
    $Packages$labeling$Package
    [1] "labeling"
    
    $Packages$labeling$Version
    [1] "0.4.2"
    
    $Packages$labeling$Source
    [1] "Repository"
    
    $Packages$labeling$Imports
    [1] "stats"    "graphics"
    
    $Packages$labeling$Repository
    [1] "CRAN"
    
    $Packages$labeling$Hash
    [1] "3d5108641f47470611a32d0bdf357a72"
    
    
    $Packages$httpuv
    $Packages$httpuv$Package
    [1] "httpuv"
    
    $Packages$httpuv$Version
    [1] "1.6.5"
    
    $Packages$httpuv$Source
    [1] "Repository"
    
    $Packages$httpuv$Depends
    [1] "R (>= 2.15.1)"
    
    $Packages$httpuv$Imports
    [1] "Rcpp (>= 1.0.7)"  "utils"            "R6"               "promises"        
    [5] "later (>= 0.8.0)"
    
    $Packages$httpuv$LinkingTo
    [1] "Rcpp"  "later"
    
    $Packages$httpuv$Repository
    [1] "CRAN"
    
    $Packages$httpuv$Hash
    [1] "97fe71f0a4a1c9890e6c2128afa04bc0"
    
    
    $Packages$biomformat
    $Packages$biomformat$Package
    [1] "biomformat"
    
    $Packages$biomformat$Version
    [1] "1.22.0"
    
    $Packages$biomformat$Source
    [1] "Bioconductor"
    
    $Packages$biomformat$Depends
    [1] "R (>= 3.2)" "methods"   
    
    $Packages$biomformat$Imports
    [1] "plyr (>= 1.8)"        "jsonlite (>= 0.9.16)" "Matrix (>= 1.2)"     
    [4] "rhdf5"               
    
    $Packages$biomformat$git_url
    [1] "https://git.bioconductor.org/packages/biomformat"
    
    $Packages$biomformat$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$biomformat$git_last_commit
    [1] "ab7c641"
    
    $Packages$biomformat$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$biomformat$Hash
    [1] "fd844de531673fbad8baac4a67b9ce48"
    
    
    $Packages$class
    $Packages$class$Package
    [1] "class"
    
    $Packages$class$Version
    [1] "7.3-20"
    
    $Packages$class$Source
    [1] "Repository"
    
    $Packages$class$Depends
    [1] "R (>= 3.0.0)" "stats"        "utils"       
    
    $Packages$class$Imports
    [1] "MASS"
    
    $Packages$class$Repository
    [1] "CRAN"
    
    $Packages$class$Hash
    [1] "da09d82223e669d270e47ed24ac8686e"
    
    
    $Packages$speedyseq
    $Packages$speedyseq$Package
    [1] "speedyseq"
    
    $Packages$speedyseq$Version
    [1] "0.5.3.9018"
    
    $Packages$speedyseq$Source
    [1] "GitHub"
    
    $Packages$speedyseq$Depends
    [1] "phyloseq"
    
    $Packages$speedyseq$Imports
     [1] "ape"              "Biostrings"       "castor"           "data.table"      
     [5] "dplyr"            "ggplot2"          "magrittr"         "methods"         
     [9] "purrr"            "rlang"            "tibble"           "tidyr"           
    [13] "scales"           "stringr"          "vctrs (>= 0.3.0)" "vegan"           
    
    $Packages$speedyseq$RemoteType
    [1] "github"
    
    $Packages$speedyseq$RemoteHost
    [1] "api.github.com"
    
    $Packages$speedyseq$RemoteUsername
    [1] "mikemc"
    
    $Packages$speedyseq$RemoteRepo
    [1] "speedyseq"
    
    $Packages$speedyseq$RemoteRef
    [1] "main"
    
    $Packages$speedyseq$RemoteSha
    [1] "ceb941fdd482fe4bf9610f80970050e24f369be9"
    
    $Packages$speedyseq$Hash
    [1] "040ffc1e11abd6bc5c9f0be6bc5ab11a"
    
    
    $Packages$hoardr
    $Packages$hoardr$Package
    [1] "hoardr"
    
    $Packages$hoardr$Version
    [1] "0.5.2"
    
    $Packages$hoardr$Source
    [1] "Repository"
    
    $Packages$hoardr$Imports
    [1] "R6 (>= 2.2.0)"       "rappdirs (>= 0.3.1)" "digest"             
    
    $Packages$hoardr$Repository
    [1] "CRAN"
    
    $Packages$hoardr$Hash
    [1] "ca8a25aa079a8cb4fe44fe12b17b0eda"
    
    
    $Packages$tokenizers
    $Packages$tokenizers$Package
    [1] "tokenizers"
    
    $Packages$tokenizers$Version
    [1] "0.2.1"
    
    $Packages$tokenizers$Source
    [1] "Repository"
    
    $Packages$tokenizers$Depends
    [1] "R (>= 3.1.3)"
    
    $Packages$tokenizers$Imports
    [1] "stringi (>= 1.0.1)"   "Rcpp (>= 0.12.3)"     "SnowballC (>= 0.5.1)"
    
    $Packages$tokenizers$LinkingTo
    [1] "Rcpp"
    
    $Packages$tokenizers$Repository
    [1] "CRAN"
    
    $Packages$tokenizers$Hash
    [1] "a064f646b3a692e62dfb5d9ea690a4ea"
    
    
    $Packages$BiocNeighbors
    $Packages$BiocNeighbors$Package
    [1] "BiocNeighbors"
    
    $Packages$BiocNeighbors$Version
    [1] "1.12.0"
    
    $Packages$BiocNeighbors$Source
    [1] "Bioconductor"
    
    $Packages$BiocNeighbors$Imports
    [1] "Rcpp"         "S4Vectors"    "BiocParallel" "stats"        "methods"     
    [6] "Matrix"      
    
    $Packages$BiocNeighbors$LinkingTo
    [1] "Rcpp"     "RcppHNSW"
    
    $Packages$BiocNeighbors$git_url
    [1] "https://git.bioconductor.org/packages/BiocNeighbors"
    
    $Packages$BiocNeighbors$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocNeighbors$git_last_commit
    [1] "3c8a290"
    
    $Packages$BiocNeighbors$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$BiocNeighbors$Hash
    [1] "5fc6686b8fcbe45e684beecb15e3f3a5"
    
    
    $Packages$contentid
    $Packages$contentid$Package
    [1] "contentid"
    
    $Packages$contentid$Version
    [1] "0.0.15"
    
    $Packages$contentid$Source
    [1] "Repository"
    
    $Packages$contentid$Depends
    [1] "R (>= 4.0)"
    
    $Packages$contentid$Imports
    [1] "openssl (>= 1.4.2)" "httr"               "curl"              
    [4] "fs"                 "tools"              "methods"           
    
    $Packages$contentid$Repository
    [1] "CRAN"
    
    $Packages$contentid$Hash
    [1] "1a6637b36573e009652ff5f98280ac97"
    
    
    $Packages$clipr
    $Packages$clipr$Package
    [1] "clipr"
    
    $Packages$clipr$Version
    [1] "0.8.0"
    
    $Packages$clipr$Source
    [1] "Repository"
    
    $Packages$clipr$Imports
    [1] "utils"
    
    $Packages$clipr$Repository
    [1] "CRAN"
    
    $Packages$clipr$Hash
    [1] "3f038e5ac7f41d4ac41ce658c85e3042"
    
    
    $Packages$cpp11
    $Packages$cpp11$Package
    [1] "cpp11"
    
    $Packages$cpp11$Version
    [1] "0.4.2"
    
    $Packages$cpp11$Source
    [1] "Repository"
    
    $Packages$cpp11$Repository
    [1] "CRAN"
    
    $Packages$cpp11$Hash
    [1] "fa53ce256cd280f468c080a58ea5ba8c"
    
    
    $Packages$annotate
    $Packages$annotate$Package
    [1] "annotate"
    
    $Packages$annotate$Version
    [1] "1.72.0"
    
    $Packages$annotate$Source
    [1] "Bioconductor"
    
    $Packages$annotate$Depends
    [1] "R (>= 2.10)"               "AnnotationDbi (>= 1.27.5)"
    [3] "XML"                      
    
    $Packages$annotate$Imports
    [1] "Biobase"                  "DBI"                     
    [3] "xtable"                   "graphics"                
    [5] "utils"                    "stats"                   
    [7] "methods"                  "BiocGenerics (>= 0.13.8)"
    [9] "httr"                    
    
    $Packages$annotate$git_url
    [1] "https://git.bioconductor.org/packages/annotate"
    
    $Packages$annotate$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$annotate$git_last_commit
    [1] "67ac76a"
    
    $Packages$annotate$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$annotate$Hash
    [1] "885c32b60eab51236e4afe8921049ff1"
    
    
    $Packages$webshot
    $Packages$webshot$Package
    [1] "webshot"
    
    $Packages$webshot$Version
    [1] "0.5.2"
    
    $Packages$webshot$Source
    [1] "Repository"
    
    $Packages$webshot$Depends
    [1] "R (>= 3.0)"
    
    $Packages$webshot$Imports
    [1] "magrittr" "jsonlite" "callr"   
    
    $Packages$webshot$Repository
    [1] "CRAN"
    
    $Packages$webshot$Hash
    [1] "e99d80ad34457a4853674e89d5e806de"
    
    
    $Packages$jsonlite
    $Packages$jsonlite$Package
    [1] "jsonlite"
    
    $Packages$jsonlite$Version
    [1] "1.8.0"
    
    $Packages$jsonlite$Source
    [1] "Repository"
    
    $Packages$jsonlite$Depends
    [1] "methods"
    
    $Packages$jsonlite$Repository
    [1] "CRAN"
    
    $Packages$jsonlite$Hash
    [1] "d07e729b27b372429d42d24d503613a0"
    
    
    $Packages$XVector
    $Packages$XVector$Package
    [1] "XVector"
    
    $Packages$XVector$Version
    [1] "0.34.0"
    
    $Packages$XVector$Source
    [1] "Bioconductor"
    
    $Packages$XVector$Depends
    [1] "R (>= 4.0.0)"             "methods"                 
    [3] "BiocGenerics (>= 0.37.0)" "S4Vectors (>= 0.27.12)"  
    [5] "IRanges (>= 2.23.9)"     
    
    $Packages$XVector$Imports
    [1] "methods"      "utils"        "tools"        "zlibbioc"     "BiocGenerics"
    [6] "S4Vectors"    "IRanges"     
    
    $Packages$XVector$LinkingTo
    [1] "S4Vectors" "IRanges"  
    
    $Packages$XVector$git_url
    [1] "https://git.bioconductor.org/packages/XVector"
    
    $Packages$XVector$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$XVector$git_last_commit
    [1] "06adb25"
    
    $Packages$XVector$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$XVector$Hash
    [1] "9830cb6fc09640591c2f6436467af3c5"
    
    
    $Packages$bit
    $Packages$bit$Package
    [1] "bit"
    
    $Packages$bit$Version
    [1] "4.0.4"
    
    $Packages$bit$Source
    [1] "Repository"
    
    $Packages$bit$Depends
    [1] "R (>= 2.9.2)"
    
    $Packages$bit$Repository
    [1] "CRAN"
    
    $Packages$bit$Hash
    [1] "f36715f14d94678eea9933af927bc15d"
    
    
    $Packages$mime
    $Packages$mime$Package
    [1] "mime"
    
    $Packages$mime$Version
    [1] "0.12"
    
    $Packages$mime$Source
    [1] "Repository"
    
    $Packages$mime$Imports
    [1] "tools"
    
    $Packages$mime$Repository
    [1] "CRAN"
    
    $Packages$mime$Hash
    [1] "18e9c28c1d3ca1560ce30658b22ce104"
    
    
    $Packages$systemfonts
    $Packages$systemfonts$Package
    [1] "systemfonts"
    
    $Packages$systemfonts$Version
    [1] "1.0.4"
    
    $Packages$systemfonts$Source
    [1] "Repository"
    
    $Packages$systemfonts$Depends
    [1] "R (>= 3.2.0)"
    
    $Packages$systemfonts$LinkingTo
    [1] "cpp11 (>= 0.2.1)"
    
    $Packages$systemfonts$Repository
    [1] "CRAN"
    
    $Packages$systemfonts$Hash
    [1] "90b28393209827327de889f49935140a"
    
    
    $Packages$gridExtra
    $Packages$gridExtra$Package
    [1] "gridExtra"
    
    $Packages$gridExtra$Version
    [1] "2.3"
    
    $Packages$gridExtra$Source
    [1] "Repository"
    
    $Packages$gridExtra$Imports
    [1] "gtable"    "grid"      "grDevices" "graphics"  "utils"    
    
    $Packages$gridExtra$Repository
    [1] "CRAN"
    
    $Packages$gridExtra$Hash
    [1] "7d7f283939f563670a697165b2cf5560"
    
    
    $Packages$ids
    $Packages$ids$Package
    [1] "ids"
    
    $Packages$ids$Version
    [1] "1.0.1"
    
    $Packages$ids$Source
    [1] "Repository"
    
    $Packages$ids$Imports
    [1] "openssl" "uuid"   
    
    $Packages$ids$Repository
    [1] "CRAN"
    
    $Packages$ids$Hash
    [1] "99df65cfef20e525ed38c3d2577f7190"
    
    
    $Packages$Rsamtools
    $Packages$Rsamtools$Package
    [1] "Rsamtools"
    
    $Packages$Rsamtools$Version
    [1] "2.10.0"
    
    $Packages$Rsamtools$Source
    [1] "Bioconductor"
    
    $Packages$Rsamtools$Depends
    [1] "methods"                   "GenomeInfoDb (>= 1.1.3)"  
    [3] "GenomicRanges (>= 1.31.8)" "Biostrings (>= 2.47.6)"   
    [5] "R (>= 3.5.0)"             
    
    $Packages$Rsamtools$Imports
    [1] "utils"                    "BiocGenerics (>= 0.25.1)"
    [3] "S4Vectors (>= 0.17.25)"   "IRanges (>= 2.13.12)"    
    [5] "XVector (>= 0.19.7)"      "zlibbioc"                
    [7] "bitops"                   "BiocParallel"            
    [9] "stats"                   
    
    $Packages$Rsamtools$LinkingTo
    [1] "Rhtslib (>= 1.17.7)" "S4Vectors"           "IRanges"            
    [4] "XVector"             "Biostrings"         
    
    $Packages$Rsamtools$git_url
    [1] "https://git.bioconductor.org/packages/Rsamtools"
    
    $Packages$Rsamtools$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$Rsamtools$git_last_commit
    [1] "b19738e"
    
    $Packages$Rsamtools$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$Rsamtools$Hash
    [1] "55e8c5ccea90dfc370e4abdbc5a01cb9"
    
    
    $Packages$stringi
    $Packages$stringi$Package
    [1] "stringi"
    
    $Packages$stringi$Version
    [1] "1.7.6"
    
    $Packages$stringi$Source
    [1] "Repository"
    
    $Packages$stringi$Depends
    [1] "R (>= 3.1)"
    
    $Packages$stringi$Imports
    [1] "tools" "utils" "stats"
    
    $Packages$stringi$Repository
    [1] "CRAN"
    
    $Packages$stringi$Hash
    [1] "bba431031d30789535745a9627ac9271"
    
    
    $Packages$taxizedb
    $Packages$taxizedb$Package
    [1] "taxizedb"
    
    $Packages$taxizedb$Version
    [1] "0.3.0"
    
    $Packages$taxizedb$Source
    [1] "Repository"
    
    $Packages$taxizedb$Imports
     [1] "curl (>= 2.4)"      "DBI (>= 0.6-1)"     "RSQLite (>= 1.1.2)"
     [4] "dplyr (>= 0.7.0)"   "tibble"             "rlang"             
     [7] "readr (>= 1.1.1)"   "dbplyr (>= 1.0.0)"  "magrittr (>= 1.5)" 
    [10] "hoardr (>= 0.1.0)" 
    
    $Packages$taxizedb$Repository
    [1] "CRAN"
    
    $Packages$taxizedb$Hash
    [1] "9303b0cde2b629f5a831e171f6f14b41"
    
    
    $Packages$processx
    $Packages$processx$Package
    [1] "processx"
    
    $Packages$processx$Version
    [1] "3.5.2"
    
    $Packages$processx$Source
    [1] "Repository"
    
    $Packages$processx$Imports
    [1] "ps (>= 1.2.0)" "R6"            "utils"        
    
    $Packages$processx$Repository
    [1] "CRAN"
    
    $Packages$processx$Hash
    [1] "0cbca2bc4d16525d009c4dbba156b37c"
    
    
    $Packages$rbibutils
    $Packages$rbibutils$Package
    [1] "rbibutils"
    
    $Packages$rbibutils$Version
    [1] "2.2.7"
    
    $Packages$rbibutils$Source
    [1] "Repository"
    
    $Packages$rbibutils$Depends
    [1] "R (>= 2.10)"
    
    $Packages$rbibutils$Imports
    [1] "utils" "tools"
    
    $Packages$rbibutils$Repository
    [1] "CRAN"
    
    $Packages$rbibutils$Hash
    [1] "69f4a840b0c2aa99f5ce473fd7fa21bd"
    
    
    $Packages$Rdpack
    $Packages$Rdpack$Package
    [1] "Rdpack"
    
    $Packages$Rdpack$Version
    [1] "2.1.4"
    
    $Packages$Rdpack$Source
    [1] "Repository"
    
    $Packages$Rdpack$Depends
    [1] "R (>= 2.15.0)" "methods"      
    
    $Packages$Rdpack$Imports
    [1] "tools"              "utils"              "rbibutils (>= 1.3)"
    
    $Packages$Rdpack$Repository
    [1] "CRAN"
    
    $Packages$Rdpack$Hash
    [1] "7befe052afcbb032eacbbbc0babca848"
    
    
    $Packages$yulab.utils
    $Packages$yulab.utils$Package
    [1] "yulab.utils"
    
    $Packages$yulab.utils$Version
    [1] "0.0.4"
    
    $Packages$yulab.utils$Source
    [1] "Repository"
    
    $Packages$yulab.utils$Imports
    [1] "utils"
    
    $Packages$yulab.utils$Repository
    [1] "CRAN"
    
    $Packages$yulab.utils$Hash
    [1] "922e11dcf40bb5dfcf3fe5e714d0dc35"
    
    
    $Packages$hardhat
    $Packages$hardhat$Package
    [1] "hardhat"
    
    $Packages$hardhat$Version
    [1] "0.2.0"
    
    $Packages$hardhat$Source
    [1] "Repository"
    
    $Packages$hardhat$Depends
    [1] "R (>= 2.10)"
    
    $Packages$hardhat$Imports
    [1] "glue"             "rlang (>= 0.4.2)" "tibble"           "vctrs (>= 0.3.0)"
    
    $Packages$hardhat$Repository
    [1] "CRAN"
    
    $Packages$hardhat$Hash
    [1] "726bc6915fd4533da5d49f4dfe6df1a6"
    
    
    $Packages$bitops
    $Packages$bitops$Package
    [1] "bitops"
    
    $Packages$bitops$Version
    [1] "1.0-7"
    
    $Packages$bitops$Source
    [1] "Repository"
    
    $Packages$bitops$Repository
    [1] "CRAN"
    
    $Packages$bitops$Hash
    [1] "b7d8d8ee39869c18d8846a184dd8a1af"
    
    
    $Packages$cli
    $Packages$cli$Package
    [1] "cli"
    
    $Packages$cli$Version
    [1] "3.2.0"
    
    $Packages$cli$Source
    [1] "Repository"
    
    $Packages$cli$Depends
    [1] "R (>= 2.10)"
    
    $Packages$cli$Imports
    [1] "glue (>= 1.6.0)" "utils"          
    
    $Packages$cli$Repository
    [1] "CRAN"
    
    $Packages$cli$Hash
    [1] "1bdb126893e9ce6aae50ad1d6fc32faf"
    
    
    $Packages$rhdf5filters
    $Packages$rhdf5filters$Package
    [1] "rhdf5filters"
    
    $Packages$rhdf5filters$Version
    [1] "1.6.0"
    
    $Packages$rhdf5filters$Source
    [1] "Bioconductor"
    
    $Packages$rhdf5filters$LinkingTo
    [1] "Rhdf5lib"
    
    $Packages$rhdf5filters$git_url
    [1] "https://git.bioconductor.org/packages/rhdf5filters"
    
    $Packages$rhdf5filters$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$rhdf5filters$git_last_commit
    [1] "5f7f3a5"
    
    $Packages$rhdf5filters$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$rhdf5filters$Hash
    [1] "429e1c984126e2a2594ca16632008027"
    
    
    $Packages$ROSE
    $Packages$ROSE$Package
    [1] "ROSE"
    
    $Packages$ROSE$Version
    [1] "0.0-4"
    
    $Packages$ROSE$Source
    [1] "Repository"
    
    $Packages$ROSE$Repository
    [1] "CRAN"
    
    $Packages$ROSE$Hash
    [1] "c6e79999925f0db3f2d88b187d014bf8"
    
    
    $Packages$RSQLite
    $Packages$RSQLite$Package
    [1] "RSQLite"
    
    $Packages$RSQLite$Version
    [1] "2.2.10"
    
    $Packages$RSQLite$Source
    [1] "Repository"
    
    $Packages$RSQLite$Depends
    [1] "R (>= 3.1.0)"
    
    $Packages$RSQLite$Imports
    [1] "bit64"           "blob (>= 1.2.0)" "DBI (>= 1.1.0)"  "memoise"        
    [5] "methods"         "pkgconfig"       "Rcpp (>= 1.0.7)"
    
    $Packages$RSQLite$LinkingTo
    [1] "plogr (>= 0.2.0)" "Rcpp"            
    
    $Packages$RSQLite$Repository
    [1] "CRAN"
    
    $Packages$RSQLite$Hash
    [1] "5d747a4aa7492f22d659be8c343e9045"
    
    
    $Packages$randomForest
    $Packages$randomForest$Package
    [1] "randomForest"
    
    $Packages$randomForest$Version
    [1] "4.7-1"
    
    $Packages$randomForest$Source
    [1] "Repository"
    
    $Packages$randomForest$Depends
    [1] "R (>= 4.1.0)" "stats"       
    
    $Packages$randomForest$Repository
    [1] "CRAN"
    
    $Packages$randomForest$Hash
    [1] "dd4f861c0353a869964024886c93b7d5"
    
    
    $Packages$mlr
    $Packages$mlr$Package
    [1] "mlr"
    
    $Packages$mlr$Version
    [1] "2.19.0"
    
    $Packages$mlr$Source
    [1] "Repository"
    
    $Packages$mlr$Depends
    [1] "ParamHelpers (>= 1.10)" "R (>= 3.0.2)"          
    
    $Packages$mlr$Imports
     [1] "backports (>= 1.1.0)"   "BBmisc (>= 1.11)"       "checkmate (>= 1.8.2)"  
     [4] "data.table (>= 1.12.4)" "ggplot2"                "methods"               
     [7] "parallelMap (>= 1.3)"   "stats"                  "stringi"               
    [10] "survival"               "utils"                  "XML"                   
    
    $Packages$mlr$Repository
    [1] "CRAN"
    
    $Packages$mlr$Hash
    [1] "cf4429d0d8bb3c06fe8b07eb85ffead9"
    
    
    $Packages$tidyr
    $Packages$tidyr$Package
    [1] "tidyr"
    
    $Packages$tidyr$Version
    [1] "1.2.0"
    
    $Packages$tidyr$Source
    [1] "Repository"
    
    $Packages$tidyr$Depends
    [1] "R (>= 3.1)"
    
    $Packages$tidyr$Imports
     [1] "dplyr (>= 1.0.0)"      "ellipsis (>= 0.1.0)"   "glue"                 
     [4] "lifecycle"             "magrittr"              "purrr"                
     [7] "rlang"                 "tibble (>= 2.1.1)"     "tidyselect (>= 1.1.0)"
    [10] "utils"                 "vctrs (>= 0.3.7)"     
    
    $Packages$tidyr$LinkingTo
    [1] "cpp11 (>= 0.4.0)"
    
    $Packages$tidyr$Repository
    [1] "CRAN"
    
    $Packages$tidyr$Hash
    [1] "d8b95b7fee945d7da6888cf7eb71a49c"
    
    
    $Packages$data.table
    $Packages$data.table$Package
    [1] "data.table"
    
    $Packages$data.table$Version
    [1] "1.14.2"
    
    $Packages$data.table$Source
    [1] "Repository"
    
    $Packages$data.table$Depends
    [1] "R (>= 3.1.0)"
    
    $Packages$data.table$Imports
    [1] "methods"
    
    $Packages$data.table$Repository
    [1] "CRAN"
    
    $Packages$data.table$Hash
    [1] "36b67b5adf57b292923f5659f5f0c853"
    
    
    $Packages$ParamHelpers
    $Packages$ParamHelpers$Package
    [1] "ParamHelpers"
    
    $Packages$ParamHelpers$Version
    [1] "1.14"
    
    $Packages$ParamHelpers$Source
    [1] "Repository"
    
    $Packages$ParamHelpers$Imports
    [1] "backports"            "BBmisc (>= 1.10)"     "checkmate (>= 1.8.2)"
    [4] "fastmatch"            "methods"             
    
    $Packages$ParamHelpers$Repository
    [1] "CRAN"
    
    $Packages$ParamHelpers$Hash
    [1] "0feaa54a0bcf59d72eb758da748e0904"
    
    
    $Packages$rstudioapi
    $Packages$rstudioapi$Package
    [1] "rstudioapi"
    
    $Packages$rstudioapi$Version
    [1] "0.13"
    
    $Packages$rstudioapi$Source
    [1] "Repository"
    
    $Packages$rstudioapi$Repository
    [1] "CRAN"
    
    $Packages$rstudioapi$Hash
    [1] "06c85365a03fdaf699966cc1d3cf53ea"
    
    
    $Packages$microbiome
    $Packages$microbiome$Package
    [1] "microbiome"
    
    $Packages$microbiome$Version
    [1] "1.16.0"
    
    $Packages$microbiome$Source
    [1] "Bioconductor"
    
    $Packages$microbiome$Depends
    [1] "R (>= 3.6.0)" "phyloseq"     "ggplot2"     
    
    $Packages$microbiome$Imports
    [1] "dplyr"    "reshape2" "Rtsne"    "scales"   "stats"    "tibble"   "tidyr"   
    [8] "utils"    "vegan"   
    
    $Packages$microbiome$git_url
    [1] "https://git.bioconductor.org/packages/microbiome"
    
    $Packages$microbiome$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$microbiome$git_last_commit
    [1] "a7b74b7"
    
    $Packages$microbiome$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$microbiome$Hash
    [1] "4fd8c7b66d1b18224348a02478c81164"
    
    
    $Packages$GenomicAlignments
    $Packages$GenomicAlignments$Package
    [1] "GenomicAlignments"
    
    $Packages$GenomicAlignments$Version
    [1] "1.30.0"
    
    $Packages$GenomicAlignments$Source
    [1] "Bioconductor"
    
    $Packages$GenomicAlignments$Depends
     [1] "R (>= 4.0.0)"                     "methods"                         
     [3] "BiocGenerics (>= 0.37.0)"         "S4Vectors (>= 0.27.12)"          
     [5] "IRanges (>= 2.23.9)"              "GenomeInfoDb (>= 1.13.1)"        
     [7] "GenomicRanges (>= 1.41.5)"        "SummarizedExperiment (>= 1.9.13)"
     [9] "Biostrings (>= 2.55.7)"           "Rsamtools (>= 1.31.2)"           
    
    $Packages$GenomicAlignments$Imports
     [1] "methods"       "utils"         "stats"         "BiocGenerics" 
     [5] "S4Vectors"     "IRanges"       "GenomicRanges" "Biostrings"   
     [9] "Rsamtools"     "BiocParallel" 
    
    $Packages$GenomicAlignments$LinkingTo
    [1] "S4Vectors" "IRanges"  
    
    $Packages$GenomicAlignments$git_url
    [1] "https://git.bioconductor.org/packages/GenomicAlignments"
    
    $Packages$GenomicAlignments$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$GenomicAlignments$git_last_commit
    [1] "9046119"
    
    $Packages$GenomicAlignments$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$GenomicAlignments$Hash
    [1] "e1e27b5fac829942cf1359555dfbed9a"
    
    
    $Packages$nlme
    $Packages$nlme$Package
    [1] "nlme"
    
    $Packages$nlme$Version
    [1] "3.1-155"
    
    $Packages$nlme$Source
    [1] "Repository"
    
    $Packages$nlme$Depends
    [1] "R (>= 3.4.0)"
    
    $Packages$nlme$Imports
    [1] "graphics" "stats"    "utils"    "lattice" 
    
    $Packages$nlme$Repository
    [1] "CRAN"
    
    $Packages$nlme$Hash
    [1] "74ad940dccc9e977189a5afe5fcdb7ba"
    
    
    $Packages$listenv
    $Packages$listenv$Package
    [1] "listenv"
    
    $Packages$listenv$Version
    [1] "0.8.0"
    
    $Packages$listenv$Source
    [1] "Repository"
    
    $Packages$listenv$Depends
    [1] "R (>= 3.1.2)"
    
    $Packages$listenv$Repository
    [1] "CRAN"
    
    $Packages$listenv$Hash
    [1] "0bde42ee282efb18c7c4e63822f5b4f7"
    
    
    $Packages$ggthemes
    $Packages$ggthemes$Package
    [1] "ggthemes"
    
    $Packages$ggthemes$Version
    [1] "4.2.4"
    
    $Packages$ggthemes$Source
    [1] "Repository"
    
    $Packages$ggthemes$Depends
    [1] "R (>= 3.3.0)"
    
    $Packages$ggthemes$Imports
    [1] "ggplot2 (>= 3.0.0)" "graphics"           "grid"              
    [4] "methods"            "purrr"              "scales"            
    [7] "stringr"            "tibble"            
    
    $Packages$ggthemes$Repository
    [1] "CRAN"
    
    $Packages$ggthemes$Hash
    [1] "617fb2c300f68e8b1ba21043fba88fd7"
    
    
    $Packages$parsnip
    $Packages$parsnip$Package
    [1] "parsnip"
    
    $Packages$parsnip$Version
    [1] "0.2.0.9000"
    
    $Packages$parsnip$Source
    [1] "GitHub"
    
    $Packages$parsnip$Depends
    [1] "R (>= 2.10)"
    
    $Packages$parsnip$Imports
     [1] "cli"                      "dplyr (>= 0.8.0.1)"      
     [3] "generics (>= 0.1.0.9000)" "ggplot2"                 
     [5] "globals"                  "glue"                    
     [7] "hardhat (>= 0.1.6.9001)"  "lifecycle"               
     [9] "magrittr"                 "prettyunits"             
    [11] "purrr"                    "rlang (>= 0.3.1)"        
    [13] "stats"                    "tibble (>= 2.1.1)"       
    [15] "tidyr (>= 1.0.0)"         "utils"                   
    [17] "vctrs (>= 0.2.0)"         "withr"                   
    
    $Packages$parsnip$RemoteType
    [1] "github"
    
    $Packages$parsnip$RemoteUsername
    [1] "tidymodels"
    
    $Packages$parsnip$RemoteRepo
    [1] "parsnip"
    
    $Packages$parsnip$RemoteRef
    [1] "HEAD"
    
    $Packages$parsnip$RemoteSha
    [1] "6029887abf7be2fb0b7d9cebe94956627bbc3308"
    
    $Packages$parsnip$RemoteHost
    [1] "api.github.com"
    
    $Packages$parsnip$Hash
    [1] "a39aaae9e0baa971d80375c4bcbe845b"
    
    
    $Packages$DECIPHER
    $Packages$DECIPHER$Package
    [1] "DECIPHER"
    
    $Packages$DECIPHER$Version
    [1] "2.22.0"
    
    $Packages$DECIPHER$Source
    [1] "Bioconductor"
    
    $Packages$DECIPHER$Depends
    [1] "R (>= 3.5.0)"           "Biostrings (>= 2.59.1)" "RSQLite (>= 1.1)"      
    [4] "stats"                  "parallel"              
    
    $Packages$DECIPHER$Imports
    [1] "methods"   "DBI"       "S4Vectors" "IRanges"   "XVector"  
    
    $Packages$DECIPHER$LinkingTo
    [1] "Biostrings" "S4Vectors"  "IRanges"    "XVector"   
    
    $Packages$DECIPHER$git_url
    [1] "https://git.bioconductor.org/packages/DECIPHER"
    
    $Packages$DECIPHER$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$DECIPHER$git_last_commit
    [1] "45da5ca"
    
    $Packages$DECIPHER$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$DECIPHER$Hash
    [1] "ae960336cac449fa17ba2b9c5cb1d092"
    
    
    $Packages$SnowballC
    $Packages$SnowballC$Package
    [1] "SnowballC"
    
    $Packages$SnowballC$Version
    [1] "0.7.0"
    
    $Packages$SnowballC$Source
    [1] "Repository"
    
    $Packages$SnowballC$Repository
    [1] "CRAN"
    
    $Packages$SnowballC$Hash
    [1] "bc26e07c0d747fd287c370fe355e7b85"
    
    
    $Packages$rsample
    $Packages$rsample$Package
    [1] "rsample"
    
    $Packages$rsample$Version
    [1] "0.1.1"
    
    $Packages$rsample$Source
    [1] "Repository"
    
    $Packages$rsample$Depends
    [1] "R (>= 3.2)"
    
    $Packages$rsample$Imports
     [1] "dplyr (>= 1.0.0)"  "ellipsis"          "furrr"            
     [4] "generics"          "lifecycle"         "methods"          
     [7] "purrr"             "rlang (>= 0.4.10)" "slider (>= 0.1.5)"
    [10] "tibble"            "tidyr"             "tidyselect"       
    [13] "vctrs (>= 0.3.0)" 
    
    $Packages$rsample$Repository
    [1] "CRAN"
    
    $Packages$rsample$Hash
    [1] "d1ac9b0eedb5dabcec1309e072b169cb"
    
    
    $Packages$gridGraphics
    $Packages$gridGraphics$Package
    [1] "gridGraphics"
    
    $Packages$gridGraphics$Version
    [1] "0.5-1"
    
    $Packages$gridGraphics$Source
    [1] "Repository"
    
    $Packages$gridGraphics$Depends
    [1] "grid"     "graphics"
    
    $Packages$gridGraphics$Imports
    [1] "grDevices"
    
    $Packages$gridGraphics$Repository
    [1] "CRAN"
    
    $Packages$gridGraphics$Hash
    [1] "5b79228594f02385d4df4979284879ae"
    
    
    $Packages$R.oo
    $Packages$R.oo$Package
    [1] "R.oo"
    
    $Packages$R.oo$Version
    [1] "1.24.0"
    
    $Packages$R.oo$Source
    [1] "Repository"
    
    $Packages$R.oo$Depends
    [1] "R (>= 2.13.0)"          "R.methodsS3 (>= 1.8.0)"
    
    $Packages$R.oo$Imports
    [1] "methods" "utils"  
    
    $Packages$R.oo$Repository
    [1] "CRAN"
    
    $Packages$R.oo$Hash
    [1] "5709328352717e2f0a9c012be8a97554"
    
    
    $Packages$optparse
    $Packages$optparse$Package
    [1] "optparse"
    
    $Packages$optparse$Version
    [1] "1.7.1"
    
    $Packages$optparse$Source
    [1] "Repository"
    
    $Packages$optparse$Depends
    [1] "R (>= 2.9.0)"
    
    $Packages$optparse$Imports
    [1] "methods"            "getopt (>= 1.20.2)"
    
    $Packages$optparse$Repository
    [1] "CRAN"
    
    $Packages$optparse$Hash
    [1] "2490600671344e847c37a7f75ee458c0"
    
    
    $Packages$infer
    $Packages$infer$Package
    [1] "infer"
    
    $Packages$infer$Version
    [1] "1.0.0"
    
    $Packages$infer$Source
    [1] "Repository"
    
    $Packages$infer$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$infer$Imports
     [1] "dplyr (>= 0.7.0)" "ggplot2"          "glue (>= 1.3.0)"  "grDevices"       
     [5] "magrittr"         "methods"          "purrr"            "rlang (>= 0.2.0)"
     [9] "tibble"           "broom"            "tidyr"            "generics"        
    [13] "patchwork"       
    
    $Packages$infer$Repository
    [1] "CRAN"
    
    $Packages$infer$Hash
    [1] "74128662081cedb158fc98d6dae73c0e"
    
    
    $Packages$dbplyr
    $Packages$dbplyr$Package
    [1] "dbplyr"
    
    $Packages$dbplyr$Version
    [1] "2.1.1"
    
    $Packages$dbplyr$Source
    [1] "Repository"
    
    $Packages$dbplyr$Depends
    [1] "R (>= 3.1)"
    
    $Packages$dbplyr$Imports
     [1] "assertthat (>= 0.2.0)" "blob (>= 1.2.0)"       "DBI (>= 1.0.0)"       
     [4] "dplyr (>= 1.0.4)"      "ellipsis"              "glue (>= 1.2.0)"      
     [7] "lifecycle (>= 1.0.0)"  "magrittr"              "methods"              
    [10] "purrr (>= 0.2.5)"      "R6 (>= 2.2.2)"         "rlang (>= 0.2.0)"     
    [13] "tibble (>= 1.4.2)"     "tidyselect (>= 0.2.4)" "utils"                
    [16] "vctrs"                 "withr"                
    
    $Packages$dbplyr$Repository
    [1] "CRAN"
    
    $Packages$dbplyr$Hash
    [1] "1f37fa4ab2f5f7eded42f78b9a887182"
    
    
    $Packages$BiocGenerics
    $Packages$BiocGenerics$Package
    [1] "BiocGenerics"
    
    $Packages$BiocGenerics$Version
    [1] "0.40.0"
    
    $Packages$BiocGenerics$Source
    [1] "Bioconductor"
    
    $Packages$BiocGenerics$Depends
    [1] "R (>= 4.0.0)" "methods"      "utils"        "graphics"     "stats"       
    
    $Packages$BiocGenerics$Imports
    [1] "methods"  "utils"    "graphics" "stats"   
    
    $Packages$BiocGenerics$git_url
    [1] "https://git.bioconductor.org/packages/BiocGenerics"
    
    $Packages$BiocGenerics$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocGenerics$git_last_commit
    [1] "0bc1e0e"
    
    $Packages$BiocGenerics$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$BiocGenerics$Hash
    [1] "ddc1d29bbe66aaef34b0d17620b61a69"
    
    
    $Packages$readxl
    $Packages$readxl$Package
    [1] "readxl"
    
    $Packages$readxl$Version
    [1] "1.3.1"
    
    $Packages$readxl$Source
    [1] "Repository"
    
    $Packages$readxl$Imports
    [1] "cellranger"        "Rcpp (>= 0.12.18)" "tibble (>= 1.3.1)"
    [4] "utils"            
    
    $Packages$readxl$LinkingTo
    [1] "progress" "Rcpp"    
    
    $Packages$readxl$Repository
    [1] "CRAN"
    
    $Packages$readxl$Hash
    [1] "63537c483c2dbec8d9e3183b3735254a"
    
    
    $Packages$lifecycle
    $Packages$lifecycle$Package
    [1] "lifecycle"
    
    $Packages$lifecycle$Version
    [1] "1.0.1"
    
    $Packages$lifecycle$Source
    [1] "Repository"
    
    $Packages$lifecycle$Depends
    [1] "R (>= 3.3)"
    
    $Packages$lifecycle$Imports
    [1] "glue"              "rlang (>= 0.4.10)"
    
    $Packages$lifecycle$Repository
    [1] "CRAN"
    
    $Packages$lifecycle$Hash
    [1] "a6b6d352e3ed897373ab19d8395c98d0"
    
    
    $Packages$timeDate
    $Packages$timeDate$Package
    [1] "timeDate"
    
    $Packages$timeDate$Version
    [1] "3043.102"
    
    $Packages$timeDate$Source
    [1] "Repository"
    
    $Packages$timeDate$Depends
    [1] "R (>= 2.15.1)" "graphics"      "utils"         "stats"        
    [5] "methods"      
    
    $Packages$timeDate$Repository
    [1] "CRAN"
    
    $Packages$timeDate$Hash
    [1] "fde4fc571f5f61978652c229d4713845"
    
    
    $Packages$ExperimentHub
    $Packages$ExperimentHub$Package
    [1] "ExperimentHub"
    
    $Packages$ExperimentHub$Version
    [1] "2.2.1"
    
    $Packages$ExperimentHub$Source
    [1] "Bioconductor"
    
    $Packages$ExperimentHub$Depends
    [1] "methods"                   "BiocGenerics (>= 0.15.10)"
    [3] "AnnotationHub (>= 2.19.3)" "BiocFileCache (>= 1.5.1)" 
    
    $Packages$ExperimentHub$Imports
    [1] "utils"       "S4Vectors"   "BiocManager" "curl"        "rappdirs"   
    
    $Packages$ExperimentHub$git_url
    [1] "https://git.bioconductor.org/packages/ExperimentHub"
    
    $Packages$ExperimentHub$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$ExperimentHub$git_last_commit
    [1] "4e10686"
    
    $Packages$ExperimentHub$git_last_commit_date
    [1] "2022-01-20"
    
    $Packages$ExperimentHub$Hash
    [1] "7baef74c130e13932918e975a8878f46"
    
    
    $Packages$commonmark
    $Packages$commonmark$Package
    [1] "commonmark"
    
    $Packages$commonmark$Version
    [1] "1.8.0"
    
    $Packages$commonmark$Source
    [1] "Repository"
    
    $Packages$commonmark$Repository
    [1] "CRAN"
    
    $Packages$commonmark$Hash
    [1] "2ba81b120c1655ab696c935ef33ea716"
    
    
    $Packages$munsell
    $Packages$munsell$Package
    [1] "munsell"
    
    $Packages$munsell$Version
    [1] "0.5.0"
    
    $Packages$munsell$Source
    [1] "Repository"
    
    $Packages$munsell$Imports
    [1] "colorspace" "methods"   
    
    $Packages$munsell$Repository
    [1] "CRAN"
    
    $Packages$munsell$Hash
    [1] "6dfe8bf774944bd5595785e3229d8771"
    
    
    $Packages$cellranger
    $Packages$cellranger$Package
    [1] "cellranger"
    
    $Packages$cellranger$Version
    [1] "1.1.0"
    
    $Packages$cellranger$Source
    [1] "Repository"
    
    $Packages$cellranger$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$cellranger$Imports
    [1] "rematch" "tibble" 
    
    $Packages$cellranger$Repository
    [1] "CRAN"
    
    $Packages$cellranger$Hash
    [1] "f61dbaec772ccd2e17705c1e872e9e7c"
    
    
    $Packages$ggsci
    $Packages$ggsci$Package
    [1] "ggsci"
    
    $Packages$ggsci$Version
    [1] "2.9"
    
    $Packages$ggsci$Source
    [1] "Repository"
    
    $Packages$ggsci$Depends
    [1] "R (>= 3.0.2)"
    
    $Packages$ggsci$Imports
    [1] "grDevices"          "scales"             "ggplot2 (>= 2.0.0)"
    
    $Packages$ggsci$Repository
    [1] "CRAN"
    
    $Packages$ggsci$Hash
    [1] "81ccb8213ed592598210afd10c3a5936"
    
    
    $Packages$R.methodsS3
    $Packages$R.methodsS3$Package
    [1] "R.methodsS3"
    
    $Packages$R.methodsS3$Version
    [1] "1.8.1"
    
    $Packages$R.methodsS3$Source
    [1] "Repository"
    
    $Packages$R.methodsS3$Depends
    [1] "R (>= 2.13.0)"
    
    $Packages$R.methodsS3$Imports
    [1] "utils"
    
    $Packages$R.methodsS3$Repository
    [1] "CRAN"
    
    $Packages$R.methodsS3$Hash
    [1] "4bf6453323755202d5909697b6f7c109"
    
    
    $Packages$hwriter
    $Packages$hwriter$Package
    [1] "hwriter"
    
    $Packages$hwriter$Version
    [1] "1.3.2"
    
    $Packages$hwriter$Source
    [1] "Repository"
    
    $Packages$hwriter$Depends
    [1] "R (>= 2.6.0)"
    
    $Packages$hwriter$Repository
    [1] "CRAN"
    
    $Packages$hwriter$Hash
    [1] "fd39c53c2fff477e78557114b5ff873b"
    
    
    $Packages$codetools
    $Packages$codetools$Package
    [1] "codetools"
    
    $Packages$codetools$Version
    [1] "0.2-18"
    
    $Packages$codetools$Source
    [1] "Repository"
    
    $Packages$codetools$Depends
    [1] "R (>= 2.1)"
    
    $Packages$codetools$Repository
    [1] "CRAN"
    
    $Packages$codetools$Hash
    [1] "019388fc48e48b3da0d3a76ff94608a8"
    
    
    $Packages$Biobase
    $Packages$Biobase$Package
    [1] "Biobase"
    
    $Packages$Biobase$Version
    [1] "2.54.0"
    
    $Packages$Biobase$Source
    [1] "Bioconductor"
    
    $Packages$Biobase$Depends
    [1] "R (>= 2.10)"              "BiocGenerics (>= 0.27.1)"
    [3] "utils"                   
    
    $Packages$Biobase$Imports
    [1] "methods"
    
    $Packages$Biobase$git_url
    [1] "https://git.bioconductor.org/packages/Biobase"
    
    $Packages$Biobase$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$Biobase$git_last_commit
    [1] "8215d76"
    
    $Packages$Biobase$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$Biobase$Hash
    [1] "767057b0e2897e2a43b59718be888ccd"
    
    
    $Packages$GenomeInfoDb
    $Packages$GenomeInfoDb$Package
    [1] "GenomeInfoDb"
    
    $Packages$GenomeInfoDb$Version
    [1] "1.30.1"
    
    $Packages$GenomeInfoDb$Source
    [1] "Bioconductor"
    
    $Packages$GenomeInfoDb$Depends
    [1] "R (>= 4.0.0)"             "methods"                 
    [3] "BiocGenerics (>= 0.37.0)" "S4Vectors (>= 0.25.12)"  
    [5] "IRanges (>= 2.13.12)"    
    
    $Packages$GenomeInfoDb$Imports
    [1] "stats"            "stats4"           "utils"            "RCurl"           
    [5] "GenomeInfoDbData"
    
    $Packages$GenomeInfoDb$git_url
    [1] "https://git.bioconductor.org/packages/GenomeInfoDb"
    
    $Packages$GenomeInfoDb$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$GenomeInfoDb$git_last_commit
    [1] "bf8b385"
    
    $Packages$GenomeInfoDb$git_last_commit_date
    [1] "2022-01-27"
    
    $Packages$GenomeInfoDb$Hash
    [1] "30dc66e49ea68fd50f5d3fe4b82f0533"
    
    
    $Packages$themis
    $Packages$themis$Package
    [1] "themis"
    
    $Packages$themis$Version
    [1] "0.1.4"
    
    $Packages$themis$Source
    [1] "Repository"
    
    $Packages$themis$Depends
    [1] "R (>= 2.10)"         "recipes (>= 0.1.15)"
    
    $Packages$themis$Imports
    [1] "dplyr"               "generics (>= 0.1.0)" "purrr"              
    [4] "RANN"                "rlang"               "ROSE"               
    [7] "tibble"              "unbalanced"          "withr"              
    
    $Packages$themis$Repository
    [1] "CRAN"
    
    $Packages$themis$Hash
    [1] "42b548d48b359286942c20fcd78c348f"
    
    
    $Packages$MultiAssayExperiment
    $Packages$MultiAssayExperiment$Package
    [1] "MultiAssayExperiment"
    
    $Packages$MultiAssayExperiment$Version
    [1] "1.20.0"
    
    $Packages$MultiAssayExperiment$Source
    [1] "Bioconductor"
    
    $Packages$MultiAssayExperiment$Depends
    [1] "R (>= 4.0.0)"                     "SummarizedExperiment (>= 1.3.81)"
    
    $Packages$MultiAssayExperiment$Imports
    [1] "methods"                    "GenomicRanges (>= 1.25.93)"
    [3] "BiocGenerics"               "S4Vectors (>= 0.23.19)"    
    [5] "IRanges"                    "Biobase"                   
    [7] "stats"                      "tidyr"                     
    [9] "utils"                     
    
    $Packages$MultiAssayExperiment$git_url
    [1] "https://git.bioconductor.org/packages/MultiAssayExperiment"
    
    $Packages$MultiAssayExperiment$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$MultiAssayExperiment$git_last_commit
    [1] "c543a7c"
    
    $Packages$MultiAssayExperiment$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$MultiAssayExperiment$Hash
    [1] "e218c8c5a55951f37fc767d4307b395c"
    
    
    $Packages$vipor
    $Packages$vipor$Package
    [1] "vipor"
    
    $Packages$vipor$Version
    [1] "0.4.5"
    
    $Packages$vipor$Source
    [1] "Repository"
    
    $Packages$vipor$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$vipor$Imports
    [1] "stats"    "graphics"
    
    $Packages$vipor$Repository
    [1] "CRAN"
    
    $Packages$vipor$Hash
    [1] "ea85683da7f2bfa63a98dc6416892591"
    
    
    $Packages$sys
    $Packages$sys$Package
    [1] "sys"
    
    $Packages$sys$Version
    [1] "3.4"
    
    $Packages$sys$Source
    [1] "Repository"
    
    $Packages$sys$Repository
    [1] "CRAN"
    
    $Packages$sys$Hash
    [1] "b227d13e29222b4574486cfcbde077fa"
    
    
    $Packages$ontologyIndex
    $Packages$ontologyIndex$Package
    [1] "ontologyIndex"
    
    $Packages$ontologyIndex$Version
    [1] "2.7"
    
    $Packages$ontologyIndex$Source
    [1] "Repository"
    
    $Packages$ontologyIndex$Depends
    [1] "R (>= 3.1.0)"
    
    $Packages$ontologyIndex$Repository
    [1] "CRAN"
    
    $Packages$ontologyIndex$Hash
    [1] "ca4b9b931c659a989d3fc8e97b1b4525"
    
    
    $Packages$xtable
    $Packages$xtable$Package
    [1] "xtable"
    
    $Packages$xtable$Version
    [1] "1.8-4"
    
    $Packages$xtable$Source
    [1] "Repository"
    
    $Packages$xtable$Depends
    [1] "R (>= 2.10.0)"
    
    $Packages$xtable$Imports
    [1] "stats" "utils"
    
    $Packages$xtable$Repository
    [1] "CRAN"
    
    $Packages$xtable$Hash
    [1] "b8acdf8af494d9ec19ccb2481a9b11c2"
    
    
    $Packages$googlesheets4
    $Packages$googlesheets4$Package
    [1] "googlesheets4"
    
    $Packages$googlesheets4$Version
    [1] "1.0.0"
    
    $Packages$googlesheets4$Source
    [1] "Repository"
    
    $Packages$googlesheets4$Depends
    [1] "R (>= 3.3)"
    
    $Packages$googlesheets4$Imports
     [1] "cellranger"             "cli (>= 3.0.0)"         "curl"                  
     [4] "gargle (>= 1.2.0)"      "glue (>= 1.3.0)"        "googledrive (>= 2.0.0)"
     [7] "httr"                   "ids"                    "magrittr"              
    [10] "methods"                "purrr"                  "rematch2"              
    [13] "rlang (>= 0.4.11)"      "tibble (>= 2.1.1)"      "utils"                 
    [16] "vctrs (>= 0.2.3)"      
    
    $Packages$googlesheets4$Repository
    [1] "CRAN"
    
    $Packages$googlesheets4$Hash
    [1] "9a6564184dc4a81daea4f1d7ce357c6a"
    
    
    $Packages$formatR
    $Packages$formatR$Package
    [1] "formatR"
    
    $Packages$formatR$Version
    [1] "1.11"
    
    $Packages$formatR$Source
    [1] "Repository"
    
    $Packages$formatR$Depends
    [1] "R (>= 3.2.3)"
    
    $Packages$formatR$Repository
    [1] "CRAN"
    
    $Packages$formatR$Hash
    [1] "2590a6a868515a69f258640a51b724c9"
    
    
    $Packages$BiocManager
    $Packages$BiocManager$Package
    [1] "BiocManager"
    
    $Packages$BiocManager$Version
    [1] "1.30.16"
    
    $Packages$BiocManager$Source
    [1] "Repository"
    
    $Packages$BiocManager$Imports
    [1] "utils"
    
    $Packages$BiocManager$Repository
    [1] "CRAN"
    
    $Packages$BiocManager$Hash
    [1] "2fdca0877debdd4668190832cdee4c31"
    
    
    $Packages$farver
    $Packages$farver$Package
    [1] "farver"
    
    $Packages$farver$Version
    [1] "2.1.0"
    
    $Packages$farver$Source
    [1] "Repository"
    
    $Packages$farver$Repository
    [1] "CRAN"
    
    $Packages$farver$Hash
    [1] "c98eb5133d9cb9e1622b8691487f11bb"
    
    
    $Packages$FNN
    $Packages$FNN$Package
    [1] "FNN"
    
    $Packages$FNN$Version
    [1] "1.1.3"
    
    $Packages$FNN$Source
    [1] "Repository"
    
    $Packages$FNN$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$FNN$Repository
    [1] "CRAN"
    
    $Packages$FNN$Hash
    [1] "b56998fff55e4a4b4860ad6e8c67e0f9"
    
    
    $Packages$parallelly
    $Packages$parallelly$Package
    [1] "parallelly"
    
    $Packages$parallelly$Version
    [1] "1.30.0"
    
    $Packages$parallelly$Source
    [1] "Repository"
    
    $Packages$parallelly$Imports
    [1] "parallel" "tools"    "utils"   
    
    $Packages$parallelly$Repository
    [1] "CRAN"
    
    $Packages$parallelly$Hash
    [1] "67db13907a9cea89c118cf82d448799f"
    
    
    $Packages$AnnotationHub
    $Packages$AnnotationHub$Package
    [1] "AnnotationHub"
    
    $Packages$AnnotationHub$Version
    [1] "3.2.2"
    
    $Packages$AnnotationHub$Source
    [1] "Bioconductor"
    
    $Packages$AnnotationHub$Depends
    [1] "BiocGenerics (>= 0.15.10)" "BiocFileCache (>= 1.5.1)" 
    
    $Packages$AnnotationHub$Imports
     [1] "utils"                      "methods"                   
     [3] "grDevices"                  "RSQLite"                   
     [5] "BiocManager"                "BiocVersion"               
     [7] "curl"                       "rappdirs"                  
     [9] "AnnotationDbi (>= 1.31.19)" "S4Vectors"                 
    [11] "interactiveDisplayBase"     "httr"                      
    [13] "yaml"                       "dplyr"                     
    
    $Packages$AnnotationHub$git_url
    [1] "https://git.bioconductor.org/packages/AnnotationHub"
    
    $Packages$AnnotationHub$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$AnnotationHub$git_last_commit
    [1] "8e761e1"
    
    $Packages$AnnotationHub$git_last_commit_date
    [1] "2022-02-28"
    
    $Packages$AnnotationHub$Hash
    [1] "7630e5d6cf99e154acd3bfe8ae3b5b8d"
    
    
    $Packages$RANN
    $Packages$RANN$Package
    [1] "RANN"
    
    $Packages$RANN$Version
    [1] "2.6.1"
    
    $Packages$RANN$Source
    [1] "Repository"
    
    $Packages$RANN$Repository
    [1] "CRAN"
    
    $Packages$RANN$Hash
    [1] "d128ea05a972d3e67c6f39de52c72bd7"
    
    
    $Packages$aplot
    $Packages$aplot$Package
    [1] "aplot"
    
    $Packages$aplot$Version
    [1] "0.1.2"
    
    $Packages$aplot$Source
    [1] "Repository"
    
    $Packages$aplot$Imports
    [1] "ggfun (>= 0.0.4)" "ggplot2"          "ggplotify"        "patchwork"       
    [5] "magrittr"         "methods"          "utils"            "yulab.utils"     
    
    $Packages$aplot$Repository
    [1] "CRAN"
    
    $Packages$aplot$Hash
    [1] "10f88b6bfa5c7ab877d1ece725aaa66d"
    
    
    $Packages$askpass
    $Packages$askpass$Package
    [1] "askpass"
    
    $Packages$askpass$Version
    [1] "1.1"
    
    $Packages$askpass$Source
    [1] "Repository"
    
    $Packages$askpass$Imports
    [1] "sys (>= 2.1)"
    
    $Packages$askpass$Repository
    [1] "CRAN"
    
    $Packages$askpass$Hash
    [1] "e8a22846fff485f0be3770c2da758713"
    
    
    $Packages$BBmisc
    $Packages$BBmisc$Package
    [1] "BBmisc"
    
    $Packages$BBmisc$Version
    [1] "1.11"
    
    $Packages$BBmisc$Source
    [1] "Repository"
    
    $Packages$BBmisc$Imports
    [1] "utils"                "methods"              "stats"               
    [4] "checkmate (>= 1.8.0)"
    
    $Packages$BBmisc$Repository
    [1] "CRAN"
    
    $Packages$BBmisc$Hash
    [1] "b2ff0878d07259998cd84d739a83b3d9"
    
    
    $Packages$duckdb
    $Packages$duckdb$Package
    [1] "duckdb"
    
    $Packages$duckdb$Version
    [1] "0.3.2-1"
    
    $Packages$duckdb$Source
    [1] "Repository"
    
    $Packages$duckdb$Depends
    [1] "DBI"          "R (>= 3.6.0)"
    
    $Packages$duckdb$Imports
    [1] "methods" "utils"  
    
    $Packages$duckdb$Repository
    [1] "CRAN"
    
    $Packages$duckdb$Hash
    [1] "29d4f41f42b2cf0feecd265f8fd8412f"
    
    
    $Packages$ggtree
    $Packages$ggtree$Package
    [1] "ggtree"
    
    $Packages$ggtree$Version
    [1] "3.2.1"
    
    $Packages$ggtree$Source
    [1] "Bioconductor"
    
    $Packages$ggtree$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$ggtree$Imports
     [1] "ape"                 "aplot (>= 0.0.4)"    "dplyr"              
     [4] "ggplot2 (>= 3.0.0)"  "grid"                "magrittr"           
     [7] "methods"             "purrr"               "rlang"              
    [10] "ggfun"               "yulab.utils"         "tidyr"              
    [13] "tidytree (>= 0.2.6)" "treeio (>= 1.8.0)"   "utils"              
    [16] "scales"             
    
    $Packages$ggtree$RemoteType
    [1] "bioconductor"
    
    $Packages$ggtree$Remotes
    [1] "GuangchuangYu/treeio"
    
    $Packages$ggtree$Remotes
    [1] "GuangchuangYu/treeio"
    
    $Packages$ggtree$git_url
    [1] "https://git.bioconductor.org/packages/ggtree"
    
    $Packages$ggtree$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$ggtree$git_last_commit
    [1] "d3747e6"
    
    $Packages$ggtree$git_last_commit_date
    [1] "2021-11-14"
    
    $Packages$ggtree$Hash
    [1] "5711c057a04e53ed1c70909939dd9ad9"
    
    
    $Packages$plogr
    $Packages$plogr$Package
    [1] "plogr"
    
    $Packages$plogr$Version
    [1] "0.2.0"
    
    $Packages$plogr$Source
    [1] "Repository"
    
    $Packages$plogr$Repository
    [1] "CRAN"
    
    $Packages$plogr$Hash
    [1] "09eb987710984fc2905c7129c7d85e65"
    
    
    $Packages$GenomicRanges
    $Packages$GenomicRanges$Package
    [1] "GenomicRanges"
    
    $Packages$GenomicRanges$Version
    [1] "1.46.1"
    
    $Packages$GenomicRanges$Source
    [1] "Bioconductor"
    
    $Packages$GenomicRanges$Depends
    [1] "R (>= 4.0.0)"             "methods"                 
    [3] "stats4"                   "BiocGenerics (>= 0.37.0)"
    [5] "S4Vectors (>= 0.27.12)"   "IRanges (>= 2.23.9)"     
    [7] "GenomeInfoDb (>= 1.15.2)"
    
    $Packages$GenomicRanges$Imports
    [1] "utils"               "stats"               "XVector (>= 0.29.2)"
    
    $Packages$GenomicRanges$LinkingTo
    [1] "S4Vectors" "IRanges"  
    
    $Packages$GenomicRanges$git_url
    [1] "https://git.bioconductor.org/packages/GenomicRanges"
    
    $Packages$GenomicRanges$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$GenomicRanges$git_last_commit
    [1] "e422642"
    
    $Packages$GenomicRanges$git_last_commit_date
    [1] "2021-11-16"
    
    $Packages$GenomicRanges$Hash
    [1] "4d9cdd0c25e0c40075a2756aa78f4960"
    
    
    $Packages$BiocIO
    $Packages$BiocIO$Package
    [1] "BiocIO"
    
    $Packages$BiocIO$Version
    [1] "1.4.0"
    
    $Packages$BiocIO$Source
    [1] "Bioconductor"
    
    $Packages$BiocIO$Depends
    [1] "R (>= 4.0)"
    
    $Packages$BiocIO$Imports
    [1] "BiocGenerics" "S4Vectors"    "methods"      "tools"       
    
    $Packages$BiocIO$git_url
    [1] "https://git.bioconductor.org/packages/BiocIO"
    
    $Packages$BiocIO$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocIO$git_last_commit
    [1] "c335932"
    
    $Packages$BiocIO$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$BiocIO$Hash
    [1] "10f50bd6daf34cd11d222befa8829635"
    
    
    $Packages$base64url
    $Packages$base64url$Package
    [1] "base64url"
    
    $Packages$base64url$Version
    [1] "1.4"
    
    $Packages$base64url$Source
    [1] "Repository"
    
    $Packages$base64url$Imports
    [1] "backports (>= 1.1.0)"
    
    $Packages$base64url$Repository
    [1] "CRAN"
    
    $Packages$base64url$Hash
    [1] "0c54cf3a08cc0e550fbd64ad33166143"
    
    
    $Packages$patchwork
    $Packages$patchwork$Package
    [1] "patchwork"
    
    $Packages$patchwork$Version
    [1] "1.1.1"
    
    $Packages$patchwork$Source
    [1] "Repository"
    
    $Packages$patchwork$Imports
    [1] "ggplot2 (>= 3.0.0)" "gtable"             "grid"              
    [4] "stats"              "grDevices"          "utils"             
    [7] "graphics"          
    
    $Packages$patchwork$Repository
    [1] "CRAN"
    
    $Packages$patchwork$Hash
    [1] "c446b30cb33ec125ff02588b60660ccb"
    
    
    $Packages$dtplyr
    $Packages$dtplyr$Package
    [1] "dtplyr"
    
    $Packages$dtplyr$Version
    [1] "1.2.1"
    
    $Packages$dtplyr$Source
    [1] "Repository"
    
    $Packages$dtplyr$Depends
    [1] "R (>= 3.3)"
    
    $Packages$dtplyr$Imports
     [1] "crayon"                 "data.table (>= 1.13.0)" "dplyr (>= 1.0.3)"      
     [4] "ellipsis"               "glue"                   "lifecycle"             
     [7] "rlang"                  "tibble"                 "tidyselect"            
    [10] "vctrs"                 
    
    $Packages$dtplyr$Repository
    [1] "CRAN"
    
    $Packages$dtplyr$Hash
    [1] "f5d195cd5fcc0a77499d9da698ef2ea3"
    
    
    $Packages$tibble
    $Packages$tibble$Package
    [1] "tibble"
    
    $Packages$tibble$Version
    [1] "3.1.6"
    
    $Packages$tibble$Source
    [1] "Repository"
    
    $Packages$tibble$Depends
    [1] "R (>= 3.1.0)"
    
    $Packages$tibble$Imports
     [1] "ellipsis (>= 0.3.2)"  "fansi (>= 0.4.0)"     "lifecycle (>= 1.0.0)"
     [4] "magrittr"             "methods"              "pillar (>= 1.6.2)"   
     [7] "pkgconfig"            "rlang (>= 0.4.3)"     "utils"               
    [10] "vctrs (>= 0.3.8)"    
    
    $Packages$tibble$Repository
    [1] "CRAN"
    
    $Packages$tibble$Hash
    [1] "8a8f02d1934dfd6431c671361510dd0b"
    
    
    $Packages$futile.options
    $Packages$futile.options$Package
    [1] "futile.options"
    
    $Packages$futile.options$Version
    [1] "1.0.1"
    
    $Packages$futile.options$Source
    [1] "Repository"
    
    $Packages$futile.options$Depends
    [1] "R (>= 2.8.0)"
    
    $Packages$futile.options$Repository
    [1] "CRAN"
    
    $Packages$futile.options$Hash
    [1] "0d9bf02413ddc2bbe8da9ce369dcdd2b"
    
    
    $Packages$cluster
    $Packages$cluster$Package
    [1] "cluster"
    
    $Packages$cluster$Version
    [1] "2.1.2"
    
    $Packages$cluster$Source
    [1] "Repository"
    
    $Packages$cluster$Depends
    [1] "R (>= 3.4.0)"
    
    $Packages$cluster$Imports
    [1] "graphics"  "grDevices" "stats"     "utils"    
    
    $Packages$cluster$Repository
    [1] "CRAN"
    
    $Packages$cluster$Hash
    [1] "ce49bfe5bc0b3ecd43a01fe1b01c2243"
    
    
    $Packages$future.apply
    $Packages$future.apply$Package
    [1] "future.apply"
    
    $Packages$future.apply$Version
    [1] "1.8.1"
    
    $Packages$future.apply$Source
    [1] "Repository"
    
    $Packages$future.apply$Depends
    [1] "R (>= 3.2.0)"       "future (>= 1.21.0)"
    
    $Packages$future.apply$Imports
    [1] "globals (>= 0.14.0)" "parallel"            "utils"              
    
    $Packages$future.apply$Repository
    [1] "CRAN"
    
    $Packages$future.apply$Hash
    [1] "f568ce73d3d59582b0f7babd0eb33d07"
    
    
    $Packages$praise
    $Packages$praise$Package
    [1] "praise"
    
    $Packages$praise$Version
    [1] "1.0.0"
    
    $Packages$praise$Source
    [1] "Repository"
    
    $Packages$praise$Repository
    [1] "CRAN"
    
    $Packages$praise$Hash
    [1] "a555924add98c99d2f411e37e7d25e9f"
    
    
    $Packages$Matrix
    $Packages$Matrix$Package
    [1] "Matrix"
    
    $Packages$Matrix$Version
    [1] "1.4-0"
    
    $Packages$Matrix$Source
    [1] "Repository"
    
    $Packages$Matrix$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$Matrix$Imports
    [1] "methods"  "graphics" "grid"     "stats"    "utils"    "lattice" 
    
    $Packages$Matrix$Repository
    [1] "CRAN"
    
    $Packages$Matrix$Hash
    [1] "130c0caba175739d98f2963c6a407cf6"
    
    
    $Packages$tidytree
    $Packages$tidytree$Package
    [1] "tidytree"
    
    $Packages$tidytree$Version
    [1] "0.3.9"
    
    $Packages$tidytree$Source
    [1] "Repository"
    
    $Packages$tidytree$Depends
    [1] "R (>= 3.4.0)"
    
    $Packages$tidytree$Imports
     [1] "ape"                    "dplyr"                  "lazyeval"              
     [4] "magrittr"               "methods"                "rlang"                 
     [7] "tibble"                 "tidyr"                  "tidyselect"            
    [10] "yulab.utils (>= 0.0.4)" "pillar"                
    
    $Packages$tidytree$Repository
    [1] "CRAN"
    
    $Packages$tidytree$Hash
    [1] "bf9d994b6c4a8e448618feac4bbff2d4"
    
    
    $Packages$ellipsis
    $Packages$ellipsis$Package
    [1] "ellipsis"
    
    $Packages$ellipsis$Version
    [1] "0.3.2"
    
    $Packages$ellipsis$Source
    [1] "Repository"
    
    $Packages$ellipsis$Depends
    [1] "R (>= 3.2)"
    
    $Packages$ellipsis$Imports
    [1] "rlang (>= 0.3.0)"
    
    $Packages$ellipsis$Repository
    [1] "CRAN"
    
    $Packages$ellipsis$Hash
    [1] "bb0eec2fe32e88d9e2836c2f73ea2077"
    
    
    $Packages$prettyunits
    $Packages$prettyunits$Package
    [1] "prettyunits"
    
    $Packages$prettyunits$Version
    [1] "1.1.1"
    
    $Packages$prettyunits$Source
    [1] "Repository"
    
    $Packages$prettyunits$Repository
    [1] "CRAN"
    
    $Packages$prettyunits$Hash
    [1] "95ef9167b75dde9d2ccc3c7528393e7e"
    
    
    $Packages$SQUAREM
    $Packages$SQUAREM$Package
    [1] "SQUAREM"
    
    $Packages$SQUAREM$Version
    [1] "2021.1"
    
    $Packages$SQUAREM$Source
    [1] "Repository"
    
    $Packages$SQUAREM$Depends
    [1] "R (>= 3.0)"
    
    $Packages$SQUAREM$Repository
    [1] "CRAN"
    
    $Packages$SQUAREM$Hash
    [1] "0cf10dab0d023d5b46a5a14387556891"
    
    
    $Packages$lubridate
    $Packages$lubridate$Package
    [1] "lubridate"
    
    $Packages$lubridate$Version
    [1] "1.8.0"
    
    $Packages$lubridate$Source
    [1] "Repository"
    
    $Packages$lubridate$Depends
    [1] "methods"    "R (>= 3.2)"
    
    $Packages$lubridate$Imports
    [1] "generics"
    
    $Packages$lubridate$LinkingTo
    [1] "cpp11 (>= 0.2.7)"
    
    $Packages$lubridate$Repository
    [1] "CRAN"
    
    $Packages$lubridate$Hash
    [1] "2ff5eedb6ee38fb1b81205c73be1be5a"
    
    
    $Packages$googledrive
    $Packages$googledrive$Package
    [1] "googledrive"
    
    $Packages$googledrive$Version
    [1] "2.0.0"
    
    $Packages$googledrive$Source
    [1] "Repository"
    
    $Packages$googledrive$Depends
    [1] "R (>= 3.3)"
    
    $Packages$googledrive$Imports
     [1] "cli (>= 3.0.0)"    "gargle (>= 1.2.0)" "glue (>= 1.4.2)"  
     [4] "httr"              "jsonlite"          "lifecycle"        
     [7] "magrittr"          "pillar"            "purrr (>= 0.2.3)" 
    [10] "rlang (>= 0.4.9)"  "tibble (>= 2.0.0)" "utils"            
    [13] "uuid"              "vctrs (>= 0.3.0)"  "withr"            
    
    $Packages$googledrive$Repository
    [1] "CRAN"
    
    $Packages$googledrive$Hash
    [1] "c3a25adbbfbb03f12e6f88c5fb1f3024"
    
    
    $Packages$reprex
    $Packages$reprex$Package
    [1] "reprex"
    
    $Packages$reprex$Version
    [1] "2.0.1"
    
    $Packages$reprex$Source
    [1] "Repository"
    
    $Packages$reprex$Depends
    [1] "R (>= 3.3)"
    
    $Packages$reprex$Imports
     [1] "callr (>= 3.6.0)" "cli (>= 2.3.1)"   "clipr (>= 0.4.0)" "fs"              
     [5] "glue"             "knitr (>= 1.23)"  "rlang (>= 0.4.0)" "rmarkdown"       
     [9] "rstudioapi"       "utils"            "withr (>= 2.3.0)"
    
    $Packages$reprex$Repository
    [1] "CRAN"
    
    $Packages$reprex$Hash
    [1] "911d101becedc0fde495bd910984bdc8"
    
    
    $Packages$rematch
    $Packages$rematch$Package
    [1] "rematch"
    
    $Packages$rematch$Version
    [1] "1.0.1"
    
    $Packages$rematch$Source
    [1] "Repository"
    
    $Packages$rematch$Repository
    [1] "CRAN"
    
    $Packages$rematch$Hash
    [1] "c66b930d20bb6d858cd18e1cebcfae5c"
    
    
    $Packages$gitcreds
    $Packages$gitcreds$Package
    [1] "gitcreds"
    
    $Packages$gitcreds$Version
    [1] "0.1.1"
    
    $Packages$gitcreds$Source
    [1] "Repository"
    
    $Packages$gitcreds$Repository
    [1] "CRAN"
    
    $Packages$gitcreds$Hash
    [1] "f3aefccc1cc50de6338146b62f115de8"
    
    
    $Packages$igraph
    $Packages$igraph$Package
    [1] "igraph"
    
    $Packages$igraph$Version
    [1] "1.2.11"
    
    $Packages$igraph$Source
    [1] "Repository"
    
    $Packages$igraph$Depends
    [1] "methods"
    
    $Packages$igraph$Imports
    [1] "graphics"             "grDevices"            "magrittr"            
    [4] "Matrix"               "pkgconfig (>= 2.0.0)" "stats"               
    [7] "utils"               
    
    $Packages$igraph$Repository
    [1] "CRAN"
    
    $Packages$igraph$Hash
    [1] "1d10cd31c2979f9c819ffe4d16b9dc2b"
    
    
    $Packages$multtest
    $Packages$multtest$Package
    [1] "multtest"
    
    $Packages$multtest$Version
    [1] "2.50.0"
    
    $Packages$multtest$Source
    [1] "Bioconductor"
    
    $Packages$multtest$Depends
    [1] "R (>= 2.10)"  "methods"      "BiocGenerics" "Biobase"     
    
    $Packages$multtest$Imports
    [1] "survival" "MASS"     "stats4"  
    
    $Packages$multtest$git_url
    [1] "https://git.bioconductor.org/packages/multtest"
    
    $Packages$multtest$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$multtest$git_last_commit
    [1] "1de9664"
    
    $Packages$multtest$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$multtest$Hash
    [1] "096b9b9a42d79a82431eb5d9ed9f0da3"
    
    
    $Packages$RcppEigen
    $Packages$RcppEigen$Package
    [1] "RcppEigen"
    
    $Packages$RcppEigen$Version
    [1] "0.3.3.9.1"
    
    $Packages$RcppEigen$Source
    [1] "Repository"
    
    $Packages$RcppEigen$Imports
    [1] "Matrix (>= 1.1-0)" "Rcpp (>= 0.11.0)"  "stats"            
    [4] "utils"            
    
    $Packages$RcppEigen$LinkingTo
    [1] "Rcpp"
    
    $Packages$RcppEigen$Repository
    [1] "CRAN"
    
    $Packages$RcppEigen$Hash
    [1] "ddfa72a87fdf4c80466a20818be91d00"
    
    
    $Packages$fgsea
    $Packages$fgsea$Package
    [1] "fgsea"
    
    $Packages$fgsea$Version
    [1] "1.20.0"
    
    $Packages$fgsea$Source
    [1] "Bioconductor"
    
    $Packages$fgsea$Depends
    [1] "R (>= 3.3)"
    
    $Packages$fgsea$Imports
     [1] "Rcpp"               "data.table"         "BiocParallel"      
     [4] "stats"              "ggplot2 (>= 2.2.0)" "gridExtra"         
     [7] "grid"               "fastmatch"          "Matrix"            
    [10] "utils"             
    
    $Packages$fgsea$LinkingTo
    [1] "Rcpp" "BH"  
    
    $Packages$fgsea$git_url
    [1] "https://git.bioconductor.org/packages/fgsea"
    
    $Packages$fgsea$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$fgsea$git_last_commit
    [1] "b704f81"
    
    $Packages$fgsea$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$fgsea$Hash
    [1] "a4913104166e73957fa32d2e595cceef"
    
    
    $Packages$unbalanced
    $Packages$unbalanced$Package
    [1] "unbalanced"
    
    $Packages$unbalanced$Version
    [1] "2.0"
    
    $Packages$unbalanced$Source
    [1] "Repository"
    
    $Packages$unbalanced$Depends
    [1] "mlr"        "foreach"    "doParallel"
    
    $Packages$unbalanced$Imports
    [1] "FNN"  "RANN"
    
    $Packages$unbalanced$Repository
    [1] "CRAN"
    
    $Packages$unbalanced$Hash
    [1] "ed45de16b03027f36b2bc2e8c84c2533"
    
    
    $Packages$dada2
    $Packages$dada2$Package
    [1] "dada2"
    
    $Packages$dada2$Version
    [1] "1.22.0"
    
    $Packages$dada2$Source
    [1] "Bioconductor"
    
    $Packages$dada2$Depends
    [1] "R (>= 3.4.0)"       "Rcpp (>= 0.12.0)"   "methods (>= 3.4.0)"
    
    $Packages$dada2$Imports
    [1] "Biostrings (>= 2.42.1)"   "ggplot2 (>= 2.1.0)"      
    [3] "reshape2 (>= 1.4.1)"      "ShortRead (>= 1.32.0)"   
    [5] "RcppParallel (>= 4.3.0)"  "parallel (>= 3.2.0)"     
    [7] "IRanges (>= 2.6.0)"       "XVector (>= 0.16.0)"     
    [9] "BiocGenerics (>= 0.22.0)"
    
    $Packages$dada2$LinkingTo
    [1] "Rcpp"         "RcppParallel"
    
    $Packages$dada2$git_url
    [1] "https://git.bioconductor.org/packages/dada2"
    
    $Packages$dada2$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$dada2$git_last_commit
    [1] "3abe06c"
    
    $Packages$dada2$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$dada2$Hash
    [1] "03bd5c0a175bcc98e5c552941c60c70e"
    
    
    $Packages$gargle
    $Packages$gargle$Package
    [1] "gargle"
    
    $Packages$gargle$Version
    [1] "1.2.0"
    
    $Packages$gargle$Source
    [1] "Repository"
    
    $Packages$gargle$Depends
    [1] "R (>= 3.3)"
    
    $Packages$gargle$Imports
     [1] "cli (>= 3.0.0)"   "fs (>= 1.3.1)"    "glue (>= 1.3.0)"  "httr (>= 1.4.0)" 
     [5] "jsonlite"         "rappdirs"         "rlang (>= 0.4.9)" "rstudioapi"      
     [9] "stats"            "utils"            "withr"           
    
    $Packages$gargle$Repository
    [1] "CRAN"
    
    $Packages$gargle$Hash
    [1] "9d234e6a87a6f8181792de6dc4a00e39"
    
    
    $Packages$testthat
    $Packages$testthat$Package
    [1] "testthat"
    
    $Packages$testthat$Version
    [1] "3.1.2"
    
    $Packages$testthat$Source
    [1] "Repository"
    
    $Packages$testthat$Depends
    [1] "R (>= 3.1)"
    
    $Packages$testthat$Imports
     [1] "brio"                "callr (>= 3.5.1)"    "cli (>= 2.2.0)"     
     [4] "crayon (>= 1.3.4)"   "desc"                "digest"             
     [7] "ellipsis (>= 0.2.0)" "evaluate"            "jsonlite"           
    [10] "lifecycle"           "magrittr"            "methods"            
    [13] "pkgload"             "praise"              "processx"           
    [16] "ps (>= 1.3.4)"       "R6 (>= 2.2.0)"       "rlang (>= 0.4.9)"   
    [19] "utils"               "waldo (>= 0.2.4)"    "withr (>= 2.4.3)"   
    
    $Packages$testthat$Repository
    [1] "CRAN"
    
    $Packages$testthat$Hash
    [1] "32454e5780e8dbe31e4b61b13d8918fe"
    
    
    $Packages$getopt
    $Packages$getopt$Package
    [1] "getopt"
    
    $Packages$getopt$Version
    [1] "1.20.3"
    
    $Packages$getopt$Source
    [1] "Repository"
    
    $Packages$getopt$Imports
    [1] "stats"
    
    $Packages$getopt$Repository
    [1] "CRAN"
    
    $Packages$getopt$Hash
    [1] "ad68e3263f0bc9269aab8c2039440117"
    
    
    $Packages$htmltools
    $Packages$htmltools$Package
    [1] "htmltools"
    
    $Packages$htmltools$Version
    [1] "0.5.2"
    
    $Packages$htmltools$Source
    [1] "Repository"
    
    $Packages$htmltools$Depends
    [1] "R (>= 2.14.1)"
    
    $Packages$htmltools$Imports
    [1] "utils"             "digest"            "grDevices"        
    [4] "base64enc"         "rlang (>= 0.4.10)" "fastmap"          
    
    $Packages$htmltools$Repository
    [1] "CRAN"
    
    $Packages$htmltools$Hash
    [1] "526c484233f42522278ab06fb185cb26"
    
    
    $Packages$BiocFileCache
    $Packages$BiocFileCache$Package
    [1] "BiocFileCache"
    
    $Packages$BiocFileCache$Version
    [1] "2.2.1"
    
    $Packages$BiocFileCache$Source
    [1] "Bioconductor"
    
    $Packages$BiocFileCache$Depends
    [1] "R (>= 3.4.0)"      "dbplyr (>= 1.0.0)"
    
    $Packages$BiocFileCache$Imports
     [1] "methods"  "stats"    "utils"    "dplyr"    "RSQLite"  "DBI"     
     [7] "rappdirs" "filelock" "curl"     "httr"    
    
    $Packages$BiocFileCache$git_url
    [1] "https://git.bioconductor.org/packages/BiocFileCache"
    
    $Packages$BiocFileCache$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocFileCache$git_last_commit
    [1] "cc91212"
    
    $Packages$BiocFileCache$git_last_commit_date
    [1] "2022-01-20"
    
    $Packages$BiocFileCache$Hash
    [1] "565eceff974df21115c8db7b4163ad9f"
    
    
    $Packages$piggyback
    $Packages$piggyback$Package
    [1] "piggyback"
    
    $Packages$piggyback$Version
    [1] "0.1.1"
    
    $Packages$piggyback$Source
    [1] "Repository"
    
    $Packages$piggyback$Imports
    [1] "gh"         "httr"       "jsonlite"   "fs"         "crayon"    
    [6] "clisymbols" "lubridate"  "memoise"   
    
    $Packages$piggyback$Repository
    [1] "CRAN"
    
    $Packages$piggyback$Hash
    [1] "0cd602e7eb4288f3e78fb4eba0507ead"
    
    
    $Packages$yaml
    $Packages$yaml$Package
    [1] "yaml"
    
    $Packages$yaml$Version
    [1] "2.3.5"
    
    $Packages$yaml$Source
    [1] "Repository"
    
    $Packages$yaml$Repository
    [1] "CRAN"
    
    $Packages$yaml$Hash
    [1] "458bb38374d73bf83b1bb85e353da200"
    
    
    $Packages$utf8
    $Packages$utf8$Package
    [1] "utf8"
    
    $Packages$utf8$Version
    [1] "1.2.2"
    
    $Packages$utf8$Source
    [1] "Repository"
    
    $Packages$utf8$Depends
    [1] "R (>= 2.10)"
    
    $Packages$utf8$Repository
    [1] "CRAN"
    
    $Packages$utf8$Hash
    [1] "c9c462b759a5cc844ae25b5942654d13"
    
    
    $Packages$interactiveDisplayBase
    $Packages$interactiveDisplayBase$Package
    [1] "interactiveDisplayBase"
    
    $Packages$interactiveDisplayBase$Version
    [1] "1.32.0"
    
    $Packages$interactiveDisplayBase$Source
    [1] "Bioconductor"
    
    $Packages$interactiveDisplayBase$Depends
    [1] "R (>= 2.10)"  "methods"      "BiocGenerics"
    
    $Packages$interactiveDisplayBase$Imports
    [1] "shiny" "DT"   
    
    $Packages$interactiveDisplayBase$git_url
    [1] "https://git.bioconductor.org/packages/interactiveDisplayBase"
    
    $Packages$interactiveDisplayBase$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$interactiveDisplayBase$git_last_commit
    [1] "0f88b2a"
    
    $Packages$interactiveDisplayBase$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$interactiveDisplayBase$Hash
    [1] "011d5724a7ccdd3c3e8788701242666d"
    
    
    $Packages$XML
    $Packages$XML$Package
    [1] "XML"
    
    $Packages$XML$Version
    [1] "3.99-0.9"
    
    $Packages$XML$Source
    [1] "Repository"
    
    $Packages$XML$Depends
    [1] "R (>= 4.0.0)" "methods"      "utils"       
    
    $Packages$XML$Repository
    [1] "CRAN"
    
    $Packages$XML$Hash
    [1] "6e75b3884423e21eda89b424a98e091d"
    
    
    $Packages$withr
    $Packages$withr$Package
    [1] "withr"
    
    $Packages$withr$Version
    [1] "2.5.0"
    
    $Packages$withr$Source
    [1] "Repository"
    
    $Packages$withr$Depends
    [1] "R (>= 3.2.0)"
    
    $Packages$withr$Imports
    [1] "graphics"  "grDevices" "stats"    
    
    $Packages$withr$Repository
    [1] "CRAN"
    
    $Packages$withr$Hash
    [1] "c0e49a9760983e81e55cdd9be92e7182"
    
    
    $Packages$scuttle
    $Packages$scuttle$Package
    [1] "scuttle"
    
    $Packages$scuttle$Version
    [1] "1.4.0"
    
    $Packages$scuttle$Source
    [1] "Bioconductor"
    
    $Packages$scuttle$Depends
    [1] "SingleCellExperiment"
    
    $Packages$scuttle$Imports
     [1] "methods"              "utils"                "stats"               
     [4] "Matrix"               "Rcpp"                 "BiocGenerics"        
     [7] "S4Vectors"            "BiocParallel"         "GenomicRanges"       
    [10] "SummarizedExperiment" "DelayedArray"         "DelayedMatrixStats"  
    [13] "beachmat"            
    
    $Packages$scuttle$LinkingTo
    [1] "Rcpp"     "beachmat"
    
    $Packages$scuttle$git_url
    [1] "https://git.bioconductor.org/packages/scuttle"
    
    $Packages$scuttle$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$scuttle$git_last_commit
    [1] "b335263"
    
    $Packages$scuttle$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$scuttle$Hash
    [1] "1c612f905434c933d27f508f294da9e0"
    
    
    $Packages$tidyverse
    $Packages$tidyverse$Package
    [1] "tidyverse"
    
    $Packages$tidyverse$Version
    [1] "1.3.1"
    
    $Packages$tidyverse$Source
    [1] "Repository"
    
    $Packages$tidyverse$Depends
    [1] "R (>= 3.3)"
    
    $Packages$tidyverse$Imports
     [1] "broom (>= 0.7.6)"         "cli (>= 2.4.0)"          
     [3] "crayon (>= 1.4.1)"        "dbplyr (>= 2.1.1)"       
     [5] "dplyr (>= 1.0.5)"         "dtplyr (>= 1.1.0)"       
     [7] "forcats (>= 0.5.1)"       "googledrive (>= 1.0.1)"  
     [9] "googlesheets4 (>= 0.3.0)" "ggplot2 (>= 3.3.3)"      
    [11] "haven (>= 2.3.1)"         "hms (>= 1.0.0)"          
    [13] "httr (>= 1.4.2)"          "jsonlite (>= 1.7.2)"     
    [15] "lubridate (>= 1.7.10)"    "magrittr (>= 2.0.1)"     
    [17] "modelr (>= 0.1.8)"        "pillar (>= 1.6.0)"       
    [19] "purrr (>= 0.3.4)"         "readr (>= 1.4.0)"        
    [21] "readxl (>= 1.3.1)"        "reprex (>= 2.0.0)"       
    [23] "rlang (>= 0.4.10)"        "rstudioapi (>= 0.13)"    
    [25] "rvest (>= 1.0.0)"         "stringr (>= 1.4.0)"      
    [27] "tibble (>= 3.1.0)"        "tidyr (>= 1.1.3)"        
    [29] "xml2 (>= 1.3.2)"         
    
    $Packages$tidyverse$Repository
    [1] "CRAN"
    
    $Packages$tidyverse$Hash
    [1] "fc4c72b6ae9bb283416bd59a3303bbab"
    
    
    $Packages$BiocParallel
    $Packages$BiocParallel$Package
    [1] "BiocParallel"
    
    $Packages$BiocParallel$Version
    [1] "1.28.3"
    
    $Packages$BiocParallel$Source
    [1] "Bioconductor"
    
    $Packages$BiocParallel$Depends
    [1] "methods"      "R (>= 3.5.0)"
    
    $Packages$BiocParallel$Imports
    [1] "stats"         "utils"         "futile.logger" "parallel"     
    [5] "snow"         
    
    $Packages$BiocParallel$LinkingTo
    [1] "BH"
    
    $Packages$BiocParallel$git_url
    [1] "https://git.bioconductor.org/packages/BiocParallel"
    
    $Packages$BiocParallel$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocParallel$git_last_commit
    [1] "2f9d88a"
    
    $Packages$BiocParallel$git_last_commit_date
    [1] "2021-12-07"
    
    $Packages$BiocParallel$Hash
    [1] "460c266560eefcf641c3d7c5169e5bb3"
    
    
    $Packages$bit64
    $Packages$bit64$Package
    [1] "bit64"
    
    $Packages$bit64$Version
    [1] "4.0.5"
    
    $Packages$bit64$Source
    [1] "Repository"
    
    $Packages$bit64$Depends
    [1] "R (>= 3.0.1)"   "bit (>= 4.0.0)" "utils"          "methods"       
    [5] "stats"         
    
    $Packages$bit64$Repository
    [1] "CRAN"
    
    $Packages$bit64$Hash
    [1] "9fe98599ca456d6552421db0d6772d8f"
    
    
    $Packages$foreach
    $Packages$foreach$Package
    [1] "foreach"
    
    $Packages$foreach$Version
    [1] "1.5.2"
    
    $Packages$foreach$Source
    [1] "Repository"
    
    $Packages$foreach$Depends
    [1] "R (>= 2.5.0)"
    
    $Packages$foreach$Imports
    [1] "codetools" "utils"     "iterators"
    
    $Packages$foreach$Repository
    [1] "CRAN"
    
    $Packages$foreach$Hash
    [1] "618609b42c9406731ead03adf5379850"
    
    
    $Packages$robustbase
    $Packages$robustbase$Package
    [1] "robustbase"
    
    $Packages$robustbase$Version
    [1] "0.93-9"
    
    $Packages$robustbase$Source
    [1] "Repository"
    
    $Packages$robustbase$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$robustbase$Imports
    [1] "stats"    "graphics" "utils"    "methods"  "DEoptimR"
    
    $Packages$robustbase$Repository
    [1] "CRAN"
    
    $Packages$robustbase$Hash
    [1] "e7b310bebca9aca51848b0a366a427cf"
    
    
    $Packages$Biostrings
    $Packages$Biostrings$Package
    [1] "Biostrings"
    
    $Packages$Biostrings$Version
    [1] "2.62.0"
    
    $Packages$Biostrings$Source
    [1] "Bioconductor"
    
    $Packages$Biostrings$Depends
    [1] "R (>= 4.0.0)"             "methods"                 
    [3] "BiocGenerics (>= 0.37.0)" "S4Vectors (>= 0.27.12)"  
    [5] "IRanges (>= 2.23.9)"      "XVector (>= 0.29.2)"     
    [7] "GenomeInfoDb"            
    
    $Packages$Biostrings$Imports
    [1] "methods"   "utils"     "grDevices" "graphics"  "stats"     "crayon"   
    
    $Packages$Biostrings$LinkingTo
    [1] "S4Vectors" "IRanges"   "XVector"  
    
    $Packages$Biostrings$git_url
    [1] "https://git.bioconductor.org/packages/Biostrings"
    
    $Packages$Biostrings$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$Biostrings$git_last_commit
    [1] "53ed287"
    
    $Packages$Biostrings$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$Biostrings$Hash
    [1] "a2b71046a5e013f4929b85937392a1f9"
    
    
    $Packages$progressr
    $Packages$progressr$Package
    [1] "progressr"
    
    $Packages$progressr$Version
    [1] "0.10.0"
    
    $Packages$progressr$Source
    [1] "Repository"
    
    $Packages$progressr$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$progressr$Imports
    [1] "digest" "utils" 
    
    $Packages$progressr$Repository
    [1] "CRAN"
    
    $Packages$progressr$Hash
    [1] "7df448f5ae46ab0a1af0fb619349d3fd"
    
    
    $Packages$curatedMetagenomicData
    $Packages$curatedMetagenomicData$Package
    [1] "curatedMetagenomicData"
    
    $Packages$curatedMetagenomicData$Version
    [1] "3.2.3"
    
    $Packages$curatedMetagenomicData$Source
    [1] "Bioconductor"
    
    $Packages$curatedMetagenomicData$Depends
    [1] "R (>= 4.1.0)"             "SummarizedExperiment"    
    [3] "TreeSummarizedExperiment"
    
    $Packages$curatedMetagenomicData$Imports
     [1] "AnnotationHub" "ExperimentHub" "S4Vectors"     "dplyr"        
     [5] "magrittr"      "mia"           "purrr"         "rlang"        
     [9] "stringr"       "tibble"        "tidyr"         "tidyselect"   
    
    $Packages$curatedMetagenomicData$git_url
    [1] "https://git.bioconductor.org/packages/curatedMetagenomicData"
    
    $Packages$curatedMetagenomicData$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$curatedMetagenomicData$git_last_commit
    [1] "7e9fe2a"
    
    $Packages$curatedMetagenomicData$git_last_commit_date
    [1] "2021-12-16"
    
    $Packages$curatedMetagenomicData$Hash
    [1] "5cd0cdfeb02aff0d45a23dc59597bf40"
    
    
    $Packages$rsvd
    $Packages$rsvd$Package
    [1] "rsvd"
    
    $Packages$rsvd$Version
    [1] "1.0.5"
    
    $Packages$rsvd$Source
    [1] "Repository"
    
    $Packages$rsvd$Depends
    [1] "R (>= 4.0.0)"
    
    $Packages$rsvd$Imports
    [1] "Matrix"
    
    $Packages$rsvd$Repository
    [1] "CRAN"
    
    $Packages$rsvd$Hash
    [1] "b462187d887abc519894874486dbd6fd"
    
    
    $Packages$ScaledMatrix
    $Packages$ScaledMatrix$Package
    [1] "ScaledMatrix"
    
    $Packages$ScaledMatrix$Version
    [1] "1.2.0"
    
    $Packages$ScaledMatrix$Source
    [1] "Bioconductor"
    
    $Packages$ScaledMatrix$Imports
    [1] "methods"      "Matrix"       "S4Vectors"    "DelayedArray"
    
    $Packages$ScaledMatrix$git_url
    [1] "https://git.bioconductor.org/packages/ScaledMatrix"
    
    $Packages$ScaledMatrix$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$ScaledMatrix$git_last_commit
    [1] "d0573e1"
    
    $Packages$ScaledMatrix$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$ScaledMatrix$Hash
    [1] "cbfb8a57fa1ee80e29bd822e8597e9a2"
    
    
    $Packages$memoise
    $Packages$memoise$Package
    [1] "memoise"
    
    $Packages$memoise$Version
    [1] "2.0.1"
    
    $Packages$memoise$Source
    [1] "Repository"
    
    $Packages$memoise$Imports
    [1] "rlang (>= 0.4.10)" "cachem"           
    
    $Packages$memoise$Repository
    [1] "CRAN"
    
    $Packages$memoise$Hash
    [1] "e2817ccf4a065c5d9d7f2cfbe7c1d78c"
    
    
    $Packages$RcppTOML
    $Packages$RcppTOML$Package
    [1] "RcppTOML"
    
    $Packages$RcppTOML$Version
    [1] "0.1.7"
    
    $Packages$RcppTOML$Source
    [1] "Repository"
    
    $Packages$RcppTOML$Depends
    [1] "R (>= 3.3.0)"
    
    $Packages$RcppTOML$Imports
    [1] "Rcpp (>= 0.11.5)"
    
    $Packages$RcppTOML$LinkingTo
    [1] "Rcpp"
    
    $Packages$RcppTOML$Repository
    [1] "CRAN"
    
    $Packages$RcppTOML$Hash
    [1] "f8a578aa91321ecec1292f1e2ffadeda"
    
    
    $Packages$evaluate
    $Packages$evaluate$Package
    [1] "evaluate"
    
    $Packages$evaluate$Version
    [1] "0.15"
    
    $Packages$evaluate$Source
    [1] "Repository"
    
    $Packages$evaluate$Depends
    [1] "R (>= 3.0.2)"
    
    $Packages$evaluate$Imports
    [1] "methods"
    
    $Packages$evaluate$Repository
    [1] "CRAN"
    
    $Packages$evaluate$Hash
    [1] "699a7a93d08c962d9f8950b2d7a227f1"
    
    
    $Packages$forcats
    $Packages$forcats$Package
    [1] "forcats"
    
    $Packages$forcats$Version
    [1] "0.5.1"
    
    $Packages$forcats$Source
    [1] "Repository"
    
    $Packages$forcats$Depends
    [1] "R (>= 3.2)"
    
    $Packages$forcats$Imports
    [1] "ellipsis" "magrittr" "rlang"    "tibble"  
    
    $Packages$forcats$Repository
    [1] "CRAN"
    
    $Packages$forcats$Hash
    [1] "81c3244cab67468aac4c60550832655d"
    
    
    $Packages$RApiSerialize
    $Packages$RApiSerialize$Package
    [1] "RApiSerialize"
    
    $Packages$RApiSerialize$Version
    [1] "0.1.0"
    
    $Packages$RApiSerialize$Source
    [1] "Repository"
    
    $Packages$RApiSerialize$Repository
    [1] "CRAN"
    
    $Packages$RApiSerialize$Hash
    [1] "96a7dca805d39693b0cb0dad099a5899"
    
    
    $Packages$tzdb
    $Packages$tzdb$Package
    [1] "tzdb"
    
    $Packages$tzdb$Version
    [1] "0.2.0"
    
    $Packages$tzdb$Source
    [1] "Repository"
    
    $Packages$tzdb$Depends
    [1] "R (>= 3.3)"
    
    $Packages$tzdb$LinkingTo
    [1] "cpp11 (>= 0.4.0)"
    
    $Packages$tzdb$Repository
    [1] "CRAN"
    
    $Packages$tzdb$Hash
    [1] "5e069fb033daf2317bd628d3100b75c5"
    
    
    $Packages$permute
    $Packages$permute$Package
    [1] "permute"
    
    $Packages$permute$Version
    [1] "0.9-7"
    
    $Packages$permute$Source
    [1] "Repository"
    
    $Packages$permute$Depends
    [1] "R (>= 2.14.0)"
    
    $Packages$permute$Imports
    [1] "stats"
    
    $Packages$permute$Repository
    [1] "CRAN"
    
    $Packages$permute$Hash
    [1] "abf0ca85c1c752e0d04f46334e635046"
    
    
    $Packages$callr
    $Packages$callr$Package
    [1] "callr"
    
    $Packages$callr$Version
    [1] "3.7.0"
    
    $Packages$callr$Source
    [1] "Repository"
    
    $Packages$callr$Imports
    [1] "processx (>= 3.5.0)" "R6"                  "utils"              
    
    $Packages$callr$Repository
    [1] "CRAN"
    
    $Packages$callr$Hash
    [1] "461aa75a11ce2400245190ef5d3995df"
    
    
    $Packages$ps
    $Packages$ps$Package
    [1] "ps"
    
    $Packages$ps$Version
    [1] "1.6.0"
    
    $Packages$ps$Source
    [1] "Repository"
    
    $Packages$ps$Depends
    [1] "R (>= 3.1)"
    
    $Packages$ps$Imports
    [1] "utils"
    
    $Packages$ps$Repository
    [1] "CRAN"
    
    $Packages$ps$Hash
    [1] "32620e2001c1dce1af49c49dccbb9420"
    
    
    $Packages$curl
    $Packages$curl$Package
    [1] "curl"
    
    $Packages$curl$Version
    [1] "4.3.2"
    
    $Packages$curl$Source
    [1] "Repository"
    
    $Packages$curl$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$curl$Repository
    [1] "CRAN"
    
    $Packages$curl$Hash
    [1] "022c42d49c28e95d69ca60446dbabf88"
    
    
    $Packages$fansi
    $Packages$fansi$Package
    [1] "fansi"
    
    $Packages$fansi$Version
    [1] "1.0.2"
    
    $Packages$fansi$Source
    [1] "Repository"
    
    $Packages$fansi$Depends
    [1] "R (>= 3.1.0)"
    
    $Packages$fansi$Imports
    [1] "grDevices" "utils"    
    
    $Packages$fansi$Repository
    [1] "CRAN"
    
    $Packages$fansi$Hash
    [1] "f28149c2d7a1342a834b314e95e67260"
    
    
    $Packages$highr
    $Packages$highr$Package
    [1] "highr"
    
    $Packages$highr$Version
    [1] "0.9"
    
    $Packages$highr$Source
    [1] "Repository"
    
    $Packages$highr$Depends
    [1] "R (>= 3.2.3)"
    
    $Packages$highr$Imports
    [1] "xfun (>= 0.18)"
    
    $Packages$highr$Repository
    [1] "CRAN"
    
    $Packages$highr$Hash
    [1] "8eb36c8125038e648e5d111c0d7b2ed4"
    
    
    $Packages$furrr
    $Packages$furrr$Package
    [1] "furrr"
    
    $Packages$furrr$Version
    [1] "0.2.3"
    
    $Packages$furrr$Source
    [1] "Repository"
    
    $Packages$furrr$Depends
    [1] "future (>= 1.19.1)" "R (>= 3.2.0)"      
    
    $Packages$furrr$Imports
    [1] "ellipsis"             "globals (>= 0.13.1)"  "lifecycle (>= 1.0.0)"
    [4] "purrr (>= 0.3.0)"     "rlang (>= 0.3.0)"     "vctrs (>= 0.3.2)"    
    
    $Packages$furrr$Repository
    [1] "CRAN"
    
    $Packages$furrr$Hash
    [1] "2aba4ab06e8707ac054c6683cb6fed56"
    
    
    $Packages$GSEABase
    $Packages$GSEABase$Package
    [1] "GSEABase"
    
    $Packages$GSEABase$Version
    [1] "1.56.0"
    
    $Packages$GSEABase$Source
    [1] "Bioconductor"
    
    $Packages$GSEABase$Depends
    [1] "R (>= 2.6.0)"             "BiocGenerics (>= 0.13.8)"
    [3] "Biobase (>= 2.17.8)"      "annotate (>= 1.45.3)"    
    [5] "methods"                  "graph (>= 1.37.2)"       
    
    $Packages$GSEABase$Imports
    [1] "AnnotationDbi" "XML"          
    
    $Packages$GSEABase$git_url
    [1] "https://git.bioconductor.org/packages/GSEABase"
    
    $Packages$GSEABase$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$GSEABase$git_last_commit
    [1] "ee7c3ca"
    
    $Packages$GSEABase$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$GSEABase$Hash
    [1] "014d9018b21520b4c335927b4866fa67"
    
    
    $Packages$arkdb
    $Packages$arkdb$Package
    [1] "arkdb"
    
    $Packages$arkdb$Version
    [1] "0.0.15"
    
    $Packages$arkdb$Source
    [1] "Repository"
    
    $Packages$arkdb$Depends
    [1] "R (>= 4.0)"
    
    $Packages$arkdb$Imports
    [1] "DBI"   "tools" "utils"
    
    $Packages$arkdb$Repository
    [1] "CRAN"
    
    $Packages$arkdb$Hash
    [1] "f95a9c88800a460544292b775216608f"
    
    
    $Packages$renv
    $Packages$renv$Package
    [1] "renv"
    
    $Packages$renv$Version
    [1] "0.15.4"
    
    $Packages$renv$Source
    [1] "Repository"
    
    $Packages$renv$Imports
    [1] "utils"
    
    $Packages$renv$Repository
    [1] "CRAN"
    
    $Packages$renv$Hash
    [1] "c1078316e1d4f70275fc1ea60c0bc431"
    
    
    $Packages$checkmate
    $Packages$checkmate$Package
    [1] "checkmate"
    
    $Packages$checkmate$Version
    [1] "2.0.0"
    
    $Packages$checkmate$Source
    [1] "Repository"
    
    $Packages$checkmate$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$checkmate$Imports
    [1] "backports (>= 1.1.0)" "utils"               
    
    $Packages$checkmate$Repository
    [1] "CRAN"
    
    $Packages$checkmate$Hash
    [1] "a667800d5f0350371bedeb8b8b950289"
    
    
    $Packages$cachem
    $Packages$cachem$Package
    [1] "cachem"
    
    $Packages$cachem$Version
    [1] "1.0.6"
    
    $Packages$cachem$Source
    [1] "Repository"
    
    $Packages$cachem$Imports
    [1] "rlang"   "fastmap"
    
    $Packages$cachem$Repository
    [1] "CRAN"
    
    $Packages$cachem$Hash
    [1] "648c5b3d71e6a37e3043617489a0a0e9"
    
    
    $Packages$desc
    $Packages$desc$Package
    [1] "desc"
    
    $Packages$desc$Version
    [1] "1.4.1"
    
    $Packages$desc$Source
    [1] "Repository"
    
    $Packages$desc$Depends
    [1] "R (>= 3.4)"
    
    $Packages$desc$Imports
    [1] "cli"       "R6"        "rprojroot" "utils"    
    
    $Packages$desc$Repository
    [1] "CRAN"
    
    $Packages$desc$Hash
    [1] "eebd27ee58fcc58714eedb7aa07d8ad1"
    
    
    $Packages$dials
    $Packages$dials$Package
    [1] "dials"
    
    $Packages$dials$Version
    [1] "0.1.0"
    
    $Packages$dials$Source
    [1] "Repository"
    
    $Packages$dials$Depends
    [1] "R (>= 2.10)" "scales"     
    
    $Packages$dials$Imports
     [1] "DiceDesign"              "dplyr (>= 0.8.5)"       
     [3] "glue"                    "hardhat (>= 0.1.6.9000)"
     [5] "lifecycle"               "purrr"                  
     [7] "rlang"                   "tibble"                 
     [9] "utils"                   "vctrs (>= 0.3.1)"       
    [11] "withr"                  
    
    $Packages$dials$Repository
    [1] "CRAN"
    
    $Packages$dials$Hash
    [1] "38d457070f3093b6b92fc2413a513f74"
    
    
    $Packages$tensorA
    $Packages$tensorA$Package
    [1] "tensorA"
    
    $Packages$tensorA$Version
    [1] "0.36.2"
    
    $Packages$tensorA$Source
    [1] "Repository"
    
    $Packages$tensorA$Depends
    [1] "R (>= 2.2.0)" "stats"       
    
    $Packages$tensorA$Repository
    [1] "CRAN"
    
    $Packages$tensorA$Hash
    [1] "fd792ceac77f96b647fa8d6e1788969a"
    
    
    $Packages$ggplot2
    $Packages$ggplot2$Package
    [1] "ggplot2"
    
    $Packages$ggplot2$Version
    [1] "3.3.5"
    
    $Packages$ggplot2$Source
    [1] "Repository"
    
    $Packages$ggplot2$Depends
    [1] "R (>= 3.3)"
    
    $Packages$ggplot2$Imports
     [1] "digest"            "glue"              "grDevices"        
     [4] "grid"              "gtable (>= 0.1.1)" "isoband"          
     [7] "MASS"              "mgcv"              "rlang (>= 0.4.10)"
    [10] "scales (>= 0.5.0)" "stats"             "tibble"           
    [13] "withr (>= 2.0.0)" 
    
    $Packages$ggplot2$Repository
    [1] "CRAN"
    
    $Packages$ggplot2$Hash
    [1] "d7566c471c7b17e095dd023b9ef155ad"
    
    
    $Packages$rlist
    $Packages$rlist$Package
    [1] "rlist"
    
    $Packages$rlist$Version
    [1] "0.4.6.2"
    
    $Packages$rlist$Source
    [1] "Repository"
    
    $Packages$rlist$Depends
    [1] "R (>= 2.15)"
    
    $Packages$rlist$Imports
    [1] "yaml"       "jsonlite"   "XML"        "data.table"
    
    $Packages$rlist$Repository
    [1] "CRAN"
    
    $Packages$rlist$Hash
    [1] "290c8ea0700d2e7258082d0025386e68"
    
    
    $Packages$ggrepel
    $Packages$ggrepel$Package
    [1] "ggrepel"
    
    $Packages$ggrepel$Version
    [1] "0.9.1"
    
    $Packages$ggrepel$Source
    [1] "Repository"
    
    $Packages$ggrepel$Depends
    [1] "R (>= 3.0.0)"       "ggplot2 (>= 2.2.0)"
    
    $Packages$ggrepel$Imports
    [1] "grid"              "Rcpp"              "rlang (>= 0.3.0)" 
    [4] "scales (>= 0.5.0)"
    
    $Packages$ggrepel$LinkingTo
    [1] "Rcpp"
    
    $Packages$ggrepel$Repository
    [1] "CRAN"
    
    $Packages$ggrepel$Hash
    [1] "08ab869f37e6a7741a64ab9069bcb67d"
    
    
    $Packages$qs
    $Packages$qs$Package
    [1] "qs"
    
    $Packages$qs$Version
    [1] "0.25.3"
    
    $Packages$qs$Source
    [1] "Repository"
    
    $Packages$qs$Depends
    [1] "R (>= 3.0.2)"
    
    $Packages$qs$Imports
    [1] "Rcpp"                   "RApiSerialize"          "stringfish (>= 0.15.1)"
    
    $Packages$qs$LinkingTo
    [1] "Rcpp"          "RApiSerialize" "stringfish"   
    
    $Packages$qs$Repository
    [1] "CRAN"
    
    $Packages$qs$Hash
    [1] "2a135bd3fcd62d0caf73cf1bc9270d8d"
    
    
    $Packages$ade4
    $Packages$ade4$Package
    [1] "ade4"
    
    $Packages$ade4$Version
    [1] "1.7-18"
    
    $Packages$ade4$Source
    [1] "Repository"
    
    $Packages$ade4$Depends
    [1] "R (>= 2.10)"
    
    $Packages$ade4$Imports
    [1] "graphics"  "grDevices" "methods"   "stats"     "utils"     "MASS"     
    [7] "pixmap"    "sp"       
    
    $Packages$ade4$Repository
    [1] "CRAN"
    
    $Packages$ade4$Hash
    [1] "c492e20c0e3e8bbb05ce6f57be525897"
    
    
    $Packages$rprojroot
    $Packages$rprojroot$Package
    [1] "rprojroot"
    
    $Packages$rprojroot$Version
    [1] "2.0.2"
    
    $Packages$rprojroot$Source
    [1] "Repository"
    
    $Packages$rprojroot$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$rprojroot$Repository
    [1] "CRAN"
    
    $Packages$rprojroot$Hash
    [1] "249d8cd1e74a8f6a26194a91b47f21d1"
    
    
    $Packages$sass
    $Packages$sass$Package
    [1] "sass"
    
    $Packages$sass$Version
    [1] "0.4.0"
    
    $Packages$sass$Source
    [1] "Repository"
    
    $Packages$sass$Imports
    [1] "fs"                   "rlang (>= 0.4.10)"    "htmltools (>= 0.5.1)"
    [4] "R6"                   "rappdirs"            
    
    $Packages$sass$Repository
    [1] "CRAN"
    
    $Packages$sass$Hash
    [1] "50cf822feb64bb3977bda0b7091be623"
    
    
    $Packages$magrittr
    $Packages$magrittr$Package
    [1] "magrittr"
    
    $Packages$magrittr$Version
    [1] "2.0.2"
    
    $Packages$magrittr$Source
    [1] "Repository"
    
    $Packages$magrittr$Repository
    [1] "CRAN"
    
    $Packages$magrittr$Hash
    [1] "cdc87ecd81934679d1557633d8e1fe51"
    
    
    $Packages$RCurl
    $Packages$RCurl$Package
    [1] "RCurl"
    
    $Packages$RCurl$Version
    [1] "1.98-1.6"
    
    $Packages$RCurl$Source
    [1] "Repository"
    
    $Packages$RCurl$Depends
    [1] "R (>= 3.4.0)" "methods"     
    
    $Packages$RCurl$Imports
    [1] "bitops"
    
    $Packages$RCurl$Repository
    [1] "CRAN"
    
    $Packages$RCurl$Hash
    [1] "53a4a5e55370f1a1a4afd92356252946"
    
    
    $Packages$tune
    $Packages$tune$Package
    [1] "tune"
    
    $Packages$tune$Version
    [1] "0.1.6"
    
    $Packages$tune$Source
    [1] "Repository"
    
    $Packages$tune$Depends
    [1] "R (>= 2.10)"
    
    $Packages$tune$Imports
     [1] "cli (>= 2.0.0)"       "dials (>= 0.0.9)"     "dplyr (>= 1.0.0)"    
     [4] "foreach"              "generics (>= 0.1.0)"  "ggplot2"             
     [7] "glue"                 "GPfit"                "hardhat (>= 0.1.6)"  
    [10] "lifecycle (>= 1.0.0)" "parsnip (>= 0.1.4)"   "purrr (>= 0.3.2)"    
    [13] "recipes (>= 0.1.15)"  "rlang (>= 0.4.0)"     "rsample (>= 0.0.9)"  
    [16] "tibble (>= 3.1.0)"    "tidyr"                "vctrs (>= 0.3.0)"    
    [19] "withr"                "workflows (>= 0.2.3)" "yardstick (>= 0.0.7)"
    
    $Packages$tune$Repository
    [1] "CRAN"
    
    $Packages$tune$Hash
    [1] "3cabba58cec834867411b81e03acfa7b"
    
    
    $Packages$ape
    $Packages$ape$Package
    [1] "ape"
    
    $Packages$ape$Version
    [1] "5.6-2"
    
    $Packages$ape$Source
    [1] "Repository"
    
    $Packages$ape$Depends
    [1] "R (>= 3.2.0)"
    
    $Packages$ape$Imports
    [1] "nlme"             "lattice"          "graphics"         "methods"         
    [5] "stats"            "tools"            "utils"            "parallel"        
    [9] "Rcpp (>= 0.12.0)"
    
    $Packages$ape$LinkingTo
    [1] "Rcpp"
    
    $Packages$ape$Repository
    [1] "CRAN"
    
    $Packages$ape$Hash
    [1] "894108412a7ec23d5de85cdcce871c8b"
    
    
    $Packages$bayesm
    $Packages$bayesm$Package
    [1] "bayesm"
    
    $Packages$bayesm$Version
    [1] "3.1-4"
    
    $Packages$bayesm$Source
    [1] "Repository"
    
    $Packages$bayesm$Depends
    [1] "R (>= 3.2.0)"
    
    $Packages$bayesm$Imports
    [1] "Rcpp (>= 0.12.0)" "utils"            "stats"            "graphics"        
    [5] "grDevices"       
    
    $Packages$bayesm$LinkingTo
    [1] "Rcpp"          "RcppArmadillo"
    
    $Packages$bayesm$Repository
    [1] "CRAN"
    
    $Packages$bayesm$Hash
    [1] "79a128bc420b4d0a3debb430b1e91f30"
    
    
    $Packages$ggplotify
    $Packages$ggplotify$Package
    [1] "ggplotify"
    
    $Packages$ggplotify$Version
    [1] "0.1.0"
    
    $Packages$ggplotify$Source
    [1] "Repository"
    
    $Packages$ggplotify$Depends
    [1] "R (>= 3.4.0)"
    
    $Packages$ggplotify$Imports
    [1] "ggplot2"      "graphics"     "grDevices"    "grid"         "gridGraphics"
    [6] "yulab.utils" 
    
    $Packages$ggplotify$Repository
    [1] "CRAN"
    
    $Packages$ggplotify$Hash
    [1] "acbcedf783cdb8710168aa0edba42ac0"
    
    
    $Packages$xml2
    $Packages$xml2$Package
    [1] "xml2"
    
    $Packages$xml2$Version
    [1] "1.3.3"
    
    $Packages$xml2$Source
    [1] "Repository"
    
    $Packages$xml2$Depends
    [1] "R (>= 3.1.0)"
    
    $Packages$xml2$Imports
    [1] "methods"
    
    $Packages$xml2$Repository
    [1] "CRAN"
    
    $Packages$xml2$Hash
    [1] "40682ed6a969ea5abfd351eb67833adc"
    
    
    $Packages$taxadb
    $Packages$taxadb$Package
    [1] "taxadb"
    
    $Packages$taxadb$Version
    [1] "0.1.4"
    
    $Packages$taxadb$Source
    [1] "Repository"
    
    $Packages$taxadb$Depends
    [1] "R (>= 4.0)"
    
    $Packages$taxadb$Imports
     [1] "dplyr"            "DBI"              "arkdb (> 0.0.10)" "tibble"          
     [5] "readr"            "rappdirs"         "rlang"            "magrittr"        
     [9] "stringi"          "progress"         "utils"            "dbplyr"          
    [13] "curl"             "jsonlite"         "duckdb"           "contentid"       
    [17] "R.utils"         
    
    $Packages$taxadb$Repository
    [1] "CRAN"
    
    $Packages$taxadb$Hash
    [1] "9586ab951460e1f40350fd68f7d37c4f"
    
    
    $Packages$httr
    $Packages$httr$Package
    [1] "httr"
    
    $Packages$httr$Version
    [1] "1.4.2"
    
    $Packages$httr$Source
    [1] "Repository"
    
    $Packages$httr$Depends
    [1] "R (>= 3.2)"
    
    $Packages$httr$Imports
    [1] "curl (>= 3.0.0)"  "jsonlite"         "mime"             "openssl (>= 0.8)"
    [5] "R6"              
    
    $Packages$httr$Repository
    [1] "CRAN"
    
    $Packages$httr$Hash
    [1] "a525aba14184fec243f9eaec62fbed43"
    
    
    $Packages$assertthat
    $Packages$assertthat$Package
    [1] "assertthat"
    
    $Packages$assertthat$Version
    [1] "0.2.1"
    
    $Packages$assertthat$Source
    [1] "Repository"
    
    $Packages$assertthat$Imports
    [1] "tools"
    
    $Packages$assertthat$Repository
    [1] "CRAN"
    
    $Packages$assertthat$Hash
    [1] "50c838a310445e954bc13f26f26a6ecf"
    
    
    $Packages$rmarkdown
    $Packages$rmarkdown$Package
    [1] "rmarkdown"
    
    $Packages$rmarkdown$Version
    [1] "2.12"
    
    $Packages$rmarkdown$Source
    [1] "Repository"
    
    $Packages$rmarkdown$Depends
    [1] "R (>= 3.0)"
    
    $Packages$rmarkdown$Imports
     [1] "bslib (>= 0.2.5.1)"   "evaluate (>= 0.13)"   "htmltools (>= 0.3.5)"
     [4] "jquerylib"            "jsonlite"             "knitr (>= 1.22)"     
     [7] "methods"              "stringr (>= 1.2.0)"   "tinytex (>= 0.31)"   
    [10] "tools"                "utils"                "xfun (>= 0.21)"      
    [13] "yaml (>= 2.1.19)"    
    
    $Packages$rmarkdown$Repository
    [1] "CRAN"
    
    $Packages$rmarkdown$Hash
    [1] "354da5088ddfdffb73c11cc952885d88"
    
    
    $Packages$globals
    $Packages$globals$Package
    [1] "globals"
    
    $Packages$globals$Version
    [1] "0.14.0"
    
    $Packages$globals$Source
    [1] "Repository"
    
    $Packages$globals$Depends
    [1] "R (>= 3.1.2)"
    
    $Packages$globals$Imports
    [1] "codetools"
    
    $Packages$globals$Repository
    [1] "CRAN"
    
    $Packages$globals$Hash
    [1] "eca8023ed5ca6372479ebb9b3207f5ae"
    
    
    $Packages$R6
    $Packages$R6$Package
    [1] "R6"
    
    $Packages$R6$Version
    [1] "2.5.1"
    
    $Packages$R6$Source
    [1] "Repository"
    
    $Packages$R6$Depends
    [1] "R (>= 3.0)"
    
    $Packages$R6$Repository
    [1] "CRAN"
    
    $Packages$R6$Hash
    [1] "470851b6d5d0ac559e9d01bb352b4021"
    
    
    $Packages$Rhdf5lib
    $Packages$Rhdf5lib$Package
    [1] "Rhdf5lib"
    
    $Packages$Rhdf5lib$Version
    [1] "1.16.0"
    
    $Packages$Rhdf5lib$Source
    [1] "Bioconductor"
    
    $Packages$Rhdf5lib$Depends
    [1] "R (>= 4.0.0)"
    
    $Packages$Rhdf5lib$git_url
    [1] "https://git.bioconductor.org/packages/Rhdf5lib"
    
    $Packages$Rhdf5lib$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$Rhdf5lib$git_last_commit
    [1] "534c497"
    
    $Packages$Rhdf5lib$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$Rhdf5lib$Hash
    [1] "720badd4eb76a481fb97d65e08dc67f0"
    
    
    $Packages$nnet
    $Packages$nnet$Package
    [1] "nnet"
    
    $Packages$nnet$Version
    [1] "7.3-17"
    
    $Packages$nnet$Source
    [1] "Repository"
    
    $Packages$nnet$Depends
    [1] "R (>= 3.0.0)" "stats"        "utils"       
    
    $Packages$nnet$Repository
    [1] "CRAN"
    
    $Packages$nnet$Hash
    [1] "cb1d8d9f300a7e536b89c8a88c53f610"
    
    
    $Packages$RcppHNSW
    $Packages$RcppHNSW$Package
    [1] "RcppHNSW"
    
    $Packages$RcppHNSW$Version
    [1] "0.3.0"
    
    $Packages$RcppHNSW$Source
    [1] "Repository"
    
    $Packages$RcppHNSW$Imports
    [1] "methods"          "Rcpp (>= 0.11.3)"
    
    $Packages$RcppHNSW$LinkingTo
    [1] "Rcpp"
    
    $Packages$RcppHNSW$Repository
    [1] "CRAN"
    
    $Packages$RcppHNSW$Hash
    [1] "ce584c040a9cfa786ffb41f938f5ed19"
    
    
    $Packages$DirichletMultinomial
    $Packages$DirichletMultinomial$Package
    [1] "DirichletMultinomial"
    
    $Packages$DirichletMultinomial$Version
    [1] "1.36.0"
    
    $Packages$DirichletMultinomial$Source
    [1] "Bioconductor"
    
    $Packages$DirichletMultinomial$Depends
    [1] "S4Vectors" "IRanges"  
    
    $Packages$DirichletMultinomial$Imports
    [1] "stats4"       "methods"      "BiocGenerics"
    
    $Packages$DirichletMultinomial$git_url
    [1] "https://git.bioconductor.org/packages/DirichletMultinomial"
    
    $Packages$DirichletMultinomial$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$DirichletMultinomial$git_last_commit
    [1] "926baff"
    
    $Packages$DirichletMultinomial$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$DirichletMultinomial$Hash
    [1] "0365ebbdb9249577b45d1c486ff70720"
    
    
    $Packages$progress
    $Packages$progress$Package
    [1] "progress"
    
    $Packages$progress$Version
    [1] "1.2.2"
    
    $Packages$progress$Source
    [1] "Repository"
    
    $Packages$progress$Imports
    [1] "hms"         "prettyunits" "R6"          "crayon"     
    
    $Packages$progress$Repository
    [1] "CRAN"
    
    $Packages$progress$Hash
    [1] "14dc9f7a3c91ebb14ec5bb9208a07061"
    
    
    $Packages$KEGGREST
    $Packages$KEGGREST$Package
    [1] "KEGGREST"
    
    $Packages$KEGGREST$Version
    [1] "1.34.0"
    
    $Packages$KEGGREST$Source
    [1] "Bioconductor"
    
    $Packages$KEGGREST$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$KEGGREST$Imports
    [1] "methods"    "httr"       "png"        "Biostrings"
    
    $Packages$KEGGREST$git_url
    [1] "https://git.bioconductor.org/packages/KEGGREST"
    
    $Packages$KEGGREST$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$KEGGREST$git_last_commit
    [1] "2056750"
    
    $Packages$KEGGREST$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$KEGGREST$Hash
    [1] "01d7ab3deeee1307b7dfa8eafae2a106"
    
    
    $Packages$treeio
    $Packages$treeio$Package
    [1] "treeio"
    
    $Packages$treeio$Version
    [1] "1.18.1"
    
    $Packages$treeio$Source
    [1] "Bioconductor"
    
    $Packages$treeio$Depends
    [1] "R (>= 3.6.0)"
    
    $Packages$treeio$Imports
    [1] "ape"                 "dplyr"               "jsonlite"           
    [4] "magrittr"            "methods"             "rlang"              
    [7] "tibble"              "tidytree (>= 0.3.0)" "utils"              
    
    $Packages$treeio$git_url
    [1] "https://git.bioconductor.org/packages/treeio"
    
    $Packages$treeio$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$treeio$git_last_commit
    [1] "a06b6b3"
    
    $Packages$treeio$git_last_commit_date
    [1] "2021-11-12"
    
    $Packages$treeio$Hash
    [1] "835f0ab27ea0cfcad448397568fdf4d1"
    
    
    $Packages$stringdist
    $Packages$stringdist$Package
    [1] "stringdist"
    
    $Packages$stringdist$Version
    [1] "0.9.8"
    
    $Packages$stringdist$Source
    [1] "Repository"
    
    $Packages$stringdist$Depends
    [1] "R (>= 2.15.3)"
    
    $Packages$stringdist$Imports
    [1] "parallel"
    
    $Packages$stringdist$Repository
    [1] "CRAN"
    
    $Packages$stringdist$Hash
    [1] "c8b6c79c3b2d7d8475f3ab35cd99737b"
    
    
    $Packages$beachmat
    $Packages$beachmat$Package
    [1] "beachmat"
    
    $Packages$beachmat$Version
    [1] "2.10.0"
    
    $Packages$beachmat$Source
    [1] "Bioconductor"
    
    $Packages$beachmat$Imports
    [1] "methods"                   "DelayedArray (>= 0.15.14)"
    [3] "BiocGenerics"              "Matrix"                   
    [5] "Rcpp"                     
    
    $Packages$beachmat$LinkingTo
    [1] "Rcpp"
    
    $Packages$beachmat$git_url
    [1] "https://git.bioconductor.org/packages/beachmat"
    
    $Packages$beachmat$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$beachmat$git_last_commit
    [1] "b7cc532"
    
    $Packages$beachmat$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$beachmat$Hash
    [1] "f58f6f48893416211092b9600ddeba67"
    
    
    $Packages$BiocVersion
    $Packages$BiocVersion$Package
    [1] "BiocVersion"
    
    $Packages$BiocVersion$Version
    [1] "3.14.0"
    
    $Packages$BiocVersion$Source
    [1] "Bioconductor"
    
    $Packages$BiocVersion$Depends
    [1] "R (>= 4.1.0)"
    
    $Packages$BiocVersion$git_url
    [1] "https://git.bioconductor.org/packages/BiocVersion"
    
    $Packages$BiocVersion$git_branch
    [1] "master"
    
    $Packages$BiocVersion$git_last_commit
    [1] "aa56d93"
    
    $Packages$BiocVersion$git_last_commit_date
    [1] "2021-05-19"
    
    $Packages$BiocVersion$Hash
    [1] "e57437f82a5c13263a9cf442b9796ba3"
    
    
    $Packages$rematch2
    $Packages$rematch2$Package
    [1] "rematch2"
    
    $Packages$rematch2$Version
    [1] "2.1.2"
    
    $Packages$rematch2$Source
    [1] "Repository"
    
    $Packages$rematch2$Imports
    [1] "tibble"
    
    $Packages$rematch2$Repository
    [1] "CRAN"
    
    $Packages$rematch2$Hash
    [1] "76c9e04c712a05848ae7a23d2f170a40"
    
    
    $Packages$HDF5Array
    $Packages$HDF5Array$Package
    [1] "HDF5Array"
    
    $Packages$HDF5Array$Version
    [1] "1.22.1"
    
    $Packages$HDF5Array$Source
    [1] "Bioconductor"
    
    $Packages$HDF5Array$Depends
    [1] "R (>= 3.4)"                "methods"                  
    [3] "DelayedArray (>= 0.15.16)" "rhdf5 (>= 2.31.6)"        
    
    $Packages$HDF5Array$Imports
    [1] "utils"                    "stats"                   
    [3] "tools"                    "Matrix"                  
    [5] "rhdf5filters"             "BiocGenerics (>= 0.31.5)"
    [7] "S4Vectors"                "IRanges"                 
    
    $Packages$HDF5Array$LinkingTo
    [1] "S4Vectors (>= 0.27.13)" "Rhdf5lib"              
    
    $Packages$HDF5Array$git_url
    [1] "https://git.bioconductor.org/packages/HDF5Array"
    
    $Packages$HDF5Array$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$HDF5Array$git_last_commit
    [1] "b3f091f"
    
    $Packages$HDF5Array$git_last_commit_date
    [1] "2021-11-13"
    
    $Packages$HDF5Array$Hash
    [1] "04b03e15910f87fc66dc76cf5f396c7b"
    
    
    $Packages$BiocSingular
    $Packages$BiocSingular$Package
    [1] "BiocSingular"
    
    $Packages$BiocSingular$Version
    [1] "1.10.0"
    
    $Packages$BiocSingular$Source
    [1] "Bioconductor"
    
    $Packages$BiocSingular$Imports
     [1] "BiocGenerics" "S4Vectors"    "Matrix"       "methods"      "utils"       
     [6] "DelayedArray" "BiocParallel" "ScaledMatrix" "irlba"        "rsvd"        
    [11] "Rcpp"         "beachmat"    
    
    $Packages$BiocSingular$LinkingTo
    [1] "Rcpp"     "beachmat"
    
    $Packages$BiocSingular$git_url
    [1] "https://git.bioconductor.org/packages/BiocSingular"
    
    $Packages$BiocSingular$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocSingular$git_last_commit
    [1] "6615ae8"
    
    $Packages$BiocSingular$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$BiocSingular$Hash
    [1] "075c2e9c722f3a98268ead15cff94b63"
    
    
    $Packages$rhdf5
    $Packages$rhdf5$Package
    [1] "rhdf5"
    
    $Packages$rhdf5$Version
    [1] "2.38.0"
    
    $Packages$rhdf5$Source
    [1] "Bioconductor"
    
    $Packages$rhdf5$Depends
    [1] "R (>= 4.0.0)" "methods"     
    
    $Packages$rhdf5$Imports
    [1] "Rhdf5lib (>= 1.13.4)" "rhdf5filters"        
    
    $Packages$rhdf5$LinkingTo
    [1] "Rhdf5lib"
    
    $Packages$rhdf5$git_url
    [1] "https://git.bioconductor.org/packages/rhdf5"
    
    $Packages$rhdf5$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$rhdf5$git_last_commit
    [1] "f6fdfa8"
    
    $Packages$rhdf5$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$rhdf5$Hash
    [1] "55f9159b07a2ec0e0f9e0dda28c8c48f"
    
    
    $Packages$ggfun
    $Packages$ggfun$Package
    [1] "ggfun"
    
    $Packages$ggfun$Version
    [1] "0.0.5"
    
    $Packages$ggfun$Source
    [1] "Repository"
    
    $Packages$ggfun$Imports
    [1] "ggplot2" "grid"    "rlang"   "utils"  
    
    $Packages$ggfun$Repository
    [1] "CRAN"
    
    $Packages$ggfun$Hash
    [1] "be6dbe9ddc9cc77fb9277611eff20a96"
    
    
    $Packages$colorspace
    $Packages$colorspace$Package
    [1] "colorspace"
    
    $Packages$colorspace$Version
    [1] "2.0-3"
    
    $Packages$colorspace$Source
    [1] "Repository"
    
    $Packages$colorspace$Depends
    [1] "R (>= 3.0.0)" "methods"     
    
    $Packages$colorspace$Imports
    [1] "graphics"  "grDevices" "stats"    
    
    $Packages$colorspace$Repository
    [1] "CRAN"
    
    $Packages$colorspace$Hash
    [1] "bb4341986bc8b914f0f0acf2e4a3f2f7"
    
    
    $Packages$generics
    $Packages$generics$Package
    [1] "generics"
    
    $Packages$generics$Version
    [1] "0.1.2"
    
    $Packages$generics$Source
    [1] "Repository"
    
    $Packages$generics$Depends
    [1] "R (>= 3.2)"
    
    $Packages$generics$Imports
    [1] "methods"
    
    $Packages$generics$Repository
    [1] "CRAN"
    
    $Packages$generics$Hash
    [1] "177475892cf4a55865868527654a7741"
    
    
    $Packages$tidytext
    $Packages$tidytext$Package
    [1] "tidytext"
    
    $Packages$tidytext$Version
    [1] "0.3.2"
    
    $Packages$tidytext$Source
    [1] "Repository"
    
    $Packages$tidytext$Depends
    [1] "R (>= 2.10)"
    
    $Packages$tidytext$Imports
     [1] "dplyr"             "generics"          "hunspell"         
     [4] "janeaustenr"       "lifecycle"         "Matrix"           
     [7] "methods"           "purrr (>= 0.1.1)"  "rlang (>= 0.4.10)"
    [10] "stringr"           "tibble"            "tokenizers"       
    [13] "vctrs"            
    
    $Packages$tidytext$Repository
    [1] "CRAN"
    
    $Packages$tidytext$Hash
    [1] "586ffb833b347003334465682edda644"
    
    
    $Packages$base64enc
    $Packages$base64enc$Package
    [1] "base64enc"
    
    $Packages$base64enc$Version
    [1] "0.1-3"
    
    $Packages$base64enc$Source
    [1] "Repository"
    
    $Packages$base64enc$Depends
    [1] "R (>= 2.9.0)"
    
    $Packages$base64enc$Repository
    [1] "CRAN"
    
    $Packages$base64enc$Hash
    [1] "543776ae6848fde2f48ff3816d0628bc"
    
    
    $Packages$compositions
    $Packages$compositions$Package
    [1] "compositions"
    
    $Packages$compositions$Version
    [1] "2.0-4"
    
    $Packages$compositions$Source
    [1] "Repository"
    
    $Packages$compositions$Depends
    [1] "R (>= 3.6)"
    
    $Packages$compositions$Imports
    [1] "methods"    "utils"      "grDevices"  "stats"      "tensorA"   
    [6] "robustbase" "bayesm"     "graphics"   "MASS"      
    
    $Packages$compositions$Repository
    [1] "CRAN"
    
    $Packages$compositions$Hash
    [1] "db2719c9c2910d9d4a3a9f2ce6c8e815"
    
    
    $Packages$pillar
    $Packages$pillar$Package
    [1] "pillar"
    
    $Packages$pillar$Version
    [1] "1.7.0"
    
    $Packages$pillar$Source
    [1] "Repository"
    
    $Packages$pillar$Imports
     [1] "cli (>= 2.3.0)"      "crayon (>= 1.3.4)"   "ellipsis (>= 0.3.2)"
     [4] "fansi"               "glue"                "lifecycle"          
     [7] "rlang (>= 0.3.0)"    "utf8 (>= 1.1.0)"     "utils"              
    [10] "vctrs (>= 0.3.8)"   
    
    $Packages$pillar$Repository
    [1] "CRAN"
    
    $Packages$pillar$Hash
    [1] "51dfc97e1b7069e9f7e6f83f3589c22e"
    
    
    $Packages$sp
    $Packages$sp$Package
    [1] "sp"
    
    $Packages$sp$Version
    [1] "1.4-6"
    
    $Packages$sp$Source
    [1] "Repository"
    
    $Packages$sp$Depends
    [1] "R (>= 3.0.0)" "methods"     
    
    $Packages$sp$Imports
    [1] "utils"     "stats"     "graphics"  "grDevices" "lattice"   "grid"     
    
    $Packages$sp$Repository
    [1] "CRAN"
    
    $Packages$sp$Hash
    [1] "ce8613f4e8c84ef4da9eba65b874ebe9"
    
    
    $Packages$uuid
    $Packages$uuid$Package
    [1] "uuid"
    
    $Packages$uuid$Version
    [1] "1.0-3"
    
    $Packages$uuid$Source
    [1] "Repository"
    
    $Packages$uuid$Depends
    [1] "R (>= 2.9.0)"
    
    $Packages$uuid$Repository
    [1] "CRAN"
    
    $Packages$uuid$Hash
    [1] "2097822ba5e4440b81a0c7525d0315ce"
    
    
    $Packages$MetBrewer
    $Packages$MetBrewer$Package
    [1] "MetBrewer"
    
    $Packages$MetBrewer$Version
    [1] "0.1.0"
    
    $Packages$MetBrewer$Source
    [1] "Repository"
    
    $Packages$MetBrewer$Repository
    [1] "CRAN"
    
    $Packages$MetBrewer$Hash
    [1] "1c71a9a7659f56389404aad6dc8ea8b0"
    
    
    $Packages$ini
    $Packages$ini$Package
    [1] "ini"
    
    $Packages$ini$Version
    [1] "0.3.1"
    
    $Packages$ini$Source
    [1] "Repository"
    
    $Packages$ini$Repository
    [1] "CRAN"
    
    $Packages$ini$Hash
    [1] "6154ec2223172bce8162d4153cda21f7"
    
    
    $Packages$GenomeInfoDbData
    $Packages$GenomeInfoDbData$Package
    [1] "GenomeInfoDbData"
    
    $Packages$GenomeInfoDbData$Version
    [1] "1.2.7"
    
    $Packages$GenomeInfoDbData$Source
    [1] "Bioconductor"
    
    $Packages$GenomeInfoDbData$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$GenomeInfoDbData$Hash
    [1] "89e8144e21da34e26b7c05945cefa3ca"
    
    
    $Packages$plyr
    $Packages$plyr$Package
    [1] "plyr"
    
    $Packages$plyr$Version
    [1] "1.8.6"
    
    $Packages$plyr$Source
    [1] "Repository"
    
    $Packages$plyr$Depends
    [1] "R (>= 3.1.0)"
    
    $Packages$plyr$Imports
    [1] "Rcpp (>= 0.11.0)"
    
    $Packages$plyr$LinkingTo
    [1] "Rcpp"
    
    $Packages$plyr$Repository
    [1] "CRAN"
    
    $Packages$plyr$Hash
    [1] "ec0e5ab4e5f851f6ef32cd1d1984957f"
    
    
    $Packages$targets
    $Packages$targets$Package
    [1] "targets"
    
    $Packages$targets$Version
    [1] "0.11.0"
    
    $Packages$targets$Source
    [1] "Repository"
    
    $Packages$targets$Depends
    [1] "R (>= 3.5.0)"
    
    $Packages$targets$Imports
     [1] "base64url (>= 1.4)"     "callr (>= 3.4.3)"       "cli (>= 2.0.2)"        
     [4] "codetools (>= 0.2.16)"  "data.table (>= 1.12.8)" "digest (>= 0.6.25)"    
     [7] "igraph (>= 1.2.5)"      "knitr (>= 1.34)"        "R6 (>= 2.4.1)"         
    [10] "rlang (>= 0.4.10)"      "stats"                  "tibble (>= 3.0.1)"     
    [13] "tidyselect (>= 1.1.0)"  "tools"                  "utils"                 
    [16] "vctrs (>= 0.2.4)"       "withr (>= 2.4.0)"       "yaml (>= 2.2.1)"       
    
    $Packages$targets$Repository
    [1] "CRAN"
    
    $Packages$targets$Hash
    [1] "bcb5125b621b728d8dc23e39f2f07428"
    
    
    $Packages$gtable
    $Packages$gtable$Package
    [1] "gtable"
    
    $Packages$gtable$Version
    [1] "0.3.0"
    
    $Packages$gtable$Source
    [1] "Repository"
    
    $Packages$gtable$Depends
    [1] "R (>= 3.0)"
    
    $Packages$gtable$Imports
    [1] "grid"
    
    $Packages$gtable$Repository
    [1] "CRAN"
    
    $Packages$gtable$Hash
    [1] "ac5c6baf7822ce8732b343f14c072c4d"
    
    
    $Packages$futile.logger
    $Packages$futile.logger$Package
    [1] "futile.logger"
    
    $Packages$futile.logger$Version
    [1] "1.4.3"
    
    $Packages$futile.logger$Source
    [1] "Repository"
    
    $Packages$futile.logger$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$futile.logger$Imports
    [1] "utils"               "lambda.r (>= 1.1.0)" "futile.options"     
    
    $Packages$futile.logger$Repository
    [1] "CRAN"
    
    $Packages$futile.logger$Hash
    [1] "99f0ace8c05ec7d3683d27083c4f1e7e"
    
    
    $Packages$rvest
    $Packages$rvest$Package
    [1] "rvest"
    
    $Packages$rvest$Version
    [1] "1.0.2"
    
    $Packages$rvest$Source
    [1] "Repository"
    
    $Packages$rvest$Depends
    [1] "R (>= 3.2)"
    
    $Packages$rvest$Imports
    [1] "httr (>= 0.5)"        "lifecycle (>= 1.0.0)" "magrittr"            
    [4] "rlang (>= 0.4.10)"    "selectr"              "tibble"              
    [7] "xml2 (>= 1.3)"       
    
    $Packages$rvest$Repository
    [1] "CRAN"
    
    $Packages$rvest$Hash
    [1] "bb099886deffecd6f9b298b7d4492943"
    
    
    $Packages$stringfish
    $Packages$stringfish$Package
    [1] "stringfish"
    
    $Packages$stringfish$Version
    [1] "0.15.5"
    
    $Packages$stringfish$Source
    [1] "Repository"
    
    $Packages$stringfish$Depends
    [1] "R (>= 3.0.2)"
    
    $Packages$stringfish$Imports
    [1] "Rcpp"         "RcppParallel"
    
    $Packages$stringfish$LinkingTo
    [1] "Rcpp (>= 0.12.18.3)" "RcppParallel"       
    
    $Packages$stringfish$Repository
    [1] "CRAN"
    
    $Packages$stringfish$Hash
    [1] "f5a78e4562f6390f61ddf14848a0a3d9"
    
    
    $Packages$knitr
    $Packages$knitr$Package
    [1] "knitr"
    
    $Packages$knitr$Version
    [1] "1.37"
    
    $Packages$knitr$Source
    [1] "Repository"
    
    $Packages$knitr$Depends
    [1] "R (>= 3.2.3)"
    
    $Packages$knitr$Imports
    [1] "evaluate (>= 0.10)" "highr"              "methods"           
    [4] "stringr (>= 0.6)"   "yaml (>= 2.1.19)"   "xfun (>= 0.27)"    
    [7] "tools"             
    
    $Packages$knitr$Repository
    [1] "CRAN"
    
    $Packages$knitr$Hash
    [1] "a4ec675eb332a33fe7b7fe26f70e1f98"
    
    
    $Packages$RcppArmadillo
    $Packages$RcppArmadillo$Package
    [1] "RcppArmadillo"
    
    $Packages$RcppArmadillo$Version
    [1] "0.10.8.1.0"
    
    $Packages$RcppArmadillo$Source
    [1] "Repository"
    
    $Packages$RcppArmadillo$Depends
    [1] "R (>= 3.3.0)"
    
    $Packages$RcppArmadillo$Imports
    [1] "Rcpp (>= 0.11.0)" "stats"            "utils"            "methods"         
    
    $Packages$RcppArmadillo$LinkingTo
    [1] "Rcpp"
    
    $Packages$RcppArmadillo$Repository
    [1] "CRAN"
    
    $Packages$RcppArmadillo$Hash
    [1] "303dcd8a3fca86087e341ed5cc360abf"
    
    
    $Packages$latticeExtra
    $Packages$latticeExtra$Package
    [1] "latticeExtra"
    
    $Packages$latticeExtra$Version
    [1] "0.6-29"
    
    $Packages$latticeExtra$Source
    [1] "Repository"
    
    $Packages$latticeExtra$Depends
    [1] "R (>= 3.6.0)" "lattice"     
    
    $Packages$latticeExtra$Imports
    [1] "grid"         "stats"        "utils"        "grDevices"    "png"         
    [6] "jpeg"         "RColorBrewer"
    
    $Packages$latticeExtra$Repository
    [1] "CRAN"
    
    $Packages$latticeExtra$Hash
    [1] "590829599d6182cf7461787af34666ee"
    
    
    $Packages$parallelMap
    $Packages$parallelMap$Package
    [1] "parallelMap"
    
    $Packages$parallelMap$Version
    [1] "1.5.1"
    
    $Packages$parallelMap$Source
    [1] "Repository"
    
    $Packages$parallelMap$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$parallelMap$Imports
    [1] "BBmisc (>= 1.8)"      "checkmate (>= 1.8.0)" "parallel"            
    [4] "stats"                "utils"               
    
    $Packages$parallelMap$Repository
    [1] "CRAN"
    
    $Packages$parallelMap$Hash
    [1] "3f3c911183a84966c793987f5d73f3e6"
    
    
    $Packages$TreeSummarizedExperiment
    $Packages$TreeSummarizedExperiment$Package
    [1] "TreeSummarizedExperiment"
    
    $Packages$TreeSummarizedExperiment$Version
    [1] "2.2.0"
    
    $Packages$TreeSummarizedExperiment$Source
    [1] "Bioconductor"
    
    $Packages$TreeSummarizedExperiment$Depends
    [1] "R(>= 3.6.0)"            "SingleCellExperiment"   "S4Vectors (>= 0.23.18)"
    [4] "Biostrings"            
    
    $Packages$TreeSummarizedExperiment$Imports
     [1] "methods"              "BiocGenerics"         "utils"               
     [4] "ape"                  "rlang"                "dplyr"               
     [7] "SummarizedExperiment" "BiocParallel"         "IRanges"             
    [10] "treeio"              
    
    $Packages$TreeSummarizedExperiment$git_url
    [1] "https://git.bioconductor.org/packages/TreeSummarizedExperiment"
    
    $Packages$TreeSummarizedExperiment$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$TreeSummarizedExperiment$git_last_commit
    [1] "3fc3f90"
    
    $Packages$TreeSummarizedExperiment$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$TreeSummarizedExperiment$Hash
    [1] "fef6f5b94c1434c05280a1bfafab7416"
    
    
    $Packages$IRanges
    $Packages$IRanges$Package
    [1] "IRanges"
    
    $Packages$IRanges$Version
    [1] "2.28.0"
    
    $Packages$IRanges$Source
    [1] "Bioconductor"
    
    $Packages$IRanges$Depends
    [1] "R (>= 4.0.0)"             "methods"                 
    [3] "utils"                    "stats"                   
    [5] "BiocGenerics (>= 0.39.2)" "S4Vectors (>= 0.29.19)"  
    
    $Packages$IRanges$Imports
    [1] "stats4"
    
    $Packages$IRanges$LinkingTo
    [1] "S4Vectors"
    
    $Packages$IRanges$git_url
    [1] "https://git.bioconductor.org/packages/IRanges"
    
    $Packages$IRanges$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$IRanges$git_last_commit
    [1] "d85ee90"
    
    $Packages$IRanges$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$IRanges$Hash
    [1] "8299b0f981c924a2b824be548f836e9d"
    
    
    $Packages$fastmap
    $Packages$fastmap$Package
    [1] "fastmap"
    
    $Packages$fastmap$Version
    [1] "1.1.0"
    
    $Packages$fastmap$Source
    [1] "Repository"
    
    $Packages$fastmap$Repository
    [1] "CRAN"
    
    $Packages$fastmap$Hash
    [1] "77bd60a6157420d4ffa93b27cf6a58b8"
    
    
    $Packages$workflowsets
    $Packages$workflowsets$Package
    [1] "workflowsets"
    
    $Packages$workflowsets$Version
    [1] "0.1.0"
    
    $Packages$workflowsets$Source
    [1] "Repository"
    
    $Packages$workflowsets$Depends
    [1] "R (>= 2.10)"
    
    $Packages$workflowsets$Imports
     [1] "cli"                  "dplyr (>= 1.0.0)"     "generics"            
     [4] "ggplot2"              "hardhat (>= 0.1.6)"   "lifecycle (>= 1.0.0)"
     [7] "prettyunits"          "purrr"                "rlang"               
    [10] "rsample (>= 0.0.9)"   "stats"                "tibble (>= 3.1.0)"   
    [13] "tidyr"                "tune (>= 0.1.3)"      "vctrs"               
    [16] "withr"                "workflows (>= 0.2.3)"
    
    $Packages$workflowsets$Repository
    [1] "CRAN"
    
    $Packages$workflowsets$Hash
    [1] "34d9ca8086d92ee78df6747831ce3ec5"
    
    
    $Packages$crosstalk
    $Packages$crosstalk$Package
    [1] "crosstalk"
    
    $Packages$crosstalk$Version
    [1] "1.2.0"
    
    $Packages$crosstalk$Source
    [1] "Repository"
    
    $Packages$crosstalk$Imports
    [1] "htmltools (>= 0.3.6)" "jsonlite"             "lazyeval"            
    [4] "R6"                  
    
    $Packages$crosstalk$Repository
    [1] "CRAN"
    
    $Packages$crosstalk$Hash
    [1] "6aa54f69598c32177e920eb3402e8293"
    
    
    $Packages$doParallel
    $Packages$doParallel$Package
    [1] "doParallel"
    
    $Packages$doParallel$Version
    [1] "1.0.17"
    
    $Packages$doParallel$Source
    [1] "Repository"
    
    $Packages$doParallel$Depends
    [1] "R (>= 2.14.0)"        "foreach (>= 1.2.0)"   "iterators (>= 1.0.0)"
    [4] "parallel"             "utils"               
    
    $Packages$doParallel$Repository
    [1] "CRAN"
    
    $Packages$doParallel$Hash
    [1] "451e5edf411987991ab6a5410c45011f"
    
    
    $Packages$AnnotationDbi
    $Packages$AnnotationDbi$Package
    [1] "AnnotationDbi"
    
    $Packages$AnnotationDbi$Version
    [1] "1.56.2"
    
    $Packages$AnnotationDbi$Source
    [1] "Bioconductor"
    
    $Packages$AnnotationDbi$Depends
    [1] "R (>= 2.7.0)"             "methods"                 
    [3] "utils"                    "stats4"                  
    [5] "BiocGenerics (>= 0.29.2)" "Biobase (>= 1.17.0)"     
    [7] "IRanges"                 
    
    $Packages$AnnotationDbi$Imports
    [1] "DBI"                   "RSQLite"               "S4Vectors (>= 0.9.25)"
    [4] "stats"                 "KEGGREST"             
    
    $Packages$AnnotationDbi$git_url
    [1] "https://git.bioconductor.org/packages/AnnotationDbi"
    
    $Packages$AnnotationDbi$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$AnnotationDbi$git_last_commit
    [1] "13fdc4a"
    
    $Packages$AnnotationDbi$git_last_commit_date
    [1] "2021-11-09"
    
    $Packages$AnnotationDbi$Hash
    [1] "ae553c2779c85946f1452919cc3e71e3"
    
    
    $Packages$broom
    $Packages$broom$Package
    [1] "broom"
    
    $Packages$broom$Version
    [1] "0.7.12"
    
    $Packages$broom$Source
    [1] "Repository"
    
    $Packages$broom$Depends
    [1] "R (>= 3.1)"
    
    $Packages$broom$Imports
     [1] "backports"           "dplyr (>= 1.0.0)"    "ellipsis"           
     [4] "generics (>= 0.0.2)" "glue"                "methods"            
     [7] "purrr"               "rlang"               "stringr"            
    [10] "tibble (>= 3.0.0)"   "tidyr (>= 1.0.0)"    "ggplot2"            
    
    $Packages$broom$Repository
    [1] "CRAN"
    
    $Packages$broom$Hash
    [1] "bfa8a039d77ae8d5413254e572c8abea"
    
    
    $Packages$openssl
    $Packages$openssl$Package
    [1] "openssl"
    
    $Packages$openssl$Version
    [1] "2.0.0"
    
    $Packages$openssl$Source
    [1] "Repository"
    
    $Packages$openssl$Imports
    [1] "askpass"
    
    $Packages$openssl$Repository
    [1] "CRAN"
    
    $Packages$openssl$Hash
    [1] "cf4329aac12c2c44089974559c18e446"
    
    
    $Packages$scales
    $Packages$scales$Package
    [1] "scales"
    
    $Packages$scales$Version
    [1] "1.1.1"
    
    $Packages$scales$Source
    [1] "Repository"
    
    $Packages$scales$Depends
    [1] "R (>= 3.2)"
    
    $Packages$scales$Imports
    [1] "farver (>= 2.0.3)" "labeling"          "lifecycle"        
    [4] "munsell (>= 0.5)"  "R6"                "RColorBrewer"     
    [7] "viridisLite"      
    
    $Packages$scales$Repository
    [1] "CRAN"
    
    $Packages$scales$Hash
    [1] "6f76f71042411426ec8df6c54f34e6dd"
    
    
    $Packages$filelock
    $Packages$filelock$Package
    [1] "filelock"
    
    $Packages$filelock$Version
    [1] "1.0.2"
    
    $Packages$filelock$Source
    [1] "Repository"
    
    $Packages$filelock$Repository
    [1] "CRAN"
    
    $Packages$filelock$Hash
    [1] "38ec653c2613bed60052ba3787bd8a2c"
    
    
    $Packages$backports
    $Packages$backports$Package
    [1] "backports"
    
    $Packages$backports$Version
    [1] "1.4.1"
    
    $Packages$backports$Source
    [1] "Repository"
    
    $Packages$backports$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$backports$Repository
    [1] "CRAN"
    
    $Packages$backports$Hash
    [1] "c39fbec8a30d23e721980b8afb31984c"
    
    
    $Packages$S4Vectors
    $Packages$S4Vectors$Package
    [1] "S4Vectors"
    
    $Packages$S4Vectors$Version
    [1] "0.32.3"
    
    $Packages$S4Vectors$Source
    [1] "Bioconductor"
    
    $Packages$S4Vectors$Depends
    [1] "R (>= 4.0.0)"             "methods"                 
    [3] "utils"                    "stats"                   
    [5] "stats4"                   "BiocGenerics (>= 0.37.0)"
    
    $Packages$S4Vectors$git_url
    [1] "https://git.bioconductor.org/packages/S4Vectors"
    
    $Packages$S4Vectors$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$S4Vectors$git_last_commit
    [1] "ad90e78"
    
    $Packages$S4Vectors$git_last_commit_date
    [1] "2021-11-18"
    
    $Packages$S4Vectors$Hash
    [1] "7c01112d4a7528159c335b486a046101"
    
    
    $Packages$ipred
    $Packages$ipred$Package
    [1] "ipred"
    
    $Packages$ipred$Version
    [1] "0.9-12"
    
    $Packages$ipred$Source
    [1] "Repository"
    
    $Packages$ipred$Depends
    [1] "R (>= 2.10)"
    
    $Packages$ipred$Imports
    [1] "rpart (>= 3.1-8)" "MASS"             "survival"         "nnet"            
    [5] "class"            "prodlim"         
    
    $Packages$ipred$Repository
    [1] "CRAN"
    
    $Packages$ipred$Hash
    [1] "8312ebd8121ad2eca1c76441040bee5d"
    
    
    $Packages$vroom
    $Packages$vroom$Package
    [1] "vroom"
    
    $Packages$vroom$Version
    [1] "1.5.7"
    
    $Packages$vroom$Source
    [1] "Repository"
    
    $Packages$vroom$Depends
    [1] "R (>= 3.1)"
    
    $Packages$vroom$Imports
     [1] "bit64"             "crayon"            "cli"              
     [4] "glue"              "hms"               "lifecycle"        
     [7] "methods"           "rlang (>= 0.4.2)"  "stats"            
    [10] "tibble (>= 2.0.0)" "tzdb (>= 0.1.1)"   "vctrs (>= 0.2.0)" 
    [13] "tidyselect"        "withr"            
    
    $Packages$vroom$LinkingTo
    [1] "progress (>= 1.2.1)" "cpp11 (>= 0.2.0)"    "tzdb (>= 0.1.1)"    
    
    $Packages$vroom$Repository
    [1] "CRAN"
    
    $Packages$vroom$Hash
    [1] "976507b5a105bc3bdf6a5a5f29e0684f"
    
    
    $Packages$decontam
    $Packages$decontam$Package
    [1] "decontam"
    
    $Packages$decontam$Version
    [1] "1.14.0"
    
    $Packages$decontam$Source
    [1] "Bioconductor"
    
    $Packages$decontam$Depends
    [1] "R (>= 3.4.1)"       "methods (>= 3.4.1)"
    
    $Packages$decontam$Imports
    [1] "ggplot2 (>= 2.1.0)"  "reshape2 (>= 1.4.1)" "stats"              
    
    $Packages$decontam$git_url
    [1] "https://git.bioconductor.org/packages/decontam"
    
    $Packages$decontam$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$decontam$git_last_commit
    [1] "b710769"
    
    $Packages$decontam$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$decontam$Hash
    [1] "d3b0675086a9a97173ea8482b15a18b9"
    
    
    $Packages$hms
    $Packages$hms$Package
    [1] "hms"
    
    $Packages$hms$Version
    [1] "1.1.1"
    
    $Packages$hms$Source
    [1] "Repository"
    
    $Packages$hms$Imports
    [1] "ellipsis (>= 0.3.2)" "lifecycle"           "methods"            
    [4] "pkgconfig"           "rlang"               "vctrs (>= 0.3.8)"   
    
    $Packages$hms$Repository
    [1] "CRAN"
    
    $Packages$hms$Hash
    [1] "5b8a2dd0fdbe2ab4f6081e6c7be6dfca"
    
    
    $Packages$Rtsne
    $Packages$Rtsne$Package
    [1] "Rtsne"
    
    $Packages$Rtsne$Version
    [1] "0.15"
    
    $Packages$Rtsne$Source
    [1] "Repository"
    
    $Packages$Rtsne$Imports
    [1] "Rcpp (>= 0.11.0)" "stats"           
    
    $Packages$Rtsne$LinkingTo
    [1] "Rcpp"
    
    $Packages$Rtsne$Repository
    [1] "CRAN"
    
    $Packages$Rtsne$Hash
    [1] "f153432c4ca15b937ccfaa40f167c892"
    
    
    $Packages$dplyr
    $Packages$dplyr$Package
    [1] "dplyr"
    
    $Packages$dplyr$Version
    [1] "1.0.8"
    
    $Packages$dplyr$Source
    [1] "Repository"
    
    $Packages$dplyr$Depends
    [1] "R (>= 3.4.0)"
    
    $Packages$dplyr$Imports
     [1] "generics"              "glue (>= 1.3.2)"       "lifecycle (>= 1.0.1)" 
     [4] "magrittr (>= 1.5)"     "methods"               "R6"                   
     [7] "rlang (>= 1.0.0)"      "tibble (>= 2.1.3)"     "tidyselect (>= 1.1.1)"
    [10] "utils"                 "vctrs (>= 0.3.5)"      "pillar (>= 1.5.1)"    
    
    $Packages$dplyr$Repository
    [1] "CRAN"
    
    $Packages$dplyr$Hash
    [1] "ef47665e64228a17609d6df877bf86f2"
    
    
    $Packages$shiny
    $Packages$shiny$Package
    [1] "shiny"
    
    $Packages$shiny$Version
    [1] "1.7.1"
    
    $Packages$shiny$Source
    [1] "Repository"
    
    $Packages$shiny$Depends
    [1] "R (>= 3.0.2)" "methods"     
    
    $Packages$shiny$Imports
     [1] "utils"                  "grDevices"              "httpuv (>= 1.5.2)"     
     [4] "mime (>= 0.3)"          "jsonlite (>= 0.9.16)"   "xtable"                
     [7] "fontawesome (>= 0.2.1)" "htmltools (>= 0.5.2)"   "R6 (>= 2.0)"           
    [10] "sourcetools"            "later (>= 1.0.0)"       "promises (>= 1.1.0)"   
    [13] "tools"                  "crayon"                 "rlang (>= 0.4.10)"     
    [16] "fastmap (>= 1.1.0)"     "withr"                  "commonmark (>= 1.7)"   
    [19] "glue (>= 1.3.2)"        "bslib (>= 0.3.0)"       "cachem"                
    [22] "ellipsis"               "lifecycle (>= 0.2.0)"  
    
    $Packages$shiny$Repository
    [1] "CRAN"
    
    $Packages$shiny$Hash
    [1] "00344c227c7bd0ab5d78052c5d736c44"
    
    
    $Packages$numDeriv
    $Packages$numDeriv$Package
    [1] "numDeriv"
    
    $Packages$numDeriv$Version
    [1] "2016.8-1.1"
    
    $Packages$numDeriv$Source
    [1] "Repository"
    
    $Packages$numDeriv$Depends
    [1] "R (>= 2.11.1)"
    
    $Packages$numDeriv$Repository
    [1] "CRAN"
    
    $Packages$numDeriv$Hash
    [1] "df58958f293b166e4ab885ebcad90e02"
    
    
    $Packages$DiceDesign
    $Packages$DiceDesign$Package
    [1] "DiceDesign"
    
    $Packages$DiceDesign$Version
    [1] "1.9"
    
    $Packages$DiceDesign$Source
    [1] "Repository"
    
    $Packages$DiceDesign$Depends
    [1] "R (>= 2.10)"
    
    $Packages$DiceDesign$Repository
    [1] "CRAN"
    
    $Packages$DiceDesign$Hash
    [1] "b7b812ae4484d4bbf0a0baac72e8fc01"
    
    
    $Packages$BH
    $Packages$BH$Package
    [1] "BH"
    
    $Packages$BH$Version
    [1] "1.78.0-0"
    
    $Packages$BH$Source
    [1] "Repository"
    
    $Packages$BH$Repository
    [1] "CRAN"
    
    $Packages$BH$Hash
    [1] "4e348572ffcaa2fb1e610e7a941f6f3a"
    
    
    $Packages$lazyeval
    $Packages$lazyeval$Package
    [1] "lazyeval"
    
    $Packages$lazyeval$Version
    [1] "0.2.2"
    
    $Packages$lazyeval$Source
    [1] "Repository"
    
    $Packages$lazyeval$Depends
    [1] "R (>= 3.1.0)"
    
    $Packages$lazyeval$Repository
    [1] "CRAN"
    
    $Packages$lazyeval$Hash
    [1] "d908914ae53b04d4c0c0fd72ecc35370"
    
    
    $Packages$crayon
    $Packages$crayon$Package
    [1] "crayon"
    
    $Packages$crayon$Version
    [1] "1.5.0"
    
    $Packages$crayon$Source
    [1] "Repository"
    
    $Packages$crayon$Imports
    [1] "grDevices" "methods"   "utils"    
    
    $Packages$crayon$Repository
    [1] "CRAN"
    
    $Packages$crayon$Hash
    [1] "741c2e098e98afe3dc26a7b0e5489f4e"
    
    
    $Packages$MASS
    $Packages$MASS$Package
    [1] "MASS"
    
    $Packages$MASS$Version
    [1] "7.3-55"
    
    $Packages$MASS$Source
    [1] "Repository"
    
    $Packages$MASS$Depends
    [1] "R (>= 3.3.0)" "grDevices"    "graphics"     "stats"        "utils"       
    
    $Packages$MASS$Imports
    [1] "methods"
    
    $Packages$MASS$Repository
    [1] "CRAN"
    
    $Packages$MASS$Hash
    [1] "c5232ffb549f6d7a04a152c34ca1353d"
    
    
    $Packages$pROC
    $Packages$pROC$Package
    [1] "pROC"
    
    $Packages$pROC$Version
    [1] "1.18.0"
    
    $Packages$pROC$Source
    [1] "Repository"
    
    $Packages$pROC$Depends
    [1] "R (>= 2.14)"
    
    $Packages$pROC$Imports
    [1] "methods"          "plyr"             "Rcpp (>= 0.11.1)"
    
    $Packages$pROC$LinkingTo
    [1] "Rcpp"
    
    $Packages$pROC$Repository
    [1] "CRAN"
    
    $Packages$pROC$Hash
    [1] "417fd0d40479932c19faf2747817c473"
    
    
    $Packages$naturalsort
    $Packages$naturalsort$Package
    [1] "naturalsort"
    
    $Packages$naturalsort$Version
    [1] "0.1.3"
    
    $Packages$naturalsort$Source
    [1] "Repository"
    
    $Packages$naturalsort$Repository
    [1] "CRAN"
    
    $Packages$naturalsort$Hash
    [1] "2e48674f7d0ab0a3ace4eba2c1145001"
    
    
    $Packages$sparseMatrixStats
    $Packages$sparseMatrixStats$Package
    [1] "sparseMatrixStats"
    
    $Packages$sparseMatrixStats$Version
    [1] "1.6.0"
    
    $Packages$sparseMatrixStats$Source
    [1] "Bioconductor"
    
    $Packages$sparseMatrixStats$Depends
    [1] "MatrixGenerics (>= 1.5.3)"
    
    $Packages$sparseMatrixStats$Imports
    [1] "Rcpp"                    "Matrix"                 
    [3] "matrixStats (>= 0.60.0)" "methods"                
    
    $Packages$sparseMatrixStats$LinkingTo
    [1] "Rcpp"
    
    $Packages$sparseMatrixStats$git_url
    [1] "https://git.bioconductor.org/packages/sparseMatrixStats"
    
    $Packages$sparseMatrixStats$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$sparseMatrixStats$git_last_commit
    [1] "78627a8"
    
    $Packages$sparseMatrixStats$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$sparseMatrixStats$Hash
    [1] "f8dd82d2581115df0ea380c91ed8b9f6"
    
    
    $Packages$viridis
    $Packages$viridis$Package
    [1] "viridis"
    
    $Packages$viridis$Version
    [1] "0.6.2"
    
    $Packages$viridis$Source
    [1] "Repository"
    
    $Packages$viridis$Depends
    [1] "R (>= 2.10)"            "viridisLite (>= 0.4.0)"
    
    $Packages$viridis$Imports
    [1] "stats"              "ggplot2 (>= 1.0.1)" "gridExtra"         
    
    $Packages$viridis$Repository
    [1] "CRAN"
    
    $Packages$viridis$Hash
    [1] "ee96aee95a7a563e5496f8991e9fde4b"
    
    
    $Packages$svglite
    $Packages$svglite$Package
    [1] "svglite"
    
    $Packages$svglite$Version
    [1] "2.1.0"
    
    $Packages$svglite$Source
    [1] "Repository"
    
    $Packages$svglite$Depends
    [1] "R (>= 3.0.0)"
    
    $Packages$svglite$Imports
    [1] "systemfonts (>= 1.0.0)"
    
    $Packages$svglite$LinkingTo
    [1] "cpp11"       "systemfonts"
    
    $Packages$svglite$Repository
    [1] "CRAN"
    
    $Packages$svglite$Hash
    [1] "68dfdf211af6aa4e5f050f064f64d401"
    
    
    $Packages$rpart
    $Packages$rpart$Package
    [1] "rpart"
    
    $Packages$rpart$Version
    [1] "4.1.16"
    
    $Packages$rpart$Source
    [1] "Repository"
    
    $Packages$rpart$Depends
    [1] "R (>= 2.15.0)" "graphics"      "stats"         "grDevices"    
    
    $Packages$rpart$Repository
    [1] "CRAN"
    
    $Packages$rpart$Hash
    [1] "ea3ca1d9473daabb3cd0f1b4f974c1ed"
    
    
    
    attr(,"class")
    [1] "renv_lockfile"
    
    $lockfile
    $R
    $R$Version
    [1] "4.1.3"
    
    $R$Repositories
    $R$Repositories$CRAN
    [1] "https://cran.rstudio.com"
    
    
    
    $Bioconductor
    $Bioconductor$Version
    [1] "3.14"
    
    
    $Python
    $Python$Version
    [1] "3.8.9"
    
    $Python$Type
    [1] "virtualenv"
    
    $Python$Name
    [1] "./renv/python/virtualenvs/renv-python-3.8"
    
    
    $Packages
    $Packages$ANCOMBC
    $Packages$ANCOMBC$Package
    [1] "ANCOMBC"
    
    $Packages$ANCOMBC$Version
    [1] "1.4.0"
    
    $Packages$ANCOMBC$Source
    [1] "Bioconductor"
    
    $Packages$ANCOMBC$git_url
    [1] "https://git.bioconductor.org/packages/ANCOMBC"
    
    $Packages$ANCOMBC$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$ANCOMBC$git_last_commit
    [1] "b9a7fb1"
    
    $Packages$ANCOMBC$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$ANCOMBC$Hash
    [1] "01d19305af91bf6cf081f8a9dfdb18d1"
    
    $Packages$ANCOMBC$Requirements
    $Packages$ANCOMBC$Requirements[[1]]
    [1] "MASS"
    
    $Packages$ANCOMBC$Requirements[[2]]
    [1] "Rdpack"
    
    $Packages$ANCOMBC$Requirements[[3]]
    [1] "microbiome"
    
    $Packages$ANCOMBC$Requirements[[4]]
    [1] "nloptr"
    
    $Packages$ANCOMBC$Requirements[[5]]
    [1] "phyloseq"
    
    
    
    $Packages$AnnotationDbi
    $Packages$AnnotationDbi$Package
    [1] "AnnotationDbi"
    
    $Packages$AnnotationDbi$Version
    [1] "1.56.2"
    
    $Packages$AnnotationDbi$Source
    [1] "Bioconductor"
    
    $Packages$AnnotationDbi$git_url
    [1] "https://git.bioconductor.org/packages/AnnotationDbi"
    
    $Packages$AnnotationDbi$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$AnnotationDbi$git_last_commit
    [1] "13fdc4a"
    
    $Packages$AnnotationDbi$git_last_commit_date
    [1] "2021-11-09"
    
    $Packages$AnnotationDbi$Hash
    [1] "ae553c2779c85946f1452919cc3e71e3"
    
    $Packages$AnnotationDbi$Requirements
    $Packages$AnnotationDbi$Requirements[[1]]
    [1] "Biobase"
    
    $Packages$AnnotationDbi$Requirements[[2]]
    [1] "BiocGenerics"
    
    $Packages$AnnotationDbi$Requirements[[3]]
    [1] "DBI"
    
    $Packages$AnnotationDbi$Requirements[[4]]
    [1] "IRanges"
    
    $Packages$AnnotationDbi$Requirements[[5]]
    [1] "KEGGREST"
    
    $Packages$AnnotationDbi$Requirements[[6]]
    [1] "RSQLite"
    
    $Packages$AnnotationDbi$Requirements[[7]]
    [1] "S4Vectors"
    
    
    
    $Packages$AnnotationHub
    $Packages$AnnotationHub$Package
    [1] "AnnotationHub"
    
    $Packages$AnnotationHub$Version
    [1] "3.2.2"
    
    $Packages$AnnotationHub$Source
    [1] "Bioconductor"
    
    $Packages$AnnotationHub$git_url
    [1] "https://git.bioconductor.org/packages/AnnotationHub"
    
    $Packages$AnnotationHub$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$AnnotationHub$git_last_commit
    [1] "8e761e1"
    
    $Packages$AnnotationHub$git_last_commit_date
    [1] "2022-02-28"
    
    $Packages$AnnotationHub$Hash
    [1] "7630e5d6cf99e154acd3bfe8ae3b5b8d"
    
    $Packages$AnnotationHub$Requirements
    $Packages$AnnotationHub$Requirements[[1]]
    [1] "AnnotationDbi"
    
    $Packages$AnnotationHub$Requirements[[2]]
    [1] "BiocFileCache"
    
    $Packages$AnnotationHub$Requirements[[3]]
    [1] "BiocGenerics"
    
    $Packages$AnnotationHub$Requirements[[4]]
    [1] "BiocManager"
    
    $Packages$AnnotationHub$Requirements[[5]]
    [1] "BiocVersion"
    
    $Packages$AnnotationHub$Requirements[[6]]
    [1] "RSQLite"
    
    $Packages$AnnotationHub$Requirements[[7]]
    [1] "S4Vectors"
    
    $Packages$AnnotationHub$Requirements[[8]]
    [1] "curl"
    
    $Packages$AnnotationHub$Requirements[[9]]
    [1] "dplyr"
    
    $Packages$AnnotationHub$Requirements[[10]]
    [1] "httr"
    
    $Packages$AnnotationHub$Requirements[[11]]
    [1] "interactiveDisplayBase"
    
    $Packages$AnnotationHub$Requirements[[12]]
    [1] "rappdirs"
    
    $Packages$AnnotationHub$Requirements[[13]]
    [1] "yaml"
    
    
    
    $Packages$BBmisc
    $Packages$BBmisc$Package
    [1] "BBmisc"
    
    $Packages$BBmisc$Version
    [1] "1.11"
    
    $Packages$BBmisc$Source
    [1] "Repository"
    
    $Packages$BBmisc$Repository
    [1] "CRAN"
    
    $Packages$BBmisc$Hash
    [1] "b2ff0878d07259998cd84d739a83b3d9"
    
    $Packages$BBmisc$Requirements
    $Packages$BBmisc$Requirements[[1]]
    [1] "checkmate"
    
    
    
    $Packages$BH
    $Packages$BH$Package
    [1] "BH"
    
    $Packages$BH$Version
    [1] "1.78.0-0"
    
    $Packages$BH$Source
    [1] "Repository"
    
    $Packages$BH$Repository
    [1] "CRAN"
    
    $Packages$BH$Hash
    [1] "4e348572ffcaa2fb1e610e7a941f6f3a"
    
    $Packages$BH$Requirements
    list()
    
    
    $Packages$Biobase
    $Packages$Biobase$Package
    [1] "Biobase"
    
    $Packages$Biobase$Version
    [1] "2.54.0"
    
    $Packages$Biobase$Source
    [1] "Bioconductor"
    
    $Packages$Biobase$git_url
    [1] "https://git.bioconductor.org/packages/Biobase"
    
    $Packages$Biobase$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$Biobase$git_last_commit
    [1] "8215d76"
    
    $Packages$Biobase$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$Biobase$Hash
    [1] "767057b0e2897e2a43b59718be888ccd"
    
    $Packages$Biobase$Requirements
    $Packages$Biobase$Requirements[[1]]
    [1] "BiocGenerics"
    
    
    
    $Packages$BiocFileCache
    $Packages$BiocFileCache$Package
    [1] "BiocFileCache"
    
    $Packages$BiocFileCache$Version
    [1] "2.2.1"
    
    $Packages$BiocFileCache$Source
    [1] "Bioconductor"
    
    $Packages$BiocFileCache$git_url
    [1] "https://git.bioconductor.org/packages/BiocFileCache"
    
    $Packages$BiocFileCache$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocFileCache$git_last_commit
    [1] "cc91212"
    
    $Packages$BiocFileCache$git_last_commit_date
    [1] "2022-01-20"
    
    $Packages$BiocFileCache$Hash
    [1] "565eceff974df21115c8db7b4163ad9f"
    
    $Packages$BiocFileCache$Requirements
    $Packages$BiocFileCache$Requirements[[1]]
    [1] "DBI"
    
    $Packages$BiocFileCache$Requirements[[2]]
    [1] "RSQLite"
    
    $Packages$BiocFileCache$Requirements[[3]]
    [1] "curl"
    
    $Packages$BiocFileCache$Requirements[[4]]
    [1] "dbplyr"
    
    $Packages$BiocFileCache$Requirements[[5]]
    [1] "dplyr"
    
    $Packages$BiocFileCache$Requirements[[6]]
    [1] "filelock"
    
    $Packages$BiocFileCache$Requirements[[7]]
    [1] "httr"
    
    $Packages$BiocFileCache$Requirements[[8]]
    [1] "rappdirs"
    
    
    
    $Packages$BiocGenerics
    $Packages$BiocGenerics$Package
    [1] "BiocGenerics"
    
    $Packages$BiocGenerics$Version
    [1] "0.40.0"
    
    $Packages$BiocGenerics$Source
    [1] "Bioconductor"
    
    $Packages$BiocGenerics$git_url
    [1] "https://git.bioconductor.org/packages/BiocGenerics"
    
    $Packages$BiocGenerics$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocGenerics$git_last_commit
    [1] "0bc1e0e"
    
    $Packages$BiocGenerics$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$BiocGenerics$Hash
    [1] "ddc1d29bbe66aaef34b0d17620b61a69"
    
    $Packages$BiocGenerics$Requirements
    list()
    
    
    $Packages$BiocIO
    $Packages$BiocIO$Package
    [1] "BiocIO"
    
    $Packages$BiocIO$Version
    [1] "1.4.0"
    
    $Packages$BiocIO$Source
    [1] "Bioconductor"
    
    $Packages$BiocIO$git_url
    [1] "https://git.bioconductor.org/packages/BiocIO"
    
    $Packages$BiocIO$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocIO$git_last_commit
    [1] "c335932"
    
    $Packages$BiocIO$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$BiocIO$Hash
    [1] "10f50bd6daf34cd11d222befa8829635"
    
    $Packages$BiocIO$Requirements
    $Packages$BiocIO$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$BiocIO$Requirements[[2]]
    [1] "S4Vectors"
    
    
    
    $Packages$BiocManager
    $Packages$BiocManager$Package
    [1] "BiocManager"
    
    $Packages$BiocManager$Version
    [1] "1.30.16"
    
    $Packages$BiocManager$Source
    [1] "Repository"
    
    $Packages$BiocManager$Repository
    [1] "CRAN"
    
    $Packages$BiocManager$Hash
    [1] "2fdca0877debdd4668190832cdee4c31"
    
    $Packages$BiocManager$Requirements
    list()
    
    
    $Packages$BiocNeighbors
    $Packages$BiocNeighbors$Package
    [1] "BiocNeighbors"
    
    $Packages$BiocNeighbors$Version
    [1] "1.12.0"
    
    $Packages$BiocNeighbors$Source
    [1] "Bioconductor"
    
    $Packages$BiocNeighbors$git_url
    [1] "https://git.bioconductor.org/packages/BiocNeighbors"
    
    $Packages$BiocNeighbors$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocNeighbors$git_last_commit
    [1] "3c8a290"
    
    $Packages$BiocNeighbors$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$BiocNeighbors$Hash
    [1] "5fc6686b8fcbe45e684beecb15e3f3a5"
    
    $Packages$BiocNeighbors$Requirements
    $Packages$BiocNeighbors$Requirements[[1]]
    [1] "BiocParallel"
    
    $Packages$BiocNeighbors$Requirements[[2]]
    [1] "Matrix"
    
    $Packages$BiocNeighbors$Requirements[[3]]
    [1] "Rcpp"
    
    $Packages$BiocNeighbors$Requirements[[4]]
    [1] "RcppHNSW"
    
    $Packages$BiocNeighbors$Requirements[[5]]
    [1] "S4Vectors"
    
    
    
    $Packages$BiocParallel
    $Packages$BiocParallel$Package
    [1] "BiocParallel"
    
    $Packages$BiocParallel$Version
    [1] "1.28.3"
    
    $Packages$BiocParallel$Source
    [1] "Bioconductor"
    
    $Packages$BiocParallel$git_url
    [1] "https://git.bioconductor.org/packages/BiocParallel"
    
    $Packages$BiocParallel$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocParallel$git_last_commit
    [1] "2f9d88a"
    
    $Packages$BiocParallel$git_last_commit_date
    [1] "2021-12-07"
    
    $Packages$BiocParallel$Hash
    [1] "460c266560eefcf641c3d7c5169e5bb3"
    
    $Packages$BiocParallel$Requirements
    $Packages$BiocParallel$Requirements[[1]]
    [1] "BH"
    
    $Packages$BiocParallel$Requirements[[2]]
    [1] "futile.logger"
    
    $Packages$BiocParallel$Requirements[[3]]
    [1] "snow"
    
    
    
    $Packages$BiocSet
    $Packages$BiocSet$Package
    [1] "BiocSet"
    
    $Packages$BiocSet$Version
    [1] "1.8.1"
    
    $Packages$BiocSet$Source
    [1] "Bioconductor"
    
    $Packages$BiocSet$git_url
    [1] "https://git.bioconductor.org/packages/BiocSet"
    
    $Packages$BiocSet$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocSet$git_last_commit
    [1] "1f7d340"
    
    $Packages$BiocSet$git_last_commit_date
    [1] "2021-11-02"
    
    $Packages$BiocSet$Hash
    [1] "1f006a9fb0988ff10b02c4378a67b78b"
    
    $Packages$BiocSet$Requirements
    $Packages$BiocSet$Requirements[[1]]
    [1] "AnnotationDbi"
    
    $Packages$BiocSet$Requirements[[2]]
    [1] "BiocIO"
    
    $Packages$BiocSet$Requirements[[3]]
    [1] "KEGGREST"
    
    $Packages$BiocSet$Requirements[[4]]
    [1] "S4Vectors"
    
    $Packages$BiocSet$Requirements[[5]]
    [1] "dplyr"
    
    $Packages$BiocSet$Requirements[[6]]
    [1] "ontologyIndex"
    
    $Packages$BiocSet$Requirements[[7]]
    [1] "plyr"
    
    $Packages$BiocSet$Requirements[[8]]
    [1] "rlang"
    
    $Packages$BiocSet$Requirements[[9]]
    [1] "tibble"
    
    $Packages$BiocSet$Requirements[[10]]
    [1] "tidyr"
    
    
    
    $Packages$BiocSingular
    $Packages$BiocSingular$Package
    [1] "BiocSingular"
    
    $Packages$BiocSingular$Version
    [1] "1.10.0"
    
    $Packages$BiocSingular$Source
    [1] "Bioconductor"
    
    $Packages$BiocSingular$git_url
    [1] "https://git.bioconductor.org/packages/BiocSingular"
    
    $Packages$BiocSingular$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$BiocSingular$git_last_commit
    [1] "6615ae8"
    
    $Packages$BiocSingular$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$BiocSingular$Hash
    [1] "075c2e9c722f3a98268ead15cff94b63"
    
    $Packages$BiocSingular$Requirements
    $Packages$BiocSingular$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$BiocSingular$Requirements[[2]]
    [1] "BiocParallel"
    
    $Packages$BiocSingular$Requirements[[3]]
    [1] "DelayedArray"
    
    $Packages$BiocSingular$Requirements[[4]]
    [1] "Matrix"
    
    $Packages$BiocSingular$Requirements[[5]]
    [1] "Rcpp"
    
    $Packages$BiocSingular$Requirements[[6]]
    [1] "S4Vectors"
    
    $Packages$BiocSingular$Requirements[[7]]
    [1] "ScaledMatrix"
    
    $Packages$BiocSingular$Requirements[[8]]
    [1] "beachmat"
    
    $Packages$BiocSingular$Requirements[[9]]
    [1] "irlba"
    
    $Packages$BiocSingular$Requirements[[10]]
    [1] "rsvd"
    
    
    
    $Packages$BiocVersion
    $Packages$BiocVersion$Package
    [1] "BiocVersion"
    
    $Packages$BiocVersion$Version
    [1] "3.14.0"
    
    $Packages$BiocVersion$Source
    [1] "Bioconductor"
    
    $Packages$BiocVersion$git_url
    [1] "https://git.bioconductor.org/packages/BiocVersion"
    
    $Packages$BiocVersion$git_branch
    [1] "master"
    
    $Packages$BiocVersion$git_last_commit
    [1] "aa56d93"
    
    $Packages$BiocVersion$git_last_commit_date
    [1] "2021-05-19"
    
    $Packages$BiocVersion$Hash
    [1] "e57437f82a5c13263a9cf442b9796ba3"
    
    $Packages$BiocVersion$Requirements
    list()
    
    
    $Packages$Biostrings
    $Packages$Biostrings$Package
    [1] "Biostrings"
    
    $Packages$Biostrings$Version
    [1] "2.62.0"
    
    $Packages$Biostrings$Source
    [1] "Bioconductor"
    
    $Packages$Biostrings$git_url
    [1] "https://git.bioconductor.org/packages/Biostrings"
    
    $Packages$Biostrings$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$Biostrings$git_last_commit
    [1] "53ed287"
    
    $Packages$Biostrings$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$Biostrings$Hash
    [1] "a2b71046a5e013f4929b85937392a1f9"
    
    $Packages$Biostrings$Requirements
    $Packages$Biostrings$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$Biostrings$Requirements[[2]]
    [1] "GenomeInfoDb"
    
    $Packages$Biostrings$Requirements[[3]]
    [1] "IRanges"
    
    $Packages$Biostrings$Requirements[[4]]
    [1] "S4Vectors"
    
    $Packages$Biostrings$Requirements[[5]]
    [1] "XVector"
    
    $Packages$Biostrings$Requirements[[6]]
    [1] "crayon"
    
    
    
    $Packages$DBI
    $Packages$DBI$Package
    [1] "DBI"
    
    $Packages$DBI$Version
    [1] "1.1.2"
    
    $Packages$DBI$Source
    [1] "Repository"
    
    $Packages$DBI$Repository
    [1] "CRAN"
    
    $Packages$DBI$Hash
    [1] "dcd1743af4336156873e3ce3c950b8b9"
    
    $Packages$DBI$Requirements
    list()
    
    
    $Packages$DECIPHER
    $Packages$DECIPHER$Package
    [1] "DECIPHER"
    
    $Packages$DECIPHER$Version
    [1] "2.22.0"
    
    $Packages$DECIPHER$Source
    [1] "Bioconductor"
    
    $Packages$DECIPHER$git_url
    [1] "https://git.bioconductor.org/packages/DECIPHER"
    
    $Packages$DECIPHER$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$DECIPHER$git_last_commit
    [1] "45da5ca"
    
    $Packages$DECIPHER$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$DECIPHER$Hash
    [1] "ae960336cac449fa17ba2b9c5cb1d092"
    
    $Packages$DECIPHER$Requirements
    $Packages$DECIPHER$Requirements[[1]]
    [1] "Biostrings"
    
    $Packages$DECIPHER$Requirements[[2]]
    [1] "DBI"
    
    $Packages$DECIPHER$Requirements[[3]]
    [1] "IRanges"
    
    $Packages$DECIPHER$Requirements[[4]]
    [1] "RSQLite"
    
    $Packages$DECIPHER$Requirements[[5]]
    [1] "S4Vectors"
    
    $Packages$DECIPHER$Requirements[[6]]
    [1] "XVector"
    
    
    
    $Packages$DEoptimR
    $Packages$DEoptimR$Package
    [1] "DEoptimR"
    
    $Packages$DEoptimR$Version
    [1] "1.0-10"
    
    $Packages$DEoptimR$Source
    [1] "Repository"
    
    $Packages$DEoptimR$Repository
    [1] "CRAN"
    
    $Packages$DEoptimR$Hash
    [1] "455af0b32e7e45047cf7f495d462e456"
    
    $Packages$DEoptimR$Requirements
    list()
    
    
    $Packages$DT
    $Packages$DT$Package
    [1] "DT"
    
    $Packages$DT$Version
    [1] "0.21"
    
    $Packages$DT$Source
    [1] "Repository"
    
    $Packages$DT$Repository
    [1] "CRAN"
    
    $Packages$DT$Hash
    [1] "45fa28dbf288cd606e13ca35d3d72437"
    
    $Packages$DT$Requirements
    $Packages$DT$Requirements[[1]]
    [1] "crosstalk"
    
    $Packages$DT$Requirements[[2]]
    [1] "htmltools"
    
    $Packages$DT$Requirements[[3]]
    [1] "htmlwidgets"
    
    $Packages$DT$Requirements[[4]]
    [1] "jquerylib"
    
    $Packages$DT$Requirements[[5]]
    [1] "jsonlite"
    
    $Packages$DT$Requirements[[6]]
    [1] "magrittr"
    
    $Packages$DT$Requirements[[7]]
    [1] "promises"
    
    
    
    $Packages$DelayedArray
    $Packages$DelayedArray$Package
    [1] "DelayedArray"
    
    $Packages$DelayedArray$Version
    [1] "0.20.0"
    
    $Packages$DelayedArray$Source
    [1] "Bioconductor"
    
    $Packages$DelayedArray$git_url
    [1] "https://git.bioconductor.org/packages/DelayedArray"
    
    $Packages$DelayedArray$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$DelayedArray$git_last_commit
    [1] "829b529"
    
    $Packages$DelayedArray$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$DelayedArray$Hash
    [1] "8400bc4ac5cf78a44eefa2395a2ed782"
    
    $Packages$DelayedArray$Requirements
    $Packages$DelayedArray$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$DelayedArray$Requirements[[2]]
    [1] "IRanges"
    
    $Packages$DelayedArray$Requirements[[3]]
    [1] "Matrix"
    
    $Packages$DelayedArray$Requirements[[4]]
    [1] "MatrixGenerics"
    
    $Packages$DelayedArray$Requirements[[5]]
    [1] "S4Vectors"
    
    
    
    $Packages$DelayedMatrixStats
    $Packages$DelayedMatrixStats$Package
    [1] "DelayedMatrixStats"
    
    $Packages$DelayedMatrixStats$Version
    [1] "1.16.0"
    
    $Packages$DelayedMatrixStats$Source
    [1] "Bioconductor"
    
    $Packages$DelayedMatrixStats$git_url
    [1] "https://git.bioconductor.org/packages/DelayedMatrixStats"
    
    $Packages$DelayedMatrixStats$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$DelayedMatrixStats$git_last_commit
    [1] "d44a3d7"
    
    $Packages$DelayedMatrixStats$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$DelayedMatrixStats$Hash
    [1] "9ae660d2827845adeaa2fbdc3e6ad468"
    
    $Packages$DelayedMatrixStats$Requirements
    $Packages$DelayedMatrixStats$Requirements[[1]]
    [1] "DelayedArray"
    
    $Packages$DelayedMatrixStats$Requirements[[2]]
    [1] "IRanges"
    
    $Packages$DelayedMatrixStats$Requirements[[3]]
    [1] "Matrix"
    
    $Packages$DelayedMatrixStats$Requirements[[4]]
    [1] "MatrixGenerics"
    
    $Packages$DelayedMatrixStats$Requirements[[5]]
    [1] "S4Vectors"
    
    $Packages$DelayedMatrixStats$Requirements[[6]]
    [1] "matrixStats"
    
    $Packages$DelayedMatrixStats$Requirements[[7]]
    [1] "sparseMatrixStats"
    
    
    
    $Packages$DiceDesign
    $Packages$DiceDesign$Package
    [1] "DiceDesign"
    
    $Packages$DiceDesign$Version
    [1] "1.9"
    
    $Packages$DiceDesign$Source
    [1] "Repository"
    
    $Packages$DiceDesign$Repository
    [1] "CRAN"
    
    $Packages$DiceDesign$Hash
    [1] "b7b812ae4484d4bbf0a0baac72e8fc01"
    
    $Packages$DiceDesign$Requirements
    list()
    
    
    $Packages$DirichletMultinomial
    $Packages$DirichletMultinomial$Package
    [1] "DirichletMultinomial"
    
    $Packages$DirichletMultinomial$Version
    [1] "1.36.0"
    
    $Packages$DirichletMultinomial$Source
    [1] "Bioconductor"
    
    $Packages$DirichletMultinomial$git_url
    [1] "https://git.bioconductor.org/packages/DirichletMultinomial"
    
    $Packages$DirichletMultinomial$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$DirichletMultinomial$git_last_commit
    [1] "926baff"
    
    $Packages$DirichletMultinomial$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$DirichletMultinomial$Hash
    [1] "0365ebbdb9249577b45d1c486ff70720"
    
    $Packages$DirichletMultinomial$Requirements
    $Packages$DirichletMultinomial$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$DirichletMultinomial$Requirements[[2]]
    [1] "IRanges"
    
    $Packages$DirichletMultinomial$Requirements[[3]]
    [1] "S4Vectors"
    
    
    
    $Packages$ExperimentHub
    $Packages$ExperimentHub$Package
    [1] "ExperimentHub"
    
    $Packages$ExperimentHub$Version
    [1] "2.2.1"
    
    $Packages$ExperimentHub$Source
    [1] "Bioconductor"
    
    $Packages$ExperimentHub$git_url
    [1] "https://git.bioconductor.org/packages/ExperimentHub"
    
    $Packages$ExperimentHub$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$ExperimentHub$git_last_commit
    [1] "4e10686"
    
    $Packages$ExperimentHub$git_last_commit_date
    [1] "2022-01-20"
    
    $Packages$ExperimentHub$Hash
    [1] "7baef74c130e13932918e975a8878f46"
    
    $Packages$ExperimentHub$Requirements
    $Packages$ExperimentHub$Requirements[[1]]
    [1] "AnnotationHub"
    
    $Packages$ExperimentHub$Requirements[[2]]
    [1] "BiocFileCache"
    
    $Packages$ExperimentHub$Requirements[[3]]
    [1] "BiocGenerics"
    
    $Packages$ExperimentHub$Requirements[[4]]
    [1] "BiocManager"
    
    $Packages$ExperimentHub$Requirements[[5]]
    [1] "S4Vectors"
    
    $Packages$ExperimentHub$Requirements[[6]]
    [1] "curl"
    
    $Packages$ExperimentHub$Requirements[[7]]
    [1] "rappdirs"
    
    
    
    $Packages$FNN
    $Packages$FNN$Package
    [1] "FNN"
    
    $Packages$FNN$Version
    [1] "1.1.3"
    
    $Packages$FNN$Source
    [1] "Repository"
    
    $Packages$FNN$Repository
    [1] "CRAN"
    
    $Packages$FNN$Hash
    [1] "b56998fff55e4a4b4860ad6e8c67e0f9"
    
    $Packages$FNN$Requirements
    list()
    
    
    $Packages$GPfit
    $Packages$GPfit$Package
    [1] "GPfit"
    
    $Packages$GPfit$Version
    [1] "1.0-8"
    
    $Packages$GPfit$Source
    [1] "Repository"
    
    $Packages$GPfit$Repository
    [1] "CRAN"
    
    $Packages$GPfit$Hash
    [1] "29a7dccade1fd037c8262c2a239775eb"
    
    $Packages$GPfit$Requirements
    $Packages$GPfit$Requirements[[1]]
    [1] "lattice"
    
    $Packages$GPfit$Requirements[[2]]
    [1] "lhs"
    
    
    
    $Packages$GSEABase
    $Packages$GSEABase$Package
    [1] "GSEABase"
    
    $Packages$GSEABase$Version
    [1] "1.56.0"
    
    $Packages$GSEABase$Source
    [1] "Bioconductor"
    
    $Packages$GSEABase$git_url
    [1] "https://git.bioconductor.org/packages/GSEABase"
    
    $Packages$GSEABase$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$GSEABase$git_last_commit
    [1] "ee7c3ca"
    
    $Packages$GSEABase$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$GSEABase$Hash
    [1] "014d9018b21520b4c335927b4866fa67"
    
    $Packages$GSEABase$Requirements
    $Packages$GSEABase$Requirements[[1]]
    [1] "AnnotationDbi"
    
    $Packages$GSEABase$Requirements[[2]]
    [1] "Biobase"
    
    $Packages$GSEABase$Requirements[[3]]
    [1] "BiocGenerics"
    
    $Packages$GSEABase$Requirements[[4]]
    [1] "XML"
    
    $Packages$GSEABase$Requirements[[5]]
    [1] "annotate"
    
    $Packages$GSEABase$Requirements[[6]]
    [1] "graph"
    
    
    
    $Packages$GSVA
    $Packages$GSVA$Package
    [1] "GSVA"
    
    $Packages$GSVA$Version
    [1] "1.42.0"
    
    $Packages$GSVA$Source
    [1] "Bioconductor"
    
    $Packages$GSVA$git_url
    [1] "https://git.bioconductor.org/packages/GSVA"
    
    $Packages$GSVA$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$GSVA$git_last_commit
    [1] "c99b10b"
    
    $Packages$GSVA$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$GSVA$Hash
    [1] "95684adcece76d9d48c78a933ae01f87"
    
    $Packages$GSVA$Requirements
    $Packages$GSVA$Requirements[[1]]
    [1] "Biobase"
    
    $Packages$GSVA$Requirements[[2]]
    [1] "BiocParallel"
    
    $Packages$GSVA$Requirements[[3]]
    [1] "BiocSingular"
    
    $Packages$GSVA$Requirements[[4]]
    [1] "DelayedArray"
    
    $Packages$GSVA$Requirements[[5]]
    [1] "DelayedMatrixStats"
    
    $Packages$GSVA$Requirements[[6]]
    [1] "GSEABase"
    
    $Packages$GSVA$Requirements[[7]]
    [1] "HDF5Array"
    
    $Packages$GSVA$Requirements[[8]]
    [1] "IRanges"
    
    $Packages$GSVA$Requirements[[9]]
    [1] "Matrix"
    
    $Packages$GSVA$Requirements[[10]]
    [1] "S4Vectors"
    
    $Packages$GSVA$Requirements[[11]]
    [1] "SingleCellExperiment"
    
    $Packages$GSVA$Requirements[[12]]
    [1] "SummarizedExperiment"
    
    $Packages$GSVA$Requirements[[13]]
    [1] "sparseMatrixStats"
    
    
    
    $Packages$GenomeInfoDb
    $Packages$GenomeInfoDb$Package
    [1] "GenomeInfoDb"
    
    $Packages$GenomeInfoDb$Version
    [1] "1.30.1"
    
    $Packages$GenomeInfoDb$Source
    [1] "Bioconductor"
    
    $Packages$GenomeInfoDb$git_url
    [1] "https://git.bioconductor.org/packages/GenomeInfoDb"
    
    $Packages$GenomeInfoDb$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$GenomeInfoDb$git_last_commit
    [1] "bf8b385"
    
    $Packages$GenomeInfoDb$git_last_commit_date
    [1] "2022-01-27"
    
    $Packages$GenomeInfoDb$Hash
    [1] "30dc66e49ea68fd50f5d3fe4b82f0533"
    
    $Packages$GenomeInfoDb$Requirements
    $Packages$GenomeInfoDb$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$GenomeInfoDb$Requirements[[2]]
    [1] "GenomeInfoDbData"
    
    $Packages$GenomeInfoDb$Requirements[[3]]
    [1] "IRanges"
    
    $Packages$GenomeInfoDb$Requirements[[4]]
    [1] "RCurl"
    
    $Packages$GenomeInfoDb$Requirements[[5]]
    [1] "S4Vectors"
    
    
    
    $Packages$GenomeInfoDbData
    $Packages$GenomeInfoDbData$Package
    [1] "GenomeInfoDbData"
    
    $Packages$GenomeInfoDbData$Version
    [1] "1.2.7"
    
    $Packages$GenomeInfoDbData$Source
    [1] "Bioconductor"
    
    $Packages$GenomeInfoDbData$Hash
    [1] "89e8144e21da34e26b7c05945cefa3ca"
    
    $Packages$GenomeInfoDbData$Requirements
    list()
    
    
    $Packages$GenomicAlignments
    $Packages$GenomicAlignments$Package
    [1] "GenomicAlignments"
    
    $Packages$GenomicAlignments$Version
    [1] "1.30.0"
    
    $Packages$GenomicAlignments$Source
    [1] "Bioconductor"
    
    $Packages$GenomicAlignments$git_url
    [1] "https://git.bioconductor.org/packages/GenomicAlignments"
    
    $Packages$GenomicAlignments$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$GenomicAlignments$git_last_commit
    [1] "9046119"
    
    $Packages$GenomicAlignments$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$GenomicAlignments$Hash
    [1] "e1e27b5fac829942cf1359555dfbed9a"
    
    $Packages$GenomicAlignments$Requirements
    $Packages$GenomicAlignments$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$GenomicAlignments$Requirements[[2]]
    [1] "BiocParallel"
    
    $Packages$GenomicAlignments$Requirements[[3]]
    [1] "Biostrings"
    
    $Packages$GenomicAlignments$Requirements[[4]]
    [1] "GenomeInfoDb"
    
    $Packages$GenomicAlignments$Requirements[[5]]
    [1] "GenomicRanges"
    
    $Packages$GenomicAlignments$Requirements[[6]]
    [1] "IRanges"
    
    $Packages$GenomicAlignments$Requirements[[7]]
    [1] "Rsamtools"
    
    $Packages$GenomicAlignments$Requirements[[8]]
    [1] "S4Vectors"
    
    $Packages$GenomicAlignments$Requirements[[9]]
    [1] "SummarizedExperiment"
    
    
    
    $Packages$GenomicRanges
    $Packages$GenomicRanges$Package
    [1] "GenomicRanges"
    
    $Packages$GenomicRanges$Version
    [1] "1.46.1"
    
    $Packages$GenomicRanges$Source
    [1] "Bioconductor"
    
    $Packages$GenomicRanges$git_url
    [1] "https://git.bioconductor.org/packages/GenomicRanges"
    
    $Packages$GenomicRanges$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$GenomicRanges$git_last_commit
    [1] "e422642"
    
    $Packages$GenomicRanges$git_last_commit_date
    [1] "2021-11-16"
    
    $Packages$GenomicRanges$Hash
    [1] "4d9cdd0c25e0c40075a2756aa78f4960"
    
    $Packages$GenomicRanges$Requirements
    $Packages$GenomicRanges$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$GenomicRanges$Requirements[[2]]
    [1] "GenomeInfoDb"
    
    $Packages$GenomicRanges$Requirements[[3]]
    [1] "IRanges"
    
    $Packages$GenomicRanges$Requirements[[4]]
    [1] "S4Vectors"
    
    $Packages$GenomicRanges$Requirements[[5]]
    [1] "XVector"
    
    
    
    $Packages$HDF5Array
    $Packages$HDF5Array$Package
    [1] "HDF5Array"
    
    $Packages$HDF5Array$Version
    [1] "1.22.1"
    
    $Packages$HDF5Array$Source
    [1] "Bioconductor"
    
    $Packages$HDF5Array$git_url
    [1] "https://git.bioconductor.org/packages/HDF5Array"
    
    $Packages$HDF5Array$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$HDF5Array$git_last_commit
    [1] "b3f091f"
    
    $Packages$HDF5Array$git_last_commit_date
    [1] "2021-11-13"
    
    $Packages$HDF5Array$Hash
    [1] "04b03e15910f87fc66dc76cf5f396c7b"
    
    $Packages$HDF5Array$Requirements
    $Packages$HDF5Array$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$HDF5Array$Requirements[[2]]
    [1] "DelayedArray"
    
    $Packages$HDF5Array$Requirements[[3]]
    [1] "IRanges"
    
    $Packages$HDF5Array$Requirements[[4]]
    [1] "Matrix"
    
    $Packages$HDF5Array$Requirements[[5]]
    [1] "Rhdf5lib"
    
    $Packages$HDF5Array$Requirements[[6]]
    [1] "S4Vectors"
    
    $Packages$HDF5Array$Requirements[[7]]
    [1] "rhdf5"
    
    $Packages$HDF5Array$Requirements[[8]]
    [1] "rhdf5filters"
    
    
    
    $Packages$HMP16SData
    $Packages$HMP16SData$Package
    [1] "HMP16SData"
    
    $Packages$HMP16SData$Version
    [1] "1.14.0"
    
    $Packages$HMP16SData$Source
    [1] "Bioconductor"
    
    $Packages$HMP16SData$git_url
    [1] "https://git.bioconductor.org/packages/HMP16SData"
    
    $Packages$HMP16SData$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$HMP16SData$git_last_commit
    [1] "557973d"
    
    $Packages$HMP16SData$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$HMP16SData$Hash
    [1] "45981b9d777b5fead36ab81ae34a1312"
    
    $Packages$HMP16SData$Requirements
    $Packages$HMP16SData$Requirements[[1]]
    [1] "AnnotationHub"
    
    $Packages$HMP16SData$Requirements[[2]]
    [1] "ExperimentHub"
    
    $Packages$HMP16SData$Requirements[[3]]
    [1] "S4Vectors"
    
    $Packages$HMP16SData$Requirements[[4]]
    [1] "SummarizedExperiment"
    
    $Packages$HMP16SData$Requirements[[5]]
    [1] "assertthat"
    
    $Packages$HMP16SData$Requirements[[6]]
    [1] "dplyr"
    
    $Packages$HMP16SData$Requirements[[7]]
    [1] "kableExtra"
    
    $Packages$HMP16SData$Requirements[[8]]
    [1] "knitr"
    
    $Packages$HMP16SData$Requirements[[9]]
    [1] "magrittr"
    
    $Packages$HMP16SData$Requirements[[10]]
    [1] "readr"
    
    $Packages$HMP16SData$Requirements[[11]]
    [1] "stringr"
    
    $Packages$HMP16SData$Requirements[[12]]
    [1] "tibble"
    
    
    
    $Packages$IRanges
    $Packages$IRanges$Package
    [1] "IRanges"
    
    $Packages$IRanges$Version
    [1] "2.28.0"
    
    $Packages$IRanges$Source
    [1] "Bioconductor"
    
    $Packages$IRanges$git_url
    [1] "https://git.bioconductor.org/packages/IRanges"
    
    $Packages$IRanges$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$IRanges$git_last_commit
    [1] "d85ee90"
    
    $Packages$IRanges$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$IRanges$Hash
    [1] "8299b0f981c924a2b824be548f836e9d"
    
    $Packages$IRanges$Requirements
    $Packages$IRanges$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$IRanges$Requirements[[2]]
    [1] "S4Vectors"
    
    
    
    $Packages$KEGGREST
    $Packages$KEGGREST$Package
    [1] "KEGGREST"
    
    $Packages$KEGGREST$Version
    [1] "1.34.0"
    
    $Packages$KEGGREST$Source
    [1] "Bioconductor"
    
    $Packages$KEGGREST$git_url
    [1] "https://git.bioconductor.org/packages/KEGGREST"
    
    $Packages$KEGGREST$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$KEGGREST$git_last_commit
    [1] "2056750"
    
    $Packages$KEGGREST$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$KEGGREST$Hash
    [1] "01d7ab3deeee1307b7dfa8eafae2a106"
    
    $Packages$KEGGREST$Requirements
    $Packages$KEGGREST$Requirements[[1]]
    [1] "Biostrings"
    
    $Packages$KEGGREST$Requirements[[2]]
    [1] "httr"
    
    $Packages$KEGGREST$Requirements[[3]]
    [1] "png"
    
    
    
    $Packages$KernSmooth
    $Packages$KernSmooth$Package
    [1] "KernSmooth"
    
    $Packages$KernSmooth$Version
    [1] "2.23-20"
    
    $Packages$KernSmooth$Source
    [1] "Repository"
    
    $Packages$KernSmooth$Repository
    [1] "CRAN"
    
    $Packages$KernSmooth$Hash
    [1] "8dcfa99b14c296bc9f1fd64d52fd3ce7"
    
    $Packages$KernSmooth$Requirements
    list()
    
    
    $Packages$MASS
    $Packages$MASS$Package
    [1] "MASS"
    
    $Packages$MASS$Version
    [1] "7.3-55"
    
    $Packages$MASS$Source
    [1] "Repository"
    
    $Packages$MASS$Repository
    [1] "CRAN"
    
    $Packages$MASS$Hash
    [1] "c5232ffb549f6d7a04a152c34ca1353d"
    
    $Packages$MASS$Requirements
    list()
    
    
    $Packages$Matrix
    $Packages$Matrix$Package
    [1] "Matrix"
    
    $Packages$Matrix$Version
    [1] "1.4-0"
    
    $Packages$Matrix$Source
    [1] "Repository"
    
    $Packages$Matrix$Repository
    [1] "CRAN"
    
    $Packages$Matrix$Hash
    [1] "130c0caba175739d98f2963c6a407cf6"
    
    $Packages$Matrix$Requirements
    $Packages$Matrix$Requirements[[1]]
    [1] "lattice"
    
    
    
    $Packages$MatrixGenerics
    $Packages$MatrixGenerics$Package
    [1] "MatrixGenerics"
    
    $Packages$MatrixGenerics$Version
    [1] "1.6.0"
    
    $Packages$MatrixGenerics$Source
    [1] "Bioconductor"
    
    $Packages$MatrixGenerics$git_url
    [1] "https://git.bioconductor.org/packages/MatrixGenerics"
    
    $Packages$MatrixGenerics$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$MatrixGenerics$git_last_commit
    [1] "4588a60"
    
    $Packages$MatrixGenerics$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$MatrixGenerics$Hash
    [1] "506b92cb263d9e014a391b41939badcf"
    
    $Packages$MatrixGenerics$Requirements
    $Packages$MatrixGenerics$Requirements[[1]]
    [1] "matrixStats"
    
    
    
    $Packages$MetBrewer
    $Packages$MetBrewer$Package
    [1] "MetBrewer"
    
    $Packages$MetBrewer$Version
    [1] "0.1.0"
    
    $Packages$MetBrewer$Source
    [1] "Repository"
    
    $Packages$MetBrewer$Repository
    [1] "CRAN"
    
    $Packages$MetBrewer$Hash
    [1] "1c71a9a7659f56389404aad6dc8ea8b0"
    
    $Packages$MetBrewer$Requirements
    list()
    
    
    $Packages$MultiAssayExperiment
    $Packages$MultiAssayExperiment$Package
    [1] "MultiAssayExperiment"
    
    $Packages$MultiAssayExperiment$Version
    [1] "1.20.0"
    
    $Packages$MultiAssayExperiment$Source
    [1] "Bioconductor"
    
    $Packages$MultiAssayExperiment$git_url
    [1] "https://git.bioconductor.org/packages/MultiAssayExperiment"
    
    $Packages$MultiAssayExperiment$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$MultiAssayExperiment$git_last_commit
    [1] "c543a7c"
    
    $Packages$MultiAssayExperiment$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$MultiAssayExperiment$Hash
    [1] "e218c8c5a55951f37fc767d4307b395c"
    
    $Packages$MultiAssayExperiment$Requirements
    $Packages$MultiAssayExperiment$Requirements[[1]]
    [1] "Biobase"
    
    $Packages$MultiAssayExperiment$Requirements[[2]]
    [1] "BiocGenerics"
    
    $Packages$MultiAssayExperiment$Requirements[[3]]
    [1] "GenomicRanges"
    
    $Packages$MultiAssayExperiment$Requirements[[4]]
    [1] "IRanges"
    
    $Packages$MultiAssayExperiment$Requirements[[5]]
    [1] "S4Vectors"
    
    $Packages$MultiAssayExperiment$Requirements[[6]]
    [1] "SummarizedExperiment"
    
    $Packages$MultiAssayExperiment$Requirements[[7]]
    [1] "tidyr"
    
    
    
    $Packages$ParamHelpers
    $Packages$ParamHelpers$Package
    [1] "ParamHelpers"
    
    $Packages$ParamHelpers$Version
    [1] "1.14"
    
    $Packages$ParamHelpers$Source
    [1] "Repository"
    
    $Packages$ParamHelpers$Repository
    [1] "CRAN"
    
    $Packages$ParamHelpers$Hash
    [1] "0feaa54a0bcf59d72eb758da748e0904"
    
    $Packages$ParamHelpers$Requirements
    $Packages$ParamHelpers$Requirements[[1]]
    [1] "BBmisc"
    
    $Packages$ParamHelpers$Requirements[[2]]
    [1] "backports"
    
    $Packages$ParamHelpers$Requirements[[3]]
    [1] "checkmate"
    
    $Packages$ParamHelpers$Requirements[[4]]
    [1] "fastmatch"
    
    
    
    $Packages$R.methodsS3
    $Packages$R.methodsS3$Package
    [1] "R.methodsS3"
    
    $Packages$R.methodsS3$Version
    [1] "1.8.1"
    
    $Packages$R.methodsS3$Source
    [1] "Repository"
    
    $Packages$R.methodsS3$Repository
    [1] "CRAN"
    
    $Packages$R.methodsS3$Hash
    [1] "4bf6453323755202d5909697b6f7c109"
    
    $Packages$R.methodsS3$Requirements
    list()
    
    
    $Packages$R.oo
    $Packages$R.oo$Package
    [1] "R.oo"
    
    $Packages$R.oo$Version
    [1] "1.24.0"
    
    $Packages$R.oo$Source
    [1] "Repository"
    
    $Packages$R.oo$Repository
    [1] "CRAN"
    
    $Packages$R.oo$Hash
    [1] "5709328352717e2f0a9c012be8a97554"
    
    $Packages$R.oo$Requirements
    $Packages$R.oo$Requirements[[1]]
    [1] "R.methodsS3"
    
    
    
    $Packages$R.utils
    $Packages$R.utils$Package
    [1] "R.utils"
    
    $Packages$R.utils$Version
    [1] "2.11.0"
    
    $Packages$R.utils$Source
    [1] "Repository"
    
    $Packages$R.utils$Repository
    [1] "CRAN"
    
    $Packages$R.utils$Hash
    [1] "a7ecb8e60815c7a18648e84cd121b23a"
    
    $Packages$R.utils$Requirements
    $Packages$R.utils$Requirements[[1]]
    [1] "R.methodsS3"
    
    $Packages$R.utils$Requirements[[2]]
    [1] "R.oo"
    
    
    
    $Packages$R6
    $Packages$R6$Package
    [1] "R6"
    
    $Packages$R6$Version
    [1] "2.5.1"
    
    $Packages$R6$Source
    [1] "Repository"
    
    $Packages$R6$Repository
    [1] "CRAN"
    
    $Packages$R6$Hash
    [1] "470851b6d5d0ac559e9d01bb352b4021"
    
    $Packages$R6$Requirements
    list()
    
    
    $Packages$RANN
    $Packages$RANN$Package
    [1] "RANN"
    
    $Packages$RANN$Version
    [1] "2.6.1"
    
    $Packages$RANN$Source
    [1] "Repository"
    
    $Packages$RANN$Repository
    [1] "CRAN"
    
    $Packages$RANN$Hash
    [1] "d128ea05a972d3e67c6f39de52c72bd7"
    
    $Packages$RANN$Requirements
    list()
    
    
    $Packages$RApiSerialize
    $Packages$RApiSerialize$Package
    [1] "RApiSerialize"
    
    $Packages$RApiSerialize$Version
    [1] "0.1.0"
    
    $Packages$RApiSerialize$Source
    [1] "Repository"
    
    $Packages$RApiSerialize$Repository
    [1] "CRAN"
    
    $Packages$RApiSerialize$Hash
    [1] "96a7dca805d39693b0cb0dad099a5899"
    
    $Packages$RApiSerialize$Requirements
    list()
    
    
    $Packages$RColorBrewer
    $Packages$RColorBrewer$Package
    [1] "RColorBrewer"
    
    $Packages$RColorBrewer$Version
    [1] "1.1-2"
    
    $Packages$RColorBrewer$Source
    [1] "Repository"
    
    $Packages$RColorBrewer$Repository
    [1] "CRAN"
    
    $Packages$RColorBrewer$Hash
    [1] "e031418365a7f7a766181ab5a41a5716"
    
    $Packages$RColorBrewer$Requirements
    list()
    
    
    $Packages$RCurl
    $Packages$RCurl$Package
    [1] "RCurl"
    
    $Packages$RCurl$Version
    [1] "1.98-1.6"
    
    $Packages$RCurl$Source
    [1] "Repository"
    
    $Packages$RCurl$Repository
    [1] "CRAN"
    
    $Packages$RCurl$Hash
    [1] "53a4a5e55370f1a1a4afd92356252946"
    
    $Packages$RCurl$Requirements
    $Packages$RCurl$Requirements[[1]]
    [1] "bitops"
    
    
    
    $Packages$ROSE
    $Packages$ROSE$Package
    [1] "ROSE"
    
    $Packages$ROSE$Version
    [1] "0.0-4"
    
    $Packages$ROSE$Source
    [1] "Repository"
    
    $Packages$ROSE$Repository
    [1] "CRAN"
    
    $Packages$ROSE$Hash
    [1] "c6e79999925f0db3f2d88b187d014bf8"
    
    $Packages$ROSE$Requirements
    list()
    
    
    $Packages$RSQLite
    $Packages$RSQLite$Package
    [1] "RSQLite"
    
    $Packages$RSQLite$Version
    [1] "2.2.10"
    
    $Packages$RSQLite$Source
    [1] "Repository"
    
    $Packages$RSQLite$Repository
    [1] "CRAN"
    
    $Packages$RSQLite$Hash
    [1] "5d747a4aa7492f22d659be8c343e9045"
    
    $Packages$RSQLite$Requirements
    $Packages$RSQLite$Requirements[[1]]
    [1] "DBI"
    
    $Packages$RSQLite$Requirements[[2]]
    [1] "Rcpp"
    
    $Packages$RSQLite$Requirements[[3]]
    [1] "bit64"
    
    $Packages$RSQLite$Requirements[[4]]
    [1] "blob"
    
    $Packages$RSQLite$Requirements[[5]]
    [1] "memoise"
    
    $Packages$RSQLite$Requirements[[6]]
    [1] "pkgconfig"
    
    $Packages$RSQLite$Requirements[[7]]
    [1] "plogr"
    
    
    
    $Packages$RSpectra
    $Packages$RSpectra$Package
    [1] "RSpectra"
    
    $Packages$RSpectra$Version
    [1] "0.16-0"
    
    $Packages$RSpectra$Source
    [1] "Repository"
    
    $Packages$RSpectra$Repository
    [1] "CRAN"
    
    $Packages$RSpectra$Hash
    [1] "a41329d24d5a98eaed2bd0159adb1b5f"
    
    $Packages$RSpectra$Requirements
    $Packages$RSpectra$Requirements[[1]]
    [1] "Matrix"
    
    $Packages$RSpectra$Requirements[[2]]
    [1] "Rcpp"
    
    $Packages$RSpectra$Requirements[[3]]
    [1] "RcppEigen"
    
    
    
    $Packages$Rcpp
    $Packages$Rcpp$Package
    [1] "Rcpp"
    
    $Packages$Rcpp$Version
    [1] "1.0.8"
    
    $Packages$Rcpp$Source
    [1] "Repository"
    
    $Packages$Rcpp$Repository
    [1] "CRAN"
    
    $Packages$Rcpp$Hash
    [1] "22b546dd7e337f6c0c58a39983a496bc"
    
    $Packages$Rcpp$Requirements
    list()
    
    
    $Packages$RcppArmadillo
    $Packages$RcppArmadillo$Package
    [1] "RcppArmadillo"
    
    $Packages$RcppArmadillo$Version
    [1] "0.10.8.1.0"
    
    $Packages$RcppArmadillo$Source
    [1] "Repository"
    
    $Packages$RcppArmadillo$Repository
    [1] "CRAN"
    
    $Packages$RcppArmadillo$Hash
    [1] "303dcd8a3fca86087e341ed5cc360abf"
    
    $Packages$RcppArmadillo$Requirements
    $Packages$RcppArmadillo$Requirements[[1]]
    [1] "Rcpp"
    
    
    
    $Packages$RcppEigen
    $Packages$RcppEigen$Package
    [1] "RcppEigen"
    
    $Packages$RcppEigen$Version
    [1] "0.3.3.9.1"
    
    $Packages$RcppEigen$Source
    [1] "Repository"
    
    $Packages$RcppEigen$Repository
    [1] "CRAN"
    
    $Packages$RcppEigen$Hash
    [1] "ddfa72a87fdf4c80466a20818be91d00"
    
    $Packages$RcppEigen$Requirements
    $Packages$RcppEigen$Requirements[[1]]
    [1] "Matrix"
    
    $Packages$RcppEigen$Requirements[[2]]
    [1] "Rcpp"
    
    
    
    $Packages$RcppHNSW
    $Packages$RcppHNSW$Package
    [1] "RcppHNSW"
    
    $Packages$RcppHNSW$Version
    [1] "0.3.0"
    
    $Packages$RcppHNSW$Source
    [1] "Repository"
    
    $Packages$RcppHNSW$Repository
    [1] "CRAN"
    
    $Packages$RcppHNSW$Hash
    [1] "ce584c040a9cfa786ffb41f938f5ed19"
    
    $Packages$RcppHNSW$Requirements
    $Packages$RcppHNSW$Requirements[[1]]
    [1] "Rcpp"
    
    
    
    $Packages$RcppParallel
    $Packages$RcppParallel$Package
    [1] "RcppParallel"
    
    $Packages$RcppParallel$Version
    [1] "5.1.5"
    
    $Packages$RcppParallel$Source
    [1] "Repository"
    
    $Packages$RcppParallel$Repository
    [1] "CRAN"
    
    $Packages$RcppParallel$Hash
    [1] "f3e94e34ff656a7c8336ce01207bc2b8"
    
    $Packages$RcppParallel$Requirements
    list()
    
    
    $Packages$RcppTOML
    $Packages$RcppTOML$Package
    [1] "RcppTOML"
    
    $Packages$RcppTOML$Version
    [1] "0.1.7"
    
    $Packages$RcppTOML$Source
    [1] "Repository"
    
    $Packages$RcppTOML$Repository
    [1] "CRAN"
    
    $Packages$RcppTOML$Hash
    [1] "f8a578aa91321ecec1292f1e2ffadeda"
    
    $Packages$RcppTOML$Requirements
    $Packages$RcppTOML$Requirements[[1]]
    [1] "Rcpp"
    
    
    
    $Packages$Rdpack
    $Packages$Rdpack$Package
    [1] "Rdpack"
    
    $Packages$Rdpack$Version
    [1] "2.1.4"
    
    $Packages$Rdpack$Source
    [1] "Repository"
    
    $Packages$Rdpack$Repository
    [1] "CRAN"
    
    $Packages$Rdpack$Hash
    [1] "7befe052afcbb032eacbbbc0babca848"
    
    $Packages$Rdpack$Requirements
    $Packages$Rdpack$Requirements[[1]]
    [1] "rbibutils"
    
    
    
    $Packages$Rhdf5lib
    $Packages$Rhdf5lib$Package
    [1] "Rhdf5lib"
    
    $Packages$Rhdf5lib$Version
    [1] "1.16.0"
    
    $Packages$Rhdf5lib$Source
    [1] "Bioconductor"
    
    $Packages$Rhdf5lib$git_url
    [1] "https://git.bioconductor.org/packages/Rhdf5lib"
    
    $Packages$Rhdf5lib$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$Rhdf5lib$git_last_commit
    [1] "534c497"
    
    $Packages$Rhdf5lib$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$Rhdf5lib$Hash
    [1] "720badd4eb76a481fb97d65e08dc67f0"
    
    $Packages$Rhdf5lib$Requirements
    list()
    
    
    $Packages$Rhtslib
    $Packages$Rhtslib$Package
    [1] "Rhtslib"
    
    $Packages$Rhtslib$Version
    [1] "1.26.0"
    
    $Packages$Rhtslib$Source
    [1] "Bioconductor"
    
    $Packages$Rhtslib$git_url
    [1] "https://git.bioconductor.org/packages/Rhtslib"
    
    $Packages$Rhtslib$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$Rhtslib$git_last_commit
    [1] "f5b20e9"
    
    $Packages$Rhtslib$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$Rhtslib$Hash
    [1] "b88f760ac9993107b94388100f136f43"
    
    $Packages$Rhtslib$Requirements
    $Packages$Rhtslib$Requirements[[1]]
    [1] "zlibbioc"
    
    
    
    $Packages$Rsamtools
    $Packages$Rsamtools$Package
    [1] "Rsamtools"
    
    $Packages$Rsamtools$Version
    [1] "2.10.0"
    
    $Packages$Rsamtools$Source
    [1] "Bioconductor"
    
    $Packages$Rsamtools$git_url
    [1] "https://git.bioconductor.org/packages/Rsamtools"
    
    $Packages$Rsamtools$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$Rsamtools$git_last_commit
    [1] "b19738e"
    
    $Packages$Rsamtools$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$Rsamtools$Hash
    [1] "55e8c5ccea90dfc370e4abdbc5a01cb9"
    
    $Packages$Rsamtools$Requirements
    $Packages$Rsamtools$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$Rsamtools$Requirements[[2]]
    [1] "BiocParallel"
    
    $Packages$Rsamtools$Requirements[[3]]
    [1] "Biostrings"
    
    $Packages$Rsamtools$Requirements[[4]]
    [1] "GenomeInfoDb"
    
    $Packages$Rsamtools$Requirements[[5]]
    [1] "GenomicRanges"
    
    $Packages$Rsamtools$Requirements[[6]]
    [1] "IRanges"
    
    $Packages$Rsamtools$Requirements[[7]]
    [1] "Rhtslib"
    
    $Packages$Rsamtools$Requirements[[8]]
    [1] "S4Vectors"
    
    $Packages$Rsamtools$Requirements[[9]]
    [1] "XVector"
    
    $Packages$Rsamtools$Requirements[[10]]
    [1] "bitops"
    
    $Packages$Rsamtools$Requirements[[11]]
    [1] "zlibbioc"
    
    
    
    $Packages$Rtsne
    $Packages$Rtsne$Package
    [1] "Rtsne"
    
    $Packages$Rtsne$Version
    [1] "0.15"
    
    $Packages$Rtsne$Source
    [1] "Repository"
    
    $Packages$Rtsne$Repository
    [1] "CRAN"
    
    $Packages$Rtsne$Hash
    [1] "f153432c4ca15b937ccfaa40f167c892"
    
    $Packages$Rtsne$Requirements
    $Packages$Rtsne$Requirements[[1]]
    [1] "Rcpp"
    
    
    
    $Packages$S4Vectors
    $Packages$S4Vectors$Package
    [1] "S4Vectors"
    
    $Packages$S4Vectors$Version
    [1] "0.32.3"
    
    $Packages$S4Vectors$Source
    [1] "Bioconductor"
    
    $Packages$S4Vectors$git_url
    [1] "https://git.bioconductor.org/packages/S4Vectors"
    
    $Packages$S4Vectors$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$S4Vectors$git_last_commit
    [1] "ad90e78"
    
    $Packages$S4Vectors$git_last_commit_date
    [1] "2021-11-18"
    
    $Packages$S4Vectors$Hash
    [1] "7c01112d4a7528159c335b486a046101"
    
    $Packages$S4Vectors$Requirements
    $Packages$S4Vectors$Requirements[[1]]
    [1] "BiocGenerics"
    
    
    
    $Packages$SQUAREM
    $Packages$SQUAREM$Package
    [1] "SQUAREM"
    
    $Packages$SQUAREM$Version
    [1] "2021.1"
    
    $Packages$SQUAREM$Source
    [1] "Repository"
    
    $Packages$SQUAREM$Repository
    [1] "CRAN"
    
    $Packages$SQUAREM$Hash
    [1] "0cf10dab0d023d5b46a5a14387556891"
    
    $Packages$SQUAREM$Requirements
    list()
    
    
    $Packages$ScaledMatrix
    $Packages$ScaledMatrix$Package
    [1] "ScaledMatrix"
    
    $Packages$ScaledMatrix$Version
    [1] "1.2.0"
    
    $Packages$ScaledMatrix$Source
    [1] "Bioconductor"
    
    $Packages$ScaledMatrix$git_url
    [1] "https://git.bioconductor.org/packages/ScaledMatrix"
    
    $Packages$ScaledMatrix$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$ScaledMatrix$git_last_commit
    [1] "d0573e1"
    
    $Packages$ScaledMatrix$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$ScaledMatrix$Hash
    [1] "cbfb8a57fa1ee80e29bd822e8597e9a2"
    
    $Packages$ScaledMatrix$Requirements
    $Packages$ScaledMatrix$Requirements[[1]]
    [1] "DelayedArray"
    
    $Packages$ScaledMatrix$Requirements[[2]]
    [1] "Matrix"
    
    $Packages$ScaledMatrix$Requirements[[3]]
    [1] "S4Vectors"
    
    
    
    $Packages$ShortRead
    $Packages$ShortRead$Package
    [1] "ShortRead"
    
    $Packages$ShortRead$Version
    [1] "1.52.0"
    
    $Packages$ShortRead$Source
    [1] "Bioconductor"
    
    $Packages$ShortRead$git_url
    [1] "https://git.bioconductor.org/packages/ShortRead"
    
    $Packages$ShortRead$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$ShortRead$git_last_commit
    [1] "4d7304d"
    
    $Packages$ShortRead$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$ShortRead$Hash
    [1] "cc412b33be7046db49d38794133ddfc4"
    
    $Packages$ShortRead$Requirements
    $Packages$ShortRead$Requirements[[1]]
    [1] "Biobase"
    
    $Packages$ShortRead$Requirements[[2]]
    [1] "BiocGenerics"
    
    $Packages$ShortRead$Requirements[[3]]
    [1] "BiocParallel"
    
    $Packages$ShortRead$Requirements[[4]]
    [1] "Biostrings"
    
    $Packages$ShortRead$Requirements[[5]]
    [1] "GenomeInfoDb"
    
    $Packages$ShortRead$Requirements[[6]]
    [1] "GenomicAlignments"
    
    $Packages$ShortRead$Requirements[[7]]
    [1] "GenomicRanges"
    
    $Packages$ShortRead$Requirements[[8]]
    [1] "IRanges"
    
    $Packages$ShortRead$Requirements[[9]]
    [1] "Rhtslib"
    
    $Packages$ShortRead$Requirements[[10]]
    [1] "Rsamtools"
    
    $Packages$ShortRead$Requirements[[11]]
    [1] "S4Vectors"
    
    $Packages$ShortRead$Requirements[[12]]
    [1] "XVector"
    
    $Packages$ShortRead$Requirements[[13]]
    [1] "hwriter"
    
    $Packages$ShortRead$Requirements[[14]]
    [1] "lattice"
    
    $Packages$ShortRead$Requirements[[15]]
    [1] "latticeExtra"
    
    $Packages$ShortRead$Requirements[[16]]
    [1] "zlibbioc"
    
    
    
    $Packages$SingleCellExperiment
    $Packages$SingleCellExperiment$Package
    [1] "SingleCellExperiment"
    
    $Packages$SingleCellExperiment$Version
    [1] "1.16.0"
    
    $Packages$SingleCellExperiment$Source
    [1] "Bioconductor"
    
    $Packages$SingleCellExperiment$git_url
    [1] "https://git.bioconductor.org/packages/SingleCellExperiment"
    
    $Packages$SingleCellExperiment$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$SingleCellExperiment$git_last_commit
    [1] "bb27609"
    
    $Packages$SingleCellExperiment$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$SingleCellExperiment$Hash
    [1] "27d907d4865fa44f6f17f757352ae7a3"
    
    $Packages$SingleCellExperiment$Requirements
    $Packages$SingleCellExperiment$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$SingleCellExperiment$Requirements[[2]]
    [1] "DelayedArray"
    
    $Packages$SingleCellExperiment$Requirements[[3]]
    [1] "GenomicRanges"
    
    $Packages$SingleCellExperiment$Requirements[[4]]
    [1] "S4Vectors"
    
    $Packages$SingleCellExperiment$Requirements[[5]]
    [1] "SummarizedExperiment"
    
    
    
    $Packages$SnowballC
    $Packages$SnowballC$Package
    [1] "SnowballC"
    
    $Packages$SnowballC$Version
    [1] "0.7.0"
    
    $Packages$SnowballC$Source
    [1] "Repository"
    
    $Packages$SnowballC$Repository
    [1] "CRAN"
    
    $Packages$SnowballC$Hash
    [1] "bc26e07c0d747fd287c370fe355e7b85"
    
    $Packages$SnowballC$Requirements
    list()
    
    
    $Packages$SummarizedExperiment
    $Packages$SummarizedExperiment$Package
    [1] "SummarizedExperiment"
    
    $Packages$SummarizedExperiment$Version
    [1] "1.24.0"
    
    $Packages$SummarizedExperiment$Source
    [1] "Bioconductor"
    
    $Packages$SummarizedExperiment$git_url
    [1] "https://git.bioconductor.org/packages/SummarizedExperiment"
    
    $Packages$SummarizedExperiment$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$SummarizedExperiment$git_last_commit
    [1] "d37f193"
    
    $Packages$SummarizedExperiment$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$SummarizedExperiment$Hash
    [1] "619feb3635b27a198477316a033b1ea8"
    
    $Packages$SummarizedExperiment$Requirements
    $Packages$SummarizedExperiment$Requirements[[1]]
    [1] "Biobase"
    
    $Packages$SummarizedExperiment$Requirements[[2]]
    [1] "BiocGenerics"
    
    $Packages$SummarizedExperiment$Requirements[[3]]
    [1] "DelayedArray"
    
    $Packages$SummarizedExperiment$Requirements[[4]]
    [1] "GenomeInfoDb"
    
    $Packages$SummarizedExperiment$Requirements[[5]]
    [1] "GenomicRanges"
    
    $Packages$SummarizedExperiment$Requirements[[6]]
    [1] "IRanges"
    
    $Packages$SummarizedExperiment$Requirements[[7]]
    [1] "Matrix"
    
    $Packages$SummarizedExperiment$Requirements[[8]]
    [1] "MatrixGenerics"
    
    $Packages$SummarizedExperiment$Requirements[[9]]
    [1] "S4Vectors"
    
    
    
    $Packages$TreeSummarizedExperiment
    $Packages$TreeSummarizedExperiment$Package
    [1] "TreeSummarizedExperiment"
    
    $Packages$TreeSummarizedExperiment$Version
    [1] "2.2.0"
    
    $Packages$TreeSummarizedExperiment$Source
    [1] "Bioconductor"
    
    $Packages$TreeSummarizedExperiment$git_url
    [1] "https://git.bioconductor.org/packages/TreeSummarizedExperiment"
    
    $Packages$TreeSummarizedExperiment$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$TreeSummarizedExperiment$git_last_commit
    [1] "3fc3f90"
    
    $Packages$TreeSummarizedExperiment$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$TreeSummarizedExperiment$Hash
    [1] "fef6f5b94c1434c05280a1bfafab7416"
    
    $Packages$TreeSummarizedExperiment$Requirements
    $Packages$TreeSummarizedExperiment$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$TreeSummarizedExperiment$Requirements[[2]]
    [1] "BiocParallel"
    
    $Packages$TreeSummarizedExperiment$Requirements[[3]]
    [1] "Biostrings"
    
    $Packages$TreeSummarizedExperiment$Requirements[[4]]
    [1] "IRanges"
    
    $Packages$TreeSummarizedExperiment$Requirements[[5]]
    [1] "S4Vectors"
    
    $Packages$TreeSummarizedExperiment$Requirements[[6]]
    [1] "SingleCellExperiment"
    
    $Packages$TreeSummarizedExperiment$Requirements[[7]]
    [1] "SummarizedExperiment"
    
    $Packages$TreeSummarizedExperiment$Requirements[[8]]
    [1] "ape"
    
    $Packages$TreeSummarizedExperiment$Requirements[[9]]
    [1] "dplyr"
    
    $Packages$TreeSummarizedExperiment$Requirements[[10]]
    [1] "rlang"
    
    $Packages$TreeSummarizedExperiment$Requirements[[11]]
    [1] "treeio"
    
    
    
    $Packages$XML
    $Packages$XML$Package
    [1] "XML"
    
    $Packages$XML$Version
    [1] "3.99-0.9"
    
    $Packages$XML$Source
    [1] "Repository"
    
    $Packages$XML$Repository
    [1] "CRAN"
    
    $Packages$XML$Hash
    [1] "6e75b3884423e21eda89b424a98e091d"
    
    $Packages$XML$Requirements
    list()
    
    
    $Packages$XVector
    $Packages$XVector$Package
    [1] "XVector"
    
    $Packages$XVector$Version
    [1] "0.34.0"
    
    $Packages$XVector$Source
    [1] "Bioconductor"
    
    $Packages$XVector$git_url
    [1] "https://git.bioconductor.org/packages/XVector"
    
    $Packages$XVector$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$XVector$git_last_commit
    [1] "06adb25"
    
    $Packages$XVector$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$XVector$Hash
    [1] "9830cb6fc09640591c2f6436467af3c5"
    
    $Packages$XVector$Requirements
    $Packages$XVector$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$XVector$Requirements[[2]]
    [1] "IRanges"
    
    $Packages$XVector$Requirements[[3]]
    [1] "S4Vectors"
    
    $Packages$XVector$Requirements[[4]]
    [1] "zlibbioc"
    
    
    
    $Packages$ade4
    $Packages$ade4$Package
    [1] "ade4"
    
    $Packages$ade4$Version
    [1] "1.7-18"
    
    $Packages$ade4$Source
    [1] "Repository"
    
    $Packages$ade4$Repository
    [1] "CRAN"
    
    $Packages$ade4$Hash
    [1] "c492e20c0e3e8bbb05ce6f57be525897"
    
    $Packages$ade4$Requirements
    $Packages$ade4$Requirements[[1]]
    [1] "MASS"
    
    $Packages$ade4$Requirements[[2]]
    [1] "pixmap"
    
    $Packages$ade4$Requirements[[3]]
    [1] "sp"
    
    
    
    $Packages$annotate
    $Packages$annotate$Package
    [1] "annotate"
    
    $Packages$annotate$Version
    [1] "1.72.0"
    
    $Packages$annotate$Source
    [1] "Bioconductor"
    
    $Packages$annotate$git_url
    [1] "https://git.bioconductor.org/packages/annotate"
    
    $Packages$annotate$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$annotate$git_last_commit
    [1] "67ac76a"
    
    $Packages$annotate$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$annotate$Hash
    [1] "885c32b60eab51236e4afe8921049ff1"
    
    $Packages$annotate$Requirements
    $Packages$annotate$Requirements[[1]]
    [1] "AnnotationDbi"
    
    $Packages$annotate$Requirements[[2]]
    [1] "Biobase"
    
    $Packages$annotate$Requirements[[3]]
    [1] "BiocGenerics"
    
    $Packages$annotate$Requirements[[4]]
    [1] "DBI"
    
    $Packages$annotate$Requirements[[5]]
    [1] "XML"
    
    $Packages$annotate$Requirements[[6]]
    [1] "httr"
    
    $Packages$annotate$Requirements[[7]]
    [1] "xtable"
    
    
    
    $Packages$ape
    $Packages$ape$Package
    [1] "ape"
    
    $Packages$ape$Version
    [1] "5.6-2"
    
    $Packages$ape$Source
    [1] "Repository"
    
    $Packages$ape$Repository
    [1] "CRAN"
    
    $Packages$ape$Hash
    [1] "894108412a7ec23d5de85cdcce871c8b"
    
    $Packages$ape$Requirements
    $Packages$ape$Requirements[[1]]
    [1] "Rcpp"
    
    $Packages$ape$Requirements[[2]]
    [1] "lattice"
    
    $Packages$ape$Requirements[[3]]
    [1] "nlme"
    
    
    
    $Packages$aplot
    $Packages$aplot$Package
    [1] "aplot"
    
    $Packages$aplot$Version
    [1] "0.1.2"
    
    $Packages$aplot$Source
    [1] "Repository"
    
    $Packages$aplot$Repository
    [1] "CRAN"
    
    $Packages$aplot$Hash
    [1] "10f88b6bfa5c7ab877d1ece725aaa66d"
    
    $Packages$aplot$Requirements
    $Packages$aplot$Requirements[[1]]
    [1] "ggfun"
    
    $Packages$aplot$Requirements[[2]]
    [1] "ggplot2"
    
    $Packages$aplot$Requirements[[3]]
    [1] "ggplotify"
    
    $Packages$aplot$Requirements[[4]]
    [1] "magrittr"
    
    $Packages$aplot$Requirements[[5]]
    [1] "patchwork"
    
    $Packages$aplot$Requirements[[6]]
    [1] "yulab.utils"
    
    
    
    $Packages$arkdb
    $Packages$arkdb$Package
    [1] "arkdb"
    
    $Packages$arkdb$Version
    [1] "0.0.15"
    
    $Packages$arkdb$Source
    [1] "Repository"
    
    $Packages$arkdb$Repository
    [1] "CRAN"
    
    $Packages$arkdb$Hash
    [1] "f95a9c88800a460544292b775216608f"
    
    $Packages$arkdb$Requirements
    $Packages$arkdb$Requirements[[1]]
    [1] "DBI"
    
    
    
    $Packages$askpass
    $Packages$askpass$Package
    [1] "askpass"
    
    $Packages$askpass$Version
    [1] "1.1"
    
    $Packages$askpass$Source
    [1] "Repository"
    
    $Packages$askpass$Repository
    [1] "CRAN"
    
    $Packages$askpass$Hash
    [1] "e8a22846fff485f0be3770c2da758713"
    
    $Packages$askpass$Requirements
    $Packages$askpass$Requirements[[1]]
    [1] "sys"
    
    
    
    $Packages$assertthat
    $Packages$assertthat$Package
    [1] "assertthat"
    
    $Packages$assertthat$Version
    [1] "0.2.1"
    
    $Packages$assertthat$Source
    [1] "Repository"
    
    $Packages$assertthat$Repository
    [1] "CRAN"
    
    $Packages$assertthat$Hash
    [1] "50c838a310445e954bc13f26f26a6ecf"
    
    $Packages$assertthat$Requirements
    list()
    
    
    $Packages$backports
    $Packages$backports$Package
    [1] "backports"
    
    $Packages$backports$Version
    [1] "1.4.1"
    
    $Packages$backports$Source
    [1] "Repository"
    
    $Packages$backports$Repository
    [1] "CRAN"
    
    $Packages$backports$Hash
    [1] "c39fbec8a30d23e721980b8afb31984c"
    
    $Packages$backports$Requirements
    list()
    
    
    $Packages$base64enc
    $Packages$base64enc$Package
    [1] "base64enc"
    
    $Packages$base64enc$Version
    [1] "0.1-3"
    
    $Packages$base64enc$Source
    [1] "Repository"
    
    $Packages$base64enc$Repository
    [1] "CRAN"
    
    $Packages$base64enc$Hash
    [1] "543776ae6848fde2f48ff3816d0628bc"
    
    $Packages$base64enc$Requirements
    list()
    
    
    $Packages$base64url
    $Packages$base64url$Package
    [1] "base64url"
    
    $Packages$base64url$Version
    [1] "1.4"
    
    $Packages$base64url$Source
    [1] "Repository"
    
    $Packages$base64url$Repository
    [1] "CRAN"
    
    $Packages$base64url$Hash
    [1] "0c54cf3a08cc0e550fbd64ad33166143"
    
    $Packages$base64url$Requirements
    $Packages$base64url$Requirements[[1]]
    [1] "backports"
    
    
    
    $Packages$bayesm
    $Packages$bayesm$Package
    [1] "bayesm"
    
    $Packages$bayesm$Version
    [1] "3.1-4"
    
    $Packages$bayesm$Source
    [1] "Repository"
    
    $Packages$bayesm$Repository
    [1] "CRAN"
    
    $Packages$bayesm$Hash
    [1] "79a128bc420b4d0a3debb430b1e91f30"
    
    $Packages$bayesm$Requirements
    $Packages$bayesm$Requirements[[1]]
    [1] "Rcpp"
    
    $Packages$bayesm$Requirements[[2]]
    [1] "RcppArmadillo"
    
    
    
    $Packages$beachmat
    $Packages$beachmat$Package
    [1] "beachmat"
    
    $Packages$beachmat$Version
    [1] "2.10.0"
    
    $Packages$beachmat$Source
    [1] "Bioconductor"
    
    $Packages$beachmat$git_url
    [1] "https://git.bioconductor.org/packages/beachmat"
    
    $Packages$beachmat$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$beachmat$git_last_commit
    [1] "b7cc532"
    
    $Packages$beachmat$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$beachmat$Hash
    [1] "f58f6f48893416211092b9600ddeba67"
    
    $Packages$beachmat$Requirements
    $Packages$beachmat$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$beachmat$Requirements[[2]]
    [1] "DelayedArray"
    
    $Packages$beachmat$Requirements[[3]]
    [1] "Matrix"
    
    $Packages$beachmat$Requirements[[4]]
    [1] "Rcpp"
    
    
    
    $Packages$beeswarm
    $Packages$beeswarm$Package
    [1] "beeswarm"
    
    $Packages$beeswarm$Version
    [1] "0.4.0"
    
    $Packages$beeswarm$Source
    [1] "Repository"
    
    $Packages$beeswarm$Repository
    [1] "CRAN"
    
    $Packages$beeswarm$Hash
    [1] "0f4e9d8caa6feaa7e409ae6c30f2ca66"
    
    $Packages$beeswarm$Requirements
    list()
    
    
    $Packages$biomformat
    $Packages$biomformat$Package
    [1] "biomformat"
    
    $Packages$biomformat$Version
    [1] "1.22.0"
    
    $Packages$biomformat$Source
    [1] "Bioconductor"
    
    $Packages$biomformat$git_url
    [1] "https://git.bioconductor.org/packages/biomformat"
    
    $Packages$biomformat$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$biomformat$git_last_commit
    [1] "ab7c641"
    
    $Packages$biomformat$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$biomformat$Hash
    [1] "fd844de531673fbad8baac4a67b9ce48"
    
    $Packages$biomformat$Requirements
    $Packages$biomformat$Requirements[[1]]
    [1] "Matrix"
    
    $Packages$biomformat$Requirements[[2]]
    [1] "jsonlite"
    
    $Packages$biomformat$Requirements[[3]]
    [1] "plyr"
    
    $Packages$biomformat$Requirements[[4]]
    [1] "rhdf5"
    
    
    
    $Packages$bit
    $Packages$bit$Package
    [1] "bit"
    
    $Packages$bit$Version
    [1] "4.0.4"
    
    $Packages$bit$Source
    [1] "Repository"
    
    $Packages$bit$Repository
    [1] "CRAN"
    
    $Packages$bit$Hash
    [1] "f36715f14d94678eea9933af927bc15d"
    
    $Packages$bit$Requirements
    list()
    
    
    $Packages$bit64
    $Packages$bit64$Package
    [1] "bit64"
    
    $Packages$bit64$Version
    [1] "4.0.5"
    
    $Packages$bit64$Source
    [1] "Repository"
    
    $Packages$bit64$Repository
    [1] "CRAN"
    
    $Packages$bit64$Hash
    [1] "9fe98599ca456d6552421db0d6772d8f"
    
    $Packages$bit64$Requirements
    $Packages$bit64$Requirements[[1]]
    [1] "bit"
    
    
    
    $Packages$bitops
    $Packages$bitops$Package
    [1] "bitops"
    
    $Packages$bitops$Version
    [1] "1.0-7"
    
    $Packages$bitops$Source
    [1] "Repository"
    
    $Packages$bitops$Repository
    [1] "CRAN"
    
    $Packages$bitops$Hash
    [1] "b7d8d8ee39869c18d8846a184dd8a1af"
    
    $Packages$bitops$Requirements
    list()
    
    
    $Packages$blob
    $Packages$blob$Package
    [1] "blob"
    
    $Packages$blob$Version
    [1] "1.2.2"
    
    $Packages$blob$Source
    [1] "Repository"
    
    $Packages$blob$Repository
    [1] "CRAN"
    
    $Packages$blob$Hash
    [1] "dc5f7a6598bb025d20d66bb758f12879"
    
    $Packages$blob$Requirements
    $Packages$blob$Requirements[[1]]
    [1] "rlang"
    
    $Packages$blob$Requirements[[2]]
    [1] "vctrs"
    
    
    
    $Packages$brio
    $Packages$brio$Package
    [1] "brio"
    
    $Packages$brio$Version
    [1] "1.1.3"
    
    $Packages$brio$Source
    [1] "Repository"
    
    $Packages$brio$Repository
    [1] "CRAN"
    
    $Packages$brio$Hash
    [1] "976cf154dfb043c012d87cddd8bca363"
    
    $Packages$brio$Requirements
    list()
    
    
    $Packages$broom
    $Packages$broom$Package
    [1] "broom"
    
    $Packages$broom$Version
    [1] "0.7.12"
    
    $Packages$broom$Source
    [1] "Repository"
    
    $Packages$broom$Repository
    [1] "CRAN"
    
    $Packages$broom$Hash
    [1] "bfa8a039d77ae8d5413254e572c8abea"
    
    $Packages$broom$Requirements
    $Packages$broom$Requirements[[1]]
    [1] "backports"
    
    $Packages$broom$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$broom$Requirements[[3]]
    [1] "ellipsis"
    
    $Packages$broom$Requirements[[4]]
    [1] "generics"
    
    $Packages$broom$Requirements[[5]]
    [1] "ggplot2"
    
    $Packages$broom$Requirements[[6]]
    [1] "glue"
    
    $Packages$broom$Requirements[[7]]
    [1] "purrr"
    
    $Packages$broom$Requirements[[8]]
    [1] "rlang"
    
    $Packages$broom$Requirements[[9]]
    [1] "stringr"
    
    $Packages$broom$Requirements[[10]]
    [1] "tibble"
    
    $Packages$broom$Requirements[[11]]
    [1] "tidyr"
    
    
    
    $Packages$bslib
    $Packages$bslib$Package
    [1] "bslib"
    
    $Packages$bslib$Version
    [1] "0.3.1"
    
    $Packages$bslib$Source
    [1] "Repository"
    
    $Packages$bslib$Repository
    [1] "CRAN"
    
    $Packages$bslib$Hash
    [1] "56ae7e1987b340186a8a5a157c2ec358"
    
    $Packages$bslib$Requirements
    $Packages$bslib$Requirements[[1]]
    [1] "htmltools"
    
    $Packages$bslib$Requirements[[2]]
    [1] "jquerylib"
    
    $Packages$bslib$Requirements[[3]]
    [1] "jsonlite"
    
    $Packages$bslib$Requirements[[4]]
    [1] "rlang"
    
    $Packages$bslib$Requirements[[5]]
    [1] "sass"
    
    
    
    $Packages$cachem
    $Packages$cachem$Package
    [1] "cachem"
    
    $Packages$cachem$Version
    [1] "1.0.6"
    
    $Packages$cachem$Source
    [1] "Repository"
    
    $Packages$cachem$Repository
    [1] "CRAN"
    
    $Packages$cachem$Hash
    [1] "648c5b3d71e6a37e3043617489a0a0e9"
    
    $Packages$cachem$Requirements
    $Packages$cachem$Requirements[[1]]
    [1] "fastmap"
    
    $Packages$cachem$Requirements[[2]]
    [1] "rlang"
    
    
    
    $Packages$callr
    $Packages$callr$Package
    [1] "callr"
    
    $Packages$callr$Version
    [1] "3.7.0"
    
    $Packages$callr$Source
    [1] "Repository"
    
    $Packages$callr$Repository
    [1] "CRAN"
    
    $Packages$callr$Hash
    [1] "461aa75a11ce2400245190ef5d3995df"
    
    $Packages$callr$Requirements
    $Packages$callr$Requirements[[1]]
    [1] "R6"
    
    $Packages$callr$Requirements[[2]]
    [1] "processx"
    
    
    
    $Packages$castor
    $Packages$castor$Package
    [1] "castor"
    
    $Packages$castor$Version
    [1] "1.7.2"
    
    $Packages$castor$Source
    [1] "Repository"
    
    $Packages$castor$Repository
    [1] "CRAN"
    
    $Packages$castor$Hash
    [1] "0280b6b1793e9b6c98b727e967b4dd37"
    
    $Packages$castor$Requirements
    $Packages$castor$Requirements[[1]]
    [1] "Matrix"
    
    $Packages$castor$Requirements[[2]]
    [1] "RSpectra"
    
    $Packages$castor$Requirements[[3]]
    [1] "Rcpp"
    
    $Packages$castor$Requirements[[4]]
    [1] "naturalsort"
    
    
    
    $Packages$cellranger
    $Packages$cellranger$Package
    [1] "cellranger"
    
    $Packages$cellranger$Version
    [1] "1.1.0"
    
    $Packages$cellranger$Source
    [1] "Repository"
    
    $Packages$cellranger$Repository
    [1] "CRAN"
    
    $Packages$cellranger$Hash
    [1] "f61dbaec772ccd2e17705c1e872e9e7c"
    
    $Packages$cellranger$Requirements
    $Packages$cellranger$Requirements[[1]]
    [1] "rematch"
    
    $Packages$cellranger$Requirements[[2]]
    [1] "tibble"
    
    
    
    $Packages$checkmate
    $Packages$checkmate$Package
    [1] "checkmate"
    
    $Packages$checkmate$Version
    [1] "2.0.0"
    
    $Packages$checkmate$Source
    [1] "Repository"
    
    $Packages$checkmate$Repository
    [1] "CRAN"
    
    $Packages$checkmate$Hash
    [1] "a667800d5f0350371bedeb8b8b950289"
    
    $Packages$checkmate$Requirements
    $Packages$checkmate$Requirements[[1]]
    [1] "backports"
    
    
    
    $Packages$class
    $Packages$class$Package
    [1] "class"
    
    $Packages$class$Version
    [1] "7.3-20"
    
    $Packages$class$Source
    [1] "Repository"
    
    $Packages$class$Repository
    [1] "CRAN"
    
    $Packages$class$Hash
    [1] "da09d82223e669d270e47ed24ac8686e"
    
    $Packages$class$Requirements
    $Packages$class$Requirements[[1]]
    [1] "MASS"
    
    
    
    $Packages$cli
    $Packages$cli$Package
    [1] "cli"
    
    $Packages$cli$Version
    [1] "3.2.0"
    
    $Packages$cli$Source
    [1] "Repository"
    
    $Packages$cli$Repository
    [1] "CRAN"
    
    $Packages$cli$Hash
    [1] "1bdb126893e9ce6aae50ad1d6fc32faf"
    
    $Packages$cli$Requirements
    $Packages$cli$Requirements[[1]]
    [1] "glue"
    
    
    
    $Packages$clipr
    $Packages$clipr$Package
    [1] "clipr"
    
    $Packages$clipr$Version
    [1] "0.8.0"
    
    $Packages$clipr$Source
    [1] "Repository"
    
    $Packages$clipr$Repository
    [1] "CRAN"
    
    $Packages$clipr$Hash
    [1] "3f038e5ac7f41d4ac41ce658c85e3042"
    
    $Packages$clipr$Requirements
    list()
    
    
    $Packages$clisymbols
    $Packages$clisymbols$Package
    [1] "clisymbols"
    
    $Packages$clisymbols$Version
    [1] "1.2.0"
    
    $Packages$clisymbols$Source
    [1] "Repository"
    
    $Packages$clisymbols$Repository
    [1] "CRAN"
    
    $Packages$clisymbols$Hash
    [1] "96c01552bfd5661b9bbdefbc762f4bcd"
    
    $Packages$clisymbols$Requirements
    list()
    
    
    $Packages$cluster
    $Packages$cluster$Package
    [1] "cluster"
    
    $Packages$cluster$Version
    [1] "2.1.2"
    
    $Packages$cluster$Source
    [1] "Repository"
    
    $Packages$cluster$Repository
    [1] "CRAN"
    
    $Packages$cluster$Hash
    [1] "ce49bfe5bc0b3ecd43a01fe1b01c2243"
    
    $Packages$cluster$Requirements
    list()
    
    
    $Packages$codetools
    $Packages$codetools$Package
    [1] "codetools"
    
    $Packages$codetools$Version
    [1] "0.2-18"
    
    $Packages$codetools$Source
    [1] "Repository"
    
    $Packages$codetools$Repository
    [1] "CRAN"
    
    $Packages$codetools$Hash
    [1] "019388fc48e48b3da0d3a76ff94608a8"
    
    $Packages$codetools$Requirements
    list()
    
    
    $Packages$colorspace
    $Packages$colorspace$Package
    [1] "colorspace"
    
    $Packages$colorspace$Version
    [1] "2.0-3"
    
    $Packages$colorspace$Source
    [1] "Repository"
    
    $Packages$colorspace$Repository
    [1] "CRAN"
    
    $Packages$colorspace$Hash
    [1] "bb4341986bc8b914f0f0acf2e4a3f2f7"
    
    $Packages$colorspace$Requirements
    list()
    
    
    $Packages$commonmark
    $Packages$commonmark$Package
    [1] "commonmark"
    
    $Packages$commonmark$Version
    [1] "1.8.0"
    
    $Packages$commonmark$Source
    [1] "Repository"
    
    $Packages$commonmark$Repository
    [1] "CRAN"
    
    $Packages$commonmark$Hash
    [1] "2ba81b120c1655ab696c935ef33ea716"
    
    $Packages$commonmark$Requirements
    list()
    
    
    $Packages$compositions
    $Packages$compositions$Package
    [1] "compositions"
    
    $Packages$compositions$Version
    [1] "2.0-4"
    
    $Packages$compositions$Source
    [1] "Repository"
    
    $Packages$compositions$Repository
    [1] "CRAN"
    
    $Packages$compositions$Hash
    [1] "db2719c9c2910d9d4a3a9f2ce6c8e815"
    
    $Packages$compositions$Requirements
    $Packages$compositions$Requirements[[1]]
    [1] "MASS"
    
    $Packages$compositions$Requirements[[2]]
    [1] "bayesm"
    
    $Packages$compositions$Requirements[[3]]
    [1] "robustbase"
    
    $Packages$compositions$Requirements[[4]]
    [1] "tensorA"
    
    
    
    $Packages$conflicted
    $Packages$conflicted$Package
    [1] "conflicted"
    
    $Packages$conflicted$Version
    [1] "1.1.0"
    
    $Packages$conflicted$Source
    [1] "Repository"
    
    $Packages$conflicted$Repository
    [1] "CRAN"
    
    $Packages$conflicted$Hash
    [1] "c6bb5e1ef58f2f1c84f238f55bd2e56a"
    
    $Packages$conflicted$Requirements
    $Packages$conflicted$Requirements[[1]]
    [1] "memoise"
    
    $Packages$conflicted$Requirements[[2]]
    [1] "rlang"
    
    
    
    $Packages$contentid
    $Packages$contentid$Package
    [1] "contentid"
    
    $Packages$contentid$Version
    [1] "0.0.15"
    
    $Packages$contentid$Source
    [1] "Repository"
    
    $Packages$contentid$Repository
    [1] "CRAN"
    
    $Packages$contentid$Hash
    [1] "1a6637b36573e009652ff5f98280ac97"
    
    $Packages$contentid$Requirements
    $Packages$contentid$Requirements[[1]]
    [1] "curl"
    
    $Packages$contentid$Requirements[[2]]
    [1] "fs"
    
    $Packages$contentid$Requirements[[3]]
    [1] "httr"
    
    $Packages$contentid$Requirements[[4]]
    [1] "openssl"
    
    
    
    $Packages$cpp11
    $Packages$cpp11$Package
    [1] "cpp11"
    
    $Packages$cpp11$Version
    [1] "0.4.2"
    
    $Packages$cpp11$Source
    [1] "Repository"
    
    $Packages$cpp11$Repository
    [1] "CRAN"
    
    $Packages$cpp11$Hash
    [1] "fa53ce256cd280f468c080a58ea5ba8c"
    
    $Packages$cpp11$Requirements
    list()
    
    
    $Packages$crayon
    $Packages$crayon$Package
    [1] "crayon"
    
    $Packages$crayon$Version
    [1] "1.5.0"
    
    $Packages$crayon$Source
    [1] "Repository"
    
    $Packages$crayon$Repository
    [1] "CRAN"
    
    $Packages$crayon$Hash
    [1] "741c2e098e98afe3dc26a7b0e5489f4e"
    
    $Packages$crayon$Requirements
    list()
    
    
    $Packages$crosstalk
    $Packages$crosstalk$Package
    [1] "crosstalk"
    
    $Packages$crosstalk$Version
    [1] "1.2.0"
    
    $Packages$crosstalk$Source
    [1] "Repository"
    
    $Packages$crosstalk$Repository
    [1] "CRAN"
    
    $Packages$crosstalk$Hash
    [1] "6aa54f69598c32177e920eb3402e8293"
    
    $Packages$crosstalk$Requirements
    $Packages$crosstalk$Requirements[[1]]
    [1] "R6"
    
    $Packages$crosstalk$Requirements[[2]]
    [1] "htmltools"
    
    $Packages$crosstalk$Requirements[[3]]
    [1] "jsonlite"
    
    $Packages$crosstalk$Requirements[[4]]
    [1] "lazyeval"
    
    
    
    $Packages$curatedMetagenomicData
    $Packages$curatedMetagenomicData$Package
    [1] "curatedMetagenomicData"
    
    $Packages$curatedMetagenomicData$Version
    [1] "3.2.3"
    
    $Packages$curatedMetagenomicData$Source
    [1] "Bioconductor"
    
    $Packages$curatedMetagenomicData$git_url
    [1] "https://git.bioconductor.org/packages/curatedMetagenomicData"
    
    $Packages$curatedMetagenomicData$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$curatedMetagenomicData$git_last_commit
    [1] "7e9fe2a"
    
    $Packages$curatedMetagenomicData$git_last_commit_date
    [1] "2021-12-16"
    
    $Packages$curatedMetagenomicData$Hash
    [1] "5cd0cdfeb02aff0d45a23dc59597bf40"
    
    $Packages$curatedMetagenomicData$Requirements
    $Packages$curatedMetagenomicData$Requirements[[1]]
    [1] "AnnotationHub"
    
    $Packages$curatedMetagenomicData$Requirements[[2]]
    [1] "ExperimentHub"
    
    $Packages$curatedMetagenomicData$Requirements[[3]]
    [1] "S4Vectors"
    
    $Packages$curatedMetagenomicData$Requirements[[4]]
    [1] "SummarizedExperiment"
    
    $Packages$curatedMetagenomicData$Requirements[[5]]
    [1] "TreeSummarizedExperiment"
    
    $Packages$curatedMetagenomicData$Requirements[[6]]
    [1] "dplyr"
    
    $Packages$curatedMetagenomicData$Requirements[[7]]
    [1] "magrittr"
    
    $Packages$curatedMetagenomicData$Requirements[[8]]
    [1] "mia"
    
    $Packages$curatedMetagenomicData$Requirements[[9]]
    [1] "purrr"
    
    $Packages$curatedMetagenomicData$Requirements[[10]]
    [1] "rlang"
    
    $Packages$curatedMetagenomicData$Requirements[[11]]
    [1] "stringr"
    
    $Packages$curatedMetagenomicData$Requirements[[12]]
    [1] "tibble"
    
    $Packages$curatedMetagenomicData$Requirements[[13]]
    [1] "tidyr"
    
    $Packages$curatedMetagenomicData$Requirements[[14]]
    [1] "tidyselect"
    
    
    
    $Packages$curl
    $Packages$curl$Package
    [1] "curl"
    
    $Packages$curl$Version
    [1] "4.3.2"
    
    $Packages$curl$Source
    [1] "Repository"
    
    $Packages$curl$Repository
    [1] "CRAN"
    
    $Packages$curl$Hash
    [1] "022c42d49c28e95d69ca60446dbabf88"
    
    $Packages$curl$Requirements
    list()
    
    
    $Packages$dada2
    $Packages$dada2$Package
    [1] "dada2"
    
    $Packages$dada2$Version
    [1] "1.22.0"
    
    $Packages$dada2$Source
    [1] "Bioconductor"
    
    $Packages$dada2$git_url
    [1] "https://git.bioconductor.org/packages/dada2"
    
    $Packages$dada2$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$dada2$git_last_commit
    [1] "3abe06c"
    
    $Packages$dada2$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$dada2$Hash
    [1] "03bd5c0a175bcc98e5c552941c60c70e"
    
    $Packages$dada2$Requirements
    $Packages$dada2$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$dada2$Requirements[[2]]
    [1] "Biostrings"
    
    $Packages$dada2$Requirements[[3]]
    [1] "IRanges"
    
    $Packages$dada2$Requirements[[4]]
    [1] "Rcpp"
    
    $Packages$dada2$Requirements[[5]]
    [1] "RcppParallel"
    
    $Packages$dada2$Requirements[[6]]
    [1] "ShortRead"
    
    $Packages$dada2$Requirements[[7]]
    [1] "XVector"
    
    $Packages$dada2$Requirements[[8]]
    [1] "ggplot2"
    
    $Packages$dada2$Requirements[[9]]
    [1] "reshape2"
    
    
    
    $Packages$data.table
    $Packages$data.table$Package
    [1] "data.table"
    
    $Packages$data.table$Version
    [1] "1.14.2"
    
    $Packages$data.table$Source
    [1] "Repository"
    
    $Packages$data.table$Repository
    [1] "CRAN"
    
    $Packages$data.table$Hash
    [1] "36b67b5adf57b292923f5659f5f0c853"
    
    $Packages$data.table$Requirements
    list()
    
    
    $Packages$dbplyr
    $Packages$dbplyr$Package
    [1] "dbplyr"
    
    $Packages$dbplyr$Version
    [1] "2.1.1"
    
    $Packages$dbplyr$Source
    [1] "Repository"
    
    $Packages$dbplyr$Repository
    [1] "CRAN"
    
    $Packages$dbplyr$Hash
    [1] "1f37fa4ab2f5f7eded42f78b9a887182"
    
    $Packages$dbplyr$Requirements
    $Packages$dbplyr$Requirements[[1]]
    [1] "DBI"
    
    $Packages$dbplyr$Requirements[[2]]
    [1] "R6"
    
    $Packages$dbplyr$Requirements[[3]]
    [1] "assertthat"
    
    $Packages$dbplyr$Requirements[[4]]
    [1] "blob"
    
    $Packages$dbplyr$Requirements[[5]]
    [1] "dplyr"
    
    $Packages$dbplyr$Requirements[[6]]
    [1] "ellipsis"
    
    $Packages$dbplyr$Requirements[[7]]
    [1] "glue"
    
    $Packages$dbplyr$Requirements[[8]]
    [1] "lifecycle"
    
    $Packages$dbplyr$Requirements[[9]]
    [1] "magrittr"
    
    $Packages$dbplyr$Requirements[[10]]
    [1] "purrr"
    
    $Packages$dbplyr$Requirements[[11]]
    [1] "rlang"
    
    $Packages$dbplyr$Requirements[[12]]
    [1] "tibble"
    
    $Packages$dbplyr$Requirements[[13]]
    [1] "tidyselect"
    
    $Packages$dbplyr$Requirements[[14]]
    [1] "vctrs"
    
    $Packages$dbplyr$Requirements[[15]]
    [1] "withr"
    
    
    
    $Packages$decontam
    $Packages$decontam$Package
    [1] "decontam"
    
    $Packages$decontam$Version
    [1] "1.14.0"
    
    $Packages$decontam$Source
    [1] "Bioconductor"
    
    $Packages$decontam$git_url
    [1] "https://git.bioconductor.org/packages/decontam"
    
    $Packages$decontam$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$decontam$git_last_commit
    [1] "b710769"
    
    $Packages$decontam$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$decontam$Hash
    [1] "d3b0675086a9a97173ea8482b15a18b9"
    
    $Packages$decontam$Requirements
    $Packages$decontam$Requirements[[1]]
    [1] "ggplot2"
    
    $Packages$decontam$Requirements[[2]]
    [1] "reshape2"
    
    
    
    $Packages$desc
    $Packages$desc$Package
    [1] "desc"
    
    $Packages$desc$Version
    [1] "1.4.1"
    
    $Packages$desc$Source
    [1] "Repository"
    
    $Packages$desc$Repository
    [1] "CRAN"
    
    $Packages$desc$Hash
    [1] "eebd27ee58fcc58714eedb7aa07d8ad1"
    
    $Packages$desc$Requirements
    $Packages$desc$Requirements[[1]]
    [1] "R6"
    
    $Packages$desc$Requirements[[2]]
    [1] "cli"
    
    $Packages$desc$Requirements[[3]]
    [1] "rprojroot"
    
    
    
    $Packages$dials
    $Packages$dials$Package
    [1] "dials"
    
    $Packages$dials$Version
    [1] "0.1.0"
    
    $Packages$dials$Source
    [1] "Repository"
    
    $Packages$dials$Repository
    [1] "CRAN"
    
    $Packages$dials$Hash
    [1] "38d457070f3093b6b92fc2413a513f74"
    
    $Packages$dials$Requirements
    $Packages$dials$Requirements[[1]]
    [1] "DiceDesign"
    
    $Packages$dials$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$dials$Requirements[[3]]
    [1] "glue"
    
    $Packages$dials$Requirements[[4]]
    [1] "hardhat"
    
    $Packages$dials$Requirements[[5]]
    [1] "lifecycle"
    
    $Packages$dials$Requirements[[6]]
    [1] "purrr"
    
    $Packages$dials$Requirements[[7]]
    [1] "rlang"
    
    $Packages$dials$Requirements[[8]]
    [1] "scales"
    
    $Packages$dials$Requirements[[9]]
    [1] "tibble"
    
    $Packages$dials$Requirements[[10]]
    [1] "vctrs"
    
    $Packages$dials$Requirements[[11]]
    [1] "withr"
    
    
    
    $Packages$diffobj
    $Packages$diffobj$Package
    [1] "diffobj"
    
    $Packages$diffobj$Version
    [1] "0.3.5"
    
    $Packages$diffobj$Source
    [1] "Repository"
    
    $Packages$diffobj$Repository
    [1] "CRAN"
    
    $Packages$diffobj$Hash
    [1] "bcaa8b95f8d7d01a5dedfd959ce88ab8"
    
    $Packages$diffobj$Requirements
    $Packages$diffobj$Requirements[[1]]
    [1] "crayon"
    
    
    
    $Packages$digest
    $Packages$digest$Package
    [1] "digest"
    
    $Packages$digest$Version
    [1] "0.6.29"
    
    $Packages$digest$Source
    [1] "Repository"
    
    $Packages$digest$Repository
    [1] "CRAN"
    
    $Packages$digest$Hash
    [1] "cf6b206a045a684728c3267ef7596190"
    
    $Packages$digest$Requirements
    list()
    
    
    $Packages$doParallel
    $Packages$doParallel$Package
    [1] "doParallel"
    
    $Packages$doParallel$Version
    [1] "1.0.17"
    
    $Packages$doParallel$Source
    [1] "Repository"
    
    $Packages$doParallel$Repository
    [1] "CRAN"
    
    $Packages$doParallel$Hash
    [1] "451e5edf411987991ab6a5410c45011f"
    
    $Packages$doParallel$Requirements
    $Packages$doParallel$Requirements[[1]]
    [1] "foreach"
    
    $Packages$doParallel$Requirements[[2]]
    [1] "iterators"
    
    
    
    $Packages$dplyr
    $Packages$dplyr$Package
    [1] "dplyr"
    
    $Packages$dplyr$Version
    [1] "1.0.8"
    
    $Packages$dplyr$Source
    [1] "Repository"
    
    $Packages$dplyr$Repository
    [1] "CRAN"
    
    $Packages$dplyr$Hash
    [1] "ef47665e64228a17609d6df877bf86f2"
    
    $Packages$dplyr$Requirements
    $Packages$dplyr$Requirements[[1]]
    [1] "R6"
    
    $Packages$dplyr$Requirements[[2]]
    [1] "generics"
    
    $Packages$dplyr$Requirements[[3]]
    [1] "glue"
    
    $Packages$dplyr$Requirements[[4]]
    [1] "lifecycle"
    
    $Packages$dplyr$Requirements[[5]]
    [1] "magrittr"
    
    $Packages$dplyr$Requirements[[6]]
    [1] "pillar"
    
    $Packages$dplyr$Requirements[[7]]
    [1] "rlang"
    
    $Packages$dplyr$Requirements[[8]]
    [1] "tibble"
    
    $Packages$dplyr$Requirements[[9]]
    [1] "tidyselect"
    
    $Packages$dplyr$Requirements[[10]]
    [1] "vctrs"
    
    
    
    $Packages$dtplyr
    $Packages$dtplyr$Package
    [1] "dtplyr"
    
    $Packages$dtplyr$Version
    [1] "1.2.1"
    
    $Packages$dtplyr$Source
    [1] "Repository"
    
    $Packages$dtplyr$Repository
    [1] "CRAN"
    
    $Packages$dtplyr$Hash
    [1] "f5d195cd5fcc0a77499d9da698ef2ea3"
    
    $Packages$dtplyr$Requirements
    $Packages$dtplyr$Requirements[[1]]
    [1] "crayon"
    
    $Packages$dtplyr$Requirements[[2]]
    [1] "data.table"
    
    $Packages$dtplyr$Requirements[[3]]
    [1] "dplyr"
    
    $Packages$dtplyr$Requirements[[4]]
    [1] "ellipsis"
    
    $Packages$dtplyr$Requirements[[5]]
    [1] "glue"
    
    $Packages$dtplyr$Requirements[[6]]
    [1] "lifecycle"
    
    $Packages$dtplyr$Requirements[[7]]
    [1] "rlang"
    
    $Packages$dtplyr$Requirements[[8]]
    [1] "tibble"
    
    $Packages$dtplyr$Requirements[[9]]
    [1] "tidyselect"
    
    $Packages$dtplyr$Requirements[[10]]
    [1] "vctrs"
    
    
    
    $Packages$duckdb
    $Packages$duckdb$Package
    [1] "duckdb"
    
    $Packages$duckdb$Version
    [1] "0.3.2-1"
    
    $Packages$duckdb$Source
    [1] "Repository"
    
    $Packages$duckdb$Repository
    [1] "CRAN"
    
    $Packages$duckdb$Hash
    [1] "29d4f41f42b2cf0feecd265f8fd8412f"
    
    $Packages$duckdb$Requirements
    $Packages$duckdb$Requirements[[1]]
    [1] "DBI"
    
    
    
    $Packages$ellipsis
    $Packages$ellipsis$Package
    [1] "ellipsis"
    
    $Packages$ellipsis$Version
    [1] "0.3.2"
    
    $Packages$ellipsis$Source
    [1] "Repository"
    
    $Packages$ellipsis$Repository
    [1] "CRAN"
    
    $Packages$ellipsis$Hash
    [1] "bb0eec2fe32e88d9e2836c2f73ea2077"
    
    $Packages$ellipsis$Requirements
    $Packages$ellipsis$Requirements[[1]]
    [1] "rlang"
    
    
    
    $Packages$evaluate
    $Packages$evaluate$Package
    [1] "evaluate"
    
    $Packages$evaluate$Version
    [1] "0.15"
    
    $Packages$evaluate$Source
    [1] "Repository"
    
    $Packages$evaluate$Repository
    [1] "CRAN"
    
    $Packages$evaluate$Hash
    [1] "699a7a93d08c962d9f8950b2d7a227f1"
    
    $Packages$evaluate$Requirements
    list()
    
    
    $Packages$fansi
    $Packages$fansi$Package
    [1] "fansi"
    
    $Packages$fansi$Version
    [1] "1.0.2"
    
    $Packages$fansi$Source
    [1] "Repository"
    
    $Packages$fansi$Repository
    [1] "CRAN"
    
    $Packages$fansi$Hash
    [1] "f28149c2d7a1342a834b314e95e67260"
    
    $Packages$fansi$Requirements
    list()
    
    
    $Packages$farver
    $Packages$farver$Package
    [1] "farver"
    
    $Packages$farver$Version
    [1] "2.1.0"
    
    $Packages$farver$Source
    [1] "Repository"
    
    $Packages$farver$Repository
    [1] "CRAN"
    
    $Packages$farver$Hash
    [1] "c98eb5133d9cb9e1622b8691487f11bb"
    
    $Packages$farver$Requirements
    list()
    
    
    $Packages$fastmap
    $Packages$fastmap$Package
    [1] "fastmap"
    
    $Packages$fastmap$Version
    [1] "1.1.0"
    
    $Packages$fastmap$Source
    [1] "Repository"
    
    $Packages$fastmap$Repository
    [1] "CRAN"
    
    $Packages$fastmap$Hash
    [1] "77bd60a6157420d4ffa93b27cf6a58b8"
    
    $Packages$fastmap$Requirements
    list()
    
    
    $Packages$fastmatch
    $Packages$fastmatch$Package
    [1] "fastmatch"
    
    $Packages$fastmatch$Version
    [1] "1.1-3"
    
    $Packages$fastmatch$Source
    [1] "Repository"
    
    $Packages$fastmatch$Repository
    [1] "CRAN"
    
    $Packages$fastmatch$Hash
    [1] "dabc225759a2c2b241e60e42bf0e8e54"
    
    $Packages$fastmatch$Requirements
    list()
    
    
    $Packages$fgsea
    $Packages$fgsea$Package
    [1] "fgsea"
    
    $Packages$fgsea$Version
    [1] "1.20.0"
    
    $Packages$fgsea$Source
    [1] "Bioconductor"
    
    $Packages$fgsea$git_url
    [1] "https://git.bioconductor.org/packages/fgsea"
    
    $Packages$fgsea$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$fgsea$git_last_commit
    [1] "b704f81"
    
    $Packages$fgsea$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$fgsea$Hash
    [1] "a4913104166e73957fa32d2e595cceef"
    
    $Packages$fgsea$Requirements
    $Packages$fgsea$Requirements[[1]]
    [1] "BH"
    
    $Packages$fgsea$Requirements[[2]]
    [1] "BiocParallel"
    
    $Packages$fgsea$Requirements[[3]]
    [1] "Matrix"
    
    $Packages$fgsea$Requirements[[4]]
    [1] "Rcpp"
    
    $Packages$fgsea$Requirements[[5]]
    [1] "data.table"
    
    $Packages$fgsea$Requirements[[6]]
    [1] "fastmatch"
    
    $Packages$fgsea$Requirements[[7]]
    [1] "ggplot2"
    
    $Packages$fgsea$Requirements[[8]]
    [1] "gridExtra"
    
    
    
    $Packages$filelock
    $Packages$filelock$Package
    [1] "filelock"
    
    $Packages$filelock$Version
    [1] "1.0.2"
    
    $Packages$filelock$Source
    [1] "Repository"
    
    $Packages$filelock$Repository
    [1] "CRAN"
    
    $Packages$filelock$Hash
    [1] "38ec653c2613bed60052ba3787bd8a2c"
    
    $Packages$filelock$Requirements
    list()
    
    
    $Packages$fontawesome
    $Packages$fontawesome$Package
    [1] "fontawesome"
    
    $Packages$fontawesome$Version
    [1] "0.2.2"
    
    $Packages$fontawesome$Source
    [1] "Repository"
    
    $Packages$fontawesome$Repository
    [1] "CRAN"
    
    $Packages$fontawesome$Hash
    [1] "55624ed409e46c5f358b2c060be87f67"
    
    $Packages$fontawesome$Requirements
    $Packages$fontawesome$Requirements[[1]]
    [1] "htmltools"
    
    $Packages$fontawesome$Requirements[[2]]
    [1] "rlang"
    
    
    
    $Packages$forcats
    $Packages$forcats$Package
    [1] "forcats"
    
    $Packages$forcats$Version
    [1] "0.5.1"
    
    $Packages$forcats$Source
    [1] "Repository"
    
    $Packages$forcats$Repository
    [1] "CRAN"
    
    $Packages$forcats$Hash
    [1] "81c3244cab67468aac4c60550832655d"
    
    $Packages$forcats$Requirements
    $Packages$forcats$Requirements[[1]]
    [1] "ellipsis"
    
    $Packages$forcats$Requirements[[2]]
    [1] "magrittr"
    
    $Packages$forcats$Requirements[[3]]
    [1] "rlang"
    
    $Packages$forcats$Requirements[[4]]
    [1] "tibble"
    
    
    
    $Packages$foreach
    $Packages$foreach$Package
    [1] "foreach"
    
    $Packages$foreach$Version
    [1] "1.5.2"
    
    $Packages$foreach$Source
    [1] "Repository"
    
    $Packages$foreach$Repository
    [1] "CRAN"
    
    $Packages$foreach$Hash
    [1] "618609b42c9406731ead03adf5379850"
    
    $Packages$foreach$Requirements
    $Packages$foreach$Requirements[[1]]
    [1] "codetools"
    
    $Packages$foreach$Requirements[[2]]
    [1] "iterators"
    
    
    
    $Packages$formatR
    $Packages$formatR$Package
    [1] "formatR"
    
    $Packages$formatR$Version
    [1] "1.11"
    
    $Packages$formatR$Source
    [1] "Repository"
    
    $Packages$formatR$Repository
    [1] "CRAN"
    
    $Packages$formatR$Hash
    [1] "2590a6a868515a69f258640a51b724c9"
    
    $Packages$formatR$Requirements
    list()
    
    
    $Packages$fs
    $Packages$fs$Package
    [1] "fs"
    
    $Packages$fs$Version
    [1] "1.5.2"
    
    $Packages$fs$Source
    [1] "Repository"
    
    $Packages$fs$Repository
    [1] "CRAN"
    
    $Packages$fs$Hash
    [1] "7c89603d81793f0d5486d91ab1fc6f1d"
    
    $Packages$fs$Requirements
    list()
    
    
    $Packages$furrr
    $Packages$furrr$Package
    [1] "furrr"
    
    $Packages$furrr$Version
    [1] "0.2.3"
    
    $Packages$furrr$Source
    [1] "Repository"
    
    $Packages$furrr$Repository
    [1] "CRAN"
    
    $Packages$furrr$Hash
    [1] "2aba4ab06e8707ac054c6683cb6fed56"
    
    $Packages$furrr$Requirements
    $Packages$furrr$Requirements[[1]]
    [1] "ellipsis"
    
    $Packages$furrr$Requirements[[2]]
    [1] "future"
    
    $Packages$furrr$Requirements[[3]]
    [1] "globals"
    
    $Packages$furrr$Requirements[[4]]
    [1] "lifecycle"
    
    $Packages$furrr$Requirements[[5]]
    [1] "purrr"
    
    $Packages$furrr$Requirements[[6]]
    [1] "rlang"
    
    $Packages$furrr$Requirements[[7]]
    [1] "vctrs"
    
    
    
    $Packages$futile.logger
    $Packages$futile.logger$Package
    [1] "futile.logger"
    
    $Packages$futile.logger$Version
    [1] "1.4.3"
    
    $Packages$futile.logger$Source
    [1] "Repository"
    
    $Packages$futile.logger$Repository
    [1] "CRAN"
    
    $Packages$futile.logger$Hash
    [1] "99f0ace8c05ec7d3683d27083c4f1e7e"
    
    $Packages$futile.logger$Requirements
    $Packages$futile.logger$Requirements[[1]]
    [1] "futile.options"
    
    $Packages$futile.logger$Requirements[[2]]
    [1] "lambda.r"
    
    
    
    $Packages$futile.options
    $Packages$futile.options$Package
    [1] "futile.options"
    
    $Packages$futile.options$Version
    [1] "1.0.1"
    
    $Packages$futile.options$Source
    [1] "Repository"
    
    $Packages$futile.options$Repository
    [1] "CRAN"
    
    $Packages$futile.options$Hash
    [1] "0d9bf02413ddc2bbe8da9ce369dcdd2b"
    
    $Packages$futile.options$Requirements
    list()
    
    
    $Packages$future
    $Packages$future$Package
    [1] "future"
    
    $Packages$future$Version
    [1] "1.24.0"
    
    $Packages$future$Source
    [1] "Repository"
    
    $Packages$future$Repository
    [1] "CRAN"
    
    $Packages$future$Hash
    [1] "5cc7addaa73372fbee0a7d06c880068e"
    
    $Packages$future$Requirements
    $Packages$future$Requirements[[1]]
    [1] "digest"
    
    $Packages$future$Requirements[[2]]
    [1] "globals"
    
    $Packages$future$Requirements[[3]]
    [1] "listenv"
    
    $Packages$future$Requirements[[4]]
    [1] "parallelly"
    
    
    
    $Packages$future.apply
    $Packages$future.apply$Package
    [1] "future.apply"
    
    $Packages$future.apply$Version
    [1] "1.8.1"
    
    $Packages$future.apply$Source
    [1] "Repository"
    
    $Packages$future.apply$Repository
    [1] "CRAN"
    
    $Packages$future.apply$Hash
    [1] "f568ce73d3d59582b0f7babd0eb33d07"
    
    $Packages$future.apply$Requirements
    $Packages$future.apply$Requirements[[1]]
    [1] "future"
    
    $Packages$future.apply$Requirements[[2]]
    [1] "globals"
    
    
    
    $Packages$gargle
    $Packages$gargle$Package
    [1] "gargle"
    
    $Packages$gargle$Version
    [1] "1.2.0"
    
    $Packages$gargle$Source
    [1] "Repository"
    
    $Packages$gargle$Repository
    [1] "CRAN"
    
    $Packages$gargle$Hash
    [1] "9d234e6a87a6f8181792de6dc4a00e39"
    
    $Packages$gargle$Requirements
    $Packages$gargle$Requirements[[1]]
    [1] "cli"
    
    $Packages$gargle$Requirements[[2]]
    [1] "fs"
    
    $Packages$gargle$Requirements[[3]]
    [1] "glue"
    
    $Packages$gargle$Requirements[[4]]
    [1] "httr"
    
    $Packages$gargle$Requirements[[5]]
    [1] "jsonlite"
    
    $Packages$gargle$Requirements[[6]]
    [1] "rappdirs"
    
    $Packages$gargle$Requirements[[7]]
    [1] "rlang"
    
    $Packages$gargle$Requirements[[8]]
    [1] "rstudioapi"
    
    $Packages$gargle$Requirements[[9]]
    [1] "withr"
    
    
    
    $Packages$generics
    $Packages$generics$Package
    [1] "generics"
    
    $Packages$generics$Version
    [1] "0.1.2"
    
    $Packages$generics$Source
    [1] "Repository"
    
    $Packages$generics$Repository
    [1] "CRAN"
    
    $Packages$generics$Hash
    [1] "177475892cf4a55865868527654a7741"
    
    $Packages$generics$Requirements
    list()
    
    
    $Packages$getopt
    $Packages$getopt$Package
    [1] "getopt"
    
    $Packages$getopt$Version
    [1] "1.20.3"
    
    $Packages$getopt$Source
    [1] "Repository"
    
    $Packages$getopt$Repository
    [1] "CRAN"
    
    $Packages$getopt$Hash
    [1] "ad68e3263f0bc9269aab8c2039440117"
    
    $Packages$getopt$Requirements
    list()
    
    
    $Packages$ggbeeswarm
    $Packages$ggbeeswarm$Package
    [1] "ggbeeswarm"
    
    $Packages$ggbeeswarm$Version
    [1] "0.6.0"
    
    $Packages$ggbeeswarm$Source
    [1] "Repository"
    
    $Packages$ggbeeswarm$Repository
    [1] "CRAN"
    
    $Packages$ggbeeswarm$Hash
    [1] "dd68b9b215b2d3119603549a794003c3"
    
    $Packages$ggbeeswarm$Requirements
    $Packages$ggbeeswarm$Requirements[[1]]
    [1] "beeswarm"
    
    $Packages$ggbeeswarm$Requirements[[2]]
    [1] "ggplot2"
    
    $Packages$ggbeeswarm$Requirements[[3]]
    [1] "vipor"
    
    
    
    $Packages$ggfun
    $Packages$ggfun$Package
    [1] "ggfun"
    
    $Packages$ggfun$Version
    [1] "0.0.5"
    
    $Packages$ggfun$Source
    [1] "Repository"
    
    $Packages$ggfun$Repository
    [1] "CRAN"
    
    $Packages$ggfun$Hash
    [1] "be6dbe9ddc9cc77fb9277611eff20a96"
    
    $Packages$ggfun$Requirements
    $Packages$ggfun$Requirements[[1]]
    [1] "ggplot2"
    
    $Packages$ggfun$Requirements[[2]]
    [1] "rlang"
    
    
    
    $Packages$ggplot2
    $Packages$ggplot2$Package
    [1] "ggplot2"
    
    $Packages$ggplot2$Version
    [1] "3.3.5"
    
    $Packages$ggplot2$Source
    [1] "Repository"
    
    $Packages$ggplot2$Repository
    [1] "CRAN"
    
    $Packages$ggplot2$Hash
    [1] "d7566c471c7b17e095dd023b9ef155ad"
    
    $Packages$ggplot2$Requirements
    $Packages$ggplot2$Requirements[[1]]
    [1] "MASS"
    
    $Packages$ggplot2$Requirements[[2]]
    [1] "digest"
    
    $Packages$ggplot2$Requirements[[3]]
    [1] "glue"
    
    $Packages$ggplot2$Requirements[[4]]
    [1] "gtable"
    
    $Packages$ggplot2$Requirements[[5]]
    [1] "isoband"
    
    $Packages$ggplot2$Requirements[[6]]
    [1] "mgcv"
    
    $Packages$ggplot2$Requirements[[7]]
    [1] "rlang"
    
    $Packages$ggplot2$Requirements[[8]]
    [1] "scales"
    
    $Packages$ggplot2$Requirements[[9]]
    [1] "tibble"
    
    $Packages$ggplot2$Requirements[[10]]
    [1] "withr"
    
    
    
    $Packages$ggplotify
    $Packages$ggplotify$Package
    [1] "ggplotify"
    
    $Packages$ggplotify$Version
    [1] "0.1.0"
    
    $Packages$ggplotify$Source
    [1] "Repository"
    
    $Packages$ggplotify$Repository
    [1] "CRAN"
    
    $Packages$ggplotify$Hash
    [1] "acbcedf783cdb8710168aa0edba42ac0"
    
    $Packages$ggplotify$Requirements
    $Packages$ggplotify$Requirements[[1]]
    [1] "ggplot2"
    
    $Packages$ggplotify$Requirements[[2]]
    [1] "gridGraphics"
    
    $Packages$ggplotify$Requirements[[3]]
    [1] "yulab.utils"
    
    
    
    $Packages$ggrepel
    $Packages$ggrepel$Package
    [1] "ggrepel"
    
    $Packages$ggrepel$Version
    [1] "0.9.1"
    
    $Packages$ggrepel$Source
    [1] "Repository"
    
    $Packages$ggrepel$Repository
    [1] "CRAN"
    
    $Packages$ggrepel$Hash
    [1] "08ab869f37e6a7741a64ab9069bcb67d"
    
    $Packages$ggrepel$Requirements
    $Packages$ggrepel$Requirements[[1]]
    [1] "Rcpp"
    
    $Packages$ggrepel$Requirements[[2]]
    [1] "ggplot2"
    
    $Packages$ggrepel$Requirements[[3]]
    [1] "rlang"
    
    $Packages$ggrepel$Requirements[[4]]
    [1] "scales"
    
    
    
    $Packages$ggsci
    $Packages$ggsci$Package
    [1] "ggsci"
    
    $Packages$ggsci$Version
    [1] "2.9"
    
    $Packages$ggsci$Source
    [1] "Repository"
    
    $Packages$ggsci$Repository
    [1] "CRAN"
    
    $Packages$ggsci$Hash
    [1] "81ccb8213ed592598210afd10c3a5936"
    
    $Packages$ggsci$Requirements
    $Packages$ggsci$Requirements[[1]]
    [1] "ggplot2"
    
    $Packages$ggsci$Requirements[[2]]
    [1] "scales"
    
    
    
    $Packages$ggthemes
    $Packages$ggthemes$Package
    [1] "ggthemes"
    
    $Packages$ggthemes$Version
    [1] "4.2.4"
    
    $Packages$ggthemes$Source
    [1] "Repository"
    
    $Packages$ggthemes$Repository
    [1] "CRAN"
    
    $Packages$ggthemes$Hash
    [1] "617fb2c300f68e8b1ba21043fba88fd7"
    
    $Packages$ggthemes$Requirements
    $Packages$ggthemes$Requirements[[1]]
    [1] "ggplot2"
    
    $Packages$ggthemes$Requirements[[2]]
    [1] "purrr"
    
    $Packages$ggthemes$Requirements[[3]]
    [1] "scales"
    
    $Packages$ggthemes$Requirements[[4]]
    [1] "stringr"
    
    $Packages$ggthemes$Requirements[[5]]
    [1] "tibble"
    
    
    
    $Packages$ggtree
    $Packages$ggtree$Package
    [1] "ggtree"
    
    $Packages$ggtree$Version
    [1] "3.2.1"
    
    $Packages$ggtree$Source
    [1] "Bioconductor"
    
    $Packages$ggtree$RemoteType
    [1] "bioconductor"
    
    $Packages$ggtree$Remotes
    [1] "GuangchuangYu/treeio"
    
    $Packages$ggtree$git_url
    [1] "https://git.bioconductor.org/packages/ggtree"
    
    $Packages$ggtree$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$ggtree$git_last_commit
    [1] "d3747e6"
    
    $Packages$ggtree$git_last_commit_date
    [1] "2021-11-14"
    
    $Packages$ggtree$Hash
    [1] "5711c057a04e53ed1c70909939dd9ad9"
    
    $Packages$ggtree$Requirements
    $Packages$ggtree$Requirements[[1]]
    [1] "ape"
    
    $Packages$ggtree$Requirements[[2]]
    [1] "aplot"
    
    $Packages$ggtree$Requirements[[3]]
    [1] "dplyr"
    
    $Packages$ggtree$Requirements[[4]]
    [1] "ggfun"
    
    $Packages$ggtree$Requirements[[5]]
    [1] "ggplot2"
    
    $Packages$ggtree$Requirements[[6]]
    [1] "magrittr"
    
    $Packages$ggtree$Requirements[[7]]
    [1] "purrr"
    
    $Packages$ggtree$Requirements[[8]]
    [1] "rlang"
    
    $Packages$ggtree$Requirements[[9]]
    [1] "scales"
    
    $Packages$ggtree$Requirements[[10]]
    [1] "tidyr"
    
    $Packages$ggtree$Requirements[[11]]
    [1] "tidytree"
    
    $Packages$ggtree$Requirements[[12]]
    [1] "treeio"
    
    $Packages$ggtree$Requirements[[13]]
    [1] "yulab.utils"
    
    
    
    $Packages$gh
    $Packages$gh$Package
    [1] "gh"
    
    $Packages$gh$Version
    [1] "1.3.0"
    
    $Packages$gh$Source
    [1] "Repository"
    
    $Packages$gh$Repository
    [1] "CRAN"
    
    $Packages$gh$Hash
    [1] "38c2580abbda249bd6afeec00d14f531"
    
    $Packages$gh$Requirements
    $Packages$gh$Requirements[[1]]
    [1] "cli"
    
    $Packages$gh$Requirements[[2]]
    [1] "gitcreds"
    
    $Packages$gh$Requirements[[3]]
    [1] "httr"
    
    $Packages$gh$Requirements[[4]]
    [1] "ini"
    
    $Packages$gh$Requirements[[5]]
    [1] "jsonlite"
    
    
    
    $Packages$gitcreds
    $Packages$gitcreds$Package
    [1] "gitcreds"
    
    $Packages$gitcreds$Version
    [1] "0.1.1"
    
    $Packages$gitcreds$Source
    [1] "Repository"
    
    $Packages$gitcreds$Repository
    [1] "CRAN"
    
    $Packages$gitcreds$Hash
    [1] "f3aefccc1cc50de6338146b62f115de8"
    
    $Packages$gitcreds$Requirements
    list()
    
    
    $Packages$globals
    $Packages$globals$Package
    [1] "globals"
    
    $Packages$globals$Version
    [1] "0.14.0"
    
    $Packages$globals$Source
    [1] "Repository"
    
    $Packages$globals$Repository
    [1] "CRAN"
    
    $Packages$globals$Hash
    [1] "eca8023ed5ca6372479ebb9b3207f5ae"
    
    $Packages$globals$Requirements
    $Packages$globals$Requirements[[1]]
    [1] "codetools"
    
    
    
    $Packages$glue
    $Packages$glue$Package
    [1] "glue"
    
    $Packages$glue$Version
    [1] "1.6.2"
    
    $Packages$glue$Source
    [1] "Repository"
    
    $Packages$glue$Repository
    [1] "CRAN"
    
    $Packages$glue$Hash
    [1] "4f2596dfb05dac67b9dc558e5c6fba2e"
    
    $Packages$glue$Requirements
    list()
    
    
    $Packages$googledrive
    $Packages$googledrive$Package
    [1] "googledrive"
    
    $Packages$googledrive$Version
    [1] "2.0.0"
    
    $Packages$googledrive$Source
    [1] "Repository"
    
    $Packages$googledrive$Repository
    [1] "CRAN"
    
    $Packages$googledrive$Hash
    [1] "c3a25adbbfbb03f12e6f88c5fb1f3024"
    
    $Packages$googledrive$Requirements
    $Packages$googledrive$Requirements[[1]]
    [1] "cli"
    
    $Packages$googledrive$Requirements[[2]]
    [1] "gargle"
    
    $Packages$googledrive$Requirements[[3]]
    [1] "glue"
    
    $Packages$googledrive$Requirements[[4]]
    [1] "httr"
    
    $Packages$googledrive$Requirements[[5]]
    [1] "jsonlite"
    
    $Packages$googledrive$Requirements[[6]]
    [1] "lifecycle"
    
    $Packages$googledrive$Requirements[[7]]
    [1] "magrittr"
    
    $Packages$googledrive$Requirements[[8]]
    [1] "pillar"
    
    $Packages$googledrive$Requirements[[9]]
    [1] "purrr"
    
    $Packages$googledrive$Requirements[[10]]
    [1] "rlang"
    
    $Packages$googledrive$Requirements[[11]]
    [1] "tibble"
    
    $Packages$googledrive$Requirements[[12]]
    [1] "uuid"
    
    $Packages$googledrive$Requirements[[13]]
    [1] "vctrs"
    
    $Packages$googledrive$Requirements[[14]]
    [1] "withr"
    
    
    
    $Packages$googlesheets4
    $Packages$googlesheets4$Package
    [1] "googlesheets4"
    
    $Packages$googlesheets4$Version
    [1] "1.0.0"
    
    $Packages$googlesheets4$Source
    [1] "Repository"
    
    $Packages$googlesheets4$Repository
    [1] "CRAN"
    
    $Packages$googlesheets4$Hash
    [1] "9a6564184dc4a81daea4f1d7ce357c6a"
    
    $Packages$googlesheets4$Requirements
    $Packages$googlesheets4$Requirements[[1]]
    [1] "cellranger"
    
    $Packages$googlesheets4$Requirements[[2]]
    [1] "cli"
    
    $Packages$googlesheets4$Requirements[[3]]
    [1] "curl"
    
    $Packages$googlesheets4$Requirements[[4]]
    [1] "gargle"
    
    $Packages$googlesheets4$Requirements[[5]]
    [1] "glue"
    
    $Packages$googlesheets4$Requirements[[6]]
    [1] "googledrive"
    
    $Packages$googlesheets4$Requirements[[7]]
    [1] "httr"
    
    $Packages$googlesheets4$Requirements[[8]]
    [1] "ids"
    
    $Packages$googlesheets4$Requirements[[9]]
    [1] "magrittr"
    
    $Packages$googlesheets4$Requirements[[10]]
    [1] "purrr"
    
    $Packages$googlesheets4$Requirements[[11]]
    [1] "rematch2"
    
    $Packages$googlesheets4$Requirements[[12]]
    [1] "rlang"
    
    $Packages$googlesheets4$Requirements[[13]]
    [1] "tibble"
    
    $Packages$googlesheets4$Requirements[[14]]
    [1] "vctrs"
    
    
    
    $Packages$gower
    $Packages$gower$Package
    [1] "gower"
    
    $Packages$gower$Version
    [1] "1.0.0"
    
    $Packages$gower$Source
    [1] "Repository"
    
    $Packages$gower$Repository
    [1] "CRAN"
    
    $Packages$gower$Hash
    [1] "e30deda901954a80e035ad2e972ba7fd"
    
    $Packages$gower$Requirements
    list()
    
    
    $Packages$graph
    $Packages$graph$Package
    [1] "graph"
    
    $Packages$graph$Version
    [1] "1.72.0"
    
    $Packages$graph$Source
    [1] "Bioconductor"
    
    $Packages$graph$git_url
    [1] "https://git.bioconductor.org/packages/graph"
    
    $Packages$graph$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$graph$git_last_commit
    [1] "7afbd26"
    
    $Packages$graph$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$graph$Hash
    [1] "008757628d8e6bacfe8f122954d8f024"
    
    $Packages$graph$Requirements
    $Packages$graph$Requirements[[1]]
    [1] "BiocGenerics"
    
    
    
    $Packages$gridExtra
    $Packages$gridExtra$Package
    [1] "gridExtra"
    
    $Packages$gridExtra$Version
    [1] "2.3"
    
    $Packages$gridExtra$Source
    [1] "Repository"
    
    $Packages$gridExtra$Repository
    [1] "CRAN"
    
    $Packages$gridExtra$Hash
    [1] "7d7f283939f563670a697165b2cf5560"
    
    $Packages$gridExtra$Requirements
    $Packages$gridExtra$Requirements[[1]]
    [1] "gtable"
    
    
    
    $Packages$gridGraphics
    $Packages$gridGraphics$Package
    [1] "gridGraphics"
    
    $Packages$gridGraphics$Version
    [1] "0.5-1"
    
    $Packages$gridGraphics$Source
    [1] "Repository"
    
    $Packages$gridGraphics$Repository
    [1] "CRAN"
    
    $Packages$gridGraphics$Hash
    [1] "5b79228594f02385d4df4979284879ae"
    
    $Packages$gridGraphics$Requirements
    list()
    
    
    $Packages$gt
    $Packages$gt$Package
    [1] "gt"
    
    $Packages$gt$Version
    [1] "0.4.0"
    
    $Packages$gt$Source
    [1] "Repository"
    
    $Packages$gt$Repository
    [1] "CRAN"
    
    $Packages$gt$Hash
    [1] "4905765870343c66704c9e495e4e3bc0"
    
    $Packages$gt$Requirements
    $Packages$gt$Requirements[[1]]
    [1] "base64enc"
    
    $Packages$gt$Requirements[[2]]
    [1] "bitops"
    
    $Packages$gt$Requirements[[3]]
    [1] "checkmate"
    
    $Packages$gt$Requirements[[4]]
    [1] "commonmark"
    
    $Packages$gt$Requirements[[5]]
    [1] "dplyr"
    
    $Packages$gt$Requirements[[6]]
    [1] "fs"
    
    $Packages$gt$Requirements[[7]]
    [1] "ggplot2"
    
    $Packages$gt$Requirements[[8]]
    [1] "glue"
    
    $Packages$gt$Requirements[[9]]
    [1] "htmltools"
    
    $Packages$gt$Requirements[[10]]
    [1] "magrittr"
    
    $Packages$gt$Requirements[[11]]
    [1] "rlang"
    
    $Packages$gt$Requirements[[12]]
    [1] "sass"
    
    $Packages$gt$Requirements[[13]]
    [1] "scales"
    
    $Packages$gt$Requirements[[14]]
    [1] "stringr"
    
    $Packages$gt$Requirements[[15]]
    [1] "tibble"
    
    $Packages$gt$Requirements[[16]]
    [1] "tidyselect"
    
    
    
    $Packages$gtable
    $Packages$gtable$Package
    [1] "gtable"
    
    $Packages$gtable$Version
    [1] "0.3.0"
    
    $Packages$gtable$Source
    [1] "Repository"
    
    $Packages$gtable$Repository
    [1] "CRAN"
    
    $Packages$gtable$Hash
    [1] "ac5c6baf7822ce8732b343f14c072c4d"
    
    $Packages$gtable$Requirements
    list()
    
    
    $Packages$hardhat
    $Packages$hardhat$Package
    [1] "hardhat"
    
    $Packages$hardhat$Version
    [1] "0.2.0"
    
    $Packages$hardhat$Source
    [1] "Repository"
    
    $Packages$hardhat$Repository
    [1] "CRAN"
    
    $Packages$hardhat$Hash
    [1] "726bc6915fd4533da5d49f4dfe6df1a6"
    
    $Packages$hardhat$Requirements
    $Packages$hardhat$Requirements[[1]]
    [1] "glue"
    
    $Packages$hardhat$Requirements[[2]]
    [1] "rlang"
    
    $Packages$hardhat$Requirements[[3]]
    [1] "tibble"
    
    $Packages$hardhat$Requirements[[4]]
    [1] "vctrs"
    
    
    
    $Packages$haven
    $Packages$haven$Package
    [1] "haven"
    
    $Packages$haven$Version
    [1] "2.4.3"
    
    $Packages$haven$Source
    [1] "Repository"
    
    $Packages$haven$Repository
    [1] "CRAN"
    
    $Packages$haven$Hash
    [1] "10bec8a8264f3eb59531e8c4c0303f96"
    
    $Packages$haven$Requirements
    $Packages$haven$Requirements[[1]]
    [1] "cpp11"
    
    $Packages$haven$Requirements[[2]]
    [1] "forcats"
    
    $Packages$haven$Requirements[[3]]
    [1] "hms"
    
    $Packages$haven$Requirements[[4]]
    [1] "readr"
    
    $Packages$haven$Requirements[[5]]
    [1] "rlang"
    
    $Packages$haven$Requirements[[6]]
    [1] "tibble"
    
    $Packages$haven$Requirements[[7]]
    [1] "tidyselect"
    
    $Packages$haven$Requirements[[8]]
    [1] "vctrs"
    
    
    
    $Packages$here
    $Packages$here$Package
    [1] "here"
    
    $Packages$here$Version
    [1] "1.0.1"
    
    $Packages$here$Source
    [1] "Repository"
    
    $Packages$here$Repository
    [1] "CRAN"
    
    $Packages$here$Hash
    [1] "24b224366f9c2e7534d2344d10d59211"
    
    $Packages$here$Requirements
    $Packages$here$Requirements[[1]]
    [1] "rprojroot"
    
    
    
    $Packages$highr
    $Packages$highr$Package
    [1] "highr"
    
    $Packages$highr$Version
    [1] "0.9"
    
    $Packages$highr$Source
    [1] "Repository"
    
    $Packages$highr$Repository
    [1] "CRAN"
    
    $Packages$highr$Hash
    [1] "8eb36c8125038e648e5d111c0d7b2ed4"
    
    $Packages$highr$Requirements
    $Packages$highr$Requirements[[1]]
    [1] "xfun"
    
    
    
    $Packages$hms
    $Packages$hms$Package
    [1] "hms"
    
    $Packages$hms$Version
    [1] "1.1.1"
    
    $Packages$hms$Source
    [1] "Repository"
    
    $Packages$hms$Repository
    [1] "CRAN"
    
    $Packages$hms$Hash
    [1] "5b8a2dd0fdbe2ab4f6081e6c7be6dfca"
    
    $Packages$hms$Requirements
    $Packages$hms$Requirements[[1]]
    [1] "ellipsis"
    
    $Packages$hms$Requirements[[2]]
    [1] "lifecycle"
    
    $Packages$hms$Requirements[[3]]
    [1] "pkgconfig"
    
    $Packages$hms$Requirements[[4]]
    [1] "rlang"
    
    $Packages$hms$Requirements[[5]]
    [1] "vctrs"
    
    
    
    $Packages$hoardr
    $Packages$hoardr$Package
    [1] "hoardr"
    
    $Packages$hoardr$Version
    [1] "0.5.2"
    
    $Packages$hoardr$Source
    [1] "Repository"
    
    $Packages$hoardr$Repository
    [1] "CRAN"
    
    $Packages$hoardr$Hash
    [1] "ca8a25aa079a8cb4fe44fe12b17b0eda"
    
    $Packages$hoardr$Requirements
    $Packages$hoardr$Requirements[[1]]
    [1] "R6"
    
    $Packages$hoardr$Requirements[[2]]
    [1] "digest"
    
    $Packages$hoardr$Requirements[[3]]
    [1] "rappdirs"
    
    
    
    $Packages$htmltools
    $Packages$htmltools$Package
    [1] "htmltools"
    
    $Packages$htmltools$Version
    [1] "0.5.2"
    
    $Packages$htmltools$Source
    [1] "Repository"
    
    $Packages$htmltools$Repository
    [1] "CRAN"
    
    $Packages$htmltools$Hash
    [1] "526c484233f42522278ab06fb185cb26"
    
    $Packages$htmltools$Requirements
    $Packages$htmltools$Requirements[[1]]
    [1] "base64enc"
    
    $Packages$htmltools$Requirements[[2]]
    [1] "digest"
    
    $Packages$htmltools$Requirements[[3]]
    [1] "fastmap"
    
    $Packages$htmltools$Requirements[[4]]
    [1] "rlang"
    
    
    
    $Packages$htmlwidgets
    $Packages$htmlwidgets$Package
    [1] "htmlwidgets"
    
    $Packages$htmlwidgets$Version
    [1] "1.5.4"
    
    $Packages$htmlwidgets$Source
    [1] "Repository"
    
    $Packages$htmlwidgets$Repository
    [1] "CRAN"
    
    $Packages$htmlwidgets$Hash
    [1] "76147821cd3fcd8c4b04e1ef0498e7fb"
    
    $Packages$htmlwidgets$Requirements
    $Packages$htmlwidgets$Requirements[[1]]
    [1] "htmltools"
    
    $Packages$htmlwidgets$Requirements[[2]]
    [1] "jsonlite"
    
    $Packages$htmlwidgets$Requirements[[3]]
    [1] "yaml"
    
    
    
    $Packages$httpuv
    $Packages$httpuv$Package
    [1] "httpuv"
    
    $Packages$httpuv$Version
    [1] "1.6.5"
    
    $Packages$httpuv$Source
    [1] "Repository"
    
    $Packages$httpuv$Repository
    [1] "CRAN"
    
    $Packages$httpuv$Hash
    [1] "97fe71f0a4a1c9890e6c2128afa04bc0"
    
    $Packages$httpuv$Requirements
    $Packages$httpuv$Requirements[[1]]
    [1] "R6"
    
    $Packages$httpuv$Requirements[[2]]
    [1] "Rcpp"
    
    $Packages$httpuv$Requirements[[3]]
    [1] "later"
    
    $Packages$httpuv$Requirements[[4]]
    [1] "promises"
    
    
    
    $Packages$httr
    $Packages$httr$Package
    [1] "httr"
    
    $Packages$httr$Version
    [1] "1.4.2"
    
    $Packages$httr$Source
    [1] "Repository"
    
    $Packages$httr$Repository
    [1] "CRAN"
    
    $Packages$httr$Hash
    [1] "a525aba14184fec243f9eaec62fbed43"
    
    $Packages$httr$Requirements
    $Packages$httr$Requirements[[1]]
    [1] "R6"
    
    $Packages$httr$Requirements[[2]]
    [1] "curl"
    
    $Packages$httr$Requirements[[3]]
    [1] "jsonlite"
    
    $Packages$httr$Requirements[[4]]
    [1] "mime"
    
    $Packages$httr$Requirements[[5]]
    [1] "openssl"
    
    
    
    $Packages$hunspell
    $Packages$hunspell$Package
    [1] "hunspell"
    
    $Packages$hunspell$Version
    [1] "3.0.1"
    
    $Packages$hunspell$Source
    [1] "Repository"
    
    $Packages$hunspell$Repository
    [1] "CRAN"
    
    $Packages$hunspell$Hash
    [1] "3987784c19192ad0f2261c456d936df1"
    
    $Packages$hunspell$Requirements
    $Packages$hunspell$Requirements[[1]]
    [1] "Rcpp"
    
    $Packages$hunspell$Requirements[[2]]
    [1] "digest"
    
    
    
    $Packages$hwriter
    $Packages$hwriter$Package
    [1] "hwriter"
    
    $Packages$hwriter$Version
    [1] "1.3.2"
    
    $Packages$hwriter$Source
    [1] "Repository"
    
    $Packages$hwriter$Repository
    [1] "CRAN"
    
    $Packages$hwriter$Hash
    [1] "fd39c53c2fff477e78557114b5ff873b"
    
    $Packages$hwriter$Requirements
    list()
    
    
    $Packages$ids
    $Packages$ids$Package
    [1] "ids"
    
    $Packages$ids$Version
    [1] "1.0.1"
    
    $Packages$ids$Source
    [1] "Repository"
    
    $Packages$ids$Repository
    [1] "CRAN"
    
    $Packages$ids$Hash
    [1] "99df65cfef20e525ed38c3d2577f7190"
    
    $Packages$ids$Requirements
    $Packages$ids$Requirements[[1]]
    [1] "openssl"
    
    $Packages$ids$Requirements[[2]]
    [1] "uuid"
    
    
    
    $Packages$igraph
    $Packages$igraph$Package
    [1] "igraph"
    
    $Packages$igraph$Version
    [1] "1.2.11"
    
    $Packages$igraph$Source
    [1] "Repository"
    
    $Packages$igraph$Repository
    [1] "CRAN"
    
    $Packages$igraph$Hash
    [1] "1d10cd31c2979f9c819ffe4d16b9dc2b"
    
    $Packages$igraph$Requirements
    $Packages$igraph$Requirements[[1]]
    [1] "Matrix"
    
    $Packages$igraph$Requirements[[2]]
    [1] "magrittr"
    
    $Packages$igraph$Requirements[[3]]
    [1] "pkgconfig"
    
    
    
    $Packages$infer
    $Packages$infer$Package
    [1] "infer"
    
    $Packages$infer$Version
    [1] "1.0.0"
    
    $Packages$infer$Source
    [1] "Repository"
    
    $Packages$infer$Repository
    [1] "CRAN"
    
    $Packages$infer$Hash
    [1] "74128662081cedb158fc98d6dae73c0e"
    
    $Packages$infer$Requirements
    $Packages$infer$Requirements[[1]]
    [1] "broom"
    
    $Packages$infer$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$infer$Requirements[[3]]
    [1] "generics"
    
    $Packages$infer$Requirements[[4]]
    [1] "ggplot2"
    
    $Packages$infer$Requirements[[5]]
    [1] "glue"
    
    $Packages$infer$Requirements[[6]]
    [1] "magrittr"
    
    $Packages$infer$Requirements[[7]]
    [1] "patchwork"
    
    $Packages$infer$Requirements[[8]]
    [1] "purrr"
    
    $Packages$infer$Requirements[[9]]
    [1] "rlang"
    
    $Packages$infer$Requirements[[10]]
    [1] "tibble"
    
    $Packages$infer$Requirements[[11]]
    [1] "tidyr"
    
    
    
    $Packages$ini
    $Packages$ini$Package
    [1] "ini"
    
    $Packages$ini$Version
    [1] "0.3.1"
    
    $Packages$ini$Source
    [1] "Repository"
    
    $Packages$ini$Repository
    [1] "CRAN"
    
    $Packages$ini$Hash
    [1] "6154ec2223172bce8162d4153cda21f7"
    
    $Packages$ini$Requirements
    list()
    
    
    $Packages$interactiveDisplayBase
    $Packages$interactiveDisplayBase$Package
    [1] "interactiveDisplayBase"
    
    $Packages$interactiveDisplayBase$Version
    [1] "1.32.0"
    
    $Packages$interactiveDisplayBase$Source
    [1] "Bioconductor"
    
    $Packages$interactiveDisplayBase$git_url
    [1] "https://git.bioconductor.org/packages/interactiveDisplayBase"
    
    $Packages$interactiveDisplayBase$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$interactiveDisplayBase$git_last_commit
    [1] "0f88b2a"
    
    $Packages$interactiveDisplayBase$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$interactiveDisplayBase$Hash
    [1] "011d5724a7ccdd3c3e8788701242666d"
    
    $Packages$interactiveDisplayBase$Requirements
    $Packages$interactiveDisplayBase$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$interactiveDisplayBase$Requirements[[2]]
    [1] "DT"
    
    $Packages$interactiveDisplayBase$Requirements[[3]]
    [1] "shiny"
    
    
    
    $Packages$ipred
    $Packages$ipred$Package
    [1] "ipred"
    
    $Packages$ipred$Version
    [1] "0.9-12"
    
    $Packages$ipred$Source
    [1] "Repository"
    
    $Packages$ipred$Repository
    [1] "CRAN"
    
    $Packages$ipred$Hash
    [1] "8312ebd8121ad2eca1c76441040bee5d"
    
    $Packages$ipred$Requirements
    $Packages$ipred$Requirements[[1]]
    [1] "MASS"
    
    $Packages$ipred$Requirements[[2]]
    [1] "class"
    
    $Packages$ipred$Requirements[[3]]
    [1] "nnet"
    
    $Packages$ipred$Requirements[[4]]
    [1] "prodlim"
    
    $Packages$ipred$Requirements[[5]]
    [1] "rpart"
    
    $Packages$ipred$Requirements[[6]]
    [1] "survival"
    
    
    
    $Packages$irlba
    $Packages$irlba$Package
    [1] "irlba"
    
    $Packages$irlba$Version
    [1] "2.3.5"
    
    $Packages$irlba$Source
    [1] "Repository"
    
    $Packages$irlba$Repository
    [1] "CRAN"
    
    $Packages$irlba$Hash
    [1] "066c11bb9bc75b343f3de1ecaf3b7ba2"
    
    $Packages$irlba$Requirements
    $Packages$irlba$Requirements[[1]]
    [1] "Matrix"
    
    
    
    $Packages$isoband
    $Packages$isoband$Package
    [1] "isoband"
    
    $Packages$isoband$Version
    [1] "0.2.5"
    
    $Packages$isoband$Source
    [1] "Repository"
    
    $Packages$isoband$Repository
    [1] "CRAN"
    
    $Packages$isoband$Hash
    [1] "7ab57a6de7f48a8dc84910d1eca42883"
    
    $Packages$isoband$Requirements
    list()
    
    
    $Packages$iterators
    $Packages$iterators$Package
    [1] "iterators"
    
    $Packages$iterators$Version
    [1] "1.0.14"
    
    $Packages$iterators$Source
    [1] "Repository"
    
    $Packages$iterators$Repository
    [1] "CRAN"
    
    $Packages$iterators$Hash
    [1] "8954069286b4b2b0d023d1b288dce978"
    
    $Packages$iterators$Requirements
    list()
    
    
    $Packages$janeaustenr
    $Packages$janeaustenr$Package
    [1] "janeaustenr"
    
    $Packages$janeaustenr$Version
    [1] "0.1.5"
    
    $Packages$janeaustenr$Source
    [1] "Repository"
    
    $Packages$janeaustenr$Repository
    [1] "CRAN"
    
    $Packages$janeaustenr$Hash
    [1] "8b07a4b9d0a0d97d9fe12de8af6d219e"
    
    $Packages$janeaustenr$Requirements
    list()
    
    
    $Packages$jpeg
    $Packages$jpeg$Package
    [1] "jpeg"
    
    $Packages$jpeg$Version
    [1] "0.1-9"
    
    $Packages$jpeg$Source
    [1] "Repository"
    
    $Packages$jpeg$Repository
    [1] "CRAN"
    
    $Packages$jpeg$Hash
    [1] "441ee36360a57b363f4fa3df0c364630"
    
    $Packages$jpeg$Requirements
    list()
    
    
    $Packages$jquerylib
    $Packages$jquerylib$Package
    [1] "jquerylib"
    
    $Packages$jquerylib$Version
    [1] "0.1.4"
    
    $Packages$jquerylib$Source
    [1] "Repository"
    
    $Packages$jquerylib$Repository
    [1] "CRAN"
    
    $Packages$jquerylib$Hash
    [1] "5aab57a3bd297eee1c1d862735972182"
    
    $Packages$jquerylib$Requirements
    $Packages$jquerylib$Requirements[[1]]
    [1] "htmltools"
    
    
    
    $Packages$jsonlite
    $Packages$jsonlite$Package
    [1] "jsonlite"
    
    $Packages$jsonlite$Version
    [1] "1.8.0"
    
    $Packages$jsonlite$Source
    [1] "Repository"
    
    $Packages$jsonlite$Repository
    [1] "CRAN"
    
    $Packages$jsonlite$Hash
    [1] "d07e729b27b372429d42d24d503613a0"
    
    $Packages$jsonlite$Requirements
    list()
    
    
    $Packages$kableExtra
    $Packages$kableExtra$Package
    [1] "kableExtra"
    
    $Packages$kableExtra$Version
    [1] "1.3.4"
    
    $Packages$kableExtra$Source
    [1] "Repository"
    
    $Packages$kableExtra$Repository
    [1] "CRAN"
    
    $Packages$kableExtra$Hash
    [1] "49b625e6aabe4c5f091f5850aba8ff78"
    
    $Packages$kableExtra$Requirements
    $Packages$kableExtra$Requirements[[1]]
    [1] "digest"
    
    $Packages$kableExtra$Requirements[[2]]
    [1] "glue"
    
    $Packages$kableExtra$Requirements[[3]]
    [1] "htmltools"
    
    $Packages$kableExtra$Requirements[[4]]
    [1] "knitr"
    
    $Packages$kableExtra$Requirements[[5]]
    [1] "magrittr"
    
    $Packages$kableExtra$Requirements[[6]]
    [1] "rmarkdown"
    
    $Packages$kableExtra$Requirements[[7]]
    [1] "rstudioapi"
    
    $Packages$kableExtra$Requirements[[8]]
    [1] "rvest"
    
    $Packages$kableExtra$Requirements[[9]]
    [1] "scales"
    
    $Packages$kableExtra$Requirements[[10]]
    [1] "stringr"
    
    $Packages$kableExtra$Requirements[[11]]
    [1] "svglite"
    
    $Packages$kableExtra$Requirements[[12]]
    [1] "viridisLite"
    
    $Packages$kableExtra$Requirements[[13]]
    [1] "webshot"
    
    $Packages$kableExtra$Requirements[[14]]
    [1] "xml2"
    
    
    
    $Packages$knitr
    $Packages$knitr$Package
    [1] "knitr"
    
    $Packages$knitr$Version
    [1] "1.37"
    
    $Packages$knitr$Source
    [1] "Repository"
    
    $Packages$knitr$Repository
    [1] "CRAN"
    
    $Packages$knitr$Hash
    [1] "a4ec675eb332a33fe7b7fe26f70e1f98"
    
    $Packages$knitr$Requirements
    $Packages$knitr$Requirements[[1]]
    [1] "evaluate"
    
    $Packages$knitr$Requirements[[2]]
    [1] "highr"
    
    $Packages$knitr$Requirements[[3]]
    [1] "stringr"
    
    $Packages$knitr$Requirements[[4]]
    [1] "xfun"
    
    $Packages$knitr$Requirements[[5]]
    [1] "yaml"
    
    
    
    $Packages$labeling
    $Packages$labeling$Package
    [1] "labeling"
    
    $Packages$labeling$Version
    [1] "0.4.2"
    
    $Packages$labeling$Source
    [1] "Repository"
    
    $Packages$labeling$Repository
    [1] "CRAN"
    
    $Packages$labeling$Hash
    [1] "3d5108641f47470611a32d0bdf357a72"
    
    $Packages$labeling$Requirements
    list()
    
    
    $Packages$lambda.r
    $Packages$lambda.r$Package
    [1] "lambda.r"
    
    $Packages$lambda.r$Version
    [1] "1.2.4"
    
    $Packages$lambda.r$Source
    [1] "Repository"
    
    $Packages$lambda.r$Repository
    [1] "CRAN"
    
    $Packages$lambda.r$Hash
    [1] "b1e925c4b9ffeb901bacf812cbe9a6ad"
    
    $Packages$lambda.r$Requirements
    $Packages$lambda.r$Requirements[[1]]
    [1] "formatR"
    
    
    
    $Packages$later
    $Packages$later$Package
    [1] "later"
    
    $Packages$later$Version
    [1] "1.3.0"
    
    $Packages$later$Source
    [1] "Repository"
    
    $Packages$later$Repository
    [1] "CRAN"
    
    $Packages$later$Hash
    [1] "7e7b457d7766bc47f2a5f21cc2984f8e"
    
    $Packages$later$Requirements
    $Packages$later$Requirements[[1]]
    [1] "Rcpp"
    
    $Packages$later$Requirements[[2]]
    [1] "rlang"
    
    
    
    $Packages$lattice
    $Packages$lattice$Package
    [1] "lattice"
    
    $Packages$lattice$Version
    [1] "0.20-45"
    
    $Packages$lattice$Source
    [1] "Repository"
    
    $Packages$lattice$Repository
    [1] "CRAN"
    
    $Packages$lattice$Hash
    [1] "b64cdbb2b340437c4ee047a1f4c4377b"
    
    $Packages$lattice$Requirements
    list()
    
    
    $Packages$latticeExtra
    $Packages$latticeExtra$Package
    [1] "latticeExtra"
    
    $Packages$latticeExtra$Version
    [1] "0.6-29"
    
    $Packages$latticeExtra$Source
    [1] "Repository"
    
    $Packages$latticeExtra$Repository
    [1] "CRAN"
    
    $Packages$latticeExtra$Hash
    [1] "590829599d6182cf7461787af34666ee"
    
    $Packages$latticeExtra$Requirements
    $Packages$latticeExtra$Requirements[[1]]
    [1] "RColorBrewer"
    
    $Packages$latticeExtra$Requirements[[2]]
    [1] "jpeg"
    
    $Packages$latticeExtra$Requirements[[3]]
    [1] "lattice"
    
    $Packages$latticeExtra$Requirements[[4]]
    [1] "png"
    
    
    
    $Packages$lava
    $Packages$lava$Package
    [1] "lava"
    
    $Packages$lava$Version
    [1] "1.6.10"
    
    $Packages$lava$Source
    [1] "Repository"
    
    $Packages$lava$Repository
    [1] "CRAN"
    
    $Packages$lava$Hash
    [1] "4c31a28528978d8689145f5274ce9058"
    
    $Packages$lava$Requirements
    $Packages$lava$Requirements[[1]]
    [1] "SQUAREM"
    
    $Packages$lava$Requirements[[2]]
    [1] "future.apply"
    
    $Packages$lava$Requirements[[3]]
    [1] "numDeriv"
    
    $Packages$lava$Requirements[[4]]
    [1] "progressr"
    
    $Packages$lava$Requirements[[5]]
    [1] "survival"
    
    
    
    $Packages$lazyeval
    $Packages$lazyeval$Package
    [1] "lazyeval"
    
    $Packages$lazyeval$Version
    [1] "0.2.2"
    
    $Packages$lazyeval$Source
    [1] "Repository"
    
    $Packages$lazyeval$Repository
    [1] "CRAN"
    
    $Packages$lazyeval$Hash
    [1] "d908914ae53b04d4c0c0fd72ecc35370"
    
    $Packages$lazyeval$Requirements
    list()
    
    
    $Packages$lhs
    $Packages$lhs$Package
    [1] "lhs"
    
    $Packages$lhs$Version
    [1] "1.1.4"
    
    $Packages$lhs$Source
    [1] "Repository"
    
    $Packages$lhs$Repository
    [1] "CRAN"
    
    $Packages$lhs$Hash
    [1] "452453e7fc05361e7c8cf98e3972d22c"
    
    $Packages$lhs$Requirements
    $Packages$lhs$Requirements[[1]]
    [1] "Rcpp"
    
    
    
    $Packages$lifecycle
    $Packages$lifecycle$Package
    [1] "lifecycle"
    
    $Packages$lifecycle$Version
    [1] "1.0.1"
    
    $Packages$lifecycle$Source
    [1] "Repository"
    
    $Packages$lifecycle$Repository
    [1] "CRAN"
    
    $Packages$lifecycle$Hash
    [1] "a6b6d352e3ed897373ab19d8395c98d0"
    
    $Packages$lifecycle$Requirements
    $Packages$lifecycle$Requirements[[1]]
    [1] "glue"
    
    $Packages$lifecycle$Requirements[[2]]
    [1] "rlang"
    
    
    
    $Packages$listenv
    $Packages$listenv$Package
    [1] "listenv"
    
    $Packages$listenv$Version
    [1] "0.8.0"
    
    $Packages$listenv$Source
    [1] "Repository"
    
    $Packages$listenv$Repository
    [1] "CRAN"
    
    $Packages$listenv$Hash
    [1] "0bde42ee282efb18c7c4e63822f5b4f7"
    
    $Packages$listenv$Requirements
    list()
    
    
    $Packages$lubridate
    $Packages$lubridate$Package
    [1] "lubridate"
    
    $Packages$lubridate$Version
    [1] "1.8.0"
    
    $Packages$lubridate$Source
    [1] "Repository"
    
    $Packages$lubridate$Repository
    [1] "CRAN"
    
    $Packages$lubridate$Hash
    [1] "2ff5eedb6ee38fb1b81205c73be1be5a"
    
    $Packages$lubridate$Requirements
    $Packages$lubridate$Requirements[[1]]
    [1] "cpp11"
    
    $Packages$lubridate$Requirements[[2]]
    [1] "generics"
    
    
    
    $Packages$magrittr
    $Packages$magrittr$Package
    [1] "magrittr"
    
    $Packages$magrittr$Version
    [1] "2.0.2"
    
    $Packages$magrittr$Source
    [1] "Repository"
    
    $Packages$magrittr$Repository
    [1] "CRAN"
    
    $Packages$magrittr$Hash
    [1] "cdc87ecd81934679d1557633d8e1fe51"
    
    $Packages$magrittr$Requirements
    list()
    
    
    $Packages$matrixStats
    $Packages$matrixStats$Package
    [1] "matrixStats"
    
    $Packages$matrixStats$Version
    [1] "0.61.0"
    
    $Packages$matrixStats$Source
    [1] "Repository"
    
    $Packages$matrixStats$Repository
    [1] "CRAN"
    
    $Packages$matrixStats$Hash
    [1] "b8e6221fc11247b12ab1b055a6f66c27"
    
    $Packages$matrixStats$Requirements
    list()
    
    
    $Packages$memoise
    $Packages$memoise$Package
    [1] "memoise"
    
    $Packages$memoise$Version
    [1] "2.0.1"
    
    $Packages$memoise$Source
    [1] "Repository"
    
    $Packages$memoise$Repository
    [1] "CRAN"
    
    $Packages$memoise$Hash
    [1] "e2817ccf4a065c5d9d7f2cfbe7c1d78c"
    
    $Packages$memoise$Requirements
    $Packages$memoise$Requirements[[1]]
    [1] "cachem"
    
    $Packages$memoise$Requirements[[2]]
    [1] "rlang"
    
    
    
    $Packages$mgcv
    $Packages$mgcv$Package
    [1] "mgcv"
    
    $Packages$mgcv$Version
    [1] "1.8-39"
    
    $Packages$mgcv$Source
    [1] "Repository"
    
    $Packages$mgcv$Repository
    [1] "CRAN"
    
    $Packages$mgcv$Hash
    [1] "055265005c238024e306fe0b600c89ff"
    
    $Packages$mgcv$Requirements
    $Packages$mgcv$Requirements[[1]]
    [1] "Matrix"
    
    $Packages$mgcv$Requirements[[2]]
    [1] "nlme"
    
    
    
    $Packages$mia
    $Packages$mia$Package
    [1] "mia"
    
    $Packages$mia$Version
    [1] "1.2.7"
    
    $Packages$mia$Source
    [1] "Bioconductor"
    
    $Packages$mia$git_url
    [1] "https://git.bioconductor.org/packages/mia"
    
    $Packages$mia$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$mia$git_last_commit
    [1] "217dc43"
    
    $Packages$mia$git_last_commit_date
    [1] "2022-02-06"
    
    $Packages$mia$Hash
    [1] "824cd6adeffd1bad73f2b6ed19f3967c"
    
    $Packages$mia$Requirements
    $Packages$mia$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$mia$Requirements[[2]]
    [1] "BiocParallel"
    
    $Packages$mia$Requirements[[3]]
    [1] "Biostrings"
    
    $Packages$mia$Requirements[[4]]
    [1] "DECIPHER"
    
    $Packages$mia$Requirements[[5]]
    [1] "DelayedArray"
    
    $Packages$mia$Requirements[[6]]
    [1] "DelayedMatrixStats"
    
    $Packages$mia$Requirements[[7]]
    [1] "DirichletMultinomial"
    
    $Packages$mia$Requirements[[8]]
    [1] "IRanges"
    
    $Packages$mia$Requirements[[9]]
    [1] "MASS"
    
    $Packages$mia$Requirements[[10]]
    [1] "MultiAssayExperiment"
    
    $Packages$mia$Requirements[[11]]
    [1] "S4Vectors"
    
    $Packages$mia$Requirements[[12]]
    [1] "SingleCellExperiment"
    
    $Packages$mia$Requirements[[13]]
    [1] "SummarizedExperiment"
    
    $Packages$mia$Requirements[[14]]
    [1] "TreeSummarizedExperiment"
    
    $Packages$mia$Requirements[[15]]
    [1] "ape"
    
    $Packages$mia$Requirements[[16]]
    [1] "decontam"
    
    $Packages$mia$Requirements[[17]]
    [1] "dplyr"
    
    $Packages$mia$Requirements[[18]]
    [1] "rlang"
    
    $Packages$mia$Requirements[[19]]
    [1] "scater"
    
    $Packages$mia$Requirements[[20]]
    [1] "scuttle"
    
    $Packages$mia$Requirements[[21]]
    [1] "tibble"
    
    $Packages$mia$Requirements[[22]]
    [1] "tidyr"
    
    $Packages$mia$Requirements[[23]]
    [1] "vegan"
    
    
    
    $Packages$microbiome
    $Packages$microbiome$Package
    [1] "microbiome"
    
    $Packages$microbiome$Version
    [1] "1.16.0"
    
    $Packages$microbiome$Source
    [1] "Bioconductor"
    
    $Packages$microbiome$git_url
    [1] "https://git.bioconductor.org/packages/microbiome"
    
    $Packages$microbiome$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$microbiome$git_last_commit
    [1] "a7b74b7"
    
    $Packages$microbiome$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$microbiome$Hash
    [1] "4fd8c7b66d1b18224348a02478c81164"
    
    $Packages$microbiome$Requirements
    $Packages$microbiome$Requirements[[1]]
    [1] "Rtsne"
    
    $Packages$microbiome$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$microbiome$Requirements[[3]]
    [1] "ggplot2"
    
    $Packages$microbiome$Requirements[[4]]
    [1] "phyloseq"
    
    $Packages$microbiome$Requirements[[5]]
    [1] "reshape2"
    
    $Packages$microbiome$Requirements[[6]]
    [1] "scales"
    
    $Packages$microbiome$Requirements[[7]]
    [1] "tibble"
    
    $Packages$microbiome$Requirements[[8]]
    [1] "tidyr"
    
    $Packages$microbiome$Requirements[[9]]
    [1] "vegan"
    
    
    
    $Packages$mime
    $Packages$mime$Package
    [1] "mime"
    
    $Packages$mime$Version
    [1] "0.12"
    
    $Packages$mime$Source
    [1] "Repository"
    
    $Packages$mime$Repository
    [1] "CRAN"
    
    $Packages$mime$Hash
    [1] "18e9c28c1d3ca1560ce30658b22ce104"
    
    $Packages$mime$Requirements
    list()
    
    
    $Packages$mlr
    $Packages$mlr$Package
    [1] "mlr"
    
    $Packages$mlr$Version
    [1] "2.19.0"
    
    $Packages$mlr$Source
    [1] "Repository"
    
    $Packages$mlr$Repository
    [1] "CRAN"
    
    $Packages$mlr$Hash
    [1] "cf4429d0d8bb3c06fe8b07eb85ffead9"
    
    $Packages$mlr$Requirements
    $Packages$mlr$Requirements[[1]]
    [1] "BBmisc"
    
    $Packages$mlr$Requirements[[2]]
    [1] "ParamHelpers"
    
    $Packages$mlr$Requirements[[3]]
    [1] "XML"
    
    $Packages$mlr$Requirements[[4]]
    [1] "backports"
    
    $Packages$mlr$Requirements[[5]]
    [1] "checkmate"
    
    $Packages$mlr$Requirements[[6]]
    [1] "data.table"
    
    $Packages$mlr$Requirements[[7]]
    [1] "ggplot2"
    
    $Packages$mlr$Requirements[[8]]
    [1] "parallelMap"
    
    $Packages$mlr$Requirements[[9]]
    [1] "stringi"
    
    $Packages$mlr$Requirements[[10]]
    [1] "survival"
    
    
    
    $Packages$modeldata
    $Packages$modeldata$Package
    [1] "modeldata"
    
    $Packages$modeldata$Version
    [1] "0.1.1"
    
    $Packages$modeldata$Source
    [1] "Repository"
    
    $Packages$modeldata$Repository
    [1] "CRAN"
    
    $Packages$modeldata$Hash
    [1] "30951e1c2022c6a8e1809d4224c0a5ab"
    
    $Packages$modeldata$Requirements
    list()
    
    
    $Packages$modelr
    $Packages$modelr$Package
    [1] "modelr"
    
    $Packages$modelr$Version
    [1] "0.1.8"
    
    $Packages$modelr$Source
    [1] "Repository"
    
    $Packages$modelr$Repository
    [1] "CRAN"
    
    $Packages$modelr$Hash
    [1] "9fd59716311ee82cba83dc2826fc5577"
    
    $Packages$modelr$Requirements
    $Packages$modelr$Requirements[[1]]
    [1] "broom"
    
    $Packages$modelr$Requirements[[2]]
    [1] "magrittr"
    
    $Packages$modelr$Requirements[[3]]
    [1] "purrr"
    
    $Packages$modelr$Requirements[[4]]
    [1] "rlang"
    
    $Packages$modelr$Requirements[[5]]
    [1] "tibble"
    
    $Packages$modelr$Requirements[[6]]
    [1] "tidyr"
    
    $Packages$modelr$Requirements[[7]]
    [1] "tidyselect"
    
    $Packages$modelr$Requirements[[8]]
    [1] "vctrs"
    
    
    
    $Packages$multtest
    $Packages$multtest$Package
    [1] "multtest"
    
    $Packages$multtest$Version
    [1] "2.50.0"
    
    $Packages$multtest$Source
    [1] "Bioconductor"
    
    $Packages$multtest$git_url
    [1] "https://git.bioconductor.org/packages/multtest"
    
    $Packages$multtest$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$multtest$git_last_commit
    [1] "1de9664"
    
    $Packages$multtest$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$multtest$Hash
    [1] "096b9b9a42d79a82431eb5d9ed9f0da3"
    
    $Packages$multtest$Requirements
    $Packages$multtest$Requirements[[1]]
    [1] "Biobase"
    
    $Packages$multtest$Requirements[[2]]
    [1] "BiocGenerics"
    
    $Packages$multtest$Requirements[[3]]
    [1] "MASS"
    
    $Packages$multtest$Requirements[[4]]
    [1] "survival"
    
    
    
    $Packages$munsell
    $Packages$munsell$Package
    [1] "munsell"
    
    $Packages$munsell$Version
    [1] "0.5.0"
    
    $Packages$munsell$Source
    [1] "Repository"
    
    $Packages$munsell$Repository
    [1] "CRAN"
    
    $Packages$munsell$Hash
    [1] "6dfe8bf774944bd5595785e3229d8771"
    
    $Packages$munsell$Requirements
    $Packages$munsell$Requirements[[1]]
    [1] "colorspace"
    
    
    
    $Packages$naturalsort
    $Packages$naturalsort$Package
    [1] "naturalsort"
    
    $Packages$naturalsort$Version
    [1] "0.1.3"
    
    $Packages$naturalsort$Source
    [1] "Repository"
    
    $Packages$naturalsort$Repository
    [1] "CRAN"
    
    $Packages$naturalsort$Hash
    [1] "2e48674f7d0ab0a3ace4eba2c1145001"
    
    $Packages$naturalsort$Requirements
    list()
    
    
    $Packages$nlme
    $Packages$nlme$Package
    [1] "nlme"
    
    $Packages$nlme$Version
    [1] "3.1-155"
    
    $Packages$nlme$Source
    [1] "Repository"
    
    $Packages$nlme$Repository
    [1] "CRAN"
    
    $Packages$nlme$Hash
    [1] "74ad940dccc9e977189a5afe5fcdb7ba"
    
    $Packages$nlme$Requirements
    $Packages$nlme$Requirements[[1]]
    [1] "lattice"
    
    
    
    $Packages$nloptr
    $Packages$nloptr$Package
    [1] "nloptr"
    
    $Packages$nloptr$Version
    [1] "2.0.0"
    
    $Packages$nloptr$Source
    [1] "Repository"
    
    $Packages$nloptr$Repository
    [1] "CRAN"
    
    $Packages$nloptr$Hash
    [1] "6e45c045fea34a9d0d1ceaa6fb7c4e91"
    
    $Packages$nloptr$Requirements
    $Packages$nloptr$Requirements[[1]]
    [1] "testthat"
    
    
    
    $Packages$nnet
    $Packages$nnet$Package
    [1] "nnet"
    
    $Packages$nnet$Version
    [1] "7.3-17"
    
    $Packages$nnet$Source
    [1] "Repository"
    
    $Packages$nnet$Repository
    [1] "CRAN"
    
    $Packages$nnet$Hash
    [1] "cb1d8d9f300a7e536b89c8a88c53f610"
    
    $Packages$nnet$Requirements
    list()
    
    
    $Packages$numDeriv
    $Packages$numDeriv$Package
    [1] "numDeriv"
    
    $Packages$numDeriv$Version
    [1] "2016.8-1.1"
    
    $Packages$numDeriv$Source
    [1] "Repository"
    
    $Packages$numDeriv$Repository
    [1] "CRAN"
    
    $Packages$numDeriv$Hash
    [1] "df58958f293b166e4ab885ebcad90e02"
    
    $Packages$numDeriv$Requirements
    list()
    
    
    $Packages$ontologyIndex
    $Packages$ontologyIndex$Package
    [1] "ontologyIndex"
    
    $Packages$ontologyIndex$Version
    [1] "2.7"
    
    $Packages$ontologyIndex$Source
    [1] "Repository"
    
    $Packages$ontologyIndex$Repository
    [1] "CRAN"
    
    $Packages$ontologyIndex$Hash
    [1] "ca4b9b931c659a989d3fc8e97b1b4525"
    
    $Packages$ontologyIndex$Requirements
    list()
    
    
    $Packages$openssl
    $Packages$openssl$Package
    [1] "openssl"
    
    $Packages$openssl$Version
    [1] "2.0.0"
    
    $Packages$openssl$Source
    [1] "Repository"
    
    $Packages$openssl$Repository
    [1] "CRAN"
    
    $Packages$openssl$Hash
    [1] "cf4329aac12c2c44089974559c18e446"
    
    $Packages$openssl$Requirements
    $Packages$openssl$Requirements[[1]]
    [1] "askpass"
    
    
    
    $Packages$optparse
    $Packages$optparse$Package
    [1] "optparse"
    
    $Packages$optparse$Version
    [1] "1.7.1"
    
    $Packages$optparse$Source
    [1] "Repository"
    
    $Packages$optparse$Repository
    [1] "CRAN"
    
    $Packages$optparse$Hash
    [1] "2490600671344e847c37a7f75ee458c0"
    
    $Packages$optparse$Requirements
    $Packages$optparse$Requirements[[1]]
    [1] "getopt"
    
    
    
    $Packages$pROC
    $Packages$pROC$Package
    [1] "pROC"
    
    $Packages$pROC$Version
    [1] "1.18.0"
    
    $Packages$pROC$Source
    [1] "Repository"
    
    $Packages$pROC$Repository
    [1] "CRAN"
    
    $Packages$pROC$Hash
    [1] "417fd0d40479932c19faf2747817c473"
    
    $Packages$pROC$Requirements
    $Packages$pROC$Requirements[[1]]
    [1] "Rcpp"
    
    $Packages$pROC$Requirements[[2]]
    [1] "plyr"
    
    
    
    $Packages$parallelMap
    $Packages$parallelMap$Package
    [1] "parallelMap"
    
    $Packages$parallelMap$Version
    [1] "1.5.1"
    
    $Packages$parallelMap$Source
    [1] "Repository"
    
    $Packages$parallelMap$Repository
    [1] "CRAN"
    
    $Packages$parallelMap$Hash
    [1] "3f3c911183a84966c793987f5d73f3e6"
    
    $Packages$parallelMap$Requirements
    $Packages$parallelMap$Requirements[[1]]
    [1] "BBmisc"
    
    $Packages$parallelMap$Requirements[[2]]
    [1] "checkmate"
    
    
    
    $Packages$parallelly
    $Packages$parallelly$Package
    [1] "parallelly"
    
    $Packages$parallelly$Version
    [1] "1.30.0"
    
    $Packages$parallelly$Source
    [1] "Repository"
    
    $Packages$parallelly$Repository
    [1] "CRAN"
    
    $Packages$parallelly$Hash
    [1] "67db13907a9cea89c118cf82d448799f"
    
    $Packages$parallelly$Requirements
    list()
    
    
    $Packages$parsnip
    $Packages$parsnip$Package
    [1] "parsnip"
    
    $Packages$parsnip$Version
    [1] "0.2.0.9000"
    
    $Packages$parsnip$Source
    [1] "GitHub"
    
    $Packages$parsnip$RemoteType
    [1] "github"
    
    $Packages$parsnip$RemoteUsername
    [1] "tidymodels"
    
    $Packages$parsnip$RemoteRepo
    [1] "parsnip"
    
    $Packages$parsnip$RemoteRef
    [1] "HEAD"
    
    $Packages$parsnip$RemoteSha
    [1] "6029887abf7be2fb0b7d9cebe94956627bbc3308"
    
    $Packages$parsnip$RemoteHost
    [1] "api.github.com"
    
    $Packages$parsnip$Hash
    [1] "a39aaae9e0baa971d80375c4bcbe845b"
    
    $Packages$parsnip$Requirements
    $Packages$parsnip$Requirements[[1]]
    [1] "cli"
    
    $Packages$parsnip$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$parsnip$Requirements[[3]]
    [1] "generics"
    
    $Packages$parsnip$Requirements[[4]]
    [1] "ggplot2"
    
    $Packages$parsnip$Requirements[[5]]
    [1] "globals"
    
    $Packages$parsnip$Requirements[[6]]
    [1] "glue"
    
    $Packages$parsnip$Requirements[[7]]
    [1] "hardhat"
    
    $Packages$parsnip$Requirements[[8]]
    [1] "lifecycle"
    
    $Packages$parsnip$Requirements[[9]]
    [1] "magrittr"
    
    $Packages$parsnip$Requirements[[10]]
    [1] "prettyunits"
    
    $Packages$parsnip$Requirements[[11]]
    [1] "purrr"
    
    $Packages$parsnip$Requirements[[12]]
    [1] "rlang"
    
    $Packages$parsnip$Requirements[[13]]
    [1] "tibble"
    
    $Packages$parsnip$Requirements[[14]]
    [1] "tidyr"
    
    $Packages$parsnip$Requirements[[15]]
    [1] "vctrs"
    
    $Packages$parsnip$Requirements[[16]]
    [1] "withr"
    
    
    
    $Packages$patchwork
    $Packages$patchwork$Package
    [1] "patchwork"
    
    $Packages$patchwork$Version
    [1] "1.1.1"
    
    $Packages$patchwork$Source
    [1] "Repository"
    
    $Packages$patchwork$Repository
    [1] "CRAN"
    
    $Packages$patchwork$Hash
    [1] "c446b30cb33ec125ff02588b60660ccb"
    
    $Packages$patchwork$Requirements
    $Packages$patchwork$Requirements[[1]]
    [1] "ggplot2"
    
    $Packages$patchwork$Requirements[[2]]
    [1] "gtable"
    
    
    
    $Packages$permute
    $Packages$permute$Package
    [1] "permute"
    
    $Packages$permute$Version
    [1] "0.9-7"
    
    $Packages$permute$Source
    [1] "Repository"
    
    $Packages$permute$Repository
    [1] "CRAN"
    
    $Packages$permute$Hash
    [1] "abf0ca85c1c752e0d04f46334e635046"
    
    $Packages$permute$Requirements
    list()
    
    
    $Packages$phyloseq
    $Packages$phyloseq$Package
    [1] "phyloseq"
    
    $Packages$phyloseq$Version
    [1] "1.38.0"
    
    $Packages$phyloseq$Source
    [1] "Bioconductor"
    
    $Packages$phyloseq$git_url
    [1] "https://git.bioconductor.org/packages/phyloseq"
    
    $Packages$phyloseq$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$phyloseq$git_last_commit
    [1] "1e2409a"
    
    $Packages$phyloseq$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$phyloseq$Hash
    [1] "cf198d6103ff6270c347f475ce325c14"
    
    $Packages$phyloseq$Requirements
    $Packages$phyloseq$Requirements[[1]]
    [1] "Biobase"
    
    $Packages$phyloseq$Requirements[[2]]
    [1] "BiocGenerics"
    
    $Packages$phyloseq$Requirements[[3]]
    [1] "Biostrings"
    
    $Packages$phyloseq$Requirements[[4]]
    [1] "ade4"
    
    $Packages$phyloseq$Requirements[[5]]
    [1] "ape"
    
    $Packages$phyloseq$Requirements[[6]]
    [1] "biomformat"
    
    $Packages$phyloseq$Requirements[[7]]
    [1] "cluster"
    
    $Packages$phyloseq$Requirements[[8]]
    [1] "data.table"
    
    $Packages$phyloseq$Requirements[[9]]
    [1] "foreach"
    
    $Packages$phyloseq$Requirements[[10]]
    [1] "ggplot2"
    
    $Packages$phyloseq$Requirements[[11]]
    [1] "igraph"
    
    $Packages$phyloseq$Requirements[[12]]
    [1] "multtest"
    
    $Packages$phyloseq$Requirements[[13]]
    [1] "plyr"
    
    $Packages$phyloseq$Requirements[[14]]
    [1] "reshape2"
    
    $Packages$phyloseq$Requirements[[15]]
    [1] "scales"
    
    $Packages$phyloseq$Requirements[[16]]
    [1] "vegan"
    
    
    
    $Packages$piggyback
    $Packages$piggyback$Package
    [1] "piggyback"
    
    $Packages$piggyback$Version
    [1] "0.1.1"
    
    $Packages$piggyback$Source
    [1] "Repository"
    
    $Packages$piggyback$Repository
    [1] "CRAN"
    
    $Packages$piggyback$Hash
    [1] "0cd602e7eb4288f3e78fb4eba0507ead"
    
    $Packages$piggyback$Requirements
    $Packages$piggyback$Requirements[[1]]
    [1] "clisymbols"
    
    $Packages$piggyback$Requirements[[2]]
    [1] "crayon"
    
    $Packages$piggyback$Requirements[[3]]
    [1] "fs"
    
    $Packages$piggyback$Requirements[[4]]
    [1] "gh"
    
    $Packages$piggyback$Requirements[[5]]
    [1] "httr"
    
    $Packages$piggyback$Requirements[[6]]
    [1] "jsonlite"
    
    $Packages$piggyback$Requirements[[7]]
    [1] "lubridate"
    
    $Packages$piggyback$Requirements[[8]]
    [1] "memoise"
    
    
    
    $Packages$pillar
    $Packages$pillar$Package
    [1] "pillar"
    
    $Packages$pillar$Version
    [1] "1.7.0"
    
    $Packages$pillar$Source
    [1] "Repository"
    
    $Packages$pillar$Repository
    [1] "CRAN"
    
    $Packages$pillar$Hash
    [1] "51dfc97e1b7069e9f7e6f83f3589c22e"
    
    $Packages$pillar$Requirements
    $Packages$pillar$Requirements[[1]]
    [1] "cli"
    
    $Packages$pillar$Requirements[[2]]
    [1] "crayon"
    
    $Packages$pillar$Requirements[[3]]
    [1] "ellipsis"
    
    $Packages$pillar$Requirements[[4]]
    [1] "fansi"
    
    $Packages$pillar$Requirements[[5]]
    [1] "glue"
    
    $Packages$pillar$Requirements[[6]]
    [1] "lifecycle"
    
    $Packages$pillar$Requirements[[7]]
    [1] "rlang"
    
    $Packages$pillar$Requirements[[8]]
    [1] "utf8"
    
    $Packages$pillar$Requirements[[9]]
    [1] "vctrs"
    
    
    
    $Packages$pixmap
    $Packages$pixmap$Package
    [1] "pixmap"
    
    $Packages$pixmap$Version
    [1] "0.4-12"
    
    $Packages$pixmap$Source
    [1] "Repository"
    
    $Packages$pixmap$Repository
    [1] "CRAN"
    
    $Packages$pixmap$Hash
    [1] "ff0e0b97265ced0db4d13d0005334d78"
    
    $Packages$pixmap$Requirements
    list()
    
    
    $Packages$pkgconfig
    $Packages$pkgconfig$Package
    [1] "pkgconfig"
    
    $Packages$pkgconfig$Version
    [1] "2.0.3"
    
    $Packages$pkgconfig$Source
    [1] "Repository"
    
    $Packages$pkgconfig$Repository
    [1] "CRAN"
    
    $Packages$pkgconfig$Hash
    [1] "01f28d4278f15c76cddbea05899c5d6f"
    
    $Packages$pkgconfig$Requirements
    list()
    
    
    $Packages$pkgload
    $Packages$pkgload$Package
    [1] "pkgload"
    
    $Packages$pkgload$Version
    [1] "1.2.4"
    
    $Packages$pkgload$Source
    [1] "Repository"
    
    $Packages$pkgload$Repository
    [1] "CRAN"
    
    $Packages$pkgload$Hash
    [1] "7533cd805940821bf23eaf3c8d4c1735"
    
    $Packages$pkgload$Requirements
    $Packages$pkgload$Requirements[[1]]
    [1] "cli"
    
    $Packages$pkgload$Requirements[[2]]
    [1] "crayon"
    
    $Packages$pkgload$Requirements[[3]]
    [1] "desc"
    
    $Packages$pkgload$Requirements[[4]]
    [1] "rlang"
    
    $Packages$pkgload$Requirements[[5]]
    [1] "rprojroot"
    
    $Packages$pkgload$Requirements[[6]]
    [1] "rstudioapi"
    
    $Packages$pkgload$Requirements[[7]]
    [1] "withr"
    
    
    
    $Packages$plogr
    $Packages$plogr$Package
    [1] "plogr"
    
    $Packages$plogr$Version
    [1] "0.2.0"
    
    $Packages$plogr$Source
    [1] "Repository"
    
    $Packages$plogr$Repository
    [1] "CRAN"
    
    $Packages$plogr$Hash
    [1] "09eb987710984fc2905c7129c7d85e65"
    
    $Packages$plogr$Requirements
    list()
    
    
    $Packages$plyr
    $Packages$plyr$Package
    [1] "plyr"
    
    $Packages$plyr$Version
    [1] "1.8.6"
    
    $Packages$plyr$Source
    [1] "Repository"
    
    $Packages$plyr$Repository
    [1] "CRAN"
    
    $Packages$plyr$Hash
    [1] "ec0e5ab4e5f851f6ef32cd1d1984957f"
    
    $Packages$plyr$Requirements
    $Packages$plyr$Requirements[[1]]
    [1] "Rcpp"
    
    
    
    $Packages$png
    $Packages$png$Package
    [1] "png"
    
    $Packages$png$Version
    [1] "0.1-7"
    
    $Packages$png$Source
    [1] "Repository"
    
    $Packages$png$Repository
    [1] "CRAN"
    
    $Packages$png$Hash
    [1] "03b7076c234cb3331288919983326c55"
    
    $Packages$png$Requirements
    list()
    
    
    $Packages$praise
    $Packages$praise$Package
    [1] "praise"
    
    $Packages$praise$Version
    [1] "1.0.0"
    
    $Packages$praise$Source
    [1] "Repository"
    
    $Packages$praise$Repository
    [1] "CRAN"
    
    $Packages$praise$Hash
    [1] "a555924add98c99d2f411e37e7d25e9f"
    
    $Packages$praise$Requirements
    list()
    
    
    $Packages$prettyunits
    $Packages$prettyunits$Package
    [1] "prettyunits"
    
    $Packages$prettyunits$Version
    [1] "1.1.1"
    
    $Packages$prettyunits$Source
    [1] "Repository"
    
    $Packages$prettyunits$Repository
    [1] "CRAN"
    
    $Packages$prettyunits$Hash
    [1] "95ef9167b75dde9d2ccc3c7528393e7e"
    
    $Packages$prettyunits$Requirements
    list()
    
    
    $Packages$processx
    $Packages$processx$Package
    [1] "processx"
    
    $Packages$processx$Version
    [1] "3.5.2"
    
    $Packages$processx$Source
    [1] "Repository"
    
    $Packages$processx$Repository
    [1] "CRAN"
    
    $Packages$processx$Hash
    [1] "0cbca2bc4d16525d009c4dbba156b37c"
    
    $Packages$processx$Requirements
    $Packages$processx$Requirements[[1]]
    [1] "R6"
    
    $Packages$processx$Requirements[[2]]
    [1] "ps"
    
    
    
    $Packages$prodlim
    $Packages$prodlim$Package
    [1] "prodlim"
    
    $Packages$prodlim$Version
    [1] "2019.11.13"
    
    $Packages$prodlim$Source
    [1] "Repository"
    
    $Packages$prodlim$Repository
    [1] "CRAN"
    
    $Packages$prodlim$Hash
    [1] "c243bf70db3a6631a0c8783152fb7db9"
    
    $Packages$prodlim$Requirements
    $Packages$prodlim$Requirements[[1]]
    [1] "KernSmooth"
    
    $Packages$prodlim$Requirements[[2]]
    [1] "Rcpp"
    
    $Packages$prodlim$Requirements[[3]]
    [1] "lava"
    
    $Packages$prodlim$Requirements[[4]]
    [1] "survival"
    
    
    
    $Packages$progress
    $Packages$progress$Package
    [1] "progress"
    
    $Packages$progress$Version
    [1] "1.2.2"
    
    $Packages$progress$Source
    [1] "Repository"
    
    $Packages$progress$Repository
    [1] "CRAN"
    
    $Packages$progress$Hash
    [1] "14dc9f7a3c91ebb14ec5bb9208a07061"
    
    $Packages$progress$Requirements
    $Packages$progress$Requirements[[1]]
    [1] "R6"
    
    $Packages$progress$Requirements[[2]]
    [1] "crayon"
    
    $Packages$progress$Requirements[[3]]
    [1] "hms"
    
    $Packages$progress$Requirements[[4]]
    [1] "prettyunits"
    
    
    
    $Packages$progressr
    $Packages$progressr$Package
    [1] "progressr"
    
    $Packages$progressr$Version
    [1] "0.10.0"
    
    $Packages$progressr$Source
    [1] "Repository"
    
    $Packages$progressr$Repository
    [1] "CRAN"
    
    $Packages$progressr$Hash
    [1] "7df448f5ae46ab0a1af0fb619349d3fd"
    
    $Packages$progressr$Requirements
    $Packages$progressr$Requirements[[1]]
    [1] "digest"
    
    
    
    $Packages$promises
    $Packages$promises$Package
    [1] "promises"
    
    $Packages$promises$Version
    [1] "1.2.0.1"
    
    $Packages$promises$Source
    [1] "Repository"
    
    $Packages$promises$Repository
    [1] "CRAN"
    
    $Packages$promises$Hash
    [1] "4ab2c43adb4d4699cf3690acd378d75d"
    
    $Packages$promises$Requirements
    $Packages$promises$Requirements[[1]]
    [1] "R6"
    
    $Packages$promises$Requirements[[2]]
    [1] "Rcpp"
    
    $Packages$promises$Requirements[[3]]
    [1] "later"
    
    $Packages$promises$Requirements[[4]]
    [1] "magrittr"
    
    $Packages$promises$Requirements[[5]]
    [1] "rlang"
    
    
    
    $Packages$ps
    $Packages$ps$Package
    [1] "ps"
    
    $Packages$ps$Version
    [1] "1.6.0"
    
    $Packages$ps$Source
    [1] "Repository"
    
    $Packages$ps$Repository
    [1] "CRAN"
    
    $Packages$ps$Hash
    [1] "32620e2001c1dce1af49c49dccbb9420"
    
    $Packages$ps$Requirements
    list()
    
    
    $Packages$purrr
    $Packages$purrr$Package
    [1] "purrr"
    
    $Packages$purrr$Version
    [1] "0.3.4"
    
    $Packages$purrr$Source
    [1] "Repository"
    
    $Packages$purrr$Repository
    [1] "CRAN"
    
    $Packages$purrr$Hash
    [1] "97def703420c8ab10d8f0e6c72101e02"
    
    $Packages$purrr$Requirements
    $Packages$purrr$Requirements[[1]]
    [1] "magrittr"
    
    $Packages$purrr$Requirements[[2]]
    [1] "rlang"
    
    
    
    $Packages$qs
    $Packages$qs$Package
    [1] "qs"
    
    $Packages$qs$Version
    [1] "0.25.3"
    
    $Packages$qs$Source
    [1] "Repository"
    
    $Packages$qs$Repository
    [1] "CRAN"
    
    $Packages$qs$Hash
    [1] "2a135bd3fcd62d0caf73cf1bc9270d8d"
    
    $Packages$qs$Requirements
    $Packages$qs$Requirements[[1]]
    [1] "RApiSerialize"
    
    $Packages$qs$Requirements[[2]]
    [1] "Rcpp"
    
    $Packages$qs$Requirements[[3]]
    [1] "stringfish"
    
    
    
    $Packages$randomForest
    $Packages$randomForest$Package
    [1] "randomForest"
    
    $Packages$randomForest$Version
    [1] "4.7-1"
    
    $Packages$randomForest$Source
    [1] "Repository"
    
    $Packages$randomForest$Repository
    [1] "CRAN"
    
    $Packages$randomForest$Hash
    [1] "dd4f861c0353a869964024886c93b7d5"
    
    $Packages$randomForest$Requirements
    list()
    
    
    $Packages$rappdirs
    $Packages$rappdirs$Package
    [1] "rappdirs"
    
    $Packages$rappdirs$Version
    [1] "0.3.3"
    
    $Packages$rappdirs$Source
    [1] "Repository"
    
    $Packages$rappdirs$Repository
    [1] "CRAN"
    
    $Packages$rappdirs$Hash
    [1] "5e3c5dc0b071b21fa128676560dbe94d"
    
    $Packages$rappdirs$Requirements
    list()
    
    
    $Packages$rbibutils
    $Packages$rbibutils$Package
    [1] "rbibutils"
    
    $Packages$rbibutils$Version
    [1] "2.2.7"
    
    $Packages$rbibutils$Source
    [1] "Repository"
    
    $Packages$rbibutils$Repository
    [1] "CRAN"
    
    $Packages$rbibutils$Hash
    [1] "69f4a840b0c2aa99f5ce473fd7fa21bd"
    
    $Packages$rbibutils$Requirements
    list()
    
    
    $Packages$readr
    $Packages$readr$Package
    [1] "readr"
    
    $Packages$readr$Version
    [1] "2.1.2"
    
    $Packages$readr$Source
    [1] "Repository"
    
    $Packages$readr$Repository
    [1] "CRAN"
    
    $Packages$readr$Hash
    [1] "9c59de1357dc209868b5feb5c9f0fe2f"
    
    $Packages$readr$Requirements
    $Packages$readr$Requirements[[1]]
    [1] "R6"
    
    $Packages$readr$Requirements[[2]]
    [1] "cli"
    
    $Packages$readr$Requirements[[3]]
    [1] "clipr"
    
    $Packages$readr$Requirements[[4]]
    [1] "cpp11"
    
    $Packages$readr$Requirements[[5]]
    [1] "crayon"
    
    $Packages$readr$Requirements[[6]]
    [1] "hms"
    
    $Packages$readr$Requirements[[7]]
    [1] "lifecycle"
    
    $Packages$readr$Requirements[[8]]
    [1] "rlang"
    
    $Packages$readr$Requirements[[9]]
    [1] "tibble"
    
    $Packages$readr$Requirements[[10]]
    [1] "tzdb"
    
    $Packages$readr$Requirements[[11]]
    [1] "vroom"
    
    
    
    $Packages$readxl
    $Packages$readxl$Package
    [1] "readxl"
    
    $Packages$readxl$Version
    [1] "1.3.1"
    
    $Packages$readxl$Source
    [1] "Repository"
    
    $Packages$readxl$Repository
    [1] "CRAN"
    
    $Packages$readxl$Hash
    [1] "63537c483c2dbec8d9e3183b3735254a"
    
    $Packages$readxl$Requirements
    $Packages$readxl$Requirements[[1]]
    [1] "Rcpp"
    
    $Packages$readxl$Requirements[[2]]
    [1] "cellranger"
    
    $Packages$readxl$Requirements[[3]]
    [1] "progress"
    
    $Packages$readxl$Requirements[[4]]
    [1] "tibble"
    
    
    
    $Packages$recipes
    $Packages$recipes$Package
    [1] "recipes"
    
    $Packages$recipes$Version
    [1] "0.2.0"
    
    $Packages$recipes$Source
    [1] "Repository"
    
    $Packages$recipes$Repository
    [1] "CRAN"
    
    $Packages$recipes$Hash
    [1] "2a0f64b148f3064f980404287b511f98"
    
    $Packages$recipes$Requirements
    $Packages$recipes$Requirements[[1]]
    [1] "Matrix"
    
    $Packages$recipes$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$recipes$Requirements[[3]]
    [1] "ellipsis"
    
    $Packages$recipes$Requirements[[4]]
    [1] "generics"
    
    $Packages$recipes$Requirements[[5]]
    [1] "glue"
    
    $Packages$recipes$Requirements[[6]]
    [1] "gower"
    
    $Packages$recipes$Requirements[[7]]
    [1] "hardhat"
    
    $Packages$recipes$Requirements[[8]]
    [1] "ipred"
    
    $Packages$recipes$Requirements[[9]]
    [1] "lifecycle"
    
    $Packages$recipes$Requirements[[10]]
    [1] "lubridate"
    
    $Packages$recipes$Requirements[[11]]
    [1] "magrittr"
    
    $Packages$recipes$Requirements[[12]]
    [1] "purrr"
    
    $Packages$recipes$Requirements[[13]]
    [1] "rlang"
    
    $Packages$recipes$Requirements[[14]]
    [1] "tibble"
    
    $Packages$recipes$Requirements[[15]]
    [1] "tidyr"
    
    $Packages$recipes$Requirements[[16]]
    [1] "tidyselect"
    
    $Packages$recipes$Requirements[[17]]
    [1] "timeDate"
    
    $Packages$recipes$Requirements[[18]]
    [1] "vctrs"
    
    $Packages$recipes$Requirements[[19]]
    [1] "withr"
    
    
    
    $Packages$rematch
    $Packages$rematch$Package
    [1] "rematch"
    
    $Packages$rematch$Version
    [1] "1.0.1"
    
    $Packages$rematch$Source
    [1] "Repository"
    
    $Packages$rematch$Repository
    [1] "CRAN"
    
    $Packages$rematch$Hash
    [1] "c66b930d20bb6d858cd18e1cebcfae5c"
    
    $Packages$rematch$Requirements
    list()
    
    
    $Packages$rematch2
    $Packages$rematch2$Package
    [1] "rematch2"
    
    $Packages$rematch2$Version
    [1] "2.1.2"
    
    $Packages$rematch2$Source
    [1] "Repository"
    
    $Packages$rematch2$Repository
    [1] "CRAN"
    
    $Packages$rematch2$Hash
    [1] "76c9e04c712a05848ae7a23d2f170a40"
    
    $Packages$rematch2$Requirements
    $Packages$rematch2$Requirements[[1]]
    [1] "tibble"
    
    
    
    $Packages$renv
    $Packages$renv$Package
    [1] "renv"
    
    $Packages$renv$Version
    [1] "0.15.4"
    
    $Packages$renv$Source
    [1] "Repository"
    
    $Packages$renv$Repository
    [1] "CRAN"
    
    $Packages$renv$Hash
    [1] "c1078316e1d4f70275fc1ea60c0bc431"
    
    $Packages$renv$Requirements
    list()
    
    
    $Packages$reprex
    $Packages$reprex$Package
    [1] "reprex"
    
    $Packages$reprex$Version
    [1] "2.0.1"
    
    $Packages$reprex$Source
    [1] "Repository"
    
    $Packages$reprex$Repository
    [1] "CRAN"
    
    $Packages$reprex$Hash
    [1] "911d101becedc0fde495bd910984bdc8"
    
    $Packages$reprex$Requirements
    $Packages$reprex$Requirements[[1]]
    [1] "callr"
    
    $Packages$reprex$Requirements[[2]]
    [1] "cli"
    
    $Packages$reprex$Requirements[[3]]
    [1] "clipr"
    
    $Packages$reprex$Requirements[[4]]
    [1] "fs"
    
    $Packages$reprex$Requirements[[5]]
    [1] "glue"
    
    $Packages$reprex$Requirements[[6]]
    [1] "knitr"
    
    $Packages$reprex$Requirements[[7]]
    [1] "rlang"
    
    $Packages$reprex$Requirements[[8]]
    [1] "rmarkdown"
    
    $Packages$reprex$Requirements[[9]]
    [1] "rstudioapi"
    
    $Packages$reprex$Requirements[[10]]
    [1] "withr"
    
    
    
    $Packages$reshape2
    $Packages$reshape2$Package
    [1] "reshape2"
    
    $Packages$reshape2$Version
    [1] "1.4.4"
    
    $Packages$reshape2$Source
    [1] "Repository"
    
    $Packages$reshape2$Repository
    [1] "CRAN"
    
    $Packages$reshape2$Hash
    [1] "bb5996d0bd962d214a11140d77589917"
    
    $Packages$reshape2$Requirements
    $Packages$reshape2$Requirements[[1]]
    [1] "Rcpp"
    
    $Packages$reshape2$Requirements[[2]]
    [1] "plyr"
    
    $Packages$reshape2$Requirements[[3]]
    [1] "stringr"
    
    
    
    $Packages$reticulate
    $Packages$reticulate$Package
    [1] "reticulate"
    
    $Packages$reticulate$Version
    [1] "1.24"
    
    $Packages$reticulate$Source
    [1] "Repository"
    
    $Packages$reticulate$Repository
    [1] "CRAN"
    
    $Packages$reticulate$Hash
    [1] "ffdf27627a3c1537478073c43b6e7980"
    
    $Packages$reticulate$Requirements
    $Packages$reticulate$Requirements[[1]]
    [1] "Matrix"
    
    $Packages$reticulate$Requirements[[2]]
    [1] "Rcpp"
    
    $Packages$reticulate$Requirements[[3]]
    [1] "RcppTOML"
    
    $Packages$reticulate$Requirements[[4]]
    [1] "here"
    
    $Packages$reticulate$Requirements[[5]]
    [1] "jsonlite"
    
    $Packages$reticulate$Requirements[[6]]
    [1] "png"
    
    $Packages$reticulate$Requirements[[7]]
    [1] "rappdirs"
    
    $Packages$reticulate$Requirements[[8]]
    [1] "withr"
    
    
    
    $Packages$rhdf5
    $Packages$rhdf5$Package
    [1] "rhdf5"
    
    $Packages$rhdf5$Version
    [1] "2.38.0"
    
    $Packages$rhdf5$Source
    [1] "Bioconductor"
    
    $Packages$rhdf5$git_url
    [1] "https://git.bioconductor.org/packages/rhdf5"
    
    $Packages$rhdf5$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$rhdf5$git_last_commit
    [1] "f6fdfa8"
    
    $Packages$rhdf5$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$rhdf5$Hash
    [1] "55f9159b07a2ec0e0f9e0dda28c8c48f"
    
    $Packages$rhdf5$Requirements
    $Packages$rhdf5$Requirements[[1]]
    [1] "Rhdf5lib"
    
    $Packages$rhdf5$Requirements[[2]]
    [1] "rhdf5filters"
    
    
    
    $Packages$rhdf5filters
    $Packages$rhdf5filters$Package
    [1] "rhdf5filters"
    
    $Packages$rhdf5filters$Version
    [1] "1.6.0"
    
    $Packages$rhdf5filters$Source
    [1] "Bioconductor"
    
    $Packages$rhdf5filters$git_url
    [1] "https://git.bioconductor.org/packages/rhdf5filters"
    
    $Packages$rhdf5filters$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$rhdf5filters$git_last_commit
    [1] "5f7f3a5"
    
    $Packages$rhdf5filters$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$rhdf5filters$Hash
    [1] "429e1c984126e2a2594ca16632008027"
    
    $Packages$rhdf5filters$Requirements
    $Packages$rhdf5filters$Requirements[[1]]
    [1] "Rhdf5lib"
    
    
    
    $Packages$rlang
    $Packages$rlang$Package
    [1] "rlang"
    
    $Packages$rlang$Version
    [1] "1.0.2"
    
    $Packages$rlang$Source
    [1] "Repository"
    
    $Packages$rlang$Repository
    [1] "CRAN"
    
    $Packages$rlang$Hash
    [1] "04884d9a75d778aca22c7154b8333ec9"
    
    $Packages$rlang$Requirements
    list()
    
    
    $Packages$rlist
    $Packages$rlist$Package
    [1] "rlist"
    
    $Packages$rlist$Version
    [1] "0.4.6.2"
    
    $Packages$rlist$Source
    [1] "Repository"
    
    $Packages$rlist$Repository
    [1] "CRAN"
    
    $Packages$rlist$Hash
    [1] "290c8ea0700d2e7258082d0025386e68"
    
    $Packages$rlist$Requirements
    $Packages$rlist$Requirements[[1]]
    [1] "XML"
    
    $Packages$rlist$Requirements[[2]]
    [1] "data.table"
    
    $Packages$rlist$Requirements[[3]]
    [1] "jsonlite"
    
    $Packages$rlist$Requirements[[4]]
    [1] "yaml"
    
    
    
    $Packages$rmarkdown
    $Packages$rmarkdown$Package
    [1] "rmarkdown"
    
    $Packages$rmarkdown$Version
    [1] "2.12"
    
    $Packages$rmarkdown$Source
    [1] "Repository"
    
    $Packages$rmarkdown$Repository
    [1] "CRAN"
    
    $Packages$rmarkdown$Hash
    [1] "354da5088ddfdffb73c11cc952885d88"
    
    $Packages$rmarkdown$Requirements
    $Packages$rmarkdown$Requirements[[1]]
    [1] "bslib"
    
    $Packages$rmarkdown$Requirements[[2]]
    [1] "evaluate"
    
    $Packages$rmarkdown$Requirements[[3]]
    [1] "htmltools"
    
    $Packages$rmarkdown$Requirements[[4]]
    [1] "jquerylib"
    
    $Packages$rmarkdown$Requirements[[5]]
    [1] "jsonlite"
    
    $Packages$rmarkdown$Requirements[[6]]
    [1] "knitr"
    
    $Packages$rmarkdown$Requirements[[7]]
    [1] "stringr"
    
    $Packages$rmarkdown$Requirements[[8]]
    [1] "tinytex"
    
    $Packages$rmarkdown$Requirements[[9]]
    [1] "xfun"
    
    $Packages$rmarkdown$Requirements[[10]]
    [1] "yaml"
    
    
    
    $Packages$robustbase
    $Packages$robustbase$Package
    [1] "robustbase"
    
    $Packages$robustbase$Version
    [1] "0.93-9"
    
    $Packages$robustbase$Source
    [1] "Repository"
    
    $Packages$robustbase$Repository
    [1] "CRAN"
    
    $Packages$robustbase$Hash
    [1] "e7b310bebca9aca51848b0a366a427cf"
    
    $Packages$robustbase$Requirements
    $Packages$robustbase$Requirements[[1]]
    [1] "DEoptimR"
    
    
    
    $Packages$rpart
    $Packages$rpart$Package
    [1] "rpart"
    
    $Packages$rpart$Version
    [1] "4.1.16"
    
    $Packages$rpart$Source
    [1] "Repository"
    
    $Packages$rpart$Repository
    [1] "CRAN"
    
    $Packages$rpart$Hash
    [1] "ea3ca1d9473daabb3cd0f1b4f974c1ed"
    
    $Packages$rpart$Requirements
    list()
    
    
    $Packages$rprojroot
    $Packages$rprojroot$Package
    [1] "rprojroot"
    
    $Packages$rprojroot$Version
    [1] "2.0.2"
    
    $Packages$rprojroot$Source
    [1] "Repository"
    
    $Packages$rprojroot$Repository
    [1] "CRAN"
    
    $Packages$rprojroot$Hash
    [1] "249d8cd1e74a8f6a26194a91b47f21d1"
    
    $Packages$rprojroot$Requirements
    list()
    
    
    $Packages$rsample
    $Packages$rsample$Package
    [1] "rsample"
    
    $Packages$rsample$Version
    [1] "0.1.1"
    
    $Packages$rsample$Source
    [1] "Repository"
    
    $Packages$rsample$Repository
    [1] "CRAN"
    
    $Packages$rsample$Hash
    [1] "d1ac9b0eedb5dabcec1309e072b169cb"
    
    $Packages$rsample$Requirements
    $Packages$rsample$Requirements[[1]]
    [1] "dplyr"
    
    $Packages$rsample$Requirements[[2]]
    [1] "ellipsis"
    
    $Packages$rsample$Requirements[[3]]
    [1] "furrr"
    
    $Packages$rsample$Requirements[[4]]
    [1] "generics"
    
    $Packages$rsample$Requirements[[5]]
    [1] "lifecycle"
    
    $Packages$rsample$Requirements[[6]]
    [1] "purrr"
    
    $Packages$rsample$Requirements[[7]]
    [1] "rlang"
    
    $Packages$rsample$Requirements[[8]]
    [1] "slider"
    
    $Packages$rsample$Requirements[[9]]
    [1] "tibble"
    
    $Packages$rsample$Requirements[[10]]
    [1] "tidyr"
    
    $Packages$rsample$Requirements[[11]]
    [1] "tidyselect"
    
    $Packages$rsample$Requirements[[12]]
    [1] "vctrs"
    
    
    
    $Packages$rstudioapi
    $Packages$rstudioapi$Package
    [1] "rstudioapi"
    
    $Packages$rstudioapi$Version
    [1] "0.13"
    
    $Packages$rstudioapi$Source
    [1] "Repository"
    
    $Packages$rstudioapi$Repository
    [1] "CRAN"
    
    $Packages$rstudioapi$Hash
    [1] "06c85365a03fdaf699966cc1d3cf53ea"
    
    $Packages$rstudioapi$Requirements
    list()
    
    
    $Packages$rsvd
    $Packages$rsvd$Package
    [1] "rsvd"
    
    $Packages$rsvd$Version
    [1] "1.0.5"
    
    $Packages$rsvd$Source
    [1] "Repository"
    
    $Packages$rsvd$Repository
    [1] "CRAN"
    
    $Packages$rsvd$Hash
    [1] "b462187d887abc519894874486dbd6fd"
    
    $Packages$rsvd$Requirements
    $Packages$rsvd$Requirements[[1]]
    [1] "Matrix"
    
    
    
    $Packages$rvest
    $Packages$rvest$Package
    [1] "rvest"
    
    $Packages$rvest$Version
    [1] "1.0.2"
    
    $Packages$rvest$Source
    [1] "Repository"
    
    $Packages$rvest$Repository
    [1] "CRAN"
    
    $Packages$rvest$Hash
    [1] "bb099886deffecd6f9b298b7d4492943"
    
    $Packages$rvest$Requirements
    $Packages$rvest$Requirements[[1]]
    [1] "httr"
    
    $Packages$rvest$Requirements[[2]]
    [1] "lifecycle"
    
    $Packages$rvest$Requirements[[3]]
    [1] "magrittr"
    
    $Packages$rvest$Requirements[[4]]
    [1] "rlang"
    
    $Packages$rvest$Requirements[[5]]
    [1] "selectr"
    
    $Packages$rvest$Requirements[[6]]
    [1] "tibble"
    
    $Packages$rvest$Requirements[[7]]
    [1] "xml2"
    
    
    
    $Packages$sass
    $Packages$sass$Package
    [1] "sass"
    
    $Packages$sass$Version
    [1] "0.4.0"
    
    $Packages$sass$Source
    [1] "Repository"
    
    $Packages$sass$Repository
    [1] "CRAN"
    
    $Packages$sass$Hash
    [1] "50cf822feb64bb3977bda0b7091be623"
    
    $Packages$sass$Requirements
    $Packages$sass$Requirements[[1]]
    [1] "R6"
    
    $Packages$sass$Requirements[[2]]
    [1] "fs"
    
    $Packages$sass$Requirements[[3]]
    [1] "htmltools"
    
    $Packages$sass$Requirements[[4]]
    [1] "rappdirs"
    
    $Packages$sass$Requirements[[5]]
    [1] "rlang"
    
    
    
    $Packages$scales
    $Packages$scales$Package
    [1] "scales"
    
    $Packages$scales$Version
    [1] "1.1.1"
    
    $Packages$scales$Source
    [1] "Repository"
    
    $Packages$scales$Repository
    [1] "CRAN"
    
    $Packages$scales$Hash
    [1] "6f76f71042411426ec8df6c54f34e6dd"
    
    $Packages$scales$Requirements
    $Packages$scales$Requirements[[1]]
    [1] "R6"
    
    $Packages$scales$Requirements[[2]]
    [1] "RColorBrewer"
    
    $Packages$scales$Requirements[[3]]
    [1] "farver"
    
    $Packages$scales$Requirements[[4]]
    [1] "labeling"
    
    $Packages$scales$Requirements[[5]]
    [1] "lifecycle"
    
    $Packages$scales$Requirements[[6]]
    [1] "munsell"
    
    $Packages$scales$Requirements[[7]]
    [1] "viridisLite"
    
    
    
    $Packages$scater
    $Packages$scater$Package
    [1] "scater"
    
    $Packages$scater$Version
    [1] "1.22.0"
    
    $Packages$scater$Source
    [1] "Bioconductor"
    
    $Packages$scater$git_url
    [1] "https://git.bioconductor.org/packages/scater"
    
    $Packages$scater$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$scater$git_last_commit
    [1] "ea2c95c"
    
    $Packages$scater$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$scater$Hash
    [1] "8e229c5126b386c658047881fe2e5299"
    
    $Packages$scater$Requirements
    $Packages$scater$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$scater$Requirements[[2]]
    [1] "BiocNeighbors"
    
    $Packages$scater$Requirements[[3]]
    [1] "BiocParallel"
    
    $Packages$scater$Requirements[[4]]
    [1] "BiocSingular"
    
    $Packages$scater$Requirements[[5]]
    [1] "DelayedArray"
    
    $Packages$scater$Requirements[[6]]
    [1] "DelayedMatrixStats"
    
    $Packages$scater$Requirements[[7]]
    [1] "Matrix"
    
    $Packages$scater$Requirements[[8]]
    [1] "RColorBrewer"
    
    $Packages$scater$Requirements[[9]]
    [1] "Rtsne"
    
    $Packages$scater$Requirements[[10]]
    [1] "S4Vectors"
    
    $Packages$scater$Requirements[[11]]
    [1] "SingleCellExperiment"
    
    $Packages$scater$Requirements[[12]]
    [1] "SummarizedExperiment"
    
    $Packages$scater$Requirements[[13]]
    [1] "beachmat"
    
    $Packages$scater$Requirements[[14]]
    [1] "ggbeeswarm"
    
    $Packages$scater$Requirements[[15]]
    [1] "ggplot2"
    
    $Packages$scater$Requirements[[16]]
    [1] "ggrepel"
    
    $Packages$scater$Requirements[[17]]
    [1] "gridExtra"
    
    $Packages$scater$Requirements[[18]]
    [1] "rlang"
    
    $Packages$scater$Requirements[[19]]
    [1] "scuttle"
    
    $Packages$scater$Requirements[[20]]
    [1] "viridis"
    
    
    
    $Packages$scuttle
    $Packages$scuttle$Package
    [1] "scuttle"
    
    $Packages$scuttle$Version
    [1] "1.4.0"
    
    $Packages$scuttle$Source
    [1] "Bioconductor"
    
    $Packages$scuttle$git_url
    [1] "https://git.bioconductor.org/packages/scuttle"
    
    $Packages$scuttle$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$scuttle$git_last_commit
    [1] "b335263"
    
    $Packages$scuttle$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$scuttle$Hash
    [1] "1c612f905434c933d27f508f294da9e0"
    
    $Packages$scuttle$Requirements
    $Packages$scuttle$Requirements[[1]]
    [1] "BiocGenerics"
    
    $Packages$scuttle$Requirements[[2]]
    [1] "BiocParallel"
    
    $Packages$scuttle$Requirements[[3]]
    [1] "DelayedArray"
    
    $Packages$scuttle$Requirements[[4]]
    [1] "DelayedMatrixStats"
    
    $Packages$scuttle$Requirements[[5]]
    [1] "GenomicRanges"
    
    $Packages$scuttle$Requirements[[6]]
    [1] "Matrix"
    
    $Packages$scuttle$Requirements[[7]]
    [1] "Rcpp"
    
    $Packages$scuttle$Requirements[[8]]
    [1] "S4Vectors"
    
    $Packages$scuttle$Requirements[[9]]
    [1] "SingleCellExperiment"
    
    $Packages$scuttle$Requirements[[10]]
    [1] "SummarizedExperiment"
    
    $Packages$scuttle$Requirements[[11]]
    [1] "beachmat"
    
    
    
    $Packages$selectr
    $Packages$selectr$Package
    [1] "selectr"
    
    $Packages$selectr$Version
    [1] "0.4-2"
    
    $Packages$selectr$Source
    [1] "Repository"
    
    $Packages$selectr$Repository
    [1] "CRAN"
    
    $Packages$selectr$Hash
    [1] "3838071b66e0c566d55cc26bd6e27bf4"
    
    $Packages$selectr$Requirements
    $Packages$selectr$Requirements[[1]]
    [1] "R6"
    
    $Packages$selectr$Requirements[[2]]
    [1] "stringr"
    
    
    
    $Packages$shiny
    $Packages$shiny$Package
    [1] "shiny"
    
    $Packages$shiny$Version
    [1] "1.7.1"
    
    $Packages$shiny$Source
    [1] "Repository"
    
    $Packages$shiny$Repository
    [1] "CRAN"
    
    $Packages$shiny$Hash
    [1] "00344c227c7bd0ab5d78052c5d736c44"
    
    $Packages$shiny$Requirements
    $Packages$shiny$Requirements[[1]]
    [1] "R6"
    
    $Packages$shiny$Requirements[[2]]
    [1] "bslib"
    
    $Packages$shiny$Requirements[[3]]
    [1] "cachem"
    
    $Packages$shiny$Requirements[[4]]
    [1] "commonmark"
    
    $Packages$shiny$Requirements[[5]]
    [1] "crayon"
    
    $Packages$shiny$Requirements[[6]]
    [1] "ellipsis"
    
    $Packages$shiny$Requirements[[7]]
    [1] "fastmap"
    
    $Packages$shiny$Requirements[[8]]
    [1] "fontawesome"
    
    $Packages$shiny$Requirements[[9]]
    [1] "glue"
    
    $Packages$shiny$Requirements[[10]]
    [1] "htmltools"
    
    $Packages$shiny$Requirements[[11]]
    [1] "httpuv"
    
    $Packages$shiny$Requirements[[12]]
    [1] "jsonlite"
    
    $Packages$shiny$Requirements[[13]]
    [1] "later"
    
    $Packages$shiny$Requirements[[14]]
    [1] "lifecycle"
    
    $Packages$shiny$Requirements[[15]]
    [1] "mime"
    
    $Packages$shiny$Requirements[[16]]
    [1] "promises"
    
    $Packages$shiny$Requirements[[17]]
    [1] "rlang"
    
    $Packages$shiny$Requirements[[18]]
    [1] "sourcetools"
    
    $Packages$shiny$Requirements[[19]]
    [1] "withr"
    
    $Packages$shiny$Requirements[[20]]
    [1] "xtable"
    
    
    
    $Packages$slider
    $Packages$slider$Package
    [1] "slider"
    
    $Packages$slider$Version
    [1] "0.2.2"
    
    $Packages$slider$Source
    [1] "Repository"
    
    $Packages$slider$Repository
    [1] "CRAN"
    
    $Packages$slider$Hash
    [1] "5237bd176dc0c4dd7eb8dcdafe514de3"
    
    $Packages$slider$Requirements
    $Packages$slider$Requirements[[1]]
    [1] "ellipsis"
    
    $Packages$slider$Requirements[[2]]
    [1] "glue"
    
    $Packages$slider$Requirements[[3]]
    [1] "rlang"
    
    $Packages$slider$Requirements[[4]]
    [1] "vctrs"
    
    $Packages$slider$Requirements[[5]]
    [1] "warp"
    
    
    
    $Packages$snow
    $Packages$snow$Package
    [1] "snow"
    
    $Packages$snow$Version
    [1] "0.4-4"
    
    $Packages$snow$Source
    [1] "Repository"
    
    $Packages$snow$Repository
    [1] "CRAN"
    
    $Packages$snow$Hash
    [1] "40b74690debd20c57d93d8c246b305d4"
    
    $Packages$snow$Requirements
    list()
    
    
    $Packages$sourcetools
    $Packages$sourcetools$Package
    [1] "sourcetools"
    
    $Packages$sourcetools$Version
    [1] "0.1.7"
    
    $Packages$sourcetools$Source
    [1] "Repository"
    
    $Packages$sourcetools$Repository
    [1] "CRAN"
    
    $Packages$sourcetools$Hash
    [1] "947e4e02a79effa5d512473e10f41797"
    
    $Packages$sourcetools$Requirements
    list()
    
    
    $Packages$sp
    $Packages$sp$Package
    [1] "sp"
    
    $Packages$sp$Version
    [1] "1.4-6"
    
    $Packages$sp$Source
    [1] "Repository"
    
    $Packages$sp$Repository
    [1] "CRAN"
    
    $Packages$sp$Hash
    [1] "ce8613f4e8c84ef4da9eba65b874ebe9"
    
    $Packages$sp$Requirements
    $Packages$sp$Requirements[[1]]
    [1] "lattice"
    
    
    
    $Packages$sparseMatrixStats
    $Packages$sparseMatrixStats$Package
    [1] "sparseMatrixStats"
    
    $Packages$sparseMatrixStats$Version
    [1] "1.6.0"
    
    $Packages$sparseMatrixStats$Source
    [1] "Bioconductor"
    
    $Packages$sparseMatrixStats$git_url
    [1] "https://git.bioconductor.org/packages/sparseMatrixStats"
    
    $Packages$sparseMatrixStats$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$sparseMatrixStats$git_last_commit
    [1] "78627a8"
    
    $Packages$sparseMatrixStats$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$sparseMatrixStats$Hash
    [1] "f8dd82d2581115df0ea380c91ed8b9f6"
    
    $Packages$sparseMatrixStats$Requirements
    $Packages$sparseMatrixStats$Requirements[[1]]
    [1] "Matrix"
    
    $Packages$sparseMatrixStats$Requirements[[2]]
    [1] "MatrixGenerics"
    
    $Packages$sparseMatrixStats$Requirements[[3]]
    [1] "Rcpp"
    
    $Packages$sparseMatrixStats$Requirements[[4]]
    [1] "matrixStats"
    
    
    
    $Packages$speedyseq
    $Packages$speedyseq$Package
    [1] "speedyseq"
    
    $Packages$speedyseq$Version
    [1] "0.5.3.9018"
    
    $Packages$speedyseq$Source
    [1] "GitHub"
    
    $Packages$speedyseq$RemoteType
    [1] "github"
    
    $Packages$speedyseq$RemoteHost
    [1] "api.github.com"
    
    $Packages$speedyseq$RemoteUsername
    [1] "mikemc"
    
    $Packages$speedyseq$RemoteRepo
    [1] "speedyseq"
    
    $Packages$speedyseq$RemoteRef
    [1] "main"
    
    $Packages$speedyseq$RemoteSha
    [1] "ceb941fdd482fe4bf9610f80970050e24f369be9"
    
    $Packages$speedyseq$Hash
    [1] "040ffc1e11abd6bc5c9f0be6bc5ab11a"
    
    $Packages$speedyseq$Requirements
    $Packages$speedyseq$Requirements[[1]]
    [1] "Biostrings"
    
    $Packages$speedyseq$Requirements[[2]]
    [1] "ape"
    
    $Packages$speedyseq$Requirements[[3]]
    [1] "castor"
    
    $Packages$speedyseq$Requirements[[4]]
    [1] "data.table"
    
    $Packages$speedyseq$Requirements[[5]]
    [1] "dplyr"
    
    $Packages$speedyseq$Requirements[[6]]
    [1] "ggplot2"
    
    $Packages$speedyseq$Requirements[[7]]
    [1] "magrittr"
    
    $Packages$speedyseq$Requirements[[8]]
    [1] "phyloseq"
    
    $Packages$speedyseq$Requirements[[9]]
    [1] "purrr"
    
    $Packages$speedyseq$Requirements[[10]]
    [1] "rlang"
    
    $Packages$speedyseq$Requirements[[11]]
    [1] "scales"
    
    $Packages$speedyseq$Requirements[[12]]
    [1] "stringr"
    
    $Packages$speedyseq$Requirements[[13]]
    [1] "tibble"
    
    $Packages$speedyseq$Requirements[[14]]
    [1] "tidyr"
    
    $Packages$speedyseq$Requirements[[15]]
    [1] "vctrs"
    
    $Packages$speedyseq$Requirements[[16]]
    [1] "vegan"
    
    
    
    $Packages$stringdist
    $Packages$stringdist$Package
    [1] "stringdist"
    
    $Packages$stringdist$Version
    [1] "0.9.8"
    
    $Packages$stringdist$Source
    [1] "Repository"
    
    $Packages$stringdist$Repository
    [1] "CRAN"
    
    $Packages$stringdist$Hash
    [1] "c8b6c79c3b2d7d8475f3ab35cd99737b"
    
    $Packages$stringdist$Requirements
    list()
    
    
    $Packages$stringfish
    $Packages$stringfish$Package
    [1] "stringfish"
    
    $Packages$stringfish$Version
    [1] "0.15.5"
    
    $Packages$stringfish$Source
    [1] "Repository"
    
    $Packages$stringfish$Repository
    [1] "CRAN"
    
    $Packages$stringfish$Hash
    [1] "f5a78e4562f6390f61ddf14848a0a3d9"
    
    $Packages$stringfish$Requirements
    $Packages$stringfish$Requirements[[1]]
    [1] "Rcpp"
    
    $Packages$stringfish$Requirements[[2]]
    [1] "RcppParallel"
    
    
    
    $Packages$stringi
    $Packages$stringi$Package
    [1] "stringi"
    
    $Packages$stringi$Version
    [1] "1.7.6"
    
    $Packages$stringi$Source
    [1] "Repository"
    
    $Packages$stringi$Repository
    [1] "CRAN"
    
    $Packages$stringi$Hash
    [1] "bba431031d30789535745a9627ac9271"
    
    $Packages$stringi$Requirements
    list()
    
    
    $Packages$stringr
    $Packages$stringr$Package
    [1] "stringr"
    
    $Packages$stringr$Version
    [1] "1.4.0"
    
    $Packages$stringr$Source
    [1] "Repository"
    
    $Packages$stringr$Repository
    [1] "CRAN"
    
    $Packages$stringr$Hash
    [1] "0759e6b6c0957edb1311028a49a35e76"
    
    $Packages$stringr$Requirements
    $Packages$stringr$Requirements[[1]]
    [1] "glue"
    
    $Packages$stringr$Requirements[[2]]
    [1] "magrittr"
    
    $Packages$stringr$Requirements[[3]]
    [1] "stringi"
    
    
    
    $Packages$survival
    $Packages$survival$Package
    [1] "survival"
    
    $Packages$survival$Version
    [1] "3.3-1"
    
    $Packages$survival$Source
    [1] "Repository"
    
    $Packages$survival$Repository
    [1] "CRAN"
    
    $Packages$survival$Hash
    [1] "f6189c70451d3d68e0d571235576e833"
    
    $Packages$survival$Requirements
    $Packages$survival$Requirements[[1]]
    [1] "Matrix"
    
    
    
    $Packages$svglite
    $Packages$svglite$Package
    [1] "svglite"
    
    $Packages$svglite$Version
    [1] "2.1.0"
    
    $Packages$svglite$Source
    [1] "Repository"
    
    $Packages$svglite$Repository
    [1] "CRAN"
    
    $Packages$svglite$Hash
    [1] "68dfdf211af6aa4e5f050f064f64d401"
    
    $Packages$svglite$Requirements
    $Packages$svglite$Requirements[[1]]
    [1] "cpp11"
    
    $Packages$svglite$Requirements[[2]]
    [1] "systemfonts"
    
    
    
    $Packages$sys
    $Packages$sys$Package
    [1] "sys"
    
    $Packages$sys$Version
    [1] "3.4"
    
    $Packages$sys$Source
    [1] "Repository"
    
    $Packages$sys$Repository
    [1] "CRAN"
    
    $Packages$sys$Hash
    [1] "b227d13e29222b4574486cfcbde077fa"
    
    $Packages$sys$Requirements
    list()
    
    
    $Packages$systemfonts
    $Packages$systemfonts$Package
    [1] "systemfonts"
    
    $Packages$systemfonts$Version
    [1] "1.0.4"
    
    $Packages$systemfonts$Source
    [1] "Repository"
    
    $Packages$systemfonts$Repository
    [1] "CRAN"
    
    $Packages$systemfonts$Hash
    [1] "90b28393209827327de889f49935140a"
    
    $Packages$systemfonts$Requirements
    $Packages$systemfonts$Requirements[[1]]
    [1] "cpp11"
    
    
    
    $Packages$tarchetypes
    $Packages$tarchetypes$Package
    [1] "tarchetypes"
    
    $Packages$tarchetypes$Version
    [1] "0.4.1"
    
    $Packages$tarchetypes$Source
    [1] "Repository"
    
    $Packages$tarchetypes$Repository
    [1] "CRAN"
    
    $Packages$tarchetypes$Hash
    [1] "24c9eb58a70cef665befbcc9891572de"
    
    $Packages$tarchetypes$Requirements
    $Packages$tarchetypes$Requirements[[1]]
    [1] "digest"
    
    $Packages$tarchetypes$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$tarchetypes$Requirements[[3]]
    [1] "fs"
    
    $Packages$tarchetypes$Requirements[[4]]
    [1] "rlang"
    
    $Packages$tarchetypes$Requirements[[5]]
    [1] "targets"
    
    $Packages$tarchetypes$Requirements[[6]]
    [1] "tibble"
    
    $Packages$tarchetypes$Requirements[[7]]
    [1] "tidyselect"
    
    $Packages$tarchetypes$Requirements[[8]]
    [1] "vctrs"
    
    $Packages$tarchetypes$Requirements[[9]]
    [1] "withr"
    
    
    
    $Packages$targets
    $Packages$targets$Package
    [1] "targets"
    
    $Packages$targets$Version
    [1] "0.11.0"
    
    $Packages$targets$Source
    [1] "Repository"
    
    $Packages$targets$Repository
    [1] "CRAN"
    
    $Packages$targets$Hash
    [1] "bcb5125b621b728d8dc23e39f2f07428"
    
    $Packages$targets$Requirements
    $Packages$targets$Requirements[[1]]
    [1] "R6"
    
    $Packages$targets$Requirements[[2]]
    [1] "base64url"
    
    $Packages$targets$Requirements[[3]]
    [1] "callr"
    
    $Packages$targets$Requirements[[4]]
    [1] "cli"
    
    $Packages$targets$Requirements[[5]]
    [1] "codetools"
    
    $Packages$targets$Requirements[[6]]
    [1] "data.table"
    
    $Packages$targets$Requirements[[7]]
    [1] "digest"
    
    $Packages$targets$Requirements[[8]]
    [1] "igraph"
    
    $Packages$targets$Requirements[[9]]
    [1] "knitr"
    
    $Packages$targets$Requirements[[10]]
    [1] "rlang"
    
    $Packages$targets$Requirements[[11]]
    [1] "tibble"
    
    $Packages$targets$Requirements[[12]]
    [1] "tidyselect"
    
    $Packages$targets$Requirements[[13]]
    [1] "vctrs"
    
    $Packages$targets$Requirements[[14]]
    [1] "withr"
    
    $Packages$targets$Requirements[[15]]
    [1] "yaml"
    
    
    
    $Packages$taxadb
    $Packages$taxadb$Package
    [1] "taxadb"
    
    $Packages$taxadb$Version
    [1] "0.1.4"
    
    $Packages$taxadb$Source
    [1] "Repository"
    
    $Packages$taxadb$Repository
    [1] "CRAN"
    
    $Packages$taxadb$Hash
    [1] "9586ab951460e1f40350fd68f7d37c4f"
    
    $Packages$taxadb$Requirements
    $Packages$taxadb$Requirements[[1]]
    [1] "DBI"
    
    $Packages$taxadb$Requirements[[2]]
    [1] "R.utils"
    
    $Packages$taxadb$Requirements[[3]]
    [1] "arkdb"
    
    $Packages$taxadb$Requirements[[4]]
    [1] "contentid"
    
    $Packages$taxadb$Requirements[[5]]
    [1] "curl"
    
    $Packages$taxadb$Requirements[[6]]
    [1] "dbplyr"
    
    $Packages$taxadb$Requirements[[7]]
    [1] "dplyr"
    
    $Packages$taxadb$Requirements[[8]]
    [1] "duckdb"
    
    $Packages$taxadb$Requirements[[9]]
    [1] "jsonlite"
    
    $Packages$taxadb$Requirements[[10]]
    [1] "magrittr"
    
    $Packages$taxadb$Requirements[[11]]
    [1] "progress"
    
    $Packages$taxadb$Requirements[[12]]
    [1] "rappdirs"
    
    $Packages$taxadb$Requirements[[13]]
    [1] "readr"
    
    $Packages$taxadb$Requirements[[14]]
    [1] "rlang"
    
    $Packages$taxadb$Requirements[[15]]
    [1] "stringi"
    
    $Packages$taxadb$Requirements[[16]]
    [1] "tibble"
    
    
    
    $Packages$taxizedb
    $Packages$taxizedb$Package
    [1] "taxizedb"
    
    $Packages$taxizedb$Version
    [1] "0.3.0"
    
    $Packages$taxizedb$Source
    [1] "Repository"
    
    $Packages$taxizedb$Repository
    [1] "CRAN"
    
    $Packages$taxizedb$Hash
    [1] "9303b0cde2b629f5a831e171f6f14b41"
    
    $Packages$taxizedb$Requirements
    $Packages$taxizedb$Requirements[[1]]
    [1] "DBI"
    
    $Packages$taxizedb$Requirements[[2]]
    [1] "RSQLite"
    
    $Packages$taxizedb$Requirements[[3]]
    [1] "curl"
    
    $Packages$taxizedb$Requirements[[4]]
    [1] "dbplyr"
    
    $Packages$taxizedb$Requirements[[5]]
    [1] "dplyr"
    
    $Packages$taxizedb$Requirements[[6]]
    [1] "hoardr"
    
    $Packages$taxizedb$Requirements[[7]]
    [1] "magrittr"
    
    $Packages$taxizedb$Requirements[[8]]
    [1] "readr"
    
    $Packages$taxizedb$Requirements[[9]]
    [1] "rlang"
    
    $Packages$taxizedb$Requirements[[10]]
    [1] "tibble"
    
    
    
    $Packages$tensorA
    $Packages$tensorA$Package
    [1] "tensorA"
    
    $Packages$tensorA$Version
    [1] "0.36.2"
    
    $Packages$tensorA$Source
    [1] "Repository"
    
    $Packages$tensorA$Repository
    [1] "CRAN"
    
    $Packages$tensorA$Hash
    [1] "fd792ceac77f96b647fa8d6e1788969a"
    
    $Packages$tensorA$Requirements
    list()
    
    
    $Packages$testthat
    $Packages$testthat$Package
    [1] "testthat"
    
    $Packages$testthat$Version
    [1] "3.1.2"
    
    $Packages$testthat$Source
    [1] "Repository"
    
    $Packages$testthat$Repository
    [1] "CRAN"
    
    $Packages$testthat$Hash
    [1] "32454e5780e8dbe31e4b61b13d8918fe"
    
    $Packages$testthat$Requirements
    $Packages$testthat$Requirements[[1]]
    [1] "R6"
    
    $Packages$testthat$Requirements[[2]]
    [1] "brio"
    
    $Packages$testthat$Requirements[[3]]
    [1] "callr"
    
    $Packages$testthat$Requirements[[4]]
    [1] "cli"
    
    $Packages$testthat$Requirements[[5]]
    [1] "crayon"
    
    $Packages$testthat$Requirements[[6]]
    [1] "desc"
    
    $Packages$testthat$Requirements[[7]]
    [1] "digest"
    
    $Packages$testthat$Requirements[[8]]
    [1] "ellipsis"
    
    $Packages$testthat$Requirements[[9]]
    [1] "evaluate"
    
    $Packages$testthat$Requirements[[10]]
    [1] "jsonlite"
    
    $Packages$testthat$Requirements[[11]]
    [1] "lifecycle"
    
    $Packages$testthat$Requirements[[12]]
    [1] "magrittr"
    
    $Packages$testthat$Requirements[[13]]
    [1] "pkgload"
    
    $Packages$testthat$Requirements[[14]]
    [1] "praise"
    
    $Packages$testthat$Requirements[[15]]
    [1] "processx"
    
    $Packages$testthat$Requirements[[16]]
    [1] "ps"
    
    $Packages$testthat$Requirements[[17]]
    [1] "rlang"
    
    $Packages$testthat$Requirements[[18]]
    [1] "waldo"
    
    $Packages$testthat$Requirements[[19]]
    [1] "withr"
    
    
    
    $Packages$themis
    $Packages$themis$Package
    [1] "themis"
    
    $Packages$themis$Version
    [1] "0.1.4"
    
    $Packages$themis$Source
    [1] "Repository"
    
    $Packages$themis$Repository
    [1] "CRAN"
    
    $Packages$themis$Hash
    [1] "42b548d48b359286942c20fcd78c348f"
    
    $Packages$themis$Requirements
    $Packages$themis$Requirements[[1]]
    [1] "RANN"
    
    $Packages$themis$Requirements[[2]]
    [1] "ROSE"
    
    $Packages$themis$Requirements[[3]]
    [1] "dplyr"
    
    $Packages$themis$Requirements[[4]]
    [1] "generics"
    
    $Packages$themis$Requirements[[5]]
    [1] "purrr"
    
    $Packages$themis$Requirements[[6]]
    [1] "recipes"
    
    $Packages$themis$Requirements[[7]]
    [1] "rlang"
    
    $Packages$themis$Requirements[[8]]
    [1] "tibble"
    
    $Packages$themis$Requirements[[9]]
    [1] "unbalanced"
    
    $Packages$themis$Requirements[[10]]
    [1] "withr"
    
    
    
    $Packages$tibble
    $Packages$tibble$Package
    [1] "tibble"
    
    $Packages$tibble$Version
    [1] "3.1.6"
    
    $Packages$tibble$Source
    [1] "Repository"
    
    $Packages$tibble$Repository
    [1] "CRAN"
    
    $Packages$tibble$Hash
    [1] "8a8f02d1934dfd6431c671361510dd0b"
    
    $Packages$tibble$Requirements
    $Packages$tibble$Requirements[[1]]
    [1] "ellipsis"
    
    $Packages$tibble$Requirements[[2]]
    [1] "fansi"
    
    $Packages$tibble$Requirements[[3]]
    [1] "lifecycle"
    
    $Packages$tibble$Requirements[[4]]
    [1] "magrittr"
    
    $Packages$tibble$Requirements[[5]]
    [1] "pillar"
    
    $Packages$tibble$Requirements[[6]]
    [1] "pkgconfig"
    
    $Packages$tibble$Requirements[[7]]
    [1] "rlang"
    
    $Packages$tibble$Requirements[[8]]
    [1] "vctrs"
    
    
    
    $Packages$tidymodels
    $Packages$tidymodels$Package
    [1] "tidymodels"
    
    $Packages$tidymodels$Version
    [1] "0.1.4"
    
    $Packages$tidymodels$Source
    [1] "Repository"
    
    $Packages$tidymodels$Repository
    [1] "CRAN"
    
    $Packages$tidymodels$Hash
    [1] "74cfd3e1de860e43453b7b36675f6480"
    
    $Packages$tidymodels$Requirements
    $Packages$tidymodels$Requirements[[1]]
    [1] "broom"
    
    $Packages$tidymodels$Requirements[[2]]
    [1] "cli"
    
    $Packages$tidymodels$Requirements[[3]]
    [1] "conflicted"
    
    $Packages$tidymodels$Requirements[[4]]
    [1] "dials"
    
    $Packages$tidymodels$Requirements[[5]]
    [1] "dplyr"
    
    $Packages$tidymodels$Requirements[[6]]
    [1] "ggplot2"
    
    $Packages$tidymodels$Requirements[[7]]
    [1] "hardhat"
    
    $Packages$tidymodels$Requirements[[8]]
    [1] "infer"
    
    $Packages$tidymodels$Requirements[[9]]
    [1] "modeldata"
    
    $Packages$tidymodels$Requirements[[10]]
    [1] "parsnip"
    
    $Packages$tidymodels$Requirements[[11]]
    [1] "purrr"
    
    $Packages$tidymodels$Requirements[[12]]
    [1] "recipes"
    
    $Packages$tidymodels$Requirements[[13]]
    [1] "rlang"
    
    $Packages$tidymodels$Requirements[[14]]
    [1] "rsample"
    
    $Packages$tidymodels$Requirements[[15]]
    [1] "rstudioapi"
    
    $Packages$tidymodels$Requirements[[16]]
    [1] "tibble"
    
    $Packages$tidymodels$Requirements[[17]]
    [1] "tidyr"
    
    $Packages$tidymodels$Requirements[[18]]
    [1] "tune"
    
    $Packages$tidymodels$Requirements[[19]]
    [1] "workflows"
    
    $Packages$tidymodels$Requirements[[20]]
    [1] "workflowsets"
    
    $Packages$tidymodels$Requirements[[21]]
    [1] "yardstick"
    
    
    
    $Packages$tidyr
    $Packages$tidyr$Package
    [1] "tidyr"
    
    $Packages$tidyr$Version
    [1] "1.2.0"
    
    $Packages$tidyr$Source
    [1] "Repository"
    
    $Packages$tidyr$Repository
    [1] "CRAN"
    
    $Packages$tidyr$Hash
    [1] "d8b95b7fee945d7da6888cf7eb71a49c"
    
    $Packages$tidyr$Requirements
    $Packages$tidyr$Requirements[[1]]
    [1] "cpp11"
    
    $Packages$tidyr$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$tidyr$Requirements[[3]]
    [1] "ellipsis"
    
    $Packages$tidyr$Requirements[[4]]
    [1] "glue"
    
    $Packages$tidyr$Requirements[[5]]
    [1] "lifecycle"
    
    $Packages$tidyr$Requirements[[6]]
    [1] "magrittr"
    
    $Packages$tidyr$Requirements[[7]]
    [1] "purrr"
    
    $Packages$tidyr$Requirements[[8]]
    [1] "rlang"
    
    $Packages$tidyr$Requirements[[9]]
    [1] "tibble"
    
    $Packages$tidyr$Requirements[[10]]
    [1] "tidyselect"
    
    $Packages$tidyr$Requirements[[11]]
    [1] "vctrs"
    
    
    
    $Packages$tidyselect
    $Packages$tidyselect$Package
    [1] "tidyselect"
    
    $Packages$tidyselect$Version
    [1] "1.1.2"
    
    $Packages$tidyselect$Source
    [1] "Repository"
    
    $Packages$tidyselect$Repository
    [1] "CRAN"
    
    $Packages$tidyselect$Hash
    [1] "17f6da8cfd7002760a859915ce7eef8f"
    
    $Packages$tidyselect$Requirements
    $Packages$tidyselect$Requirements[[1]]
    [1] "ellipsis"
    
    $Packages$tidyselect$Requirements[[2]]
    [1] "glue"
    
    $Packages$tidyselect$Requirements[[3]]
    [1] "purrr"
    
    $Packages$tidyselect$Requirements[[4]]
    [1] "rlang"
    
    $Packages$tidyselect$Requirements[[5]]
    [1] "vctrs"
    
    
    
    $Packages$tidytext
    $Packages$tidytext$Package
    [1] "tidytext"
    
    $Packages$tidytext$Version
    [1] "0.3.2"
    
    $Packages$tidytext$Source
    [1] "Repository"
    
    $Packages$tidytext$Repository
    [1] "CRAN"
    
    $Packages$tidytext$Hash
    [1] "586ffb833b347003334465682edda644"
    
    $Packages$tidytext$Requirements
    $Packages$tidytext$Requirements[[1]]
    [1] "Matrix"
    
    $Packages$tidytext$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$tidytext$Requirements[[3]]
    [1] "generics"
    
    $Packages$tidytext$Requirements[[4]]
    [1] "hunspell"
    
    $Packages$tidytext$Requirements[[5]]
    [1] "janeaustenr"
    
    $Packages$tidytext$Requirements[[6]]
    [1] "lifecycle"
    
    $Packages$tidytext$Requirements[[7]]
    [1] "purrr"
    
    $Packages$tidytext$Requirements[[8]]
    [1] "rlang"
    
    $Packages$tidytext$Requirements[[9]]
    [1] "stringr"
    
    $Packages$tidytext$Requirements[[10]]
    [1] "tibble"
    
    $Packages$tidytext$Requirements[[11]]
    [1] "tokenizers"
    
    $Packages$tidytext$Requirements[[12]]
    [1] "vctrs"
    
    
    
    $Packages$tidytree
    $Packages$tidytree$Package
    [1] "tidytree"
    
    $Packages$tidytree$Version
    [1] "0.3.9"
    
    $Packages$tidytree$Source
    [1] "Repository"
    
    $Packages$tidytree$Repository
    [1] "CRAN"
    
    $Packages$tidytree$Hash
    [1] "bf9d994b6c4a8e448618feac4bbff2d4"
    
    $Packages$tidytree$Requirements
    $Packages$tidytree$Requirements[[1]]
    [1] "ape"
    
    $Packages$tidytree$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$tidytree$Requirements[[3]]
    [1] "lazyeval"
    
    $Packages$tidytree$Requirements[[4]]
    [1] "magrittr"
    
    $Packages$tidytree$Requirements[[5]]
    [1] "pillar"
    
    $Packages$tidytree$Requirements[[6]]
    [1] "rlang"
    
    $Packages$tidytree$Requirements[[7]]
    [1] "tibble"
    
    $Packages$tidytree$Requirements[[8]]
    [1] "tidyr"
    
    $Packages$tidytree$Requirements[[9]]
    [1] "tidyselect"
    
    $Packages$tidytree$Requirements[[10]]
    [1] "yulab.utils"
    
    
    
    $Packages$tidyverse
    $Packages$tidyverse$Package
    [1] "tidyverse"
    
    $Packages$tidyverse$Version
    [1] "1.3.1"
    
    $Packages$tidyverse$Source
    [1] "Repository"
    
    $Packages$tidyverse$Repository
    [1] "CRAN"
    
    $Packages$tidyverse$Hash
    [1] "fc4c72b6ae9bb283416bd59a3303bbab"
    
    $Packages$tidyverse$Requirements
    $Packages$tidyverse$Requirements[[1]]
    [1] "broom"
    
    $Packages$tidyverse$Requirements[[2]]
    [1] "cli"
    
    $Packages$tidyverse$Requirements[[3]]
    [1] "crayon"
    
    $Packages$tidyverse$Requirements[[4]]
    [1] "dbplyr"
    
    $Packages$tidyverse$Requirements[[5]]
    [1] "dplyr"
    
    $Packages$tidyverse$Requirements[[6]]
    [1] "dtplyr"
    
    $Packages$tidyverse$Requirements[[7]]
    [1] "forcats"
    
    $Packages$tidyverse$Requirements[[8]]
    [1] "ggplot2"
    
    $Packages$tidyverse$Requirements[[9]]
    [1] "googledrive"
    
    $Packages$tidyverse$Requirements[[10]]
    [1] "googlesheets4"
    
    $Packages$tidyverse$Requirements[[11]]
    [1] "haven"
    
    $Packages$tidyverse$Requirements[[12]]
    [1] "hms"
    
    $Packages$tidyverse$Requirements[[13]]
    [1] "httr"
    
    $Packages$tidyverse$Requirements[[14]]
    [1] "jsonlite"
    
    $Packages$tidyverse$Requirements[[15]]
    [1] "lubridate"
    
    $Packages$tidyverse$Requirements[[16]]
    [1] "magrittr"
    
    $Packages$tidyverse$Requirements[[17]]
    [1] "modelr"
    
    $Packages$tidyverse$Requirements[[18]]
    [1] "pillar"
    
    $Packages$tidyverse$Requirements[[19]]
    [1] "purrr"
    
    $Packages$tidyverse$Requirements[[20]]
    [1] "readr"
    
    $Packages$tidyverse$Requirements[[21]]
    [1] "readxl"
    
    $Packages$tidyverse$Requirements[[22]]
    [1] "reprex"
    
    $Packages$tidyverse$Requirements[[23]]
    [1] "rlang"
    
    $Packages$tidyverse$Requirements[[24]]
    [1] "rstudioapi"
    
    $Packages$tidyverse$Requirements[[25]]
    [1] "rvest"
    
    $Packages$tidyverse$Requirements[[26]]
    [1] "stringr"
    
    $Packages$tidyverse$Requirements[[27]]
    [1] "tibble"
    
    $Packages$tidyverse$Requirements[[28]]
    [1] "tidyr"
    
    $Packages$tidyverse$Requirements[[29]]
    [1] "xml2"
    
    
    
    $Packages$timeDate
    $Packages$timeDate$Package
    [1] "timeDate"
    
    $Packages$timeDate$Version
    [1] "3043.102"
    
    $Packages$timeDate$Source
    [1] "Repository"
    
    $Packages$timeDate$Repository
    [1] "CRAN"
    
    $Packages$timeDate$Hash
    [1] "fde4fc571f5f61978652c229d4713845"
    
    $Packages$timeDate$Requirements
    list()
    
    
    $Packages$tinytex
    $Packages$tinytex$Package
    [1] "tinytex"
    
    $Packages$tinytex$Version
    [1] "0.37"
    
    $Packages$tinytex$Source
    [1] "Repository"
    
    $Packages$tinytex$Repository
    [1] "CRAN"
    
    $Packages$tinytex$Hash
    [1] "a80abeb527a977e4bef21873d29222dd"
    
    $Packages$tinytex$Requirements
    $Packages$tinytex$Requirements[[1]]
    [1] "xfun"
    
    
    
    $Packages$tokenizers
    $Packages$tokenizers$Package
    [1] "tokenizers"
    
    $Packages$tokenizers$Version
    [1] "0.2.1"
    
    $Packages$tokenizers$Source
    [1] "Repository"
    
    $Packages$tokenizers$Repository
    [1] "CRAN"
    
    $Packages$tokenizers$Hash
    [1] "a064f646b3a692e62dfb5d9ea690a4ea"
    
    $Packages$tokenizers$Requirements
    $Packages$tokenizers$Requirements[[1]]
    [1] "Rcpp"
    
    $Packages$tokenizers$Requirements[[2]]
    [1] "SnowballC"
    
    $Packages$tokenizers$Requirements[[3]]
    [1] "stringi"
    
    
    
    $Packages$treeio
    $Packages$treeio$Package
    [1] "treeio"
    
    $Packages$treeio$Version
    [1] "1.18.1"
    
    $Packages$treeio$Source
    [1] "Bioconductor"
    
    $Packages$treeio$git_url
    [1] "https://git.bioconductor.org/packages/treeio"
    
    $Packages$treeio$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$treeio$git_last_commit
    [1] "a06b6b3"
    
    $Packages$treeio$git_last_commit_date
    [1] "2021-11-12"
    
    $Packages$treeio$Hash
    [1] "835f0ab27ea0cfcad448397568fdf4d1"
    
    $Packages$treeio$Requirements
    $Packages$treeio$Requirements[[1]]
    [1] "ape"
    
    $Packages$treeio$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$treeio$Requirements[[3]]
    [1] "jsonlite"
    
    $Packages$treeio$Requirements[[4]]
    [1] "magrittr"
    
    $Packages$treeio$Requirements[[5]]
    [1] "rlang"
    
    $Packages$treeio$Requirements[[6]]
    [1] "tibble"
    
    $Packages$treeio$Requirements[[7]]
    [1] "tidytree"
    
    
    
    $Packages$tune
    $Packages$tune$Package
    [1] "tune"
    
    $Packages$tune$Version
    [1] "0.1.6"
    
    $Packages$tune$Source
    [1] "Repository"
    
    $Packages$tune$Repository
    [1] "CRAN"
    
    $Packages$tune$Hash
    [1] "3cabba58cec834867411b81e03acfa7b"
    
    $Packages$tune$Requirements
    $Packages$tune$Requirements[[1]]
    [1] "GPfit"
    
    $Packages$tune$Requirements[[2]]
    [1] "cli"
    
    $Packages$tune$Requirements[[3]]
    [1] "dials"
    
    $Packages$tune$Requirements[[4]]
    [1] "dplyr"
    
    $Packages$tune$Requirements[[5]]
    [1] "foreach"
    
    $Packages$tune$Requirements[[6]]
    [1] "generics"
    
    $Packages$tune$Requirements[[7]]
    [1] "ggplot2"
    
    $Packages$tune$Requirements[[8]]
    [1] "glue"
    
    $Packages$tune$Requirements[[9]]
    [1] "hardhat"
    
    $Packages$tune$Requirements[[10]]
    [1] "lifecycle"
    
    $Packages$tune$Requirements[[11]]
    [1] "parsnip"
    
    $Packages$tune$Requirements[[12]]
    [1] "purrr"
    
    $Packages$tune$Requirements[[13]]
    [1] "recipes"
    
    $Packages$tune$Requirements[[14]]
    [1] "rlang"
    
    $Packages$tune$Requirements[[15]]
    [1] "rsample"
    
    $Packages$tune$Requirements[[16]]
    [1] "tibble"
    
    $Packages$tune$Requirements[[17]]
    [1] "tidyr"
    
    $Packages$tune$Requirements[[18]]
    [1] "vctrs"
    
    $Packages$tune$Requirements[[19]]
    [1] "withr"
    
    $Packages$tune$Requirements[[20]]
    [1] "workflows"
    
    $Packages$tune$Requirements[[21]]
    [1] "yardstick"
    
    
    
    $Packages$tzdb
    $Packages$tzdb$Package
    [1] "tzdb"
    
    $Packages$tzdb$Version
    [1] "0.2.0"
    
    $Packages$tzdb$Source
    [1] "Repository"
    
    $Packages$tzdb$Repository
    [1] "CRAN"
    
    $Packages$tzdb$Hash
    [1] "5e069fb033daf2317bd628d3100b75c5"
    
    $Packages$tzdb$Requirements
    $Packages$tzdb$Requirements[[1]]
    [1] "cpp11"
    
    
    
    $Packages$unbalanced
    $Packages$unbalanced$Package
    [1] "unbalanced"
    
    $Packages$unbalanced$Version
    [1] "2.0"
    
    $Packages$unbalanced$Source
    [1] "Repository"
    
    $Packages$unbalanced$Repository
    [1] "CRAN"
    
    $Packages$unbalanced$Hash
    [1] "ed45de16b03027f36b2bc2e8c84c2533"
    
    $Packages$unbalanced$Requirements
    $Packages$unbalanced$Requirements[[1]]
    [1] "FNN"
    
    $Packages$unbalanced$Requirements[[2]]
    [1] "RANN"
    
    $Packages$unbalanced$Requirements[[3]]
    [1] "doParallel"
    
    $Packages$unbalanced$Requirements[[4]]
    [1] "foreach"
    
    $Packages$unbalanced$Requirements[[5]]
    [1] "mlr"
    
    
    
    $Packages$utf8
    $Packages$utf8$Package
    [1] "utf8"
    
    $Packages$utf8$Version
    [1] "1.2.2"
    
    $Packages$utf8$Source
    [1] "Repository"
    
    $Packages$utf8$Repository
    [1] "CRAN"
    
    $Packages$utf8$Hash
    [1] "c9c462b759a5cc844ae25b5942654d13"
    
    $Packages$utf8$Requirements
    list()
    
    
    $Packages$uuid
    $Packages$uuid$Package
    [1] "uuid"
    
    $Packages$uuid$Version
    [1] "1.0-3"
    
    $Packages$uuid$Source
    [1] "Repository"
    
    $Packages$uuid$Repository
    [1] "CRAN"
    
    $Packages$uuid$Hash
    [1] "2097822ba5e4440b81a0c7525d0315ce"
    
    $Packages$uuid$Requirements
    list()
    
    
    $Packages$vctrs
    $Packages$vctrs$Package
    [1] "vctrs"
    
    $Packages$vctrs$Version
    [1] "0.3.8"
    
    $Packages$vctrs$Source
    [1] "Repository"
    
    $Packages$vctrs$Repository
    [1] "CRAN"
    
    $Packages$vctrs$Hash
    [1] "ecf749a1b39ea72bd9b51b76292261f1"
    
    $Packages$vctrs$Requirements
    $Packages$vctrs$Requirements[[1]]
    [1] "ellipsis"
    
    $Packages$vctrs$Requirements[[2]]
    [1] "glue"
    
    $Packages$vctrs$Requirements[[3]]
    [1] "rlang"
    
    
    
    $Packages$vegan
    $Packages$vegan$Package
    [1] "vegan"
    
    $Packages$vegan$Version
    [1] "2.5-7"
    
    $Packages$vegan$Source
    [1] "Repository"
    
    $Packages$vegan$Repository
    [1] "CRAN"
    
    $Packages$vegan$Hash
    [1] "01771d8de354fa30c87c68e599e2429b"
    
    $Packages$vegan$Requirements
    $Packages$vegan$Requirements[[1]]
    [1] "MASS"
    
    $Packages$vegan$Requirements[[2]]
    [1] "cluster"
    
    $Packages$vegan$Requirements[[3]]
    [1] "lattice"
    
    $Packages$vegan$Requirements[[4]]
    [1] "mgcv"
    
    $Packages$vegan$Requirements[[5]]
    [1] "permute"
    
    
    
    $Packages$vipor
    $Packages$vipor$Package
    [1] "vipor"
    
    $Packages$vipor$Version
    [1] "0.4.5"
    
    $Packages$vipor$Source
    [1] "Repository"
    
    $Packages$vipor$Repository
    [1] "CRAN"
    
    $Packages$vipor$Hash
    [1] "ea85683da7f2bfa63a98dc6416892591"
    
    $Packages$vipor$Requirements
    list()
    
    
    $Packages$viridis
    $Packages$viridis$Package
    [1] "viridis"
    
    $Packages$viridis$Version
    [1] "0.6.2"
    
    $Packages$viridis$Source
    [1] "Repository"
    
    $Packages$viridis$Repository
    [1] "CRAN"
    
    $Packages$viridis$Hash
    [1] "ee96aee95a7a563e5496f8991e9fde4b"
    
    $Packages$viridis$Requirements
    $Packages$viridis$Requirements[[1]]
    [1] "ggplot2"
    
    $Packages$viridis$Requirements[[2]]
    [1] "gridExtra"
    
    $Packages$viridis$Requirements[[3]]
    [1] "viridisLite"
    
    
    
    $Packages$viridisLite
    $Packages$viridisLite$Package
    [1] "viridisLite"
    
    $Packages$viridisLite$Version
    [1] "0.4.0"
    
    $Packages$viridisLite$Source
    [1] "Repository"
    
    $Packages$viridisLite$Repository
    [1] "CRAN"
    
    $Packages$viridisLite$Hash
    [1] "55e157e2aa88161bdb0754218470d204"
    
    $Packages$viridisLite$Requirements
    list()
    
    
    $Packages$vroom
    $Packages$vroom$Package
    [1] "vroom"
    
    $Packages$vroom$Version
    [1] "1.5.7"
    
    $Packages$vroom$Source
    [1] "Repository"
    
    $Packages$vroom$Repository
    [1] "CRAN"
    
    $Packages$vroom$Hash
    [1] "976507b5a105bc3bdf6a5a5f29e0684f"
    
    $Packages$vroom$Requirements
    $Packages$vroom$Requirements[[1]]
    [1] "bit64"
    
    $Packages$vroom$Requirements[[2]]
    [1] "cli"
    
    $Packages$vroom$Requirements[[3]]
    [1] "cpp11"
    
    $Packages$vroom$Requirements[[4]]
    [1] "crayon"
    
    $Packages$vroom$Requirements[[5]]
    [1] "glue"
    
    $Packages$vroom$Requirements[[6]]
    [1] "hms"
    
    $Packages$vroom$Requirements[[7]]
    [1] "lifecycle"
    
    $Packages$vroom$Requirements[[8]]
    [1] "progress"
    
    $Packages$vroom$Requirements[[9]]
    [1] "rlang"
    
    $Packages$vroom$Requirements[[10]]
    [1] "tibble"
    
    $Packages$vroom$Requirements[[11]]
    [1] "tidyselect"
    
    $Packages$vroom$Requirements[[12]]
    [1] "tzdb"
    
    $Packages$vroom$Requirements[[13]]
    [1] "vctrs"
    
    $Packages$vroom$Requirements[[14]]
    [1] "withr"
    
    
    
    $Packages$waldo
    $Packages$waldo$Package
    [1] "waldo"
    
    $Packages$waldo$Version
    [1] "0.3.1"
    
    $Packages$waldo$Source
    [1] "Repository"
    
    $Packages$waldo$Repository
    [1] "CRAN"
    
    $Packages$waldo$Hash
    [1] "ad8cfff5694ac5b3c354f8f2044bd976"
    
    $Packages$waldo$Requirements
    $Packages$waldo$Requirements[[1]]
    [1] "cli"
    
    $Packages$waldo$Requirements[[2]]
    [1] "diffobj"
    
    $Packages$waldo$Requirements[[3]]
    [1] "fansi"
    
    $Packages$waldo$Requirements[[4]]
    [1] "glue"
    
    $Packages$waldo$Requirements[[5]]
    [1] "rematch2"
    
    $Packages$waldo$Requirements[[6]]
    [1] "rlang"
    
    $Packages$waldo$Requirements[[7]]
    [1] "tibble"
    
    
    
    $Packages$warp
    $Packages$warp$Package
    [1] "warp"
    
    $Packages$warp$Version
    [1] "0.2.0"
    
    $Packages$warp$Source
    [1] "Repository"
    
    $Packages$warp$Repository
    [1] "CRAN"
    
    $Packages$warp$Hash
    [1] "2982481615756e24e79fee95bdc95daa"
    
    $Packages$warp$Requirements
    list()
    
    
    $Packages$webshot
    $Packages$webshot$Package
    [1] "webshot"
    
    $Packages$webshot$Version
    [1] "0.5.2"
    
    $Packages$webshot$Source
    [1] "Repository"
    
    $Packages$webshot$Repository
    [1] "CRAN"
    
    $Packages$webshot$Hash
    [1] "e99d80ad34457a4853674e89d5e806de"
    
    $Packages$webshot$Requirements
    $Packages$webshot$Requirements[[1]]
    [1] "callr"
    
    $Packages$webshot$Requirements[[2]]
    [1] "jsonlite"
    
    $Packages$webshot$Requirements[[3]]
    [1] "magrittr"
    
    
    
    $Packages$withr
    $Packages$withr$Package
    [1] "withr"
    
    $Packages$withr$Version
    [1] "2.5.0"
    
    $Packages$withr$Source
    [1] "Repository"
    
    $Packages$withr$Repository
    [1] "CRAN"
    
    $Packages$withr$Hash
    [1] "c0e49a9760983e81e55cdd9be92e7182"
    
    $Packages$withr$Requirements
    list()
    
    
    $Packages$workflows
    $Packages$workflows$Package
    [1] "workflows"
    
    $Packages$workflows$Version
    [1] "0.2.4"
    
    $Packages$workflows$Source
    [1] "Repository"
    
    $Packages$workflows$Repository
    [1] "CRAN"
    
    $Packages$workflows$Hash
    [1] "a087bd9bcc22f572730a8949a2ab4a47"
    
    $Packages$workflows$Requirements
    $Packages$workflows$Requirements[[1]]
    [1] "cli"
    
    $Packages$workflows$Requirements[[2]]
    [1] "ellipsis"
    
    $Packages$workflows$Requirements[[3]]
    [1] "generics"
    
    $Packages$workflows$Requirements[[4]]
    [1] "glue"
    
    $Packages$workflows$Requirements[[5]]
    [1] "hardhat"
    
    $Packages$workflows$Requirements[[6]]
    [1] "lifecycle"
    
    $Packages$workflows$Requirements[[7]]
    [1] "parsnip"
    
    $Packages$workflows$Requirements[[8]]
    [1] "rlang"
    
    $Packages$workflows$Requirements[[9]]
    [1] "tidyselect"
    
    $Packages$workflows$Requirements[[10]]
    [1] "vctrs"
    
    
    
    $Packages$workflowsets
    $Packages$workflowsets$Package
    [1] "workflowsets"
    
    $Packages$workflowsets$Version
    [1] "0.1.0"
    
    $Packages$workflowsets$Source
    [1] "Repository"
    
    $Packages$workflowsets$Repository
    [1] "CRAN"
    
    $Packages$workflowsets$Hash
    [1] "34d9ca8086d92ee78df6747831ce3ec5"
    
    $Packages$workflowsets$Requirements
    $Packages$workflowsets$Requirements[[1]]
    [1] "cli"
    
    $Packages$workflowsets$Requirements[[2]]
    [1] "dplyr"
    
    $Packages$workflowsets$Requirements[[3]]
    [1] "generics"
    
    $Packages$workflowsets$Requirements[[4]]
    [1] "ggplot2"
    
    $Packages$workflowsets$Requirements[[5]]
    [1] "hardhat"
    
    $Packages$workflowsets$Requirements[[6]]
    [1] "lifecycle"
    
    $Packages$workflowsets$Requirements[[7]]
    [1] "prettyunits"
    
    $Packages$workflowsets$Requirements[[8]]
    [1] "purrr"
    
    $Packages$workflowsets$Requirements[[9]]
    [1] "rlang"
    
    $Packages$workflowsets$Requirements[[10]]
    [1] "rsample"
    
    $Packages$workflowsets$Requirements[[11]]
    [1] "tibble"
    
    $Packages$workflowsets$Requirements[[12]]
    [1] "tidyr"
    
    $Packages$workflowsets$Requirements[[13]]
    [1] "tune"
    
    $Packages$workflowsets$Requirements[[14]]
    [1] "vctrs"
    
    $Packages$workflowsets$Requirements[[15]]
    [1] "withr"
    
    $Packages$workflowsets$Requirements[[16]]
    [1] "workflows"
    
    
    
    $Packages$xfun
    $Packages$xfun$Package
    [1] "xfun"
    
    $Packages$xfun$Version
    [1] "0.30"
    
    $Packages$xfun$Source
    [1] "Repository"
    
    $Packages$xfun$Repository
    [1] "CRAN"
    
    $Packages$xfun$Hash
    [1] "e83f48136b041845e50a6658feffb197"
    
    $Packages$xfun$Requirements
    list()
    
    
    $Packages$xml2
    $Packages$xml2$Package
    [1] "xml2"
    
    $Packages$xml2$Version
    [1] "1.3.3"
    
    $Packages$xml2$Source
    [1] "Repository"
    
    $Packages$xml2$Repository
    [1] "CRAN"
    
    $Packages$xml2$Hash
    [1] "40682ed6a969ea5abfd351eb67833adc"
    
    $Packages$xml2$Requirements
    list()
    
    
    $Packages$xtable
    $Packages$xtable$Package
    [1] "xtable"
    
    $Packages$xtable$Version
    [1] "1.8-4"
    
    $Packages$xtable$Source
    [1] "Repository"
    
    $Packages$xtable$Repository
    [1] "CRAN"
    
    $Packages$xtable$Hash
    [1] "b8acdf8af494d9ec19ccb2481a9b11c2"
    
    $Packages$xtable$Requirements
    list()
    
    
    $Packages$yaml
    $Packages$yaml$Package
    [1] "yaml"
    
    $Packages$yaml$Version
    [1] "2.3.5"
    
    $Packages$yaml$Source
    [1] "Repository"
    
    $Packages$yaml$Repository
    [1] "CRAN"
    
    $Packages$yaml$Hash
    [1] "458bb38374d73bf83b1bb85e353da200"
    
    $Packages$yaml$Requirements
    list()
    
    
    $Packages$yardstick
    $Packages$yardstick$Package
    [1] "yardstick"
    
    $Packages$yardstick$Version
    [1] "0.0.9"
    
    $Packages$yardstick$Source
    [1] "Repository"
    
    $Packages$yardstick$Repository
    [1] "CRAN"
    
    $Packages$yardstick$Hash
    [1] "fd1588dbcbb85aacd0d9b9009351c79d"
    
    $Packages$yardstick$Requirements
    $Packages$yardstick$Requirements[[1]]
    [1] "dplyr"
    
    $Packages$yardstick$Requirements[[2]]
    [1] "generics"
    
    $Packages$yardstick$Requirements[[3]]
    [1] "pROC"
    
    $Packages$yardstick$Requirements[[4]]
    [1] "rlang"
    
    $Packages$yardstick$Requirements[[5]]
    [1] "tidyselect"
    
    $Packages$yardstick$Requirements[[6]]
    [1] "vctrs"
    
    
    
    $Packages$yulab.utils
    $Packages$yulab.utils$Package
    [1] "yulab.utils"
    
    $Packages$yulab.utils$Version
    [1] "0.0.4"
    
    $Packages$yulab.utils$Source
    [1] "Repository"
    
    $Packages$yulab.utils$Repository
    [1] "CRAN"
    
    $Packages$yulab.utils$Hash
    [1] "922e11dcf40bb5dfcf3fe5e714d0dc35"
    
    $Packages$yulab.utils$Requirements
    list()
    
    
    $Packages$zlibbioc
    $Packages$zlibbioc$Package
    [1] "zlibbioc"
    
    $Packages$zlibbioc$Version
    [1] "1.40.0"
    
    $Packages$zlibbioc$Source
    [1] "Bioconductor"
    
    $Packages$zlibbioc$git_url
    [1] "https://git.bioconductor.org/packages/zlibbioc"
    
    $Packages$zlibbioc$git_branch
    [1] "RELEASE_3_14"
    
    $Packages$zlibbioc$git_last_commit
    [1] "3f116b3"
    
    $Packages$zlibbioc$git_last_commit_date
    [1] "2021-10-26"
    
    $Packages$zlibbioc$Hash
    [1] "3598bf6c766b0d2c15d1eefefa26d9ef"
    
    $Packages$zlibbioc$Requirements
    list()
    
    
    
    attr(,"class")
    [1] "renv_lockfile"
    
    $synchronized
    [1] TRUE
    



```R
library(dplyr)
```

    
    Attaching package: ‘dplyr’
    
    
    The following objects are masked from ‘package:data.table’:
    
        between, first, last
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    



```R

```

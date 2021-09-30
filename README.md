# Reproducible analyses for untitled manuscript on trait-based taxon sets.  
# Quang Nguyen   

Analyses performed using the `targets` package in R. Access `run.R` to run individual pipelines.    

Individual pipelines include:   
* `script-dada2.R` defines targets for processing data from Gevers et al. using the `dada2` manuscript.  


## Installing the required python packages  

This project uses `reticulate` package to call python functions in R. This is linked to a virtual environment that is managed by `renv` using `requirements.txt`.  To ensure the best usage, please edit `.Renviron` file to force `RETICULATE_PYTHON` to be at the correct virtual environment. This will also correct the behavior of `renv` to point to the correct python binaries.   

In order to conveniently access the MetaCyc pathway class hierarchy, the package `pythoncyc` is used. Since this package is not distributed via `pip` or `conda`, please create a new conda or virtual environment, and then use the [manual installation instructions](https://github.com/ecocyc/PythonCyc) when your desired virtual environment is activated.  
 


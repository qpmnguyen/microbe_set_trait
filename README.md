# Reproducible analyses for untitled manuscript on trait-based taxon sets.  
# Quang Nguyen   

## R pipelines  
Most of the analyses was performed using R and the `targets` package. The entire pipeline was broken down into smaller "steps" for ease of processing. The file `run.R` contains code to choose between different pipelines to run. 

Individual mini-pipelines include:   
* `script-dada2.R` defines targets for processing raw data from Gevers et al. into ASV tables using the `dada2` method by Callahan et al. 2016.  


## Installing the required python packages  

Certain parts of this project requires using python functions. Inter-operability with R is maintained using the `reticulate` package to call python functions in the R environment. Due to certain issues with `reticulate` and pre-specified conda environments, we suggest to create your own conda environment using specifications listed in `environment.yml`, and explicitly declare the environment variable `RETICULATE_PYTHON` to point towards the Python 3 executable associated with the conda environment. This will also correct the behavior of `renv` to point to the correct python version and associated library.  

In order to conveniently access the MetaCyc pathway class hierarchy, the package `pythoncyc` is used. Since this package is not distributed via `pip` or `conda`, please create a new conda or virtual environment, and then use the [manual installation instructions](https://github.com/ecocyc/PythonCyc) when your desired conda environment is activated. Due to the perks of the `pythoncyc` library, users need to also install the associated `PathwayTools` program (v. 25.1) and run it with the python API open (see instructions in the link above). This step doesn't need to be repeated (hence not included in the major pipelines) and the product of this step is saved as a binary file called `metacyc_parse.rds` in the `databases` folder.  
 


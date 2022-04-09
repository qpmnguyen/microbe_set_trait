# Reproducible analyses for untitled manuscript on trait-based taxon sets.  
# Quang Nguyen   

## Dependency management 

All R analyses were performed under R version 4.1.2 and the corresponding Bioconductor version 3.14. Specific package versions are managed via `renv` and details can be found in the corresponding `renv.lock` file.  

All python analyses are performed under Python version 3.8.9. Specific package versions are managed via `Miniconda`. The specifications for the primary virtual environment is in the `microbe_trait_env.yml` file. Since we also use `QIIME2` and `picrust2`, separate conda environments were utilized for each software respectively due to dependency conflicts. Specifications for the `QIIME2` 2022.2 environment can be found in the `qiime2-2022.2-py38-linux-conda.yml` file (for different OS, please refer to the main `QIIME2` documentation for instructions to install version 2022.2 in your machine).  

## Analyses  

All analyses are in Jupyter Notebooks located in the `analysis` folder.   
-  Placeholder  

## Pipelines outside of jupyter notebooks  

Most analyses are located within reproducible jupyter notebooks. However, there are some computations that require more resources and therefore are ran as stand-alone scripts or as part of more reproducible workflows.    

#### R 
Intensive computing steps in R were performed using the `targets` package. The entire pipeline was broken down into smaller "steps" for ease of processing. The file `run.R` contains code to choose between different pipelines to run. 

Individual mini-pipelines include:   
- Placeholder

## Function files    

#### R
- `metaphlan_db.R` defines functions to process the MetaPhlAn3 marker information files.  
- `db_preprocess.R` defines functions to pre-process the database into `BiocSet` format.  
- `hmp_preprocess.R` defines functions to pre-process the HMP data data files (newest release).   

#### Python
- Placeholder 

## Additional notes on reproducibility:      

In order to conveniently access the MetaCyc pathway class hierarchy, the package `pythoncyc` is used. Since this package is not distributed via `pip` or `conda`, please create a new conda or virtual environment, and then use the [manual installation instructions](https://github.com/ecocyc/PythonCyc) when your desired conda environment is activated. Due to the perks of the `pythoncyc` library, users need to also install the associated `PathwayTools` program (v. 25.1) and run it with the python API open (see instructions in the link above). This step doesn't need to be repeated (hence not included in the major pipelines) and the product of this step is saved as a binary file called `metacyc_parse.rds` in the `databases` folder.  
 


# GWAS_QA_QC_pipeline

### Description
----------------
This suite has two stand-alone scripts to automate general quality control and data clean up in genotyping data used in GWAS.  The first program, located in the data_visualization_stats folder with the executable called basic_stats.py, is used to help the user determine what parameters should be used in the actual filtering pipeline, located in the folder filtering_pipeline. The executable called sample_snp_qc.py, can take in a set of parameters and will actually filter out samples and variants that do not meet the specified thresholds input by the user.  A set of default parameters is provided if the user does not specify any parameters.  Please see README.md in the filtering_pipeline directory for more information.

### General Overview
-------------
![Alt text](https://github.com/tbrunetti/GWAS_QA_QC_pipeline/blob/master/GWAS_pipeline_image_workflow.jpg)


### Software requirements
-------------------------
For specific software requirements for each program please see their corresponding README.md files.  If the purpose is to use all modules in the GWAS_QA_QC_pipeline the following software requirements must be met:
* Python version 2.7
* PLINK version 1.9

__Python Dependencies__
All python dependencies can be installed with pip or from source at the html links provided
* matplotlib, numpy, pandas (https://www.scipy.org/install.html)
* chunkypipes (https://pypi.python.org/pypi/ChunkyPipes)
* seaborn (https://pypi.python.org/pypi/seaborn)
* statistics (https://pypi.python.org/pypi/statistics)

### User Generated/User Provided File Requirements
--------------------------------------------------


### Installation and Configuration
-----------------------------------



### Running the Pipeline
-------------------------


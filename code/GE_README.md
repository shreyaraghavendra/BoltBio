# Code Description

## Overview

The code downloads data from TCGA using the ```gdc-client``` API and creates the expression table.

First run 
```
conda env create -f env.yml
```
to setup the environment 


I am following the paper ```https://www.mdpi.com/2072-6694/12/12/3799#``` [Code : https://github.com/fvalle1/topicTCGA/tree/master/lung]. 

*Note Through primarily I should be using LUAD data, the model in the paper classifies a given sample as LUAD or LUSC postive and hence I am downloading LUSC data as well. That does not alter the relevant analysis and can be changed easily*

The ```getDATA.ipynb``` creates the gdc-client complaint manifest file according to the parameters, such as Workflow_type, data_format, case_id, etc. required. Once the ```manifest.txt``` file has been created, run the following command in the terminal 

```
mkdir -p data && mv manifest.txt data/. && cd data
gdc-client download -m manifest.txt
```
This will create the data folder and then download the files

Now run through the ```getTable.ipynb``` notebook to create the expression table saved as ```mainTable_all.csv```

Then run the ```preprocess.ipynb``` file. In the context of the paper, the selection of highly variable genes is a preprocessing step and hence in our context this is not exactly preprocessing. In this notebook all the related files are assimilated, including the downloaded files metadata and a final table is created. The top 3000 variable genes are then sampled with scanpy. in addition to operations for the paper, I also run a few differential expression tests on the dataset and the results can be seen in the jupyter notebook. The top 25 differentialy expressed genes are stored and their scores and ids can be accessed from the created anndata object. This notebook contains most of the code relevant to our pipeline.

Note : Remeber to change the path to files according to your working directory. The variables have been decalred at the beginning of each notebook






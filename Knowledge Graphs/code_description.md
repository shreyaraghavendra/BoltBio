# Code Description

## Overview

The code performs Link Prediction on Hetionet.

The ```Hetionet-Json``` file must provided in the ```data``` subdirectory. The user settings such as the ```BATCH_SIZE``` or ```HYPERPARAMETERS``` in the ```config.yaml``` file.
The directory must be structured as follows:

```
\
|_ data\
        |_ hetionet-v1.0.json
|_ train.py
|_config.yaml
```

The user simply needs to run the following command to process the graph and obtain results
```
python3 train.py -c config.yaml
```

The user can specify if they want to use trackers like ```wandb``` to track results.

The ```graph.py``` file contains code to process the Hetionet json file and returns the graph and subgraphs along with the ids of every single node in the graph.

```models.py``` contains the model architeture being used.




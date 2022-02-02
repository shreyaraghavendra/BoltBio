# BoltBio

### Problem Statement

Diseases with the highest fatality rates such as heart disease, diabetes, cancer, asthma, etc. and many untreatable diseases such as Alzheimer’s are highly heterogeneous, with many etiologies for a single phenotype. This makes biomarker discovery very challenging. In the last decade, the application of different individual omic studies (e.g., genomics, epigenomics, transcriptomics, proteomics, metagenomics) to understand a particular problem in human disease have been successful to a great extent. However, Multi-omics data generated for the same set of samples can also represent the flow of biological information at multiple levels and thus provide a holistic system view of the disease in relation to the biological system of interest. Bioinformatics methods and Artificial Intelligence models such as Graph Convolutional Networks can be very useful to integrate datasets with very different data modalities originating from varied assay types and increased dimensionality. Systems biology and network analysis can leverage interaction networks in trans-omics to significantly replicate the disease ecosystem.

<p align="center" >
  <img src="https://user-images.githubusercontent.com/30191888/151766060-ca512a82-cef5-4eda-aaa6-e742f6ecbd30.png" style="width:50%;height:50%;"/>
</p>


### Introduction

BoltBio is a package for gene target prioritization in diseases using multi-omics data. The package will consist of individual omic analysis pipelines for mutation (genomics), gene expression (transcriptomics), DNA Methylation (epigenomics), miRNA (epigenomics) and pathway (interactomics) data and an integrated pipeline for multi-omics analysis. Each model will prioritize a list of gene biomarkers, the results of which will be compared and correlated. The user can provide one or any combination of aforementioned datasets and use the gene rankings from the corresponding pipelines to choose an optimum target. Chosen targets can then be tested in wetlabs.





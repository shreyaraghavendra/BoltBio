### **Summary of Analysis Workflow**

- pyGNA aims to identify two things: whether the genes witin a geneset form a pathway and if there's associaton between genes in two different genesets using gene network topology tests(GNA)  and gene network assosication (GNA) tests respectively.
- It uses three models--direct, shortest path(SP), Random Walk with Restart(RWR), in the increasing order of robustness.
- Weights and direction aren't taken into consideration in this analysis, but can be implemented using READDATA abstract class.
- It takes in input genesets in GMT format, and has commands to convert CSV, deseq table for GMT conversion, and networks are imported as TSV file. The RWR and SP matrices are generated using the *build* function and guring GNA and GNT anaysis, genes in genesets are mapped to the network and the assoiation is calculated. In order to make results comparable, observed test statistics are transformed in normalised z-scores. 
- Diagnostic distribution plot can be obtained by passing -d <diagnost folder/> which gives us the null empirical distribution (blue) and observed value (red bar) for each test. A null distribution can be computed either by sampling two random genesets or by sampling only one of the two in case of GNA--the latter is more conservative and is recommended when checking for association with known pathways.
- Visualisation is done through bar plots, dot plots, heat maps, and volcano plots and top 5 terms can b annotated. 
- GNT is done between the genes in a geneset an GNA does the same among genes in two differet geneets. In case there is no second geneset file, set names can be specified or it's run within the geneets in the file. _paint-comparison-matrix_ computes heatmap by default of the table that has the GNA resuts. 
- This can be carried out using Jupyter Notebook or Snakemake workflow. In case  of Snakemake workflow, there are three important files that can be tweaked accordingly to obtain results for different dataset: the config file, snakemake file, and the rules file. Rules file contains the inormation to generate the results(something similar to modules), and config file contains the informtion about what the parameters are and where to find the input files,and snakefile is where the results are actully generated. In this case, the csv file file for each cancer type is stored in seperate folders in the results folder and the folders are specified in a list and passed as parameter for the generation of a gmt file containing all the genesets. Depening on the results one would like to generate, the rules can be customised. There is an R script to generate the csv file. This can be done for any TCGA tisue type, but for other datasets, a different script has to be written to generate the data for creating the CSV file which is subsequently used to create the gmt file.
 
### **Algorithms**

These three interaction models capture different topological properties. Direct models provide information about the neighborhood of a gene and its observed links. However, they might not be sufficiently powered to detect mid- and long-range interactions, thus statistics defined under these models are usually sensitive to missing links. Conversely, modelling gene interactions using shortest path provides a simple analytical framework to include local and global awareness of the connectivity. However, this approach is also sensitive to missing links and small-world effects, which is common in biological networks and could lead to false positives. The RWR model is more robust than the shortest path model, because it effectively adjusts interaction effects for network structure; it rewards nodes connected with many shortest paths, and penalizes those that are connected only by path going through high degree nodes.

The the strength of interaction between genes in the geneset (geneset network topology, GNT) and with genes in another geneset (geneset network association, GNA) are evaluated based on statistics using these model. 

#### GNT Statistics:

Let S=s1,...,sn be a geneset of n genes, each mapped to a node in G=(V,E). The aim is to find out whether the strength of interaction between nodes of the geneset is higher than expected by chance for a geneset of the same size.

-   Total degree statistic: The importance of a geneset S can be quantified as the number of edges connecting each node in S to any other node in the network

    ![](https://lh4.googleusercontent.com/TFBxV9JakouiKx5a1UpJZcm30tI2cJot38B2D9jIsryM76_82NAG16gkJ9yus5BbxBAi_VydD-Eogi0HfhofssBRaLnHhkD57hdpny5CqauiMDCJ381CtquDXBwMFILPF66XiDjl)

-   Internal degree: With the direct interaction model, the strength of interaction for a geneset S can be quantified as the number of edges connecting each node in S to any other node in the geneset; we refer to this quantity as the internal degree of the node
    ![](https://lh3.googleusercontent.com/sWv5p9S3ZWHEqAN23L-VEj6UQSqHK7n-9TrZuktzsuorEgVMoyZ1F1axdYG3D5CxppMKx60vVnmFUWrs62iOdGpTGLY0g_Qe8LGnGzwTP-ntX640Kq_6UTkIPTwKPB4oLCeRScWr)
-   Shortest Path: The main concern regarding direct interaction methods is that they could fail in presence of missing links, which can be circumvented by shortest path method which takes into account the distance between nodes. Shortest path is the average of the minimum distance between each gene and the rest of those in S:

  ![](https://lh6.googleusercontent.com/0aAzF-yJwTV_GyQaCUti3iXS9tC9PSpB5zIukWBo-suYxs6AcOqVjpeggTJAr4y2pNtQRR3H41l9zfKQGUNzEMurSf5hMyNjQy7Ndi8kCQDKyFkepq2WJemiNO5vxfzxC7mUElIQ)

-   Under a RWR model, we can consider hij∈H as the heat transferred from node i to node j, which can be used as a measure of interaction strength between the nodes in the geneset S, as follows:

   ![](https://lh3.googleusercontent.com/R29kXbSK4gGjzcTDGOrv5vZX9vPOlv5LNjqBDh7cAyNZt4JJm2vpi6PgCiCP89SyjSg3EDjQWnfVm36MOoM2AqKzbNfRMK0zV6JxI9hknWgqDrUxkNLs0LwNmDOlDXP3KchLOPrh)

#### Geneset network association statistics

Let S1 and S2 be two geneset with n and m genes respectively, we want to estimate the association between S1 and S2 as a function of the strength of interaction between their nodes.

Under a shortest path model, the association statistics USP is defined as follows:

Whereas, under a RWR model, we measure association as a function of the heat, UH, transferred between the two genesets as follows:

![](https://lh4.googleusercontent.com/pNFrWarxqb1BktJC4aWQCBtDERjy0DshwL03xuSeeanyFcNeGjuixsXbQtnyoZZsTqtqVOvz4zqf_CQaKYuzwlRczG9XFADGZ1bXnUXgN7xcxBC5Zga6wC9FccWHZuoPGsrMFTKQ)

Where the heat withheld by a gene is considered when there are overlapping genes between S1 and S2.
![](https://lh6.googleusercontent.com/9WUqYLEtoYJjwofcsQPPfIPv0xlfGsmODfhA7ETtZmtQ0ZewVP_oOnmwBmmTaQd2a0lTchhJtSuDcqoQ2nL1--UqAe866plxRMfNHI88Eg1KOBQIbPTznpyi-Hw4Q48OgiwGxYzl)

#### Hypothesis Testing

The topological and association statistics are ultimately used for hypothesis testing. A calibrated null distribution is needed to estimate whether the observed statistics are more extreme than what expected by chance using which, a bootstrap procedure is used to estimate null distributions of the test statistics, conditioned on the geneset size.

Thus, w.l.o.g, let Q be the null distribution of the test statistic q estimated for a geneset of size n, and ![](https://lh4.googleusercontent.com/XJeK67eNboaq0tH-m5YCWtbET8SYrlDwEE559PatIXDWxggXY3WwwJPgagrYaRvQECK60wSYx5sPIO_PTjRAtNy9dvx51BEQ1Kaa-B8ObW4yU3agKvjWESljd9pArBnCAGtuEEdv) the observed value. It is possible to derive an empirical p-value as follows:

 ![](https://lh6.googleusercontent.com/dU4jhvcZn0hQBFLBkB50iFBS_OVim_dfuZz8pqUPAMrejNFxl7lqcK3ques3HUgqlcPVcUb8soPwDMlQwZwXuxjNKdzf1axvkfcCwEVbrBTqdqBYjRDt5fF-n7fqTBM0ULZ7Enj8)

where I is the indicator function returning 1 if and only if the evaluated condition i

#### Benchmarking

**Stochastic block models (SBM) was used for both GNA and GNT.** 

SBM framework is used to simulate a network with k blocks, with a baseline probability of interaction within and between blocks, p0.

GNT:  k+<k blocks from the SBM matrix are selected  and set their within probability of connection M+ii = αp0, where α>1 is a scaling factor controlling the strength of interaction of the genes within block i compared to the rest of the genes in any other block.

GNA: k+  blocks are selected at random and set their within block connection probability to M+ii=αp0 and their between blocks connection to M+ij=γp0 for i≠j and α,γ>1. γ  was reparameterized as a function α, in order to control the relationship between the within and between block connection probability. 

Let β=γ/(α-1), the between block connection probability can be set as M+ij=p0+βp0(α-1). With this parametrization, 3 different scenarios can be simulated:

1)if β=0⇒M+ij=p0, the connection probability between blocks is equal to the baseline, thus genes in a block are highly connected.

2)if 0<β<1⇒p0<M+ij<M+ii, then the connection probability between the blocks is higher than the baseline, and thus we obtain assortative genesets.

3)if β>1⇒M+ij>M+ii, then we have non assortative genesets, thus we expect them to be detected by a GNA test.

After building a network, genesets can be generated by selecting two distinct blocks, i, j, with m nodes each, and add π×m nodes from block i and (1-π)×m nodes from block j; for simplicity, we picked genes from blocks containing the same number of genes. 

By varying the size of highly connected blocks and their interaction probability, along with the geneset composition, it is possible to assess the true positive rate (TPR) and false positive rate (FPR) of both GNT and GNA tests.

**High degree nodes model for GNT benchmarking**

The high degree nodes (HDN) model generates networks with a controllable number of hubs, nhd, whose probability of connection with another node, phd, is higher than the baseline probability p0 assigned to any other node in the network. The model is fully specified by four parameters, namely the number of nodes in the network, n, the number of HDN nodes, nhd, the baseline connection probability, p0, and the HDN connection probability, phd>p0.

In order to benchmark GNT tests in presence of HDN nodes, we created geneset as a mixture of HDNs and non HDN nodes; we denoted these genesets as extended genesets. Specifically, each geneset is made of πhd×nhd nodes, with πhd∈(0,1], and ρπhdnhd random high degree nodes, where ρ is the ratio between high degree nodes and other nodes in the network (see Additional file 1 for a graphical representation).

With the HDN model, we can replicate a common scenario where the tested geneset is made of a few master regulators and many, possibly, unrelated genes. Here, the idea is that a robust GNT test should have a low false positive rate, even when observed statistics might be skewed by few highly connected nodes.


### **Results**

**SBM and HDN results**
TID, TH and TSP are the statistics with the best overall performances with TPR >70% for all instances, whereas TM and TTD were able to detect a network effects only for highly connected genesets. In general, all tests are robust to false positives (FPR<10% for all tests), with TSP being the most conservative.
HDN model:Low FPR (<10%) regardless of geneset composition and network structure. Interestingly, while TSP was the most robust on SBM networks, it is the most sensitive to HDN in the networks, with FPR as high as 20% even for genesets with only 3 HDNs In this case TID is the most robust test (FPR<10%), while TH has FPR>0.2 when the number of HDN increases.
Taken together, the TID statistic is the one achieving the best performances and it is faster to compute respect to the other best performer, TH, which requires the computation of a random walk matrix. Nonetheless, for exploratory analyses, using the TH test is recommended, which is confirmed to be well powered to detect network effects and has a low FPR, and might less sensitive to missing links. 

**GNT:** Only differentially expressed genes in lung and lymphoid cancers show a significant network effect, albeit this is detected by all tests for lung cancer and only by TH and TID lymphoid neoplasm. 

Network effect wasn't observed for the other cancers; this could be explained by the fact these cancers might be controlled not by one highly connected network, but by multiple distinct ones.

PyGNA diagnostic plot generated by the GNT analysis detected network effect.

**GNA:** While most of them do not seem to show a consistent network effect, GNA can be used to test whether each set of differentially expressed genes are more connected with each other than expected by chance.Using either UH and USP tests, a significant association between breast, bladder, and prostate carcinomas, and between leukemia and lymphoid neoplasms was observed which is consistent with other gene expression analyses, that have shown that anatomically related cancers or with similar histopathology share similar changes in gene expression.Significant association between lung and lymphoid neoplasms observed in the study might be explained by the fact that lungs contain a vast lymphatic network, which might also be dysregulated in lung tumors.

**Conclusion:** PyGNA enables network analysis of RNA sequencing datasets and provide useful biological insights about the association between the genes.


### **Observations**

- The computation of SP matrix takes over 3-4 hours and doesn't get completed (and 2hrs on Jupyter Notebook but the matrix always gives error when used).
- Snakemake workflow seems seamless which could be due to the default number of permutations being very small(in this case 2). Analyses that do not involve rwr or SP matrix (topology module, internal degree, total degree) do not take a lot of time and those that do, take a very log time, especially GNA analysis for the same low number of permutations. Although they have mentioned in the paper the analysis can be carried out in a low-memory computer, due to the very large size of the matrices (SP matrix around 2.39 GB), the analyses do require rather powerful compters and maynot be ideal for low memory computers.       

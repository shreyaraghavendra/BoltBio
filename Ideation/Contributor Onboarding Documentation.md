#### AIM: Create a pipeline for target identification by integrating the different methods of analysis

Types of Analysis:
1. Differential gene analyis
   1. What is differential gene analysis?
   2. Which kinof data does it require?
   3. Which data formats are used to store the data used in the analysis?
   4. Python packages used in the analysis
   5. Algorithms and tests required for the analysis
   6. Conclusions that can be derived from the analysis
3. Mutation Analysis
5. Pathway Analysis
   - What is PA?
   - Why is it important?
   - Overview of existing methods and workflow
   - Datatypes, formats, platforms
   - Python packages
   - Algorithms employed 
   - Conclusions
   - Useful resources
   - References
7. Literary Analysis 
8. Knowledge Graphs

**OVERVIEW 

High-throughput technologies have made it possible to measure gene expression levels of tens of thousands of genes and to compare the gene expression profiles between two phenotypes(disease vs. control, drug A vs. drug B). Various statistical approaches are used to identify the genes which are differentially expressed (DE) between these phenotypes, such as t test, Z-score, and ANOVA. Although such lists of genes provide valuable information regarding the changes across phenotypes, and play important roles in the downstream analysis, further analysis is required as they don’t explain the complex mechanisms that are involved in the given condition[1].  

Pathway analysis and gene set enrichment analysis(GSEA) are two very important approaches that aim to take gene expression levels and leverage existing knowledge about the given organism in order to identify the underlying biological processes and mechanisms. Pathways are models describing the interactions of genes, proteins, or metabolites within cells, tissues, or organisms, and the interactions are represented as a network with multiple, interrelated nodes and edges and not simple lists of genes. Gene sets, on the other hand, are unordered and unstructured collection of genes associated with a specific biological process (e.g. cell cycle), location (e.g. on chromosome 1), disease (e.g. breast cancer), or even the set of genes that are present in a given pathway (e.g. the set of 128 genes involved in the KEGG cell cycle pathway). While there are several categories by which gene sets can be classified, some of them can be arbitrary.  A very well rounded review detailing the difference between the two can be found here, but the key takeaways can be summarised as follows:
You’ll need to perform pathway analysis when: 
You care about how genes are known to interact
You want to take full advantage of the sizes and directions of measured expression changes
You want to account for the type and direction of interactions on a pathway
You want to predict or explain downstream or pathway-level effects
when you are looking for mechanisms that are specifically affected in your experiment
when you want your results to be based on the most recent knowledge
You’ll need GSEA when:
You are looking for “quick and dirty” answers 
You have arbitrarily defined gene sets.  

Pathway analysis approaches use available pathway databases and the given gene expression data to identify the pathways which are significantly impacted in a given condition[2].It has helped in the  identification of the biological roles of candidate genes to only target the cancer cells and in determining the similarity and dissimilarity between cell lines and tumor samples, and helped in examining the biological function of a set of genes with putative connection between them which was yet to be validated[3].  

The integrated analysis of omics data with pathway analysis has become very popular and  is employed in a number of large-scale projects like The Cancer Genome Atlas (TCGA) and Clinical Proteome Tumor Analysis Consortium (CPTAC). Numerous approaches are available for interpreting single gene lists. For example, the GSEA algorithm can detect upregulated and downregulated pathways in gene expression datasets Web-based methods such as Panther, ToppCluster detect significantly enriched pathways amongst ranked or unranked gene lists and are generally applicable to genes  and proteins from various analyses[4].

Due to the evident importance of this type of analysis, more than 70 pathway analysis methods have been proposed so far. The review by Nguyen, TM. et al., presents the comparison of the performances of 13 representative pathways analysis methods and the key takeaways are discussed briefly here and in detail in the Algorithms and Tests section. 

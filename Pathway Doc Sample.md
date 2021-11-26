PATHWAY ANALYSIS

OVERVIEW 

High-throughput technologies have made it possible to measure gene expression levels of tens of thousands of genes and to compare the gene expression profiles between two phenotypes(disease vs. control, drug A vs. drug B). Various statistical approaches are used to identify the genes which are differentially expressed (DE) between these phenotypes, such as t test, Z-score, and ANOVA. Although such lists of genes provide valuable information regarding the changes across phenotypes, and play important roles in the downstream analysis, further analysis is required as they don’t explain the complex mechanisms that are involved in the given condition[1].

Pathway analysis and gene set enrichment analysis(GSEA) are two very important approaches that aim to take gene expression levels and leverage existing knowledge about the given organism in order to identify the underlying biological processes and mechanisms. Pathways are models describing the interactions of genes, proteins, or metabolites within cells, tissues, or organisms, and the interactions are represented as a network with multiple, interrelated nodes and edges and not simple lists of genes. Gene sets, on the other hand, are unordered and unstructured collection of genes associated with a specific biological process (e.g. cell cycle), location (e.g. on chromosome 1), disease (e.g. breast cancer), or even the set of genes that are present in a given pathway (e.g. the set of 128 genes involved in the KEGG cell cycle pathway). While there are several categories by which gene sets can be classified, some of them can be arbitrary.  A very well rounded review detailing the difference between the two can be found here, but the key takeaways can be summarised as follows:
You’ll need to perform pathway analysis when: 

1.You care about how genes are known to interact
2.You want to take full advantage of the sizes and directions of measured expression changes
3.You want to account for the type and direction of interactions on a pathway
4.You want to predict or explain downstream or pathway-level effects
5.when you are looking for mechanisms that are specifically affected in your experiment
6.when you want your results to be based on the most recent knowledge

You’ll need GSEA when:

You are looking for “quick and dirty” answers 
You have arbitrarily defined gene sets.
  
  
Pathway analysis approaches use available pathway databases and the given gene expression data to identify the pathways which are significantly impacted in a given condition[2].It has helped in the  identification of the biological roles of candidate genes to only target the cancer cells and in determining the similarity and dissimilarity between cell lines and tumor samples, and helped in examining the biological function of a set of genes with putative connection between them which was yet to be validated[3].

The integrated analysis of omics data with pathway analysis has become very popular and  is employed in a number of large-scale projects like The Cancer Genome Atlas (TCGA) and Clinical Proteome Tumor Analysis Consortium (CPTAC). Numerous approaches are available for interpreting single gene lists. For example, the GSEA algorithm can detect upregulated and downregulated pathways in gene expression datasets Web-based methods such as Panther, ToppCluster detect significantly enriched pathways amongst ranked or unranked gene lists and are generally applicable to genes  and proteins from various analyses[4].

Due to the evident importance of this type of analysis, more than 70 pathway analysis methods have been proposed so far. The review by Nguyen, TM. et al., presents the comparison of the performances of 13 representative pathways analysis methods and the key takeaways are discussed briefly here and in detail in the Algorithms and Tests section. 
These pathway analysis methods can be divided into two different categories. The first category includes “non-topology-based” methods (non-TB methods, also known as gene set analysis methods), i.e., methods that do not take advantage of the existing knowledge regarding the positions and roles of the genes within the pathways, the directions and types of the signals transmitted from one gene to another, etc. The first generation in the non-TB category is the over-representation analysis (ORA).Some of the tools that leverage this method are GeneMAPP, GoMiner, DAVID, WebGestalt. This approach takes a list of DE genes as input and identifies the pathways in which the DE genes are over- or underrepresented. The second generation of non-TB approaches includes functional class scoring methods (FCS). The hypothesis behind this approach is that small but coordinated changes in sets of functionally related genes may also be important. This approach eliminates the dependency on the gene selection criteria by taking all gene expressions into consideration.Some of the popular FCS approaches are GSEA, Catmap, GlobalTest, sigPathway. 
 
In considering the pathways as simple un-ordered and unstructured collection of genes—as the non-TB methods do— we discard a substantial amount of knowledge about the biological processes described by these pathways.All the dependencies and interactions between genes that are meant to capture and describe the biological phenomenon are completely ignored. Topology-based methods (TB) have been developed in an attempt to include all this additional knowledge in the analysis. The impact analysis was the first such approach . This was followed by a plethora of over 30 tools and methods that fall in this category including Pathway-Express, SPIA, NetGSA, etc.
The data formats and availability, web-based tools and softwares,  Python packages, and algorithms are discussed in detail in the respective sections. 

DATA TYPES AND FORMATS

Pathway analysis approach typically requires two types of input: experimental data and known biological networks. Most pathway analysis methods analyze data from high-throughput experiments, such as microarrays, next-generation sequencing, or proteomics. They accept either a list of gene IDs or a list of such gene IDs associated with measured changes. These changes could be measured with different technologies and therefore can serve as proxies for different biochemical entities.Biological networks or pathways are often represented in the form of graphs that capture our current knowledge about the interactions of genes, proteins, metabolites, or compounds in an organism.

Input Data Types and File Formats: 

Different analysis methods use different input formats and support (i) a list of differentially expressed (DE) genes, (ii) genes and their fold changes and (iii) an expression matrix. The first two input types can be directly entered on the website or uploaded from users’ local machine as a .txt or .tsv file, in which each row represents a gene. For expression matrix input, a dataset can be represented by two .csv files (command-separated)—one for expression matrix and one for sample grouping. The sample grouping file has two columns in which the first column includes samples and the second column are their corresponding groups (e.g. control or disease). 

Gene ID Conversion: Typically tools are in-built to convert arbitrary gene IDs to Enterz ID

Data Availability: The data can be experimental data or data from NCBI GEO, Cancer Genome Atlas, etc. 

The implementation of analysis methods constrains the software to accept a specific input pathway data format, while the underlying graph models in the methods are independent of the input format. Regardless of the pathway format, this information must be parsed into a computer readable graph data structure before being processed. The implementation may incorporate a parser,or this may be up to the user. For instance, SPIA accepts any signaling pathway or network if it can be transformed into an adjacency matrix representing a directed graph where all nodes are components and all edges are interactions. NetGSA is similarly flexible with regard to signaling and metabolic pathways. SPIA provides KEGG signaling pathways as a set of pre-parsed adjacency matrices. The methods described in this chapter may be restricted to only one pathway database, or may accept several. Analysis methods and input types are listed in the table below from this article.

Output Formats

Pathway ontologies are the definitions of “pathway” used by each pathway database.Unifying ontologies across these databases is accomplished through the use of pathway standard languages. These are standard formats that seek to facilitate the exchange of pathway data between the databases and pathway analysis tools. A number of pathway data are based on the Extensible Markup Language (.xml) or in plain text (.txt) formats. Encoding pathways in such formats makes them readable for both humans and machines. Examples of these standard languages are: the Systems Biology Markup Language (SBML;) the Systems Biology Graphical Notation (SBGN), or the Biological Pathway Exchange (BioPAX)[3]. Gene Matrix Transposed (GMT) format, which is used in gene set enrichment analyses and Simple Interaction Format (SIF) and extended SIF with additional fields, which are useful for network analysis and visualization are also commonly provided for pathway information download[6].

GSEA Data Types and Format: The supporting file formats are listed below:
1. Expression Data Formats
1.1 GCT: Gene Cluster Text file format (*.gct)
1.2 RES: ExpRESsion (with P and A calls) file format (*.res)
1.3 PCL: Stanford cDNA file format (*.pcl)
1.4 TXT: Text file format for expression dataset (*.txt)
2. Phenotype Data Formats
2.1 CLS: Categorical (e.g tumor vs normal) class file format (*.cls)
2.2 CLS: Continuous (e.g time-series or gene profile) file format (*.cls)
3. Gene Set Database Formats
3.1 GMX: Gene MatriX file format (*.gmx)
3.2 GMT: Gene Matrix Transposed file format (*.gmt)
3.3 GRP: Gene set file format (*.grp)
3.4 XML: Molecular signature database file format (msigdb_*.xml)
4. Microarray Chip Annotation Formats
4.1 CHIP: Chip file format (*.chip)
5. Ranked Gene Lists
5.1 RNK: Ranked list file format (*.rnk)
 
DATABASES, WEB-BASED TOOLS, AND SOFTWARES

KEGG:KEGG PATHWAY is a collection of manually drawn pathway maps representing our knowledge of the molecular interaction, reaction and relation networks for:

1. Metabolism

2. Genetic Information Processing

3. Environmental Information Processing

4. Cellular Processes

5. Organismal Systems

6. Human Diseases

7. Drug Development
GO: The Gene Ontology (GO) knowledgebase is the world’s largest source of information on the functions of genes. This knowledge is both human-readable and machine-readable, and is a foundation for computational analysis of large-scale molecular biology and genetics experiments in biomedical research.

REACTOME:Reactome is a free, open-source, curated and peer-reviewed pathway database that provides tools for the visualization, interpretation and analysis of pathway knowledge to support basic research, genome analysis, modeling, systems biology and education.

STRING:STRING is a database of known and predicted protein-protein interactions. The interactions include direct (physical) and indirect (functional) associations; they stem from computational prediction, from knowledge transfer between organisms, and from interactions aggregated from other (primary) databases.

DAVID: Provides a comprehensive set of functional annotation tools for investigators to understand biological meaning behind large list of genes. 

Pathway Commons:

Pathvisio:Free open-source pathway analysis and drawing software which allows drawing, editing, and analyzing biological pathways

Webgestalt(WEB-based Gene SeT AnaLysis Toolkit):  Functional enrichment analysis web tool that supports three well-established and complementary methods for enrichment analysis, including Over-Representation Analysis (ORA), Gene Set Enrichment Analysis (GSEA), and Network Topology-based Analysis (NTA). 

GO-Elite: Application designed to identify a non-redundant set of ontology terms, gene sets and pathways to describe a particular set of genes or metabolites. Can calculate advanced over-representation analysis (ORA) statistics from user gene or metabolite lists, determine the minimal set of biologically distinct ontology terms and pathways from these results and summarize these results at multiple levels

QIAGEN Ingenuity Pathway Analysis (QIAGEN IPA):Quickly visualize and understand complex ‘omics data and perform insightful data analysis and interpretation by placing experimental results within the context of biological systems

Metacore:Delivers high-quality biological systems content in the context of comprehensive, validated biological pathways, providing the essential data and analytics to avoid failed drug development and to accelerate scientific research

GSEA: Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether an a priori defined set of genes shows statistically significant, concordant differences between two biological states (UCSF Web-tool)
Knowledge Engine for Genomics (KnowEnG):  A free-to-use computational system for analysis of genomics data sets, designed to accelerate biomedical discovery. It includes tools for popular bioinformatics tasks such as gene prioritization, sample clustering, gene set analysis, and expression signature analysis. 
Analysis of Pan-omics Data in Human Interactome Network (APODHIN): Platform for integrative analysis of transcriptomics, proteomics, genomics, and metabolomics data for identification of key molecular players and their interconnections exemplified in cancer scenario. APODHIN works on a meta-interactome network consisting of human protein–protein interactions (PPIs), miRNA-target gene regulatory interactions, and transcription factor-target gene regulatory relationships. 

PETAL: A Python tool that automatically explores and detects the most relevant nodes within a KEGG pathway, scanning and performing an in-depth search. PETAL can contribute to discovering novel therapeutic targets or biomarkers that are potentially hidden and not considered in the network under study.

PaintOmics 3: A web-based resource for the integrated visualization of multiple omic data types onto KEGG pathway diagrams. PaintOmics 3 combines server-end capabilities for data analysis with the potential of modern web resources for data visualization, providing researchers with a powerful framework for interactive exploration of their multi-omics information

Consensus pathway analysis (CPA): Aids in conducting pathway analyses using multiple datasets to discover consistent biological mechanisms behind a condition.

Additional Resources: 
Cloud-based tools for pathway analysis:https://academic.oup.com/bib/article/22/1/66/5813255?login=true
Statistical methods and user-friendly software tools to analyze these multi-omics data:http://conesalab.org/

PYTHON PACKAGES 

Python Geneset Network Analysis (PyGNA):A tool for network-aware geneset analysis. PyGNA can either be readily used and easily integrated into existing high-performance data analysis pipelines or as a Python package to implement new tests and analyses. With the increasing availability of population-scale omic data, PyGNA provides a viable approach for large scale geneset network analysis[3]

GOATOOLS: A Python-based library, makes it more efficient to stay current with the latest ontologies and annotations, perform gene ontology enrichment analyses to determine over- and under-represented terms, and organize results for greater clarity and easier interpretation using a novel GOATOOLS GO grouping method. 

PyPathway: An extensible free and open source Python package for functional enrichment analysis, network modeling, and network visualization. The network process module supports various interaction network and pathway databases such as Reactome, WikiPathway, STRING, and BioGRID. The network analysis module implements overrepresentation analysis, gene set enrichment analysis, network-based enrichment, and de novo network modeling. 

SharePathway: Python package for KEGG pathway enrichment analysis with multiple gene lists. Aims at providing users a simple and easy-to-use tool for enrichment analysis on multiple lists of genes simultaneously, which may help gain insight into the underlying biological background of these lists of genes

GSEApy:  Used for RNA-seq, ChIP-seq, Microarray data analysis. It can be used for convenient GO enrichment and to produce publication quality figures in python.

pyMultiOmics: Python package for multi-omics data integration and analysis. It uses the Reactome database to map entities (genes, transcripts, proteins, compounds) to their reactions and pathways. The results is then shown as a network graph. Various analyses such as differential analysis, pathway activity analysis can be performed on this network graph, with the results overlaid on the graph

PyIOmica: Open-source Python package focusing on integrating longitudinal multiple omics datasets, characterizing and categorizing temporal trends. The package includes multiple bioinformatics tools including data normalization, annotation, categorization, visualization and enrichment analysis for gene ontology terms and pathways. Additionally, the package includes an implementation of visibility graphs to visualize time series as networks.

iPANDA:Owned by Insilico Medicine Inc. iPANDA is a pathway analysis method that aids interpretation of large-scale gene expression data. It includes coexpression analysis and gene importance estimation to robustly identify relevant pathways and biomarkers for patient stratification, drug discovery and other applications.


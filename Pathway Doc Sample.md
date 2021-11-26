PATHWAY ANALYSIS

OVERVIEW

High-throughput technologies have made it possible to measure gene expression levels of tens of thousands of genes and to compare the gene expression profiles between two phenotypes(disease vs. control, drug A vs. drug B). Various statistical approaches are used to identify the genes which are differentially expressed (DE) between these phenotypes, such as t test, Z-score, and ANOVA. Although such lists of genes provide valuable information regarding the changes across phenotypes, and play important roles in the downstream analysis, further analysis is required as they don't explain the complex mechanisms that are involved in the given condition[1].

Pathway analysis and gene set enrichment analysis(GSEA) are two very important approaches that aim to take gene expression levels and leverage existing knowledge about the given organism in order to identify the underlying biological processes and mechanisms. Pathways are models describing the interactions of genes, proteins, or metabolites within cells, tissues, or organisms, and the interactions are represented as a network with multiple, interrelated nodes and edges and not simple lists of genes. Gene sets, on the other hand, are unordered and unstructured collection of genes associated with a specific biological process (e.g. cell cycle), location (e.g. on chromosome 1), disease (e.g. breast cancer), or even the set of genes that are present in a given pathway (e.g. the set of 128 genes involved in the KEGG cell cycle pathway). While there are several categories by which gene sets can be classified, some of them can be arbitrary.  A very well rounded review detailing the difference between the two can be found [here](https://advaitabio.com/ipathwayguide/pathway-analysis-vs-gene-set-analysis/), but the key takeaways can be summarised as follows:

You'll need to perform pathway analysis when: 

-   You care about how genes are known to interact

-   You want to take full advantage of the sizes and directions of measured expression changes

-   You want to account for the type and direction of interactions on a pathway

-   You want to predict or explain downstream or pathway-level effects

-   when you are looking for mechanisms that are specifically affected in your experiment

-   when you want your results to be based on the most recent knowledge

You'll need GSEA when:

-   You are looking for "quick and dirty" answers 

-   You have arbitrarily defined gene sets.

Pathway analysis approaches use available pathway databases and the given gene expression data to identify the pathways which are significantly impacted in a given condition[2].It has helped in the  identification of the biological roles of candidate genes to only target the cancer cells and in determining the similarity and dissimilarity between cell lines and tumor samples, and helped in examining the biological function of a set of genes with putative connection between them which was yet to be validated[3].

The integrated analysis of omics data with pathway analysis has become very popular and  is employed in a number of large-scale projects like The Cancer Genome Atlas (TCGA) and Clinical Proteome Tumor Analysis Consortium (CPTAC). Numerous approaches are available for interpreting single gene lists. For example, the GSEA algorithm can detect upregulated and downregulated pathways in gene expression datasets  Web-based methods such as Panther, ToppCluster  detect significantly enriched pathways amongst ranked or unranked gene lists and are generally applicable to genes  and proteins from various analyses[4].

Due to the evident importance of this type of analysis, more than 70 pathway analysis methods have been proposed so far. The [review ](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1790-4)by Nguyen, TM. et al., presents the comparison of the performances of 13 representative pathways analysis methods and the key takeaways are discussed briefly here and in detail in the Algorithms and Tests section. 

These pathway analysis methods can be divided into two different categories. The first category includes "non-topology-based" methods (non-TB methods, also known as gene set analysis methods), i.e., methods that do not take advantage of the existing knowledge regarding the positions and roles of the genes within the pathways, the directions and types of the signals transmitted from one gene to another, etc. The first generation in the non-TB category is the over-representation analysis (ORA).Some of the tools that leverage this method are GeneMAPP, GoMiner, DAVID, WebGestalt. This approach takes a list of DE genes as input and identifies the pathways in which the DE genes are over- or underrepresented. The second generation of non-TB approaches includes functional class scoring methods (FCS). The hypothesis behind this approach is that small but coordinated changes in sets of functionally related genes may also be important. This approach eliminates the dependency on the gene selection criteria by taking all gene expressions into consideration.Some of the popular FCS approaches are GSEA, Catmap, GlobalTest, sigPathway.

In considering the pathways as simple un-ordered and unstructured collection of genes---as the non-TB methods do--- we discard a substantial amount of knowledge about the biological processes described by these pathways.All the dependencies and interactions between genes that are meant to capture and describe the biological phenomenon are completely ignored. Topology-based methods (TB) have been developed in an attempt to include all this additional knowledge in the analysis. The impact analysis was the first such approach . This was followed by a plethora of over 30 tools and methods that fall in this category including Pathway-Express, SPIA, NetGSA, etc.

The data formats and availability, web-based tools and softwares,  Python packages, and algorithms are discussed in detail in the respective sections.

DATA TYPES AND FORMATS

Pathway analysis approach typically requires two types of input: experimental data and known biological networks. Most pathway analysis methods analyze data from high-throughput experiments, such as microarrays, next-generation sequencing, or proteomics. They accept either a list of gene IDs or a list of such gene IDs associated with measured changes. These changes could be measured with different technologies and therefore can serve as proxies for different biochemical entities.Biological networks or pathways are often represented in the form of graphs that capture our current knowledge about the interactions of genes, proteins, metabolites, or compounds in an organism.

Input Data Types and File Formats: Different analysis methods use different input formats and support (i) a list of differentially expressed (DE) genes, (ii) genes and their fold changes and (iii) an expression matrix. The first two input types can be directly entered on the website or uploaded from users' local machine as a .txt or .tsv file, in which each row represents a gene. For expression matrix input, a dataset can be represented by two .csv files (command-separated)---one for expression matrix and one for sample grouping. The sample grouping file has two columns in which the first column includes samples and the second column are their corresponding groups (e.g. control or disease).

Gene ID Conversion: Typically tools are in-built to convert arbitrary gene IDs to Enterz ID

Data Availability: The data can be experimental data or data from NCBI GEO, Cancer Genome Atlas, etc.

The implementation of analysis methods constrains the software to accept a specific input pathway data format, while the underlying graph models in the methods are independent of the input format. Regardless of the pathway format, this information must be parsed into a computer readable graph data structure before being processed. The implementation may incorporate a parser,or this may be up to the user. For instance, SPIA accepts any signaling pathway or network if it can be transformed into an adjacency matrix representing a directed graph where all nodes are components and all edges are interactions. NetGSA is similarly flexible with regard to signaling and metabolic pathways. SPIA provides KEGG signaling pathways as a set of pre-parsed adjacency matrices. The methods described in this chapter may be restricted to only one pathway database, or may accept several. Analysis methods and input types are listed in the table below from [this](https://bioinformatics.cse.unr.edu/publications/4966.pdf) article.

![](https://lh4.googleusercontent.com/ahMKPx8Omd_zKiZJrx-j7hWNZFngscab5qbBaDA-0kZEpMote2-IsJyX8SeDmRWbNcdEBVSY-nTm0jPFa0FZNxBaL7TzQaTDXFiDhvxOO03mDNpndv2xkMI2ohQAGPHppvggGeaC)

Output Formats

Pathway ontologies are the definitions of "pathway" used by each pathway database.Unifying ontologies across these databases is accomplished through the use of pathway standard languages. These are standard formats that seek to facilitate the exchange of pathway data between the databases and pathway analysis tools. A number of pathway data are based on the Extensible Markup Language (.xml) or in plain text (.txt) formats. Encoding pathways in such formats makes them readable for both humans and machines. Examples of these standard languages are: the Systems Biology Markup Language (SBML;) the Systems Biology Graphical Notation (SBGN), or the Biological Pathway Exchange (BioPAX)[3]. Gene Matrix Transposed (GMT) format, which is used in gene set enrichment analyses and Simple Interaction Format (SIF) and extended SIF with additional fields, which are useful for network analysis and visualization are also commonly provided for pathway information download[6].

GSEA Data Types and Format: The supporting file formats are listed below:

[1\. Expression Data Formats](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Expression_Data_Formats)

[1.1 GCT: Gene Cluster Text file format (*.gct)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GCT:_Gene_Cluster_Text_file_format_.28.2A.gct.29)

[1.2 RES: ExpRESsion (with P and A calls) file format (*.res)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RES:_ExpRESsion_.28with_P_and_A_calls.29_file_format_.28.2A.res.29)

[1.3 PCL: Stanford cDNA file format (*.pcl)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#PCL:_Stanford_cDNA_file_format_.28.2A.pcl.29)

[1.4 TXT: Text file format for expression dataset (*.txt)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#TXT:_Text_file_format_for_expression_dataset_.28.2A.txt.29)

[2\. Phenotype Data Formats](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Phenotype_Data_Formats)

[2.1 CLS: Categorical (e.g tumor vs normal) class file format (*.cls)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Categorical_.28e.g_tumor_vs_normal.29_class_file_format_.28.2A.cls.29)

[2.2 CLS: Continuous (e.g time-series or gene profile) file format (*.cls)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Continuous_.28e.g_time-series_or_gene_profile.29_file_format_.28.2A.cls.29)

[3\. Gene Set Database Formats](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats)

[3.1 GMX: Gene MatriX file format (*.gmx)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29)

[3.2 GMT: Gene Matrix Transposed file format (*.gmt)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)

[3.3 GRP: Gene set file format (*.grp)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GRP:_Gene_set_file_format_.28.2A.grp.29)

[3.4 XML: Molecular signature database file format (msigdb_*.xml)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#XML:_Molecular_signature_database_file_format_.28msigdb_.2A.xml.29)

[4\. Microarray Chip Annotation Formats](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Microarray_Chip_Annotation_Formats)

[4.1 CHIP: Chip file format (*.chip)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CHIP:_Chip_file_format_.28.2A.chip.29)

[5\. Ranked Gene Lists](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Ranked_Gene_Lists)

[5.1 RNK: Ranked list file format (*.rnk)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29)

DATABASES, WEB-BASED TOOLS, AND SOFTWARES

[KEGG](https://www.genome.jp/kegg/pathway.html):KEGG PATHWAY is a collection of manually drawn pathway maps representing our knowledge of the molecular interaction, reaction and relation networks for:

1\. Metabolism

2\. Genetic Information Processing

3\. Environmental Information Processing

4\. Cellular Processes

5\. Organismal Systems

6\. Human Diseases

7\. Drug Development

[GO:](http://geneontology.org/docs/ontology-documentation/) The Gene Ontology (GO) knowledgebase is the world's largest source of information on the functions of genes. This knowledge is both human-readable and machine-readable, and is a foundation for computational analysis of large-scale molecular biology and genetics experiments in biomedical research.

[REACTOME](https://reactome.org/):Reactome is a free, open-source, curated and peer-reviewed pathway database that provides tools for the visualization, interpretation and analysis of pathway knowledge to support basic research, genome analysis, modeling, systems biology and education.

[STRING](https://string-db.org/cgi/about):STRING is a database of known and predicted protein-protein interactions. The interactions include direct (physical) and indirect (functional) associations; they stem from computational prediction, from knowledge transfer between organisms, and from interactions aggregated from other (primary) databases.

[DAVID](https://david.ncifcrf.gov/): Provides a comprehensive set of functional annotation tools for investigators to understand biological meaning behind large list of genes. 

Pathway Commons:

[Pathvisio](https://pathvisio.org/):Free open-source pathway analysis and drawing software which allows drawing, editing, and analyzing biological pathways

[Webgestalt](http://www.webgestalt.org/)(WEB-based Gene SeT AnaLysis Toolkit):  Functional enrichment analysis web tool that supports three well-established and complementary methods for enrichment analysis, including Over-Representation Analysis (ORA), Gene Set Enrichment Analysis (GSEA), and Network Topology-based Analysis (NTA). 

[GO-Elite](http://www.genmapp.org/go_elite/help_main.htm): Application designed to identify a non-redundant set of ontology terms, gene sets and pathways to describe a particular set of genes or metabolites. Can calculate advanced over-representation analysis (ORA) statistics from user gene or metabolite lists, determine the minimal set of biologically distinct ontology terms and pathways from these results and summarize these results at multiple levels

[QIAGEN Ingenuity Pathway Analysis](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa/) (QIAGEN IPA):Quickly visualize and understand complex 'omics data and perform insightful data analysis and interpretation by placing experimental results within the context of biological systems

[Metacore](https://clarivate.com/cortellis/solutions/early-research-intelligence-solutions/):Delivers high-quality biological systems content in the context of comprehensive, validated biological pathways, providing the essential data and analytics to avoid failed drug development and to accelerate scientific research

[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp): Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether an a priori defined set of genes shows statistically significant, concordant differences between two biological states (UCSF Web-tool)

[Knowledge Engine for Genomics (KnowEnG)](https://knoweng.org/analyze/):  A free-to-use computational system for analysis of genomics data sets, designed to accelerate biomedical discovery. It includes tools for popular bioinformatics tasks such as gene prioritization, sample clustering, gene set analysis, and expression signature analysis. 

[Analysis of Pan-omics Data in Human Interactome Network (APODHIN)](http://www.hpppi.iicb.res.in/APODHIN/home.html): Platform for integrative analysis of transcriptomics, proteomics, genomics, and metabolomics data for identification of key molecular players and their interconnections exemplified in cancer scenario. APODHIN works on a meta-interactome network consisting of human protein--protein interactions (PPIs), miRNA-target gene regulatory interactions, and transcription factor-target gene regulatory relationships. 

[PETAL](https://github.com/Pex2892/PETAL): A Python tool that automatically explores and detects the most relevant nodes within a KEGG pathway, scanning and performing an in-depth search. PETAL can contribute to discovering novel therapeutic targets or biomarkers that are potentially hidden and not considered in the network under study.

PaintOmics 3: A web-based resource for the integrated visualization of multiple omic data types onto KEGG pathway diagrams. PaintOmics 3 combines server-end capabilities for data analysis with the potential of modern web resources for data visualization, providing researchers with a powerful framework for interactive exploration of their multi-omics information

[Consensus pathway analysis (CPA)](https://bioinformatics.cse.unr.edu/software/cpa/): Aids in conducting pathway analyses using multiple datasets to discover consistent biological mechanisms behind a condition.

Additional Resources: 

1.  Cloud-based tools for pathway analysis:<https://academic.oup.com/bib/article/22/1/66/5813255?login=true>

2.  Statistical methods and user-friendly software tools to analyze these multi-omics data:<http://conesalab.org/>

PYTHON PACKAGES 

[Python Geneset Network Analysis (PyGNA)](http://github.com/stracquadaniolab/pygna):A tool for network-aware geneset analysis. PyGNA can either be readily used and easily integrated into existing high-performance data analysis pipelines or as a Python package to implement new tests and analyses. With the increasing availability of population-scale omic data, PyGNA provides a viable approach for large scale geneset network analysis[3]

[GOATOOLS](https://github.com/tanghaibao/goatools): A Python-based library, makes it more efficient to stay current with the latest ontologies and annotations, perform gene ontology enrichment analyses to determine over- and under-represented terms, and organize results for greater clarity and easier interpretation using a novel GOATOOLS GO grouping method. 

[PyPathway](https://github.com/iseekwonderful/PyPathway): An extensible free and open source Python package for functional enrichment analysis, network modeling, and network visualization. The network process module supports various interaction network and pathway databases such as Reactome, WikiPathway, STRING, and BioGRID. The network analysis module implements overrepresentation analysis, gene set enrichment analysis, network-based enrichment, and de novo network modeling. 

[SharePathway:](https://github.com/GuipengLi/SharePathway) Python package for KEGG pathway enrichment analysis with multiple gene lists. Aims at providing users a simple and easy-to-use tool for enrichment analysis on multiple lists of genes simultaneously, which may help gain insight into the underlying biological background of these lists of genes

[GSEApy](https://pypi.org/project/gseapy/):  Used for RNA-seq, ChIP-seq, Microarray data analysis. It can be used for convenient GO enrichment and to produce publication quality figures in python.

[pyMultiOmics](https://github.com/glasgowcompbio/pyMultiOmics): Python package for multi-omics data integration and analysis. It uses the Reactome database to map entities (genes, transcripts, proteins, compounds) to their reactions and pathways. The results is then shown as a network graph. Various analyses such as differential analysis, pathway activity analysis can be performed on this network graph, with the results overlaid on the graph

[PyIOmica](https://pyiomica.readthedocs.io/en/latest/): Open-source Python package focusing on integrating longitudinal multiple omics datasets, characterizing and categorizing temporal trends. The package includes multiple bioinformatics tools including data normalization, annotation, categorization, visualization and enrichment analysis for gene ontology terms and pathways. Additionally, the package includes an implementation of visibility graphs to visualize time series as networks.

[iPANDA:](https://github.com/varnivey/ipanda)Owned by Insilico Medicine Inc. iPANDA is a pathway analysis method that aids interpretation of large-scale gene expression data. It includes coexpression analysis and gene importance estimation to robustly identify relevant pathways and biomarkers for patient stratification, drug discovery and other applications.
<table class="tg">
<thead>
  <tr>
    <th class="tg-1wig"><span style="font-weight:700;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Python Packages</span></th>
    <th class="tg-1wig"><span style="font-weight:700;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Types of Analysis Performed</span></th>
    <th class="tg-1wig"><span style="font-weight:700;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Input Data Types and Data Formats</span></th>
    <th class="tg-1wig"><span style="font-weight:700;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Output Data Types and Formats</span></th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">PyGNA</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Geneset network topology (GNT) and geneset network analysis (GNA) tests. </span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Genset(.gmt or txt with gene IDs), network(.tsv, graph pickle), gene tables(.csv,)</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Network tables-.csv, .gmt and.tsv</span></td>
  </tr>
  <tr>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">GOATOOLS</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">GOEA:Test the overrepresentation of gene ontology terms in a list of genes or gene products in order to understand their biological significance</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">A copy of the ontology, which describes terms and relationships among them(.obo, JASON), and a set of annotations, which associates the GO terms to specific gene products(GAF, GPAD formats). </span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Directed acyclic graph (DAG)--EXCEL spreadsheet, tab-separated text file, JSON file, or Python variable containing a list of results with the GO results grouped by function </span></td>
  </tr>
  <tr>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">PyPathway </span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Functional set based and network based enrichment analysis algorithms implemented: ORA, GSEA and SPIA</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">ORA: GMT file and for specific ORA with KEGG, REACTOME and GO--specific classes and preloaded datsets</span><br><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">GSEA:genesets, cls files, </span><br><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Gene expression file-pandas dataframe(excel/csv)</span><br><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">GO analysis-python dictionary/TXT file</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Table and graph visualised by respective methods in the </span><a href="https://nadp.me/PyPathway/enrichment/"><span style="font-weight:400;font-style:normal;text-decoration:underline;color:#15C;background-color:transparent">documentation</span></a></td>
  </tr>
  <tr>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">SharePathway</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Enrichment analysis on multiple lists of genes(gene lists of different sample groups) simultaneously for better biological context as it is usually done separately</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Summary file containing path of all the gene list files(.txt, one path per line)</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Result file displayed in .html format</span></td>
  </tr>
  <tr>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">GSEApy</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Convenient GO enrichments and produce publication-quality figures from python</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Input requires a txt file(FPKM, Expected Counts, TPM, et.al), a cls file, and gene_sets file in gmt format</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Pandas dataframe and the plots can be saved as pdf</span></td>
  </tr>
  <tr>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">pyMultiOmics</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">genes, transcripts, proteins, compounds mapped to reactome database</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Transcriptomics, proteomics, metabolomics data(csv, tsv, any table data that can be read by pandas)</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Mapper objects whose queries can be visualised or saved as .csv, .tsv</span></td>
  </tr>
  <tr>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">PyIomica</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Enrichment Analysis:</span><br><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">KEGG Analysis(ORA)</span><br><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">GO Analysis(ORA)</span><br><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Reactome Analysis</span><br><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">(Reactome POST-GET-style analysis.)</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Pandasdataframe/list</span></td>
    <td class="tg-0lax"><span style="font-weight:400;font-style:normal;text-decoration:none;color:#000;background-color:transparent">Generates “.xlsx” report</span></td>
  </tr>
</tbody>
</table>
ALGORITHMS AND TESTS

### Non-TopologyBased (TB) pathway analysis methods

1.  Fisher's exact (FE) test 

2.  Kolmogorov-Smirnov (KS) test 

3.  Wilcoxon rank sum (WRS)

4.  GSEA

5.  GSA 

6.  PADOG

### TB pathway analysis methods

1.  Impact analysis 

2.  CePaGSA and CePaORA 

3.  PathNet 

Details about each of these methods can be found in the paper:<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1790-4#ref-CR67>. 

Additional resources 

1.  <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002375>

2.  <https://en.wikipedia.org/wiki/Pathway_analysis>

3.  <https://academic.oup.com/bib/article/20/2/690/4986938?login=true>

4.  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4681784/>

5.  <https://www.nature.com/articles/s41596-018-0103-9>

CASE STUDY 1

Loss of ASH1L in developing brains causes autistic-like behaviors in a mouse model 

BACKGROUND

Autism spectrum disorder (ASD) is a neurodevelopmental disease associated with various gene mutations. Recent genetic and clinical studies report that mutations of the epigenetic gene ASH1L are highly associated with human ASD and intellectual disability (ID). However, the causal link between ASH1L mutations and ASD/ID remains undetermined. 

HYPOTHESIS/PROBLEM STATEMENT 

ASH1L is an epigenetic factor that regulates gene expression during development. To identify the genes regulated by ASH1L in developing brains, we performed RNA-sequencing (RNA-seq) analyses to examine differential gene expression between wild-type and Ash1L-deleted neural cells.

WORKFLOW

The NPCs (neural progenitor cell)were isolated from the subventricular zone (SVZ) of brains and maintained in serum-free NPC culture medium. The deletion of Ash1L gene in the established NPCs was induced by 4-hydrotamoxifen (4OH-TAM) added in the medium for 10 days(mouse model).

To identify the differentially expressed genes in early NPC differentiation, we performed RNA-seq analyses 0, 12 and 24 hours after induced differentiation.

Results

DGE:The results identified total 2,475 upregulated and 2,808 downregulated genes during induced differentiation (cutoff: fold changes > 1.5, p < 0.01)

Gene ontology (GO) enrichment analyses showed the upregulated genes had enriched GO terms involving nervous system development, while the downregulated genes involved metabolic processes and cell cycle regulation (cutoff: FDR < 0.05) (Figure S3J-K, Table S1, S2). Among all 2,475 genes upregulated in the wild-type NPC differentiation, 70 genes were found to have significantly reduced expression in the Ash1L-KO cells (cutoff: fold changes > 1.5, p < 0.01) These genes downregulated in the Ash1L-KO cells had enriched GO terms involving telencephalon development, regulation of cell communication, brain development, and central nervous development (cutoff: FDR < 0.05) 

Reactome pathway enrichment analysis: Revealed involvement in postsynaptic signal transmission, including NMDA receptor activation (cutoff: FDR < 0.05). Notably, multiple genes, such as Emx2, Dbx2, Pcdh10, and Foxg1, that were previously reported to play critical roles in normal brain development and related neurodevelopmental diseases were found to have significantly reduced expression in the differentiating Ash1L-KO cells (cutoff: fold changes > 1.5, p < 0.01) (Figure 4D), suggesting that loss of ASH1L's function in activating genes critical for normal brain development might be a possible molecular mechanism linking ASH1L mutations to neurodevelopmental defects of ASD/ID.

INFERENCE

In this study, we used an animal model to demonstrate that the loss of Ash1L gene alone in developing mouse brains is sufficient to cause autistic-like behaviors and ID-like deficits in adult mice, which strongly suggests ASH1L gene mutations found in patients are likely to be the causative drivers leading to clinical ASD/ID.

Data Formats and Platform:

The RNA-seq data presented in this study has been deposited to the Gene Expression Omnibus database, GEO accession: GSE173262. Other data have been disclosed in the sections above, or are available from the corresponding author upon reasonable request.

CASE STUDY 2

PyGNA: a unified framework for geneset network analysis

Aim

Analysis of RNA sequencing data generated by The Cancer Genome Atlas (TCGA) project for 6 different types of cancer: 4 epithelial tumors, including 2 from urogenital tissues (BLCA and PRAD), 1 from breast (BRCA) and 1 from lung (LUSC), and 2 from liquid cancers (LAML and DLBC)

To  find whether differentially expressed genes in each cancer show a network effect, and whether they are similar to any other cancer analysed. It is possible to address these questions using the GNT and GNA tests implemented in PyGNA.

Workflow/methodology

TCGA data was retrieved  and differential expression analysis (DEA) was performed using the TCGABiolinks package. 

Data Platforms:  Gene expression data from the Genotype-Tissue EXpression (GTEX) project as control-no control for  LUSC, LAML, and DLBC, and the TCGA tumor data processed by the Recount2 project to avoid biases introduced by different RNA quantification pipelines..Taken together, we retrieved 6 datasets providing mRNA abundance for ≈15000 genes for each tumor and performed differential expression analysis. For each dataset, all genes with FDR<0.01 and logFC|>3 were considered significant.

Analysis

PyGNA was used to perform GNT analysis and GNA analysis between all cancer datasets, using the BioGRID interaction network.

For each test, PyGNA returns the results as a CSV file, which includes descriptive statistics and the parameters of the null distribution used for hypothesis testing.  PyGNA plotting tool (paint-summary-gnt) twas used o visualize a summary of the GNT results for all datasets,  where the test statistic was reported as a z-score, to make them comparable across different tests.GNA analysis between all differentially expressed genesets was performed using the command (paint-comparison-matrix) in PyGNA.

RESULTS

GNT: Only differentially expressed genes in lung and lymphoid cancers show a significant network effect, albeit this is detected by all tests for lung cancer and only by TH and TID lymphoid neoplasm. 

Network effect wasn't observed for the other cancers; this could be explained by the fact these cancers might be controlled not by one highly connected network, but by multiple distinct ones.

PyGNA diagnostic plot generated by the GNT analysis-network effect detected

GNA: While most of them do not seem to show a consistent network effect, GNA can be used to test whether each set of differentially expressed genes are more connected with each other than expected by chance.Using either UH and USP tests, a significant association between breast, bladder, and prostate carcinomas, and between leukemia and lymphoid neoplasms was observed which is consistent with other gene expression analyses, that have shown that anatomically related cancers or with similar histopathology share similar changes in gene expression.

Significant association between lung and lymphoid neoplasms observed in the study might be explained by the fact that lungs contain a vast lymphatic network, which might also be dysregulated in lung tumors.

Conclusion: PyGNA enables network analysis of RNA sequencing datasets and provide useful biological insights about the association between the genes.

CASE STUDY 3

In silico Pathway Activation Network Decomposition Analysis (iPANDA) as a method for biomarker development

BACKGROUND

The inherent complexity of gene network interactions and high diversity of experimental platforms and inconsistency of the data coming from the various types of equipment pose a major challenge in analysing transcriptomics data posing a need for the development of the new large-scale analytical methodologies that infer complex transcriptomic changes more accurately into the network of biologically relevant signalling axes.

Breast cancer data was chosen for the analysis as one of the most challenging in several ways. Since breast cancer has a high degree of intertumour and intratumoural heterogeneity, this cancer type is one of the most difficult in terms of outcome and treatment response prediction. This is especially true for a group of tumours with poor prognosis and fewer number of effective treatments. Thus, traditional methods for transcriptomic data analysis may not be sufficient in this particular case.

The aim is to show that iPANDA is capable of producing highly robust sets of pathway markers, which can be further used for stratification of samples into responder and non-responder groups using neoadjuvant therapy pretreatment breast cancer data with known treatment outcome and receptor status (estrogen receptor and HER2). 

Data Availability and Platform

Six data sets containing gene expression data related to breast cancer patients treated with paclitaxel and three data sets containing transcriptomic data from normal cancer-free breast tissue which were used as a reference from GEO.The MicroArray Quality Control (MAQC) data set (GEO identifier 5350) was obtained from the GEO database.

Platforms: Affymetrix, Agilent

Algorithm

-   The topological weight of each gene is proportional to the number of independent paths through the pathway gene network represented as a directed graph.

-   The computation of topological coefficients for a set of coexpressed genes is inefficient, unless a group of coexpressed genes is being considered as a single unit. To circumvent this challenge, gene modules reflecting the coexpression of genes are introduced in the iPANDA algorithm. 

-   The contribution of gene units (including gene modules and individual genes) to pathway activation is computed as a product of their fold changes in logarithmic scale, topological and statistical weights. Then the contributions are multiplied by a discrete coefficient which equals to +1 or -1 in the case of pathway activation or suppression by the particular unit, respectively. Finally, the activation scores, which we refer to as iPANDA values, are obtained as a linear combination of the scores calculated for gene units that contribute to the pathway activation/suppression. 

### Biomarker identification and relevance: Methodology and Results

-   Measure receiver operating characteristics area under curve (AUC) values to assess the capability of transcriptomic pathway markers to distinguish between two groups of samples: 

-   Gene expression data sets from breast cancer patients with measured response to paclitaxel treatment

-   Pathway activation scores for each sample obtained

-   Lists of the top 30 paclitaxel treatment sensitivity pathway markers obtained for the estrogen receptor negative (ERN) HER2-positive (HER2P) and ERN HER2-negative (HER2N) breast cancer types were obtained.

-   Four and five independent data sets were used for comparison of ERN HER2P and ERN HER2N cancer types

-   Ranked based on AUC values

-   The lists of markers differ significantly between cancer types--complies with the observation that the mechanisms of paclitaxel treatment resistance depend on the breast cancer subtype

-   Produces consistent results for different data sets obtained independently for the same biological case--iPANDA values for 19 and 8 pathways for ERN HER2P and ERN HER2N breast cancer types, respectively, can be utilized as paclitaxel response classifiers with AUC values higher than 0.7 for all data sets examined

-   The common marker pathway (CMP) index to estimate the robustness of the biomarker lists: the combined implementation of the gene modules along with the topology-based coefficients serves as an effective way of noise reduction in gene expression data and allows one to obtain stable pathway activation scores for a set of independent data

-   Pathway activation scores obtained using iPANDA were applied to the identification of paclitaxel neoadjuvant therapy sensitivity in breast cancer: the models developed using normalized iPANDA scores distinguished paclitaxel treatment responders from non-responders with high accuracy

Conclusions

Application of the pathway activation measurement implemented in iPANDA leads to significant noise reduction in the input data to produce highly consistent sets of biologically relevant biomarkers acquired on multiple transcriptomic data sets. More robust since it takes both gene expression and pathway activation into account since they both have their unique insights to offer from the context of a disease.

REFERENCES

1.  <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1790-4>

2.  <https://advaitabio.com/ipathwayguide/pathway-analysis-vs-gene-set-analysis/>

3.  <https://www.frontiersin.org/articles/10.3389/fphys.2015.00383/full>

4.  <https://www.nature.com/articles/s41467-019-13983-9>

5.  <https://academic.oup.com/nar/article/49/W1/W114/6284183>

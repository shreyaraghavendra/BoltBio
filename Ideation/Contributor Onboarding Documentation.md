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
4. Pathway Analysis 
5. Literary Analysis 
6. Knowledge Graphs
7. Biological Network Datasets  

    These contain graphs with gene/protein/disease/drug interactions:  
    http://snap.stanford.edu/biodata/index.html  
    https://www.cs.ucdavis.edu/~filkov/classes/289l-W11/NetworkResources.htm  
    http://snap.stanford.edu/deepnetbio-ismb/  
    http://ncg.kcl.ac.uk/download.php  
    https://thebiogrid.org/  
    http://www.ndexbio.org/#/networkset/e8ebbdde-86dc-11e7-a10d-0ac135e8bacf?accesskey=7fbd23635b798321954e66c63526c46397a3f45b40298cf43f22d07d4feed0fa  

8. Methylation Analysis  

   _What is Methylation, CpG islands? Role in Disease?_  
      https://www.youtube.com/watch?v=bc3wtVXyAXo&ab_channel=NeuralAcademy  
      https://www.slideshare.net/ibadali14/dna-methylation-ppt  

   _What diseases is Methylation linked to?_  
      https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6147084/  
      
   _Epigenetic deregulation:_  
	Most human protein-coding genes harbor multiple transcription starting sites. Promoter regions near these sites regulate transcription and transcript isoforms (multiple combinations of exons within introns translate to different proteins). Tumor-causing methylation in promoter regions can use this mechanism to activate oncogenes/ reduce tumor suppressor genes. Similar to CNA, methylation alteration at gene promoters typically does not alter the coding sequences of 
genes, but contributes to cancer by influencing gene expression. 

	DNAme is also thought to be an important disease-defining feature in many cancers, pointing to the cancer's cell of origin and predictive of the outcome. DNAme changes in cancer have been described along two principal axes: global hypomethylation affecting retroviral elements and genome stability, and focal hypermethylation at promoters of tumor suppressor genes.

	Promoter hypermethylation of TSGs has been surveyed across cancer in TCGA and other studies : DNA repair (MLH1, RBBP8), cell cycle (CDKN2A, CDKN2B), p53 network (CDKN2A, TP73), apoptosis (WIF1, SFRP1), Ras signaling (RASSF1), Wnt signaling (SOX17), and tyrosine kinase cascades (SOCS3).

	Measuring genome-wide DNAm in large numbers of specimens typically uses microarray-based technologies such as the Illumina HumanMethylation450 (450 K) and HumanMethylationEPIC (850 K) [14] arrays, which yield an approximation to the proportion of DNA copies that are methylated at each specific cytosine locus, and are reported as beta values.


   a) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5925605/  
   b) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3620319/
   c) https://github.com/Christensen-Lab-Dartmouth/PyMethylProcess  
   d) https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02220-y  
   e) https://cancerdiscovery.aacrjournals.org/content/11/9/2266  
   	GitHub: https://github.com/HengPan2007/MethSig

	**AIM**: Identify methylation driver genes from passenger genes
MethSig accounts for varying stochastic and disordered hypermethylation rates across genome which results in heterogeneity in the tumor. This is a big challenge in distinguishing driver from passenger genes.
Widely used statistical methods produce thousands of candidate promoter hypermethylation sites, most of which are passengers, very few drivers.

	**DHcR**:  
- Promoter hypermethylation can be measured using differentially hypermethylated cytosine ratio (DHcR) = hypermethylated cytosines (HC) / total number of CpGs profiled in promoters.   
- Uniform background models use the average DHcR which implies assuming a constant HC rate, when in reality the HC rate varies across the genome. This leads to many passenger genes being identified as drivers.   
- Our model uses variable DHcR. This means, MethSig compares a calculated expected methylation rate with promoter specific methylation rate (tumor vs expected).
- This expected value is calculated using a beta regression model and below covariates.  
- How do we know if the BETA model is a good estimate of tumor DHcR value? The drivers were cross-patient aggregated from a large population, not a particular subgroup, as a driver would affect large number of patients.   

	**Covariates included**:  
	(i) Silenced gene expression   
	(ii) Replication time  
	(iii) Proportion of discordant reads (PDR)- If all the CpGs on a specific read are methylated or unmethylated, the read is classified as concordant. Otherwise, it is classified as discordant. At each CpG, the PDR is equal to the number of discordant reads divided by the total number of reads that cover that location. Promoter PDR is given by averaging the values of individual CpGs, as calculated for all CpGs within the promoter of interest that are covered by a minimum of 10 reads that contain at least 4 CpGs. The normal PDR was calculated by averaging PDR of all the normal samples.
Samples with High PDR- low reliability as drivers and vice-versa.

	**Benchmarking against other models:**  
	(i) Q-Q plots with p-values

	(ii) Robustness across datasets- MethSig was applied to individual datasets of 3 tumor types in CLL RRBS. MethSig resulted in a well-calibrated Q-Q plot and a deviation factor closer to 1, compared with benchmarked methods.
	(iii) Accuracy of identified drivers-  MethSig drivers had large overlap with top selections from other studies with similar p value cutoffs. MethSig achieved higher performance compared with benchmarked methods in identifying DNAme drivers associated with gene silencing.

	(iv) MethSig-array designed for methylation arrays generated from Illumina also performed better against other methods with q-q plots and deviation closer to 1.

	**Biological validation:**   
	(i) Diagnosis: Using pathway enrichment analysis (link), established TSGs in many identified drivers, many genes part of core cell pathways associated with cancer, and the whole section under “Candidate CLL DNAme Drivers Include Established TSGs and Were Functionally Validated to Enhance Cancer Cell Fitness”

	(ii) Prognosis: Refer section under “MethSig-nominated CLL DNAme Drivers Provide Independent Prognostic Information and Are Associated with Adverse Outcomes”.

	(iii) Evaluating the clustering of methylated CpG positions- Promoters with a nonrandom distribution of methylated CpGs are more likely to exert repression on corresponding genes and result in a substantial phenotypic impact. MethSig observed higher clustering of methylated positions in hypermethylated drivers. 

   f) https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3443-8  
   g) https://www.sciencedirect.com/science/article/abs/pii/S0888754319309449?dgcid=rss_sd_all  
   h) https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0226461  

   





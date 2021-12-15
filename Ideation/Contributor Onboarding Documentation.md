#### AIM: Create a pipeline for target identification by integrating the different methods of analysis

Types of Analysis:  

I) Differential gene analyis    
   1. What is differential gene analysis?
   2. Which kind of data does it require?
   3. Which data formats are used to store the data used in the analysis?
   4. Python packages used in the analysis
   5. Algorithms and tests required for the analysis
   6. Conclusions that can be derived from the analysis

II) Mutation Analysis

i) What are mutations?
      
A mutation is a change in a DNA sequence. Mutations can result from DNA copying mistakes made during cell division, exposure to ionizing radiation, exposure to chemicals called mutagens, or infection by viruses. Germ line mutations occur in the eggs and sperm and can be passed on to offspring, while somatic mutations occur in body cells and are not passed on.
      
Types of DNA Mutations and Their Impact:
    <table>
    <thead>
        <tr>
            <th>Class of Mutation</th>
            <th>Type of Mutation</th>
            <th>Description</th>
           <th> Human Disease(s) Linked to This Mutation </th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td rowspan=3>Point mutation</td>
           <td> Substitution</td> 
            <td> One base is incorrectly added during replication and replaces the pair in the corresponding position on the complementary strand</td>
            <td>Sickle-cell anemia</td>
        </tr>
        <tr>
            <td>Insertion</td>
           <td> One or more extra nucleotides are inserted into replicating DNA, often resulting in a frameshift </td> 
           <td> One form of beta-thalassemia </td>
        </tr>
        <tr>
            <td>Deletion</td>
            <td>One or more nucleotides is "skipped" during replication or otherwise excised, often resulting in a frameshift</td>
           <td> Cystic fibrosis </td>
        </tr>
        <tr>
            <td rowspan=4>Chromosomal mutation</td>
           <td>Inversion </td>
           <td> One region of a chromosome is flipped and reinserted</td>
           <td> Opitz-Kaveggia syndrome</td>
        </tr>
       <tr>
          <td> Deletion</td>
          <td> A region of a chromosome is lost, resulting in the absence of all the genes in that area</td>
          <td> Cri du chat syndrome</td>
       </tr>
       <tr>
          <td> Duplication</td>
          <td> A region of a chromosome is repeated, resulting in an increase in dosage from the genes in that region </td>
          <td> Some cancers</td>
       </tr>
       <tr>
          <td>Translocation </td>
          <td> A region from one chromosome is aberrantly attached to another chromosome</td>
          <td> One form of leukemia</td>
       </tr>
   <tr>
      <td rowspan=2> Copy number variation</td>
      <td> Gene amplification</td>
      <td> The number of tandem copies of a locus is increased</td>
      <td>Some breast cancers </td>
   </tr>
   <tr>
      <td> Expanding trinucleotide repeat</td>
      <td> The normal number of repeated trinucleotide sequences is expanded</td>
      <td> Fragile X syndrome, Huntington's disease</td>
   </tr>
  </tbody>
</table>
   
ii) What is mutation analysis? 
   
Mutational Analysis is testing for the presence of a specific mutation or set of mutations, as opposed to    complete gene sequencing or mutation scanning, which detect most, if not all, mutations in the tested region. 
   - Sequence homology approach: Helps us to understand what are the significant positions of the protein sequence when it comes to the mutation and whether they are conserved or not.
   - Sequence Structural approach: Analyses the structure and the sequence of the protein
   - Sequence similarity approach: Carries out a sequence similarity matrix in order to understand what the mutations are
   - Conservation based approach  

iii) Packages used in the analysis

a) [DeepVarient](https://github.com/google/deepvariant#how-deepvariant-works) [Python and C++]: An analysis pipeline that uses a deep neural network to call genetic variants from next-generation DNA sequencing data.  

*Pipeline:* 
1. Make examples: creates tf.Example protos for training/calling
2. Call Variants: calling variants with a trained DeepVariant model 
3. Post process variants: Postprocess output from call_variants to produce a VCF file

*Input:* 
1.  A reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format)
    format and its corresponding
    [.fai index file](http://www.htslib.org/doc/faidx.html) generated using the
    `samtools faidx` command.

2.  An aligned reads file in [BAM](http://genome.sph.umich.edu/wiki/BAM) or [CRAM] format
    and its corresponding index file (.bai). The reads must be aligned to the
    reference genome described above.
    
*Different input types and their corresponding models:*  

Aligned reads from:  
*   NGS (Illumina) data for either a
    [whole genome](docs/deepvariant-case-study.md) (WGS model) or [whole exome](docs/deepvariant-exome-case-study.md) (WES model).
*   PacBio HiFi data (PACBIO model) , see the
    [PacBio case study](docs/deepvariant-pacbio-model-case-study.md).
*   Hybrid PacBio HiFi + Illumina WGS (HYBRID model), see the
    [hybrid case study](docs/deepvariant-hybrid-case-study.md).
*   Oxford Nanopore long-read data by using
    [PEPPER-DeepVariant](https://github.com/kishwarshafin/pepper).
*   GenapSys data, by using a
    [model retrained by GenapSys](https://github.com/GenapsysInc/genapsys_deepvariant/blob/master/docs/GenapSys_DeepVariant_WES_Model.md).

*Output:*  
1. A list of all variant calls in VCF format or the [gVCF format](https://github.com/google/deepvariant/blob/r1.2/docs/deepvariant-gvcf-support.md#creating-gvcf-output-with-deepvariant) . 
2. An [VCF stats report](https://github.com/google/deepvariant/blob/r1.2/docs/deepvariant-vcf-stats-report.md) as a .html file

b) [DeepTrio](https://github.com/google/deepvariant/blob/r1.2/docs/deeptrio-details.md): DeepTrio is built on top of DeepVariant. It is intended for variant calling of trios or duos. The main advantage of DeepTrio is that genetic inheritance is considered by a neural network for calling variants in trio samples. Also, variant candidates are generated from all samples at once, which ensures a genotype call is made for any position in the trio with a variant.  

*Pipeline (Same as DeepVariant):*
1. Make examples: creates tf.Example protos for training/calling
2. Call Variants: calling variants with a trained DeepVariant model 
3. Post process variants: Postprocess output from call_variants to produce a VCF file

*Input:*  
1. A reference genome in FASTA format and its corresponding .fai index file generated using the samtools faidx command.
2. An aligned reads files for child and one or two parents in BAM format and its corresponding index file (.bai). The reads must be aligned to the reference genome described above.

*Different input data types and their corresponding models:*

*  Illumina whole genome data (WGS).
*  Illumina whole exome data (WES).
*  PacBio HiFi whole genome data (PacBio WGS).

*Output:*
A set of variants in VCF format representing the child and one or two parents.

c) [BayesPI-BAR2](https://junbaiw.github.io/BayesPI-BAR2/) [Python]: BayesPI-BAR2 is a package designed to predict how non-coding somatic mutations in cancer samples affect protein-DNA binding at the mutated place. Changes in binding of transcription factors to mutated regulatory sequences can lead to disrupted gene regulation, which may promote tumorigenesis. BayesPI-BAR2 takes into account the possibility for several nearby mutations to affect binding of the same protein. The predicted effects are tested for significance in the given patient cohort, and only those that appear in patient samples more frequently than expected by chance are reported.

d) [seqfam](https://seqfam.readthedocs.io/en/latest/tutorial_and_api.html): Seqfam package is primarily designed for analysing next generation sequencing (NGS) DNA data from families with known pedigree information in order to identify rare variants that are potentially causal of a disease/trait of interest.

e) [MutaNET](https://sourceforge.net/projects/mutanet/) [Python]: It includes a next generation sequencing (NGS) pipeline that calls mutations based on paired-end NGS reads, an automated analysis tool and various file converters and mergers. The mutation analysis feature considers the coding region, protein domains, regulation and transcription factor binding site information, and can be used to analyse the potential impact of mutations on genes of interest.  

f) [Orchid](https://academic.oup.com/bioinformatics/article/34/6/936/4587584): Orchid is a python based software package for the management, annotation and machine learning of cancer mutations.

- Annotation python packages

a) [Varcode](https://pypi.org/project/varcode/)  
b) [Viral-ngs](https://viral-ngs.readthedocs.io/en/latest/index.html)  
c) [Gvanno](https://github.com/sigven/gvanno)  
d) [Renovo](https://github.com/mazzalab-ieo/renovo)  
e) [VAPr](https://sci-hub.mksa.top/10.1093/bioinformatics/bty192)


iv) Existing mutation analysis tools and services:

a) [Genewiz](https://www.genewiz.com/Public/Services/Molecular-Genetics/Mutation-Analysis?sc_device=Mobile): GENEWIZ’s Mutation Analysis service helps scientists ramp up mutation detection in coding exons, enabling scientists to quickly analyze and identify mutations that may affect the function of their gene of interest.  
*Input*: Purified genomic DNA or biosafety level 1 (BSL1) and 2 (BSL2) supplied sources from which genomic DNA can be extracted  
*Output*: Raw sequence data files and  report identifying mutations compared to the provided reference sequence  

 b) [ExPecto](https://hb.flatironinstitute.org/expecto): ExPecto is a framework for ab initio sequence-based prediction of mutation gene expression effects and disease risks , [this](https://hb.flatironinstitute.org/expecto) is an explorer of tissue-specific expression effect predictions. 
 
c) [DeepSea](https://hb.flatironinstitute.org/deepsea/): DeepSEA is a deep learning-based algorithmic framework for predicting the chromatin effects of sequence alterations with single nucleotide sensitivity. DeepSEA can accurately predict the epigenetic state of a sequence, including transcription factors binding, DNase I sensitivities and histone marks in multiple cell types, and further utilize this capability to predict the chromatin effects of sequence variants and prioritize regulatory variants.

v) Annotation Tools 

a) [VarAFT](https://varaft.eu/)(Annotation and filtering): VarAFT provides experiments’ quality, annotates, and allows the filtration of VCF files. Data from multiple samples may be combined to address different Mendelian Inherited Disorders, Population Genetics or Cancers.  
b) [Annovar](https://annovar.openbioinformatics.org/en/latest/): ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes (including human genome hg18, hg19, hg38, as well as mouse, worm, fly, yeast and many others). Annovar performs gene-based, region-based and filter-based annotation. [Research paper explaining Annovar](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2938201/)

vi) Databases 

a) [Protein Data Bank](https://www.rcsb.org/#Category-welcome)  
b) [GenBank](https://www.ncbi.nlm.nih.gov/genbank/)  
c) [Online Mendelian Inheritance in Man (OMIM)](https://www.omim.org/#)  
d) [The Human Gene Mutation Database (HGMD)](http://www.hgmd.cf.ac.uk/ac/index.php)  
e) [The International Genome Sample Resource (IGSR)](https://www.internationalgenome.org/)
f) [dbSNP](https://www.ncbi.nlm.nih.gov/snp/)
g) [ENCODE](https://www.genome.gov/Funded-Programs-Projects/ENCODE-Project-ENCyclopedia-Of-DNA-Elements)

vii) Annotation Databases:

a) [Ensembl site](http://www.ensembl.org)  
b) [NCBI](http://www.ncbi.nlm.nih.gov)  
c) [UCSC Genome browser](http://genome.cse.ucsc.edu)  
d) [gnomAD](https://gnomad.broadinstitute.org/)

viii) File formats used in mutational analysis:

a) [BED](https://www.youtube.com/watch?v=uZZ2pHsATCo): The BED (Browser Extensible Data) format is a text file format used to store genomic regions as coordinates and associated annotations. The data are presented in the form of columns separated by spaces or tabs.  
b) [SAM/BAM](https://www.youtube.com/watch?v=tADZk2GsEaE&t=1s): SAM stands for Sequence Alignment/Map format. It is a TAB-delimited text format consisting of a header section, which is optional, and an alignment section.  
c) [Fastq](https://www.youtube.com/watch?v=SZ6suqu-eLA): FASTQ format is a text-based format for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores. Both the sequence letter and quality score are each encoded with a single ASCII character for brevity.  
d) [Fasta](https://www.youtube.com/watch?v=cJm_BGpjnWg): FASTA format is a text-based format for representing either nucleotide sequences or peptide sequences, in which base pairs or amino acids are represented using single-letter codes. A sequence in FASTA format begins with a single-line description, followed by lines of sequence data.  
e) [VCF](https://www.youtube.com/watch?v=Qgb4Ja5VnUQ): VCF stands for Variant Call Format. It is a standardized text file format used for representing SNP, indel, and structural variation calls.
For more information refer [here](http://samtools.github.io/hts-specs/).

Additional Sources:
* [Commonly used biology terms and their definations](https://sg.idtdna.com/pages/education/decoded/article/genotyping-terms-to-know)
* [Computational Approaches in Comparative Genomics book](https://www.ncbi.nlm.nih.gov/books/NBK20260)
* [A survey of tools for variant analysis of next-generation genome sequencing data](https://academic.oup.com/bib/article/15/2/256/210976)

Related Research papers:  
* [KRAS mutation analysis: a comparison between primary tumours and matched liver metastases in 305 colorectal cancer patients](https://www.nature.com/articles/bjc201126)
* [Integrative whole-genome sequence analysis reveals roles of regulatory mutations in BCL6 and BCL2 in follicular lymphoma](https://www.nature.com/articles/s41598-017-07226-4)
* [Detection of Pathogenic Variants With Germline Genetic Testing Using Deep Learning vs Standard Methods in Patients With Prostate Cancer and Melanoma](https://jamanetwork.com/journals/jama/article-abstract/2772962)
* [A systematic genome-wide mapping of oncogenic mutation selection during CRISPR-Cas9 genome editing](https://www.nature.com/articles/s41467-021-26788-6)
* [A framework for detecting noncoding rare variant associations of large-scale whole-genome sequencing studies](https://www.biorxiv.org/content/10.1101/2021.11.05.467531v1.full.pdf)
* [An active learning framework improves tumor variant interpretation](https://www.biorxiv.org/content/10.1101/2021.11.08.467747v1.full.pdf)
* [Population aware deepVariant model](https://www.biorxiv.org/content/10.1101/2021.01.06.425550v2.full)
* [Mutational analysis of known ALS genes in an Italian population-based cohort](https://sci-hub.mksa.top/10.1212/WNL.0000000000011209)
* [Population genetics for target identification](https://sci-hub.mksa.top/10.1016/j.ddtec.2004.08.012)
* [Computational Compensatory Mutation Discovery Approach: Predicting a PARP1 Variant Rescue Mutation](https://www.biorxiv.org/content/10.1101/2021.11.21.469407v1.full.pdf)
* [A Complete Pedigree-Based Graph Workflow for Rare Candidate Variant Analysis](https://www.biorxiv.org/content/10.1101/2021.11.24.469912v1.full.pdf)
* [A general framework for estimating the relative pathogenicity of human genetic variants](https://sci-hub.mksa.top/10.1038/ng.2892)
* [Variant interpretation using population databases: lessons from gnomAD](https://arxiv.org/pdf/2107.11458.pdf)
* [A survey of tools for variant analysis of next-generation genome sequencing data](https://academic.oup.com/bib/article/15/2/256/210976)
* [HugeSeq](https://sci-hub.mksa.top/10.1038/nbt.2134)
* [Assembly of a pan-genome from deep sequencing of 910 humans of African descent](https://www.nature.com/articles/s41588-018-0273-y)


III) Pathway Analysis  
IV) Literary Analysis  
V) Knowledge Graphs  
VI) Biological Network Datasets    

  These contain graphs with gene/protein/disease/drug interactions:  
  http://snap.stanford.edu/biodata/index.html  
  https://www.cs.ucdavis.edu/~filkov/classes/289l-W11/NetworkResources.htm  
  http://snap.stanford.edu/deepnetbio-ismb/  
  http://ncg.kcl.ac.uk/download.php  
  https://thebiogrid.org/  
  http://www.ndexbio.org/#/networkset/e8ebbdde-86dc-11e7-a10d-0ac135e8bacf?accesskey=7fbd23635b798321954e66c63526c46397a3f45b40298cf43f22d07d4feed0fa  

VII) Methylation Analysis  

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

#### AIM: Create a pipeline for target identification by integrating the different methods of analysis

Types of Analysis:  

I) Differential gene analyis</dt>     
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

a) [DeepVarient](https://github.com/google/deepvariant#how-deepvariant-works) [Python]: An analysis pipeline that uses a deep neural network to call genetic variants from next-generation DNA sequencing data.  
*Input:* Aligned reads  
*Input data formats:* .fasta, .bam and .cram files  
*Output:*  
*Output data formats:* .vcf and .gvcf
b) [MutaNET](https://sourceforge.net/projects/mutanet/)[Python]: It includes a next generation sequencing (NGS) pipeline that calls mutations based on paired-end NGS reads, an automated analysis tool and various file converters and mergers. The mutation analysis feature considers the coding region, protein domains, regulation and transcription factor binding site information, and can be used to analyse the potential impact of mutations on genes of interest.  

iv) Existing mutation analysis services:

a) [Genewiz](https://www.genewiz.com/Public/Services/Molecular-Genetics/Mutation-Analysis?sc_device=Mobile): . GENEWIZ’s Mutation Analysis service helps scientists ramp up mutation detection in coding exons, enabling scientists to quickly analyze and identify mutations that may affect the function of their gene of interest.  
*Input*: Purified genomic DNA or biosafety level 1 (BSL1) and 2 (BSL2) supplied sources from which genomic DNA can be extracted  
*Output*: Raw sequence data files and  report identifying mutations compared to the provided reference sequence  
 b) [ExPecto](https://hb.flatironinstitute.org/expecto): ExPecto is a framework for ab initio sequence-based prediction of mutation gene expression effects and disease risks. , [this](https://hb.flatironinstitute.org/expecto) an explorer of tissue-specific expression effect predictions.  
c) [DeepSea](https://hb.flatironinstitute.org/deepsea/): DeepSEA is a deep learning-based algorithmic framework for predicting the chromatin effects of sequence alterations with single nucleotide sensitivity. DeepSEA can accurately predict the epigenetic state of a sequence, including transcription factors binding, DNase I sensitivities and histone marks in multiple cell types, and further utilize this capability to predict the chromatin effects of sequence variants and prioritize regulatory variants.

III) Pathway Analysis 
    
  

IV) Literary Analysis 

V) Knowledge Graphs

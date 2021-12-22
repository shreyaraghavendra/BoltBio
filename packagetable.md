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



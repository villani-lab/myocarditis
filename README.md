[![DOI](https://zenodo.org/badge/478265605.svg)](https://zenodo.org/doi/10.5281/zenodo.11519192)

Welcome! Here, we present our study of immunotherapy-related myocarditis.

## Read the paper :mortar_board: 

Please read and cite our original research article:

- Blum SM, Zlotoff DA, Smith NP, Kernin IJ, Ramesh S, Zubiri L, et al. [Immune Responses in Checkpoint Myocarditis Across Heart, Blood, and Tumor.](https://doi.org/10.1101/2023.09.15.557794) _bioRxiv_. 2023. p. 2023.09.15.557794. doi:10.1101/2023.09.15.557794

<details>
  <summary>BibTeX</summary>
<pre>
@ARTICLE{Blum2023-dg,
  title    = "{Immune Responses in Checkpoint Myocarditis Across Heart, Blood,
              and Tumor}",
  author   = "Blum, Steven M and Zlotoff, Daniel A and Smith, Neal P and Kernin,
              Isabela J and Ramesh, Swetha and Zubiri, Leyre and Caplin, Joshua
              and Tirard, Alice and Sen, Prith and Song, Yuhui and Barth, Jaimie
              L and Slowikowski, Kamil and Nasrallah, Mazen and Tantivit,
              Jessica and Manakongtreecheep, Kasidet and Arnold, Benjamin Y and
              McGuire, John and Pinto, Christopher J and McLoughlin, Daniel and
              Jackson, Monica and Chan, Puiyee and Lawless, Aleigha and Sharova,
              Tatyana and Nieman, Linda T and Gainor, Justin F and Juric, Dejan
              and Mino-Kenudsen, Mari and Sullivan, Ryan J and Boland, Genevieve
              M and Stone, James R and Thomas, Molly F and Neilan, Tomas G and
              Reynolds, Kerry L and Villani, Alexandra-Chloe",
  journal  = "bioRxiv",
  pages    = "2023.09.15.557794",
  abstract = "Immune checkpoint inhibitors (ICIs) are widely used anti-cancer
              therapies that can cause morbid and potentially fatal
              immune-related adverse events (irAEs). ICI-related myocarditis
              (irMyocarditis) is uncommon but has the highest mortality of any
              irAE. The pathogenesis of irMyocarditis and its relationship to
              anti-tumor immunity remain poorly understood. We sought to define
              immune responses in heart, tumor, and blood during irMyocarditis
              and identify biomarkers of clinical severity by leveraging
              single-cell (sc)RNA-seq coupled with T cell receptor (TCR)
              sequencing, microscopy, and proteomics analysis of 28
              irMyocarditis patients and 23 controls. Our analysis of 284,360
              cells from heart and blood specimens identified cytotoxic T cells,
              inflammatory macrophages, conventional dendritic cells (cDCs), and
              fibroblasts enriched in irMyocarditis heart tissue. Additionally,
              potentially targetable, pro-inflammatory transcriptional programs
              were upregulated across multiple cell types. TCR clones enriched
              in heart and paired tumor tissue were largely non-overlapping,
              suggesting distinct T cell responses within these tissues. We also
              identify the presence of cardiac-expanded TCRs in a circulating,
              cycling CD8 T cell population as a novel peripheral biomarker of
              fatality. Collectively, these findings highlight critical biology
              driving irMyocarditis and putative biomarkers for therapeutic
              intervention. \#\#\# Competing Interest Statement S.M.B has been a
              paid consultant to Two River Consulting and Third Rock Ventures.
              He has equity positions in Kronos Bio, 76Bio, and Allogene
              Therapeutics. D.A.Z. has been a paid consultant to Bristol Myers
              Squibb, Freeline Therapeutics, and Intrinsic Imaging. L.Z. has
              received consulting fees from Bristol Myers Squibb and Merck.
              R.J.S has been a paid consultant to Bristol Myers Squibb, Merck,
              Pfizer, Marengo Therapeutics, Novartis, Eisai, Iovance, OncoSec,
              and AstraZeneca and has received research funding from Merck.
              T.G.N has been a paid consultant to Bristol Myers Squibb,
              Genentech, CRC Oncology, Roche, Sanofi and Parexel Imaging
              Pharmaceuticals and has received grant funding from Astra Zeneca
              and Bristol Myers Squibb related to the cardiac effects of immune
              checkpoint inhibitors. K.L.R has served as an advisory board to
              SAGA Diagnostics and received speakers fees from CMEOutfitters and
              Medscape as well as research funding from Bristol Myers Squibb.
              A.C.V. has been a paid consultant to Bristol Myers Squibb.",
  month    =  sep,
  year     =  2023,
  doi      = "10.1101/2023.09.15.557794",
  language = "en"
}
</pre>
</details>

## Explore the data :microscope: 

<table>
<tr>
<td width="33%">
<a href="https://villani.mgh.harvard.edu/myocarditis/app/?ds=tissue_global&gene=IFI27">
<img src="https://github.com/villani-lab/myocarditis/assets/209714/3397fa69-aeea-4d4d-bba2-769492a8a3c0"></img>
</a>
</td>
<td>
<b>Cell Clusters</b>

Metadata variables and gene expression in two-dimensional embeddings.

Heart tissue cells:
- T cells, Myeloid cells, Non-immune cells, Fibroblasts

Blood immune cells:
- Mononuclear phagocytes, CD4 T cells, CD8 T/NK cells, B cells

<table><tr><td><a href="https://villani.mgh.harvard.edu/myocarditis/app/?ds=tissue_global&gene=IFI27">View Cell Clusters</a> :microscope:</a></td></tr></table>
</td>
</tr>
</table>

## Read the source code &#x1F4BB;

This repository includes two main folders:

[functions](functions)
- Functions used in the analysis scripts.

[figure_scripts](figure_scripts)
- The scripts for our analyses.

## Download the data &#x1F4BE;

The raw and processed scRNA-seq gene expression files are available at NCBI GEO [GSE228597].

Sequencing reads will be available at dbGAP accession **PLACEHOLDER**.

[GSE228597]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228597



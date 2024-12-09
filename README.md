Integration of single-nucleus and spatial transcriptomics reveals the
molecular landscape of the human hippocampus
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

Welcome to the `spatial_HPC` project! In this study, we generated
spatially-resolved transcriptomics (SRT) and single-nucleus
RNA-sequencing (snRNA-seq) data from adjacent tissue sections of the
anterior human hippocampus across ten adult neurotypical donors. SRT
data was generated using [10x Genomics
**Visium**](https://www.10xgenomics.com/products/spatial-gene-expression)
(n=36 capture areas) and [10x Genomics **Visium Spatial Proteogenomics**
(SPG)](https://www.10xgenomics.com/products/spatial-proteogenomics) (n=8
capture areas). snRNA-seq data was generated using [10x Genomics
**Chromium**](https://www.10xgenomics.com/products/single-cell-gene-expression)
(n=26 total snRNA-seq libraries).

If you tweet about this website, the data or the R package please use
the <code>\#spatial_HPC</code> hashtag. You can find previous tweets
that way as shown
<a href="https://twitter.com/search?q=%23spatialDLPFC&src=typed_query">here</a>.

Thank you for your interest in our work!

## Study design

<img src="https://research.libd.org/spatial_hpc/img/Copy%20of%20HPC%20figure%201.png" width="1000px" align="left" />

 Experimental design to generate paired single-nucleus RNA-sequencing (snRNA-seq) and spatially-resolved transcriptomics (SRT) data in the human hippocampus. 
(A) Postmortem human tissue blocks containing the anterior hippocampus were dissected from 10 adult neurotypical brain donors. 
(B) Tissue blocks were scored and cryosectioned for snRNA-seq assays (gold), and placement on Visium slides (Visium-H&E, blue). 
(C) Top: Tissue sections (2-4 100μm cryosections per donor) from all ten donors were collected from the same tissue blocks for measurement with the 10x Genomics Chromium 3’ gene expression platform. For each donor, two samples were generated, one sorted based on propidium iodide (PI, purple) and the second sorted based on PI+ and NeuN+ (green). Replicate samples were collected from three donors for a total of n=26 total snRNA-seq libraries. Bottom: 10μm tissue sections from all ten donors were placed onto 2-5 capture areas to include the extent of the HPC (n=36 total capture areas), for measurement with the 10x Genomics Visium-H&E platform. Orientation was verified based on expression of known marker genes. 
(D) Canonical marker genes were identified as spatially variable genes using nnSVG (29).
(E) SRT data was clustered using PRECAST (30) with k=18 and clusters were annotated (columns) based on expression of known marker genes (rows). Cluster groupings indicated at the top of the heatmap define which clusters contributed to the broad domains of Neuron, Neuropil, white matter (WM), and vascular/ cerebrospinal fluid cell-enriched (Vasc/CSF). RHP: retrohippocampus, SUB: subiculum, CA2.4: cornu ammonis (CA) regions 2 through 4 (CA2, CA3, CA4), GCL: dentate gyrus granule cell layer, ML: dentate gyrus molecular layer, SL: stratum lucidum, SR: stratum radiatum, SLM: stratum lacunosum-moleculare, SGZ: dentate gyrus subgranular zone.. (This figure was created with
[Biorender](https://biorender.com))

## Interactive Websites

All of these interactive websites are powered by open source software,
namely:

- 🔍 [`samui`](http://dx.doi.org/10.1017/S2633903X2300017X)
- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)

We provide the following interactive websites, organized by dataset with
software labeled by emojis:

- Visium (n = 44)
  - 👀 <https://libd.shinyapps.io/pseudobulk_HPC/>
    - Provides tools for visualization of pseudobulked Visium data.
  - 🔍 [HPC Samui
    browser](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=Br3942&s=Br8325&s=Br2720&s=Br2743&s=Br3942-VSPG&s=Br6423&s=Br6432&s=Br6471&s=Br6522&s=Br8325-VSPG&s=Br8492&s=Br8667)
    - Provides interactive spot-level visualization of Visium data.
- snRNA-seq (n = 26)
  - 👀 <https://libd.shinyapps.io/HPC_snRNAseq_data/>
    - Provides tools for visualization of snRNA-seq data.

## Data Access

All data, including raw FASTQ files and `SpaceRanger`/`CellRanger`
processed data outputs, can be accessed via Gene Expression Omnibus
(GEO) under accessions
[GSE264692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264692)
(SRT) and
[GSE264624](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264624)
(snRNA-seq).

The spatially-resolved transcriptomics (SRT) and single-nucleus RNA-sequencing (snRNA-seq) data can also be accessed through the bioconductor package at [humanHippocampus2024](https://bioconductor.org/packages/devel/data/experiment/html/humanHippocampus2024.html)
## Contact

We value public questions, as they allow other users to learn from the
answers. If you have any questions, please ask them at
[LieberInstitute/spatial_hpc/issues](https://github.com/LieberInstitute/spatial_hpc/issues)
and refrain from emailing us. Thank you again for your interest in our
work!

## Internal

- JHPCE location:
  `/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc`

### Files:

- `code`: Scripts for running all analyses.
- `plots`: plots generated by analysis scripts.
- `processed-data`
  - `Images`: images used for running `SpaceRanger` and `VistoSeg`.
  - `spaceranger`: `SpaceRanger` output files.
- `raw-data`
  - `sample_info`: metadata about samples.
- `snRNAseq_HPC`: code, plots, and data for snRNA-seq analyses.
- Code for running
  [`GraphST`](https://doi.org/10.1038/s41467-023-36796-3) clustering
  pipeline can be found here:
  <https://github.com/JianingYao/SpatialHPC_graphST_multipleSample>

This GitHub repository is organized along the [*R/Bioconductor-powered
Team Data Science* group
guidelines](https://lcolladotor.github.io/bioc_team_ds/organizing-your-work.html#.Yaf9fPHMIdk).
It follows the
[LieberInstitute/template_project](https://github.com/LieberInstitute/template_project)
structure.

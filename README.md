# NBT
### Network-based transfer (NBT) of pan-cancer immunotherapy responses to guide breast cancer prognosis 
![graphical abstract](./image/Graphical_abstract.jpg)
<p style="text-align: justify;line-height: 1.5;">A new perspective on transferring pan-cancer immunotherapy responses to guide breast cancer prognosis was introduced, providing deeper insights into tumor homogeneity. The protein-protein interaction network was utilized to bridge the genes and clinical phenotypes. A network-based method was proposed to support the identification of gene signatures for breast cancer prognosis based on immunotherapy responses. </p>

<p style="text-align: justify;text-indent: 2em;line-height: 1.5;">This repository provides the code and data for network-based transfer across clinical phenotypes. The figure illustrates the detailed workflow for the identification of gene signatures. The dotted-line module represents the input data originating from various cancer studies.</p>

![workflow](./image/workflow.jpg)

| File                     | Description                                               |
|--------------------------|-----------------------------------------------------------|
| **data**                 | Data for filtering anchor genes and identifying gene signatures |
| **EnrichGeneSignature.R** | Enrichment analysis for gene signatures                   |
| **FilterAnchorGenes.R**   | Filter anchor genes from the pan-cancer immunotherapy cohort |
| **IdentifyGeneSignature.R** | Identify gene signatures from the breast cancer cohort    |
| **main.Rmd**             | R Markdown test file for various modules                  |
| **NetworkPropagation.R** | Network propagation of anchor genes on protein-protein interactions (PPIs) |
| **utils.R**              | Utility functions                                         |


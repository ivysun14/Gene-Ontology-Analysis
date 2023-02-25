# Gene-Ontology-Analysis

This repository contains code in R for a gene ontology analysis workflow. The workflow includes:

- With dataframe manipulation, extract target genes from the input csv
- Conversion of gene symbols to Entrez ID (done with AnnotationDbi package)
- Perform gene ontology analysis using clusterProfiler package
- Visualize ontology terms with package ggplot2

Gene ontology analysis is used to analyze genes associated with breast cancer-related phenotypes obtained from PathFX results on cardiovascular drugs and endocrine drugs. The analysis is part of a small study expansion of the PathFX algorithm, and was completed in the class Computational and Systems Biology 185 at UCLA during Winter quarter 2022.

# Gene-Ontology-Analysis

This repository contains code in R for a gene ontology analysis workflow, which was completed as part of a small study expansion of the [PathFX algorithm](https://github.com/jenwilson521/PathFX). The original PathFX paper can be accessed [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006614). The workflow includes:

- With dataframe manipulation, extract target genes from the input csv
- Conversion of gene symbols to Entrez ID (done with AnnotationDbi package)
- Perform gene ontology analysis using clusterProfiler package
- Visualize ontology terms with package ggplot2

Gene ontology analysis is used to analyze genes associated with breast cancer-related phenotypes obtained from PathFX results on cardiovascular drugs and endocrine drugs. The analysis was completed in the class Computational and Systems Biology 185 at UCLA during Winter quarter 2022. The complete project presentation slides can also be found in the repository.

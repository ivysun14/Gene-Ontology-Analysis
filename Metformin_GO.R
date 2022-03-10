library(ggplot2)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DT)

#Read in PathFX result file.
DB00331 <- read.csv("/Users/junweisun/Documents/C&S BIO 185/Metformin_DB00331.csv")
DB00331 <- as.data.frame(DB00331) %>% column_to_rownames("X")

#Filter for phenotypes involving keyword "breast".
DB00331_assocBreastCancer <- DB00331[grep('[bB]reast', DB00331$phenotype), ]

#Extract genes associated with those phenotypes.
DB00331_assocBreastCancer_geneID <- as.data.frame(DB00331_assocBreastCancer$genes)
colnames(DB00331_assocBreastCancer_geneID)[1] <- "gene"

DB00331_assocBreastCancer_geneID <- unlist(strsplit(DB00331_assocBreastCancer_geneID$gene, ","))
DB00331_assocBreastCancer_geneID <- unique(DB00331_assocBreastCancer_geneID)

#Converting gene names to Entrez ID.
genome <- org.Hs.eg.db
symbols <- DB00331_assocBreastCancer_geneID
geneLst <- select(genome, 
                  keys = symbols,
                  columns = c("ENTREZID", "SYMBOL"),
                  keytype = "SYMBOL")
geneLst <- geneLst$ENTREZID

#GO over-representation analysis using package clusterProfiler.
go <- enrichGO(geneLst, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "BP", readable = TRUE)

#Visualize result.
head(go, 25) %>% dplyr::select(ID, Description, pvalue, p.adjust) %>% datatable(elementId = "goEle")

#Export GO result.
DB00331_assocBreastCancer_GO <- as.data.frame(go@result)
write.csv(DB00331_assocBreastCancer_GO, "DB00331_GO.csv")

#Bar plot.
topGO <- DB00331_assocBreastCancer_GO[1:10, ]
ggplot(data=topGO, aes(x=reorder(Description, -pvalue), y=pvalue)) +
  geom_bar(stat="identity", width=0.5, fill="steelblue") +
  theme_minimal() +
  ggtitle("Pathways enriched for gene variants associated with breast cancer in Metformin") +
  xlab("Biological Processes") +
  theme(plot.title = element_text(hjust = 1.0)) +
  coord_flip()
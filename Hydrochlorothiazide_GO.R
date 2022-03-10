library(ggplot2)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DT)

#Read in PathFX result file.
DB00999 <- read.csv("/Users/junweisun/Documents/C&S BIO 185/Hydrochlorothiazide_DB00999.csv")
DB00999 <- as.data.frame(DB00999) %>% column_to_rownames("X")

#Filter for phenotypes involving keyword "breast".
DB00999_assocBreastCancer <- DB00999[grep('[bB]reast', DB00999$phenotype), ]

#Extract genes associated with those phenotypes.
DB00999_assocBreastCancer_geneID <- as.data.frame(DB00999_assocBreastCancer$genes)
colnames(DB00999_assocBreastCancer_geneID)[1] <- "gene"

DB00999_assocBreastCancer_geneID <- unlist(strsplit(DB00999_assocBreastCancer_geneID$gene, ","))
DB00999_assocBreastCancer_geneID <- unique(DB00999_assocBreastCancer_geneID)

#Converting gene names to Entrez ID.
genome <- org.Hs.eg.db
symbols <- DB00999_assocBreastCancer_geneID
geneLst <- select(genome, 
                  keys = symbols,
                  columns = c("ENTREZID", "SYMBOL"),
                  keytype = "SYMBOL")
geneLst <- geneLst$ENTREZID

#GO over-representation analysis using package clusterProfiler.
go <- enrichGO(geneLst, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "ALL", readable = TRUE)

#Visualize result.
head(go, 25) %>% dplyr::select(ID, Description, pvalue, p.adjust) %>% datatable(elementId = "goEle")

#Export GO result.
DB00999_assocBreastCancer_GO <- as.data.frame(go@result)
write.csv(DB00999_assocBreastCancer_GO, "DB00999_GO.csv")

#Bar plot.
topGO <- DB00999_assocBreastCancer_GO[1:7, ]
ggplot(data=topGO, aes(x=reorder(Description, -pvalue), y=pvalue)) +
  geom_bar(stat="identity", width=0.5, fill="steelblue") +
  theme_minimal() +
  ggtitle("Cellular Functions associated with gene variants in Hydrochlorothiazide") +
  xlab("Cellular Components") +
  theme(plot.title = element_text(hjust = 1.0)) +
  coord_flip()
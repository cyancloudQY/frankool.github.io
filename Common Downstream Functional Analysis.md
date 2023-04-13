# Common Downstream Functional Analysis
## Differentially expressed gene
Use R Deseq2

Just imagine we already had a expression matrix with gene symbol as row names, and sample names as column names. 


```r
######### Group information ##############
library(tidyverse)
library(DESeq2)

counts <- raw_counts

## Could add multiple group info at this step to include covariates
group_info <- factor(c(rep('SLE', 3),rep('C', 3)), levels = c('SLE', 'C'))

######### DESeq2 #############################
DESeq_object <- DESeqDataSetFromMatrix(counts, colData = DataFrame(group_info), design= ~group_info)
DESeq_object <- DESeq(DESeq_object)
DESeq_results <- as.data.frame(results(DESeq_object))


######## Filter ##############################
fold_change <- 1
p_value <- 0.05

DESeq_results_gene$type <- ifelse((DESeq_results_gene$log2FoldChange > fold_change) & (DESeq_results_gene$pvalue < p_value), 'UP',
                             ifelse((DESeq_results_gene$log2FoldChange < -fold_change) & (DESeq_results_gene$pvalue < p_value), 'DOWN','NOT'))
DESeq_diff <- filter(DESeq_results_gene, (DESeq_results_gene$log2FoldChange < -fold_change | DESeq_results_gene$log2FoldChange > fold_change )
                     & DESeq_results_gene$pvalue < p_value)
gene_list <- rownames(DESeq_diff)
gene_list = bitr(gene_list, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")


```

## GO & KEGG
R clusterProfiler
```r
#################################################
############### GO ##############################
#################################################

library(clusterProfiler)
library(org.Hs.eg.db)

ego_all <- enrichGO(gene = gene_list$SYMBOL, 
                    OrgDb = 'org.Hs.eg.db', 
                    ont = "ALL", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE,
                    keyType = 'SYMBOL',
)

go_results <- ego_all@results


#################### Visualization ######################
BP <- go_results[go_results$ONTOLOGY == 'BP',][1:10,]
MF <- go_results[go_results$ONTOLOGY == 'MF',][1:10,]
CC <- go_results[go_results$ONTOLOGY == 'CC',][1:10,]


go_enrich_df<-data.frame(ID=c(BP$ID, CC$ID, MF$ID),
                         Description=c(BP$Description, CC$Description, MF$Description),
                         GeneNumber=c(BP$Count, CC$Count, MF$Count),
                         type=factor(c(rep("biological process", 10), rep("cellular component", 10),rep("molecular function",10)),levels=c("molecular function", "cellular component", "biological process")))


go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))


labels = go_enrich_df$Description
names(labels) = rev(1:nrow(go_enrich_df))
## Draw
CPCOLS <- c("#6495ED", "#8FBC8F", "#F4A460")

ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity") + coord_flip() + 
  scale_fill_manual(values = CPCOLS) +  theme_classic() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray30")) +
  labs(title = "The Most Enriched GO Terms") +
  scale_y_continuous(expand = c(0,0))


########################################################
############### KEGG ###################################
########################################################

kegg_all <- enrichKEGG(gene = gene_list$ENTREZID,
                       organism = "hsa",
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
)

kegg_all <- setReadable(kegg_all, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
View(kegg_all@result)



```

## GSEA
1. R clusterRrofiler. Assume we have DESeq results. 
2. Remeber, the results should not filtered by any method, since we want all the gene lists. 

```r
################## GO ######################
gsea_list <- read.table('./DESeq_results.txt')

gsea_file <- rownames_to_column(gsea_list, 'gene_name')
gsea_file <- dplyr::select(gsea_file, gene_name, log2FoldChange)
gsea_file <- gsea_file[order(gsea_file$log2FoldChange, decreasing = TRUE),]

bitr <- bitr(gsea_file$gene_name,
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = 'org.Hs.eg.db')


gsea_file <- left_join(gsea_file, bitr, by = c('gene_name' = 'SYMBOL'))

gsea_file <- gsea_file[complete.cases(gsea_file), ]

geneList <- gsea_file$log2FoldChange 
names(geneList) <- gsea_file$gene_name

gsea_go <- gseGO(geneList = geneList,
                 OrgDb = org.Hs.eg.db,
                 ont = 'ALL',
                 verbose = FALSE,
                 seed = FALSE,
                 by = 'fgsea',
                 keyType = 'SYMBOL')

### Visualization
gseaplot2(gsea_go, geneSetID = c('GO:0005509', 'GO:0071277', 'GO:0016339'), title = 'ABCDEFG')


################## KEGG #############################
names(geneList) <- gsea_file$ENTREZID
gsea_kegg_DMSO_OSI <- gseKEGG(geneList = geneList,
                     organism = 'hsa',
                     verbose = FALSE,
                     seed = FALSE,
                     by = 'fgsea',
                     pvalueCutoff = 0.5)

gsea_kegg_DMSO_OSI <- setReadable(gsea_kegg_DMSO_OSI, OrgDb = org.Hs.eg.db, keyType="ENTREZID")



################ Special Term #######################
# Need to download the terms from website
# https://www.gsea-msigdb.org/gsea/msigdb/

gmt_H1 <- read.gmt('./gmt/h.all.v2023.1.Hs.symbols.gmt')

names(geneList) <- gsea_file$Symbol
Hallmark <- GSEA(geneList,TERM2GENE = gmt_H1, pvalueCutoff = 0.5)
View(Hallmark@result)

```


## GSVA
1. Need have a expression matrix with gene symbol as row names, and sample names as column names. 
2. Need to adjust kcdf method according to the expression type. (Normalized to Gaussian, Raw counts to Poisson)

```
library(GSVA)
library(pheatmap)


gmt_H1 <- split(gmt_H1$gene, gmt_H1$term, drop = TRUE)
gsva_HM <- gsva(expr = as.matrix(expression_matrix), gmt_H1,
                  kcdf = "Poisson",
                  verbose = F, method = 'gsva', mx.diff=FALSE)

group_info <- factor(c(rep('SLE',3), rep('C', 3)), levels = c('SLE', 'C'))

design_limma <- model.matrix(~group_info)

fit <- lmFit(gsva_HM,design_limma)
fit <- eBayes(fit)

limma_results <- topTable(fit, n = Inf)


pheatmap(gsva_HM, cluster_cols = FALSE, cluster_rows = TRUE)

```






















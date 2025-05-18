library(DESeq2)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(FactoMineR)
library(tibble)
DEseq_custom <- function(RPKM, Count, CTR, TRE){
  cond <- colnames(Count) %>% strsplit('_')
  Condition <- c()
  for (i in cond){
    Condition <- append(Condition, i[1])
  }
  Condition <- factor(Condition)
  ColData <- data.frame(row.names=colnames(Count), Condition)
  dds <- DESeqDataSetFromMatrix(Count, colData = ColData, design = ~ Condition)
  dds$Condition <- relevel(dds$Condition,ref=as.character(Condition[1]))
  dds <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
  res <- results(dds, contrast=c("Condition",TRE,CTR))
  res_de <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res_de <- rownames_to_column(res_de,"Gene_name")
  RPKM <- RPKM %>% filter(rowMeans(.)> 0.4) %>% rownames_to_column("Gene_name")
  mutate(res_de, direction = if_else(abs(log2FoldChange) < 1 | pvalue > 0.05, 'ns', 
                                     if_else(log2FoldChange >= 1, 'up', 'down'))) %>% 
    left_join(RPKM,by = "Gene_name") %>% na.omit() %>% 
    select(Gene_name, contains(c(CTR,TRE)),log2FoldChange, pvalue, direction) %>% 
    group_by(direction) %>% arrange(desc(abs(log2FoldChange)),.by_group = TRUE) %>% summarise(counts = n()) %>% print()
  
  output <- mutate(res_de, direction = if_else(abs(log2FoldChange) < 1 | pvalue > 0.05, 'ns', 
                                               if_else(log2FoldChange >= 1, 'up', 'down'))) %>% 
    left_join(RPKM,by = "Gene_name") %>% na.omit() %>%
    select(Gene_name, contains(c(CTR,TRE)),log2FoldChange, pvalue, direction) %>% 
    group_by(direction) %>% arrange(desc(abs(log2FoldChange)),.by_group = TRUE)
  name <- paste(CTR,'_vs_',TRE,'.DEresults.txt', sep = '')
  write.table(output,name,quote = F, 
              sep = "\t", row.names = F, col.names = T)
  return(output)
}
RPKM <- read.delim("~/workspace/RNA-seq/wenlong/mapping/MOLM13-FTO-mRNA-seq.FTO.RPKM.txt",row.names = 1)
Counts <- read.delim("~/workspace/RNA-seq/wenlong/mapping/MOLM13-FTO-mRNA-seq.FTO.Counts.txt",row.names = 1)
RPKM1 <- filter(RPKM, rowMeans(select(RPKM,1:6))>0.5)

###PCA
#RPKM1 <- column_to_rownames(RPKM, "Gene_name")
gene.pca <- PCA(t(RPKM1), ncp = 2, scale.unit = TRUE, graph = FALSE)
plot(gene.pca)

###DEseq
Drug_DMSO_DE <- DEseq_custom(RPKM, Counts, 'DMSO','FP54')




### GO
DMSO_vs_FP54.DEresults <- read.delim("~/workspace/RNA-seq/wenlong/MOLM13/R/DMSO_vs_FP54.DEresults.txt")


ggo <- enrichGO(gene     = Gene_symbol_to_ENTREZID(down_gene),
                OrgDb    = org.Hs.eg.db,
                ont      = "BP",
                readable = TRUE)

head(ggo,20)
write.table(ggo,"DMSO_vs_FP54_down_GO.txt",quote = F, 
            sep = "\t", row.names = F, col.names = T)


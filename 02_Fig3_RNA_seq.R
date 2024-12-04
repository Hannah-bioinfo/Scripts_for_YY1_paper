setwd("/mnt/hwt2_data1/chuhan/0_co_work/2_Liyang/2_RNAseq_YY1/3_anqi/3_counts/")

library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(dplyr)

cts <- read.table('/lustre/chuhan/wkDIR/')
rownames(cts) <- ko2$Geneid

gene.df <- bitr(str_split(rownames(cts) , "\\.",simplify = T)[,1],
                fromType = "ENSEMBL",
                toType = c("SYMBOL"),
                OrgDb = org.Mm.eg.db)

cts_expr <- inner_join(
  x = cts, 
  y = gene.df, by = "ENSEMBL"
)

metaData <- data.frame(condition=c(rep('ko', 6), rep('wt', 6)), 
                       sample=colnames(cts), 
                       type=c(paste0('assay', 0:5), paste0('assay', 0:5)))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = metaData,
                              design = ~ type + condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","ko","wt"))
summary(res)
head(res[order(res$pvalue),])

library(ggplot2)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p1 <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape = type)) +
  geom_point(size=3) +
  xlim(-12, 12) +
  ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$type)
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


DEG_up <- (res %>% as.data.frame() %>% filter(padj < 0.05 & log2FoldChange > 0.5))
DEG_up$ENSEMBL <- rownames(DEG_up)
DEG_up$ENSEMBL <- unlist(lapply(strsplit(DEG_up$ENSEMBL, "." ,fixed=TRUE),head,1))
gene.df <- bitr(DEG_up$ENSEMBL,
                fromType = "ENSEMBL",
                toType = c("SYMBOL"),
                OrgDb = org.Mm.eg.db)
DEG_up<-merge(DEG_up,gene.df,by.y="ENSEMBL")

DEG_down <- (res %>% as.data.frame() %>% filter(padj < 0.075 & log2FoldChange < -0.5))
DEG_down$ENSEMBL <- rownames(DEG_down)
DEG_down$ENSEMBL <- unlist(lapply(strsplit(DEG_down$ENSEMBL, "." ,fixed=TRUE),head,1))
gene.df <- bitr(DEG_down$ENSEMBL,
                fromType = "ENSEMBL",
                toType = c("SYMBOL"),
                OrgDb = org.Mm.eg.db)
DEG_down<-merge(DEG_down,gene.df,by.y="ENSEMBL")
write.csv(DEG_up[order(-DEG_up$log2FoldChange),], file='./DEG_up_full_5.csv', row.names=FALSE, quote=FALSE, col.names=TRUE)
write.csv(DEG_down[order(DEG_down$log2FoldChange),], file='./DEG_down_full_5.csv', row.names=FALSE, quote=FALSE, col.names=TRUE)

up_david <- read.csv('./DAVID_DEG_up.csv', head=TRUE)
up_david$Category <- str_split(up_david$Term,"\\~",simplify = T)[,1]
up_david$Term <- str_split(up_david$Term,"\\~",simplify = T)[,2]
up_david_data <- data.frame(ID=up_david$Category, Description=up_david$Term, 
                            GeneRatio=up_david[,"X."], BgRatio=up_david$Fold.Enrichment, 
                            pvalue=up_david$PValue, p.adjust=up_david$Bonferroni, 
                            qvalue=up_david$FDR, geneID=up_david$Genes,Count=up_david$Count )
up_david_data$GeneRatio <- paste0(as.character(up_david_data$Count),"/180")
up_david_y=new("enrichResult", result=up_david_data)
pdf(file="up_david.pdf",width = 10,height = 8)
dotplot(up_david_y, showCategory=15,title="EnrichmentGO_up_david")
dev.off()

down_david <- read.csv('./DAVID_DEG_down.csv', head=TRUE)
down_david$Category <- str_split(down_david$Term,"\\~",simplify = T)[,1]
down_david$Term <- str_split(down_david$Term,"\\~",simplify = T)[,2]
down_david_data <- data.frame(ID=down_david$Category, Description=down_david$Term, 
                              GeneRatio=down_david[,"X."], BgRatio=down_david$Fold.Enrichment, 
                              pvalue=down_david$PValue, p.adjust=down_david$Bonferroni, 
                              qvalue=down_david$FDR, geneID=down_david$Genes,Count=down_david$Count )
down_david_data$GeneRatio <- paste0(as.character(down_david_data$Count),"/1106")
down_david_y=new("enrichResult", result=down_david_data)
pdf(file="down_david.pdf",width = 10,height = 8)
dotplot(down_david_y, showCategory=15,title="EnrichmentGO_down_david")
dev.off()

de <- res %>% as.data.frame() %>% filter(!is.na(padj))
de$diffexpressed <- ifelse(de$log2FoldChange > 0.5 & de$padj < 0.05, "Up-regulated",
                           ifelse(de$log2FoldChange < (-0.5) & de$padj < 0.05, "Down-regulated", "n.s."))

pdf(file="1_RNA_DEG_volcano.pdf",width = 10,height = 8)
ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values=c("blue", "#767171", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), linetype='longdash')
dev.off()

TSSg <- FigR::mm10TSSRanges
names(TSSg) <- as.character(TSSg$gene_name)
TSSflank <- GenomicRanges::flank(TSSg,
                                   width = 1000,
                                   both = TRUE)
TSSflank_DEG <- TSSflank[TSSflank$gene_name %in% c(DEG_up_symbol$SYMBOL, DEG_down_symbol$SYMBOL)]
write.table(as.data.frame(TSSflank_DEG), file='./2_DEG_promoter_2k.bed', row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
DEGs_YY1_bound <- read.table('./3_DEG_promoter_overlapped_YY1.bed')[, 7]
DEG_up_symbol %>% as.data.frame() %>% filter(SYMBOL %in% DEGs_YY1_bound)

DEGs_YY1_bound_p001 <- read.table('./3.2_DEG_promoter_overlapped_YY1.bed')[, 7]
write.csv(DEG_up_symbol %>% as.data.frame() %>% filter(SYMBOL %in% DEGs_YY1_bound_p001), 
          file='./4.2_DEG_up_YY1_bound_p001.csv', row.names=FALSE, quote=FALSE, col.names=TRUE)
write.csv(DEG_down_symbol %>% as.data.frame() %>% filter(SYMBOL %in% DEGs_YY1_bound_p001), 
          file='./4.2_DEG_down_YY1_bound_p001.csv', row.names=FALSE, quote=FALSE, col.names=TRUE)

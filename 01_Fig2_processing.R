library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(Signac)
library(monocle)


seurat_integrated <- readRDS(file= '/lustre/chuhan/wkDIR/seurat_integrated.rds')
Idents(seurat_integrated) <- factor(Idents(seurat_integrated), levels= rev(c('MuSC', 'Myocyte', "EC", 'FAP', "Cycling basal cell", "Peri", "SMMC", 'Teno', 
                                                                             "Schwann", 'BC', 'DC', 'TC',"Neutro", 'Mono/Mphage')))
cell_markers <- c( "Pax7", "Myod1", "Myf5", "Myh1", "Acta1", "Pecam1", "Cdh5", "Cd34", "Pdgfra", "Mxd3", 'Rgs5', 'Abcc9', "Myl9",
                   "Scx", "Ptn", "Tnmd", "Mpz","Sox10", "Cd79a", "Cd74", 'Flt3',"Cd3d", "S100a9","Cxcr2", "Ccl6", "C1qa",  "Mrc1")
DotPlot(seurat_integrated,features = cell_markers)+ggtitle("Dot plot of representative cell type-specific markers")+theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 15),
                                                                                                                          axis.text.y = element_text(size = 15),
                                                                                                                          legend.text = element_text(size= 15),legend.title= element_text(size= 15),
                                                                                                                          plot.title = element_text(hjust = 0.5,size = 18))+
  scale_colour_gradientn(colours = c( '#d9d9d9', '#d1c7dc', '#d1b3d8', '#d89dce', '#e184bb', '#df73ae', '#dd60a1', '#da4c92', '#cc3f8f', '#bd328c','#ad258a', '#9d1787'))


seu.wt <- subset(seurat_integrated, subset = (orig.ident == "wt"))
seu.ko  <- subset(seurat_integrated, subset = (orig.ident == "ko"))

DefaultAssay(seu.wt) <- 'RNA'
DefaultAssay(seu.ko) <- 'RNA'

FAPcells <- subset(seurat_integrated, idents='FAP')
FAPcells <- SCTransform(FAPcells) %>% RunPCA(npcs = 30) %>% RunPCA(npcs = 30) %>% FindNeighbors(dims = 1:30)
FAPcells <- RunUMAP(FAPcells, dims = 1:30)
FAPcells <- FindClusters(FAPcells, resolution = c(0.35, 0.3))
markers_all_FAP <- FindAllMarkers(object = FAPcells, 
                                  group.by='integrated_snn_res.0.35',
                                  only.pos = TRUE,
                                  logfc.threshaged = 0.5)

FAPcells.wt <- subset(x = FAPcells, subset= (orig.ident == 'wt'))
FAPcells.ko <- subset(x = FAPcells, subset= (orig.ident == 'ko'))

# the subtypes of WT-FAP
DefaultAssay(FAPcells_reSCT) <- 'SCT'
markers_ko_vs_wt <- FindMarkers(object = FAPcells_reSCT,
                                ident.1="ko",
                                ident.2='wt',
                                group.by="sample",
                                logfc.threshold = 0.25)

markers_all_FAP.wt <- FindAllMarkers(object = FAPcells.wt, 
                                     group.by='SCT_snn_res.0.35',
                                     only.pos = TRUE,
                                     logfc.threshold = 0.5)
markers_all_FAP.ko <- FindAllMarkers(object = FAPcells.ko, 
                                     group.by='SCT_snn_res.0.35',
                                     only.pos = TRUE,
                                     logfc.threshold = 0.5)

counts.data <- as(as.matrix(FAPcells.wt@assays$RNA@data), 'sparseMatrix')
pheno.data <- new('AnnotatedDataFrame', data = FAPcells.wt@meta.data)
feature.data <- data.frame(gene_short_name = row.names(counts.data), row.names = row.names(counts.data))
feature.data <- new('AnnotatedDataFrame', data = feature.data)

cds.wt <- newCellDataSet(counts.data, phenoData = pheno.data, featureData = feature.data, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
cds.wt <- estimateSizeFactors(cds.wt)
cds.wt <- estimateDispersions(cds.wt)
cds.wt <- detectGenes(cds.wt, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds.wt),
                                    num_cells_expressed >= 10))
cds.wt$subtype <- cds.wt$integrated_snn_res.0.35
diff_test_res_1 <- differentialGeneTest(cds.wt[expressed_genes,],
                                        fullModelFormulaStr = "~subtype")

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds.wt.copy <- cds.wt

library(clusterProfiler)
library(org.Mm.eg.db)

# exclude genes related with cell cycle
cell_cycle <- read.table('/lustre/chuhan/reference/cell_cycle_markers/')
gene.df <- bitr(cell_cycle$V2,
                fromType = "ENSEMBL",
                toType = c("SYMBOL"),
                OrgDb = org.Mm.eg.db)

cds.wt <- monocle::setOrderingFilter(cds.wt.copy, setdiff(markers_all_FAP.wt$gene[markers_all_FAP.wt$avg_log2FC > 0.7], gene.df$SYMBOL))

cds.wt <- monocle::reduceDimension(cds.wt, method = 'DDRTree')
cds.wt <- monocle::orderCells(cds.wt)
plot_cell_trajectory(cds.wt, color_by = "subtype")
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$subtype)[,"2"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds.wt <- monocle::orderCells(cds.wt, root_state = GM_state(cds.wt))
plot_cell_trajectory(cds.wt, color_by = "Pseudotime",
                     show_branch_points=FALSE) + theme(text = element_text(size = 15)) + scale_color_gradientn(colours = c("#fde0dd", "#fcc5c0", 
                                                                                                                           "#fa9fb5", "#f768a1",
                                                                                                                           "#dd3497",  "#7a0177", 
                                                                                                                           "#49006a"))
plot_cell_trajectory(cds.wt, color_by = "subtype",show_branch_points=FALSE) + theme(text = element_text(size = 15))

cds.ko <- monocle::setOrderingFilter(cds.ko.copy, markers_all_FAP.ko$gene[markers_all_FAP.ko$avg_log2FC > 0.6]) 
cds.ko <- monocle::reduceDimension(cds.ko, method = 'DDRTree')
cds.ko <- monocle::orderCells(cds.ko)
cds.ko <- monocle::orderCells(cds.ko, root_state = GM_state(cds.ko))
cds.ko <- monocle::setOrderingFilter(cds.ko.copy, markers_all_FAP.ko$gene[markers_all_FAP.ko$cluster == '0'| markers_all_FAP.ko$avg_log2FC
                                                                          > 0.65])
cds.ko <- monocle::reduceDimension(cds.ko, method = 'DDRTree')
cds.ko <- monocle::orderCells(cds.ko)
cds.ko <- monocle::orderCells(cds.ko, root_state = GM_state(cds.ko))
plot_cell_trajectory(cds.ko, color_by = "subtype",show_branch_points=FALSE) + theme(text = element_text(size = 15))

plot_cell_trajectory(cds.ko, color_by = "Pseudotime",
                     show_branch_points=FALSE) + theme(text = element_text(size = 15)) + scale_color_gradientn(colours = c("#fde0dd", "#fcc5c0",
                                                                                                                           "#fa9fb5", "#f768a1",
                                                                                                                           "#dd3497",  "#7a0177",
                                                                                                                           "#49006a"))
FAPcells.wt <- RenameIdents(
  object = FAPcells.wt,
  '0' = 'Fibrogenic',
  '1' = 'Stressed',
  '2' = 'Activated',
  '3' = 'Adipogenic',
  '4' = 'Transitional')

FAPcells.ko <- RenameIdents(
  object = FAPcells.ko,
  '0' = 'Activated',
  '1' = 'Fibrogenic',
  '2' = 'Stressed',
  '3' = 'Adipogenic',
  '4' = 'Inflammantory')


plot_cell_trajectory(cds.ko, color_by = "Pseudotime",
                     show_branch_points=FALSE) + theme(text = element_text(size = 15)) + scale_color_gradientn(colours = c("#fde0dd", "#fcc5c0",
                                                                                                                           "#fa9fb5", "#f768a1",
                                                                                                                           "#dd3497",  "#7a0177",
                                                                                                                           "#49006a")) 
plot_cell_trajectory(cds.ko, color_by = "subtype", show_branch_points=FALSE) + theme(text = element_text(size = 15)) + scale_color_manual(values=c('#00B478',  '#F65144',  '#B044F2',  '#0096D2',  '#F09937'))

immune_genes <- read.table('/lustre/chuhan/wkDIR/immune_genes.txt')[,1]
anti_apop_genes <- read.table('/lustre/chuhan/wkDIR/anti_apop_genes.txt')[,1]
pro_apop_genes <- read.table('/lustre/chuhan/wkDIR/pro_apop_genes.txt')[,1]
x <- AddModuleScore(
  object = FAPcells,
  features = list(immune_genes, anti_apop_genes, pro_apop_genes),
  enrich.name = 'immune_ModScore'
)

dittoRidgePlot(FAPcells, 'immune_ModScore1', group.by='sample', color.panel=rev(c('#BCBBBC', '#7898BF')),ridgeplot.lineweight = 0.5, 
               ridgeplot.scale=1.8, add.line=c(-0.21, -0.05), line.color=c('#5D5B5D', '#2F4763'))+ theme(
                 axis.title.x = element_text(size = 22),
                 axis.text.x = element_text(size = 22, angle=0, hjust= 0.5),
                 axis.title.y = element_text(size = 20,  hjust= 0.5),
                 axis.text.y = element_text(size = 20),
                 legend.title = element_text(size=20), legend.text = element_text(size=20))
dittoRidgePlot(FAPcells, 'pro_apoptosis_ModScore', group.by='sample', color.panel=rev(c('#BCBBBC', '#7898BF')),ridgeplot.lineweight = 0.5,
               ridgeplot.scale=1.8, add.line=c(0.01369707, -0.04346278), line.color=c('#5D5B5D','#2F4763'), sub='n.s.')+ theme(
                 axis.title.x = element_text(size = 22),
                 axis.text.x = element_text(size = 22, angle=0, hjust= 0.5),
                 axis.title.y = element_text(size = 20,  hjust= 0.5),
                 axis.text.y = element_text(size = 20),
                 legend.title = element_text(size=20), legend.text = element_text(size=20))



genes_ko <- c('Cxcl5', 'Cxcl3', 'Ccl7', 'Ccl2', 
              'Cxcl14', 'Lum', 'Smoc2', 'Podn',  #Lpl
              'Hspb1', 'Hspd1', 'Hspe1', 'Fos', 'Atf3','Jun',
              'Pi16', 'Dpp4', 'Igfbp5', 
              'Cd74','Apoe','C1qa','C1qb','H2-Eb1','H2-Ab1')

genes_ct <- c('Cxcl5', 'Cxcl3', 'Ccl7', 'Ccl2', 
              'Cxcl14', 'Lum', 'Smoc2', 'Podn', 
              'Hspb1', 'Hspd1', 'Hspe1', 'Fos', 'Atf3','Jun',
              'Pi16', 'Dpp4', 'Igfbp5',   
              'Apod', 'Ptx3','Muoc', 'Mt1')

library(dittoSeq)
library(export)
dittoHeatmap(FAPcells.ko, intersect(genes_ko, rownames(FAPcells.ko@assays$SCT)), annot.by = c("subtype"),
             annot.colors = c('#00B478', '#F65144',  '#B044F2',  '#0096D2',  '#F09937'), 
             assay = "SCT",
             show_colnames = FALSE,
             show_rownames = TRUE,
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             complex=FALSE, heatmap.colors =rev(c('#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf', '#e0f3f8',
                                                  '#abd9e9', '#74add1', '#4575b4', '#313695')), fontsize_row=18)
color_bar <- rev(c('#DA0034', '#E05952', '#F6815C', '#FDC263', "#FED972", '#ffffbf', '#e0f3f8',
                   '#abd9e9', '#74add1', '#4575b4', '#313695'))
c('#F8766D',   '#B044F2',  '#00B478',  '#0096D2',  '#A0A200')
dittoHeatmap(FAPcells.wt, intersect(genes_ct, rownames(FAPcells.wt@assays$SCT)), annot.by = c("subtype"),
             annot.colors = c('#00B478',   '#F8766D',  '#B044F2',  '#0096D2',  '#A0A200'), 
             assay = "SCT",
             show_colnames = FALSE,
             show_rownames = TRUE,
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             complex=FALSE, heatmap.colors =color_bar, fontsize_row=18)
dev.off()

dittoHeatmap(FAPcells.ko, genes,
             annot.by = c("subtype"),
             assay = "SCT",
             show_colnames = FALSE,
             show_rownames = TRUE,
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             complex=FALSE, heatmap.colors = c("#006837","#238443", "#41ab5d",
                                               "#78c679","#addd8e", "#d9f0a3",
                                               "#fde0dd", "#fcc5c0", "#fa9fb5", "#f768a1",
                                               "#dd3497", "#ae017e", "#7a0177", "#49006a"), fontsize_row=18)
c("#004529","#006837","#238443", "#41ab5d",
  "#78c679","#addd8e", "#d9f0a3", "#f7fcb9",
  "#fde0dd", "#fcc5c0", "#fa9fb5", "#f768a1",
  "#dd3497", "#ae017e", "#7a0177", "#49006a")
c("#004529","#006837","#238443", "#41ab5d",
  "#78c679","#addd8e",  "#fcc5c0", "#fa9fb5", "#f768a1",
  "#dd3497", "#ae017e", "#7a0177", "#49006a")
c('#f7fcb9', '#fff1bd', '#ffe8c7', '#ffe3d3', '#fde0dd')

gene_ct_MP <- c('Cxcl3','Vcan', 'Chil3',
                'Spp1','Fabp5', 'Cd36',
                'C1qa','C1qb', 'C1qc',
                'H2-Aa','H2-Eb1', 'H2-Ab1',
                'Ccna2', 'Ccnb2','Cdk1','Cdc20','Cdca3','Cdca8',
)
gene_ko_MP <- c('Cxcl3','Vcan', 'Chil3',
                'Spp1','Fabp5', 'Cd36',
                'C1qa','C1qb', 'C1qc',
                'H2-Aa','H2-Eb1', 'H2-Ab1',
                'Ccl7', 'Ccl4', 'Ccl2', 'Ccl8', 'Ccl6', 'Ccl24', 'Ccl3')


dittoHeatmap(MPcells.ko, intersect(gene_ko_MP, rownames(MPcells.ko@assays$SCT)), annot.by = c("subtype"),
             annot.colors = c('#00B478', '#B044F2',  '#0096D2', '#F8766D',  '#F09937'), 
             assay = "SCT",
             show_colnames = FALSE,
             show_rownames = TRUE,
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             complex=FALSE, heatmap.colors =color_bar, fontsize_row=18)
dittoHeatmap(MPcells.wt, intersect(gene_ct_MP, rownames(MPcells.wt@assays$SCT)), annot.by = c("subtype"),
             annot.colors = c('#00B478', '#B044F2',  '#0096D2', '#F8766D',  '#A0A200'), 
             assay = "SCT",
             show_colnames = FALSE,
             show_rownames = TRUE,
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             complex=FALSE, heatmap.colors =color_bar, fontsize_row=18)

library(methods)
library(optparse)
library(ggplot2)
library(qqman)
library(dplyr)
library(data.table)
library(stringr)
#bioclite install
#library(ensembldb)
#library(EnsDb.Hsapiens.v79)
#library(EnsDb.Hsapiens.v86)

library("DESeq2")
library("RColorBrewer")
library("gplots")#
library("genefilter")
library("optparse")
library("ggplot2")
library("pheatmap")
library(Biobase)
library(goseq)
library(GO.db)
library(tidyr)
library(ggrepel)
library(windowscanr)

plotPCA2<-function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:4]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
                                                                                                            100), "% variance")) + coord_fixed()
}


# import genes of interest ------------------------------------------------

#This whole script looks at just a subsample of genes
#first take those 29 vrom venn diagram
df.de.venn <- read.csv("output/venn_diagram/overlap_29_genes_all_three_cats_fdr_0.05.csv")
df.de.t <- read.csv("output/venn_diagram/t_implant_de_genes_only_fdr_0.05.csv")
de <- df.de.t$gene_name
df.sex.null <- read.csv("output/venn_diagram/M_more_F_more_shoulder_same_pheno.csv")
df.t_ssp_overlap_exclude_sex <- read.csv("output/venn_diagram/overlap_genes_t_treat_subspecies_More_Lore_exclude_sex.csv")
# first is just SP --------------------------------------------------------


sample_id <- dir(file.path("data/star_DoCounts"))
kal_dirs <- file.path("data/star_DoCounts", sample_id)
directory <- "data/star_DoCounts/"
i.s2c <- read.table("output/sample_table_RSEM.csv", header = TRUE, stringsAsFactors = TRUE, sep = ",")
s2c <- i.s2c

s2c <- dplyr::filter(s2c, samplePart == "SP") %>% 
  #filter(samplePopulation == "lorentzi" | samplePopulation == "moretoni") %>% 
  filter(!grepl("M", sampleName) & sampleName!="47745-moretoni-Chest" & Paired!="N")
#remove apparent outlier - noted as very long feathers
#s2c<-dplyr::filter(s2c,sampleName!="97513-SP-before") #outlier inspections indiv genes in Timplant comparison. Looks ok otherwise though

s2c <- s2c %>% mutate(path = paste("data/star_DoCounts/",fileName, sep = ""))

# Restrict to females and remove T treatment individuals
s2c <- s2c %>% dplyr::rename(sample = sampleName)
#s2c <- s2c %>% filter(Sex != "M" & Treatment != "After")



all.hts<-DESeqDataSetFromHTSeqCount(sampleTable=s2c, directory=directory, design=~PopPartColor)

all.dds<-DESeq(all.hts)
resultsNames(all.dds)
all.res<-results(all.dds)
resultsNames(all.dds)

#QC using PCA
all.vsd <- vst(all.dds, blind=TRUE)
pcaData<-plotPCA2(all.vsd, intgroup = c("PopPart"),returnData=TRUE)
pcaData <- pcaData %>% dplyr::select(-PopPart)
pca2<-cbind(pcaData,s2c)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pca2, aes(PC1, PC2,color = sampleColor, shape = samplePopulation)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggplot(pca2, aes(PC2, PC3,color = sampleColor, shape = samplePopulation)) +
  geom_point(size=6) +
  xlab(paste0("PC2: ",percentVar[3],"% variance")) +
  ylab(paste0("PC3: ",percentVar[4],"% variance")) + 
  coord_fixed()

# get DE genes ------------------------------------------------------------
#s2c.t <- s2c %>% filter(Treatment!="")
#s2c.t$ind <- as.factor(s2c.t$ind)
#t.hts <- DESeqDataSetFromHTSeqCount(sampleTable=s2c.t, directory=directory, design=~ind + Treatment)
#
#t.dds<-DESeq(t.hts)
#resultsNames(t.dds)
#t.res<-results(t.dds, name = "Treatment_Before_vs_After")
#de <- rownames(t.res[t.res$padj<0.1 & !is.na(t.res$padj) & t.res$log2FoldChange < 0, ])
#de <- rownames(t.res[t.res$padj<0.1 & !is.na(t.res$padj), ])

# -------------------------------------------------------------------------
#https://www.biostars.org/p/371379/

sampleDists <- dist(t(assay(all.vsd[de,])))


library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(all.vsd$samplePopulation, all.vsd$ind, sep="-")
#rownames(sampleDistMatrix) <- all.vsd$ind
colnames(sampleDistMatrix) <- NULL
#colors <- colorRampPalette( brewer.pal(9, "Greys") )(255)
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)

dir.create("output/heatmaps_2020", showWarnings = F)

pdf("output/heatmaps_2020/shoulder_only_all_t_de_genes.pdf", width = 12, height = 10)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         cutree_cols = 3,
         treeheight_row = 0)
dev.off()


# all chest and shoulder --------------------------------------------------

all.parts <- i.s2c %>% filter(samplePart != "Crown") %>% 
  filter(!grepl("M", sampleName) & sampleName!="47745-moretoni-Chest" & Paired!="N") %>% 
  mutate(path = paste("data/star_DoCounts/",fileName, sep = "")) %>% 
  dplyr::rename(sample = sampleName)

all.parts.hts<-DESeqDataSetFromHTSeqCount(sampleTable=all.parts, directory=directory, design=~PopPartColor)

all.parts.dds<-DESeq(all.parts.hts)
all.parts.res<-results(all.parts.dds)

#QC using PCA
all.parts.vsd <- vst(all.parts.dds, blind=TRUE)

all.parts.sampleDists <- dist(t(assay(all.parts.vsd[de,])))

all.parts.sampleDistMatrix <- as.matrix(all.parts.sampleDists)
rownames(all.parts.sampleDistMatrix) <- paste(all.parts.vsd$PopPartColor, all.parts.vsd$ind, sep="-")
#rownames(sampleDistMatrix) <- all.vsd$ind
colnames(all.parts.sampleDistMatrix) <- NULL
#colors <- colorRampPalette( brewer.pal(9, "Greys") )(255)
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)



pdf("output/heatmaps_2020/allparts_only_all_t_de_genes.pdf", width = 12, height = 10)
pheatmap(all.parts.sampleDistMatrix,
         clustering_distance_rows=all.parts.sampleDists,
         clustering_distance_cols=all.parts.sampleDists,
         col=colors,
         cutree_cols = 3,
         treeheight_row = 0)
dev.off()

pcaData<-plotPCA2(all.parts.vsd, intgroup = c("PopPart"),returnData=TRUE)
pcaData <- pcaData %>% dplyr::select(-PopPart)
pca2<-cbind(pcaData,all.parts)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pca2, aes(PC1, PC2,color = sampleColor, shape = samplePopulation)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggplot(pca2, aes(PC2, PC3,color = sampleColor, shape = samplePopulation)) +
  geom_point(size=6) +
  xlab(paste0("PC2: ",percentVar[3],"% variance")) +
  ylab(paste0("PC3: ",percentVar[4],"% variance")) + 
  coord_fixed()



# what about males? shoulders-------------------------------------------------------

sex <- i.s2c %>% filter(samplePart == "SP") %>% 
  filter(sampleName!="47745-moretoni-Chest" & Paired!="N") %>% 
  mutate(path = paste("data/star_DoCounts/",fileName, sep = "")) %>% 
  dplyr::rename(sample = sampleName)

sex.hts<-DESeqDataSetFromHTSeqCount(sampleTable=sex, directory=directory, design=~PopPartColor)

sex.dds<-DESeq(sex.hts)
sex.res<-results(sex.dds)

#QC using PCA
sex.vsd <- vst(sex.dds, blind=TRUE)

sex.sampleDists <- dist(t(assay(sex.vsd[as.character(df.de.venn$gene_name),])))
#sex.sampleDists <- dist(t(assay(sex.vsd)))

sex.sampleDistMatrix <- as.matrix(sex.sampleDists)
rownames(sex.sampleDistMatrix) <- paste(sex.vsd$Sex, sex.vsd$PopPartColor, sex.vsd$ind, sep="-")
#rownames(sampleDistMatrix) <- all.vsd$ind
colnames(sex.sampleDistMatrix) <- NULL
#colors <- colorRampPalette( brewer.pal(9, "Greys") )(255)
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)

pdf("output/heatmaps_2020/sex_t_venn_genes.pdf", width = 12, height = 10)
pheatmap(sex.sampleDistMatrix,
         clustering_distance_rows=sex.sampleDists,
         clustering_distance_cols=sex.sampleDists,
         col=colors,
         cutree_cols = 2,
         treeheight_row = 0)
dev.off()

#check T implant genes (not those that overlap in venn diagram)
sex.sampleDists2 <- dist(t(assay(sex.vsd[as.character(de),])))
#sex.sampleDists <- dist(t(assay(sex.vsd)))

sex.sampleDistMatrix2 <- as.matrix(sex.sampleDists2)
rownames(sex.sampleDistMatrix2) <- paste(sex.vsd$Sex, sex.vsd$PopPartColor, sex.vsd$ind, sep="-")
#rownames(sampleDistMatrix) <- all.vsd$ind
colnames(sex.sampleDistMatrix2) <- NULL
pheatmap(sex.sampleDistMatrix2,
         clustering_distance_rows=sex.sampleDists2,
         clustering_distance_cols=sex.sampleDists2,
         col=colors,
         cutree_cols = 2,
         treeheight_row = 0)

#check only T implant and ssp genes, that arnt overlapping with sex (bottom right of venn)
sex.sampleDists3 <- dist(t(assay(sex.vsd[as.character(df.t_ssp_overlap_exclude_sex$gene_name),])))
#sex.sampleDists <- dist(t(assay(sex.vsd)))

sex.sampleDistMatrix3 <- as.matrix(sex.sampleDists3)
rownames(sex.sampleDistMatrix3) <- paste(sex.vsd$Sex, sex.vsd$PopPartColor, sex.vsd$ind, sep="-")
#rownames(sampleDistMatrix) <- all.vsd$ind
colnames(sex.sampleDistMatrix3) <- NULL
pheatmap(sex.sampleDistMatrix3,
         clustering_distance_rows=sex.sampleDists3,
         clustering_distance_cols=sex.sampleDists3,
         col=colors,
         cutree_cols = 2,
         treeheight_row = 0)

#only M v F in more genes
sex.sampleDists4 <- dist(t(assay(sex.vsd[as.character(df.sex.null$gene_name),])))
#sex.sampleDists <- dist(t(assay(sex.vsd)))

sex.sampleDistMatrix4 <- as.matrix(sex.sampleDists4)
rownames(sex.sampleDistMatrix4) <- paste(sex.vsd$Sex, sex.vsd$PopPartColor, sex.vsd$ind, sep="-")
#rownames(sampleDistMatrix) <- all.vsd$ind
colnames(sex.sampleDistMatrix4) <- NULL
pdf("output/heatmaps_2020/sex_M_more_F_more_shoulder_genes.pdf", width = 12, height = 10)
pheatmap(sex.sampleDistMatrix4,
         clustering_distance_rows=sex.sampleDists4,
         clustering_distance_cols=sex.sampleDists4,
         col=colors,
         cutree_cols = 2,
         treeheight_row = 0)
dev.off()
# all chest and shoulder include Males --------------------------------------------------

all.parts <- i.s2c %>% filter(samplePart != "Crown") %>% 
  filter(sampleName!="47745-moretoni-Chest" & Paired!="N") %>% 
  mutate(path = paste("data/star_DoCounts/",fileName, sep = "")) %>% 
  dplyr::rename(sample = sampleName)
all.parts$ind <- as.factor(all.parts$ind)
all.parts.hts<-DESeqDataSetFromHTSeqCount(sampleTable=all.parts, directory=directory, design=~PopPartColor)

all.parts.dds<-DESeq(all.parts.hts)
all.parts.res<-results(all.parts.dds)
plotDispEsts(all.parts.dds)

#QC using PCA
all.parts.vsd <- vst(all.parts.dds, blind=TRUE)

#all.parts.sampleDists <- dist(t(assay(all.parts.vsd[de,])))
all.parts.sampleDists <- dist(t(assay(all.parts.vsd[as.character(df.de.venn$gene_name),]))) #using venn diagram overlaps

all.parts.sampleDistMatrix <- as.matrix(all.parts.sampleDists)
rownames(all.parts.sampleDistMatrix) <- paste(all.parts.vsd$Sex, all.parts.vsd$PopPartColor, all.parts.vsd$ind, sep="-")
#colnames(all.parts.sampleDistMatrix) <- paste(all.parts.vsd$Sex, all.parts.vsd$PopPartColor, all.parts.vsd$ind, sep="-")
colnames(all.parts.sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)
pheatmap(all.parts.sampleDistMatrix,
         clustering_distance_rows=all.parts.sampleDists,
         clustering_distance_cols=all.parts.sampleDists,
         col=colors,
         cutree_cols = 3,
         treeheight_row = 0)



#rownames(sampleDistMatrix) <- all.vsd$ind
#colnames(all.parts.sampleDistMatrix) <- NULL
#colors <- colorRampPalette( brewer.pal(9, "Greys") )(255)
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)

rld <- rlog(all.parts.dds, blind=T)
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
rld_cor <- cor(rld_mat)    ## cor() is a base R function

pdf("output/heatmaps_2020/allparts_only_all_t_de_genes_with_males_rld.pdf", width = 12, height = 10)
pheatmap(rld_cor) #this would be all genes, not just the T implant ones
dev.off()


sex.sampleDists <- dist(t(assay(all.parts.vsd[de,])))

sex.sampleDistMatrix <- as.matrix(sex.sampleDists)
rownames(sex.sampleDistMatrix) <- paste(all.parts.vsd$Sex, all.parts.vsd$PopPartColor, all.parts.vsd$ind, sep="-")
#rownames(sampleDistMatrix) <- all.vsd$ind
colnames(sex.sampleDistMatrix) <- NULL
#colors <- colorRampPalette( brewer.pal(9, "Greys") )(255)
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)

pdf("output/heatmaps_2020/sex_chest_sp_all_t_de_genes.pdf", width = 12, height = 10)
pheatmap(all.parts.sampleDistMatrix,
         clustering_distance_rows=all.parts.sampleDistMatrix,
         clustering_distance_cols=all.parts.sampleDistMatrix,
         col=colors,
         cutree_cols = 3,
         treeheight_row = 0)
dev.off()







##test only sp 
#keep.these <- grep("SP",colnames(all.parts.sampleDistMatrix), value = T)
#m1 <- all.parts.sampleDistMatrix[, colnames(all.parts.sampleDistMatrix) %in% keep.these]
#m2 <- m1[rownames(m1) %in% keep.these, ]
#m2
#
#colnames(m2) <- NULL
#
#keep.these.2 <- grep("SP", colnames(all.parts.sampleDists), value = T)
#
#pdf("output/heatmaps_2020/allparts_only_all_t_de_genes_with_males.pdf", width = 12, height = 10)
#pheatmap(m2,
#         clustering_distance_rows=all.parts.sampleDists,
#         clustering_distance_cols=all.parts.sampleDists,
#         col=colors,
#         cutree_cols = 3,
#         treeheight_row = 0)
#dev.off()
#
#
#
## tweek params ------------------------------------------------------------
##none of these seem to make much of a diff
#
#sex <- i.s2c %>% filter(samplePart == "SP") %>% 
#  filter(sampleName!="47745-moretoni-Chest" & Paired!="N") %>% 
#  mutate(path = paste("data/star_DoCounts/",fileName, sep = "")) %>% 
#  dplyr::rename(sample = sampleName)
#
#sex.hts<-DESeqDataSetFromHTSeqCount(sampleTable=sex, directory=directory, design=~PopPartColor)
#
#sex.dds<-DESeq(sex.hts)
#sex.res<-results(sex.dds)
#
##QC using PCA
#sex.vsd <- vst(sex.dds, blind=TRUE) #true or false blind doesnt effect the heatmap 
#rld <- rlog(sex.dds, blind=TRUE)
#
#sex.sampleDists <- dist(t(assay(sex.vsd[de,])))
#
#sex.sampleDistMatrix <- as.matrix(sex.sampleDists)
#rownames(sex.sampleDistMatrix) <- paste(sex.vsd$PopPartColor, sex.vsd$ind, sep="-")
##rownames(sampleDistMatrix) <- all.vsd$ind
#colnames(sex.sampleDistMatrix) <- NULL
##colors <- colorRampPalette( brewer.pal(9, "Greys") )(255)
#colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)
#
##pdf("output/heatmaps_2020/sex_only_all_t_de_genes.pdf", width = 12, height = 10)
#pheatmap(sex.sampleDistMatrix,
#         clustering_distance_rows=sex.sampleDists,
#         clustering_distance_cols=sex.sampleDists,
#         col=colors,
#         cutree_cols = 3,
#         treeheight_row = 0)
##dev.off()
#
#
#sex.sampleDists <- dist(t(assay(rld[de,])))
#
#sex.sampleDistMatrix <- as.matrix(sex.sampleDists)
#rownames(sex.sampleDistMatrix) <- paste(rld$PopPartColor, rld$ind, sep="-")
##rownames(sampleDistMatrix) <- all.vsd$ind
#colnames(sex.sampleDistMatrix) <- NULL
##colors <- colorRampPalette( brewer.pal(9, "Greys") )(255)
#colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)
#
##pdf("output/heatmaps_2020/sex_only_all_t_de_genes.pdf", width = 12, height = 10)
#pheatmap(sex.sampleDistMatrix,
#         clustering_distance_rows=sex.sampleDists,
#         clustering_distance_cols=sex.sampleDists,
#         col=colors,
#         cutree_cols = 3,
#         treeheight_row = 0)
##dev.off()
#
#
#ntd <- normTransform(sex.dds)
#library("vsn")
#meanSdPlot(assay(ntd))
#meanSdPlot(assay(sex.vsd))
#meanSdPlot(assay(rld))
#
#
#
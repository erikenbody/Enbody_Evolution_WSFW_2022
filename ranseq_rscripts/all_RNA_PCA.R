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



sample_id <- dir(file.path("data/star_DoCounts"))
kal_dirs <- file.path("data/star_DoCounts", sample_id)
directory <- "data/star_DoCounts/"
i.s2c <- read.table("output/sample_table_RSEM.csv", header = TRUE, stringsAsFactors = TRUE, sep = ",")
s2c <- i.s2c
s2c<-dplyr::filter(s2c,sampleName!="47745-moretoni-Chest") #remove apparent outlier - noted as very long feathers
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

all.vsd <- getVarianceStabilizedData(all.dds)
all.rv <- rowVars(all.vsd)
q95_wpn <- quantile( rowVars(all.vsd), .95) 
expr_normalized <- all.vsd[ all.rv > q95_wpn, ]



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

ggplot(pca2, aes(PC3, PC4,color=sampleColor,shape=samplePopulation)) +
  geom_point(size=6) +
  xlab(paste0("PC3: ",percentVar[3],"% variance")) +
  ylab(paste0("PC4: ",percentVar[4],"% variance")) + 
  coord_fixed()

all.vsd <- vst(all.dds, blind=TRUE)
#all.rld <- rlog(all.dds, blind=TRUE)


pcaData<-plotPCA2(all.vsd, intgroup = c("samplePopulation","samplePart"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, label=name,color=samplePopulation,shape=samplePart)) +
  geom_label_repel() +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#QC using PCA
all.vsd <- vst(all.dds, blind=TRUE)
pcaData<-plotPCA(all.vsd, intgroup = c("sampleColor"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, label=name,color=sampleColor)) +
  geom_label_repel() +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

plotCounts(all.dds,gene="WSFW005676",intgroup="PopPart")
plotCounts(all.dds, gene ="WSFW014970", intgroup = "PopPart", main = "CAMK4")
ESD.out <- plotCounts(all.dds, gene ="WSFW013094", intgroup = "PopPart", main = "ESD", returnData = T) 
write.csv(ESD.out, "output/esd_counts.csv", row.names = F)


kitlg.out <- plotCounts(all.dds,gene="WSFW009166",intgroup="PopPart", main = "KITLG", returnData = TRUE)
write.csv(kitlg.out, "output/kitlg_counts.csv", row.names = F)
egr1.out <- plotCounts(all.dds,gene="WSFW014426",intgroup="PopPart", main = "EGR1", returnData = TRUE)
write.csv(egr1.out, "output/egr1_counts.csv", row.names = F)




par(mfrow=c(1,3))
plotCounts(all.dds,gene="WSFW003303",intgroup="PopPart",main="androgen receptor")#, returnData = T) %>% View()
plotCounts(all.dds,gene="WSFW002320",intgroup="PopPart",main="aromatase")
plotCounts(all.dds,gene="WSFW001639",intgroup="PopPart",main="estrogen receptor")#, returnData = T) %>% View()


#some QC
res<-results(all.dds)

plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)

addmargins(table(filtering=(res$padj < .1)))
summary(res)
###
#
keep <- rowSums(counts(all.dds)) > 5
filtdds <- all.dds[keep,]
nrow(filtdds)


# heatmap -----------------------------------------------------------------
mypalette <- brewer.pal(11,"RdYlBu")

morecols <- colorRampPalette(mypalette)


vsd.df <- as.data.frame(assay(all.vsd) )
vsd.df$gene <- row.names(vsd.df)
t.vsd <- vsd.df %>% dplyr::select(ends_with("after"), gene)
t.vsd[timplat.de$gene_name, ] %>% View()

down.reg.t <- timplat.de %>% filter(log2FoldChange < 0)

subset.vsd <- assay(all.vsd)[down.reg.t$gene_name, ]
#subset.vsd <- assay(all.vsd)[timplat.de$gene_name, ]
#subset.vsd.df <- as.data.frame(subset.vsd)

heatmap.2(subset.vsd, colRow=NULL, scale="row", Colv = FALSE, density.info = "none",
          trace="none", dendrogram="column", margins=c(5, 10),key = TRUE,
          col=rev(morecols(50)),labRow=FALSE)

heatmap.2(assay(all.vsd), colRow=NULL,scale="row",Colv = FALSE,density.info = "none",
          trace="none", dendrogram="column", margins=c(5, 10),key = TRUE, ColSideColors = as.character(as.numeric(all.vsd$samplePopulation)),
          col=rev(morecols(50)),labCol=all.vsd$samplePart,labRow=FALSE)


subset.vsd <- assay(all.vsd)[timplat.overlap.dr$gene_name, ]


pheatmap::pheatmap(subset.vsd), cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, scale = "row", annotation_col = df,show_colnames = FALSE, border_color = NA)



my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)

#col=rev(morecols(50))

subset_names <- vsd.df %>% dplyr::select(contains("lorentzi-SP")) %>% names()
test <- subset.vsd[,subset_names]
test.col <- as.character(as.numeric(col.names(test)))

pheatmap(test)
heatmap.2(test, 
          colRow=NULL,
          scale="row",
          Colv = TRUE,
          density.info = "histogram",
          trace="none", 
          dendrogram="column", 
          margins=c(5, 10),
          key = TRUE,
          ColSideColors = c("1","1","1","2","2","2","2"),
          col=redgreen(75),
          labCol=as.factor(c("before","before","before","after","after","after","after")),
          labRow=FALSE)

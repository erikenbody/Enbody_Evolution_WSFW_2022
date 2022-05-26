suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(reshape)))
suppressWarnings(suppressMessages(library(ggrepel)))

#########FOR RUNNING THIS ON COMMAND LINE
option_list <- list(make_option(c('-t','--type'), action='store', type='character', default=NULL, help='Color with all parts or between parts: AllParts, Chest, SP, Crown'),
                    make_option(c('-c','--colorA'), action='store', type='character', default=NULL, help='colorA: White, Black, Brown'),
                    make_option(c('-d','--colorB'), action='store', type='character', default=NULL, help='colorB: White, Black, Brown'),
                    make_option(c('-p','--popA'), action='store', type='character', default=NA, help='popA: aida, naimii, lorentzi, moretoni'),
                    make_option(c('-q','--popB'), action='store', type='character', default=NA, help='popB: aida, naimii, lorentzi, moretoni')
)
opt <- parse_args(OptionParser(option_list = option_list))

#################################################################################
suppressWarnings(source("annotation_import.r")) #generic import of annotations for this type of thing I put together

#uncomment for testing. access variables from opt object
#opt$colorA <- "White"
#opt$colorB <- "Brown"
#opt$type   <- "SP"
#opt$popA   <- "lorentzi"
#opt$popB   <- "lorentzi"
#consider levels (direction of DE)

outdir <- "output/deseq_DE_output_fdr_0.05/"
if(!file.exists(outdir)){ 
  dir.create(outdir)
  print("Ouptput directory created.")
} else {
  print("Ouptput directory exists, appending and adding new files.")
}

#pull in functional information
annie.final <- annie.final %>% select(-dup, -dup2, -combined)

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
s2c <- s2c %>% filter(Sex != "M" & Treatment != "After")

if (opt$type != "AllParts"){
  s2c <- s2c %>% filter(samplePart == opt$type)
}

#set up conditions for t treatment analysi
if (!is.na(opt$popA) & opt$popA == "lorentzi" & opt$popB == "lorentzi" & opt$type == "SP"){
  s2c <- i.s2c %>% filter(samplePopulation == "lorentzi" & samplePart == "SP" & Paired == "Y") %>% 
    mutate(ind = as.factor(ind))
  s2c <- arrange(s2c, ind, Treatment)
}  

#if (!is.null(opt$popA)){
#  s2c <- s2c %>% filter(samplePopulation == opt$popA | samplePopulation == opt$popB)
#} 

# Set up color contrast desired
s2c <- s2c %>% filter(sampleColor == opt$colorA | sampleColor == opt$colorB)
#setup factor name for deseq
comp <- paste("PopPart", opt$popA, opt$type, opt$popB, opt$type, sep = "_")
if (!is.na(opt$popA)){
  out_prefix <- paste(outdir, opt$popA, "_", opt$popB, "_", opt$type, "_", opt$colorA, "_", opt$colorB, sep = "")
} else {
  out_prefix <- paste(outdir, opt$type, "_", opt$colorA, "_", opt$colorB, sep = "")
}

# run DESEQ2 --------------------------------------------------------------

if (is.na(opt$popA)){ #if no population defined, run the comparison between parts/color
  deseq_hts<-DESeqDataSetFromHTSeqCount(sampleTable=s2c, directory=directory, design=~PartColor)
  deseq_dds<-DESeq(deseq_hts)
  resultsNames(deseq_dds)
  deseq_res<-results(deseq_dds,name=resultsNames(deseq_dds)[2], alpha = 0.05)
} else if (opt$popA == "lorentzi" & opt$popB == "lorentzi" & opt$type == "SP"){ #if the =
  deseq_hts<-DESeqDataSetFromHTSeqCount(sampleTable=s2c, directory=directory, design=~ind + PartColor)
  deseq_dds<-DESeq(deseq_hts)
  resultsNames(deseq_dds)
  deseq_res<-results(deseq_dds,name=resultsNames(deseq_dds)[4], alpha = 0.05)
} else {
  deseq_hts<-DESeqDataSetFromHTSeqCount(sampleTable=s2c, directory=directory, design=~PopPart)
  deseq_dds<-DESeq(deseq_hts)
  
  #pull out comparison of interest
  compA <- paste(opt$popA, opt$type, sep = "_")
  compB <- paste(opt$popB, opt$type, sep = "_")
  
  deseq_res<-results(deseq_dds, c("PopPart", compA, compB), alpha = 0.05)
}

#plotCounts(deseq_dds, gene = "WSFW009166", intgroup="PartColor")

#plotMA(deseq_res)
# output results ----------------------------------------------------------
# 
#first output raw counts
deseq_sf <- estimateSizeFactors(deseq_dds)
deseq_counts<-counts(deseq_sf, normalized=TRUE)
write.csv(deseq_counts,paste(out_prefix, "counts.csv", sep = "_"),row.names=TRUE)

#DEGs
#no filtering on pvalue
deseq_res$gene_name<-row.names(deseq_res)
deseq_res.annie<-merge(as.data.frame(deseq_res),annie.final,by.x="gene_name",by.y="gene_name",all.x=TRUE,all.y=FALSE)
deseq_res.annie$blast_name[is.na(deseq_res.annie$blast_name)]<-"unknown"
fwrite(deseq_res.annie,paste(out_prefix,"DE_genes.txt", sep = "_"), quote = FALSE, sep = "\t", row.names = FALSE)

#basic QC
png(filename = paste(out_prefix, "QC.png", sep = "_"))
plot(metadata(deseq_res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
abline(v=metadata(deseq_res)$filterTheta)
lines(metadata(deseq_res)$lo.fit, col="red")
dev.off()

#simple PCA
if (is.na(opt$popA) | opt$popA == "lorentzi" & opt$popB == "lorentzi" & opt$type == "SP"){
  pca.comp <- "PartColor"
} else {
  pca.comp <- "PopPart"
}
deseq_vsd <- vst(deseq_dds, blind=FALSE)
pcaData<-plotPCA(deseq_vsd, intgroup = c(pca.comp),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

if (!is.na(opt$popA)){
  pcaData <- pcaData %>% filter(grepl(opt$popA, name) | grepl(opt$popB, name))
} 
png(filename = paste(out_prefix, "PCA.png", sep = "_"))
ggplot(data = pcaData, aes_string("PC1", "PC2", label="name",color=pca.comp)) +
    geom_label_repel() +
    geom_point(size=5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() 
dev.off()

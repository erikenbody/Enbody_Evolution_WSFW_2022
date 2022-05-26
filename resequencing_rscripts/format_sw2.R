#!/usr/bin/env Rscript
library(dplyr)
args <- commandArgs(trailingOnly = TRUE) #apparently trailing is critial here
input<-args[1]

df.mafs <- read.table(input, header = T)

prefix <- gsub(".mafs.gz", "", basename(input))

df.sf <- df.mafs %>% 
  mutate(x = round(unknownEM * nInd,0)) %>% 
  mutate(folded = ifelse(anc == "N", 1, 0)) %>% 
  filter(!(x==0 & folded == 0)) %>% #remove polarized monomorphic sites
  select(position, x, nInd, folded, n = nInd)

write.table(df.sf, paste0(prefix, "_sf_input.txt"),
            col.names = F, quote = F, row.names = F, sep = "\t")


for (contig in unique(df.mafs$chromo)){
  
  dir.create(prefix, showWarnings = F)
  df.mafs.f <- df.mafs %>%  filter(chromo == contig)
  df.mafs.f %>% 
    mutate(x = round(unknownEM * nInd,0)) %>% 
    mutate(folded = ifelse(anc == "N", 1, 0)) %>% 
    filter(!(x==0 & folded == 0)) %>% #remove polarized monomorphic sites
    select(position, x, nInd, folded, n = nInd) %>% 
    write.table(paste0(prefix,"/",prefix,"_",contig, "_sf_input.txt"),
              col.names = F, quote = F, row.names = F, sep = "\t")
}

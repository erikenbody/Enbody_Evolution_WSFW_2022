# satsuma alignment -------------------------------------------------------

#wsfw.zefi.chain.ens <- fread("data/satsuma_chain_file/ensembl/satsuma_summary.chained.out_sorted", sep = "\t", fill = T)
#
#names(wsfw.zefi.chain.ens) <- c("qName", "qStart", "qStop", "zefi_chr", "refStart", "refStop", "V7", "dir")
#wsfw.zefi.chain.ens.out <- wsfw.zefi.chain.ens %>% 
#  filter(zefi_chr!="") %>% 
#  separate(zefi_chr, c("A"), remove = FALSE, sep = "_", extra = "drop") %>% 
#  mutate(refName = paste("chromosome", A, sep = "_")) %>% 
#  select(-A) %>% 
#  mutate(qLength = qStop - qStart)
#
##get the maxmatch summary
#wsfw.chain.df <- ddply(wsfw.zefi.chain.ens.out,
#                .(qName,refName),
#                summarize,
#                MaxMatch=sum(qLength),
#                refStart=min(refStart))
#
#wsfw.chain <- arrange(wsfw.chain.df,refName,refStart,MaxMatch) %>% 
#  mutate(refName = ifelse(MaxMatch < 5000, "NA", refName))
#
##select only the place where the whole contig matches best
#wsfw.chain <- ddply(wsfw.chain,.(qName),function(e){                
#  a <- subset(e,refName==e$refName[e$MaxMatch==max(MaxMatch)])
#  if(nrow(a)==1){
#    a
#  }
#})
#
##unique(wsfw.chain$refName)
#####wsfw.chain$refName <- gsub("chromosome_","", wsfw.chain$refName)
#
#write.csv(wsfw.chain, "data/satsuma_chain_file/ensembl/FORMATTED_satsuma_summary.chained_ensembl_V2.out", row.names = F)
#
wsfw.chain<- read.csv("data/satsuma_chain_file/ensembl/FORMATTED_satsuma_summary.chained_ensembl_V2.out")

chr_order <- c("chromosome_1", "chromosome_1A", "chromosome_2", "chromosome_3", "chromosome_4", 
               "chromosome_4A", "chromosome_5", "chromosome_6", "chromosome_7", "chromosome_8", 
               "chromosome_9", "chromosome_10", "chromosome_11", "chromosome_12", "chromosome_13", 
               "chromosome_14", "chromosome_15", "chromosome_16", "chromosome_17", "chromosome_18", 
               "chromosome_19", "chromosome_20", "chromosome_21", "chromosome_22", "chromosome_23", 
               "chromosome_24", "chromosome_25", "chromosome_26", "chromosome_27", "chromosome_28", 
               "chromosome_LG2", "chromosome_LG5","chromosome_LGE22", "chromosome_Z", "chromosome_MT", 
               "chromosome_1B", "chromosome_Un")

#could this possible be improved with incorporating the orientation of each scaffold?

simple_labels <- c("chromosome_1", "chromosome_1A", "chromosome_2", "chromosome_3", "chromosome_4", 
                   "chromosome_4A", "chromosome_5", "chromosome_6", "chromosome_7", "chromosome_8", 
                   "chromosome_9", "chromosome_10", "chromosome_11", "chromosome_12", "chromosome_13", 
                   "chromosome_14", "chromosome_15", "chromosome_16", "chromosome_17", "chromosome_18", 
                   "chromosome_19", "chromosome_20", "chromosome_Z", "chromosome_Un")


# inport fst --------------------------------------------------------------

input_fst<-function(x,y){
  pops<-x
  colnames(pops)<-c("region","scaff","midPos","Nsites","fst")
  zpops<-y
  colnames(zpops)<-c("region","scaff","midPos","Nsites","fst")
  
  pops_df <- rbind(pops,zpops)
  
  sl50<-read.table("data/autosome_scaffolds_gr_50000.txt",header=TRUE)
  pops_df<-merge(pops_df,sl50,all.x=FALSE,all.y=FALSE,by="scaff")
  out<-strsplit(as.character(pops_df$scaff), "_")
  pops_df$scaff_num<-as.numeric(lapply(out,function(ss){ss[2]}))
  #remove those with < 2 windows and < 10 SNPs
  pops_df<-pops_df[as.numeric(ave(pops_df$scaff_num, pops_df$scaff_num, FUN=length)) >= 2, ]
  pops_df<-dplyr::filter(pops_df, Nsites>10)
  #remove any NA
  pops_df<-pops_df[complete.cases(pops_df),]
  
  pops.chr <- suppressWarnings(left_join(pops_df, wsfw.chain, by = c("scaff" = "qName")))
  
  pops.chr <- pops.chr %>% mutate(refName = ifelse(!refName %in% chr_order, "chromosome_Un", as.character(refName))) %>% 
    dplyr::rename(chr = refName)
  #chr = paste("chromosome_", refName2, sep = "")) #%>% 
  #select(-refName) %>% dplyr::rename(refName = refName2)
  
  pops.chr <- pops.chr %>% 
    mutate(chr_ordered = factor(chr, levels = chr_order)) %>% 
    arrange(chr, refStart) %>% 
    dplyr::rename(zefi_start = refStart)
  
  #then do zfst seperately
  zdf <- pops.chr %>% filter(chr == "Z") %>% mutate(zfst = (fst - mean(fst)) / sd(fst))
  adf <- pops.chr %>% filter(chr != "Z" | is.na(chr)) %>% mutate(zfst = (fst - mean(fst)) / sd(fst))
  
  chr_pops_comb <- rbind(adf, zdf)
  
  return(chr_pops_comb)
  
}


# plottingm anhattan----------------------------------------------------------------

manc.pbs<-function(fst){
  
  fst <- fst %>% arrange(chr_ordered, zefi_start)
  fst$row<-1:nrow(fst)
  fst$pbsrollmean <- zoo::rollmean(fst$pbs,50,fill=NA)
  
  chr_breaks <- fst %>% group_by(chr_ordered) %>% 
    dplyr::summarise(chr_breaks = mean(row)) %>% 
    mutate(chr_ordered = as.factor(chr_ordered),
           chr_labels = ifelse(chr_ordered %in% simple_labels, as.character(chr_ordered), ".")) %>% 
    mutate(chr_labels = gsub("chromosome_", "", chr_labels))
  
  olf<-filter(fst,chr_ordered != "chromosome_Z" & chr_ordered!="NA")
  olf3<-filter(olf,pbs>=quantile(pbs,.999,na.rm=T))
  
  olz<-filter(fst,chr_ordered=="chromosome_Z")
  olz3<-filter(olz,pbs>=quantile(pbs,.999,na.rm=T))
  
  #FST
  fst_main <- fst %>% 
    ggplot(aes(x = row, y = pbs, col = chr_ordered)) + theme_bw()+
    scale_color_manual(values=rep(c("grey30","grey80"), length(levels(factor(fst$chr_ordered)))/2+1))+
    #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1), minor_breaks = NULL)+
    geom_point(size=0.4,shape=21,stroke=0.2) +
    labs(y="PBS") +
    #scale_x_continuous(breaks = chr_breaks$chr_breaks, 
    #                   labels = chr_breaks$chr_labels) +
    scale_x_continuous(breaks = chr_breaks$chr_breaks, #https://stackoverflow.com/questions/50399838/how-to-alternate-a-new-line-for-overlapping-x-axis-labels
                       labels = function(labels) {
                         sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                       }) +
    theme(legend.position="none",
          panel.border=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x = element_text(color = "black", vjust = 0.99, size = 12),
          panel.grid = element_blank(),
          panel.grid.major.y=element_line(color="grey60",size=0.2),
          panel.grid.minor.y=element_line(color="grey60",size=0.1),
          axis.title.y = element_text(size=12),
          axis.text = element_text(size=7),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2))
  
  fst_out <- fst_main +
    geom_point(data=olf3,shape=21,stroke=0.4,size=0.5,col="red")+
    geom_point(data=olz3,shape=21,stroke=0.4,size=0.5,col="red")
  
  
  return(fst_out)
}


# manhattan with gene labels ----------------------------------------------

manc.labels.pbs<-function(fst, genes){
  fst <- fst %>% arrange(chr_ordered, zefi_start)
  fst$row <- 1:nrow(fst)
  fst$pbsrollmean <- zoo::rollmean(fst$pbs,50,fill=NA)
  
  chr_breaks <- fst %>% group_by(chr_ordered) %>% 
    dplyr::summarise(chr_breaks = mean(row)) %>% 
    mutate(chr_ordered = as.factor(chr_ordered),
           chr_labels = ifelse(chr_ordered %in% simple_labels, as.character(chr_ordered), ".")) %>% 
    mutate(chr_labels = gsub("chromosome_", "", chr_labels))
  
  olf <- filter(fst,chr_ordered != "chromosome_Z" & chr_ordered!="NA")
  olf3 <- filter(olf,pbs >= quantile(pbs,.999, na.rm=T))
  
  olz <- filter(fst,chr_ordered == "chromosome_Z")
  olz3 <- filter(olz,pbs >= quantile(pbs,.999, na.rm=T))
  
  genes.interest <- genes %>% distinct(blast_name, .keep_all = TRUE) %>% 
    mutate(blast_name = ifelse(is.na(blast_name), as.character(Name), as.character(blast_name)))
  
  gr1.x <- GRanges(genes.interest$seqid, ranges = IRanges(genes.interest$start, width = 1))
  gr2.x <- GRanges(fst$scaff, ranges = IRanges(fst$midPos, width = 1))
  
  fo <- suppressWarnings(findOverlaps(gr1.x, gr2.x, maxgap = 50000))
  
  genes.fst <- cbind(genes.interest[queryHits(fo), ], fst[subjectHits(fo), c("row", "pbs")])
  genes.fst <- genes.fst %>% select(seqid, start, chr, blast_name, Name, row, pbs) %>% 
    group_by(Name) %>% filter(pbs == max(pbs))
  
  genes.fst <- genes.fst %>% mutate(pbs = ifelse(blast_name == "Lgr6", 0.482, pbs))
  
  better_label <- genes.fst %>% 
    group_by(seqid) %>% 
    dplyr::summarise(blast_name = paste(blast_name, collapse = " \n "),
                     row = round(mean(row), 0),
                     pbs = max(pbs))
  
  better_label <- better_label %>% mutate(pbs = ifelse(pbs < quantile(fst$pbs,.999, na.rm=T), max(pbs), pbs))
  
  #FST
  fst_main <- fst %>% 
    ggplot(aes(x = row, y = pbs, col = chr_ordered)) + theme_bw()+
    scale_color_manual(values=rep(c("grey30","grey80"), length(levels(factor(fst$chr_ordered)))/2+1))+
    #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1), minor_breaks = NULL)+
    geom_point(size=1,shape=19,stroke=0.2) +
    labs(y="PBS") +
    #scale_x_continuous(breaks = chr_breaks$chr_breaks, 
    #                   labels = chr_breaks$chr_labels) +
    scale_x_continuous(breaks = chr_breaks$chr_breaks, #https://stackoverflow.com/questions/50399838/how-to-alternate-a-new-line-for-overlapping-x-axis-labels
                       labels = function(labels) {
                         sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                       }) +
    theme(legend.position="none",
          panel.border=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x = element_text(color = "black", vjust = 0.99, size = 18),
          panel.grid = element_blank(),
          panel.grid.major.y=element_line(color="grey60",size=0.2),
          panel.grid.minor.y=element_line(color="grey60",size=0.1),
          axis.title.y = element_text(size=12),
          axis.text = element_text(size=20),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    geom_text_repel(data = better_label, aes(x = row, y = pbs, label = blast_name), 
                    color = "black", 
                    nudge_y = 0.5, 
                    #force=1, point.padding=unit(1,'lines'),
                    hjust=1,
                    direction='both',
                    #nudge_x = 15,
                    segment.size=0.2,
                    size = 6) 
  
  fst_out <- fst_main +
    geom_point(data=olf3,shape=19,stroke=0.4,size=1,col="red")+
    geom_point(data=olz3,shape=19,stroke=0.4,size=1,col="red")
  
  return(fst_out)
}

#dfA <- input_fst(bc.f, bc.z)
#manc(dfA) + ylim(0,0.75)
#manc.labels(dfA, "output/chest_genes_of_interest.csv") + ylim(0, 0.75)
#


# manhattan for a single chrom --------------------------------------------


PBS.filt_manc<-function(fst, chr_idx){
  
  fst <- fst %>% filter(chr == chr_idx)
  
  fst$chr <- ifelse(!grepl("chromosome",fst$chr), "chromosome_unknown", fst$chr)
  
  #specific to the vertebrate genomes unpublished ZEFI assembly
  chr_order <- c("chromosome_1", "chromosome_1A", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_4A", "chromosome_5", "chromosome_6", "chromosome_7", "chromosome_8", "chromosome_9", "chromosome_10", "chromosome_11", "chromosome_12", "chromosome_13", "chromosome_14", "chromosome_15", "chromosome_16", "chromosome_17", "chromosome_18", "chromosome_19", "chromosome_20", "chromosome_21", "chromosome_22", "chromosome_23", "chromosome_24", "chromosome_25", "chromosome_26", "chromosome_27", "chromosome_28", "chromosome_29", "chromosome_30", "chromosome_Z", "chromosome_unknown")
  
  fst$chr_ordered <- factor(fst$chr, levels = chr_order)
  
  fst<-fst%>%arrange(chr_ordered, zefi_start)
  fst$row<-1:nrow(fst)
  fst$pbsrollmean <- zoo::rollmean(fst$pbs,50,fill=NA)
  
  fst$chr_labels <- gsub("chromosome_", "", fst$chr_ordered)
  chr_breaks <- fst %>% group_by(chr_labels) %>% dplyr::summarise(chr_breaks = mean(row))
  
  olf<-filter(fst,chr!="chromosome_Z")
  olf2<-filter(olf,pbs>=quantile(pbs,.99,na.rm=T))
  olf3<-filter(olf,pbs>=quantile(pbs,.999,na.rm=T))
  
  olz<-filter(fst,chr=="chromosome_Z")
  olz2<-filter(olz,pbs>=quantile(pbs,.99,na.rm=T))
  olz3<-filter(olz,pbs>=quantile(pbs,.999,na.rm=T))
  
  #FST
  fst_main <- fst %>% 
    ggplot(aes(x=row,y=fst,col=chr_ordered))+theme_bw()+
    theme(legend.position="none",
          panel.border=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, color = "black"),
          panel.grid = element_blank(),
          panel.grid.major.y=element_line(color="grey60",size=0.2),
          panel.grid.minor.y=element_line(color="grey60",size=0.1),
          axis.title.y = element_text(size=12),
          axis.text = element_text(size=7),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    scale_color_manual(values=rep(c("grey30","grey80"),length(levels(factor(fst$chr_ordered)))/2+1))+
    #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL)+
    geom_point(size=0.4,shape=21,stroke=0.2)+
    geom_line(aes(y=fstrollmean),lwd=0.15,col="black") +
    labs(y="PBS") +
    scale_x_continuous(breaks=chr_breaks$chr_breaks, 
                       labels = chr_breaks$chr_labels)
  
  fst_out <- fst_main +
    geom_point(data=olf2,shape=21,stroke=0.4,size=0.5,col="orange")+
    geom_point(data=olz2,shape=21,stroke=0.4,size=0.5,col="orange")+
    geom_point(data=olf3,shape=21,stroke=0.4,size=0.5,col="red")+
    geom_point(data=olz3,shape=21,stroke=0.4,size=0.5,col="red")
  
  return(fst_out)
}


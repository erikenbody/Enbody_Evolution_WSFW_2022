library(stringr)
library(tidyr)
library(reshape)

#includes low confidence calls from gff
gene.df<-dget("data/gene_only.db")
fil.gene.df<-dplyr::select(gene.df,seqid,start,end,strand,Name,Note)
fil.gene.df$gff_gene_name<-word(fil.gene.df$Note,3)
fil.gene.df$gff_gene_name<-gsub(":","",fil.gene.df$gff_gene_name,fixed = TRUE)
fil.gene.df$gene_length<-fil.gene.df$end - fil.gene.df$start
gene_length<-dplyr::select(fil.gene.df,Name,gene_length)

#Melanin genes
mel<-read.csv("data/melano_genes.csv",header=TRUE)
mel$gene_name<-as.character(mel$gene_name)
mel<-mel[complete.cases(mel$gene_name),]

#go terms from IPR + Annie
wsfw_goterms<-read.delim("data/WSFW_maker_only_ipr.annie",header=FALSE)
wsfw_goterms$V2<-NULL
wsfw_goterms$go2<-gsub("|","_",wsfw_goterms$V3,fixed = TRUE)
#sep.go<-read.table(text=wsfw_goterms$go2,sep="_")

wsfw_goterms<-separate(wsfw_goterms, 'go2', paste("go", 1:12, sep="_"), sep="_", extra="drop")
wsfw_goterms$V3<-NULL
wsfw_goterms<-separate(wsfw_goterms, 'V1', paste("name", 1:2, sep="_"),sep="-", extra="drop")
wsfw_goterms$name_2<-NULL
wsfw_goterms.m<-melt(wsfw_goterms, id.var = 'name_1', variable.name = 'go')
wsfw_goterms.m$variable<-NULL
wsfw_goterms.m$go<-wsfw_goterms.m$value ; wsfw_goterms.m$value<-NULL
wsfw_goterms_final<-wsfw_goterms.m[!is.na(wsfw_goterms.m$go),]

#Three versions of Annie blast output, by order of increasing strictness. Last one is that advocated by MAKER developers
#annie_filt<-read.delim("data/annotation_files/WSFW_maker_blast_0.01.annie",header=FALSE)
#annie_filt<-read.delim("data/annotation_files/WSFW_maker_blast_1e-6.annie",header=FALSE)
annie_filt<-read.delim("data/WSFW_maker_blast_sprot_np_STRICT.annie",header=FALSE)
colnames(annie_filt)<-c("gene_name","type","blast")


#format for use as annotation
annie_filt.gene<-dplyr::filter(annie_filt,type=="name")
filtannie<-annie_filt.gene %>% tidyr::separate(blast, c("blast_name", "transcript_num"),"_")
filtannie$combined<-paste(filtannie$gene_name,filtannie$blast_name,sep="_")
#remove duplicates
filtannie.dups<-filtannie[filtannie$combined %in% filtannie$combined[!duplicated(filtannie$combined)],]
filtannie.drm<-filtannie[!duplicated(filtannie$combined), ]
#need to mark as duplicate to note later that this gene shows up twice. also for easy filtering
#https://stackoverflow.com/questions/38315683/how-do-i-mark-duplicates-in-a-new-column
filtannie.drm$dup <- ave(filtannie.drm$blast_name, filtannie.drm$blast_name, FUN = length) > 1L
filtannie.drm$gene_name<-as.character(filtannie.drm$gene_name)
filtannie.drm$dup2 <- ave(filtannie.drm$gene_name, filtannie.drm$gene_name, FUN = length) > 1L

annie.final <- filtannie.drm

#make a transcript index file
tx_annie_filt<-annie_filt %>% dplyr::rename(tx_name=gene_name) %>% 
  tidyr::separate(tx_name, c("gene_name","ext"),"-",remove=FALSE) %>%
  dplyr::select(gene_name,tx_name) %>%
  dplyr::left_join(annie.final,by=c("gene_name")) %>%
  dplyr::distinct(tx_name,.keep_all = TRUE) %>%
  dplyr::select(tx_name,gene_name,blast_name) %>%
  dplyr::rename(target_id=tx_name,ext_gene=blast_name)


#feather growth
#fg<-read.csv("data/annotation_files/base_tip_diff_exp.csv")
#genes2<-fg$ID
#fg.conv<-genes(EnsDb.Hsapiens.v79, filter=list(GenenameFilter(genes2),GeneIdFilter("ENSG", "startsWith")),
#return.type="data.frame", columns=c("gene_id"))

##BELOW ARE FROM NG ET AL 2015 BMC
#from feather study appendix
feathers_just_body_bot_top<-read.csv("data/feather_growth_gene_ID.csv")

#genes specifically mentioned in text related to feather development
deg_body_bottom_top<-c("PLXNA1", "NRP1", "DPYSL3", "MAPK1", "CDK5","SP3", "NRP1", "EPAS1", "TP63", "BMPR-II", "BMPR1A","SMAD5", "MADH2", "BMPR1A", "NRP1", "BMPR-II", "EGFR","BMPRIA","BMPRIB","EGF")
deg_body_bottom_top<-as.data.frame(deg_body_bottom_top)
deg_body_bottom_top$part<-c("deg_body_bottom_top")
deg_body_bottom_top$gene_id<-deg_body_bottom_top$deg_body_bottom_top ; deg_body_bottom_top$deg_body_bottom_top<-NULL

deg_body_flight<-c("PRKAR1A", "LMO4", "TP63", "TWSG1", "JAG1", "FGFR2", "ACVR1", "CA2","SMAD5", "JAG1", "FGFR2", "HIF1A", "ACVR1", "WNT7B","JAG1","FGF","FGFR1","FGFR2","FGFR3")
deg_body_flight<-as.data.frame(deg_body_flight)
deg_body_flight$part<-c("deg_body_flight")
deg_body_flight$gene_id<-deg_body_flight$deg_body_flight ; deg_body_flight$deg_body_flight <-NULL

deg_flight_bottom_top<-c("JUN", "DKK3", "WNT5A", "CDH3", "CD44", "CDH5", "PPARD", "WNT2B", "SFRP4", "SOX14", "SFRP2", "PPP2R2B", "SOX7", "FRZB", "RARB","JUN", "INHBA", "SMAD2", "SMAD1", "RUNX3", "MAPK11", "INHBB","MGP", "GLI1", "WWOX", "IGFBP5", "GJA5", "SHH", "SOX14", "SMAD1", "CBFB","WNT5A","WNT5B","WNT6","LMO4", "LIPA", "SHH", "EDNRA", "RARB", "GJA5", "CRH", "GLI1","L-CAM", "N-CAM","ADIPOQ", "COL1A2", "COL3A1", "COL4A1", "COL6A1", "COL6A2", "COL8A1","DKK2","CD44","TIMP3","CRISP1","FSTL1","ADRA1A", "ITPR3", "PLA2G4A", "ACTG2", "ARHGEF12", "PLCB4", "PLA2G12A", "ITPR2", "PLA2G10", "RAMP2","PLA2G4C")
deg_flight_bottom_top<-as.data.frame(deg_flight_bottom_top)
deg_flight_bottom_top$part<-c("deg_flight_bottom_top")
deg_flight_bottom_top$gene_id<-deg_flight_bottom_top$deg_flight_bottom_top ; deg_flight_bottom_top$deg_flight_bottom_top<-NULL

deg_flight_calamus<-c("INHBA", "RUNX3", "PMEPA1", "RUNX2", "INHBB","TUBB3", "TUBA1B", "CLDN4", "TJP3", "JAM3", "ACTN1", "MTMR2","TUBB3", "TUBA1B", "CLDN4", "TJP3", "JAM3", "ACTN1", "MTMR2")
deg_flight_calamus<-as.data.frame(deg_flight_calamus)
deg_flight_calamus$part<-"deg_flight_calamus"
deg_flight_calamus$gene_id<-deg_flight_calamus$deg_flight_calamus ; deg_flight_calamus$deg_flight_calamus <-NULL

axon<-c("MYL4", "CDK5",
        "SEMA4B", "PRKAR1A", "NFATC3", "PLXNA1", "PLXNB2", "ARPC4", "NRP1", "MAPK1", "WNT6", "PRKCI", "RASA1", "ECE2", "SEMA5A")
axon<-as.data.frame(axon)
axon$part<-c("axon")
axon$gene_id<-axon$axon ; axon$axon <-NULL

fe_cand<-rbind(deg_body_bottom_top,deg_body_flight,deg_flight_bottom_top,deg_flight_calamus,axon)

###ANDROGENS
###Paper: Using the canary genome to decipher the evolution of hormone-sensitive gene regulation in seasonal singing birds
andrC<-read.csv("data/canary_T_genes.csv",header=TRUE)
andrE<-read.csv("data/euro_robin_T_genes.csv",header=TRUE)
cangened_andr<-rbind(andrC,andrE)


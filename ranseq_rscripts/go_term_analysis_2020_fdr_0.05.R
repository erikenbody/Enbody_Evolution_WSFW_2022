library(annotables)
library(clusterProfiler)
library(org.Gg.eg.db) #chicken dataframe with gene links to ens
library(org.Hs.eg.db) #human df with gene links to ens
library(tidyverse)
library(enrichplot)
library(goseq)
library(GO.db)
#inspiration:
#https://github.com/hbctraining/DGE_workshop/blob/master/lessons/09_functional_analysis.md

df.de.venn <- read.csv("output/venn_diagram/overlap_29_genes_all_three_cats_fdr_0.05.csv") %>% mutate(symbol = toupper(blast_name)) #uppercase for matching chicken
df.de.t <- read.csv("output/venn_diagram/t_implant_de_genes_only_fdr_0.05.csv") %>% mutate(symbol = toupper(blast_name)) 
df.de.ssp <- read.csv("output/venn_diagram/overlap_genes_t_treat_subspecies_More_Lore_fdr_0.05.csv") %>% mutate(symbol = toupper(blast_name)) 
df.background <- read.csv("output/go_terms_2020/all_T_genes_for_background_fdr_0.05.csv") %>% mutate(symbol = toupper(blast_name)) 

#df.de.t %>% filter(!symbol %in% galgal5$symbol) %>% nrow() #56 genes not found in galgal5
#df.de.t %>% filter(!symbol %in% grch37$symbol) %>% nrow() #29 genes not found in human
#df.background %>% filter(!symbol %in% galgal5$symbol) %>% nrow() #6318
#df.background %>% filter(!symbol %in% grch37$symbol) %>% nrow() #3052

wsfw2esn <- function(df.gene.list, ens.sp){

  if(ens.sp == "chicken"){ 
    ens.sp.db <- galgal5
  } else if(ens.sp == "human") {
    ens.sp.db <- grch37
  }
  
  df.ens.sp <- ens.sp.db %>% filter(symbol %in% df.gene.list$symbol)
  df.ens.sp.nd <- which(duplicated(df.ens.sp$symbol) == FALSE)
  df.ens.sp <- df.ens.sp[df.ens.sp.nd, ] 
  df.gene.list.ens <- inner_join(df.gene.list, df.ens.sp, by=c("symbol"="symbol")) %>% 
    distinct(symbol, .keep_all = T) #beacuse sometimes same blast hit multiple genes (expansions, possiblye)
  return(df.gene.list.ens)
}

wsfw.go <- function(df.ens.de.genes, ens.sp){

  if(ens.sp == "chicken"){ 
    go.sp.db <- org.Gg.eg.db
    df.background.go <- df.background.chick
  } else if(ens.sp == "human") {
    go.sp.db <- org.Hs.eg.db
    df.background.go <- df.background.human
  }
  ego.raw <- enrichGO(gene = df.ens.de.genes$ensgene, 
                          universe = df.background.go$ensgene,
                          keyType = "ENSEMBL",
                          OrgDb = go.sp.db, 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05, 
                          readable = TRUE)
  return(ego.raw)
  #cluster_summary <- data.frame(ego.raw)
  #return(cluster_summary)
}



# setup files -------------------------------------------------------------


df.de.venn.chick <- wsfw2esn(df.de.venn, "chicken")
df.de.venn.human <- wsfw2esn(df.de.venn, "human")
df.de.t.chick <- wsfw2esn(df.de.t, "chicken")
df.de.t.human <- wsfw2esn(df.de.t, "human")

df.de.ssp.chick <- wsfw2esn(df.de.ssp, "chicken")
df.de.ssp.human<- wsfw2esn(df.de.ssp, "human")

df.background.chick <- wsfw2esn(df.background, "chicken")
df.background.human <- wsfw2esn(df.background, "human")



# go analysis -------------------------------------------------------------

df.go.venn.chick <- wsfw.go(df.de.venn.chick, "chicken")
df.go.venn.human <- wsfw.go(df.de.venn.human, "human")

df.go.t.chick <- wsfw.go(df.de.t.chick, "chicken")
df.go.t.human <- wsfw.go(df.de.t.human, "human")

df.go.ssp.chick <- wsfw.go(df.de.ssp.chick, "chicken")
df.go.ssp.human <- wsfw.go(df.de.ssp.human, "human")


# plots -------------------------------------------------------------------

png("output/go_terms_2020/go_venn_human_fdr_0.05.png", height = 8, width = 15, units = "in", res = 300)
clusterProfiler::dotplot(df.go.venn.human)
dev.off()

png("output/go_terms_2020/go_t_treatment_human_fdr_0.05.png", height = 8, width = 15, units = "in", res = 300)
clusterProfiler::dotplot(df.go.t.human)
dev.off()

png("output/go_terms_2020/go_t_treatment&subspecies_human_fdr_0.05.png", height = 8, width = 15, units = "in", res = 300)
clusterProfiler::dotplot(df.go.ssp.human)
dev.off()


# goseq -------------------------------------------------------------------

# shared t implant and white brown ----------------------------------------

df.in.all <- read.csv("output/all_comparison_DE_output_fdr_0.05.csv")
unique(df.in.all$comparison)
df.X1 <- df.in.all %>% filter(comparison == "lorentzi_lorentzi_SP_White_Brown_DE_genes" & padj < 0.05)
df.X2 <- df.in.all %>% filter(comparison == "SP_White_Brown_DE_genes" & padj < 0.05)

#some values for text
df.X1 %>% nrow() #number of DE genes from T
df.X2 %>% nrow() #number of DE genes between all white and brown SP feathers
df.X2 %>% filter(grepl("KRF", blast_name)) %>% nrow() #number of KRF genes
df.X2 %>% filter(gene_name %in% df.X1$gene_name) %>% nrow() #overlapping
df.X2 %>% filter(gene_name %in% df.X1$gene_name) %>% filter(log2FoldChange > 0) %>% nrow() #overlapping up
df.X2 %>% filter(gene_name %in% df.X1$gene_name) %>% filter(log2FoldChange < 0) %>% nrow() #overlapping up


#shared SP: brown and white with t implant
df.shared.t.sp <- df.X2 %>% filter(gene_name %in% df.X1$gene_name) %>% 
  filter(blast_name!="unknown") %>% distinct(blast_name, .keep_all = T) %>% 
  mutate(sig = "Y")

#unique background
df.BG.shared.t.sp <- df.in.all %>% filter(!blast_name %in% df.shared.t.sp$blast_name) %>% 
  filter(blast_name!="unknown") %>% distinct(blast_name, .keep_all = T) %>% 
  mutate(sig = "N")

goseq.in.shared.t.sp <- rbind(df.shared.t.sp, df.BG.shared.t.sp)

genes = as.integer(goseq.in.shared.t.sp$sig == "N")

not_na = !is.na(genes)
names(genes) = goseq.in.shared.t.sp$blast_name
genes = genes[not_na]

pwf=nullp(genes,"hg19","geneSymbol")
head(pwf)
GO.wall=goseq(pwf,"hg19","geneSymbol")
head(GO.wall)



GO.wall.sig<-GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]

GO.wall.out <- GO.wall %>% filter(category %in% GO.wall.sig) %>% 
  dplyr::select(category, ontology, numDEInCat, numInCat, term, GOterm = category, description = term)

GO.wall.sig.lg<-length(GO.wall.sig)

GO.wall.sig.lg.godsec<-capture.output(for(go in GO.wall.sig[1:GO.wall.sig.lg]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
})

ft <- flextable(GO.wall.out)
ft <- autofit(ft)
save_as_docx(ft, path = "output/go_terms_2020/SP.white.brown.overlap.testosterone.goterms.docx")


#x<-df.de.t
#contrast.res.func<-x
#
##get rid of any genes that matched to more than one WSFW gene model that werent significant diferentially expressed in this comparison. This is just for background data, so probably fine to do this
##to include "NA" columns as not DEGs, easiest to call them 999
#contrast.res.func$padj[is.na(contrast.res.func$padj)]<-999
##contrast_nonsig_notdup<-dplyr::filter(contrast.res.func,padj>0.05)
#contrast_nonsig_notdup<-df.background
#contrast_nonsig_notdup <- contrast_nonsig_notdup %>% filter(!blast_name %in% contrast.res.func$blast_name) %>% 
#  distinct(blast_name, .keep_all = T)
##next, setup a dataframe for parsing the significant genes
#contrast_sig<-dplyr::filter(contrast.res.func,padj<0.05) %>% filter(blast_name!="unknown")
##collapse any genes that match to two different RNAseq found genes. This doesnt remove that gene, nor does it fix the problem that it identifies with two gene models of WSFW, but it is a bandaid.  
#contrast_sig<-contrast_sig[!duplicated(contrast_sig$blast_name),]
##any gene models that are unknown function are NA right now and need a placeholder name
#contrast_sig <- contrast_sig %>% filter(blast_name!="unknown")
##contrast_sig$blast_name[is.na(contrast_sig$blast_name)]<-"unknown"
##put nonsig and sig values back together
#contrast.res.func2<-rbind(contrast_nonsig_notdup,contrast_sig)
#
##WITH HUMAN
#genes = as.integer(contrast.res.func2$padj < 0.05)
#not_na = !is.na(genes)
#names(genes) = contrast.res.func2$blast_name
#genes = genes[not_na]
#
#pwf=nullp(genes,"hg19","geneSymbol")
#head(pwf)
#GO.wall=goseq(pwf,"hg19","geneSymbol")
#head(GO.wall)
#



# chest ----------------------------------------
#code is not working for unclear reasons

unique(df.in.all$comparison)
df.chest.bw <- df.in.all %>% filter(comparison == "Chest_White_Black_DE_genes" & padj < 0.05)

df.chest.bw %>%  filter(log2FoldChange > 0) %>% nrow() #overlapping up
df.chest.bw %>%  filter(log2FoldChange < 0) %>% nrow() #overlapping up
df.chest.bw %>% filter(melano.list == "Y" & log2FoldChange < 0)

df.chest.bw<-df.chest.bw %>% 
  filter(blast_name!="unknown") %>% distinct(blast_name, .keep_all = T) %>% 
  mutate(sig = "Y")

#unique background
df.BG.chest.bw <- df.in.all %>% filter(comparison == "Chest_White_Black_DE_genes") %>% 
  filter(!blast_name %in% df.chest.bw$blast_name) %>% 
  filter(blast_name!="unknown") %>% distinct(blast_name, .keep_all = T) %>% 
  mutate(sig = "N")

goseq.in.chest.bw <- rbind(df.chest.bw, df.BG.chest.bw)

goseq.in.chest.bw %>% group_by(blast_name) %>% filter(blast_name == "unknown")

genes = as.integer(goseq.in.chest.bw$sig == "N")

not_na = !is.na(genes)
names(genes) = goseq.in.chest.bw$blast_name
genes = genes[not_na]

df.pwf=nullp(genes,"hg19","geneSymbol")
head(df.pwf)
class(df.pwf)
pwf <- df.pwf %>% filter(!is.na(pwf))
GO.wall=goseq(pwf,"hg19","geneSymbol")
head(GO.wall)



GO.wall.sig<-GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]

GO.wall.out <- GO.wall %>% filter(category %in% GO.wall.sig) %>% 
  dplyr::select(category, ontology, numDEInCat, numInCat, term, GOterm = category, description = term)

GO.wall.sig.lg<-length(GO.wall.sig)

GO.wall.sig.lg.godsec<-capture.output(for(go in GO.wall.sig[1:GO.wall.sig.lg]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
})

ft <- flextable(GO.wall.out)
ft <- autofit(ft)
save_as_docx(ft, path = "output/go_terms_2020/SP.white.brown.overlap.testosterone.goterms.docx")

# dump --------------------------------------------------------------------

### Return the IDs for the gene symbols in the DE results
#
#df.t.chicken <- galgal5 %>% filter(symbol %in% df.de.t$symbol)
#df.t.human <- grch37 %>% filter(symbol %in% df.de.t$symbol)
#df.background.chicken <- galgal5 %>% filter(symbol %in% df.background$symbol)
#df.background.human <- grch37 %>% filter(symbol %in% df.background$symbol)
#
#
### The gene names can map to more than one Ensembl ID (some genes change ID over time), 
### so we need to remove duplicate IDs prior to assessing enriched GO terms
#df.t.chicken.nd <- which(duplicated(df.t.chicken$symbol) == FALSE)
#df.t.human.nd <- which(duplicated(df.t.human$symbol) == FALSE)
#df.background.human.nd <- which(duplicated(df.background.human$symbol) == FALSE)
#df.background.chicken.nd <- which(duplicated(df.background.chicken$symbol) == FALSE)
#
#df.t.chicken <- df.t.chicken[df.t.chicken.nd, ] 
#df.t.human <- df.t.human[df.t.human.nd, ] 
#df.background.human <- df.background.human[df.background.human.nd, ]
#df.background.chicken <- df.background.chicken[df.background.chicken.nd, ]
#
### Merge the IDs with the results 
#df.de.t.chicken <- inner_join(df.de.t, df.t.chicken, by=c("symbol"="symbol"))
#df.de.t.human<- inner_join(df.de.t, df.t.human, by=c("symbol"="symbol"))
#df.background.chiken.merge <- inner_join(df.background, df.background.chicken, by=c("symbol"="symbol"))
#df.background.human.merge <- inner_join(df.background, df.background.human, by=c("symbol"="symbol"))



library(data.table)
library(tidyverse)


# load data ---------------------------------------------------------------

#load annotation
suppressWarnings(source("annotation_import.r")) #generic import of annotations for this type of thing I put together
#load sample table for summaries
s2c <- read.table("output/sample_table_RSEM.csv", header = TRUE, stringsAsFactors = TRUE, sep = ",")

#excludes outlier in moretoni
n.comps <- read.table("output/sample_table_RSEM_sample_size.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

#load genes in outlier windows from Fst scans

chest.fst<-dget("data/genome_outliers/blackchest_whitechest_autosome_nocov_NR_bycolor_50000_OL-Win.fst.txt_Genes_in_99.9q.db") %>% 
  rbind(dget("data/genome_outliers/blackchest_whitechest_zchr_nocov_NR_bycolor_50000_OL-Win.fst.txt_Genes_in_99.9q.db")) %>%
  select(Name, scaff, Nsites, fst)

sp.fst<-dget("data/genome_outliers/whitesp_lorentzi_autosome_nocov_NR_bycolor_50000_OL-Win.fst.txt_Genes_in_99.9q.db") %>% 
  rbind(dget("data/genome_outliers/whitesp_lorentzi_zchr_nocov_NR_bycolor_50000_OL-Win.fst.txt_Genes_in_99.9q.db")) %>% 
  select(Name, scaff, Nsites, fst)

#define paths
path <- "output/deseq_DE_output"
files <- dir(path, pattern = '\\DE_genes.txt', full.names = TRUE)
tables <- lapply(files, fread, header=TRUE)

#read in all comparisons
all.comps = tibble(comparison = files) %>%
  mutate(Data = lapply(comparison, fread, header = T)) %>%
  unnest(Data) %>%
  mutate(comparison = gsub(".txt","",comparison)) %>%
  mutate(comparison = gsub("output/deseq_DE_output/","",comparison)) %>% 
  #separate(comparison, c("V1", "type2"), extra = "merge", sep = "_", remove = FALSE) %>% 
  #select(-type2) %>% 
  mutate(type = ifelse(grepl("AllParts", comparison), "AllParts", 
                       ifelse(grepl("aida", comparison) | 
                                grepl("moretoni", comparison) | 
                                grepl("naimii", comparison) | 
                                grepl("lorentzi", comparison), "Population", "Part"))) %>% 
  mutate(melano.list = ifelse(blast_name %in% mel$gene_name, "Y", "N")) %>% 
  mutate(keratin.follicle = ifelse(grepl("KRF", blast_name), "Y", "N")) %>% 
  mutate(chest.outlier = ifelse(gene_name %in% chest.fst$Name, "Y", "N")) %>% 
  mutate(sp.outlier = ifelse(gene_name %in% sp.fst$Name, "Y", "N"))

all.comps %>% distinct(gene_name) %>% nrow()

#NOTE THAT CANNOT DO DIFFERENT P VALUE CUTOFF WITHOUT EDITING run_DESEQ2.R
#note this is FDR not p value

#t implant and sp contrasts are reversed in direction (because After is before Before alphabetically)
write_csv(all.comps, "output/all_comparison_DE_output_fdr_0.05.csv")

# set up male female analysis ---------------------------------------------

unique(all.comps$comparison)
df.subspecies <- all.comps %>% filter(comparison == "moretoni_lorentzi_SP_White_Brown_DE_genes", padj < 0.05 & !is.na(padj))
df.timplant <- all.comps %>% filter(comparison == "lorentzi_lorentzi_SP_White_Brown_DE_genes", padj < 0.05 & !is.na(padj))
df.timplant.overlap.white.brown.sp <- all.comps %>% filter(comparison == "SP_White_Brown_DE_genes" & padj < 0.05 & !is.na(padj) & gene_name %in% df.timplant$gene_name)

all.comps %>% filter(comparison == "lorentzi_lorentzi_SP_White_Brown_DE_genes" & padj > 0.05 | comparison == "lorentzi_lorentzi_SP_White_Brown_DE_genes" & is.na(padj)) %>% write_csv("output/go_terms_2020/all_T_genes_for_background_fdr_0.05.csv")
df.sex <- all.comps %>% filter(comparison == "M_moretoni_vs_F_lorentzi_Before_DE_genes", padj < 0.05 & !is.na(padj))
df.sex.null <- all.comps %>% filter(comparison == "M_moretoni_vs_F_moretoni_DE_genes", padj < 0.05 & !is.na(padj))

#this overlap should take into account expression direction! they should all be

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(df.sex$gene_name,df.timplant$gene_name,  df.subspecies$gene_name),
  category.names = c("Sexual dichromatism", "implant" , "Subspecies color" ),
  filename = 'output/venn_diagram/venn_sex_comps_0.05.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 900 , 
  width = 900 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .8,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

df.subspecies %>% filter(gene_name %in% df.timplant$gene_name & gene_name %in% df.sex$gene_name) %>% 
  write_csv("output/venn_diagram/overlap_29_genes_all_three_cats_fdr_0.05.csv")

#omit those diff between M F more
df.subspecies %>% filter(gene_name %in% df.timplant$gene_name & gene_name %in% df.sex$gene_name & !gene_name %in% df.sex.null$gene_name) %>% 
  write_csv("output/venn_diagram/overlap_genes_all_three_cats_null_omitted_fdr_0.05.csv")

df.subspecies %>% filter(gene_name %in% df.timplant$gene_name) %>% #similar to what was in original ms
  write_csv("output/venn_diagram/overlap_genes_t_treat_subspecies_More_Lore_fdr_0.05.csv")
df.sex.null %>% write_csv("output/venn_diagram/M_more_F_more_shoulder_same_pheno_fdr_0.05.csv")
df.subspecies %>% filter(gene_name %in% df.timplant$gene_name & !gene_name %in% df.sex$gene_name) %>% #similar to what was in original ms
  write_csv("output/venn_diagram/overlap_genes_t_treat_subspecies_More_Lore_exclude_sex_fdr_0.05.csv")
df.timplant.overlap.white.brown.sp %>% write_csv("output/go_terms_2020/SP_white_brown_overlap_t_treatment_fdr_0.05.csv")
  
df.timplant %>% write_csv("output/venn_diagram/t_implant_de_genes_only_fdr_0.05.csv")



# output file for supplemental --------------------------------------------
unique(all.comps$comparison)

all.comps %>% filter(padj < 0.05) %>% 
  filter(!grepl("Crown", comparison)) %>% 
  filter(!grepl("AllParts", comparison)) %>% 
  mutate(idx = paste(comparison, gene_name, sep = "_")) %>% 
  distinct(idx, .keep_all = T) %>% #remove duplicated genes (they got two blast hits)
  select(-idx) %>% 
  select(comparison, type, gene_name, blast_name, melano.list, baseMean, log2FoldChange, pvalue, padj) %>% 
  write_csv("output/all_de_comparisons.csv")

# paper values ------------------------------------------------------------

#df.sex %>% distinct(gene_name) %>% 
df.sex %>% filter(gene_name %in% df.timplant$gene_name) %>% distinct(gene_name) %>% nrow()

#middle value
df.sex %>% filter(gene_name %in% df.timplant$gene_name & gene_name %in% df.subspecies$gene_name) %>% 
  select(gene_name, blast_name) %>% write_csv("output/venn_diagram/14_middle_genes_fdr_0.05.csv")


#how many are female specific?
df.implant.female <- df.timplant %>% distinct(gene_name, .keep_all = T) %>% 
  filter(!gene_name %in% df.sex$gene_name) %>% nrow() 
df.subspecies.female <- df.subspecies %>%  distinct(gene_name, .keep_all = T) %>% 
  filter(!gene_name %in% df.sex$gene_name) %>% nrow()

#how many are female specific?
(1-(196+14+1) / ((df.implant.female+df.subspecies.female) + (196+14+1))) * 100

#male female same color
df.sex.null %>% nrow()


#implant overlap with subspecies
df.implant.female.overlap.ssp <- df.timplant %>% distinct(gene_name, .keep_all = T) %>% 
  filter(gene_name %in% df.subspecies$gene_name) %>% nrow()

#the paper value should be how many genes are shared between T and between subspecies comparisons, not shared with male:
14 / df.implant.female.overlap.ssp

df.timplant %>% distinct(gene_name, .keep_all = T) %>% nrow()

df.timplant %>% distinct(gene_name, .keep_all = T) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(gene_name %in% df.subspecies$gene_name) %>% nrow() 

df.timplant %>% distinct(gene_name, .keep_all = T) %>% 
  filter(log2FoldChange < 0) %>% 
  filter(gene_name %in% df.subspecies$gene_name) %>% nrow() 

# Summary -----------------------------------------------------------------

summary.all.comps <- all.comps %>% filter(padj < 0.05 & !is.na(padj)) %>%
  group_by(type, comparison) %>% 
  summarise(n = n(),
            n.mel = sum(melano.list == "Y"),
            n.ker = sum(keratin.follicle == "Y"),
            chest.outlier = sum(chest.outlier == "Y"),
            sp.outlier = sum(sp.outlier == "Y"),
            mel.chest.outlier = sum(chest.outlier == "Y" & melano.list == "Y"),
            mel.sp.outlier = sum(sp.outlier == "Y" & melano.list == "Y"),
            ker.sp.outlier = sum(sp.outlier == "Y" & keratin.follicle == "Y"),
            ker.chest.outlier = sum(chest.outlier == "Y" & keratin.follicle == "Y")) %>% 
  left_join(n.comps, by = c("type", "comparison"))
write_csv(summary.all.comps, "output/de_summary_table_0.05.csv")

  #how many up and down reg
dr <- all.comps %>% filter(padj < 0.1 & !is.na(padj)) %>%
    group_by(type, comparison) %>% 
  filter(log2FoldChange < 0) %>% 
  summarize(down_regulated = n ())

ur <- all.comps %>% filter(padj < 0.1 & !is.na(padj)) %>%
  group_by(type, comparison) %>% 
  filter(log2FoldChange > 0) %>% 
  summarize(up_regulated = n ())


dr.ur <- left_join(dr, ur, by = c("type", "comparison")) 

summary.all.comps.full <- left_join(summary.all.comps, dr.ur, by = c("type", "comparison"))
write_csv(summary.all.comps.full, "output/full_de_summary_table_0.05.csv")

all_de_genes <- all.comps %>% filter(padj < 0.1 & !is.na(padj)) %>%
  dplyr::select(comparison, gene_name, padj)
write_csv(all_de_genes, "output/all_de_genes_0.05.csv")

all_de_genes2 <- all.comps %>% filter(padj < 0.1 & !is.na(padj))
  
timplat.de <- all_de_genes2 %>% filter(comparison == "lorentzi_lorentzi_SP_White_Brown_DE_genes")
sp.white.brown.de <- all_de_genes2 %>% filter(comparison == "SP_White_Brown_DE_genes")

#how many are also in between pop comparisons?
timplat.de %>% filter(timplat.de$gene_name %in% sp.white.brown.de$gene_name & log2FoldChange > 0) %>% nrow()
timplat.de %>% filter(timplat.de$gene_name %in% sp.white.brown.de$gene_name & log2FoldChange < 0) %>% nrow()

timplat.overlap.dr <- timplat.de %>% filter(timplat.de$gene_name %in% sp.white.brown.de$gene_name & log2FoldChange < 0)



# playground --------------------------------------------------------------

cwb <- all.comps %>% filter(padj < 0.1 & !is.na(padj)) %>% filter(comparison == "Chest_White_Black_DE_genes")
#(HPGDS, FRZB, MLANA, PMEL, SLC24A4, TYR) all these are still DEed in this

swb <- all.comps %>% filter(padj < 0.1 & !is.na(padj)) %>% filter(comparison == "SP_White_Brown_DE_genes")

all.comps %>% filter(padj < 0.1 & !is.na(padj)) %>% 
  filter(comparison == "lorentzi_lorentzi_SP_White_Brown_DE_genes" & sp.outlier == "Y")

all.comps %>% filter(padj < 0.1 & !is.na(padj)) %>% 
  filter(comparison == "Chest_White_Black_DE_genes" | comparison == "SP_White_Brown_DE_genes" | comparison == "lorentzi_lorentzi_SP_White_Brown_DE_genes") %>% 
  filter(chest.outlier == "Y" | sp.outlier == "Y") #%>% View()

all.comps %>% filter(gene_name == "WSFW009166")  %>% select(comparison, log2FoldChange, padj) %>% arrange(log2FoldChange)
  ggplot() + geom_point(aes(x = comparison, y = log2FoldChange))
  
chest.white.black <- all.comps %>% filter(comparison == "Chest_White_Black_DE_genes")
chest.white.black.sig <- all.comps %>% filter(comparison == "Chest_White_Black_DE_genes" & padj < 0.1)

ggplot() + 
  geom_point(data = chest.white.black, aes(y = log2FoldChange, x = baseMean)) +
  geom_point(data = chest.white.black.sig, aes(y = log2FoldChange, x = baseMean), color = "red") 
  xlim(0,1000) + ylim(-2,2)
  

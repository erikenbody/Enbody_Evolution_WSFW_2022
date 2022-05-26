library(tidyverse)
library(data.table)
df.genos <- fread("data/genotypes/WSFW_scaffold_142.1.2.1.5.pos.genos") %>% select(-V42) %>% 
  mutate(full.pos = paste(CHROM, POS, sep = "_"))
df.rbfw <- fread("data/genotypes/RBFW_sc142_inv.pos.genos") %>% select(-V6, -REF, -ALT, A_1_M.melanocephalus_A = RBFW_220bp_GCTCAT) %>% 
  mutate(full.pos = paste(CHROM, POS, sep = "_")) %>% select(-CHROM, -POS)

df.genos <- df.genos %>% left_join(df.rbfw, by = "full.pos")
df.genos <- df.genos %>% select(-full.pos)

df.outliers <- read.csv("data/persite_fst/kitlg_outlier_pos.csv")

df.genos <- df.genos %>% filter(POS %in% df.outliers$V2)

nrow(df.outliers) #210
#remove tri allelic
df.genos <- df.genos %>% filter(!grepl(",", ALT))

long.df.genos <- df.genos %>% 
  pivot_longer(cols = -c("CHROM", "POS", "ALT", "REF"), names_to = "sample", values_to = "genos") %>% 
  mutate(genos = ifelse(genos == "./." | is.na(genos), "NA", genos))



df.pops <- tibble(pop = c("aida","moretoni", "naimii", "lorentzi", "M.melanocephalu"),
       pop.num = c(1,2,3,4,5))

long.df.genos <- long.df.genos %>% separate(sample, into = c(NA,"band","pop",NA), remove = F, sep = "_") %>% 
  left_join(df.pops, by = "pop") %>% 
  mutate(label = paste(pop, band, sep = "-")) %>% 
  mutate(Genotype = as.factor(genos)) %>% 
  mutate(Genotype = ifelse(Genotype %in% c("1/2","0/2","0/4","0/3","2/2","2/3"), NA, as.character(Genotype))) #only applies to RBFW but removes triallelic sites

#colors <- c('#cd853f','black',"#3f87cd", "white")
#colors <- c('#e9a3c9','black','#ffffbf','white')
colors <- c('#ef8a62','black','#67a9cf', "white")

long.df.genos %>% ggplot(aes(x = as.factor(POS), y = fct_reorder(label, pop.num))) + 
  geom_tile(aes(fill = Genotype), color = "white") +
  scale_fill_manual(values=colors)+
  theme(axis.text.x = element_text(angle = 90, size = 3),
        axis.text.y = element_text(size = 20)) + 
  xlab(element_blank()) + ylab(element_blank())

ggsave("output/haplotype_heatmap/kitlg_region_heatmap.png", height = 10, width = 10)

#kitlg range
df.outliers %>% filter(V2 > 1352258 & V2 < 1380411) %>% select(V2) %>% max()
df.outliers %>% filter(V2 > 1352258 & V2 < 1380411) %>% select(V2) %>% min()
#TMTC3
df.outliers %>% filter(V2 > 1239858 & V2 < 1269494) %>% select(V2) %>% max()
df.outliers %>% filter(V2 > 1239858 & V2 < 1269494) %>% select(V2) %>% min()

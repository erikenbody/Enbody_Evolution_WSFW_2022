library(tidyverse)
library(data.table)
library(ggrepel)
bc_persite <- fread("data/persite_fst/nocov_NR_bycolor/Chest_persite_fst_scaffold_142.txt")
bc_persite$fst <- bc_persite$V3 / bc_persite$V4

Cgene_scaf<-142
Cgene_start<-1352258
Cgene_end<-1380411
Cgene_name <- "KITLG"

gene.subset <- gene.ranges %>% filter(seqid == "scaffold_142" & start > 1.1e6 & end < 1.5e6) %>% 
  mutate(start.kb = start/1000, end.kb = end/1000)

#pdf("output/persite_fst-KITLG.pdf",height=4,width=12)
png("output/persite_fst/persite_fst-KITLG.png", width = 1200, height = 900)
#pdf("output/persite_fst/persite_fst-KITLG.pdf", width = 1200, height = 900)
bc_persite %>% filter(V2 > 1.1e6 & V2 < 1.5e6) %>% 
  ggplot() + geom_point(aes(x = V2/1000, y = fst), size = 3) + 
  theme_classic() + 
  theme(legend.position = "none", text = element_text(size=50)) +
  #scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
  geom_segment(data = gene.subset, aes(x = start.kb, y = 0.95, xend = end.kb, yend = 0.95), size = 2) +
  geom_label_repel(data = gene.subset, aes(x = start.kb, y = 1, label = blast_name), size = 12) +
  #xlab(paste0("scaffold_", Cgene_scaf,"/", "kb")) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = seq(1000, 1500, by = 100))
dev.off()

bc_persite %>% #filter(V2 > 0.0e6 & V2 < 2.5e6) %>% 
  filter(fst > 0.3) %>% 
  ggplot() + geom_point(aes(x = V2/1000, y = fst)) + 
  theme_classic() + theme(legend.position = "none", text = element_text(size=20)) +
  #scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
  geom_segment(data = gene.subset, aes(x = start.kb, y = 0.95, xend = end.kb, yend = 0.95)) +
  geom_label_repel(data = gene.subset, aes(x = start.kb, y = 1, label = blast_name), size = 8) +
  xlab(Cgene_scaf)  
  

df.outliers <- bc_persite %>% filter(fst > 0.7 & V2> 1.2e6 & V2 < 1.5e6)
write.csv(df.outliers, "data/persite_fst/kitlg_outlier_pos.csv")
# aida naimii -------------------------------------------------------------

a_n_persite <- fread("data/persite_fst/persite_fst_aida_naimii_scaffold_142.txt")
a_n_persite$fst <- a_n_persite$V3 / a_n_persite$V4
 

#pdf("output/persite_fst-KITLG.pdf",height=4,width=12)
png("output/persite_fst_aida_naimii-KITLG.png")
a_n_persite %>% filter(V2 > 0.5e6 & V2 < 1.8e6) %>% 
  ggplot() + geom_point(aes(x = V2, y = fst), pch='.') + 
  theme_classic() + theme(legend.position = "none", text = element_text(size=20)) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
  geom_segment(data = gene.subset, aes(x = start, y = 0.95, xend = end, yend = 0.95)) +
  geom_label_repel(data = gene.subset, aes(x = start, y = 1.25, label = blast_name)) +
  xlab(Cgene_scaf)
dev.off()


# -------------------------------------------------------------------------
#include sweepfinder


df.sw1 <- read.table("data/sweepfinder/moretoni_sf2_200_wSFS.txt", header = T) %>% 
  mutate(pop = "moretoni")
df.sw2 <- read.table("data/sweepfinder/aida_sf2_200_wSFS.txt", header = T) %>% 
  mutate(pop = "aida")
df.sw3 <- read.table("data/sweepfinder/naimii_sf2_200_wSFS.txt", header = T) %>% 
  mutate(pop = "naimii")
df.sw4 <- read.table("data/sweepfinder/lorentzi_sf2_200_wSFS.txt", header = T) %>% 
  mutate(pop = "lorentzi")

df.sw <- rbind(df.sw1, df.sw2, df.sw3, df.sw4)

pfst <- bc_persite %>% 
  filter(V2 > 0.7e6 & V2 < 1.7e6) %>% 
  ggplot() + 
    geom_point(aes(x = V2, y = fst)) + 
    theme_bw() + 
  theme(legend.position = "none", text = element_text(size=20)) +
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
    geom_segment(data = gene.subset, aes(x = start, y = 0.95, xend = end, yend = 0.95)) +
    geom_label_repel(data = gene.subset, aes(x = start, y = 1.25, label = blast_name)) +
    xlab(NULL) + ylab("Fst") 

df.sw <- df.sw %>% dplyr::rename(subspecies = pop)

p.clr <- ggplot() + 
  geom_line(data = df.sw, aes(x = location, y = LR, color = subspecies)) +
  theme_bw() + 
  theme(legend.position = c(0.8, 0.8), text = element_text(size=20)) +
  labs(y = "CLR", x = "Scaffold 142 position (bp)") +
  lims(x = c(0.7e6, 1.7e6))
  

png("output/persite_fst_aida_naimii-KITLG_with_CLR.png", width = 1200, height = 800)
pfst / p.clr
dev.off()

ap<-ggplot() + 
  geom_point(data = df.sw, aes(x = location, y = LR, color = subspecies)) +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size=20)) +
  labs(y = "CLR", x = "Scaffold 142 position (bp)") +
  lims(x = c(1.1e6, 1.3e6))

bp <- bc_persite %>% 
  #filter(V2 > 1.1e6 & V2 < 1.3e6) %>% 
  ggplot() + 
  geom_point(aes(x = V2, y = fst)) + 
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size=20)) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
  #geom_segment(data = gene.subset, aes(x = start, y = 0.95, xend = end, yend = 0.95)) +
  #geom_label_repel(data = gene.subset, aes(x = start, y = 1.25, label = blast_name)) +
  xlab(NULL) + ylab("Fst") 
ap/bp


df.sw2 <- read.table("data/sweepfinder/naimii_autosome_scaffolds_scaffold_142_sw2_5kb.txt", header = T)
head(df.sw2)
b3 <- df.sw2 %>% ggplot() + geom_point(aes(x = location, y = LR))
b3/bp

library(tidyverse)
library(patchwork)
library(RcppCNPy)
library(ggthemes)
library(forcats)

location.order <- c("lorentzi","moretoni","naimii","aida")

# NGSadmix method ---------------------------------------------------------
#files say super 1 but are whole genome

S.admix <- read.table("data/ngsadmix/ngs_admix_k4.qopt") # Reads results from selection scan
pca.names <- read.table("data/ngsadmix/ngsadmix_names.txt", stringsAsFactors = F) #custom formatted 

df.admix.k4 <- as.data.frame(S.admix)
df.admix.k4$sampleID <- pca.names$V1
df.admix.k4$loc <- pca.names$V2

df.admix.k4 <- df.admix.k4 %>% pivot_longer(names_to = "popGroup", values_to = "prob", cols = -c(sampleID, loc)) %>% 
  mutate(popGroup = gsub("V", "", popGroup))

df.admix.k4$Location <- factor(df.admix.k4$loc, levels = location.order)

#https://luisdva.github.io/rstats/model-cluster-plots/

p.ngsadmix.k4 <- ggplot(df.admix.k4, aes(x = factor(sampleID), y = prob, fill = factor(popGroup))) +
  geom_col(aes(color = popGroup), size = 0.1) +
  #facet_grid(~loc, switch = "x", scales = "free", space = "free") +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=4", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE) + xlab(NULL) +
  scale_fill_manual(values = c("black", "#CD853F","grey80",  "grey40"), guide = F) +
  scale_color_manual(values = c("black", "#CD853F","grey80", "grey40"), guide = F)

p.ngsadmix.k4

ggsave("output/ngsadmix/wsfw_k4_ngsadmix.pdf", width = 10, height = 4)



# with pcangsd ------------------------------------------------------------
pop<-read.table("popfile.txt")

adm <- npyLoad("data/pcangsd/WSFW_pcangsd_admix.admix.Q.npy") 
df.pcangsd.adx <- as.data.frame(adm)
df.pcangsd.adx$sampleID <- pop$V1
df.pcangsd.adx$loc <- pop$V2


df.pcangsd.adx <- df.pcangsd.adx %>% pivot_longer(names_to = "popGroup", values_to = "prob", cols = -c(sampleID, loc)) %>% 
  mutate(popGroup = gsub("V", "", popGroup))

df.pcangsd.adx$Location <- factor(df.pcangsd.adx$loc, levels = location.order)

#https://luisdva.github.io/rstats/model-cluster-plots/

p.pcangsd.k4 <- ggplot(df.pcangsd.adx, aes(x = factor(sampleID), y = prob, fill = factor(popGroup))) +
  geom_col(aes(color = popGroup), size = 0.1) +
  #facet_grid(~loc, switch = "x", scales = "free", space = "free") +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=4", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE) + xlab(NULL) +
  scale_fill_manual(values = c("#CD853F", "grey40", "black", "grey80"), guide = F) +
  scale_color_manual(values =  c("#CD853F", "grey40", "black", "grey80"), guide = F)

p.pcangsd.k4
ggsave("output/ngsadmix/wsfw_k4_pcangsd.pdf", width = 10, height = 4)

library(patchwork)
p.ngsadmix.k4 / p.pcangsd.k4

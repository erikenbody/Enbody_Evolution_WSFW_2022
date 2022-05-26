library(tidyverse)
library(patchwork)

pca.wsfw <- read.table("data/pcangsd/WSFW_2.cov")
pca.names<-read.table("popfile.txt")
names(pca.names) <- c("sample","pop")

#pca.wsfw$pop <- pca.names$pop

#pca.names <- pca.names %>% left_join(pop.names, by = "pop")

eigen.wsfw <- eigen(pca.wsfw, symmetric = T)
eigenvectors.wsfw <- as.data.frame(eigen.wsfw$vectors)
eigenvalues.wsfw <- eigen.wsfw$values

eigenvectors.wsfw$sample <- pca.names$sample

eigenvectors.wsfw$pop <- pca.names$pop

p1 <- ggplot() + 
  geom_point(data = eigenvectors.wsfw, aes(x = V2, y = V1, fill = pop), color = "black", size = 6,shape = 21) +
  ylab(paste0("PC1 ", "(", round(eigenvalues.wsfw[1],1),"%)")) +
  xlab(paste0("PC2 ", "(", round(eigenvalues.wsfw[2], 1),"%)")) +
  theme_bw() +
  theme(text = element_text(size = 24))+
  scale_fill_manual(values = c("black", "#CD853F","black","white"), guide = F)

p1

#pdf("output/pcangsd/wsfw_pca_plot.pdf", width = 1800, height = 1800)
#p1
#dev.off()


ggsave("output/pcangsd/wsfw_pca_plot.pdf", width = 10, height = 10)


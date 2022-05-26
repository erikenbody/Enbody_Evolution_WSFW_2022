library(tidyverse)

de.fd.1 <- read.csv("output/all_comparison_DE_output.csv")
de.fd.05 <- read.csv("output/all_comparison_DE_output_fdr_0.05.csv")

de.fd.1.s <- de.fd.1 %>% filter(padj<0.1)
de.fd.05.s <- de.fd.05 %>% filter(padj<0.05)
nrow(de.fd.1.s)
nrow(de.fd.05.s)

de.fd.1.s.f <- de.fd.1.s %>% select("comparison","gene_name","blast_name","baseMean","log2FoldChange", "padj")
de.fd.05.s.f <- de.fd.05.s %>% select("comparison","gene_name","blast_name","baseMean","log2FoldChange", "padj")

de.compare.thresh <- full_join(de.fd.1.s.f, de.fd.05.s.f, by = c("comparison","gene_name","blast_name"))

de.sum.1 <- read.csv("output/full_de_summary_table.csv")
de.sum.05 <- read.csv("output/full_de_summary_table_0.05.csv")

comp.sum <- left_join(de.sum.1[,c(1,2,3)], de.sum.05[,c(1,2,3)], by = c("type","comparison"))

ggplot(comp.sum) + geom_point(aes(x = n.x, y = n.y)) +
  labs(x = "P < 0.1", y = "P < 0.05") + xlim(0,2000) + ylim(0,2000)

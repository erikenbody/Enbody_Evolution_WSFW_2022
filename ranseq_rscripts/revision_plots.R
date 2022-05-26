library(tidyverse)
library(gghighlight)
library(ggrepel)

de.fd.1 <- read.csv("output/all_comparison_DE_output.csv") %>% 
  mutate(padj.plot = ifelse(log2FoldChange < 0, log10(padj) * -1, log10(padj))) %>% 
  mutate(sig = ifelse(padj < 0.1, "Y", "N"))
de.fd.05 <- read.csv("output/all_comparison_DE_output_fdr_0.05.csv") %>% 
  mutate(padj.plot = ifelse(log2FoldChange < 0, log10(padj) * -1, log10(padj))) %>% 
  mutate(sig = ifelse(padj < 0.05, "Y", "N"))

unique(de.fd.05$comparison)

df.SP.white.brown <- de.fd.05 %>% filter(comparison == "SP_White_Brown_DE_genes") %>% 
  select("gene_name","blast_name","padj.plot", "sig") 
df.SP.more.lore <- de.fd.05 %>% filter(comparison == "moretoni_lorentzi_SP_White_Brown_DE_genes") %>% 
  select("gene_name","blast_name","padj.plot", "sig") 
df.t.treat <- de.fd.05 %>% filter(comparison == "lorentzi_lorentzi_SP_White_Brown_DE_genes") %>% 
  select("gene_name","blast_name","padj.plot", "sig") 


# fig 3 -------------------------------------------------------------------

df.SP_t <- left_join(df.SP.white.brown, df.t.treat, by = c("gene_name","blast_name")) %>% 
  mutate(sig =  ifelse(sig.x == "Y" & sig.y == "Y", "Y", "N"))
df.more.lore_t <- left_join(df.SP.more.lore, df.t.treat, by = c("gene_name","blast_name")) %>% 
  mutate(sig =  ifelse(sig.x == "Y" & sig.y == "Y", "Y", "N"))
  
 ggplot(df.SP_t, aes(x = -padj.plot.x, y = -padj.plot.y))+ geom_point(aes(color = sig), size = 1) +
  theme_bw() + scale_color_manual(values = c("black",  "red")) +
   geom_hline(yintercept=0) +
   geom_vline(xintercept=0) +
  geom_text_repel(data = subset(df.SP_t, -padj.plot.x > 0 & sig=="Y"), 
                  aes(label = blast_name), color = "red") +
  geom_text_repel(data = subset(df.SP_t, -padj.plot.x < -3.2 & -padj.plot.y < -2.2  & sig=="Y"), 
                  aes(label = blast_name), color = "red") +
   theme(legend.position="none",
         axis.text = element_text(size = 14),
         axis.title = element_blank()) +
   xlim(-20,20) +
   ylim(-16,16)
 ggsave("output/revision_plots/SP_white_brown_against_t_treat.png",width = 9, height = 9)

  
 ggplot(df.more.lore_t, aes(x = -padj.plot.x, y = -padj.plot.y))+ geom_point(aes(color = sig), size = 1) +
   theme_bw() + scale_color_manual(values = c("black",  "red")) +
   geom_hline(yintercept=0) +
   geom_vline(xintercept=0) +
   geom_text_repel(data = subset(df.more.lore_t, -padj.plot.x > 0 & sig=="Y"), 
                   aes(label = blast_name), color = "red") +
   geom_text_repel(data = subset(df.more.lore_t, -padj.plot.x < -3 & -padj.plot.y < -2  & sig=="Y"), 
                   aes(label = blast_name), color = "red") +
   theme(legend.position="none",
         axis.text = element_text(size = 14),
         axis.title = element_blank()) +
   xlim(-20,20) +
   ylim(-12,12)
 ggsave("output/revision_plots/SP_more_lore_against_t_treat.png",width = 9, height = 9)


# fig S3 ------------------------------------------------------------------
 df.Chest_White_Black <- de.fd.05 %>% filter(comparison == "Chest_White_Black_DE_genes") %>% 
   dplyr::select("gene_name","blast_name","padj","padj.plot","log2FoldChange","padj.plot", "sig", "melano.list") #%>% 
   #mutate(melano.list = ifelse(melano.list == "Y","melano","not")) %>% 
   #mutate(melano.list = factor(melano.list, levels = c("melano","not")))
 #df.Chest_White_Black$melano.list
 
 ggplot(df.Chest_White_Black, aes(x = -log2FoldChange, y = -log10(padj))) + 
   geom_point(data = df.Chest_White_Black, aes(color = sig), alpha = 0.8, size = 2) +
   geom_text_repel(data = subset(df.Chest_White_Black, -log10(padj) > 2 & sig=="Y"), 
                   aes(label = blast_name, color = melano.list), alpha = 0.8) +
   #geom_text_repel(data = subset(df.Chest_White_Black, -log10(padj) > 2 & sig=="Y" & melano.list == "Y"), 
   #                aes(label = blast_name), color = "black") +
   scale_color_manual(values = c("black",  "red")) +
   theme_bw() +
   theme(legend.position="none",
         axis.text.x = element_text(size = 14),
         axis.text.y = element_text(size = 14),
         axis.title = element_blank()) 
 ggsave("output/revision_plots/chest_volcano_plot.png",width = 9, height = 8)
 
 
  
 
 
library(tidyverse)

df.kitlg <- read.csv("output/kitlg_counts_edited.csv") %>% 
  filter(!grepl("Crown", PopPart))

df.kitlg$PopPart <- as.factor(df.kitlg$PopPart)

df.kitlg$Color.numeric <- as.numeric(df.kitlg$Feather.Color)

df.sum <- df.kitlg %>% group_by(PopPart) %>% 
  summarise(mean_count = mean(count),
            min_count = min(count))

df.kitlg <- df.kitlg %>% left_join(df.sum, by = "PopPart")

ggplot(df.kitlg) + 
  geom_point(aes(x = fct_reorder(PopPart, mean_count), y = count, fill = Feather.Color), size = 5, pch = 21) +
  scale_fill_manual(values = c("black", "#CD853F", "white")) +
  theme_bw() + theme(axis.title.x = element_blank(),
                     legend.position = c(0.1, 0.8)) +
  ylab("Normalized counts")

ggsave("output/kitlg_counts_nice.pdf", width = 13, height = 6, useDingbats=FALSE)


df.kitlg <- df.kitlg %>% separate(PopPart, into = c("pop","part"), sep = "_", remove = F)

ggplot(df.kitlg) + 
  geom_point(aes(x = fct_reorder(PopPart, part), y = count, fill = Feather.Color), size = 5, pch = 21) +
  scale_fill_manual(values = c("black", "#CD853F", "white")) +
  theme_bw() + theme(axis.title.x = element_blank(),
                     legend.position = c(0.1, 0.8)) +
  ylab("Normalized counts")

ggsave("output/kitlg_counts_nice_V2.pdf", width = 13, height = 6, useDingbats=FALSE)


# ESD ---------------------------------------------------------------------


library(tidyverse)

df.esd <- read.csv("output/esd_counts_edited.csv") %>% 
  filter(!grepl("Crown", PopPart))

df.esd$PopPart <- as.factor(df.esd$PopPart)

df.esd$Color.numeric <- as.numeric(df.esd$Feather.Color)

df.sum <- df.esd %>% group_by(PopPart) %>% 
  summarise(mean_count = mean(count),
            min_count = min(count))

df.esd <- df.esd %>% left_join(df.sum, by = "PopPart")

ggplot(df.esd) + 
  geom_point(aes(x = fct_reorder(PopPart, mean_count), y = count, fill = Feather.Color), size = 5, pch = 21) +
  scale_fill_manual(values = c("black", "#CD853F", "white")) +
  theme_bw() + theme(axis.title.x = element_blank(),
                     legend.position = c(0.1, 0.8)) +
  ylab("Normalized counts")

ggsave("output/esd_counts_nice.pdf", width = 13, height = 6, useDingbats=FALSE)

df.esd <- df.esd %>% separate(PopPart, into = c("pop","part"), sep = "_", remove = F)


ggplot(df.esd) + 
  geom_point(aes(x = fct_reorder(PopPart, part), y = count, fill = Feather.Color), size = 5, pch = 21) +
  scale_fill_manual(values = c("black", "#CD853F", "white")) +
  theme_bw() + theme(axis.title.x = element_blank(),
                     legend.position = c(0.1, 0.8)) +
  ylab("Normalized counts")

ggsave("output/esd_counts_nice_V2.pdf", width = 13, height = 6, useDingbats=FALSE)

# EGR1 ---------------------------------------------------------------------

df.EGR1 <- read.csv("output/egr1_counts.csv")

df.EGR1$PopPart <- as.factor(df.EGR1$PopPart)
df.EGR1$Feather.Color <- df.esd$Feather.Color
df.EGR1$Color.numeric <- as.numeric(df.EGR1$Feather.Color)

df.sum <- df.EGR1 %>% group_by(PopPart) %>% 
  summarise(mean_count = mean(count),
            min_count = min(count))

df.EGR1 <- df.EGR1 %>% left_join(df.sum, by = "PopPart")

ggplot(df.EGR1) + 
  geom_point(aes(x = fct_reorder(PopPart, mean_count), y = count, fill = Feather.Color), size = 5, pch = 21) +
  scale_fill_manual(values = c("black", "#CD853F", "white")) +
  theme_bw() + theme(axis.title.x = element_blank(),
                     legend.position = c(0.1, 0.8)) +
  ylab("Normalized counts")

ggsave("output/EGR1_counts_nice.pdf", width = 13, height = 6, useDingbats=FALSE)

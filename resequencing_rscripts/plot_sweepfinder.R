library(tidyverse)
library(patchwork)

df.mafs <- read.table("data/sweepfinder/aida.mafs.gz", header = T)
head(df.mafs)

df.sf <- df.mafs %>% 
  mutate(x = round(unknownEM * nInd,0)) %>% 
  mutate(folded = ifelse(anc == "N", 1, 0)) %>% 
  filter(!(x==0 & folded == 0)) %>% #remove polarized monomorphic sites
  select(position, x, nInd, folded, n = nInd)
write_tsv(df.sf, "data/sweepfinder/aida.sf")
  
hist(df.sf$x/df.sf$n)



# -------------------------------------------------------------------------

df.sw <- read.table("data/sweepfinder/aida_sw", header = T)
ggplot(df.sw) + geom_point(aes(x = location, y = LR)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411)


df.sw <- read.table("data/sweepfinder/moretoni_sf2.txt", header = T)
ggplot(df.sw) + geom_point(aes(x = location, y = LR)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411)

# -------------------------------------------------------------------------

df.aida <- read.table("data/sweepfinder/aida.sf", header = T)
df.naim <- read.table("data/sweepfinder/naimii_sf_input.txt", header = T)
df.more <- read.table("data/sweepfinder/moretoni_sf_input.txt", header = T)
df.lore <- read.table("data/sweepfinder/lorentzi_sf_input.txt", header = T)

df.outliers <- read.csv("data/persite_fst/kitlg_outlier_pos.csv")

df.ai.na <- left_join(df.aida, df.naim, by = "position")
df.ai.na %>% 
  mutate(af.a = (x.x/n.x),
         af.b = (x.y/n.y),
         daf = abs(af.a - af.b)) %>% 
  ggplot() + geom_point(aes(x = position, y = daf)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411)


df.mo.na <- left_join(df.more, df.naim, by = "position")
df.mo.na %>% 
  mutate(af.a = (x.x/n.x),
         af.b = (x.y/n.y),
         daf = abs(af.a - af.b)) %>% 
  ggplot() + geom_point(aes(x = position, y = daf)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411)

df.mo.lo <- left_join(df.more, df.lore, by = "position")
df.mo.lo %>% 
  mutate(af.a = (x.x/n.x),
         af.b = (x.y/n.y),
         daf = abs(af.a - af.b)) %>% 
  ggplot() + geom_point(aes(x = position, y = daf)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411)

df.mo.ai <- left_join(df.more, df.aida, by = "position")
df.mo.ai %>% 
  mutate(af.a = (x.x/n.x),
         af.b = (x.y/n.y),
         daf = abs(af.a - af.b)) %>% 
  ggplot() + geom_point(aes(x = position, y = daf)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411)

# -------------------------------------------------------------------------


df.mo.lo %>% 
  filter(position %in% df.outliers$V2)
  mutate(af.a = (x.x/n.x),
         af.b = (x.y/n.y),
         N = (af.a - af.b)^2,
         D = (af.a + af.b - 2*af.a*af.b),
         FST = N/D) %>% 
  ggplot() + geom_point(aes(x = position, y = FST)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411)


# -------------------------------------------------------------------------

df.mafs.M <- read.table("data/sweepfinder/moretoni.mafs.gz", header = T)
df.mafs.L <- read.table("data/sweepfinder/lorentzi.mafs.gz", header = T)
df.mafs.N <- read.table("data/sweepfinder/naimii.mafs.gz", header = T)
df.mafs.A <- read.table("data/sweepfinder/aida.mafs.gz", header = T)

df.maf.j <- inner_join(df.mafs.M, df.mafs.L, by = "position")

df.maf.j <- df.maf.j %>% 
  #filter(unknownEM.x > 0.1 | unknownEM.y > 0.1) %>% 
  mutate(N = (unknownEM.x  - unknownEM.y)^2,
         D = (unknownEM.x + unknownEM.y - 2*unknownEM.x*unknownEM.y),
         FST = N/D)
ggplot(df.maf.j) + geom_point(aes(x = position, y = FST)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411)

df.maf.AN <- inner_join(df.mafs.A, df.mafs.N, by = "position")

df.maf.AN <- df.maf.AN %>% 
  #filter(unknownEM.x > 0.1 | unknownEM.y > 0.1) %>% 
  mutate(N = (unknownEM.x  - unknownEM.y)^2,
         D = (unknownEM.x + unknownEM.y - 2*unknownEM.x*unknownEM.y),
         FST = N/D)
ggplot(df.maf.AN) + geom_point(aes(x = position, y = FST)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411)


# -------------------------------------------------------------------------


df.sw1 <- read.table("data/sweepfinder/moretoni_sf2_200.txt", header = T) %>% 
  mutate(pop = "moretoni")
df.sw2 <- read.table("data/sweepfinder/aida_sf2_200.txt", header = T) %>% 
  mutate(pop = "aida")
df.sw3 <- read.table("data/sweepfinder/naimii_sf2_200.txt", header = T) %>% 
  mutate(pop = "naimii")
df.sw4 <- read.table("data/sweepfinder/lorentzi_sf2_200.txt", header = T) %>% 
  mutate(pop = "lorentzi")

df.sw <- rbind(df.sw1, df.sw2, df.sw3, df.sw4)
min(df.outliers$V2)
ggplot(df.sw) + geom_point(aes(x = location, y = LR, color = pop)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411) +
  geom_vline(xintercept = min(df.outliers$V2), color = "red") +
  geom_vline(xintercept = max(df.outliers$V2), color = "red") +
  xlim(500000, 2000000)


df.sw1 <- read.table("data/sweepfinder/moretoni_sf2_5k.txt", header = T) %>% 
  mutate(pop = "moretoni")
df.sw2 <- read.table("data/sweepfinder/aida_sf2_5k.txt", header = T) %>% 
  mutate(pop = "aida")
df.sw3 <- read.table("data/sweepfinder/naimii_sf2_5k.txt", header = T) %>% 
  mutate(pop = "naimii")
df.sw4 <- read.table("data/sweepfinder/lorentzi_sf2_5k.txt", header = T) %>% 
  mutate(pop = "lorentzi")

df.sw <- rbind(df.sw1, df.sw2, df.sw3, df.sw4)
min(df.outliers$V2)
ggplot(df.sw) + geom_point(aes(x = location, y = LR, color = pop)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411) +
  geom_vline(xintercept = min(df.outliers$V2), color = "red") +
  geom_vline(xintercept = max(df.outliers$V2), color = "red") +
  xlim(500000, 2000000)



df.sw1 <- read.table("data/sweepfinder/moretoni_sf2_1000.txt", header = T) %>% 
  mutate(pop = "moretoni")
df.sw2 <- read.table("data/sweepfinder/aida_sf2_1000.txt", header = T) %>% 
  mutate(pop = "aida")
df.sw3 <- read.table("data/sweepfinder/naimii_sf2_1000.txt", header = T) %>% 
  mutate(pop = "naimii")
df.sw4 <- read.table("data/sweepfinder/lorentzi_sf2_1000.txt", header = T) %>% 
  mutate(pop = "lorentzi")

df.sw <- rbind(df.sw1, df.sw2, df.sw3, df.sw4)
min(df.outliers$V2)
ggplot(df.sw) + geom_point(aes(x = location, y = LR, color = pop)) +
  geom_vline(xintercept = 1352258) +
  geom_vline(xintercept = 1380411) +
  geom_vline(xintercept = min(df.outliers$V2), color = "red") +
  geom_vline(xintercept = max(df.outliers$V2), color = "red") +
  xlim(500000, 2000000)


df.sw1 <- read.table("data/sweepfinder/moretoni_sf2_200_wSFS.txt", header = T) %>% 
  mutate(pop = "moretoni")
df.sw2 <- read.table("data/sweepfinder/aida_sf2_200_wSFS.txt", header = T) %>% 
  mutate(pop = "aida")
df.sw3 <- read.table("data/sweepfinder/naimii_sf2_200_wSFS.txt", header = T) %>% 
  mutate(pop = "naimii")
df.sw4 <- read.table("data/sweepfinder/lorentzi_sf2_200_wSFS.txt", header = T) %>% 
  mutate(pop = "lorentzi")

df.sw <- rbind(df.sw1, df.sw2, df.sw3, df.sw4)

bc_persite.f <- bc_persite %>% filter(fst > 0.2)

pfst <- ggplot(bc_persite) +
  geom_point(data = bc_persite.f, aes(x = V2, y  = fst)) +
  theme_bw() 

p.clr <- ggplot() + 
  geom_line(data = df.sw, aes(x = location, y = LR, color = pop)) +
  geom_vline(xintercept = min(df.outliers$V2), color = "red") +
  geom_vline(xintercept = max(df.outliers$V2), color = "red") +
  theme_bw() 
pfst / p.clr

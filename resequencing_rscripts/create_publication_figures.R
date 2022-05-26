
#manhattan plots
#fst and pbs
source("summarize_outlier_windows.R")
#knit("custom_pbs_V2.Rmd") #just have to go run these
#knit("chr_by_color_fstplots_V2.Rmd") #just have to go run these

library(patchwork)
lore.label.plot <- manc.labels.pbs(lore.pbs, lore.ann)
plot.bc_fst.l <- manc.labels(bc_fst, "output/chest_genes_of_interest.csv") + ylim(0, 0.8)

plot.bc_fst.l / lore.label.plot
#ggsave("output/manhattan_figure.pdf", height = 8, width = 20) 
ggsave("output/manhattan_figure.png", height = 10, width = 20) 



# all pop comps -----------------------------------------------------------


source("plot_fst_perpop_V2.R")

ai_lo.man / ai_mo.man / ai_na.man / lo_mo.man / na_mo.man / na_lo.man
#ggsave("output/all_pops_fst_figure.pdf", height = 25, width = 35) 
ggsave("output/all_pops_fst_figure.png", height = 16, width = 18) 

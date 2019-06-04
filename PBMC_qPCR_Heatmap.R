# 1. Tell R which packages you need to use

library(ggplot2)
library(ggdendro)
library(plyr)
library(reshape2)
library(scales)
library(viridis)
library(viridisLite)
library(readr)
library(heatmaply)
library(pvclust)

# 2. Upload data *replace ###### with path to file "PBMC_miRNA_qPCR.txt"

targets <- read.delim("~/XXXX/PBMC_miRNA_corr.txt")
View(targets)

# 3. Heatmap of Pearson correlation coefficients for all qPCR data.
# (Figure 6B)

targets.a <- targets[,c(1:10)]
targets.l <- (targets.a[,2:10])
targets.x <- (targets.a[,2:10])
row.names(targets.l) <- targets.a$X
colnames(targets.l) <- targets.a$X
targets.l$gene <- targets.a$X
targets.x$order <- c(1,2,3,4,5,6,7,8,9)
targets.l$gene <- factor(targets.l$gene, levels=(targets.l$gene)[order(targets.x$order)])

mir_all <- melt(targets.l)
p <- ggplot(mir_all, aes(gene, variable)) + 
  geom_tile(aes(fill = value), colour = "white", size = 1) + theme_grey(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar")


png("miRNAfigure6b.png", height = 8500, width = 10000, res = 1200)
p
dev.off()


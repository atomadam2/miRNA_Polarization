# 1. Tell R which packages you need to use

library(ggplot2)

# 2. Upload data *replace ###### with path to file "PBMC_miRNA_qPCR.txt"

PBMC_miRNA <- read.delim("~/XXXX/PBMC_miRNA_qPCR.txt")
View(PBMC_miRNA)

# 3. Scatterplots of log transformed qPCR data for 2 miRNAs.
# (Figure 6A - miR-221-5p vs miR-222-5p)

png("miRNAfigure6a.png", height = 9000, width = 10000, res = 1200)

ggplot(PBMC_miRNA, aes(x = log2(miR221), y = log2(miR222), color=pheno)) + 
  labs(x = "log2(miR-221-5p/SNORD68)", y = "log2(miR-222-5p/SNORD68)") +
  geom_point(size = 3, aes(colour=as.factor(pheno))) +
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_color_manual(values = c("red", "black")) +
  scale_fill_manual(values = c("red", "black")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=18),
        axis.text.y = element_text(face="bold", color="#000000", size=18),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "miR-221 vs miR-222")) +
  theme(legend.position=c(.8, .2)) 

dev.off()

# (Figure 6C - miR-125a-5p vs miR-99b-5p)

png("miRNAfigure6c.png", height = 9000, width = 10000, res = 1200)

ggplot(PBMC_miRNA, aes(x = log2(mir125a5p), y = log2(miR99b5p), color=pheno)) + 
  labs(x = "log2(miR-125a-5p/SNORD68)", y = "log2(miR-99b-5p/SNORD68)") +
  geom_point(size = 3, aes(colour=as.factor(pheno))) +
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_color_manual(values = c("red", "black")) +
  scale_fill_manual(values = c("red", "black")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=18),
        axis.text.y = element_text(face="bold", color="#000000", size=18),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "miR-125a-5p vs miR-99b-5p")) +
  theme(legend.position=c(.8, .8)) 

dev.off()


# DIY - Draw scatterplot of any 2 miRNAs of your choice
# Replace XXXX and YYYY with miRNA names

# miRNA names: miR2145p, miR2143p, miR199a3p, miR99b5p, miR125a3p, miR125a5p, let7c, miR221, miR222

ggplot(PBMC_miRNA, aes(x = log2(XXXX), y = log2(YYYY), color=pheno)) + 
  labs(x = "XXXX/SNORD68", y = "YYYY/SNORD68") +
  geom_point(size = 3, aes(colour=as.factor(pheno))) +
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_color_manual(values = c("red", "black")) +
  scale_fill_manual(values = c("red", "black")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "XXXX vs YYYY")) +
  theme(legend.position=c(.8, .2)) 




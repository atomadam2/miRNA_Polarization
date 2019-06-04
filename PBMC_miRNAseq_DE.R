# Tell R which packages you need to use

library(edgeR)
library(limma)
library(WGCNA)
library(ggplot2)
library(readr)
library(rafalib)
library(devtools)

##########################
# EdgeR

# 1. Upload Raw Counts Data
#    *replace XXXXXX with path to file "All_Counts.txt"

miRNA_counts <- read_delim("~/XXXXXX/All_Counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# 2. Remove unrelated data

miRNA_counts <- miRNA_counts[-which(apply(miRNA_counts[,5:16], 1, function(x)sum(x=="0"))>2),]

# 3. Preliminary MDS plot of data
y <- DGEList(counts=miRNA_counts[,5:16], genes=miRNA_counts[,1:4])
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]

nrow(y)
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
y$samples
plotMDS(y, col=c(rep("black",1), rep("red",1)))

# 4. Create Design Matrix

Drug <- factor(c("B","B","B","L","L","L","B","B","B","L","L","L"))
Sample <- factor(c("E","E","E","E","E","E","P","P","P","P","P","P"))
data.frame(Sample=colnames(y),Sample,Drug)
design <- model.matrix(~Drug+Sample)
rownames(design) <- colnames(y)
design

# 5. Estimate Dispersion

y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion
plotBCV(y)

# 6. Calculate Differentially Expressed Genes

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
colnames(design)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(de <- decideTestsDGE(lrt))
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1,1), col="blue")

# 7. Write file with all differentially expressed genes

write.table(topTags(lrt, n=25000), file="PBMC_miRNA_DE")

# 10. Volcano plot
# (Figure 1)

df_allhits <- data.frame(topTags(lrt, n=15000))
df_allhits$threshold <- as.factor(abs( df_allhits$FDR < 0.2))

volcano1 <- ggplot(df_allhits, aes(logFC, -log10(PValue), colour=threshold))

p <- volcano1 + geom_point(size=1)+
  #xlim(-10, 10) +                          # x axis
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  #  geom_vline(aes(x = 0))+
  scale_color_manual(values = c("#0000CC", "#FF3300" )) +
  scale_fill_manual(values = c("#0000CC", "#FF3300")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "Volcano Plot of Gene Expression Changes", x = "log2 fold change", y = "-log10 p-value"))

png("miRNAFig1.png", height =3600, width = 4000, res = 600)
p
dev.off()






















#First, PCA of all macrophage miRNA samples

targets <- read_delim("~/Desktop/miRNAPolar/AlmostAllCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

targets <- targets[-which(apply(targets[,2:23], 1, function(x)sum(x=="0"))>=16),]
targets[,2:23] <- t(1000000*t(targets[,2:23])/rowSums(t(targets[,2:23])))  #Transpose twice to divide each element by column sum
targets[,2:23] <- log2(targets[,2:23] + 1)
correcting <- (targets[,2:23])
row.names(correcting) <- targets$Gene
batch1 = c("A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","C","C","C","C","C","C")
corrected <- removeBatchEffect(correcting,  batch1)

write.table(corrected, na="NA", file="AlmostAllCounts_TPM_Batched.txt")

boxplot(as.data.frame(correcting),main="Original")
boxplot(as.data.frame(corrected),main="Batch corrected")

# NO LOG, all data

pcadata <- data.frame(corrected)
#pca <- princomp(pcadata,cor=TRUE)
pcat <- prcomp(pcadata, scale. = FALSE)
summary(pcat)
loadings(pcat)
plot(pcat,type="lines")
#pca$scores
biplot(pcat)

#library(psych)

pcadata2 <- t(pcadata)
pcadata2 <- pcadata2[,colSums(pcadata2)>0]
pcaOut <- prcomp(pcadata2, scale = T, center = F)
summary(pcaOut)

explained <- as.list((pcaOut$sdev)^2 / sum(pcaOut$sdev^2))
names(explained) <- colnames(pcaOut$x)
print(plot(as.numeric(explained[1:10]), main = "baseName"))
PC <- as.data.frame(pcaOut$x)

# Label

PC$label <- rownames(pcadata2)
PC$otherLabel <- c("B","B","B","B","B","B","P","P","P","P","P","P","P","P","P","P","K","K","K","K","K","K")
PC$sample <- c("C","C","C","L","L","L","1","1","2a","2a","2c","2c","o","o","p","p","E","E","E","P","P","P")

#Plot

ggplot(PC, aes(x = PC1, y = PC2, label=sample)) + 
  geom_point(size = 2, aes(colour=as.factor(otherLabel))) + 
  geom_text( nudge_y = 0.3) +
  xlab(paste("PC1", " (", round(as.numeric(explained["PC1"]) * 100, 2), "%)", sep = "")) +
  ylab(paste("PC2", " (", round(as.numeric(explained["PC2"]) * 100, 2), "%)", sep = "")) +
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  #scale_color_manual(values = c("#0000CC", "#FF3300" )) +
  #scale_fill_manual(values = c("#0000CC", "#FF3300")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "PCA of samples")) +
  theme(legend.position=c(.2, .8))


summary(pcat)




# NO LOG - just TM
pcadata <- data.frame(targets[,20:37])
#pca <- princomp(pcadata,cor=TRUE)
pcat <- prcomp(pcadata, scale. = FALSE)
summary(pcat)
loadings(pcat)
plot(pcat,type="lines")
#pca$scores
biplot(pcat)

#library(psych)

pcadata2 <- t(pcadata)
pcadata2 <- pcadata2[,colSums(pcadata2)>0]
pcaOut <- prcomp(pcadata2, scale = T, center = F)
summary(pcaOut)

explained <- as.list((pcaOut$sdev)^2 / sum(pcaOut$sdev^2))
names(explained) <- colnames(pcaOut$x)
print(plot(as.numeric(explained[1:10]), main = "baseName"))
PC <- as.data.frame(pcaOut$x)

# Label

PC$label <- rownames(pcadata2)
PC$otherLabel <- c("0","0","0","1","1","1","0","0","0","1","1","1","0","0","0","1","1","1")
PC$sample <- c("G","G","G","H","H","H","I","I","I","J","J","J","K","K","K","L","L","L")

#Plot

ggplot(PC, aes(x = PC1, y = PC2, label=sample)) + 
  geom_point(size = 2, aes(colour=as.factor(otherLabel))) + 
  geom_text( nudge_y = 0.5) +
  xlab(paste("PC1", " (", round(as.numeric(explained["PC1"]) * 100, 2), "%)", sep = "")) +
  ylab(paste("PC2", " (", round(as.numeric(explained["PC2"]) * 100, 2), "%)", sep = "")) +
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_color_manual(values = c("#0000CC", "#FF3300" )) +
  scale_fill_manual(values = c("#0000CC", "#FF3300")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "PCA of samples before and after CQ no gametocyte")) +
  theme(legend.position=c(.9, .9))


summary(pcat)



# NO LOG - just CA
pcadata <- data.frame(targets[,2:19])
#pca <- princomp(pcadata,cor=TRUE)
pcat <- prcomp(pcadata, scale. = FALSE)
summary(pcat)
loadings(pcat)
plot(pcat,type="lines")
#pca$scores
biplot(pcat)

#library(psych)

pcadata2 <- t(pcadata)
pcadata2 <- pcadata2[,colSums(pcadata2)>0]
pcaOut <- prcomp(pcadata2, scale = T, center = F)
summary(pcaOut)

explained <- as.list((pcaOut$sdev)^2 / sum(pcaOut$sdev^2))
names(explained) <- colnames(pcaOut$x)
print(plot(as.numeric(explained[1:10]), main = "baseName"))
PC <- as.data.frame(pcaOut$x)

# Label

PC$label <- rownames(pcadata2)
PC$otherLabel <- c("0","0","0","1","1","1","0","0","0","1","1","1","0","0","0","1","1","1")
PC$sample <- c("A","A","A","B","B","B","C","C","C","D","D","D","E","E","E","F","F","F")

#Plot

ggplot(PC, aes(x = PC1, y = PC2, label=sample)) + 
  geom_point(size = 2, aes(colour=as.factor(otherLabel))) + 
  geom_text( nudge_y = 0.5) +
  xlab(paste("PC1", " (", round(as.numeric(explained["PC1"]) * 100, 2), "%)", sep = "")) +
  ylab(paste("PC2", " (", round(as.numeric(explained["PC2"]) * 100, 2), "%)", sep = "")) +
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_color_manual(values = c("#0000CC", "#FF3300" )) +
  scale_fill_manual(values = c("#0000CC", "#FF3300")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "PCA of samples before and after CQ no gametocyte")) +
  theme(legend.position=c(.9, .9))


summary(pcat)




#TPM for all everything

targets <- read_delim("~/Desktop/miRNAPolar/AllCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
targets <- targets[-which(apply(targets[,2:41], 1, function(x)sum(x=="0"))>=30),]
targets[,2:41] <- t(1000000*t(targets[,2:41])/rowSums(t(targets[,2:41])))  #Transpose twice to divide each element by column sum
targets[,2:41] <- log2(targets[,2:41] + 1)
correcting <- (targets[,2:41])
row.names(correcting) <- targets$Gene
batch1 = c("A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C")
corrected <- removeBatchEffect(correcting,  batch1)

write.table(corrected, na="NA", file="AllCounts_TPM_Batched.txt")

boxplot(as.data.frame(correcting),main="Original")
boxplot(as.data.frame(corrected),main="Batch corrected")





#Fischer

Mphagepolar <- matrix(c(9, 32, 65, 669),
                    nrow = 2,
                    dimnames =
                      list(c("DE", "Not DE"), c("Polar", "Not Polar")))

fisher.test(Mphagepolar, alternative = "greater")













library(ggplot2)
library(ggdendro)
library(plyr)
library(reshape2)
library(scales)
library(viridis)
library(viridisLite)
library(readr)
library(heatmaply)


targets <- read_delim("~/Desktop/Phospho_2012.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

targets.a <- targets[,c(2:37)]
targets.l <- log(targets.a[,1:36] +1)

row.names(targets.l) <- targets$GeneID

###
#dendogram data
x <- as.matrix(scale(targets.l))
dd.col <- as.dendrogram(hclust(dist(targets.l)))
dd.row <- as.dendrogram(hclust(dist(t(targets.l))))
dx <- dendro_data(dd.row)
dy <- dendro_data(dd.col)

# helper function for creating dendograms
ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

# x/y dendograms
px <- ggdend(dx$segments)
py <- ggdend(dy$segments) + coord_flip()

# heatmap
col.ord <- order.dendrogram(dd.col)
row.ord <- order.dendrogram(dd.row)
xx <- scale(targets.l)[col.ord, row.ord]
xxb <- data.matrix(targets.l[col.ord, row.ord])
xx_names <- attr(xx, "dimnames")
df <- as.data.frame(xxb)
colnames(df) <- xx_names[[2]]
df$gene <- xx_names[[1]]
df$gene <- with(df, factor(gene, levels=gene, ordered=TRUE))
mdf <- reshape2::melt(df, id.vars="gene")
#p <- ggplot(mdf, aes(x = variable, y = gene)) + geom_tile(aes(fill = value)) + 
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_gradient(low="blue", high="red")

p <- ggplot(mdf, aes(x = variable, y = gene)) + geom_tile(aes(fill = value)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_gradientn(colours = c("white", "lightpink",  "red"), values = c(-1,0,1))

p

# Scatter plot

mat <- matrix(unlist(dplyr::select(df,-gene)),nrow=nrow(df))
colnames(mat) <- colnames(df)[1:ncol(df)-1]
rownames(mat) <- rownames(df)

# hide axis ticks and grid lines
eaxis <- list(
  showticklabels = FALSE,
  showgrid = FALSE,
  zeroline = FALSE
)

p_empty <- plot_ly(filename="r-docs/dendrogram") %>%
  # note that margin applies to entire plot, so we can
  # add it here to make tick labels more readable
  layout(margin = list(l = 200),
         xaxis = eaxis,
         yaxis = eaxis)

subplot(px, p_empty, p, py, nrows = 2, margin = 0.01)
##differential expression 
library(edgeR)
rm(list = ls())
#setwd("path")
gene = read.csv("gene_data.csv",header = T,row.names = 1)
y = DGEList(counts = gene)
y = calcNormFactors(y)
patient = factor(c(8,8,33,33,51,51))
tissue = factor(c("N","T","N","T","N","T"))
design = model.matrix(~patient + tissue)
keep = filterByExpr(y, design, min.total.count = 200, min.count = 200)
y = y[keep,]
y = estimateDisp(y, design)
fit = glmFit(y, design)
lrt = glmLRT(fit)
de_df = lrt$table
gene_order = order(lrt$table$PValue)
de_df$DE = as.vector(decideTestsDGE(lrt))
de_df = de_df[gene_order,]
de_df$FDR = p.adjust(de_df$PValue,method = "fdr")
write.csv(de_df,"edgeR_pval.csv",quote = F)

##plot
library(ggplot2)
ggplot()+
  geom_point(data = de_df, aes(x = logFC, y = logCPM, color = as.character(DE)), size = 2)+
  theme_bw()+
  theme(aspect.ratio = 1, text = element_text(size = 22))+
  scale_color_manual(values = c("blue","grey","red"))
ggplot()+
  geom_point(data = de_df, aes(x = logFC, y = -log(PValue), color = as.character(DE)), size = 2)+
  theme_bw()+
  theme(aspect.ratio = 1, text = element_text(size = 22))+
  scale_color_manual(values = c("blue","grey","red"))
  
  
  
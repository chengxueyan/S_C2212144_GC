
library(tidyverse)
library(ggpubr)
library(reshape2)#install.packages("reshape2")
library(openxlsx)

#######批量villion_44例所有患者######

real<-read.csv(file = "data/merged.signature_edited.csv")

###################1.pCR vs. non_pCR#####################
my_comparisons = list(c("pCR", "non_pCR"),ordered = TRUE)

#单个指标循环
for (var_col in colnames(real[5:8])){
  p <-  ggviolin(real, x = "Group1", y = var_col, fill = "Group1",
                 palette = c("#FFE4C4", "#B0C4DE"),
                 add = "boxplot", add.params = list(fill = "white"))+
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
    stat_compare_means(label.y = NULL)
  
  pdf(paste0("pCR_", var_col, ".pdf"), width=4, height=4,onefile=T)
  #tiff(paste0("PD-L1_mIHC_", var_col, ".tiff"), compression="lzw",units="in",res=600,pointsize=8,width=3, height=4)
  print(p)
  dev.off()
}


library(ggsci)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggstatsplot)
library(grid)
RS_cli_marker<-read.csv("data/merged.signature_edited_barplot.csv",row.names = 1)
dir.create("Barplot")
for(idx in colnames(RS_cli_marker[,19:29])){
  if(class(RS_cli_marker[,idx])=="character")
  {
    real<-RS_cli_marker%>%drop_na(idx)
    run_cmd_str <- paste('p <- ggbarstats(real, ',  idx, ' , "Group1" ,package = "RColorBrewer",
                       palette = "Pastel2",legend.title = "Group",label.text.size = 8,bar.outline.color = "white",
                        bar.proptest = FALSE)')
    eval(parse(text=run_cmd_str))
    pdf(paste("Barplot/pCR_", idx, ".pdf"),width=4, height=4)
    grid.draw(p)
    dev.off()
  }
}
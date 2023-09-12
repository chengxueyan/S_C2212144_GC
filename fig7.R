require("getopt",quiet=T)
spec = matrix(c(
        'help','h',0,'logical','this help',
        'matrix_burden','m',1,'character','cnv burden matrix',
        'group_file','g',1,'character','group information',
        'prefix_out','p','1','character','output file(pdf & png)'
),byrow=TRUE,ncol=5);
opt=getopt(spec)
if(!is.null(opt$help) || is.null(opt$matrix_burden) || is.null(opt$group_file) || is.null(opt$prefix_out)){
        cat(paste(getopt(spec,usage = T)))
        q(status=1)
}
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
mat=read.table(opt$matrix_burden,header=TRUE,stringsAsFactors=FALSE,sep="\t")
mat<-as.data.frame(mat)
mat<-gather(
       data=mat,
       key="type",
       value="number",
       colnames(mat)[2]:colnames(mat)[ncol(mat)]
)
group_data<-as.data.frame(read.table(opt$group_file,header=1,sep="\t"))
group_data<-group_data[order(group_data$Efficacy),]
match.id1<-match(mat$sampleid,group_data$Tumor_Sample_Barcode,nomatch=0)
mat$Group<-group_data[match.id1,"Efficacy"]
col<-c("#00AFBB","#E7B800","#FC4E07","seagreen","orchid","gray66","maroon","violet","purple","royalblue","#000000","#CC0033","#FFCCCC","#CC6600","#006600")
col_t<-col[1:length(unique(mat$Group))]
y_lim = max(as.numeric(as.matrix(mat[,"number"])))+median(as.numeric(as.matrix(mat[,"number"])))
number<-length(unique(mat$Group))
if (number == 2 ){
   p<-ggplot(mat,aes(x=type,y=number,fill = Group))+ geom_boxplot()+stat_boxplot(geom='errorbar') + labs(y="score",x="",family="Times")+theme(plot.title=element_text(hjust=0.5))+ scale_fill_manual(values=c("#007FC4","#E5201F"))  +stat_compare_means(method ="wilcox.test",label = "p.format") +theme(panel.border=element_blank())+theme(panel.background=element_blank())+theme(axis.line=element_line(colour="black"))+ theme(axis.text.x=element_text(size=12,vjust=0.5,face="bold",angle=45))+theme(axis.text.y=element_text(size=12))+theme(axis.title.y=element_text(size=14))+theme(legend.text=element_text(size=12),legend.background = element_rect(fill="white",size=0.5,color="white"),legend.title=element_text(size=14),legend.key.height = unit(12, "pt"),legend.key.width = unit(12, "pt"))
}
if (number == 3 ){
   p<-ggplot(mat,aes(x=type,y=number,fill = Group))+ geom_boxplot()+stat_boxplot(geom='errorbar') + labs(y="score",x="",family="Times")+theme(plot.title=element_text(hjust=0.5))+ scale_fill_manual(values=c("#007FC4","#E5201F","#6AC239"))  +stat_compare_means(method ="kruskal.test",label = "p.format") +theme(panel.border=element_blank())+theme(panel.background=element_blank())+theme(axis.line=element_line(colour="black"))+ theme(axis.text.x=element_text(size=12,vjust=0.5,face="bold",angle=45))+theme(axis.text.y=element_text(size=12))+theme(axis.title.y=element_text(size=14))+theme(legend.text=element_text(size=12),legend.background = element_rect(fill="white",size=0.5,color="white"),legend.title=element_text(size=14),legend.key.height = unit(12, "pt"),legend.key.width = unit(12, "pt"))
}

if (number == 4 ){   
   p<-ggplot(mat,aes(x=type,y=number,fill = Group))+ geom_boxplot()+stat_boxplot(geom='errorbar') + labs(y="score",x="",family="Times")+theme(plot.title=element_text(hjust=0.5))+ scale_fill_manual(values=c("#007FC4","#E5201F","#6AC239","#F9DE8A"))  +stat_compare_means(method ="kruskal.test",label = "p.format") +theme(panel.border=element_blank())+theme(panel.background=element_blank())+theme(axis.line=element_line(colour="black"))+ theme(axis.text.x=element_text(size=12,vjust=0.5,face="bold",angle=45))+theme(axis.text.y=element_text(size=12))+theme(axis.title.y=element_text(size=14))+theme(legend.text=element_text(size=12),legend.background = element_rect(fill="white",size=0.5,color="white"),legend.title=element_text(size=14),legend.key.height = unit(12, "pt"),legend.key.width = unit(12, "pt"))
}
if (number == 5 ){  
   p<-ggplot(mat,aes(x=type,y=number,fill = Group))+ geom_boxplot()+stat_boxplot(geom='errorbar') + labs(y="score",x="",family="Times")+theme(plot.title=element_text(hjust=0.5))+ scale_fill_manual(values=c("#007FC4","#E5201F","#6AC239","#F9DE8A","#BF57B8")) +stat_compare_means(method ="kruskal.test",label = "p.format") +theme(panel.border=element_blank())+theme(panel.background=element_blank())+theme(axis.line=element_line(colour="black"))+ theme(axis.text.x=element_text(size=12,vjust=0.5,face="bold",angle=45))+theme(axis.text.y=element_text(size=12))+theme(axis.title.y=element_text(size=14))+theme(legend.text=element_text(size=12),legend.background = element_rect(fill="white",size=0.5,color="white"),legend.title=element_text(size=14),legend.key.height = unit(12, "pt"),legend.key.width = unit(12, "pt"))
}
out1_pdf<-paste(opt$prefix_out,".boxplot.pdf",sep="")
pdf(out1_pdf,height = 6,width = 8)
print(p)
dev.off()


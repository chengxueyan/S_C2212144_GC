require("getopt",quiet=T)
spec = matrix(c(
        'help','h',0,'logical','this help',
        'matrix_cni','m',1,'character','snp matrix input',
        'group_file','g',1,'character','group information',
        'prefix_out','p','1','character','output file prefix'
),byrow=TRUE,ncol=5);
opt=getopt(spec)
if(!is.null(opt$help) || is.null(opt$matrix_cni) || is.null(opt$group_file) || is.null(opt$prefix_out)){
        cat(paste(getopt(spec,usage = T)))
        q(status=1)
}
library(ggpubr)
library(ggplot2)
mat=read.table(opt$matrix_cni,header=TRUE,stringsAsFactors=FALSE,sep="\t")
mat<-as.data.frame(mat)
group_data<-as.data.frame(read.table(opt$group_file,header=1,sep="\t"))
group_data<-group_data[order(group_data$Efficacy),]
match.id1<-match(mat$sample,group_data$Tumor_Sample_Barcode,nomatch=0)
mat$Group<-group_data[match.id1,"Efficacy"]
col<-c("#00AFBB","#E7B800","#FC4E07","seagreen","orchid","gray66","maroon","violet","purple","royalblue","#000000","#CC0033","#FFCCCC","#CC6600","#006600")
col_t<-col[1:length(unique(mat$Group))]
my_comparisons<-unique(group_data$Efficacy)
y_lim=max(as.numeric(mat[,"CNI.Normal.2std"]))+1000
out1_pdf<-paste(opt$prefix_out,".boxplot.pdf",sep="")
out2_pdf<-paste(opt$prefix_out,".violin.pdf",sep="")
p1<-ggboxplot(mat,x="Group",y="CNI.Normal.2std",color="Group",palette=col_t,add="none",ylab="CNI")+ stat_compare_means(method = "wilcox.test",label = "p.format",comparisons = list(my_comparisons))+ ylab("CNI")+ggtitle("CNI Boxplot")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))+theme(axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))
pdf(out1_pdf,height = 6,width = 8)
print(p1)
dev.off()
anum <- mat[which(mat$Group == my_comparisons[1]),"CNI.Normal.2std"]
bnum <- mat[which(mat$Group == my_comparisons[2]),"CNI.Normal.2std"]
pvalue <- wilcox.test(anum,bnum)$p.value
pvalue <- round(pvalue,2)
p2 <- ggplot(data = mat, aes(x = Group, y =as.numeric(as.vector(CNI.Normal.2std)),fill=Group))+geom_violin()+
geom_boxplot(width=0.05,fill="white")+scale_fill_manual(values=c("#00AFBB","#E7B800"))+ylab("CNI")+xlab("Group")+ggtitle("CNI Boxplot")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5))+theme(axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))+ theme(legend.title=element_text(size=16),legend.text=element_text(size=14))+annotate("text",x=1.5,y=y_lim,label=paste("pvalue: ",pvalue,sep=""),size=4)
pdf(out2_pdf,height = 6,width = 8)
print(p2)
dev.off()
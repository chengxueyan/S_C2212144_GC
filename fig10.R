library(getopt)
opt_spec <- matrix(c(
        'help','h', 0, 'logical', 'help manual',
        'input', 'i', 1, 'character', 'the matrix file',
        'group', 'g', 1, 'character', 'the group file',
        'output', 'o', 1, 'character','the directory of output file'
),byrow=TRUE, ncol=5)

opt = getopt(opt_spec, commandArgs(TRUE))
library(ggplot2)
library(reshape2)
library(ggpubr)
if (!is.null(opt$help) || is.null(opt$i) || is.null(opt$g) || is.null(opt$o)){
         cat(getopt(opt_spec, usage=TRUE))
         q(save='no', status=1)
}
dat <- read.table(opt$input,header=T,row.names =1,encoding='UTF-8',check.names=F,stringsAsFactors = FALSE,sep="\t")
dat <- as.data.frame(t(dat),check.names = FASLE)
group <-read.table(opt$group,header=T,encoding='UTF-8',check.names=F,sep="\t",stringsAsFactors = FALSE)
my_comparisons <- unique(group$groups)
dat <- dat[,which(colnames(dat) %in% group$sample)]
group <- group[which(group$sample %in% colnames(dat)),]
dat$cell_type <-rownames(dat)
dat2 <- melt(dat)
matchinfo<-match(dat2$variable,group$sample,nomatch=0)
dat2$Group<-group[matchinfo,"groups"]
colnames(dat2) <- c("Gene_id","sample","value","Group")
dat2<-dat2[order(dat2$Group,decreasing = FALSE),]
pdf(paste(opt$output,"/",my_comparisons[1],"vs",my_comparisons[2],"_diff_boxplot.pdf",sep=""), width = 12, height = 6,onefile=FALSE)
if (grepl("paired",opt$group)){
        dat2$patient <- group[matchinfo,"types"]
        dat2 <- dat2[order(dat2$Group,dat2$patient),]
        p <- ggboxplot(dat2,x="Group",y="value",color = "Group",id ="Type",palette = "lancet", line.color="gray",add="jitter", line.size=0.4,facet.by = "Gene_id",short.panel.labs = FALSE)+stat_compare_means(method ="wilcox.test",label = "p.format", label.x = 1.5,paired = T)+facet_wrap(~Gene_id,scales="free",ncol=5)+ theme(axis.text.x = element_text(size = 8,color="black"),axis.text.y = element_text(size = 8,color="black")) + theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8))+theme(strip.text = element_text(size=8))+ geom_line(aes(group=patient),colour = "grey",linetype = 1,size = 0.2)
}else{
        p <- ggboxplot(dat2,x="Group",y="value",color = "Group",id ="Type",palette = "lancet", line.color="gray",add="jitter", line.size=0.4,facet.by = "Gene_id",short.panel.labs = FALSE)+ stat_compare_means(method ="wilcox.test",label = "p.format", label.x = 1.5)+ facet_wrap(~Gene_id, scales="free",ncol=5)+ theme(axis.text.x = element_text(size = 8, color="black"), axis.text.y = element_text(size = 8,color="black")) + theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8))+theme(strip.text = element_text(size=8))
}
print(p)
dev.off()


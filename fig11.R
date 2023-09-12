library(getopt)
opt_spec <- matrix(c(
        'help','h', 0, 'logical', 'help manual',
        'input', 'i', 1, 'character', 'the matrix file',
        'group', 'g', 1, 'character', 'the group file',
        'output', 'o', 1, 'character','prefix'
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
colnames(dat2) <- c("cell_type","sample","value","Group")
dat2<-dat2[order(dat2$Group,decreasing = FALSE),]
pdf(paste(opt$output,".boxplot.pdf",sep=""), width = 12, height = 10,onefile=FALSE)
p <- ggboxplot(dat2,x="Group",y="value",color = "Group",id ="Type",palette = c("#ADB6B6FF","#00468BFF","#ED0000FF"), line.color="gray",add="jitter", line.size=0.4,facet.by = "cell_type",short.panel.labs = FALSE)+facet_wrap(~cell_type,scales="free")+stat_compare_means(comparisons = combn(my_comparisons,2,simplify = FALSE),method = "wilcox.test",label = "p.format")+stat_compare_means(method ="kruskal.test",label.x = 1.5)
print(p)
dev.off()


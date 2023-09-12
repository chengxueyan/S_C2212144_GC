
require("getopt",quiet=T)
library(ggplot2)
spec = matrix(c(
	'help','h',0,'logical','this help',
	'snp_matrix','s',1,'character','snp matrix input',
	'group_file','g',1,'character','group information',
	'prefix_out','p','1','character','output file prefix'
),byrow=TRUE,ncol=5);
opt=getopt(spec)
if(!is.null(opt$h) || is.null(opt$s) || is.null(opt$g) || is.null(opt$p)){
	cat(paste(getopt(spec,usage = T)))
	q(status=1)
}
require(ComplexHeatmap,quietly=T)
mat=read.table(opt$s,fill=T,header=TRUE,stringsAsFactors=FALSE,sep="\t",row.names=1,encoding='UTF-8',check.names=F)
if (nrow(mat)>0){
        mat<-as.data.frame(mat)
        number=dim(mat)[2]
        distance=dim(mat)[1]
        group_data<-as.data.frame(read.table(opt$g,header=1,sep="\t",encoding='UTF-8',check.names=F))
        group_data<-group_data[order(group_data$Efficacy),]
        match.id1<-match(rownames(mat),group_data$Tumor_Sample_Barcode,nomatch=0)
        mat$Group<-group_data[match.id1,"Efficacy",drop=FALSE]
        mat<-mat[order(mat$Group),]
        mat.bak<-mat
        mat<-mat[,1:number,drop=FALSE]
        mat[is.na(mat)] = ""
        mat = t(as.matrix(mat))
        muttype <- c("nonsynonymous_SNV","frameshift_deletion","nonframeshift_deletion", "frameshift_insertion", "nonframeshift_insertion", "frameshift_substitution", "nonframeshift_substitution", "stopgain", "stoploss", 'UTR', 'splicing',"UNKNOWN")
        col = RColorBrewer::brewer.pal(n = length(muttype), name = 'Paired')
        col <- col[c(2, 1, 7, 4, 3, 6, 5, 8, 9, 10,11,12)]
        names(col) = muttype
        alter_fun_my = function(x, y, w, h, v) {
                n = sum(v)
                w = w * 0.95
                h = h * 0.9
                grid.rect(x, y, w*0.95, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA))
                if(n){
                grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h,gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
                }
        }
        color_n = length(unique(group_data$Efficacy))
        if(color_n==1){
                col_nw = "darkorange"
        }else if (color_n == 2){
                col_nw = c("#003399","#CC0033")
        }else{
                col_nw = RColorBrewer::brewer.pal(n=color_n, name="Set3")
        }
        group.assign<-setNames(col_nw,unique(group_data$Efficacy))
        match.id<-match(colnames(mat),group_data$Tumor_Sample_Barcode,nomatch=0)
        df<-group_data[match.id,"Efficacy",drop=FALSE]
        ha<-HeatmapAnnotation(df=df,col=list(Efficacy=group.assign))
        out_pdf<-paste(opt$p,".complexlandscape.pdf",sep="")
        ht<-oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],show_column_names=TRUE,column_title = "SNV Heatmap",column_names_gp = gpar(fontsize=6),pct_gp = gpar(fontsize=8),row_names_gp = gpar(fontsize=8), bottom_annotation=ha,alter_fun = alter_fun_my, col = col,heatmap_legend_param = list(title = "Mutation type",bottom_annotation_height=1),remove_empty_columns=FALSE)
        ht = draw(ht)
        or <- column_order(ht)
        re_order <- data.frame(or_freq = colnames(mat)[or], Group = NA, stringsAsFactors = F)
        matchid<-match(re_order$or_freq,group_data$Tumor_Sample_Barcode,nomatch=0)
        re_order$Group<-group_data[matchid,"Efficacy"]
        re_order <- re_order[order(re_order$Group, decreasing = F),]
        or <- row_order(ht)
        row.re_order <- data.frame(or_freq = rownames(mat)[or], stringsAsFactors = F)
        mat<-as.data.frame(mat)
        p<-oncoPrint(mat[row.re_order$or_freq[1:30],re_order$or_freq,drop=FALSE], get_type = function(x) strsplit(x, ";")[[1]], show_column_names=TRUE, row_order=row.re_order$or_freq[1:30], column_order=re_order$or_freq, column_title = "SNV Heatmap",column_names_gp = gpar(fontsize=6),pct_gp = gpar(fontsize=8),row_names_gp = gpar(fontsize=8), bottom_annotation=ha,alter_fun = alter_fun_my, col = col,heatmap_legend_param = list(title = "Mutation type",bottom_annotation_height=1),remove_empty_columns=FALSE)
        pdf(out_pdf,height = 8,width = 12)
        print(p)
        dev.off()
}else{
        out_pdf<-paste(opt$prefix_out,".complexlandscape.pdf",sep="")
        pdf(out_pdf)
        plot(1:5, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
        dev.off()
}
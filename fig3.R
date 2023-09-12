
require("getopt",quiet=T)
spec = matrix(c(
	'help','h',0,'logical','this help',
	'snp_matrix','s',1,'character','snp matrix input',
	'group_info','g',1,'character','Efficacy information',
	'prefix_out','p','1','character','output file prefix',
	'anno','a','1','character','annotation file',
        'kegg','k','1','character','formated pathway annotation file'
),byrow=TRUE,ncol=5);
opt=getopt(spec)
if(!is.null(opt$h) || is.null(opt$s) || is.null(opt$g) || is.null(opt$p) || is.null(opt$a) || is.null(opt$k) ){
        cat(paste(getopt(spec,usage = T)))
        q(status=1)
}

library(ggplot2)
library(ComplexHeatmap)

mat=read.table(opt$s,fill=T,header=TRUE,stringsAsFactors=FALSE,sep="\t",row.names=1,encoding='UTF-8',check.names=F)
number=dim(mat)[2]
distance=dim(mat)[1]
Efficacy_data<-as.data.frame(read.table(opt$g,header=1,sep="\t",encoding='UTF-8',check.names=F))
Efficacy_data<-Efficacy_data[order(Efficacy_data$Efficacy),]
match.id1<-match(rownames(mat),Efficacy_data$Tumor_Sample_Barcode,nomatch=0)
mat$Efficacy<-Efficacy_data[match.id1,"Efficacy",drop=FALSE]

mat<-mat[order(mat$Efficacy),]
mat<-mat[,1:number,drop=FALSE]
mat[is.na(mat)] = ""
mat = t(as.matrix(mat))

pathway=read.table(opt$a,fill=T,header=TRUE,stringsAsFactors=FALSE,sep="\t",encoding='UTF-8',check.names=F)
match1<-match(pathway$gene,rownames(mat),nomatch=0)
newdata<-as.data.frame(mat[match1,,drop=F])
match.id2<-match(rownames(newdata),pathway$gene,nomatch=0)
if(length(match.id2)>0){
        number<-nrow(newdata)
        newdata$KEGG<-pathway[match.id2,"pathway",drop=FALSE]
        drawdata<-as.data.frame(newdata[order(newdata$KEGG),1:(dim(newdata)[2]-1)])
        muttype <- c("nonsynonymous_SNV","frameshift_deletion","nonframeshift_deletion", "frameshift_insertion", "nonframeshift_insertion", "frameshift_substitution", "nonframeshift_substitution", "stopgain","stoploss",'UTR','splicing',"UNKNOWN")
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
        color_n = length(unique(Efficacy_data$Efficacy))
        if(color_n==1){
                col_nw = "darkorange"
        }else if (color_n == 2){
                col_nw = c("#003399","#CC0033")
        }else{
                col_nw = RColorBrewer::brewer.pal(n=color_n, name="Set3")
        }
        Efficacy.assign<-setNames(col_nw,unique(Efficacy_data$Efficacy))
        match.id<-match(colnames(mat),Efficacy_data$Tumor_Sample_Barcode,nomatch=0)
        df<-Efficacy_data[match.id,"Efficacy",drop=FALSE]
        ha<-HeatmapAnnotation(df=df,col=list(Efficacy=Efficacy.assign),annotation_height = unit(1, "cm"))
        out_pdf<-paste(opt$prefix_out,".pathway.anno.complexheatmap.pdf",sep="")
        pdf(out_pdf,width=18,height=8,onefile = FALSE)
        ht<-oncoPrint(drawdata, get_type = function(x) strsplit(x, ";")[[1]],show_column_names=TRUE,column_title = "SNV Heatmap",column_names_gp = gpar(fontsize=6),pct_gp = gpar(fontsize=8),row_names_gp = gpar(fontsize=8), bottom_annotation=ha,alter_fun = alter_fun_my, col = col,heatmap_legend_param = list(title = "Mutation type",bottom_annotation_height=1),remove_empty_columns=FALSE)
        ht = draw(ht)
        or <- column_order(ht)
        re_order <- data.frame(or_freq = colnames(drawdata)[or], Efficacy = NA, stringsAsFactors = F)
        matchid<-match(re_order$or_freq,Efficacy_data$Tumor_Sample_Barcode,nomatch=0)
        re_order$Efficacy<-Efficacy_data[matchid,"Efficacy"]
        re_order <- re_order[order(re_order$Efficacy, decreasing = F),]
        or <- row_order(ht)
        row.re_order <- data.frame(or_freq = rownames(drawdata)[or], stringsAsFactors = F)
        p<-oncoPrint( drawdata, get_type = function(x) strsplit(x, ";")[[1]],show_column_names=TRUE,row_order=row.re_order$or_freq,column_order=re_order$or_freq,column_title = "SNV Heatmap",column_names_gp = gpar(fontsize=6),pct_gp = gpar(fontsize=8),row_names_gp = gpar(fontsize=8), bottom_annotation=ha,alter_fun = alter_fun_my, col = col,heatmap_legend_param = list(title = "Mutation type",bottom_annotation_height=1),remove_empty_columns=FALSE)
        formatp<-read.table(opt$k,header=T,row.names=1,sep="\t")
        matchp<-match(rownames(drawdata),rownames(formatp),nomatch=0)
        newp<-formatp[matchp,]
        finalp<-newp[,which(colSums(newp) > 0),drop=FALSE]
        ha_cn = HeatmapAnnotation(cn = anno_text(colnames(finalp), rot = -45, just = "left", location = unit(1, "npc") + unit(0, "mm"), gp = gpar(fontsize = 8)), annotation_height = unit(3, "cm"))
        hk=Heatmap(finalp, col = c("0" = "white", "1" = "purple"), rect_gp = gpar(col = "grey"), show_row_names = FALSE, cluster_columns = TRUE,show_column_dend = FALSE, bottom_annotation = ha_cn, show_column_names = FALSE,show_heatmap_legend = FALSE, width = unit(7, "cm"), column_title = "KEGG pathway")
        draw(p+hk)
        dev.off()
}else{
        out_pdf<-paste(opt$prefix_out,".complexheatmap.pdf",sep="")
        pdf(out_pdf,width=12,height=8)
        plot(1:5, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
        dev.off()
}
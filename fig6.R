require("getopt",quiet=T)
spec = matrix(c(
	'help','h',0,'logical','this help',
	'snp_matrix','s',1,'character','snp matrix input',
	'group_file','g',1,'character','group information',
	'prefix_out','p',1,'character','output file',
	'anno','a',1,'character','annotation file',
        'kegg','k',1,'character','formated pathway annotation file'
),byrow=TRUE,ncol=5);
opt=getopt(spec)
if(is.null(opt$s) || is.null(opt$g) || is.null(opt$p) || is.null(opt$a) || is.null(opt$k) ){
        cat(paste(getopt(spec,usage = T)))
        q(status=1)
}

library(ggplot2)
library(ComplexHeatmap)

mat=read.table(opt$s,fill=T,header=TRUE,stringsAsFactors=FALSE,sep="\t",row.names=1,encoding='UTF-8',check.names=F)
mat<-as.data.frame(t(mat))
number=dim(mat)[2]
distance=dim(mat)[1]
group_data<-as.data.frame(read.table(opt$g,header=1,sep="\t",encoding='UTF-8',check.names=F))
group_data<-group_data[order(group_data$Efficacy),]
match.id1<-match(rownames(mat),group_data$Tumor_Sample_Barcode,nomatch=0)
mat$Group<-group_data[match.id1,"Efficacy"]
mat<-mat[order(mat$Group),]
mat<-mat[,1:number,drop=FALSE]

mat[is.na(mat)] = ""
mat = t(as.matrix(mat))
pathway=read.table(opt$a,fill=T,header=TRUE,stringsAsFactors=FALSE,sep="\t",encoding='UTF-8',check.names=F)
match1<-match(pathway$gene,rownames(mat),nomatch=0)
newdata<-as.data.frame(mat[match1,,drop=F])
match.id2<-match(rownames(newdata),pathway$gene,nomatch=0)
if(length(match.id2)>0){
        number<-nrow(newdata)
        newdata$KEGG<-pathway[match.id2,"pathway"]
        drawdata<-as.data.frame(newdata[order(newdata$KEGG),1:(dim(newdata)[2]-1)])
        muttype <- c("gain","loss")
        col <- c("tomato","steelblue")
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
        ha<-HeatmapAnnotation(df=df,col=list(Efficacy=group.assign),annotation_height = unit(1, "cm"))
        out_pdf<-paste(opt$prefix_out,".pathway.anno.complexlandscape.pdf",sep="")
        ht<-oncoPrint(drawdata, get_type = function(x) strsplit(x, ";")[[1]],
                show_column_names=TRUE,column_title = "CNV Heatmap",
                column_names_gp = gpar(fontsize=8),pct_gp = gpar(fontsize=8),
                row_names_gp = gpar(fontsize=10),bottom_annotation=ha,
                alter_fun = alter_fun_my, col = col,
                heatmap_legend_param = list(title = "Mutation type",bottom_annotation_height=1),remove_empty_columns=FALSE)
        ht = draw(ht)
        or <- column_order(ht)
        re_order <- data.frame(or_freq = colnames(drawdata)[or], Group = NA, stringsAsFactors = F)
        matchid<-match(re_order$or_freq,group_data$Tumor_Sample_Barcode,nomatch=0)
        re_order$Group<-group_data[matchid,"Efficacy"]
        re_order <- re_order[order(re_order$Group, decreasing = F),]
        or <- row_order(ht)
        row.re_order <- data.frame(or_freq = rownames(drawdata)[or], Group = NA, stringsAsFactors = F)
        p<-oncoPrint( drawdata, get_type = function(x) strsplit(x, ";")[[1]],
                show_column_names=TRUE,row_order=row.re_order$or_freq,column_order=re_order$or_freq,
                column_title = "CNV Heatmap",column_names_gp = gpar(fontsize=8),
                pct_gp = gpar(fontsize=8),row_names_gp = gpar(fontsize=10),
                bottom_annotation=ha,alter_fun = alter_fun_my, col = col,
                heatmap_legend_param = list(title = "Mutation type", bottom_annotation_height=1), remove_empty_columns=FALSE)

        formatp<-read.table(opt$k,header=T,row.names=1,sep="\t")
        matchp<-match(rownames(drawdata),rownames(formatp),nomatch=0)
        newp<-formatp[matchp,]

        finalp<-newp[,which(colSums(newp) > 0),drop=FALSE]
        ha_cn = HeatmapAnnotation(cn = anno_text(colnames(finalp), rot = -45, just = "left", 
                location = unit(1, "npc") + unit(0, "mm"), gp = gpar(fontsize = 8)), 
                annotation_height = unit(3, "cm"))
        hk=Heatmap(finalp, col = c("0" = "white", "1" = "purple"), rect_gp = gpar(col = "grey"), 
                show_row_names = FALSE, cluster_columns = TRUE,show_column_dend = FALSE, 
                bottom_annotation = ha_cn, show_column_names = FALSE,show_heatmap_legend = FALSE, 
                width = unit(7, "cm"), column_title = "KEGG pathway")
        pdf(out_pdf,width=16,height=8)
                draw(p+hk)
        dev.off()
}else{
        out_pdf<-paste(opt$prefix_out,".complexlandscape.pdf",sep="")
        pdf(out_pdf,width=12,height=8)
                plot(1:5, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
        dev.off()
}


# S_C2212144_GC
#DNApanel Fig1
#snv landscape
Rscript fig1.R -s data/snv.landscape.xls -g data/Group1_condition.txt -p test_snv

#DNApanel Fig3
#snv landscape add pathway annotation
Rscript fig3.R -s data/snv.landscape.xls -g data/Group1_condition.txt -p test_snv -a data/new_pathwayinfo.xls -k data/formated_pathway.xls

#DNApanel Fig5
##cnv landscape
Rscript fig5.R -s data/all.merged.cnv.landscape.xls -g data/Group1_condition.txt -p test_cnv

#DNApanel Fig6
#cnv landscape add pathway annotation
Rscript fig6.R -s data/all.merged.cnv.landscape.xls -g data/Group1_condition.txt -p test_cnv -a data/new_pathwayinfo.xls -k data/formated_pathway.xls

#DNApanel Fig7
#cnv burden boxplot
Rscript fig7.R -m data/cnv_burden.xls -g data/Group1_condition.txt -p cnv_burden

#DNApanel Fig8
##CNI boxplot
Rscript fig8.R -m data/Group1.CNI.txt -g data/Group1_condition.txt -p CNI

#DNApanel Fig11
#barplot
#demo data:data/merged.signature_edited_barplot.csv
Rscript fig11_2.R 
merged.signature_edited_barplot.csv
demo data: data/merged.signature_edited.csv
#violin
#demo data: data/merged.signature_edited.csv
Rscript fig11_1.R

#IHC Page10
#boxplot
Rscript fig10.R -i data/New.Cytokine.data.xls -g data/Group5_condition.txt -o ./

#IHC Page24
Rscript fig10.R -i data/New.Cytokine.data.xls -g data/Group8_paired_condition.txt -o ./

#IHC Page80
Rscript fig11.R -i data/New.Cytokine.data.xls -g data/Group3_condition.txt -o test

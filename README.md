# S_C2212144_GC
## 1. DNApanel Fig1    
### snv landscape  

```shell
Rscript fig1.R -s data/snv.landscape.xls -g data/Group1_condition.txt -p test_snv    
```
## 2. DNApanel Fig3    
### snv landscape add pathway annotation    

```shell
Rscript fig3.R -s data/snv.landscape.xls -g data/Group1_condition.txt -p test_snv -a data/new_pathwayinfo.xls -k data/formated_pathway.xls
```
## 3. DNApanel Fig5
### cnv landscape
```shell
Rscript fig5.R -s data/all.merged.cnv.landscape.xls -g data/Group1_condition.txt -p test_cnv
```

## 4. DNApanel Fig6
### cnv landscape add pathway annotation
```shell
Rscript fig6.R -s data/all.merged.cnv.landscape.xls -g data/Group1_condition.txt -p test_cnv -a data/new_pathwayinfo.xls -k data/formated_pathway.xls
```
## 5. DNApanel Fig7
### cnv burden boxplot
```shell
Rscript fig7.R -m data/cnv_burden.xls -g data/Group1_condition.txt -p cnv_burden
```
## 6. DNApanel Fig8
### CNI boxplot
```shell
Rscript fig8.R -m data/Group1.CNI.txt -g data/Group1_condition.txt -p CNI
```
## 7. DNApanel Fig11
### 7.1 barplot
### demo data:data/merged.signature_edited_barplot.csv
```shell
Rscript fig11_2.R
```
### 7.2 violin
### demo data: data/merged.signature_edited.csv
```shell
Rscript fig11_1.R
```
## 9. IHC Page10
### boxplot
```shell
Rscript fig10.R -i data/New.Cytokine.data.xls -g data/Group5_condition.txt -o ./
```
## 10. IHC Page24
```shell
Rscript fig10.R -i data/New.Cytokine.data.xls -g data/Group8_paired_condition.txt -o ./
```
## 11. IHC Page80
```shell
Rscript fig11.R -i data/New.Cytokine.data.xls -g data/Group3_condition.txt -o test
```

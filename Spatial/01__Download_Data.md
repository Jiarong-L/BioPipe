

下载 [GSM5833536](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5833536) human glioma 样本的空间数据，是 SpaceRanger 标准输出

```bash
## Download GSM5833536 -- SpaceRanger Output
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833536/suppl/GSM5833536%5FGBM4%5Fspaceranger%5Fout.tar.gz
tar -xzf GSM5833536_GBM4_spaceranger_out.tar.gz
tree GBM4_spaceranger_out
```




下载 human glioma 单细胞数据，来自 [GSE138794](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138794) 中的 scRNA-Seq 数据 (GSM4119531-GSM4119539)

```bash
mkdir scRNA
cd scRNA
for dd in 1 2 3 4 5 6 7 8 9; do curl https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM411953${dd}  | grep ftp | grep href | cut -d '=' -f3 | cut -d '"' -f2 | while read uu; do echo -e "wget $uu" ;  done;  done > url.sh
sh url.sh

for dd in 1 2 3 4 5 6 7 8 9; do mkdir GSM411953${dd}; mv GSM411953${dd}_*.gz GSM411953${dd}; done

for dd in 1 2 3 4 5 6 7 8 9; do cd GSM411953${dd}; mv *barcodes.tsv.gz barcodes.tsv.gz; mv *features.tsv.gz features.tsv.gz;  mv *matrix.mtx.gz matrix.mtx.gz; cd ..; done
```

不要忘记下载 scRNA_Seq_cell_types
```bash
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138794/suppl/GSE138794%5FscRNA%5FSeq%5Fcell%5Ftypes.txt.gz
## GSE138794_scRNA_Seq_cell_types.txt.gz
```

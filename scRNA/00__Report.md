


## 01 数据

同spatial练习，或[Scanpy笔记(pbcm3k)](https://jiarong-l.github.io/notes/Bioinfo/Blocks/Scanpy/)

```bash
pbmc3k
├── barcodes.tsv
├── genes.tsv
└── matrix.mtx
```


## 02 预处理

开始前，先进行基本处理：Filtering + SCTransform() + PCA/Clustering

a. 去除Doublet（UMI异常高的细胞），但需要注意（免疫细胞+组织细胞）的情况    
b. Imputation进行降噪    
c. 去除Contamination（破碎细胞飘入液滴）    



![DoubletFinder(2411 Singlet)](./img/02a_1.png)









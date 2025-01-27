
[流程示例一](https://djhcod.github.io/r-notes/single_cell/sc_supplementary/DecontX.html)

[Omicverse: Bulk/sc 一站式流程](https://omicverse.readthedocs.io/en/latest/index.html)


## 01 数据

同Spatial练习，或[Scanpy笔记(pbcm3k)](https://jiarong-l.github.io/notes/Bioinfo/Blocks/Scanpy/)

```bash
pbmc3k
├── barcodes.tsv
├── genes.tsv
└── matrix.mtx
```


## 02 预处理

以下操作之后（之前），还需进行基本处理：Filtering(MT/..) + SCTransform() + PCA/Clustering + 批次合并(Spatial_04a/b) + 细胞周期

建议综合考虑多个软件的结果（此练习只是单纯的分别运行代码、没有投票整合），不过也可以只用SOTA方法

### a. 去除Doublet
* 说明：（UMI异常高的细胞），但需要注意（免疫细胞+组织细胞）的情况
* 工具：Scrublet(Py), DoubletFinder(R), scDblFinder(R)

查看Doublet的分布，随后从数据中将它们去除 
![DoubletFinder(均匀分布可能是ok的？)](./img/02a_1.png)

其中采用 Doublet Simulation 的方法，对于相同细胞类型构成的 Doublet 不敏感 --- 可能因为聚类时会聚到一起，影响了 Score 的计算？？

### b. Imputation进行降噪  
* 说明：（对抗 Dropout event），当基因表达量过低时候（正常情况中位数约1500）
* 工具：[SAVER/SAVERX](https://singlecell.wharton.upenn.edu/saver-x/)，[MAGIC](https://cloud.tencent.com/developer/article/1803724)，scImpute，...
* 原理：基因并不是单独存在的，表达通路使它们关联，或可依据高表达量者推断其它基因的表达

预训练一个模型进行数据恢复（贝叶斯/Markov/AE/...），一般同时也具有降噪效果（但很难说模型会不会引入bias）

### c. 去除Contamination
* 说明：（破碎细胞的基因飘入液滴），当UMI曲线末段下降幅度不够陡峭时，当聚类图像显得有些模糊毛躁时
* 工具： [SoupX](https://github.com/constantAmateur/SoupX)）/ [DecontX](https://bioc.r-universe.dev/decontX/doc/manual.html)
* 原理：利用背景游离的RNA（需要空液滴数据）/ 使用Beyesian模型，随后从捕获的细胞数据中去除污染部分


## 03 细胞注释

自动化的SingleR注释并不一定准确，可以手动注释、使用各种DL工具

[Score_of_Celltype_i = 交集(Input_Marker,DB_CellMarker) / 并集(Input_Marker,DB_CellMarker)](https://blog.csdn.net/m0_72224305/article/details/127921124)


![SingleR注释，此例中质量似乎不太行](./img/03_1.png)

其它：[Omicverse 自动注释](https://zhuanlan.zhihu.com/p/653391043)，可选用 cellmarker/panglaodb/... 



## 04 细胞轨迹

也称拟时序分析。定义细胞类型后（且需要自定义初始状态的细胞），选定xx功能相关的基因（生物学背景/HVG/PAGA-DiffusionMap），构建xx功能相关的细胞轨迹（即 branching tree），绘制基因在轨迹上的表达（e.g.发育过程中基因开关）


绘制时，（点图+graph）主要关心其layout（e.g.[PHATE](https://zhuanlan.zhihu.com/p/143266371)），热图自己按顺序定义path中包含的node即可


* 工具：注意，并不是所有细胞都有发育关系 (e.g.体细胞-免疫细胞)
    - [scvelo 速率分析](https://www.jianshu.com/p/bfff8a4cf611) (细胞状态信息 -- unspliced/spliced mRNA，需fastq)，建议选用动态模型 ```scv.tl.velocity(adata,mode='dynamics')```
    - PAGA (细胞类型在特征空间中的距离，不一定符合生物学意义，**不建议使用**)
    - [monocle3](https://www.jianshu.com/p/c402b6588e17) (选取隔断的HVG)，方法
    - CytoTRACE (推断分化起点，细胞表达的基因数目越多 越接近干细胞) --其它方法需要指定发育起点
    - URD

* 选择输入基因集
    - 文献：与研究的生物过程、细胞类型或疾病相关的基因
    - HVG / 不同时间点或处理条件下的差异基因 / 回归、GRN得知最具有预测能力的基因 / 时间序列上有动态变化的基因
    - 注：有时降维后图像过于离散/隔断（e.g.用monocle自选的HVG），但发育过程应该是连续的，此时或许可优化基因集（e.g.使用Seurat的HVG），或者修改降维方法（推荐 [Diffusion Map](https://www.bilibili.com/video/BV1et411k7Yn/)）


## 05 细胞通讯

参考 Spatial 部分，使用Ligand-Receptor数据库（膜蛋白/细胞因子/外泌体/..?）进行注释推断

可自行从PPI/STRING中总结(CellTalkDB)，或使用各种特制工具 CellPhoneDB(人)/CellChat(人/鼠)/[NicheNet(人/鼠) ](https://www.jianshu.com/p/30c6e8a24415)

以上只是关心细胞外的 配体-受体 相互作用，NicheNet/[FlowSig](https://github.com/axelalmet/flowsig)进一步关联了通信引起的细胞内部Pathway（Q：或许可用从 细胞通讯 推断下一个时间节点的命运？）


## 06 其它

* sc-ATAC 提示候选 TFBS，但如果只有单细胞数据也可以推断 GRN，见 [GRN相关笔记](https://jiarong-l.github.io/notes/Readings/GRN/)
* inferCNV 比较目标细胞与参考细胞的基因表达模式，提示拷贝数显著的增加/降低(i.e.缺失)  ---- 评估肿瘤细胞

























## 环境

Conda中：R version 4.4.1

```bash
# chooseCRANmirror(graphics=F)
# chooseBioCmirror(graphics=F)

conda activate scenv
conda install R -y
conda install r-base -y
conda install r-essentials -y 
conda install r-seurat -y
conda install conda-forge::r-ggplot2 -y

sudo apt-get install libgit2-dev
conda install libgit2 -y
conda install r-devtools -y
# install.packages("BiocManager")


# BiocManager::install("SingleCellExperiment")    ## Don't update curl !!
# BiocManager::install("ComplexHeatmap")


conda install r-hdf5r -y
# install.packages("remotes")
# remotes::install_github("mojaveazure/seurat-disk")


# devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
# BiocManager::install('glmGamPoi')


# BiocManager::install("scDblFinder") 


conda install r-circlize -y
conda install r-enrichr -y
# BiocManager::install("celda")
conda install r-rstan
# devtools::install_github("campbio/decontX")


# options(timeout=10000)
# install.packages('SoupX')


# BiocManager::install("SingleR")
conda install r-v8
# BiocManager::install("celldex")

```




其它
```R
> packageVersion(c("Seurat"))
[1] ‘5.1.0’

> BiocManager::version()
[1] ‘3.20’


## 每次加装新包之后matrixStats就有error，需降低版本，否则SCT出错
## wget https://cran.r-project.org/src/contrib/Archive/matrixStats/matrixStats_1.1.0.tar.gz
install.packages("matrixStats_1.1.0.tar.gz",repos=NULL) 
```




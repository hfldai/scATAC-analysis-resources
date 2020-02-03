# Optimization-of-scATAC-seq-analysis
Improving dimensionality reduction, clustering, visualization and motif analysis for single-cell ATAC-seq

Public datasets of scATAC-seq:
* [10x](https://www.10xgenomics.com/resources/datasets/)
  * [5k 1:1 mixture of fresh frozen human (GM12878) and mouse (A20) cells](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_v1_hgmm_5k)
  * [10k 1:1 mixture of fresh frozen human (GM12878) and mouse (A20) cells](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_v1_hgmm_10k)
  * [5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_v1_pbmc_5k)
  * [10k Peripheral blood mononuclear cells (PBMCs) from a healthy donor](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_v1_pbmc_10k)
* [Buenrostro_2018 (GM12878/HEK293T?)](https://github.com/pinellolab/scATAC-benchmarking/tree/master/Real_Data/Buenrostro_2018_bulkpeaks/input)
* [Cusanovich_2018 (GM12878/HL-60?)](https://github.com/pinellolab/scATAC-benchmarking/tree/master/Real_Data/Cusanovich_2018/input)

* [Dataset used in SCALE](https://cloud.tsinghua.edu.cn/d/eb4371c556bc46ef8516/)
  * [Splenocyte (mouse)](https://cloud.tsinghua.edu.cn/d/eb4371c556bc46ef8516/?p=%2FSplenocyte&mode=list) |[original paper](https://www.nature.com/articles/s41467-018-07771-0) (k=12, n=3166, #peaks=77453, FACS)
  * [Forebrain (mouse)](https://cloud.tsinghua.edu.cn/d/eb4371c556bc46ef8516/?p=%2FForebrain&mode=list) (dense matrix, k=8, n=2088, , #peaks=11285, cellular indexing)
  * [Mouse Atlas (mouse)](https://cloud.tsinghua.edu.cn/d/eb4371c556bc46ef8516/?p=%2Fmouse_atlas&mode=list) (sparse matrix, k=30, n=~80,000)
  * [Leukemia](https://cloud.tsinghua.edu.cn/d/eb4371c556bc46ef8516/?p=%2FLeukemia&mode=list) (k=6, n=319, #peaks=7602, Fluidigm C1)
  * [InSilico](https://cloud.tsinghua.edu.cn/d/eb4371c556bc46ef8516/?p=%2FInSilico&mode=list) (k=6, n=828, #peaks=13668, Fluidigm C1)
  * [Breast tumor](https://cloud.tsinghua.edu.cn/d/eb4371c556bc46ef8516/?p=%2FBreast_Tumor&mode=list) (k=2, n=384, #peaks=27884, FACS)
  * [GM12878vsHEK](https://cloud.tsinghua.edu.cn/d/eb4371c556bc46ef8516/?p=%2FGM12878vsHEK&mode=list) (k=2, n=526, #peaks=12938, cellular indexing)
  * [GM12878vsHL](https://cloud.tsinghua.edu.cn/d/eb4371c556bc46ef8516/?p=%2FGM12878vsHL&mode=list) (k=2, n=597, #peaks=10431, cellular indexing)

[Benchmarking](https://github.com/pinellolab/scATAC-benchmarking)

Existing pipelines for scATAC-seq analysis
* SCALE - [source](https://github.com/jsxlei/SCALE) | [paper](https://www.nature.com/articles/s41467-019-12630-7)
* cisTopic - [source](https://github.com/aertslab/cisTopic) | [paper](https://www.nature.com/articles/s41592-019-0367-1)
* snapATAC - [source](https://github.com/r3fang/SnapATAC) | [paper](https://www.biorxiv.org/content/10.1101/615179v2)
* scABC - [source](https://github.com/SUwonglab/scABC) | [paper](https://www.nature.com/articles/s41467-018-04629-3)

[Batch effects control](http://bioconductor.org/packages/devel/bioc/vignettes/batchelor/inst/doc/correction.html)

[scRNA-seq tutorial](https://github.com/yuchaojiang/ISMB2020_SingleCellTutorial)

Labs working on sc analysis
http://buenrostrolab.com/

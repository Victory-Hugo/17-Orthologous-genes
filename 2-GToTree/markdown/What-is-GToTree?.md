
系统发育组学是一种尝试在比单个基因树更接近基因组水平上推断进化关系的实践（如果想要了解更多信息，请参阅[此处](https://astrobiomike.github.io/genomics/phylogenomics#concepts)；或者[这段视频](https://merenlab.org/2020/06/27/seminar-series-on-microbial-omics/#august-12-2020)讲解系统发育学与系统发育组学的区别）。它正成为许多生物学家工作中越来越重要的步骤。

[GToTree](https://github.com/AstrobioMike/GToTree/wiki)是一个旨在让更多研究人员能够生成系统发育基因组树以帮助指导其工作的程序。从本质上讲，它只是接收基因组，并根据指定的[单拷贝基因集](https://github.com/AstrobioMike/GToTree/wiki/SCG-sets)输出比对结果和系统发育基因组树。但我认为它的价值部分来自于：
1. 其在输入基因组格式方面的灵活性（支持核苷酸/氨基酸fasta文件、GenBank文件和/或NCBI accession），以及它将如何为我们检索和处理我们想要包含的每个单个参考基因组（而不是必须自己查找和下载每个基因组）
2. 自动化工具之间的必要任务，例如按基因长度过滤命中结果，过滤掉与目标基因命中过少的基因组，在合并之前对所有目标进行单独比对和修剪，以及将基因组标签替换为更有用的内容（即所需分类等级的谱系信息，而不仅仅是accession，和/或附加我们关心的特征标识符），以便我们可以更轻松地导航和探索最终的比对和树
3. 其可扩展性 - **GToTree可以在标准笔记本电脑上在5分钟内将200个基因组转换为树，使系统发育基因组树的迭代变得轻而易举 :)**

GToTree还附带了[15个单拷贝基因集](https://github.com/AstrobioMike/GToTree/wiki/SCG-sets)，适用于不同主要类群的系统发育组学分析。

> **开放获取的Bioinformatics Journal出版物可在[此处](https://doi.org/10.1093/bioinformatics/btz188)获得。**

---

**可以像这样运行快速[conda安装](https://github.com/AstrobioMike/GToTree/wiki/installation#conda-quickstart)：**

```
# 如果需要首先安装mamba（用于更快的conda安装）
conda install -n base -c conda-forge mamba
mamba create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda -c defaults gtotree
```

---

## 概述
上面展示了GToTree处理过程的概览图，下面有更多详细信息。有关GToTree的实用方法，请查看[示例用法页面](https://github.com/AstrobioMike/GToTree/wiki/example-usage)。有关使用GToTree的详细信息，请参阅[用户指南](https://github.com/AstrobioMike/GToTree/wiki/user-guide)。

### 输入文件 - 任何组合的fasta文件（核苷酸或氨基酸）、GenBank文件和/或NCBI组装accession
* **fasta文件** - 如果是核苷酸，将使用[prodigal](https://github.com/hyattpd/Prodigal)识别编码序列(CDS)
* **GenBank文件** - 如果在GenBank文件中注释了CDS，将提取它们，如果没有，则使用[prodigal](https://github.com/hyattpd/Prodigal)识别它们
* **NCBI组装accession** - 下载NCBI组装摘要文件，构建ftp链接以下载适当的组装，尝试下载CDS的氨基酸(AA)序列（如果存在注释），如果没有，则以fasta格式下载组装并使用[prodigal](https://github.com/hyattpd/Prodigal)识别CDS - 从NCBI通过网站或命令行生成此accession文件的示例，以及使用GToTree辅助程序搜索和使用出色的[基因组分类数据库](https://gtdb.ecogenomic.org/)的示例，在[示例页面](https://github.com/AstrobioMike/GToTree/wiki/example-usage#genomes)中展示

### 识别目标基因
* GToTree然后使用[HMMER3](http://hmmer.org/)搜索每个基因组，寻找由提供的HMM文件指定的目标基因
  * 其中15个随软件一起提供，列在[SCG-sets页面](https://github.com/AstrobioMike/GToTree/wiki/SCG-sets)上

### 估计基因组完整性/冗余度
* 使用来自HMM搜索的信息，报告每个基因组的%完整性和冗余度估计，还输出每个基因组每个目标基因的命中结果表（这是为了粗略概述，现在有更好的方法来估计基因组质量，比如[CheckM2](https://github.com/chklovski/CheckM2)）

### 可选地识别额外感兴趣的目标基因
* GToTree使用[HMMER3](http://hmmer.org/)搜索每个基因组中基于[PFam](https://pfam.xfam.org/) accession或[KO](https://www.genome.jp/kegg/ko.html) accession的任何额外感兴趣的目标基因（利用[KOFamScan](https://github.com/takaram/kofam_scan)），并生成一个[交互式生命树](https://itol.embl.de/)-兼容的映射文件，以便于可视化（请参阅[此处的示例](https://github.com/AstrobioMike/GToTree/wiki/example-usage#visualization-of-gene-presenceabsence-across-the-bacterial-domain)）

### 过滤基因命中结果和基因组
* 基于长度过滤基因 - 获取该集中所有基因的中位数，过滤掉长度不在中位数长度特定范围内的基因（默认20%）
* 如果基因组没有至少达到搜索的总基因的一定比例的命中，则过滤掉这些基因组（默认50%）

### 添加所需的间隙序列
* 为分析中保留的基因组中缺失的目标基因添加适当大小的间隙序列

### 比对、修剪、合并
* 使用[Muscle](https://www.drive5.com/muscle/)比对每个基因集
* 使用[Trimal](http://trimal.cgenomics.org/)执行自动修剪
* 合并在一起形成完整比对

### 添加分类学信息
* 将NCBI或GTDB分类学信息添加到将在树上显示的标签中

### 向标签添加额外信息
* 可以提供两列或三列的制表符分隔映射文件，第一列包含NCBI accession或输入文件名（取决于输入源），第二列包含所需的基因组标签，和/或第三列包含要附加到标签的文本（不需要提供所有输入基因组）

### 构建树
* 使用[FastTree](http://www.microbesonline.org/fasttree/)或[IQ-TREE](http://www.iqtree.org/)构建树

### 输出
主要输出包括：
* 树文件和比对文件
* 基因组摘要表，将所有修改后的标签映射到原始基因组ID，完整性/冗余度估计，以及任何可用的分类学信息
* 显示每个基因组每个单拷贝基因的命中次数的表格，以及如果有任何额外搜索的基因的表格
* 报告在哪个步骤过滤了什么（如果有的话）
* 一个[引用文件](https://github.com/AstrobioMike/GToTree#citation-information)，包含特定于运行的所有GToTree使用的程序的引用信息，以帮助正确引用👍

## 引用信息

GToTree依赖于许多优秀的程序。除了所有其他输出外，它还将生成一个`citations.txt`文件，其中包含针对每次运行的引用信息，说明了它所依赖的所有程序。请确保适当引用开发者🙂

以下是运行生成的`citations.txt`文件示例，以及我在方法中如何引用它：

```
GToTree v1.6.31
Lee MD. GToTree: a user-friendly workflow for phylogenomics. Bioinformatics. 2019; (March):1-3. doi:10.1093/bioinformatics/btz188

Prodigal v2.6.3
Hyatt, D. et al. Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics. 2010; 28, 2223–2230. doi.org/10.1186/1471-2105-11-119

HMMER3 v3.3.2
Eddy SR. Accelerated profile HMM searches. PLoS Comput. Biol. 2011; (7)10. doi:10.1371/journal.pcbi.1002195

Muscle v5.1
Edgar RC. MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping. bioRxiv. 2021. doi.org/10.1101/2021.06.20.449169

TrimAl v1.4.rev15
Gutierrez SC. et al. TrimAl: a Tool for automatic alignment trimming. Bioinformatics. 2009; 25, 1972–1973. doi:10.1093/bioinformatics/btp348

TaxonKit v0.9.0
Shen W and Ren H. TaxonKit: a practical and efficient NCBI Taxonomy toolkit. Journal of Genetics and Genomics. 2021. doi.org/10.1016/j.jgg.2021.03.006

FastTree 2 v2.1.11
Price MN et al. FastTree 2 - approximately maximum-likelihood trees for large alignments. PLoS One. 2010; 5. doi:10.1371/journal.pone.0009490
```

**基于上述引用输出的示例方法文本（请根据您的运行情况进行修改）**
> *古细菌系统发育基因组树使用GToTree v1.6.31 (Lee 2019)生成，使用预先打包的古细菌单拷贝基因集（76个目标基因）。简要地说，prodigal v2.6.3 (Hyatt et al. 2010)用于预测以fasta文件形式提供的输入基因组上的基因。使用HMMER3 v3.2.2 (Eddy 2011)识别目标基因，使用muscle v5.1 (Edgar 2021)进行单独比对，使用trimal v1.4.rev15 (Capella-Gutiérrez et al. 2009)进行修剪，并在使用FastTree2 v2.1.11 (Price et al. 2010)进行系统发育估计之前合并。TaxonKit (Shen and Ren 2021)用于将完整谱系连接到分类ID。*
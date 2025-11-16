

# 欢迎来到GToTree Wiki！

[GToTree](https://github.com/AstrobioMike/GToTree)是一个用户友好的系统发育组学工作流，旨在让更多研究人员能够创建系统发育基因组树，并使系统发育基因组树的迭代过程更容易。开放获取的Bioinformatics Journal出版物可在[此处](https://doi.org/10.1093/bioinformatics/btz188)获得。GToTree可以安装并在Mac或Linux机器上运行，也可以在Windows的Windows Subsystem for Linux环境中运行👍

**有关这一切的概述，请参阅["什么是GToTree？"页面](what-is-gtotree%3F)**。或者，要了解GToTree的一些实用方法，可以浏览[示例用法页面](example-usage)。

---

**可以像这样运行快速[conda安装](https://github.com/AstrobioMike/GToTree/wiki/installation#conda-quickstart)：**

```
conda create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda gtotree
```

---

# Wiki内容

* [**什么是GToTree？**](what-is-gtotree%3F)
* [**安装**](installation)
  * [Conda快速开始！](installation#conda-quickstart)
  * [更新到新版本](installation#updating-to-a-newer-version)
  * [无需Conda的安装](installation#installation-without-conda)
* [**示例用法**](example-usage)
  * [Alteromonas示例](example-usage#alteromonas-example)
  * [将GToTree与基因组分类数据库(GTDB)结合使用](example-usage#using-gtotree-with-the-genome-taxonomy-database-gtdb)
  * [细菌间基因存在/缺失的可视化](example-usage#visualization-of-gene-presenceabsence-across-the-bacterial-domain)
  * [生命树示例](example-usage#tol-example)
  * [将比对和分区文件与其他程序一起使用](example-usage#using-the-alignment-and-partitions-file-with-another-program)
* [**用户指南**](user-guide)
  * [必需输入](user-guide#required-inputs)
  * [输出](user-guide#outputs)
  * [可选参数](user-guide#optional-arguments-and-parameters)
  * [运行程序的选项设置](user-guide#options-set-for-programs-run)
  * [基因组完整性和冗余度估计](user-guide#genome-completeness-and-redundancy-estimations)
  * [所有使用程序的引用信息](user-guide#citation-information)
* [**单拷贝基因集**](scg-sets)
  * [代码](scg-sets#code)
* [**注意事项**](things-to-consider)
  * [关于系统发育组学工作流程概念的重要警告](things-to-consider#an-important-caveat-on-the-idea-of-a-workflow-for-phylogenomics)
  * [何时使用和何时不使用GToTree](things-to-consider#when-to-use-gtotree-and-when-not)
    * [GToTree对分配分类学有用吗？](things-to-consider#is-gtotree-useful-for-assigning-taxonomy)
    * [应该使用GToTree来估计基因组/MAG/bin质量吗？](things-to-consider#should-gtotree-be-used-for-estimating-genomemagbin-quality)
  * [考虑使用"代表性"基因组](things-to-consider#consider-using-representative-genomes)
  * [处理许多基因组](things-to-consider#working-with-many-genomes)
  * [按长度筛选基因](things-to-consider#filtering-hits-by-gene-length)
  * [按目标命中率筛选基因组](things-to-consider#filtering-genomes-by-fraction-of-hits-to-targets)
  * [最佳命中模式](things-to-consider#best-hit-mode)
  * [GToTree适用于真核生物吗？](things-to-consider#are-eukaryotic-genomes-appropriate-for-use-with-gtotree)

---

<p align="center">
<a href="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GToTree-Overview-main.png"><img src="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GToTree-Overview-main.png"></a>
</p>
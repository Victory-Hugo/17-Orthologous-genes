---

这个页面展示了几个完整的示例。大多数示例的数据文件、运行日志和结果（除了大型基因/存在性可视化的示例）都随GToTree一起打包提供。

---

# 示例

### [1. ***Alteromonas* 示例**](example-usage#alteromonas-example)
- 在这里，我们假设我们有一个新获得的感兴趣的基因组，并且我们想看看它与已知的亲缘关系最近的基因组在进化图谱中处于什么位置。
### [2. **将 GToTree 与基因组分类数据库 (GTDB) 结合使用**](example-usage#using-gtotree-with-the-genome-taxonomy-database-gtdb)
- 我们可以通过提供 `-D` 选项来添加 [GTDB](https://gtdb.ecogenomic.org/) 分类信息，而不是 NCBI 分类信息。
- 我们还可以使用 `gtt-get-accessions-from-GTDB` 程序，根据 GTDB 分类法搜索我们想要包含在系统发育树中的基因组，如 [这个例子](example-usage#using-gtotree-with-the-genome-taxonomy-database-gtdb) 所示。
### [3. **细菌域中基因存在/缺失的可视化**](example-usage#visualization-of-gene-presenceabsence-across-the-bacterial-domain)
- 有时，可视化一个基因或性状在感兴趣的进化支或整个域中的存在/缺失情况是很有用的，因为这可能会揭示其进化分布的有趣模式。
### [4. **将比对和分区文件与另一个程序一起使用**](example-usage#using-the-alignment-and-partitions-file-with-another-program)
- 这展示了我们可以利用 GToTree 生成的比对和分区文件，使用出色的 [IQ-TREE](http://www.iqtree.org/) 程序创建混合模型树的一种方法。
### [5. **生命树示例**](example-usage#tol-example)
- 在这里，我们生成了 [GToTree 原始论文](https://doi.org/10.1093/bioinformatics/btz188) 中展示的三域树。

---
---

# *Alteromonas* 示例

在这里，我们假设我们有一个新获得的感兴趣的基因组，并且我们想看看它与已知的亲缘关系最近的基因组在进化图谱中处于什么位置。这个示例最初完成时的确切文件可在 [此处](https://figshare.com/articles/dataset/GToTree_Alteromonas_example_run_data/19372313) 获得，但也可以按照下面的步骤进行操作（不过由于有更多的基因组可用，情况会略有不同）。

**将我们的新基因组（以 fasta 格式提供）与 Alteromonas 参考基因组（以 NCBI accession 提供）一起构建成系统发育树大约需要 6 分钟。**

## 设置

在我的博士期间，我有幸研究了一种非常酷的、能够固氮的海洋蓝细菌，名为 *Trichodesmium*。这种生物似乎只与其他微生物群永久共生——没有纯培养的（无菌的）*Trichodesmium*，没有它的伙伴，它似乎就不开心。这项工作中出现的一种高度保守的生物（这里指的是在所有 *Trichodesmium* 样本中一致存在，但在非 *Trichodesmium*“对照”样本中不存在）是一个 *Alteromonas* **宏基因组组装基因组**（MAG）——部分有趣的工作发表在 [这篇论文](https://www.nature.com/articles/ismej201749) 中。

**在这里，我们将使用 GToTree 生成一个新的系统发育树，包含这个“新的”*Alteromonas* MAG 和 NCBI RefSeq 数据库中可用的所有 [完整](https://www.ncbi.nlm.nih.gov/assembly/help/#glossary) *Alteromonas* 基因组。**

## 生成输入

### 基因组

本示例的输入将是：1) 我们的“新”MAG 的 fasta 文件；2) NCBI RefSeq 中 *Alteromonas* 基因组的 accession 列表；3) 一个用作根的 α-变形菌。

**1)** 我们可以从 NCBI 下载 MAG 的 fasta 文件并解压缩，如下所示：

```bash
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/271/865/GCA_002271865.1_ASM227186v1/GCA_002271865.1_ASM227186v1_genomic.fna.gz | gunzip - > GCA_002271865.1.fa
```

<a id="accessions"></a>
**2)** 我们可以在他们的网站上搜索 [NCBI 的组装数据库](https://www.ncbi.nlm.nih.gov/assembly/)，使用以下搜索字符串获取所有 RefSeq *Alteromonas* 基因组的 accession：`Alteromonas[ORGN] AND "latest refseq"[filter] AND "complete genome"[filter]`（当这篇文章在 2019 年 1 月 1 日整理时，返回了 31 个结果）。您可以通过选择右上角的“Send to:”并设置如下所示的选项来下载摘要文件：

<img src="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GToTree-alteromonas-search.png">

点击创建文件将把这些下载为“assembly_result.txt”。在这里，我们使用的是 RefSeq 组装 accession（那些以“GCF_...”开头的，但 GToTree 也处理 GenBank 组装 accession（那些以“GCA_...”开头的）。在我们这里，我们将从第三列获取 RefSeq accession，但如果您按照此步骤操作并希望使用 GenBank 中可用但 RefSeq 中不可用的基因组，您应该获取第一列。

在我们的例子中，我们可以使用以下命令获取 RefSeq accession：

```bash
tail -n +2 assembly_result.txt | cut -f 3 > alteromonas_refseq_accessions.txt
```

如果您不熟悉这行代码，并且想要更好地了解如何在命令行工作，一个很好的起点是 [这里](https://astrobiomike.github.io/unix/unix-intro) :)

---
> **注意**：这个过程的这一部分（生成我们想要的参考基因组的 accession 列表）实际上可以在命令行完全完成。[EDirect](https://dataguide.nlm.nih.gov/edirect/documentation.html) 是一个用于访问 NCBI 数据库的命令行工具。它可能有一定的学习曲线，但如果您经常使用 NCBI 数据库的信息，绝对值得学习。如果使用 conda 安装，EDirect 会与 GToTree 一起安装。使用 EDirect，可以在命令行创建这个 accession 文件，如下所示：`esearch -query 'Alteromonas[ORGN] AND "latest refseq"[filter] AND "complete genome"[filter] AND (latest[filter] AND all[filter] NOT anomalous[filter])' -db assembly | esummary | xtract -pattern DocumentSummary -def "NA" -element AssemblyAccession > alteromonas_refseq_accessions.txt`。如果您感兴趣，我在 [这个页面](https://astrobiomike.github.io/unix/ncbi_eutils) 上有更多使用 EDirect 的例子 🙂

---

**3)** 为了有一个外群来给树定根，并且为了合并一个 GToTree 也可以作为输入处理的 GenBank 文件，我们将选取一种 α-变形菌：

```bash
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/365/GCF_000011365.1_ASM1136v1/GCF_000011365.1_ASM1136v1_genomic.gbff.gz | gunzip - > GCF_000011365.1.gbff
```

### 为特定基因组贴标签的映射文件

通常，在树中为特定基因组设置特定标签是很有帮助的。在这种情况下，我们可能希望将我们的“新”MAG 标记为“Our_Alteromonas_MAG”，而不仅仅是“GCA_002271865.1”，我们可能希望将我们的根标记为“Alpha_root”，而不仅仅是“GCF_000011365.1”。GToTree 使用 [TaxonKit](https://bioinf.shenwei.me/taxonkit/) 向任何具有相关分类信息的基因组添加谱系信息（无论是作为 NCBI accession 还是 GenBank 文件提供），但我们也可以交换我们知道自己关心的特定基因组的标签，以便更容易找到它们。

要做到这一点，我们只需要提供一个 2 列的、制表符分隔的文件，其中第一列是初始基因组 ID（这将是 NCBI accession 或文件名，取决于基因组是如何提供的）。可以在任何地方创建，但这里是为这个例子快速创建一个：

```bash
printf "GCA_002271865.1.fa\tOur_Alteromonas_MAG\nGCF_000011365.1.gbff\tGCF_000011365.1_Alpha_Outgroup\n" > genome_to_id_map.tsv
```

它看起来像这样：

```bash
cat genome_to_id_map.tsv
```
```
GCA_002271865.1.fa	Our_Alteromonas_MAG
GCF_000011365.1.gbff	GCF_000011365.1_Alpha_Outgroup
```
---
> **注意**：给映射文件中列出的基因组提供的用户自定义标签（通过 `-m` 选项传递给程序）将始终优先于任何自动谱系交换。

---

## 运行 GToTree

accessions 文件可以按原样提供，但为了告诉 GToTree 要处理哪些 fasta 和 genbank 文件，我们需要将它们的名称（或 [路径](https://astrobiomike.github.io/unix/getting-started#the-unix-file-system-structure)）放入文件中。在这种情况下，这样做就可以了：

```bash
ls *.fa > fasta_files.txt
ls *.gbff > genbank_files.txt
```

现在我们已经准备好运行 GToTree 了 :)

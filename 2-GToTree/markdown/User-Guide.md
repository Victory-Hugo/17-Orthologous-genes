

# 用户指南目录

* [**必填输入**](user-guide#required-inputs)
  * [**输入基因组**](user-guide#input-genomes)
  * [**指定要使用的单拷贝基因集**](user-guide#specifying-which-single-copy-gene-set-to-use)
* [**输出**](user-guide#outputs)
  * [**主要输出文件**](user-guide#primary-output-files)
  * [**报告输出文件**](user-guide#report-output-files)
* [**可选参数**](user-guide#optional-arguments-and-parameters)
* [**运行程序的选项设置**](user-guide#options-set-for-programs-run)
* [**基因组完整性和冗余度估计**](user-guide#genome-completeness-and-redundancy-estimations)
* [**所有使用程序的引用信息**](user-guide#citation-information)

---

> **注意：** 不带任何参数运行 `GToTree` 将提供帮助菜单。

---

# 必填输入
GToTree 所需的最低输入是指定您要整合的基因组（通过 NCBI Accessions、GenBank 文件和/或核苷酸或氨基酸 fasta 文件的任意组合提供），以及指定要使用的单拷贝基因集。

## 输入基因组
输入基因组可以指定为 NCBI 组装 accession、GenBank 文件和/或 fasta 文件的任意组合。

### NCBI Accessions
您可以通过向 **`-a`** 参数提供一个包含 NCBI 组装 accession 的单列文件，来指定您想要整合的哪些 NCBI 存档基因组。这个文件可以通过搜索 NCBI 网站并下载结果表来"手动"创建，也可以通过使用 [Entrez-Direct](https://dataguide.nlm.nih.gov/edirect/install.html) 在命令行生成 - 两种方法的示例都在 [此处的示例页面](https://github.com/AstrobioMike/GToTree/wiki/example-usage#accessions) 中提供。

* 提供的 accession 可以有版本号（accession 中 "." 后面的内容，例如 GCF_000153765.1），也可以没有版本号（例如 GCF_000153765）。在没有提供版本的情况下，GToTree 会自动采用该 accession 的最新发布版本。
* 如果提供的任何 accession 在 NCBI 上找不到，它们将在运行开始时打印到屏幕上，并在输出目录中的 "NCBI_accessions_not_found.txt" 文件中报告。
* 一个示例输入 accession 文件可以在 GToTree 子目录中找到：`GToTree/test_data/ncbi_accessions.txt`。

### GenBank 文件
要指定要包含哪些 GenBank 文件，您需要提供一个单列文件，其中包含每个您想要整合的 GenBank 文件的文件名（或 [路径](https://astrobiomike.github.io/unix/basics#absolute-vs-relative-path)）。这会传递给 **`-g`** 参数。
* 一个示例文件可以在 GToTree 子目录中找到：`GToTree/test_data/genbank_files.txt`。

### Fasta 文件
核苷酸 fasta 文件的提供方式与 GenBank 文件类似，但会传递给 **`-f`** 参数。您需要提供一个单列文件，其中包含每个您想要整合的 fasta 文件的文件名（或 [路径](https://astrobiomike.github.io/unix/basics#absolute-vs-relative-path)）。
* 一个示例文件可以在 GToTree 子目录中找到：`GToTree/test_data/fasta_files.txt`。

### 氨基酸文件
氨基酸文件的提供方式与核苷酸 fasta 文件类似，但会传递给 **`-A`** 参数。您需要提供一个单列文件，其中包含每个您想要整合的氨基酸 fasta 文件的文件名（或 [路径](https://astrobiomike.github.io/unix/basics#absolute-vs-relative-path)）。
* 一个示例文件可以在 GToTree 子目录中找到：`GToTree/test_data/amino_acid_files.txt`。

## 指定要使用的单拷贝基因集
GToTree 还需要知道要使用哪个 SCG 集 - 通过 **`-H`** 参数传递。程序随附了 14 个 SCG 集，存储在 `hmm_sets` 子目录中（在 [此处](https://github.com/AstrobioMike/GToTree/wiki/scg-sets) 有更详细的讨论）。如果您按照 [conda 快速开始安装说明](https://github.com/AstrobioMike/GToTree/wiki/installation#conda-quickstart) 进行安装，或者自己按照 [此处](https://github.com/AstrobioMike/GToTree/wiki/installation#add-path-to-included-hmm-files) 所述设置了适当的环境变量，您可以通过单独运行 **`gtt-hmms`** 来查看可用的 HMM 文件（您不需要指定 HMM 文件的完整路径，只需指定 `gtt-hmms` 打印的名称，例如 `-H Bacterial.hmm`）。

# 输出
每次 GToTree 运行都会创建一个输出目录来保存所有输出文件。默认名为 "GToTree_output"，但可以通过 **`-o`** 参数指定，下面包含 "GToTree_output" 的文件名将相应更改。

### 主要输出文件

#### 树文件
* **Aligned_SCGs.tre**
  * 最终树文件，采用 newick 格式。
  * FastTree 报告 "本地支持值"，作为内部节点上的标签，以估计树中每个分支的可靠性。您可以在他们的用户页面 [此处](http://www.microbesonline.org/fasttree/#Support) 找到更多信息。
  * IQ-TREE 报告超快 bootstrap (UFBoot) 支持值。他们的 [帮助页面](http://www.iqtree.org/doc/Frequently-Asked-Questions) 指出，95% 的值表示该分支有 95% 的概率是真实的。
  * 如果使用 `-N` 选项运行，则不会生成树，只会生成比对结果。

#### 比对文件
* **Aligned_SCGs.faa**
  * 比对文件，采用 fasta 格式。
* **Aligned_SCGs_mod_names.faa**（如果 [TaxonKit](https://bioinf.shenwei.me/taxonkit/) 被用来向标签添加谱系信息 - 通过 `-t` 参数指定）

#### 分区文件
* **Partitions.txt**
  * 一个与能够为每个基因使用不同模型的树程序兼容的分区文件。例如，参见 [iqtree 此处的信息](http://www.iqtree.org/doc/Advanced-Tutorial)。

#### 基因组摘要信息
* **Genomes_summary_info.tsv**
  * 每个基因组的摘要信息的制表符分隔表，包括以下列：

|列号|名称|内容|
|:-----:|:---:|---|
| 1 | assembly_id | 输入组装 ID（根据输入源，是 accession 或基础文件名） |
| 2 | label | 输出树文件中分配给该基因组的标签 |
| 3 | label_source | 标签的来源 |
| 4 | taxid | 如果基因组是通过 NCBI accession 或带有 taxid 信息的 GenBank 提供的，则为 NCBI taxid |
| 5 | num_SCG_hits | 对目标 HMM 的基因命中次数 |
| 6 | uniq_SCG_hits | 对目标 HMM 的独特基因命中次数 |
| 7 | perc_comp | 基于目标 HMM 的估计完成百分比 |
| 8 | perc_redund | 基于目标 HMM 的估计冗余百分比 |
| 9 | num_SCG_hits_after_len_filt | 长度过滤后对目标 HMM 的基因命中次数 |
| 10 | in_final_tree | 是或否，该基因组是否最终出现在最终树中 |

> 根据是否添加了分类学信息，可能会有更多包含分类学信息的列。

#### 每个基因组的 SCG 命中计数
* **All_genomes_SCG_hit_counts.tsv**
  * 一个制表符分隔的文件，其中第一列包含每个基因组 ID，其余列包含每个基因组对每个目标基因的命中次数。

### 报告输出文件
只有在需要时才会写入报告文件。例如，如果一个基因组由于对目标基因的命中太少而被从分析中删除，文件 "Genomes_removed_for_too_few_hits.tsv" 将被创建。但如果没有基因组因此原因被删除，该文件将不会生成。因此，您不应该期望在任何特定运行后找到所有这些文件。这些文件将包含在输出子目录 `run_files/` 中。

**Redundant_input_accessions.txt**
  * 如果输入 NCBI accession 文件中有重复的 accession，它们将在此处报告。

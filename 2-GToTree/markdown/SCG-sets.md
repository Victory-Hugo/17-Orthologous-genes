GToTree 包含 15 个单拷贝基因 (SCG) 集合的 HMM，涵盖多个分类等级。其中 14 个是按照以下描述新创建的，还有一个包含所有 3 个域的集合（来自 [Hug et al. 2016](https://www.nature.com/articles/nmicrobiol201648)）。

我们可以通过单独运行 **`gtt-hmms`** 来查看可用的 HMM 文件。

以下是当前可用的 HMM 集合：

|HMM 集合|基因目标数量|
|----|:----:|
|Actinobacteria|138|
|Alphaproteobacteria|117|
|Archaea|76|
|Bacteria|74|
|Bacteria_and_Archaea|25|
|Bacteroidetes|90|
|Betaproteobacteria|203|
|Chlamydiae|286|
|Cyanobacteria|251|
|Epsilonproteobacteria|260|
|Firmicutes|119|
|Gammaproteobacteria|172|
|Proteobacteria|119|
|Tenericutes|99|
|Universal (Hug et al.) |16|

有关每个集合的更多信息，请查看 [此处](https://github.com/AstrobioMike/GToTree/blob/master/hmm_sets/hmm-sources-and-info.tsv) 的 "hmm-sources-and-info.tsv" 文件，以及 GToTree 安装位置（可通过运行 `gtt-hmms` 查看）。

以下是生成新集合的大致思路，并在下面展示了生成细菌集合的精确代码。

1) 下载了来自版本 32.0 的所有 17,929 个 Pfam。过滤掉那些平均覆盖不到 50% 构建该 Pfam 的 HMM 所依据的潜在蛋白质序列的 Pfam（这确保不会包含来自同一源蛋白质的两个 Pfam）。剩下 8,924 个 Pfam。

2) 从 NCBI RefSeq 数据库中下载了所有的基因组（在 2018 年 4 月 23 日）。根据 GTDB 的分类，每个属中选择了不超过 5 个成员，并且只保留了 "完成" 状态的基因组。这为真细菌留下了 2,000 个基因组，为古菌留下了 172 个基因组。

3) 对于所有选定的基因组，将 Pfam HMMs 与这些基因组的预测蛋白质序列进行比对。

4) 对于每个 Pfam，将 Pfam 的 "拷贝数"（即每个基因组与该 Pfam 的 HMM 有至少一个命中的次数）与分类谱系进行关联。然后选择那些在每个分类等级的 1000 个基因组中平均有超过 0.95 个拷贝且少于 1.05 个拷贝的 Pfam。

5) 然后，对所有分类等级重复此过程，为每个分类等级留下了高度保守的单拷贝基因集。

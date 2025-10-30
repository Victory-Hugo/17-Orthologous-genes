# 查询基因

在 [NCBI](https://www.ncbi.nlm.nih.gov/nuccore?itool=toolbar) 查询感兴趣的基因名称。以 `LolD` 基因为例，这是一个与脂蛋白转运相关的基因。
按照如下操作下载表格：
![image.png|500](https://picturerealm.oss-cn-chengdu.aliyuncs.com/obsidian/20251023202648973.png)
![image.png|500](https://picturerealm.oss-cn-chengdu.aliyuncs.com/obsidian/20251023202704904.png)

# 获得表格
获得的表格格式如下所示：
![image.png|500](https://picturerealm.oss-cn-chengdu.aliyuncs.com/obsidian/20251023202831992.png)
其中，`GeneID` 是该基因的**唯一性编号**。尽管上图展示的都是 `LolD` 基因，但是却有不同的 `GeneID`，这是因为：
1. **来源于不同物种**
	- 不同细菌物种中的`LolD`基因会被分配不同的GeneID
	- 即使基因名称相同，但物种不同，数据库会给予独立的编号
2. **同一物种的不同菌株**
	- 同一物种的不同菌株可能携带序列略有差异的`LolD`基因
	- 每个菌株的基因会获得独立的GeneID
3. **基因的不同变异体**
	- 同一基因可能存在不同的等位基因变异
	- 序列差异导致数据库分配不同编号
4. **数据库收录的独立性**
	- NCBI等数据库按照提交的独立序列分配GeneID
	- 即使功能相同，每个独立提交的序列都会获得唯一编号

# 提取 `GeneID`
使用 `Excel` 或者命令行工具提取 `GeneID` 到一个单独的 `txt` 文件。使用 `1-数据获取/script/0-GeneID→FASTA.sh` 下载对应的 `fasta` 文件。
# GTDB-Tk 工具配置和使用指南

## 一、环境要求

GTDB-Tk 是一个用于微生物基因组分类的工具，需要以下依赖：
- Python 3.8+
- HMMER
- Pplacer
- FastTree
- Muscle (可选)

## 二、安装步骤

### 1. 安装 Conda 环境

```bash
# 创建 gtdbtk 虚拟环境
conda create -n gtdbtk python=3.10 -y

# 激活环境
conda activate gtdbtk

# 安装必要的依赖
conda install -c bioconda gtdbtk hmmer pplacer fasttree -y
```

### 2. 解压数据库文件

```bash
cd /mnt/d/3-GTDB-Database/data

# 查看文件大小
ls -lh gtdbtk_data.tar.gz

# 解压数据库（需要较长时间）
tar -xzf gtdbtk_data.tar.gz

# 解压后会生成 release### 目录（如：release220 等）
ls -la
```

### 3. 配置环境变量

将以下内容添加到 `~/.bashrc` 或 `~/.bash_profile`：

```bash
# 设置 GTDB-Tk 数据目录
export GTDBTK_DATA_PATH=/mnt/d/3-GTDB-Database/data/release226

# 其中 release### 是解压后生成的目录名称
# 例如：export GTDBTK_DATA_PATH=/mnt/d/3-GTDB-Database/data/release220
```

然后执行：
```bash
source ~/.bashrc
```

### 4. 验证安装

```bash
# 激活环境
conda activate gtdbtk

# 检查 GTDB-Tk 安装
gtdbtk --version

# 设置数据路径（如果未配置环境变量）
export GTDBTK_DATA_PATH=/mnt/d/3-GTDB-Database/data/release###

# 验证数据库
gtdbtk check_install_dir 
gtdbtk test --out_dir gtdbtk_test  # 测试结果存储在gtdbtk_test目录

```

## 三、数据库信息

GTDB-Tk 数据库包含：
- **参考基因组集合** - 用于分类的标准菌株
- **多序列比对** - 16S rRNA 和蛋白编码基因的比对
- **系统树** - 预构建的进化树
- **分类信息** - 完整的分类学注释

## 四、基本使用命令

### 1. 完整分类流程

```bash
# 激活环境
conda activate gtdbtk


# 运行完整分类
gtdbtk classify_wf \
    --genome_dir /path/to/genomes \
    --out_dir ./gtdbtk_output \
    --cpus 4 \
    --extension fasta
```

### 2. 只进行识别（Identify）

```bash
gtdbtk identify \
    --genome_dir /path/to/genomes \
    --out_dir ./gtdbtk_output \
    --cpus 4 \
    --extension fasta
```

### 3. 只进行分类（Classify）

```bash
gtdbtk classify \
    --genome_dir /path/to/genomes \
    --out_dir ./gtdbtk_output \
    --cpus 4 \
    --extension fasta
```

### 4. 构建进化树（Align + Tree）

```bash
gtdbtk align \
    --identify_dir ./gtdbtk_output \
    --out_dir ./gtdbtk_output \
    --cpus 4

gtdbtk infer \
    --msa_dir ./gtdbtk_output/align \
    --out_dir ./gtdbtk_output \
    --cpus 4
```

## 五、输出文件说明

```
gtdbtk_output/
├── identify/
│   ├── intermediate_results/
│   ├── *.summary.tsv          # 识别结果摘要
│   └── *.markers_summary.tsv   # 标记基因摘要
├── classify/
│   ├── *.summary.tsv          # 分类结果摘要
│   └── tree/
│       ├── *.tree             # 进化树文件
│       └── *.decorated.tree   # 装饰后的树文件
└── align/
    ├── *.msa.fasta            # 多序列比对
    └── *.user_msa.fasta       # 用户序列的比对
```

## 六、常用参数说明

| 参数 | 说明 |
|------|------|
| `--genome_dir` | 输入基因组目录 |
| `--out_dir` | 输出目录 |
| `--extension` | 基因组文件扩展名（fasta、fa、fna 等） |
| `--cpus` | 使用的 CPU 数量 |
| `--min_perc_aa` | 最小蛋白质百分比阈值（默认 10%） |
| `--prot_model` | 蛋白质进化模型（LG 或 WAG，默认 LG） |
| `--force` | 强制覆盖已有结果 |

## 七、性能优化建议

1. **使用多 CPU**：根据你的硬件调整 `--cpus` 参数
2. **批量处理**：一次处理多个基因组会比分别处理更快
3. **中间结果保留**：保留 identify 结果，可以直接进行 classify 而无需重新识别
4. **临时文件清理**：处理完后删除 `intermediate_results` 文件夹节省空间

## 八、常见问题

### Q1: 错误 - "GTDBTK_DATA_PATH not set"
**解决方案**：
```bash
export GTDBTK_DATA_PATH=/mnt/d/3-GTDB-Database/data/release###
```

### Q2: 基因组识别失败（0个标记基因）
**原因**：基因组质量差或者不是原核生物
**解决方案**：检查基因组 completeness 和 contamination，或使用 `--min_perc_aa 5`

### Q3: 运行时间过长
**解决方案**：
- 增加 `--cpus` 参数
- 分开运行 identify 和 classify
- 确保硬盘有足够空间和较好的 I/O 性能

### Q4: 内存不足
**解决方案**：
- 减少 `--cpus` 参数
- 减少同时处理的基因组数量
- 运行 identify 而不是完整流程

## 九、推荐工作流程

```bash
# 1. 数据准备
mkdir -p ~/gtdbtk_work/genomes ~/gtdbtk_work/results
cp /path/to/your/genomes/*.fasta ~/gtdbtk_work/genomes/

# 2. 激活环境
conda activate gtdbtk
export GTDBTK_DATA_PATH=/mnt/d/3-GTDB-Database/data/release###

# 3. 分类
gtdbtk classify_wf \
    --genome_dir ~/gtdbtk_work/genomes \
    --out_dir ~/gtdbtk_work/results \
    --cpus 8 \
    --extension fasta

# 4. 查看结果
cat ~/gtdbtk_work/results/classify/classify_wf.summary.tsv
cat ~/gtdbtk_work/results/classify/classify_wf.markers_summary.tsv
```

## 十、参考资源

- 官方文档：https://ecogenomics.org/gtdbtk/
- GitHub：https://github.com/Ecogenomics/GTDBTk
- GTDB 数据库：https://gtdb.ecogenomics.org/

---
**最后更新**：2025年10月
**适用版本**：GTDB-Tk R220+

#!/bin/bash

#* ============================================================
#* GTDB-Tk 常用工作流程整理说明
#* ============================================================
#
# 工作流概览
# - classify_wf：依据 GTDB 参考系统发育树对基因组进行分类
#   步骤：ani_screening → identify → align → classify
#   含义：先做 ANI 筛选 → 识别标记基因 → 多序列比对 → 分类定位
#
# - de_novo_wf：推断一棵新的系统发育树，并用 GTDB 分类信息标注
#   步骤：identify → align → infer → root → decorate
#   含义：识别标记基因 → 比对 → 构树 → 定根 → 加分类标签
#
# 说明
# - 支持输入：.fna / .fa / .fasta（每个文件 = 一个基因组，由 contigs 组成）
# - 未替换任何变量/路径/选项；仅规范化注释格式。
#
#* ============================================================
#* 1) classify_wf：基因组分类（挂到 GTDB 官方参考树）
#* ------------------------------------------------------------
# 常用参数说明：
# - --genome_dir：包含待分类基因组的目录（每个 FASTA 为一个基因组）
# - --out_dir：输出目录
# - --extension：输入文件后缀（如 fasta / fa / fna）
# - --prefix：结果文件前缀
# - --cpus：总 CPU 数
# - --pplacer_cpus：分配给 pplacer 的 CPU 数（默认 1）
# - --min_af：基因组与参考的最小 alignment fraction（如 0.5）
# - --mash_db：MASH 数据库路径；首次运行会自动创建，可复用
# - --full_tree：使用完整未裁剪的细菌参考树（GTDB-Tk v1 模式，需 >320 GB RAM）
# - --skip_ani_screen：跳过 ANI 筛选（默认否）
# - --force：单个基因组出错时继续
#
# 结果说明：
# - 输出包含多级分类（界-门-纲-目-科-属-种）、进化树文件与统计报告
# - 可直接用于后续系统发育分析或多样性研究
#* ============================================================

gtdbtk classify_wf \
  --genome_dir /home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/input/1-classify/ \
  --out_dir /home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/output/1-classify \
  --min_af 0.5 \
  --extension fasta \
  --prefix classify_Result \
  --cpus 50 \
  --pplacer_cpus 50 \
  --mash_db /data_ssd3/7-luolintao-ssd/MASH_DB

#* ============================================================
#* 1.1) classify（仅在 align 完成后使用的精简分类命令）
#* ------------------------------------------------------------
# usage:
# gtdbtk classify (--genome_dir GENOME_DIR | --batchfile BATCHFILE)
#                 --align_dir ALIGN_DIR --out_dir OUT_DIR
#                 [--skip_ani_screen] [-x EXTENSION] [--prefix PREFIX]
#                 [--cpus CPUS] [--pplacer_cpus PPLACER_CPUS]
#                 [--scratch_dir SCRATCH_DIR] [--genes] [-f]
#                 [--min_af MIN_AF] [--tmpdir TMPDIR] [--debug] [-h]
#
# 示例（保留为注释，不执行）：
# gtdbtk classify \
#   --batchfile /home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/conf/example.tsv \
#   --align_dir /home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/output/2_classify/alignments \
#   --out_dir /home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/output/2_classify \
#   --cpus 50

#* ============================================================
#* 2) de_novo_wf：推断新的系统发育树（含你的基因组 + 参考菌株）
#* ------------------------------------------------------------
# 使用场景：
# - classify_wf：把你的菌“挂”到 GTDB 官方大树上，告诉你它属于哪里
# - de_novo_wf ：自己“长”一棵树，研究自有样本与参考的系统发育关系
#
# 关键参数：
# - --genome_dir <my_genomes>：待分析的基因组目录（FASTA）
# - --outgroup_taxon <outgroup>：外群分类名（如 p__Firmicutes），用于定根
# - --bacteria / --archaea：选择标记基因集
# - --out_dir <output_dir>：输出目录
# - --cpus：CPU 数
# - 其他常见筛选参数（按需）：
#   --min_perc_aa：比对中最低氨基酸覆盖比例（默认 50%）
#   --taxa_filter：仅保留指定门/纲的基因组参与构树
#   --prot_model：蛋白替换模型（WAG 或 LG，默认 WAG）
#* ============================================================

gtdbtk de_novo_wf \
  --genome_dir <my_genomes> \
  --outgroup_taxon <outgroup> \
  --out_dir <output_dir> \
  --cpus 50 \
  --bacteria  # 或 --archaea

#* ============================================================
#* 3) decorate：给树添加分类注释
#* ------------------------------------------------------------
# 可选参数：
# - --gtdbtk_classification_file：指定 GTDB-Tk 分类结果（如 gtdbtk.bac120.summary.tsv），
#   用于给树节点标注分类信息（通常来自 classify 步骤）
# - --custom_taxonomy_file：自定义分类表（可为自有基因组加名或指定外群）
#   格式示例（制表符分隔）：
#   MyGenome1  d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus_casei
#   OutgroupX  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae
# - --tmpdir：临时目录（默认 /tmp）
# - --debug：输出调试文件
#* ============================================================

gtdbtk decorate \
  --input_tree /home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/output/1-classify-wf/classify/classify_Result.backbone.bac120.classify.tree \
  --output_tree /home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/output/1-classify-wf/classify/bac120.decorated.tree \
  --gtdbtk_classification_file /home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/output/1-classify-wf/classify/classify_Result.bac120.summary.tsv

#* ============================================================
#* 4) convert_to_itol：转换为 iTOL 友好格式（可选，用处有限）
#* ------------------------------------------------------------

gtdbtk convert_to_itol \
  --input_tree /home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/output/1-classify-wf/classify/classify_Result.bac120.classify.tree.8.tree \
  --output_tree /home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/output/1-classify-wf/classify/bac120.itol.tree

#* ============================================================
#* 5) root：为树定根（基于指定外群）
#* ------------------------------------------------------------

gtdbtk root \
  --input_tree ./output/archaea.tree \
  --outgroup_taxon p__Nanoarchaeota \
  --output_tree ./output/archaea.rooted.tree

#* ============================================================
#* 6) export_msa：导出参考数据中的未裁剪 MSA（多序列比对）
#* ------------------------------------------------------------
# 功能说明：
# - 用于导出 GTDB-Tk 参考数据库中，未经过修剪的细菌或古菌多序列比对（MSA）。
# - 主要用于检查、比较或自定义后续系统发育分析。
#
# 用法：
# gtdbtk export_msa --domain {arc,bac} --output OUTPUT [--debug] [-h]
#
# 必需参数：
# - --domain：指定要导出的领域（domain）
#     可选项：
#       arc → 古菌（Archaea）
#       bac → 细菌（Bacteria）
# - --output：输出文件路径（导出的 MSA 将保存至此）
#
# 可选参数：
# - --debug：生成调试用中间文件（默认关闭）
#
# 示例：
#* ------------------------------------------------------------
# 输入命令：
# gtdbtk export_msa --domain arc --output /tmp/msa.faa
#
# 输出日志：
# [2020-04-13 10:03:05] INFO: GTDB-Tk v1.1.0
# [2020-04-13 10:03:05] INFO: gtdbtk export_msa --domain arc --output /tmp/msa.faa
# [2020-04-13 10:03:05] INFO: Using GTDB-Tk reference data version r89: /release89
# [2020-04-13 10:03:05] INFO: Done.
#* ============================================================

gtdbtk export_msa \
  --domain arc \
  --output /data_ssd3/7-luolintao-ssd/0-GTDB-Database/GTDB_arc_MSA_aln.faa

#!/usr/bin/env bash
# 作者: BigLin
# 依赖: bash>=4.2, GNU parallel, EMBOSS transeq, coreutils, util-linux

set -euo pipefail  # 遇错退出、未定义变量报错、管道错误传播
IFS=$'\n\t'        # 设置内部字段分隔符为换行与制表符

#=============================
# 参数定义（带中文注释）
#=============================
# | 表号（Table ID） | 名称（NCBI 名称）                                                              | 起始密码子差异           | 使用生物示例        | 备注                  |
# | ------------ | ------------------------------------------------------------------------ | ----------------- | ------------- | ------------------- |
# | **1**        | Standard Code（标准密码）                                                 | ATG (Met)         | 大多数细菌、真核核基因组  | 默认通用表               |
# | **2**        | Vertebrate Mitochondrial Code（脊椎动物线粒体）                           | ATA/ATG → Met     | 人类、小鼠线粒体      | UGA→Trp, AUA→Met    |
# | **3**        | Yeast Mitochondrial Code（酵母线粒体）                                    | ATA/ATG → Met     | 酿酒酵母 mtDNA    | CUA→Thr, AGA/AGG 停止 |
# | **4**        | Mold, Protozoan, Coelenterate Mitochondrial, Mycoplasma/Spiroplasma Code | ATA→Met           | 原虫、霉菌、支原体     | TGA→Trp             |
# | **5**        | Invertebrate Mitochondrial Code（无脊椎线粒体）                            | ATA→Met           | 果蝇、蚊等昆虫 mtDNA | AGA/AGG→Ser         |
# | **6**        | Ciliate, Dasycladacean, Hexamita Nuclear Code                            | TAA/TAG→Gln       | 草履虫等原生生物      | 罕见核编码系统             |
# | **9**        | Echinoderm and Flatworm Mitochondrial Code                               | ATA→Met           | 海胆、扁形动物 mtDNA | AAA→Asn             |
# | **10**       | Euplotid Nuclear Code                                                    | TGA→Cys           | Euplotes 属纤毛虫 | 特殊原生生物核基因组          |
# | **11**       | Bacterial, Archaeal and Plant Plastid Code（细菌、古菌及植物质体代码）      | ATG/GTG/TTG → Met | 细菌、蓝藻、叶绿体     |本脚本使用的表          |
# | **12**       | Alternative Yeast Nuclear Code                                           | CTG→Ser           | Candida 属酵母   | 真菌特有                |
# | **13**       | Ascidian Mitochondrial Code                                              | ATA→Met           | 被囊类（海鞘）mtDNA  | UGA→Trp             |
# | **14**       | Flatworm Mitochondrial Code                                              | ATA→Met           | 吸虫类 mtDNA     | AGA/AGG→Ser         |
# | **16**       | Chlorophycean Mitochondrial Code                                         | TGA→Trp           | 绿藻类 mtDNA     | 植物样 mtDNA           |
# | **21**       | Trematode Mitochondrial Code                                             | ATA→Met           | 吸虫类           | AGA/AGG→Ser         |
# | **22**       | Scenedesmus obliquus Mitochondrial Code                                  | TCA→STOP          | 绿藻            | 极少见                 |
# | **23**       | Thraustochytrium Mitochondrial Code                                      | TTA→STOP          | 海生真菌          | 特殊线粒体密码             |

INPUT_LIST="/data_ssd3/7-luolintao-ssd/0-GTDB-Database/GTDB_ALL.genome.list.txt"   # 输入文件列表
OUTPUT_DIR="/data_ssd3/7-luolintao-ssd/0-GTDB-Database/all_faa"                    # 输出目录
TEMP_DIR="/data_ssd3/7-luolintao-ssd/0-GTDB-Database/temp"                        #? 临时目录（用于解压缩文件）,如果原文件不是压缩文件，直接运行即可
PARALLEL_JOBS="8"                                                                  # 并行任务数
CHECKPOINT_FILE="/data_ssd3/7-luolintao-ssd/0-GTDB-Database/.completed_paths.log"  # 已完成任务记录文件
LOCK_FILE="/data_ssd3/7-luolintao-ssd/0-GTDB-Database/.completed_paths.lock"        # 文件锁，防止并发写入
JOBLOG_FILE="/data_ssd3/7-luolintao-ssd/0-GTDB-Database/.parallel_joblog.tsv"       # GNU parallel 的日志
TRANSEQ_BIN="transeq"                                                              # EMBOSS 的翻译工具
TRANSLATION_TABLE="11"                                                             # 翻译表编号（11为标准细菌遗传密码表）
TMP_SUFFIX=".inprogress"                                                           # 临时文件后缀
FINAL_EXTENSION=".faa"                                                             # 最终输出文件后缀
SIGNAL_CAUGHT="0"                                                                  # 信号捕获标志
COLOR_INFO="\033[1;34m"                                                            # 信息颜色（蓝）
COLOR_WARN="\033[1;33m"                                                            # 警告颜色（黄）
COLOR_ERR="\033[1;31m"                                                             # 错误颜色（红）
COLOR_OK="\033[1;32m"                                                              # 成功颜色（绿）
COLOR_RESET="\033[0m"                                                              # 颜色重置

#=============================
# 日志函数
#=============================

log_info() {
	printf '%b[%s]%b %s\n' "${COLOR_INFO}" "INFO" "${COLOR_RESET}" "$*"
}

log_warn() {
	printf '%b[%s]%b %s\n' "${COLOR_WARN}" "WARN" "${COLOR_RESET}" "$*" >&2
}

log_err() {
	printf '%b[%s]%b %s\n' "${COLOR_ERR}" "ERROR" "${COLOR_RESET}" "$*" >&2
}

log_ok() {
	printf '%b[%s]%b %s\n' "${COLOR_OK}" "DONE" "${COLOR_RESET}" "$*"
}

#=============================
# 清理函数：删除未完成的临时输出
#=============================
cleanup_partial() {
	if [[ -d "${OUTPUT_DIR}" ]]; then
		find "${OUTPUT_DIR}" -type f -name "*${TMP_SUFFIX}" -print0 | xargs -0r rm -f
	fi
	if [[ -d "${TEMP_DIR}" ]]; then
		find "${TEMP_DIR}" -type f -name "*.fna" -print0 | xargs -0r rm -f
	fi
}

#=============================
# 退出时清理逻辑
#=============================
exit_trap() {
	local exit_code
	exit_code="$?"
	trap - EXIT
	cleanup_partial
	if [[ "${exit_code}" -ne 0 ]]; then
		if [[ "${SIGNAL_CAUGHT}" -eq 1 ]]; then
			log_warn "脚本被中断，已清理未完成输出。"
		else
			log_err "脚本执行失败，退出码 ${exit_code}。"
		fi
	else
		log_ok "全部翻译任务已成功完成。"
	fi
}

#=============================
# 信号捕获（如 Ctrl+C）
#=============================
signal_trap() {
	SIGNAL_CAUGHT="1"
	log_warn "检测到终止信号，正在清理。"
	exit 130
}

trap exit_trap EXIT
trap signal_trap SIGINT SIGTERM

#=============================
# 依赖性检查
#=============================
if ! command -v "${TRANSEQ_BIN}" >/dev/null 2>&1; then
	log_err "缺少依赖: EMBOSS transeq。"
	exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
	log_err "缺少依赖: GNU parallel。"
	exit 1
fi

if ! command -v flock >/dev/null 2>&1; then
	log_err "缺少依赖: flock (util-linux)。"
	exit 1
fi

#=============================
# 初始化目录与检查点文件
#=============================
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TEMP_DIR}"
touch "${CHECKPOINT_FILE}"

#=============================
# 读取已完成的任务
#=============================
declare -A COMPLETED_MAP=()
if [[ -s "${CHECKPOINT_FILE}" ]]; then
	while IFS= read -r completed_path; do
		if [[ -n "${completed_path}" ]]; then
			COMPLETED_MAP["${completed_path}"]=1
		fi
	done < "${CHECKPOINT_FILE}"
fi

#=============================
# 生成任务列表
#=============================
declare -a TASKS=()
TOTAL_INPUTS="0"

if [[ ! -f "${INPUT_LIST}" ]]; then
	log_err "未找到输入列表文件: ${INPUT_LIST}"
	exit 1
fi

while IFS= read -r raw_path; do
	if [[ -z "${raw_path}" ]]; then
		continue
	fi
	if [[ "${raw_path}" == \#* ]]; then
		continue
	fi
	input_path="${raw_path%$'\r'}"
	input_path="${input_path#${input_path%%[![:space:]]*}}"
	input_path="${input_path%${input_path##*[![:space:]]}}"
	if [[ -z "${input_path}" ]]; then
		continue
	fi
	TOTAL_INPUTS=$((TOTAL_INPUTS + 1))
	if [[ ! -f "${input_path}" ]]; then
		log_warn "文件不存在，跳过: ${input_path}"
		continue
	fi
	if [[ -n "${COMPLETED_MAP["$input_path"]+x}" ]]; then
		continue
	fi
	output_name="$(basename "${input_path}")${FINAL_EXTENSION}"
	final_path="${OUTPUT_DIR}/${output_name}"
	if [[ -f "${final_path}" ]]; then
		if [[ -z "${COMPLETED_MAP["$input_path"]+x}" ]]; then
			{
				flock -x 200
				printf '%s\n' "${input_path}" >> "${CHECKPOINT_FILE}"
			} 200>"${LOCK_FILE}"
			COMPLETED_MAP["${input_path}"]=1
		fi
		continue
	fi
	TASKS+=("${input_path}|${final_path}")
done < "${INPUT_LIST}"

if [[ "${#TASKS[@]}" -eq 0 ]]; then
	log_info "无待处理任务，所有翻译均已完成。"
	exit 0
fi

log_info "扫描到输入文件总数: ${TOTAL_INPUTS}"
log_info "待翻译任务数: ${#TASKS[@]}"

#=============================
# 导出环境变量供 parallel 使用
#=============================
export OUTPUT_DIR TEMP_DIR TMP_SUFFIX FINAL_EXTENSION CHECKPOINT_FILE LOCK_FILE TRANSEQ_BIN TRANSLATION_TABLE COLOR_OK COLOR_ERR COLOR_RESET

#=============================
# 核心处理函数：执行 transeq 翻译
#=============================
process_entry() {
	local work_item input_path final_path temp_path base_dir actual_input temp_unzip
	work_item="$1"
	input_path="${work_item%%|*}"
	final_path="${work_item#*|}"
	temp_path="${final_path}${TMP_SUFFIX}"
	base_dir="$(dirname "${final_path}")"
	mkdir -p "${base_dir}"
	rm -f "${temp_path}"

	# 处理压缩文件（如果输入是 .gz 格式，先解压到临时目录）
	actual_input="${input_path}"
	temp_unzip=""
	if [[ "${input_path}" == *.gz ]]; then
		local basename_file
		basename_file="$(basename "${input_path}" .gz)"
		temp_unzip="${TEMP_DIR}/${basename_file}"
		if ! gunzip -c "${input_path}" > "${temp_unzip}" 2>/dev/null; then
			rm -f "${temp_unzip}"
			printf '%b[%s]%b %s (解压失败)\n' "${COLOR_ERR}" "FAIL" "${COLOR_RESET}" "${input_path}" >&2
			return 1
		fi
		actual_input="${temp_unzip}"
	fi

	# 执行翻译并检查结果
	if "${TRANSEQ_BIN}" -sequence "${actual_input}" -outseq "${temp_path}" -table "${TRANSLATION_TABLE}" -clean >/dev/null 2>&1; then
		mv "${temp_path}" "${final_path}"
		[[ -n "${temp_unzip}" ]] && rm -f "${temp_unzip}"
		{
			flock -x 200
			printf '%s\n' "${input_path}" >> "${CHECKPOINT_FILE}"
		} 200>"${LOCK_FILE}"
		printf '%b[%s]%b %s\n' "${COLOR_OK}" "OK" "${COLOR_RESET}" "${input_path}"
	else
		rm -f "${temp_path}"
		[[ -n "${temp_unzip}" ]] && rm -f "${temp_unzip}"
		printf '%b[%s]%b %s\n' "${COLOR_ERR}" "FAIL" "${COLOR_RESET}" "${input_path}" >&2
		return 1
	fi
}

export -f process_entry

#=============================
# 并行执行所有任务
#=============================
parallel --jobs "${PARALLEL_JOBS}" --bar --halt soon,fail=1 --joblog "${JOBLOG_FILE}" process_entry ::: "${TASKS[@]}"

log_ok "全部处理完成。"

python3 /home/luolintao/test_mail.py "0-GTDB数据库/script/5-fasta→faa.sh任务完成通知" "<p>0-GTDB数据库/script/5-fasta→faa.sh分析已完成，请查看结果目录。</p>"

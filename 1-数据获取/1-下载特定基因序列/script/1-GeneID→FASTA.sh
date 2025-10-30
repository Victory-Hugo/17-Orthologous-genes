#!/bin/bash
# 从基因ID列表下载蛋白质FASTA序列
# 支持并行处理、NCBI API密钥、断点续跑、重试机制
# 日期：2025年10月5日

unset http_proxy
unset https_proxy

# 默认配置
INPUT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/1-下载特定基因序列/conf/lolA_id.txt" #? 每行一个GeneID
OUTDIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/1-下载特定基因序列/data/"
LOG_FILE="$OUTDIR/gene_download_log.txt"
CSV_FILE="$OUTDIR/gene_download_map.csv"
DEFAULT_PARALLEL_JOBS=8  # 基因查询相对较轻，可以使用更多并发
FORMAT="FAA"  # 下载格式: FAA (蛋白质) 或 FNA (核苷酸)，或 AUTO (优先蛋白质，无则下载核苷酸)

# NCBI API配置
NCBI_API_KEY="29b326d54e7a21fc6c8b9afe7d71f441d809" # 请替换为您自己的API密钥
export NCBI_API_KEY

# 并行任务数和下载格式（可通过命令行参数调整）
PARALLEL_JOBS=${1:-$DEFAULT_PARALLEL_JOBS}
FORMAT=${2:-$FORMAT}

# 验证下载格式
FORMAT=$(echo "$FORMAT" | tr '[:lower:]' '[:upper:]')
if [[ "$FORMAT" != "FAA" ]] && [[ "$FORMAT" != "FNA" ]] && [[ "$FORMAT" != "AUTO" ]]; then
    echo "错误: 下载格式必须是 FAA (蛋白质) / FNA (核苷酸) / AUTO (自动选择)"
    echo "使用方法: $0 [并行任务数] [下载格式]"
    echo "示例: $0 10 FAA    # 使用10个并行任务，只下载蛋白质序列"
    echo "示例: $0 8 FNA     # 使用8个并行任务，只下载核苷酸序列"
    echo "示例: $0 12 AUTO   # 使用12个并行任务，自动选择（优先蛋白质）"
    echo "示例: $0           # 使用默认设置 (并行数=$DEFAULT_PARALLEL_JOBS, 格式=$FORMAT)"
    exit 1
fi

# 验证并行任务数
if [[ ! "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || [ "$PARALLEL_JOBS" -lt 1 ] || [ "$PARALLEL_JOBS" -gt 20 ]; then
    echo "错误: 并行任务数必须是1-20之间的整数"
    echo "使用方法: $0 [并行任务数] [下载格式]"
    echo "示例: $0 10 FAA    # 使用10个并行任务，只下载蛋白质序列"
    echo "示例: $0 8 FNA     # 使用8个并行任务，只下载核苷酸序列"
    echo "示例: $0 12 AUTO   # 使用12个并行任务，自动选择（优先蛋白质）"
    exit 1
fi

mkdir -p "$OUTDIR"

if [[ ! -f "$INPUT" ]]; then
  echo "ERROR: 输入文件不存在: $INPUT" >&2
  exit 1
fi

# 去重、清理输入
mapfile -t GENE_IDS < <(
  sed -e 's/\r$//' "$INPUT" \
    | awk '{ gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0); if($0!="" && $0 !~ /^#/) print $0 }' \
    | sort -u
)

if (( ${#GENE_IDS[@]} == 0 )); then
  echo "WARN: 未从 $INPUT 读取到任何基因ID" >&2
  exit 0
fi

echo "开始处理 ${#GENE_IDS[@]} 个基因ID (使用 $PARALLEL_JOBS 个并行任务, 下载格式: $FORMAT)"

# 初始化日志文件
init_log_file() {
    if [ ! -f "$LOG_FILE" ]; then
        echo "# Gene->Protein FASTA Download Log - Created on $(date)" > "$LOG_FILE"
        echo "# Format: TIMESTAMP | GENE_ID | STATUS | OUTPUT_FILE | FILE_SIZE | SEQUENCE_COUNT | NOTES" >> "$LOG_FILE"
        echo "# STATUS: SUCCESS/FAILED/SKIPPED" >> "$LOG_FILE"
        echo "" >> "$LOG_FILE"
    fi
}

# 记录日志的函数
log_result() {
    local gene_id="$1"
    local status="$2"
    local output_file="$3"
    local file_size="$4"
    local protein_count="$5"
    local notes="$6"
    
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "$timestamp | $gene_id | $status | $output_file | $file_size | $protein_count | $notes" >> "$LOG_FILE"
}

# 检查基因ID是否已经成功处理
is_gene_completed() {
    local gene_id="$1"
    if [ -f "$LOG_FILE" ]; then
        grep -q "| $gene_id | SUCCESS |" "$LOG_FILE"
        return $?
    fi
    return 1
}

# 获取失败的基因ID列表
get_failed_genes() {
    if [ -f "$LOG_FILE" ]; then
        grep "| FAILED |" "$LOG_FILE" | cut -d'|' -f2 | tr -d ' ' | sort -u
    fi
}

# 初始化文件
init_log_file

# 如果 CSV 不存在就初始化表头
if [[ ! -f "$CSV_FILE" ]]; then
  echo "Gene_ID,Output_File,File_Size,Sequence_Count,Type" > "$CSV_FILE"
fi

# 重试函数
retry_esearch() {
  local query="$1"
  local db="$2"
  local max_retries=3
  local retry_count=0
  local result=""
  
  while (( retry_count < max_retries )); do
    result=$(esearch -db "$db" -query "$query" 2>/dev/null)
    if [[ -n "$result" ]]; then
      echo "$result"
      return 0
    fi
    ((retry_count++))
    echo "    重试 $retry_count/$max_retries ..." >&2
    sleep $((retry_count * 2))  # 递增等待时间
  done
  return 1
}

# 检查基因是否有关联的蛋白质
check_gene_proteins() {
  local gene_id="$1"
  local protein_count
  
  protein_count=$(retry_esearch "$gene_id" "gene" | elink -target protein | grep -o '<Count>[0-9]*</Count>' | sed 's/<Count>\|<\/Count>//g' 2>/dev/null)
  
  if [[ -n "$protein_count" ]] && [[ "$protein_count" -gt 0 ]]; then
    echo "$protein_count"
    return 0
  else
    echo "0"
    return 1
  fi
}

# 获取基因信息用于诊断
get_gene_info() {
  local gene_id="$1"
  local gene_info
  
  gene_info=$(retry_esearch "$gene_id" "gene" | efetch -format docsum 2>/dev/null | grep -E '<Name>|<Description>|<ScientificName>' | head -3 | sed 's/<[^>]*>//g' | tr '\n' '|' | sed 's/|$//')
  
  echo "$gene_info"
}

# 安全下载蛋白质FASTA
safe_download_protein_fasta() {
  local gene_id="$1"
  local output_file="$2"
  local max_retries=3
  local retry_count=0
  
  echo "    基因 $gene_id 开始下载蛋白质序列..." >&2
  
  while (( retry_count < max_retries )); do
    # 清空文件以防止部分下载
    > "$output_file"
    
    # 执行下载命令
    if retry_esearch "$gene_id" "gene" | elink -target protein | efetch -format fasta > "$output_file" 2>/dev/null; then
      # 检查文件是否有实际内容（大于 100 字节且包含FASTA格式）
      if [[ -s "$output_file" ]] && (( $(stat -c%s "$output_file") > 100 )) && grep -q "^>" "$output_file"; then
        return 0
      fi
    fi
    
    ((retry_count++))
    echo "    蛋白质FASTA下载重试 $retry_count/$max_retries ..." >&2
    sleep $((retry_count * 3))  # 递增等待时间
  done
  
  # 删除失败的空文件
  rm -f "$output_file"
  return 1
}

# 安全下载核苷酸FASTA
safe_download_nucleotide_fasta() {
  local gene_id="$1"
  local output_file="$2"
  local max_retries=3
  local retry_count=0
  
  echo "    基因 $gene_id 开始下载核苷酸序列..." >&2
  
  while (( retry_count < max_retries )); do
    # 清空文件以防止部分下载
    > "$output_file"
    
    # 执行下载命令 - 直接从基因数据库获取核苷酸序列
    if retry_esearch "$gene_id" "gene" | elink -target nuccore | efetch -format fasta > "$output_file" 2>/dev/null; then
      # 检查文件是否有实际内容（大于 100 字节且包含FASTA格式）
      if [[ -s "$output_file" ]] && (( $(stat -c%s "$output_file") > 100 )) && grep -q "^>" "$output_file"; then
        return 0
      fi
    fi
    
    ((retry_count++))
    echo "    核苷酸FASTA下载重试 $retry_count/$max_retries ..." >&2
    sleep $((retry_count * 3))  # 递增等待时间
  done
  
  # 删除失败的空文件
  rm -f "$output_file"
  return 1
}

# 计算FASTA文件中的蛋白质数量
count_proteins() {
    local fasta_file="$1"
    if [[ -f "$fasta_file" ]]; then
        grep -c "^>" "$fasta_file" 2>/dev/null || echo "0"
    else
        echo "0"
    fi
}

# 单个任务处理函数
process_gene() {
  local gene_id="$1"
  local current="$2"
  local total="$3"

  # 去除可能的空白字符
  gene_id=$(echo "$gene_id" | tr -d '[:space:]')
  
  # 跳过空的gene_id
  if [ -z "$gene_id" ]; then
      echo "跳过空的gene_id"
      return
  fi

  echo "[$current/$total] 处理基因ID: $gene_id (PID: $$)"

  # 检查是否已经成功处理过
  if is_gene_completed "$gene_id"; then
      echo "  [$$] ✓ 基因 $gene_id 已在之前成功处理，跳过"
      log_result "$gene_id" "SKIPPED" "PREVIOUSLY_COMPLETED" "EXISTING" "EXISTING" "Previously completed successfully"
      return
  fi

  # 定义输出文件（默认蛋白质格式）
  local output_file="${OUTDIR}/Gene_ID_${gene_id}.faa"
  local nucleotide_file="${OUTDIR}/Gene_ID_${gene_id}.fna"
  
  # 检查是否已有蛋白质或核苷酸文件存在
  local existing_file=""
  local existing_type=""
  
  if [[ -f "$output_file" ]] && [[ -s "$output_file" ]] && grep -q "^>" "$output_file"; then
      existing_file="$output_file"
      existing_type="蛋白质"
  elif [[ -f "$nucleotide_file" ]] && [[ -s "$nucleotide_file" ]] && grep -q "^>" "$nucleotide_file"; then
      existing_file="$nucleotide_file"
      existing_type="核苷酸"
  fi
  
  if [[ -n "$existing_file" ]]; then
      local file_size=$(stat -c%s "$existing_file")
      local sequence_count=$(count_proteins "$existing_file")
      echo "  [$$] ✓ ${existing_type}文件已存在且有效: $gene_id ($sequence_count 个序列)"
      log_result "$gene_id" "SUCCESS" "$existing_file" "$file_size" "$sequence_count" "File already exists and valid ($existing_type)"
      
      # 更新CSV
      if ! grep -q "^${gene_id}," "$CSV_FILE"; then
          local file_type=$(if [[ "$existing_type" == "蛋白质" ]]; then echo "PROTEIN"; else echo "NUCLEOTIDE"; fi)
          echo "$gene_id,$existing_file,$file_size,$sequence_count,$file_type" >> "$CSV_FILE"
      fi
      return 0
  fi

  echo "  [$$] 正在下载基因 $gene_id 的序列..."
  
  # 根据 FORMAT 参数决定下载策略
  case "$FORMAT" in
      "FAA")
          # 仅下载蛋白质序列
          echo "  [$$] 格式: FAA (仅下载蛋白质序列)"
          if safe_download_protein_fasta "$gene_id" "$output_file"; then
              local file_size=$(stat -c%s "$output_file")
              local sequence_count=$(count_proteins "$output_file")
              
              echo "  [$$] ✓ 蛋白质序列下载成功: $gene_id ($sequence_count 个序列, ${file_size} 字节)"
              log_result "$gene_id" "SUCCESS" "$output_file" "$file_size" "$sequence_count" "Protein sequences downloaded successfully"
              
              # 更新CSV
              if ! grep -q "^${gene_id}," "$CSV_FILE"; then
                  echo "$gene_id,$output_file,$file_size,$sequence_count,PROTEIN" >> "$CSV_FILE"
              fi
          else
              echo "  [$$] ✗ 蛋白质序列下载失败: $gene_id"
              log_result "$gene_id" "FAILED" "" "0" "0" "Protein download failed after retries"
              return 1
          fi
          ;;
          
      "FNA")
          # 仅下载核苷酸序列
          echo "  [$$] 格式: FNA (仅下载核苷酸序列)"
          local nucleotide_output_file="${output_file%.faa}.fna"
          if safe_download_nucleotide_fasta "$gene_id" "$nucleotide_output_file"; then
              local file_size=$(stat -c%s "$nucleotide_output_file")
              local sequence_count=$(count_proteins "$nucleotide_output_file")
              
              echo "  [$$] ✓ 核苷酸序列下载成功: $gene_id ($sequence_count 个序列, ${file_size} 字节)"
              log_result "$gene_id" "SUCCESS" "$nucleotide_output_file" "$file_size" "$sequence_count" "Nucleotide sequences downloaded successfully"
              
              # 更新CSV
              if ! grep -q "^${gene_id}," "$CSV_FILE"; then
                  echo "$gene_id,$nucleotide_output_file,$file_size,$sequence_count,NUCLEOTIDE" >> "$CSV_FILE"
              fi
          else
              echo "  [$$] ✗ 核苷酸序列下载失败: $gene_id"
              log_result "$gene_id" "FAILED" "" "0" "0" "Nucleotide download failed after retries"
              return 1
          fi
          ;;
          
      "AUTO")
          # 自动选择：优先蛋白质，无则下载核苷酸
          echo "  [$$] 格式: AUTO (自动选择: 优先蛋白质，无则下载核苷酸)"
          local protein_count
          protein_count=$(check_gene_proteins "$gene_id")
          
          if [[ "$protein_count" -gt 0 ]]; then
              # 有蛋白质序列，下载蛋白质FASTA
              echo "  [$$] 基因 $gene_id 有 $protein_count 个蛋白质序列，下载蛋白质..."
              
              if safe_download_protein_fasta "$gene_id" "$output_file"; then
                  local file_size=$(stat -c%s "$output_file")
                  local sequence_count=$(count_proteins "$output_file")
                  
                  echo "  [$$] ✓ 蛋白质序列下载成功: $gene_id ($sequence_count 个序列, ${file_size} 字节)"
                  log_result "$gene_id" "SUCCESS" "$output_file" "$file_size" "$sequence_count" "Protein sequences downloaded successfully"
                  
                  # 更新CSV
                  if ! grep -q "^${gene_id}," "$CSV_FILE"; then
                      echo "$gene_id,$output_file,$file_size,$sequence_count,PROTEIN" >> "$CSV_FILE"
                  fi
              else
                  echo "  [$$] ✗ 蛋白质序列下载失败: $gene_id"
                  log_result "$gene_id" "FAILED" "" "0" "0" "Protein download failed after retries"
                  return 1
              fi
          else
              # 没有蛋白质序列，尝试下载核苷酸序列
              echo "  [$$] 基因 $gene_id 无蛋白质序列，尝试下载核苷酸序列..."
              
              # 改变文件扩展名为.fna
              local nucleotide_output_file="${output_file%.faa}.fna"
              
              if safe_download_nucleotide_fasta "$gene_id" "$nucleotide_output_file"; then
                  local file_size=$(stat -c%s "$nucleotide_output_file")
                  local sequence_count=$(count_proteins "$nucleotide_output_file")
                  
                  echo "  [$$] ✓ 核苷酸序列下载成功: $gene_id ($sequence_count 个序列, ${file_size} 字节)"
                  log_result "$gene_id" "SUCCESS" "$nucleotide_output_file" "$file_size" "$sequence_count" "Nucleotide sequences downloaded successfully"
                  
                  # 更新CSV
                  if ! grep -q "^${gene_id}," "$CSV_FILE"; then
                      echo "$gene_id,$nucleotide_output_file,$file_size,$sequence_count,NUCLEOTIDE" >> "$CSV_FILE"
                  fi
              else
                  echo "  [$$] ✗ 核苷酸序列下载也失败: $gene_id"
                  
                  # 获取基因信息用于诊断
                  local gene_info
                  gene_info=$(get_gene_info "$gene_id")
                  echo "  [$$] 基因信息: $gene_info"
                  
                  log_result "$gene_id" "FAILED" "" "0" "0" "Both protein and nucleotide download failed: $gene_info"
                  return 1
              fi
          fi
          ;;
  esac
}

# 导出函数供并行处理使用
export -f process_gene
export -f is_gene_completed
export -f log_result
export -f retry_esearch
export -f safe_download_protein_fasta
export -f safe_download_nucleotide_fasta
export -f count_proteins
export -f check_gene_proteins
export -f get_gene_info
export CSV_FILE OUTDIR LOG_FILE NCBI_API_KEY FORMAT

# 使用并行处理
echo "使用 $PARALLEL_JOBS 个并行任务处理 ${#GENE_IDS[@]} 个基因ID..."
printf '%s\n' "${GENE_IDS[@]}" | nl -nln | xargs -n1 -P"$PARALLEL_JOBS" -I{} bash -c '
    line="{}"
    current=$(echo "$line" | cut -f1)
    gene_id=$(echo "$line" | cut -f2)
    total='"${#GENE_IDS[@]}"'
    process_gene "$gene_id" "$current" "$total"
    sleep 0.1
'

echo ""
echo "全部任务完成！"
echo "FASTA文件目录: $OUTDIR"
echo "CSV对照表: $CSV_FILE"

# 显示统计信息
protein_count=$(find "$OUTDIR" -name "*.faa" 2>/dev/null | wc -l)
nucleotide_count=$(find "$OUTDIR" -name "*.fna" 2>/dev/null | wc -l)
total_files=$((protein_count + nucleotide_count))
total_sequences=0
total_size=0

# 计算总序列数和文件大小（蛋白质文件）
for fasta_file in "$OUTDIR"/*.faa; do
    if [[ -f "$fasta_file" ]]; then
        sequences=$(count_proteins "$fasta_file")
        size=$(stat -c%s "$fasta_file" 2>/dev/null || echo "0")
        total_sequences=$((total_sequences + sequences))
        total_size=$((total_size + size))
    fi
done

# 计算总序列数和文件大小（核苷酸文件）
for fasta_file in "$OUTDIR"/*.fna; do
    if [[ -f "$fasta_file" ]]; then
        sequences=$(count_proteins "$fasta_file")  # 这个函数同样适用于核苷酸序列计数
        size=$(stat -c%s "$fasta_file" 2>/dev/null || echo "0")
        total_sequences=$((total_sequences + sequences))
        total_size=$((total_size + size))
    fi
done

echo ""
echo "统计信息："
echo "- 成功下载的FASTA文件总数: $total_files 个"
echo "  * 蛋白质文件 (.faa): $protein_count 个"
echo "  * 核苷酸文件 (.fna): $nucleotide_count 个"
echo "- 总序列数: $total_sequences 个"
echo "- 总文件大小: $(( total_size / 1024 )) KB"

# 显示日志统计
if [ -f "$LOG_FILE" ]; then
    echo ""
    echo "详细统计（基于日志文件）："
    success_count=$(grep -c "| SUCCESS |" "$LOG_FILE" 2>/dev/null || echo "0")
    failed_count=$(grep -c "| FAILED |" "$LOG_FILE" 2>/dev/null || echo "0")
    skipped_count=$(grep -c "| SKIPPED |" "$LOG_FILE" 2>/dev/null || echo "0")
    
    echo "- 成功下载: $success_count 个"
    echo "- 下载失败: $failed_count 个"
    echo "- 已跳过: $skipped_count 个"
    
    # 显示失败的基因ID
    failed_genes=$(get_failed_genes)
    if [ -n "$failed_genes" ]; then
        echo ""
        echo "需要重新处理的失败基因ID："
        echo "$failed_genes"
        
        # 创建失败基因ID的重试文件
        retry_file="$OUTDIR/retry_gene_ids.txt"
        echo "$failed_genes" > "$retry_file"
        echo ""
        echo "失败基因ID已保存到: $retry_file"
        echo "可以使用以下命令重新处理失败的基因ID："
        echo "  sed -i 's|INPUT=.*|INPUT=\"$retry_file\"|' $0"
        echo "  $0"
    fi
fi

echo ""
echo "断点续跑功能："
echo "- 脚本会记录每个基因ID的处理状态到日志文件"
echo "- 重新运行脚本时会自动跳过已成功处理的基因ID"
echo "- 失败的基因ID会记录详细错误信息，方便后续重试"
echo "- 日志文件位置: $LOG_FILE"
echo ""
echo "脚本支持并行处理和下载格式选择，可通过命令行参数调整："
echo "  $0                  # 使用默认设置 (并行数=$DEFAULT_PARALLEL_JOBS, 格式=$FORMAT)"
echo "  $0 10               # 使用10个并行任务，默认格式"
echo "  $0 1 FAA            # 单线程处理，仅下载蛋白质序列"
echo "  $0 12 FNA           # 12个并行任务，仅下载核苷酸序列"
echo "  $0 15 AUTO          # 15个并行任务，自动选择（优先蛋白质）"
echo ""
echo "下载格式说明："
echo "  FAA   - 仅下载蛋白质序列（.faa 文件）"
echo "  FNA   - 仅下载核苷酸序列（.fna 文件）"
echo "  AUTO  - 自动选择（默认）：优先下载蛋白质序列，无则下载核苷酸序列"
echo ""
echo "API密钥优势："
echo "- 10次/秒的请求限制（vs 普通3次/秒）"
echo "- 支持更高的并行数（8-15个任务）"
echo "- 更稳定的服务质量"
echo ""
echo "注意事项："
echo "- 基因查询相对较轻，可使用较高的并行数（8-15个）"
echo "- 如果网络不稳定，可降低并行数到5-8个"
echo "- 每个基因可能对应多个序列"
echo "- AUTO 模式：优先下载蛋白质序列，无蛋白质时下载核苷酸序列"
echo "- FAA 模式：仅下载蛋白质序列，无蛋白质则处理失败"
echo "- FNA 模式：仅下载核苷酸序列"
echo "- 输出文件命名格式: Gene_ID_[基因ID].faa (蛋白质) 或 Gene_ID_[基因ID].fna (核苷酸)"
echo ""
echo "文件说明："
echo "- .faa文件: 包含该基因对应的蛋白质序列（FORMAT=FAA 或 AUTO 模式下）"
echo "- .fna文件: 包含该基因对应的核苷酸序列（FORMAT=FNA 模式下，或 AUTO 模式下无蛋白质时）"
echo "- CSV对照表: 记录基因ID、文件路径、文件大小、序列数量、类型等信息"
echo "- 日志文件: 记录详细的下载过程和错误信息"
echo ""
echo "如需修改输入文件路径，请编辑脚本开头的 INPUT 变量"
echo "如需修改输出目录，请编辑脚本开头的 OUTDIR 变量"

#!/bin/bash
# 从 GCF/GCA accession 获取 Assembly 的 BioSample (SAMN) 并下载 XML
# 自动补全标准 Accession（带版本号）
# 支持并行处理、API密钥、断点续跑
# 输出 XML 文件 + CSV 对照表（Input_ID ↔ Standard_ID ↔ SAMN ↔ Organism）
# 作者：Lintao Luo  
# 日期：2025年8月19日

unset http_proxy
unset https_proxy

# 默认配置
INFILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/conf/balanced_accessions_ID.txt"
XML_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/meta"
CSV_FILE="$XML_DIR/assembly_biosample_map.csv"
LOG_FILE="$XML_DIR/assembly_biosample_log.txt"
DEFAULT_PARALLEL_JOBS=5  # 默认并行任务数（Assembly查询较重，用较少并发）

# NCBI API配置
NCBI_API_KEY="29b326d54e7a21fc6c8b9afe7d71f441d809" #!请自己在NCBI申请API密钥
export NCBI_API_KEY

# 并行任务数（可通过命令行参数调整）
PARALLEL_JOBS=${1:-$DEFAULT_PARALLEL_JOBS}

# 验证并行任务数
if [[ ! "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || [ "$PARALLEL_JOBS" -lt 1 ] || [ "$PARALLEL_JOBS" -gt 15 ]; then
    echo "错误: 并行任务数必须是1-15之间的整数"
    echo "使用方法: $0 [并行任务数]"
    echo "示例: $0 6  # 使用6个并行任务"
    exit 1
fi

mkdir -p "$XML_DIR"

if [[ ! -f "$INFILE" ]]; then
  echo "ERROR: 输入文件不存在: $INFILE" >&2
  exit 1
fi

# 去重、清理输入
mapfile -t ASMS < <(
  sed -e 's/\r$//' "$INFILE" \
    | awk '{ gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0); if($0!="" && $0 !~ /^#/) print $0 }' \
    | sort -u
)

if (( ${#ASMS[@]} == 0 )); then
  echo "WARN: 未从 $INFILE 读取到任何 Assembly ID" >&2
  exit 0
fi

echo "开始处理 ${#ASMS[@]} 个 Assembly ID (使用 $PARALLEL_JOBS 个并行任务)"

# 初始化日志文件
init_log_file() {
    if [ ! -f "$LOG_FILE" ]; then
        echo "# Assembly->BioSample Download Log - Created on $(date)" > "$LOG_FILE"
        echo "# Format: TIMESTAMP | INPUT_ID | STATUS | STANDARD_ID | BIOSAMPLE_ID | XML_FILE | NOTES" >> "$LOG_FILE"
        echo "# STATUS: SUCCESS/FAILED/SKIPPED" >> "$LOG_FILE"
        echo "" >> "$LOG_FILE"
    fi
}

# 记录日志的函数
log_result() {
    local input_id="$1"
    local status="$2"
    local standard_id="$3"
    local biosample_id="$4"
    local xml_file="$5"
    local notes="$6"
    
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "$timestamp | $input_id | $status | $standard_id | $biosample_id | $xml_file | $notes" >> "$LOG_FILE"
}

# 检查Assembly是否已经成功处理
is_assembly_completed() {
    local input_id="$1"
    if [ -f "$LOG_FILE" ]; then
        grep -q "| $input_id | SUCCESS |" "$LOG_FILE"
        return $?
    fi
    return 1
}

# 获取失败的Assembly列表
get_failed_assemblies() {
    if [ -f "$LOG_FILE" ]; then
        grep "| FAILED |" "$LOG_FILE" | cut -d'|' -f2 | tr -d ' ' | sort -u
    fi
}

# 初始化文件
init_log_file

# 如果 CSV 不存在就初始化表头
if [[ ! -f "$CSV_FILE" ]]; then
  echo "Input_ID,Standard_ID,BioSample,Organism" > "$CSV_FILE"
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

# 安全下载 XML - 使用 BioSample ID 下载，而不是 Assembly ID
safe_download_xml() {
  local samn_id="$1"
  local output_file="$2"
  local max_retries=3
  local retry_count=0
  
  if [[ -z "$samn_id" ]] || [[ "$samn_id" == "SAMN" ]]; then
    echo "    跳过下载：无效的 BioSample ID" >&2
    return 1
  fi
  
  while (( retry_count < max_retries )); do
    # 清空文件以防止部分下载
    > "$output_file"
    
    # 使用 BioSample ID 下载 XML，而不是 Assembly ID
    if esearch -db biosample -query "$samn_id" | efetch -format xml > "$output_file" 2>/dev/null; then
      # 检查文件是否有实际内容（大于 50 字节）
      if [[ -s "$output_file" ]] && (( $(stat -c%s "$output_file") > 50 )); then
        return 0
      fi
    fi
    
    ((retry_count++))
    echo "    XML下载重试 $retry_count/$max_retries ..." >&2
    sleep $((retry_count * 3))  # 递增等待时间
  done
  
  # 删除失败的空文件
  rm -f "$output_file"
  return 1
}

# 单个任务
process_asm() {
  local input_id="$1"
  local current="$2"
  local total="$3"

  # 去除可能的空白字符
  input_id=$(echo "$input_id" | tr -d '[:space:]')
  
  # 跳过空的assembly_id
  if [ -z "$input_id" ]; then
      echo "跳过空的assembly_id"
      return
  fi

  echo "[$current/$total] 处理 Input: $input_id (PID: $$)"

  # 检查是否已经成功处理过
  if is_assembly_completed "$input_id"; then
      echo "  [$$] ✓ Assembly $input_id 已在之前成功处理，跳过"
      log_result "$input_id" "SKIPPED" "PREVIOUSLY_COMPLETED" "EXISTING" "EXISTING" "Previously completed successfully"
      return
  fi

  # 获取标准化 Accession（带版本号）
  echo "  [$$] 正在获取标准化 Accession..."
  local full_acc
  full_acc=$(retry_esearch "$input_id" "assembly" | efetch -format docsum 2>/dev/null | \
             xtract -pattern DocumentSummary -element AssemblyAccession | head -n1)

  if [[ -z "$full_acc" ]]; then
    echo "  [$$] ✗ 错误: 未找到标准 Accession for $input_id"
    log_result "$input_id" "FAILED" "NOT_FOUND" "" "" "Standard accession not found"
    if ! grep -q "^${input_id}," "$CSV_FILE"; then
      echo "$input_id,,," >> "$CSV_FILE"
    fi
    return 1
  fi

  echo "  [$$] ✓ 标准化 Accession: $full_acc"

  # 从 docsum 提取信息（使用正确的字段名）
  echo "  [$$] 正在获取 BioSample 信息..."
  local info
  info=$(retry_esearch "$full_acc" "assembly" | efetch -format docsum 2>/dev/null | \
         xtract -pattern DocumentSummary -element AssemblyAccession BioSampleAccn Organism | head -n1)

  local std_acc samn org
  std_acc=$(echo "$info" | cut -f1)
  samn=$(echo "$info" | cut -f2)
  org=$(echo "$info" | cut -f3-)

  if [[ -z "$samn" ]]; then
    echo "  [$$] ⚠ 警告: 未找到 BioSample for $full_acc"
    samn=""
  else
    echo "  [$$] ✓ 找到 BioSample: $samn"
  fi

  # 如果 CSV 已有该条目，跳过
  if grep -q "^${input_id},${full_acc}," "$CSV_FILE"; then
    echo "  [$$] ✓ 跳过（CSV已存在）: $input_id → $full_acc → $samn"
    log_result "$input_id" "SUCCESS" "$full_acc" "$samn" "EXISTING" "Previously processed in CSV"
    return 0
  fi

  # 下载 BioSample XML（如果有 SAMN ID）
  local xml_file=""
  if [[ -n "$samn" ]]; then
    xml_file="$XML_DIR/${samn}.xml"
    
    if [[ ! -f "$xml_file" ]] || [[ ! -s "$xml_file" ]] || (( $(stat -c%s "$xml_file") <= 50 )); then
      echo "  [$$] 正在下载 BioSample XML: $samn"
      if safe_download_xml "$samn" "$xml_file"; then
        echo "  [$$] ✓ XML下载成功: $xml_file"
      else
        echo "  [$$] ✗ 警告: BioSample XML下载失败 $samn"
        log_result "$input_id" "FAILED" "$std_acc" "$samn" "" "XML download failed"
        return 1
      fi
    else
      echo "  [$$] ✓ BioSample XML已存在: $samn"
    fi
  fi

  # 写入 CSV（避免重复）
  if ! grep -q "^${input_id},${full_acc}," "$CSV_FILE"; then
    echo "$input_id,$std_acc,$samn,\"$org\"" >> "$CSV_FILE"
  fi

  echo "  [$$] ✓ 完成: $input_id → $std_acc → $samn"
  log_result "$input_id" "SUCCESS" "$std_acc" "$samn" "$xml_file" "Completed successfully"
}

# 导出函数供并行处理使用
export -f process_asm
export -f is_assembly_completed
export -f log_result
export -f retry_esearch
export -f safe_download_xml
export CSV_FILE XML_DIR LOG_FILE NCBI_API_KEY

# 使用并行处理
echo "使用 $PARALLEL_JOBS 个并行任务处理 ${#ASMS[@]} 个 Assembly ID..."
printf '%s\n' "${ASMS[@]}" | nl -nln | xargs -n1 -P"$PARALLEL_JOBS" -I{} bash -c '
    line="{}"
    current=$(echo "$line" | cut -f1)
    asm=$(echo "$line" | cut -f2)
    total='"${#ASMS[@]}"'
    process_asm "$asm" "$current" "$total"
    sleep 0.2
'

echo ""
echo "全部任务完成！"
echo "CSV 对照表: $CSV_FILE"
echo "XML 文件目录: $XML_DIR"

# 显示统计信息
xml_count=$(find "$XML_DIR" -name "*.xml" 2>/dev/null | wc -l)
echo ""
echo "统计信息："
echo "- 成功获取的BioSample XML: $xml_count 个"

# 显示日志统计
if [ -f "$LOG_FILE" ]; then
    echo ""
    echo "详细统计（基于日志文件）："
    success_count=$(grep -c "| SUCCESS |" "$LOG_FILE" 2>/dev/null || echo "0")
    failed_count=$(grep -c "| FAILED |" "$LOG_FILE" 2>/dev/null || echo "0")
    skipped_count=$(grep -c "| SKIPPED |" "$LOG_FILE" 2>/dev/null || echo "0")
    
    echo "- 成功处理: $success_count 个"
    echo "- 处理失败: $failed_count 个"
    echo "- 已跳过: $skipped_count 个"
    
    # 显示失败的Assembly
    failed_assemblies=$(get_failed_assemblies)
    if [ -n "$failed_assemblies" ]; then
        echo ""
        echo "需要重新处理的失败Assembly："
        echo "$failed_assemblies"
        
        # 创建失败Assembly的重试文件
        retry_file="$XML_DIR/retry_assemblies.txt"
        echo "$failed_assemblies" > "$retry_file"
        echo ""
        echo "失败Assembly已保存到: $retry_file"
        echo "可以使用以下命令重新处理失败的Assembly："
        echo "  sed -i 's|INFILE=.*|INFILE=\"$retry_file\"|' $0"
        echo "  $0"
    fi
fi

echo ""
echo "断点续跑功能："
echo "- 脚本会记录每个Assembly的处理状态到日志文件"
echo "- 重新运行脚本时会自动跳过已成功处理的Assembly"
echo "- 失败的Assembly会记录详细错误信息，方便后续重试"
echo "- 日志文件位置: $LOG_FILE"
echo ""
echo "脚本支持并行处理，可通过命令行参数调整并行数："
echo "  $0           # 使用默认并行数 ($DEFAULT_PARALLEL_JOBS) - API密钥优化"
echo "  $0 6         # 使用6个并行任务"
echo "  $0 1         # 单线程处理（最安全，但最慢）"
echo "  $0 10        # 使用10个并行任务（API密钥支持高并发）"
echo ""
echo "API密钥优势："
echo "- 10次/秒的请求限制（vs 普通3次/秒）"
echo "- 支持更高的并行数（5-10个任务）"
echo "- 更稳定的服务质量"
echo ""
echo "注意事项："
echo "- Assembly查询相对较重，建议使用适中的并行数（6-10个）"
echo "- 如果网络不稳定，可降低并行数到3-5个"

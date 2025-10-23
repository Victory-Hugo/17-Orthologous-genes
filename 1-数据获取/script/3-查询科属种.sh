#!/bin/bash
# 从物种名称列表查询科属种分类信息
# 支持并行处理、NCBI API密钥、断点续跑、重试机制
# 日期：2025年10月6日
# 检查物种是否已经成功处理
is_species_completed() {
    local species_name="$1"
    if [ -f "$LOG_FILE" ]; then
        grep -q "| $species_name | SUCCESS |" "$LOG_FILE"
        return $?
    fi
    return 1
}

# 检查记录是否已经完成处理（基于完整输入行）
is_record_completed() {
    local input_line="$1"
    IFS=$'\t' read -r gene_id description species_name <<< "$input_line"
    if [ -f "$LOG_FILE" ]; then
        grep -q "| $species_name | SUCCESS |.*ID: $gene_id" "$LOG_FILE"
        return $?
    fi
    return 1
}

# 检查输出文件中是否已存在该记录（基于基因ID）
is_record_in_output_file() {
    local gene_id="$1"
    if [ -f "$OUT" ]; then
        grep -q "^$gene_id\t" "$OUT"
        return $?
    fi
    return 1
}

# 检查输出文件中是否已存在该物种（保持向后兼容）
is_in_output_file() {
    local species_name="$1"
    if [ -f "$OUT" ]; then
        grep -q "^$species_name\t" "$OUT"
        return $?
    fi
    return 1
}

unset http_proxy
unset https_proxy

# 默认配置
IN="/mnt/c/Users/Administrator/Desktop/10.1.5.40/202510221601/TamB_merge.ID.txt"
OUT="/mnt/c/Users/Administrator/Desktop/10.1.5.40/202510221601/TamBtaxonomy.tsv"
LOG_FILE="/mnt/c/Users/Administrator/Desktop/10.1.5.40/202510221601/TAMB_query_log.txt"
DEFAULT_PARALLEL_JOBS=5  # 分类查询相对较重，使用适中的并发数

# NCBI API配置
NCBI_API_KEY="29b326d54e7a21fc6c8b9afe7d71f441d809" # 请替换为您自己的API密钥
export NCBI_API_KEY

# 并行任务数（可通过命令行参数调整）
PARALLEL_JOBS=${1:-$DEFAULT_PARALLEL_JOBS}

# 验证并行任务数
if [[ ! "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || [ "$PARALLEL_JOBS" -lt 1 ] || [ "$PARALLEL_JOBS" -gt 15 ]; then
    echo "错误: 并行任务数必须是1-15之间的整数"
    echo "使用方法: $0 [并行任务数]"
    echo "示例: $0 8  # 使用8个并行任务"
    exit 1
fi

mkdir -p "$(dirname "$OUT")"

if [[ ! -f "$IN" ]]; then
  echo "ERROR: 输入文件不存在: $IN" >&2
  exit 1
fi

# 读取完整的输入文件内容，保留所有原始信息
mapfile -t INPUT_LINES < <(
  cat "$IN" | tr -d '\r' | sed '/^\s*$/d' | \
  awk '{ gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0); if($0!="" && $0 !~ /^#/) print $0 }'
)

if (( ${#INPUT_LINES[@]} == 0 )); then
  echo "WARN: 未从 $IN 读取到任何有效数据行" >&2
  exit 0
fi

echo "开始处理 ${#INPUT_LINES[@]} 行数据 (使用 $PARALLEL_JOBS 个并行任务)"

# 初始化输出文件（只在第一次运行时）
if [ ! -f "$OUT" ]; then
    printf "id\tdescription\tspecies_name\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n" > "$OUT"
fi

# 初始化日志文件
init_log_file() {
    if [ ! -f "$LOG_FILE" ]; then
        echo "# Taxonomy Query Log - Created on $(date)" > "$LOG_FILE"
        echo "# Format: TIMESTAMP | SPECIES_NAME | STATUS | KINGDOM | PHYLUM | CLASS | ORDER | FAMILY | GENUS | SPECIES | NOTES" >> "$LOG_FILE"
        echo "# STATUS: SUCCESS/FAILED/SKIPPED" >> "$LOG_FILE"
        echo "" >> "$LOG_FILE"
    fi
}

# 记录日志的函数
log_result() {
    local species_name="$1"
    local status="$2"
    local kingdom="$3"
    local phylum="$4"
    local class="$5"
    local order="$6"
    local family="$7"
    local genus="$8"
    local species="$9"
    local notes="${10}"
    
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "$timestamp | $species_name | $status | $kingdom | $phylum | $class | $order | $family | $genus | $species | $notes" >> "$LOG_FILE"
}

# 检查物种是否已经成功处理
is_species_completed() {
    local species_name="$1"
    if [ -f "$LOG_FILE" ]; then
        grep -q "| $species_name | SUCCESS |" "$LOG_FILE"
        return $?
    fi
    return 1
}

# 检查输出文件中是否已存在该物种
is_in_output_file() {
    local species_name="$1"
    if [ -f "$OUT" ]; then
        grep -q "^$species_name	" "$OUT"
        return $?
    fi
    return 1
}

# 获取失败的物种列表
get_failed_species() {
    if [ -f "$LOG_FILE" ]; then
        grep "| FAILED |" "$LOG_FILE" | cut -d'|' -f2 | tr -d ' ' | sort -u
    fi
}

# 初始化文件
init_log_file

# 重试函数
retry_esearch() {
  local query="$1"
  local db="$2"
  local max_retries=3
  local retry_count=0
  local result=""
  
  while (( retry_count < max_retries )); do
    result=$(esearch -db "$db" -query "$query" 2>/dev/null || true)
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

# 安全获取分类信息
safe_get_taxonomy() {
  local species_name="$1"
  local max_retries=3
  local retry_count=0
  
  while (( retry_count < max_retries )); do
    # 获取TaxID
    local txid=$(retry_esearch "$species_name" "taxonomy" | efetch -format uid | head -n 1 2>/dev/null || true)
    
    if [[ -n "$txid" ]] && [[ "$txid" =~ ^[0-9]+$ ]]; then
      # 获取分类信息：既从LineageEx取各级分类，也读取当前节点获取species
      local fields=$(
        efetch -db taxonomy -id "$txid" -format xml 2>/dev/null \
        | xtract -pattern Taxon -element Rank,ScientificName -block LineageEx/Taxon -element Rank,ScientificName \
        | awk -F'\t' '{
            # 处理当前节点（前两个字段）
            current_rank=$1; current_name=$2;
            if(current_rank=="species") species_name=current_name;
            
            # 处理LineageEx中的祖先节点（从第3个字段开始）
            for(i=3; i<=NF; i+=2) {
              if($i=="kingdom") k=$(i+1);
              if($i=="phylum") p=$(i+1);
              if($i=="class") c=$(i+1);
              if($i=="order") o=$(i+1);
              if($i=="family") f=$(i+1);
              if($i=="genus") g=$(i+1);
              if($i=="species") s=$(i+1)  # 虽然通常不在LineageEx中，但以防万一
            }
            
            # 如果当前节点是species级别，优先使用当前节点的名称
            if(species_name!="") s=species_name;
          } END { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", (k?k:""), (p?p:""), (c?c:""), (o?o:""), (f?f:""), (g?g:""), (s?s:"") }' 2>/dev/null || true
      )
      
      # 检查是否成功获取到分类信息（至少有kingdom、phylum、class、order、family或genus中的一个）
      if [[ -n "$fields" ]] && [[ "$fields" != $'\t\t\t\t\t\t\t' ]]; then
        echo "$fields"
        return 0
      fi
    fi
    
    ((retry_count++))
    echo "    分类查询重试 $retry_count/$max_retries (物种: $species_name) ..." >&2
    sleep $((retry_count * 3))  # 递增等待时间
  done
  
  return 1
}

# 单个任务处理函数
process_species() {
  local input_line="$1"
  local current="$2"
  local total="$3"

  # 解析输入行的三个字段：ID、描述、物种名称
  IFS=$'\t' read -r gene_id description species_name <<< "$input_line"
  
  # 去除可能的空白字符
  gene_id=$(echo "$gene_id" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
  description=$(echo "$description" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
  species_name=$(echo "$species_name" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
  
  # 清理物种名称中的特殊字符（NCBI查询用）
  # 移除方括号 [] 和多余的空格
  species_query=$(echo "$species_name" | sed 's/\[//g; s/\]//g; s/^[[:space:]]*//; s/[[:space:]]*$//')
  
  # 跳过空的记录
  if [ -z "$gene_id" ] || [ -z "$species_name" ]; then
      echo "跳过空的记录行"
      return
  fi

  echo "[$current/$total] 处理记录: ID=$gene_id, 物种=$species_name (PID: $$)"

  # 检查是否已经成功处理过（基于完整的输入行）
  if is_record_completed "$input_line" || is_record_in_output_file "$gene_id"; then
      echo "  [$$] ✓ 记录 $gene_id 已在之前成功处理，跳过"
      log_result "$species_name" "SKIPPED" "EXISTING" "EXISTING" "EXISTING" "EXISTING" "EXISTING" "EXISTING" "EXISTING" "Previously completed successfully"
      return
  fi

  echo "  [$$] 正在查询物种 $species_name 的分类信息..."
  
  # 查询分类信息（使用清理后的物种名称）
  if taxonomy_info=$(safe_get_taxonomy "$species_query"); then
      echo "  [$$] ✓ 查询成功: $species_name"
      
      # 解析结果
      IFS=$'\t' read -r kingdom phylum class order family genus species <<< "$taxonomy_info"
      
      # 写入输出文件（使用锁机制避免并发写入冲突）
      {
          flock -x 200
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$gene_id" "$description" "$species_name" "$kingdom" "$phylum" "$class" "$order" "$family" "$genus" "$species" >> "$OUT"
      } 200>>"$OUT.lock"
      
      log_result "$species_name" "SUCCESS" "$kingdom" "$phylum" "$class" "$order" "$family" "$genus" "$species" "Query completed successfully for ID: $gene_id"
  else
      echo "  [$$] ✗ 查询失败: $species_name"
      
      # 写入空结果到输出文件
      {
          flock -x 200
          printf "%s\t%s\t%s\t\t\t\t\t\t\t\n" "$gene_id" "$description" "$species_name" >> "$OUT"
      } 200>>"$OUT.lock"
      
      log_result "$species_name" "FAILED" "" "" "" "" "" "" "" "Query failed after retries for ID: $gene_id"
      return 1
  fi
  
  # 添加延迟以避免过于频繁的请求
  sleep 0.2
}

# 导出函数供并行处理使用
export -f process_species
export -f is_species_completed
export -f is_record_completed
export -f is_record_in_output_file
export -f is_in_output_file
export -f log_result
export -f retry_esearch
export -f safe_get_taxonomy
export OUT LOG_FILE NCBI_API_KEY

# 使用并行处理
echo "使用 $PARALLEL_JOBS 个并行任务处理 ${#INPUT_LINES[@]} 行数据..."
{
    i=0
    for input_line in "${INPUT_LINES[@]}"; do
        ((i++))
        echo -e "$i\t$input_line"
    done
} | xargs -P"$PARALLEL_JOBS" -I{} bash -c '
    line="{}"
    current=$(echo "$line" | cut -f1)
    input_line=$(echo "$line" | cut -f2-)
    total='"${#INPUT_LINES[@]}"'
    process_species "$input_line" "$current" "$total"
    sleep 0.1
'

echo ""
echo "全部任务完成！"
echo "分类信息文件: $OUT"

# 显示统计信息
if [ -f "$OUT" ]; then
    total_records=$(tail -n +2 "$OUT" | wc -l)  # 跳过标题行
    success_records=$(tail -n +2 "$OUT" | awk -F'\t' '$4!="" || $5!="" || $6!="" || $7!="" || $8!="" || $9!="" || $10!=""' | wc -l)
    failed_records=$((total_records - success_records))
    
    echo ""
    echo "统计信息："
    echo "- 总处理记录数: $total_records 个"
    echo "- 成功获取分类信息: $success_records 个"
    echo "- 未找到分类信息: $failed_records 个"
fi

# 显示日志统计
if [ -f "$LOG_FILE" ]; then
    echo ""
    echo "详细统计（基于日志文件）："
    success_count=$(grep -c "| SUCCESS |" "$LOG_FILE" 2>/dev/null || echo "0")
    failed_count=$(grep -c "| FAILED |" "$LOG_FILE" 2>/dev/null || echo "0")
    skipped_count=$(grep -c "| SKIPPED |" "$LOG_FILE" 2>/dev/null || echo "0")
    
    echo "- 成功查询: $success_count 个"
    echo "- 查询失败: $failed_count 个"
    echo "- 已跳过: $skipped_count 个"
    
    # 显示失败的物种
    failed_species=$(get_failed_species)
    if [ -n "$failed_species" ]; then
        echo ""
        echo "需要重新处理的失败物种："
        echo "$failed_species"
        
        # 创建失败物种的重试文件
        retry_file="$(dirname "$OUT")/retry_species.txt"
        echo "$failed_species" > "$retry_file"
        echo ""
        echo "失败物种已保存到: $retry_file"
    fi
fi

# 清理锁文件
rm -f "$OUT.lock"

echo ""
echo "断点续跑功能："
echo "- 脚本会记录每个物种的处理状态到日志文件"
echo "- 重新运行脚本时会自动跳过已成功处理的物种"
echo "- 失败的物种会记录详细错误信息，方便后续重试"
echo "- 日志文件位置: $LOG_FILE"
echo ""
echo "脚本支持并行处理，可通过命令行参数调整并行数："
echo "  $0           # 使用默认并行数 ($DEFAULT_PARALLEL_JOBS) - 分类查询优化"
echo "  $0 8         # 使用8个并行任务"
echo "  $0 1         # 单线程处理（最安全，但最慢）"
echo "  $0 10        # 使用10个并行任务（API密钥支持高并发）"
echo ""
echo "API密钥优势："
echo "- 10次/秒的请求限制（vs 普通3次/秒）"
echo "- 支持更高的并行数（5-10个任务）"
echo "- 更稳定的服务质量"
echo ""
echo "注意事项："
echo "- 分类查询相对较重，推荐使用5-8个并行任务"
echo "- 如果网络不稳定，可降低并行数到3-5个"
echo "- 某些物种名称可能在NCBI分类数据库中不存在"
echo "- 输出格式: ID\t描述\t物种名称\t界\t门\t纲\t目\t科\t属\t种"
echo "- 输出文件包含所有输入文件的原始信息，不进行去重处理"

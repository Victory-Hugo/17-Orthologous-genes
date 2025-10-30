#!/usr/bin/env python3
# coding: utf-8
"""
批量下载 NCBI Datasets 数据（已配置默认邮箱和 API Key）
支持断点续跑、彩色输出、进度显示和清理未完成文件

使用方法示例：
  export NCBI_API_KEY="..."      # 可选：覆盖脚本内默认
  export NCBI_EMAIL="..."        # 可选：覆盖脚本内默认
  python3 download_datasets.py -b /mnt/c/Users/Administrator/Desktop/ -f /path/to/下载NCBI.txt -w 5 --resume
"""
import os
import argparse
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import time
import signal
import sys
import glob
import json
from datetime import datetime
from urllib.parse import urlparse, parse_qs, urlencode, urlunparse
from pathlib import Path

# 彩色输出支持
try:
    from colorama import init, Fore, Style, Back
    init(autoreset=True)
    HAS_COLORAMA = True
except ImportError:
    HAS_COLORAMA = False
    
# 进度条支持
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

# -------------------- 默认配置（可通过环境变量或命令行覆盖） --------------------
DEFAULT_NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "giantlinlinlin@gmail.com")
DEFAULT_NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "29b326d54e7a21fc6c8b9afe7d71f441d809")
# -------------------------------------------------------------------------------

# 全局速率控制（会在 main 中根据是否有 key 调整 allowed_rps）
RATE_LOCK = threading.Lock()
LAST_REQUEST_TIME = 0.0
MIN_INTERVAL = 0.0  # 每次请求最小间隔秒，之后在 main() 中计算

# 全局变量用于清理
TEMP_FILES = set()
TEMP_FILES_LOCK = threading.Lock()
SHUTDOWN_EVENT = threading.Event()

# 彩色输出函数
def colored_print(message, color="white", style="normal"):
    """打印彩色消息"""
    if not HAS_COLORAMA:
        print(message)
        return
    
    color_map = {
        "red": Fore.RED,
        "green": Fore.GREEN,
        "yellow": Fore.YELLOW,
        "blue": Fore.BLUE,
        "magenta": Fore.MAGENTA,
        "cyan": Fore.CYAN,
        "white": Fore.WHITE,
    }
    
    style_map = {
        "normal": Style.NORMAL,
        "bright": Style.BRIGHT,
        "dim": Style.DIM,
    }
    
    color_code = color_map.get(color, Fore.WHITE)
    style_code = style_map.get(style, Style.NORMAL)
    
    print(f"{style_code}{color_code}{message}{Style.RESET_ALL}")

def print_success(message):
    """打印成功消息"""
    colored_print(f"✓ {message}", "green", "bright")

def print_error(message):
    """打印错误消息"""
    colored_print(f"✗ {message}", "red", "bright")

def print_warning(message):
    """打印警告消息"""
    colored_print(f"⚠ {message}", "yellow", "bright")

def print_info(message):
    """打印信息消息"""
    colored_print(f"ℹ {message}", "blue", "bright")

def print_progress(message):
    """打印进度消息"""
    colored_print(f"⚡ {message}", "cyan", "bright")

# 临时文件管理
def add_temp_file(filepath):
    """添加临时文件到跟踪列表"""
    with TEMP_FILES_LOCK:
        TEMP_FILES.add(filepath)

def remove_temp_file(filepath):
    """从跟踪列表移除临时文件"""
    with TEMP_FILES_LOCK:
        TEMP_FILES.discard(filepath)

def cleanup_temp_files():
    """清理所有临时文件"""
    with TEMP_FILES_LOCK:
        if TEMP_FILES:
            print_warning(f"正在清理 {len(TEMP_FILES)} 个未完成的下载文件...")
            for temp_file in TEMP_FILES.copy():
                try:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
                        print_info(f"已删除: {temp_file}")
                except Exception as e:
                    print_error(f"删除文件失败 {temp_file}: {e}")
            TEMP_FILES.clear()
            print_success("临时文件清理完成")

def signal_handler(signum, frame):
    """信号处理器，用于优雅关闭"""
    print_warning("\n收到中断信号，正在优雅关闭...")
    SHUTDOWN_EVENT.set()
    cleanup_temp_files()
    print_info("程序已终止")
    sys.exit(0)

# 断点续跑相关函数
def load_progress_state(base_path):
    """加载之前的下载进度"""
    state_file = os.path.join(base_path, "download_progress.json")
    if os.path.exists(state_file):
        try:
            with open(state_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            print_warning(f"无法加载进度文件: {e}")
    return {"completed": set(), "failed": set()}

def save_progress_state(base_path, completed, failed):
    """保存下载进度"""
    state_file = os.path.join(base_path, "download_progress.json")
    state = {
        "completed": list(completed),
        "failed": list(failed),
        "last_update": datetime.now().isoformat()
    }
    try:
        with open(state_file, 'w', encoding='utf-8') as f:
            json.dump(state, f, indent=2, ensure_ascii=False)
    except Exception as e:
        print_error(f"保存进度失败: {e}")

def check_existing_files(base_path, accession_numbers):
    """检查已存在的文件，返回已完成的accession列表"""
    completed = set()
    for acc in accession_numbers:
        zip_file = os.path.join(base_path, f"{acc}_downloaded.zip")
        if os.path.exists(zip_file) and os.path.getsize(zip_file) > 0:
            completed.add(acc)
    return completed

def parse_args():
    parser = argparse.ArgumentParser(
        description="批量下载 NCBI Datasets 数据",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 首次下载
  python3 %(prog)s -b /path/to/output/ -f accessions.txt -w 5
  
  # 断点续跑
  python3 %(prog)s -b /path/to/output/ -f accessions.txt -w 5 --resume
  
  # 强制重新下载所有文件
  python3 %(prog)s -b /path/to/output/ -f accessions.txt -w 5 --force-redownload
        """
    )
    parser.add_argument(
        "--base-path", "-b",
        required=True,
        help="输出结果和下载文件的基础路径，例如 /mnt/c/Users/Administrator/Desktop/"
    )
    parser.add_argument(
        "--file-path", "-f",
        required=True,
        help="存放 accession 列表的文本文件路径，例如 /path/to/下载NCBI.txt"
    )
    parser.add_argument(
        "--email", "-e",
        default=DEFAULT_NCBI_EMAIL,
        help=f"联系邮箱（默认: {DEFAULT_NCBI_EMAIL}）"
    )
    parser.add_argument(
        "--api-key", "-k",
        default=DEFAULT_NCBI_API_KEY,
        help=f"NCBI API Key（默认已内置）"
    )
    parser.add_argument(
        "--workers", "-w",
        type=int,
        default=5,
        help="并发线程数（默认5）"
    )
    parser.add_argument(
        "--resume", "-r",
        action="store_true",
        help="启用断点续跑，跳过已下载的文件"
    )
    parser.add_argument(
        "--force-redownload",
        action="store_true",
        help="强制重新下载所有文件，忽略已存在的文件"
    )
    parser.add_argument(
        "--no-color",
        action="store_true",
        help="禁用彩色输出"
    )
    return parser.parse_args()

def read_accession_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return [accession.strip() for accession in file if accession.strip()]

def generate_urls(accession_numbers):
    URL_TEMPLATE = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/download"
        "?include_annotation_type=GENOME_FASTA"
        "&include_annotation_type=GENOME_GFF"
        "&include_annotation_type=CDS_FASTA"
        "&include_annotation_type=PROT_FASTA"
    )
    return [URL_TEMPLATE.format(acc) for acc in accession_numbers]

def append_api_key_to_url(url, api_key):
    """
    如果 api_key 非空且 URL 中还没有 api_key 参数，则追加 api_key 参数并返回新的 URL。
    """
    if not api_key:
        return url
    parsed = urlparse(url)
    qs = parse_qs(parsed.query)
    if "api_key" in qs and qs["api_key"]:
        return url  # 已有 api_key，直接返回
    qs["api_key"] = [api_key]
    new_query = urlencode(qs, doseq=True)
    new_parsed = parsed._replace(query=new_query)
    return urlunparse(new_parsed)

def rate_limit_wait():
    """
    全局速率限制：确保线程间的请求间隔至少为 MIN_INTERVAL 秒
    """
    global LAST_REQUEST_TIME
    with RATE_LOCK:
        now = time.time()
        elapsed = now - LAST_REQUEST_TIME
        if elapsed < MIN_INTERVAL:
            to_sleep = MIN_INTERVAL - elapsed
            time.sleep(to_sleep)
            LAST_REQUEST_TIME = time.time()
        else:
            LAST_REQUEST_TIME = now

def create_session(email):
    """
    创建共享的 requests.Session，带 retry 策略和常用 headers（包含 email）
    """
    session = requests.Session()
    retry_strategy = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("https://", adapter)
    # 将 email 放入 User-Agent 和 From header（不少机构希望看到 contact）
    ua = f"ncbi-datasets-downloader ({email})"
    session.headers.update({
        "User-Agent": ua,
        "From": email
    })
    return session

def download_file(session, url, output_filename, api_key, accession):
    """
    下载单个文件（stream），在每次请求前做全局速率限制，并返回结果字典
    """
    if SHUTDOWN_EVENT.is_set():
        return {"status": "cancelled", "accession": accession, "message": "下载被用户取消"}
    
    try:
        final_url = append_api_key_to_url(url, api_key)
        # 全局速率限制
        rate_limit_wait()
        
        if SHUTDOWN_EVENT.is_set():
            return {"status": "cancelled", "accession": accession, "message": "下载被用户取消"}
        
        # 发起请求（流式）
        with session.get(final_url, stream=True, timeout=60) as response:
            if response.status_code == 200:
                # 临时写入文件（避免半截文件被误判为成功）
                tmp_out = output_filename + ".part"
                add_temp_file(tmp_out)  # 添加到临时文件跟踪
                
                try:
                    # 获取文件大小用于进度显示
                    total_size = int(response.headers.get('content-length', 0))
                    
                    with open(tmp_out, 'wb') as f:
                        downloaded = 0
                        for chunk in response.iter_content(chunk_size=1024*16):
                            if SHUTDOWN_EVENT.is_set():
                                return {"status": "cancelled", "accession": accession, "message": "下载被用户取消"}
                            
                            if chunk:
                                f.write(chunk)
                                downloaded += len(chunk)
                    
                    # 下载成功后重命名并从临时文件列表移除
                    os.replace(tmp_out, output_filename)
                    remove_temp_file(tmp_out)
                    
                    # 格式化文件大小
                    size_mb = downloaded / (1024 * 1024)
                    return {
                        "status": "success",
                        "accession": accession,
                        "message": f"成功下载 {accession} ({size_mb:.2f} MB)",
                        "filename": output_filename,
                        "size": downloaded
                    }
                    
                except Exception as e:
                    # 下载过程中出错，清理临时文件
                    if os.path.exists(tmp_out):
                        try:
                            os.remove(tmp_out)
                        except:
                            pass
                    remove_temp_file(tmp_out)
                    raise e
                    
            else:
                return {
                    "status": "failed",
                    "accession": accession,
                    "message": f"HTTP错误 {response.status_code}: {url}",
                    "url": url
                }
                
    except Exception as e:
        return {
            "status": "failed",
            "accession": accession,
            "message": f"下载失败 {accession}: {str(e)}",
            "url": url
        }

def concurrent_download(urls_and_accessions, base_path, success_file, failure_file, api_key, workers, email, resume=False):
    """
    并发下载文件，支持断点续跑
    urls_and_accessions: [(url, accession), ...] 元组列表
    """
    os.makedirs(base_path, exist_ok=True)
    session = create_session(email)
    
    # 统计信息
    stats = {
        "total": len(urls_and_accessions),
        "completed": 0,
        "failed": 0,
        "skipped": 0,
        "cancelled": 0,
        "start_time": time.time(),
        "total_size": 0
    }
    
    # 如果启用断点续跑，检查已完成的文件
    completed_accessions = set()
    if resume:
        # 从进度文件加载状态
        progress_state = load_progress_state(base_path)
        completed_accessions = set(progress_state.get("completed", []))
        
        # 检查实际文件存在性
        existing_files = check_existing_files(base_path, [acc for _, acc in urls_and_accessions])
        completed_accessions.update(existing_files)
        
        if completed_accessions:
            print_info(f"发现 {len(completed_accessions)} 个已完成的下载，将跳过")
    
    # 过滤需要下载的项目
    to_download = [(url, acc) for url, acc in urls_and_accessions if acc not in completed_accessions]
    stats["skipped"] = len(urls_and_accessions) - len(to_download)
    
    if stats["skipped"] > 0:
        print_progress(f"跳过 {stats['skipped']} 个已完成的文件")
    
    if not to_download:
        print_success("所有文件都已下载完成!")
        return stats
    
    print_info(f"准备下载 {len(to_download)} 个文件，使用 {workers} 个并发线程")
    
    # 创建进度条 (如果可用)
    pbar = None
    if HAS_TQDM and len(to_download) > 1:
        pbar = tqdm(total=len(to_download), desc="下载进度", 
                   bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]")
    
    with ThreadPoolExecutor(max_workers=workers) as executor, \
         open(success_file, 'a', encoding='utf-8') as sf, \
         open(failure_file, 'a', encoding='utf-8') as ff:
        
        future_to_info = {}
        
        # 提交所有下载任务
        for url, accession in to_download:
            if SHUTDOWN_EVENT.is_set():
                break
                
            out_fname = os.path.join(base_path, f"{accession}_downloaded.zip")
            future = executor.submit(download_file, session, url, out_fname, api_key, accession)
            future_to_info[future] = {"url": url, "accession": accession}
        
        # 处理完成的任务
        completed_in_session = set()
        failed_in_session = set()
        
        try:
            for future in as_completed(future_to_info):
                if SHUTDOWN_EVENT.is_set():
                    break
                    
                info = future_to_info[future]
                result = future.result()
                
                if result["status"] == "success":
                    print_success(result["message"])
                    sf.write(f"{datetime.now().isoformat()}: {result['message']}\n")
                    sf.flush()
                    stats["completed"] += 1
                    stats["total_size"] += result.get("size", 0)
                    completed_in_session.add(result["accession"])
                    
                elif result["status"] == "cancelled":
                    print_warning(result["message"])
                    stats["cancelled"] += 1
                    
                else:  # failed
                    print_error(result["message"])
                    ff.write(f"{datetime.now().isoformat()}: {result['message']}\n")
                    ff.flush()
                    stats["failed"] += 1
                    failed_in_session.add(result["accession"])
                
                # 更新进度条
                if pbar:
                    pbar.update(1)
                    pbar.set_postfix({
                        "成功": stats["completed"],
                        "失败": stats["failed"],
                        "取消": stats["cancelled"]
                    })
                
                # 定期保存进度
                if (stats["completed"] + stats["failed"]) % 10 == 0:
                    all_completed = completed_accessions | completed_in_session
                    save_progress_state(base_path, all_completed, failed_in_session)
        
        except KeyboardInterrupt:
            print_warning("\n检测到键盘中断，正在停止下载...")
            SHUTDOWN_EVENT.set()
            
        finally:
            if pbar:
                pbar.close()
            
            # 最终保存进度
            all_completed = completed_accessions | completed_in_session
            save_progress_state(base_path, all_completed, failed_in_session)
    
    # 打印最终统计信息
    elapsed = time.time() - stats["start_time"]
    print_info("\n" + "="*50)
    print_info("下载完成统计:")
    print_success(f"  总计: {stats['total']} 个文件")
    print_success(f"  成功: {stats['completed']} 个")
    if stats["skipped"] > 0:
        print_info(f"  跳过: {stats['skipped']} 个 (已存在)")
    if stats["failed"] > 0:
        print_error(f"  失败: {stats['failed']} 个")
    if stats["cancelled"] > 0:
        print_warning(f"  取消: {stats['cancelled']} 个")
    print_info(f"  用时: {elapsed:.1f} 秒")
    if stats["total_size"] > 0:
        size_mb = stats["total_size"] / (1024 * 1024)
        print_info(f"  下载大小: {size_mb:.2f} MB")
    print_info("="*50)
    
    return stats

def main():
    global MIN_INTERVAL, HAS_COLORAMA
    args = parse_args()

    # 禁用彩色输出如果用户指定
    if args.no_color:
        HAS_COLORAMA = False

    # 设置信号处理器
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    BASE_PATH = args.base_path.rstrip("/") + "/"
    FILE_PATH = args.file_path
    EMAIL = args.email
    API_KEY = args.api_key
    WORKERS = max(1, args.workers)
    RESUME = args.resume
    FORCE_REDOWNLOAD = args.force_redownload

    # 如果强制重新下载，则禁用断点续跑
    if FORCE_REDOWNLOAD:
        RESUME = False
        print_warning("启用强制重新下载模式，将覆盖所有现有文件")

    OUTPUT_FILE  = os.path.join(BASE_PATH, "generated_urls.txt")
    SUCCESS_FILE = os.path.join(BASE_PATH, "success.txt")
    FAILURE_FILE = os.path.join(BASE_PATH, "failure.txt")

    os.makedirs(BASE_PATH, exist_ok=True)

    # 打印启动信息
    print_info("="*60)
    print_info("NCBI Assembly 批量下载工具 (增强版)")
    print_info("="*60)
    print_info(f"输出目录: {BASE_PATH}")
    print_info(f"输入文件: {FILE_PATH}")
    print_info(f"并发线程: {WORKERS}")
    print_info(f"断点续跑: {'启用' if RESUME else '禁用'}")
    print_info(f"彩色输出: {'启用' if HAS_COLORAMA else '禁用'}")
    if HAS_TQDM:
        print_info("进度条: 启用")
    print_info("-"*60)

    # 根据是否有 API Key 设置允许的全局 rps
    if API_KEY:
        ALLOWED_RPS = 10.0   # 带 key 的建议速率
        print_info("检测到 API Key，使用较高请求频率")
    else:
        ALLOWED_RPS = 3.0    # 无 key 更保守
        print_warning("未提供 API Key，使用保守请求频率")
    
    # 计算最小间隔（多个并发线程共享）
    MIN_INTERVAL = float(WORKERS) / float(ALLOWED_RPS)
    print_info(f"请求间隔: {MIN_INTERVAL:.2f}秒 (目标频率: {ALLOWED_RPS} req/s)")

    try:
        # 1. 读取 accession 列表，并生成 URL
        print_progress("正在读取 accession 列表...")
        if not os.path.exists(FILE_PATH):
            print_error(f"输入文件不存在: {FILE_PATH}")
            return 1
            
        accession_numbers = read_accession_file(FILE_PATH)
        if not accession_numbers:
            print_error("输入文件为空或没有有效的 accession")
            return 1
            
        print_success(f"成功读取 {len(accession_numbers)} 个 accession")
        
        urls = generate_urls(accession_numbers)
        urls_and_accessions = list(zip(urls, accession_numbers))

        # 2. 保存 URL 列表
        print_progress("正在保存下载链接...")
        with open(OUTPUT_FILE, 'w', encoding='utf-8') as outf:
            outf.write(f"# Generated on {datetime.now().isoformat()}\n")
            outf.write(f"# Total URLs: {len(urls)}\n")
            for u in urls:
                outf.write(u + "\n")
        print_success(f"下载链接已保存至: {OUTPUT_FILE}")

        # 3. 清理之前的临时文件 (如果不是断点续跑模式)
        if not RESUME or FORCE_REDOWNLOAD:
            print_progress("正在清理旧的临时文件...")
            temp_pattern = os.path.join(BASE_PATH, "*.part")
            temp_files = glob.glob(temp_pattern)
            for temp_file in temp_files:
                try:
                    os.remove(temp_file)
                except:
                    pass
            if temp_files:
                print_info(f"清理了 {len(temp_files)} 个临时文件")

        # 4. 并发下载
        print_progress("开始下载...")
        stats = concurrent_download(
            urls_and_accessions, BASE_PATH, SUCCESS_FILE, FAILURE_FILE, 
            API_KEY, WORKERS, EMAIL, RESUME
        )
        
        # 5. 最终清理和总结
        cleanup_temp_files()
        
        if stats["failed"] == 0 and stats["cancelled"] == 0:
            print_success("🎉 所有文件下载完成!")
        elif stats["completed"] > 0:
            print_warning(f"⚠️  部分文件下载完成: {stats['completed']}/{stats['total']}")
        else:
            print_error("❌ 没有文件成功下载")
            
        print_info(f"详细日志: 成功={SUCCESS_FILE}, 失败={FAILURE_FILE}")
        return 0 if stats["failed"] == 0 else 1
        
    except KeyboardInterrupt:
        print_warning("\n用户中断操作")
        cleanup_temp_files()
        return 130
    except Exception as e:
        print_error(f"程序异常: {str(e)}")
        cleanup_temp_files()
        return 1

if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except Exception as e:
        print_error(f"致命错误: {str(e)}")
        cleanup_temp_files()
        sys.exit(1)

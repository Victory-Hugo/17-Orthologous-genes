[GToTree](https://github.com/AstrobioMike/GToTree/wiki) 在类 Unix 的命令行环境中运行。这意味着它可以在 Mac 和 Linux 计算机上使用它们附带的标准终端程序工作。要在 Windows 计算机上使用 GToTree，我建议安装 [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install)，然后在 WSL 终端中安装 [miniconda](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda) 的 Linux 版本。然后按照下面的说明使用 `conda` 安装将在 WSL 环境中工作 👍 

## Conda 快速开始！
如果你还没有这个出色的包管理器 [conda](https://conda.io/docs/)，我**强烈**建议你安装它。这真的不是讨论它为什么如此有用的场合，但它确实很有用，我保证 🙂

要让 conda 正常运行（这非常快），你可以按照适合你的系统的说明安装 miniconda（一个轻量级版本），从 [这里](https://conda.io/en/latest/miniconda.html) 开始。你需要 python 3.X 版本，而且很可能是 64 位版本。如果你想在某个时候了解更多关于 conda 的信息，我有 [一个介绍页面](https://astrobiomike.github.io/unix/conda-intro) 🙂

---

以下行将创建一个 gtotree conda 环境并安装 GToTree，你需要在基础 conda 环境中运行这些命令：

```bash
conda create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda gtotree
```

**完成！**

现在你应该可以使用 `conda activate gtotree` 和 `conda deactivate gtotree` 进入和退出环境。如果你进入环境并运行以下命令：

```bash
gtt-hmms
```

它将打印出 GToTree 默认 HMM 目录的位置，并列出可用的预构建 HMM。如果你输入 `GToTree` 且不带任何参数，你可以看到帮助菜单。

## 测试运行
你可以运行一个大约需要 3 分钟的测试：

```bash
gtt-test.sh
```

标准输出的结尾应该如下所示：

```
#################################################################################
####                                 Done!!                                  ####
#################################################################################

  Overall, 12 genomes of the input 14 were retained (see notes below).

    Tree written to:
        GToTree-test-output/GToTree-test-output.tre

    Alignment written to:
        GToTree-test-output/Aligned_SCGs_mod_names.faa

    Main genomes summary table written to:
        GToTree-test-output/Genomes_summary_info.tsv
```

## 更新到新版本
如果想更新到最新的 GToTree 版本，最好删除之前的 conda 环境并重新安装。这可以通过以下方式完成：

```bash
# 从 gtotree conda 环境之外（假设环境名称与上面的安装相同）
conda env remove -n gtotree

# 然后按照上面的方式在新环境中重新安装
# 如果需要，首先安装 mamba（用于更快的 conda 安装）
conda install -n base -c conda-forge mamba
mamba create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda -c defaults gtotree
```

然后可以使用 `conda activate gtotree` 激活新环境。

---

## 无需 conda 的安装（不推荐）

再次强调，**强烈**推荐使用 conda 安装，因为它在不同系统上更加健壮。但要尝试不使用 conda 安装，可以将 GToTree 下载并解压缩到系统上你想要的任何位置（确保将下面的版本更改为 [这里](https://github.com/AstrobioMike/GToTree/releases/latest) 找到的最新版本）：

```bash
curl -L https://github.com/AstrobioMike/GToTree/archive/v1.5.22.tar.gz -o GToTree-v1.5.22.tar.gz
tar -xzvf GToTree-v1.5.22.tar.gz
```

### 将 bin 添加到你的 PATH
现在我们需要将 "bin" 目录添加到我们的 PATH 中（如果你不熟悉 PATH 是什么，想要了解更多信息，请参阅 [这里](https://astrobiomike.github.io/bash/modifying_your_path)）。
我们可以这样做：将目录切换到 bin，然后在 `echo` 命令中使用 `pwd` 将完整路径放入我们的 PATH：

```bash
cd GToTree-1.5.22/bin # 确保你在这个 bin 目录中
echo "export PATH=\"$(pwd):\$PATH\"" >> ~/.bash_profile
```

### 添加包含的 HMM 文件的路径
如果你想更容易地使用包含的单拷贝基因 HMM 配置文件，你还可以在你的 bash profile 中添加一个变量，这样你就不需要在使用它们时提供它们的完整路径。如果你将目录切换到 "hmm_sets" 目录，可以按照与上面类似的方式进行：

```bash
cd ../hmm_sets/ # 从我们上面所在的位置开始
echo "export GToTree_HMM_dir=\"$(pwd)/\"" >> ~/.bash_profile
```

最后要做的是 `source` 我们刚刚修改的 ~/.bash_profile，以便这些更改在我们当前的会话中生效：

```bash
source ~/.bash_profile
```

你可以不带任何参数运行 `gtt-hmms`，以确保默认 HMM 目录已设置，并查看当前可用的 HMM 文件可以更具体地针对哪些分类单元。

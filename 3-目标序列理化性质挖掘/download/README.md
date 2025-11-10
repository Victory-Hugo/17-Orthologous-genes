# Signalp4
SignalP 4.1 是一个用于预测蛋白质信号肽的工具。信号肽是位于蛋白质N端的一段短序列，指导蛋白质进入分泌途径。SignalP 4.1 结合了神经网络和隐马尔可夫模型（HMM）来提高预测的准确性。

# 安装
```sh
conda install predector::signalp4
```

然后：
```sh
signalp4-register signalp-4.1g.Linux.tar.gz
```

# TMHMM2
TMHMM2 是一个用于预测跨膜蛋白质的工具。它基于隐马尔可夫模型（HMM），能够准确地识别蛋白质中的跨膜区域。

# 安装

```sh
conda install predector::tmhmm
tmhmm2-register tmhmm-2.0c.Linux.tar.gz
```

# 在WSL中使用的问题
**为什么 SignalP-4.1 在 WSL2 下会崩溃（segfault）**，以及 **如何修复**。
---

## 要点总结（一句话版）

SignalP-4.1 捆绑了一个很老的辅助二进制（文中叫 `nnhowplayer` / helper），它使用旧式的 x86 `vsyscall` 接口。WSL2 出于安全默认把 `vsyscall` 设为 `none`，于是每次那个老二进制触发该接口就会**访问非法地址并崩溃（segfault）**。解决办法是让 WSL2 内核**模拟（emulate）vsyscall**，或改用真正的 Linux 虚拟机／新版 SignalP（v5/6）或其它工具。

---

## 背后技术说明（不必全部记住）

* `vsyscall`：一种旧的“虚拟系统调用”接口，早期 x86/Linux 用它做快速系统调用。因为安全问题，现代内核通常禁用或限制它（用 `vsyscall=none`）。
* 那个 SignalP-4.1 的 helper 是编译得很老、还会用到 `vsyscall` 的二进制。WSL2 默认禁用 `vsyscall`，导致 helper 每次尝试访问对应地址就触发“fatal signal 11”（segfault）。
* 解决办法一是允许内核对 `vsyscall` 做“模拟”（`vsyscall=emulate`），这样老二进制就不会崩溃；另一种是用不依赖老 helper 的环境（例如真正的 Linux VM、Docker 容器，或升级到 SignalP v5/6 或用 alternative 工具）。

---

## 如果你想修复（在 WSL2 中允许模拟 vsyscall）

> 在 Windows 上做的改动 —— 请按照下面步骤操作

1. 在 Windows 的用户目录下创建或编辑文件 `%UserProfile%\.wslconfig`（就是 `C:\Users\<你的用户名>\.wslconfig`）。可以用记事本打开：
   打开 PowerShell（或 cmd），运行：

   ```powershell
   notepad $env:USERPROFILE\.wslconfig
   ```

2. 在文件中加入（或修改为）下面内容：

   ```
   [wsl2]
   kernelCommandLine = vsyscall=emulate
   ```

   保存并退出记事本。

3. 在 Windows 中运行（PowerShell 或 cmd）：

   ```powershell
   wsl --shutdown
   ```

   这会停止所有 WSL2 实例并在下一次启动时使用新的内核命令行参数。

4. 重新打开你的 WSL 发行版（例如 Ubuntu），再运行 SignalP-4.1（你之前的命令）。那时 helper 不再因 vsyscall 而崩溃，程序应能正常运行（或至少跑到下一步报别的错误）。

---

## 风险与注意事项

* `vsyscall=emulate` 允许内核模拟旧接口，理论上降低了一点安全性（因为开启了对旧接口的兼容），但对典型桌面/开发机风险很小。如果在生产/受控环境（或对安全要求很高的服务器）上要谨慎评估。
* 如果你不希望改内核参数，替代方案更加稳妥（见下）。

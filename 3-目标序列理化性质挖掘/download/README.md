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
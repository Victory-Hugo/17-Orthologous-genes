## 必须扫描的 PfamA 域（核心证据：NBD + exporter 型 TMD）

这些 PfamA 域能回答序列是否为 ABC 转运 ATPase，以及是否具有 exporter 的跨膜域：

### NBD（ATPase 域）
- **ABC_tran**：最核心的 ABC transporter NBD
- **ABC_N**：常见 NBD 相关域，有时与 ABC_tran 互补
- **ABC_ATPase**：部分 Pfam 版本提供的不同颗粒度 ATPase 域

> 实操建议：通常 ABC_tran + ABC_N 就能稳妥识别 NBD。

### TMD（跨膜域，exporter 侧重点）
- **ABC2_membrane**（包括 ABC2_membrane_2/3/4/5/6/7）
- **ABC_membrane**（包括 ABC_membrane_2 / ABC_membrane_3）
- **ABC_export**：对“exporter 型”识别很有帮助

> 实操建议：至少包含 ABC2_membrane、ABC2_membrane_2、ABC2_membrane_3、ABC_membrane、ABC_export。  
> 其余 ABC2_membrane_4~7 可根据命中情况补充（部分模型为同类细分/补充）。

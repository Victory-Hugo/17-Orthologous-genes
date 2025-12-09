# 下载Pfam-A数据库

Pfam-A数据库是一个包含蛋白质家族的数据库，广泛用于生物信息学研究。

# 下载步骤

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
```

# 解压

```bash
gunzip Pfam-A.hmm.gz
```
# 建立 HMM 索引
```bash
hmmpress Pfam-A.hmm
```
我准备了一个备份在OSS：

```bash
https://s1-software.oss-cn-chengdu.aliyuncs.com/Pfam-A.zip
```

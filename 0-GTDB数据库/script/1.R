library(ggtree)
library(treeio)
library(tidyverse)
class_tree <- read.tree("/home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/output/1-classify/classify/classify_Result.backbone.bac120.classify.tree")

ggtree(class_tree) -> p1
ggsave(p1, filename = "/home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/output/1-classify/class_tree.pdf")


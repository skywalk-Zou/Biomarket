# 加载必要的库
```R
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2", "cowplot"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(readxl)
library(dplyr)
library(cowplot)
```



# 读取 Excel 文件
```R
digest <- read_excel("D:/Users/skywa/Desktop/Network/Network/Digest.xlsx")
scu_target <- read_excel("D:/Users/skywa/Desktop/Network/Network/SCU-target.xlsx")
```




# 提取基因列表
```R
genes_digest <- digest$GENE
genes_scu_target <- scu_target$GENE
```



# 找到交集基因
```R
intersecting_genes <- intersect(genes_digest, genes_scu_target)
lapply(unique(intersecting_genes), function(x) gsub("\\s+", "", x))
print(intersecting_genes)
```



# 将基因符号转换为 Entrez ID
```R
gene_ids <- bitr(intersecting_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
```



# 进行 KEGG 富集分析，展示前 20 条通路
```R
kegg_enrichment <- enrichKEGG(gene = gene_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05, qvalueCutoff = 0.2)
kegg_enrichment <- kegg_enrichment %>% head(n = 20)
```



# 进行 GO 富集分析，分别分析 BP、MF 和 CC，并展示前 10 条通路
```R
go_enrichment_BP <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05) %>% head(n = 10)
go_enrichment_MF <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = 0.05) %>% head(n = 10)
go_enrichment_CC <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", pvalueCutoff = 0.05) %>% head(n = 10)
```



# 自定义点图函数
```R
custom_dotplot <- function(enrichment_result, title) {
  ggplot(enrichment_result, aes(x = Count, y = reorder(Description, Count), size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),  # 去掉背景网格
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Gene Count",
      y = "Pathway/GO Category",
      size = "Gene Count",
      color = "-log10(p.adjust)",
      title = title
    )
}
```



# 对富集结果进行排序，综合考虑count和-log10(p.adjust)值
```R
go_enrichment_BP_sorted <- go_enrichment_BP[order(go_enrichment_BP$Count, -log10(go_enrichment_BP$p.adjust)), ]
go_enrichment_MF_sorted <- go_enrichment_MF[order(go_enrichment_MF$Count, -log10(go_enrichment_MF$p.adjust)), ]
go_enrichment_CC_sorted <- go_enrichment_CC[order(go_enrichment_CC$Count, -log10(go_enrichment_CC$p.adjust)), ]
```



# 生成 GO 富集结果的点图，分别为 BP、MF 和 CC
```
go_plot_BP <- custom_dotplot(go_enrichment_BP_sorted, "GO Enrichment - Biological Process")
go_plot_MF <- custom_dotplot(go_enrichment_MF_sorted, "GO Enrichment - Molecular Function")
go_plot_CC <- custom_dotplot(go_enrichment_CC_sorted, "GO Enrichment - Cellular Component")
```



# 对KEGG富集结果进行排序，综合考虑count和-log10(p.adjust)值
kegg_enrichment_sorted <- kegg_enrichment[order(kegg_enrichment$Count, -log10(kegg_enrichment$p.adjust)), ]

# 生成 KEGG 富集结果的点图
kegg_plot <- custom_dotplot(kegg_enrichment_sorted, "KEGG Pathway Enrichment")

# 拼接 GO 分析的三张图片
go_combined_plot <- cowplot::plot_grid(go_plot_BP, go_plot_MF, go_plot_CC, labels = c("BP", "MF", "CC"), ncol = 1, align = "v")

# 保存图像
ggsave("kegg_enrichment.png", kegg_plot, width = 12, height = 8)
ggsave("go_enrichment_combined.png", go_combined_plot, width = 12, height = 24)

# 打印图像到控制台
print(kegg_plot)
print(go_combined_plot)
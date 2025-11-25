# 相关性聚类热图


import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import os

# 1. 读取 h5ad 文件
adata = sc.read_h5ad(r"result.h5ad")  

if 'rank_genes_groups' not in adata.uns:
    if adata.raw is not None:
        sc.tl.rank_genes_groups(adata, groupby='leiden', use_raw=True, method='wilcoxon')
    else:
        sc.pp.log1p(adata) 
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# 2. 提取每个 cluster 的 top10 marker
markers = []
for group in adata.uns['rank_genes_groups']['names'].dtype.names:
    top_genes = adata.uns['rank_genes_groups']['names'][group][:10]
    top_lfc   = adata.uns['rank_genes_groups']['logfoldchanges'][group][:10]
    for g, lfc in zip(top_genes, top_lfc):
        markers.append({'cluster': group, 'gene': g, 'avg_log2FC': lfc})

markers_df = pd.DataFrame(markers)

# 3. 取出所有出现过的基因
all_marker_genes = markers_df['gene'].drop_duplicates().tolist()

# 4. 提取表达矩阵
if adata.raw is not None:
    expr = adata.raw.to_adata()[:, all_marker_genes].to_df()
else:
    expr = adata[:, all_marker_genes].to_df()

# 5. 计算相关性
corr = expr.corr(method='pearson')

# 6. 基因排序 —— 按 cluster + log2FC 排，让同一个 cluster 的基因聚在一起
# order = markers_df.sort_values(['cluster', 'avg_log2FC'], ascending=[True, False])['gene'].tolist()
# # 注意：可能有重复基因，这里只取第一次出现的位置，保证长度一致
# order = list(dict.fromkeys(order))  # 去重
# corr_ordered = corr.loc[order, order]

# 6. 让同一个 cluster 的基因连续排列（聚成色块）
# 把你原来的这两行全部删掉：
# order = markers_df.sort_values(['cluster', 'avg_log2FC'], ascending=[True, False])['gene'].tolist()
# order = list(dict.fromkeys(order))

markers_df = markers_df.sort_values(['cluster', 'avg_log2FC'], ascending=[True, False])
# 先去重（保留第一次出现的 cluster 归属），再取基因顺序
order = markers_df.drop_duplicates('gene', keep='first')['gene'].tolist()
corr_ordered = corr.loc[order, order]



# 7. 每个基因对应其所属 cluster 的颜色
gene_to_cluster = {}
for _, row in markers_df.iterrows():
    if row['gene'] not in gene_to_cluster:
        gene_to_cluster[row['gene']] = row['cluster']

# 生成安全的颜色映射
# unique_clusters = sorted(markers_df['cluster'].unique(), key=lambda x: int(x))
unique_clusters = sorted(markers_df['cluster'].unique())
n_clusters = len(unique_clusters)
palette = sns.color_palette("husl", n_clusters)  # husl 可以无限扩展颜色
cluster_color_map = dict(zip(unique_clusters, palette))

row_colors = [cluster_color_map[gene_to_cluster.get(g, 'unknown')] for g in order]

gene_to_cluster = dict(zip(markers_df.drop_duplicates('gene', keep='first')['gene'],
                          markers_df.drop_duplicates('gene', keep='first')['cluster']))

row_colors = [cluster_color_map[gene_to_cluster[g]] for g in order]


# 8. 画图
g = sns.clustermap(
    corr_ordered,
    cmap='RdBu_r', center=0, vmin=-1, vmax=1,
    # figsize=(max(10, len(order)/4), max(10, len(order)/4)),
    figsize=(13,13),
    # annot=True, fmt='.2f', annot_kws={'size': 7},
    annot=False, # 热图中的相关性数字不展示
    linewidths=0.5, linecolor='lightgray',
    row_colors=row_colors,
    col_colors=row_colors,
    dendrogram_ratio=(0.02, 0.05),
    cbar_pos=(0.02, 0.83, 0.03, 0.13),
    tree_kws=dict(linewidths=1.5)
    # row_cluster=True,            # 保留左侧树
    # col_cluster=False  # 上侧不保留，但顺序会乱
)
g.fig.subplots_adjust(right=0.78)

handles = [Patch(facecolor=cluster_color_map[c], 
                 edgecolor='black', 
                 linewidth=1.2) 
           for c in unique_clusters]
unique_clusters = sorted(markers_df['cluster'].unique(), key=lambda x: int(x))
labels  = [f'Cluster {c}' for c in unique_clusters]

g.fig.legend(handles, labels,
             title='Cell Cluster',
             title_fontsize=14,
             fontsize=12,
             loc='center left',
             bbox_to_anchor=(0.85, 0.68),
             frameon=True,
             fancybox=False,
             edgecolor='black',
             ncol=1)     

g.fig.subplots_adjust(top=0.94)
plt.suptitle('Top5 Marker Genes Correlation (Pearson)', 
             fontsize=18, fontweight='bold', y=0.97, x=0.45)

plt.savefig('AllClusters_top10_markers_correlation_fixed.png', dpi=300, bbox_inches='tight')
plt.savefig('AllClusters_top10_markers_correlation_fixed.pdf', dpi=300, bbox_inches='tight')
plt.show()

print(f"finish! total {len(unique_clusters)} cluster，{len(order)}  marker gene")

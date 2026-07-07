# -*- coding: utf-8 -*-
"""
Created on 2025/7/11
@author: ljzzz
"""
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib.ticker import PercentFormatter
from matplotlib.font_manager import FontProperties
import warnings
import matplotlib.pyplot as plt
import matplotlib as mpl
warnings.filterwarnings("ignore")
# 设置全局字体和文本样式
mpl.rcParams["font.family"] = "Times New Roman"  # 设置字体为 Times New Roman
mpl.rcParams["axes.labelsize"] = 24             # 坐标轴标签（xlabel, ylabel）字体大小
mpl.rcParams["axes.labelweight"] = "bold"       # 坐标轴标签加粗
mpl.rcParams["xtick.labelsize"] = 12            # x轴刻度数字字体大小
mpl.rcParams["ytick.labelsize"] = 12            # y轴刻度数字字体大小
# mpl.rcParams["xtick.labelweight"] = "bold"    # 如果需要刻度也加粗，可以取消注释
# mpl.rcParams["ytick.labelweight"] = "bold"    # 同上
mpl.rcParams["figure.titlesize"] = 24           # 整个 Figure 的标题（suptitle）字体大小
mpl.rcParams["axes.titlesize"] = 24             # 每个子图的标题（plt.title() 或 ax.set_title()）字体大小


"""制作OXA变异子饼图及随时间变化分布图"""
df = pd.read_excel(r"dataset_4886_ecoli_blaOXA_metadata.xlsx")
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(211)
colors = list(plt.get_cmap('Pastel1').colors)
print(len(colors))

df_OXA = df["variants"]
print(df_OXA.value_counts())
print(df_OXA.value_counts(normalize=True))
OXA_name = df_OXA.value_counts().index  #获得OXA_label
color_dict = dict(zip(OXA_name,colors))
print(color_dict)
OXA_counts = df_OXA.value_counts().values
ax1.pie(OXA_counts,
        colors = [color_dict[i] for i in OXA_name],
        labels = OXA_name,
        autopct = '%1.2f%%',
        # textprops = {'fontstyle': 'italic'}

        )

# 创建斜体字体属性
italic_font = FontProperties(style='italic')
# 绘制图例，并设置字体为斜体
ax1.legend(prop=italic_font,bbox_to_anchor=(1.05,1))

# 绘制OXA检出数目变化图
df_OXA_year = df[["variants","Year"]]
df_OXA_year = df_OXA_year[df_OXA_year["Year"] != "Unknown"]
# df_OXA_year = df_OXA_year[df_OXA_year["Year"] != "2025"]

df_OXA_year["Year"] = df_OXA_year["Year"].astype("int64")

df_OXA_year.loc[df_OXA_year["Year"] <= 2011,"Year"] = 2011

OXA_counts_by_variant = df_OXA_year.groupby(["Year", "variants"]).size().unstack(fill_value=0)
print(OXA_counts_by_variant)
years = OXA_counts_by_variant.index
variants = OXA_counts_by_variant.columns
OXA_counts = OXA_counts_by_variant.values
print(OXA_counts)
x = range(len(years))
width = 0.8  # 每个条形的宽度
# 绘制堆叠条形图
bottom = [0] * len(years)
ax1 = fig.add_subplot(223)
for i, variant in enumerate(variants):
    color = color_dict[variant]
    ax1.bar(x, OXA_counts[:, i], width, bottom=bottom, color=color, label=variant)
    bottom = [b + p for b, p in zip(bottom, OXA_counts[:, i])]
ax1.set_xticks(x)
ax1.set_xticklabels(years)
for xticklabel in ax1.get_xticklabels():
        xticklabel.set_rotation(45)
        xticklabel.set_weight("bold")
        xticklabel.set_size(14)
for yticklabel in ax1.get_yticklabels():
        yticklabel.set_weight("bold")
        yticklabel.set_size(14)

#堆积条形图优化代码
ax1.tick_params(axis="both",which="minor",length=0)
ax1.grid(False)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.set_ylim(0,700)

ax1.set_xlabel("Years")
ax1.set_ylabel("Number")


# 绘制百分比图
ax2 = fig.add_subplot(224)
# 计算每年总检出数
year_total = OXA_counts_by_variant.sum(axis=1)
# 计算比例（每行/每年）
variant_proportions = OXA_counts_by_variant.div(year_total, axis=0)

proportions = variant_proportions.values
print("proportions")
print(proportions)
# 绘制堆叠条形图
bottom = [0] * len(years)
# 绘制分组条形图
for i, variant in enumerate(variants):
    color = color_dict[variant]
    ax2.bar(x, proportions[:, i], width, bottom=bottom, color=color, label=variant)
    bottom = [b + p for b, p in zip(bottom, proportions[:, i])]
ax2.set_xticks(x)
ax2.set_xticklabels(years)
for xticklabel in ax2.get_xticklabels():
        xticklabel.set_rotation(45)
        xticklabel.set_weight("bold")
        xticklabel.set_size(14)
for yticklabel in ax2.get_yticklabels():
        yticklabel.set_weight("bold")
        yticklabel.set_size(14)
ax2.tick_params(axis="both",which="minor",length=0)
ax2.grid(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.set_ylim(0,1)
ax2.yaxis.set_major_formatter(PercentFormatter(xmax=1))


# 调整布局
plt.tight_layout()
plt.savefig("fig1_OXA_variant_year.svg",
            dpi=600,
            bbox_inches="tight",
            transparent=True)



plt.show()
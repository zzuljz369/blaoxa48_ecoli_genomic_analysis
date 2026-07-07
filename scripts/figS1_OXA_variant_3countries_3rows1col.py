# -*- coding: utf-8 -*-
"""
Created on 2025/7/11
@author: ljzzz
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.ticker import PercentFormatter
from matplotlib.font_manager import FontProperties
import warnings
import matplotlib as mpl

warnings.filterwarnings("ignore")

# ===================== 全局字体样式配置 =====================
mpl.rcParams['svg.hashsalt'] = 'optimize'
mpl.rcParams['svg.fonttype'] = 'none'  # SVG字体为矢量文字（可编辑）
mpl.rcParams['pdf.fonttype'] = 42  # PDF字体为TrueType（AI可编辑/替换）
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["font.family"] = "Times New Roman"
mpl.rcParams["axes.labelsize"] = 24
mpl.rcParams["axes.labelweight"] = "bold"
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12
mpl.rcParams["figure.titlesize"] = 24
mpl.rcParams["axes.titlesize"] = 24

italic_font = FontProperties(style='italic')

# ===================== 读取数据、筛选目标三国 =====================
df = pd.read_excel(r"dataset_4886_ecoli_blaOXA_metadata.xlsx")
target_countries = ["France", "USA", "United Kingdom","Spain",'Australia','Ireland']
df_all = df[df["Country"].isin(target_countries)].copy()

# 全局配色（所有子图统一配色）
all_variants = df_all["variants"].value_counts().index
colors = list(plt.get_cmap('Pastel1').colors)
color_dict = dict(zip(all_variants, colors))

# ===================== 3行1列画布 =====================
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(14, 16), tight_layout=True)
axes = axes.flatten()
country_axes = dict(zip(target_countries, axes))

# 循环每个国家分别绘图
for country in target_countries:
    ax = country_axes[country]
    # 筛选当前国家数据
    df_sub = df_all[df_all["Country"] == country].copy()
    df_year = df_sub[["variants", "Year"]].copy()

    # 年份清洗
    df_year = df_year[df_year["Year"] != "Unknown"]
    df_year["Year"] = df_year["Year"].astype("int64")
    df_year.loc[df_year["Year"] <= 2011, "Year"] = 2011

    # 年份-变异子透视表
    pivot = df_year.groupby(["Year", "variants"]).size().unstack(fill_value=0)
    years = pivot.index
    variants = pivot.columns
    x = list(range(len(years)))
    width = 0.8

    # 计算每年占比
    year_total = pivot.sum(axis=1)
    prop = pivot.div(year_total, axis=0).values

    # 堆叠柱状图
    bottom = np.zeros(len(years))
    for i, var in enumerate(variants):
        ax.bar(x, prop[:, i], width, bottom=bottom, color=color_dict[var], label=var)
        bottom += prop[:, i]

    # 坐标轴美化
    ax.set_xticks(x)
    ax.set_xticklabels(years)
    for lab in ax.get_xticklabels():
        lab.set_rotation(45)
        lab.set_weight("bold")
        lab.set_size(14)
    for lab in ax.get_yticklabels():
        lab.set_weight("bold")
        lab.set_size(14)

    ax.tick_params(axis="both", which="minor", length=0)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))

    # 标题、坐标轴标签
    ax.set_title(f"{country}", fontsize=24, weight="bold")
    # ax.set_xlabel("Year", fontsize=24, weight="bold")
    ax.set_ylabel("Proportion", fontsize=24, weight="bold")

# 保存图片
plt.savefig("OXA_variant_3countries_3rows1col.svg", dpi=600, bbox_inches="tight", transparent=True)
plt.show()
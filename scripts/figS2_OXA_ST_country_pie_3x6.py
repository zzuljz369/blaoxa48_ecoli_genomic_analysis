# coding:utf-8
# Author:Ljz
# date:2025/10/24 16:53
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Patch
import matplotlib as mpl
# 全局字体设置
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
# 读取数据
df = pd.read_csv(r"./dataset_4886_ecoli_blaOXA_metadata.csv")

# 基础过滤
# df = df[df["MLST"] != "Unknown"]
# df = df[df["Year"] != "2025"]
# df = df[df["Year"] != "Unknown"]

# 分组定义
OXA_variants = ["blaOXA-48", "blaOXA-181", "blaOXA-244"]
target_countries = ["France", "USA", "United Kingdom", "Spain", "Australia", "Ireland"]
STs = ["ST10", "ST1284", "ST131", "ST167", "ST354", "ST361", "ST38", "ST405", "ST410", "ST69", "ST940"]

# MLST归类
df['MLST'] = "ST" + df["MLST"]
df['MLST'] = np.where(df['MLST'].isin(STs), df['MLST'], 'Other STs')
ST_order = STs + ["Other STs"]

# 完整配色（保留你原有）
color_dict = {
    'blaOXA-48': (0.984313725490196, 0.7058823529411765, 0.6823529411764706),
    'blaOXA-181': (0.7019607843137254, 0.803921568627451, 0.8901960784313725),
    'blaOXA-244': (0.8, 0.9215686274509803, 0.7725490196078432),
    'blaOXA-484': (0.8705882352941177, 0.796078431372549, 0.8941176470588236),
    'blaOXA-232': (0.996078431372549, 0.8509803921568627, 0.6509803921568628),
    'blaOXA-204': (1.0, 1.0, 0.8),
    'blaOXA-162': (0.8980392156862745, 0.8470588235294118, 0.7411764705882353),
    'blaOXA-245': (0.9921568627450981, 0.8549019607843137, 0.9254901960784314),
    'blaOXA-163': (0.9490196078431372, 0.9490196078431372, 0.9490196078431372),
    'USA': '#E995E6', 'France': '#BAEE8C', 'Others': '#D2D84C', 'United Kingdom': '#5D05E8',
    'Singapore': '#89CFB1', 'Unknown': '#DCDCDC', 'Germany': '#5F632D', 'Norway': '#8639BD',
    'Netherlands': '#2B39D0', 'Spain': '#155325', 'Switzerland': '#9BE795', 'Canada': '#2A8095',
    'Australia': '#562F8A', 'India': '#73F78E', 'Qatar': '#587984', 'Lebanon': '#3EE78A',
    'Slovenia': '#B4CAD0', 'New Zealand': '#B6754C', 'Ireland': '#2FCC5E', 'Poland': '#C7AC5A',
    'Other STs': '#E6E6FA',
    'ST410': (0.12156862745098039, 0.47058823529411764, 0.7058823529411765, 1.0),
    'ST361': (0.984313725490196, 0.6039215686274509, 0.6, 1.0),
    'ST448': '#4CE60C',
    'ST69': (0.6509803921568628, 0.807843137254902, 0.8901960784313725, 1.0),
    'ST38': (0.2, 0.6274509803921569, 0.17254901960784313, 1.0),
    'ST1722': (0.9921568627450981, 0.7490196078431373, 0.43529411764705883, 1.0),
    'ST131': (0.792156862745098, 0.6980392156862745, 0.8392156862745098, 1.0),
    'ST405': (0.6980392156862745, 0.8745098039215686, 0.5411764705882353, 1.0),
    'ST617': '#D51F6E',
    'ST167': (1.0, 0.4980392156862745, 0.0, 1.0),
    'ST10': (0.41568627450980394, 0.23921568627450981, 0.6039215686274509, 1.0),
    'ST648': '#1B4612',
    'ST1284': '#5D7799',
    'ST940': '#95AC90',
    'ST354': (0.8901960784313725, 0.10196078431372549, 0.10980392156862745, 1.0),
    'ST205': '#D09ED8',
    'ST58': '#3ADEB3',
    'ST127': '#4972BE'
}

# 创建画布 3行6列
fig, axes = plt.subplots(nrows=3, ncols=6, figsize=(22, 10), dpi=300)
plt.subplots_adjust(wspace=0.3, hspace=0.35)

# 循环绘图：行=OXA亚型，列=国家
for row_idx, oxa in enumerate(OXA_variants):
    for col_idx, country in enumerate(target_countries):
        ax = axes[row_idx, col_idx]
        # 筛选当前亚型+当前国家
        subdf = df[(df["variants"] == oxa) & (df["Country"] == country)]
        st_counts = subdf["MLST"].value_counts()
        total_n = len(subdf)

        if total_n == 0:
            # 无样本
            ax.text(0.5, 0.5, "No isolates", ha="center", va="center", transform=ax.transAxes, fontsize=10)
            ax.set_title(f"{country}", fontsize=18)
            ax.axis("off")
            continue

        # 提取对应ST颜色
        pie_colors = [color_dict[st] for st in st_counts.index]
        ax.pie(st_counts.values, labels=None, colors=pie_colors, wedgeprops={"edgecolor":"white", "linewidth":0.4})
        ax.set_title(f"{country}\n(n={total_n})", fontsize=18)
        ax.axis("equal")

plt.tight_layout(rect=[0, 0.08, 1, 0.96])
plt.savefig("OXA_ST_country_pie_3x6.svg", dpi=600)
plt.show()
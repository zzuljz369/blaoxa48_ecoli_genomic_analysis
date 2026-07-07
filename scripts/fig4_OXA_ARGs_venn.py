# coding:utf-8
# Author:Ljz
# date:2025/10/30 10:18
import pandas as pd
df = pd.read_csv(r"ARGs_result.csv")
ARGs_df = df.loc[:,"aac(2')-IIa":"tet(X)"]
non_zero_ratio = (ARGs_df != 0).mean().sort_values(ascending=False)
print(non_zero_ratio.head(20))
non_zero_counts = (ARGs_df != 0).sum().sort_values(ascending=False)
print(non_zero_counts.head(20))

#统计粘菌素耐药基因
mcr_ARGs_df = ARGs_df.filter(like="mcr")
mcr_ARGs_sum = mcr_ARGs_df.sum(axis=1)
print("mcr非0基因组计数",(mcr_ARGs_sum != 0).sum())
non_zero_counts = (mcr_ARGs_df != 0).sum().sort_values(ascending=False)
print("mcr非0统计\n",non_zero_counts)

non_zero_ratio = (mcr_ARGs_df != 0).mean().sort_values(ascending=False)
print("mcr非0统计比例\n",non_zero_ratio)

# 统计三个OXA变异子耐药情况
major_OXA_variants = ["blaOXA-48","blaOXA-181","blaOXA-244"]
major_ARGs_df = df[df["variants"].isin(major_OXA_variants)]

major_ARGs_df = major_ARGs_df.loc[:,"variants":"tet(X)"]
major_ARGs_df = major_ARGs_df.loc[:,major_ARGs_df.sum()!=0]
major_ARGs_df.set_index("variants",inplace=True)
major_ARGs_df.to_csv(r"major_ARGs.csv")
print("三个变异子共有耐药基因数目",len(major_ARGs_df.columns))

#绘制三个变异子耐药基因韦恩图
from venn import venn
import matplotlib.pyplot as plt
# from proplot import rc
from matplotlib.colors import ListedColormap

# rc["font.family"] = "Times New Roman"
# rc["axes.labelsize"] = 24
# rc["axes.labelweight"] = "bold"
# rc["tick.labelsize"] = 24
# rc["suptitle.size"] = 24
# rc["title.size"] = 24
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
OXA48_ARGs_df = major_ARGs_df[major_ARGs_df.index == "blaOXA-48"]
OXA48_ARGs_df = OXA48_ARGs_df.loc[:,OXA48_ARGs_df.sum() != 0]
print("OXA48耐药基因组数目",len(OXA48_ARGs_df.columns))
OXA181_ARGs_df = major_ARGs_df[major_ARGs_df.index == "blaOXA-181"]
OXA181_ARGs_df = OXA181_ARGs_df.loc[:,OXA181_ARGs_df.sum() != 0]
print("OXA181耐药基因组数目",len(OXA181_ARGs_df.columns))
OXA244_ARGs_df = major_ARGs_df[major_ARGs_df.index == "blaOXA-244"]
OXA244_ARGs_df = OXA244_ARGs_df.loc[:,OXA244_ARGs_df.sum() != 0]
print("OXA244耐药基因组数目",len(OXA244_ARGs_df.columns))
venn_dict3 = dict(
set_OXA_48 = set(OXA48_ARGs_df.columns),
set_OXA_181 = set(OXA181_ARGs_df.columns),
set_OXA_244 = set(OXA244_ARGs_df.columns)
)
fig,ax = plt.subplots(figsize=(5,4),dpi=600)
color_dict = {'blaOXA-48': (0.984313725490196, 0.7058823529411765, 0.6823529411764706), 'blaOXA-181': (0.7019607843137254, 0.803921568627451, 0.8901960784313725), 'blaOXA-244': (0.8, 0.9215686274509803, 0.7725490196078432), 'blaOXA-484': (0.8705882352941177, 0.796078431372549, 0.8941176470588236), 'blaOXA-232': (0.996078431372549, 0.8509803921568627, 0.6509803921568628), 'blaOXA-204': (1.0, 1.0, 0.8), 'blaOXA-162': (0.8980392156862745, 0.8470588235294118, 0.7411764705882353), 'blaOXA-245': (0.9921568627450981, 0.8549019607843137, 0.9254901960784314), 'blaOXA-163': (0.9490196078431372, 0.9490196078431372, 0.9490196078431372)}
colors = [color_dict[i] for i in major_OXA_variants]
cmap = ListedColormap(colors)
venn_figure = venn(venn_dict3,
               ax=ax,
               cmap=cmap,
               alpha=0.8)
plt.tight_layout()
plt.savefig(r"fig4_OXA_ARGs_venn.svg",dpi=600)
plt.show()

# 更专业的统计学标注版本
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statannotations.Annotator import Annotator
# from proplot import rc


OXA48_ARGs_sum = OXA48_ARGs_df.sum(axis=1)
OXA181_ARGs_sum = OXA181_ARGs_df.sum(axis=1)
OXA244_ARGs_sum = OXA244_ARGs_df.sum(axis=1)

# 准备数据（长格式）
def prepare_long_format_data(OXA48_ARGs_sum,OXA181_ARGs_sum,OXA244_ARGs_sum):
    """将数据转换为长格式"""
    data_list = []

    # OXA-48
    for value in OXA48_ARGs_sum.values:
        data_list.append({'Group': 'blaOXA-48', 'Value': value})

    # OXA-181
    for value in OXA181_ARGs_sum.values:
        data_list.append({'Group': 'blaOXA-181', 'Value': value})

    # OXA-244
    for value in OXA244_ARGs_sum.values:
        data_list.append({'Group': 'blaOXA-244', 'Value': value})

    return pd.DataFrame(data_list)


# 创建图形
fig, ax = plt.subplots(figsize=(6, 5), dpi=600)

# 准备数据
df = prepare_long_format_data(OXA48_ARGs_sum,OXA181_ARGs_sum,OXA244_ARGs_sum)

# 小提琴图
sns.violinplot(data=df, x='Group', y='Value',
               palette=[color_dict["blaOXA-48"], color_dict["blaOXA-181"], color_dict["blaOXA-244"]],
               ax=ax)

# # 添加数据点
# sns.stripplot(data=df, x='Group', y='Value', color='black', alpha=0.6, size=3, ax=ax)

# 设置比较对
pairs = [("blaOXA-48", "blaOXA-181"),
         ("blaOXA-48", "blaOXA-244"),
         ("blaOXA-181", "blaOXA-244")]

# 使用statannotations添加显著性标注
try:
    from statannotations.Annotator import Annotator

    annotator = Annotator(ax, pairs, data=df, x='Group', y='Value')
    annotator.configure(test='t-test_ind', text_format='star', loc='inside')
    annotator.apply_and_annotate()

except ImportError:
    print("statannotations库未安装，使用简单标注")
    # 手动添加标注
    add_statistical_annotations_simple(ax, df)

# 美化图形
ax.set_xlabel("")
ax.set_ylim(0,70)
ax.set_yticks(np.arange(0, 71, 10))
ax.set_ylabel("Number of ARGs", fontsize=20, fontweight='bold')
ax.grid(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig(r"OXA耐药基因数目小提琴图_专业版.svg", dpi=600)
plt.show()


def add_statistical_annotations_simple(ax, df):
    """简单版本的统计学标注"""
    groups = df['Group'].unique()
    y_max = df['Value'].max()
    y_step = y_max * 0.08

    from itertools import combinations
    from scipy.stats import ttest_ind

    for idx, (group1, group2) in enumerate(combinations(groups, 2)):
        data1 = df[df['Group'] == group1]['Value']
        data2 = df[df['Group'] == group2]['Value']

        t_stat, p_val = ttest_ind(data1, data2)

        if p_val < 0.05:
            # 找到组的位置
            x1 = list(groups).index(group1)
            x2 = list(groups).index(group2)

            # 绘制线条
            y_line = y_max + (idx + 1) * y_step
            ax.plot([x1, x2], [y_line, y_line], color='black', linewidth=2)
            ax.plot([x1, x1], [y_line - y_step * 0.3, y_line], color='black', linewidth=2)
            ax.plot([x2, x2], [y_line - y_step * 0.3, y_line], color='black', linewidth=2)

            # 添加星号
            if p_val < 0.001:
                stars = '**'
            elif p_val < 0.05:
                stars = '*'
            else:
                stars = '*'

            ax.text((x1 + x2) / 2, y_line + y_step * 0.2, stars,
                    ha='center', va='bottom', fontsize=16, fontweight='bold')





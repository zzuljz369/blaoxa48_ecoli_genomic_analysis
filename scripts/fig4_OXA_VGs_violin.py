# coding:utf-8
# Author:Ljz
# date:2026/07/06
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statannotations.Annotator import Annotator
try:
    from proplot import rc
except ImportError:
    import matplotlib as mpl

# ===================== 分层抽样函数：每个OXA亚型内按国家抽样，每国最多20株 =====================
def sample_by_country_per_variant(sub_df, max_n=20, seed=42):
    np.random.seed(seed)
    sample_list = []
    for _, cnt_group in sub_df.groupby("country"):
        if len(cnt_group) > max_n:
            samp = cnt_group.sample(n=max_n, random_state=seed)
        else:
            samp = cnt_group.copy()
        sample_list.append(samp)
    return pd.concat(sample_list, ignore_index=True)

# ===================== 读取数据 + 分亚型重采样 =====================
df = pd.read_csv(r"VFs.csv")
major_OXA_variants = ["blaOXA-48","blaOXA-181","blaOXA-244"]

# 对三种OXA分别分层抽样
sampled_subsets = []
for var in major_OXA_variants:
    var_df = df[df["variants"] == var].copy()
    var_samp = sample_by_country_per_variant(var_df, max_n=20, seed=42)
    sampled_subsets.append(var_samp)
df_sampled = pd.concat(sampled_subsets, ignore_index=True)
print(f"重采样后总样本量：{df_sampled.shape[0]}")

# 全部毒力基因矩阵
VFs_all = df_sampled.loc[:, "aap/aspU":"ykgK/ecpR"]
non_zero_ratio = (VFs_all != 0).mean().sort_values(ascending=False)
print("毒力基因检出比例TOP20：\n", non_zero_ratio.head(20))
non_zero_counts = (VFs_all != 0).sum().sort_values(ascending=False)
print("毒力基因检出数量TOP20：\n", non_zero_counts.head(20))

# 筛选目标OXA，拆分标签与纯数值毒力矩阵
major_all_df = df_sampled[df_sampled["variants"].isin(major_OXA_variants)].copy()
VFs_only = major_all_df.loc[:, "aap/aspU":"ykgK/ecpR"]
# 删除全部为0的毒力基因列
VFs_only = VFs_only.loc[:, VFs_only.sum() != 0]
# 拼接分型标签，设置索引
major_VFs_df = pd.concat([major_all_df[["variants"]], VFs_only], axis=1)
major_VFs_df.set_index("variants", inplace=True)
major_VFs_df.to_csv(r"major_VFs_sampled.csv")
print("三种OXA共有的非零毒力基因数目：", len(major_VFs_df.columns))

# 分亚型提取毒力矩阵并计算单株毒力基因总数
OXA48_VF_df = major_VFs_df[major_VFs_df.index == "blaOXA-48"]
OXA48_VF_df = OXA48_VF_df.loc[:, OXA48_VF_df.sum() != 0]
OXA48_VF_sum = OXA48_VF_df.sum(axis=1)
print("blaOXA-48携带非零毒力基因总数：", len(OXA48_VF_df.columns))

OXA181_VF_df = major_VFs_df[major_VFs_df.index == "blaOXA-181"]
OXA181_VF_df = OXA181_VF_df.loc[:, OXA181_VF_df.sum() != 0]
OXA181_VF_sum = OXA181_VF_df.sum(axis=1)
print("blaOXA-181携带非零毒力基因总数：", len(OXA181_VF_df.columns))

OXA244_VF_df = major_VFs_df[major_VFs_df.index == "blaOXA-244"]
OXA244_VF_df = OXA244_VF_df.loc[:, OXA244_VF_df.sum() != 0]
OXA244_VF_sum = OXA244_VF_df.sum(axis=1)
print("blaOXA-244携带非零毒力基因总数：", len(OXA244_VF_df.columns))

# 配色
color_dict = {
    'blaOXA-48': (0.984313725490196, 0.7058823529411765, 0.6823529411764706),
    'blaOXA-181': (0.7019607843137254, 0.803921568627451, 0.8901960784313725),
    'blaOXA-244': (0.8, 0.9215686274509803, 0.7725490196078432)
}

# 字体全局配置
try:
    rc["font.family"] = "Times New Roman"
    rc["axes.labelsize"] = 24
    rc["axes.labelweight"] = "bold"
    rc["tick.labelsize"] = 24
except:
    mpl.rcParams["font.family"] = "Times New Roman"
    mpl.rcParams["axes.labelsize"] = 24
    mpl.rcParams["axes.labelweight"] = "bold"
    mpl.rcParams["xtick.labelsize"] = 12
    mpl.rcParams["ytick.labelsize"] = 12
    mpl.rcParams['svg.fonttype'] = 'none'
    mpl.rcParams['pdf.fonttype'] = 42

# 转换长格式绘图数据
def prepare_long_format_data(sum48, sum181, sum244):
    data_list = []
    for val in sum48.values:
        data_list.append({"Group": "blaOXA-48", "Value": val})
    for val in sum181.values:
        data_list.append({"Group": "blaOXA-181", "Value": val})
    for val in sum244.values:
        data_list.append({"Group": "blaOXA-244", "Value": val})
    return pd.DataFrame(data_list)

# 无statannotations时备用显著性标注函数
def add_statistical_annotations_simple(ax, df_plot):
    groups = df_plot["Group"].unique()
    y_max = df_plot["Value"].max()
    y_step = y_max * 0.08
    from itertools import combinations
    for idx, (g1, g2) in enumerate(combinations(groups, 2)):
        d1 = df_plot[df_plot["Group"] == g1]["Value"]
        d2 = df_plot[df_plot["Group"] == g2]["Value"]
        _, p = stats.ttest_ind(d1, d2)
        if p < 0.05:
            x1 = list(groups).index(g1)
            x2 = list(groups).index(g2)
            y_line = y_max + (idx + 1) * y_step
            ax.plot([x1, x2], [y_line, y_line], c="black", lw=2)
            ax.plot([x1, x1], [y_line - y_step*0.3, y_line], c="black", lw=2)
            ax.plot([x2, x2], [y_line - y_step*0.3, y_line], c="black", lw=2)
            if p < 0.001:
                star = "***"
            elif p < 0.01:
                star = "**"
            else:
                star = "*"
            ax.text((x1+x2)/2, y_line + y_step*0.2, star, ha="center", va="bottom", fontsize=16, weight="bold")

# ===================== 绘制毒力基因数量小提琴图 =====================
fig, ax = plt.subplots(figsize=(6, 5), dpi=600)
df_plot = prepare_long_format_data(OXA48_VF_sum, OXA181_VF_sum, OXA244_VF_sum)

sns.violinplot(
    data=df_plot,
    x="Group",
    y="Value",
    palette=[color_dict["blaOXA-48"], color_dict["blaOXA-181"], color_dict["blaOXA-244"]],
    ax=ax
)

# 三组两两比较
pairs = [
    ("blaOXA-48", "blaOXA-181"),
    ("blaOXA-48", "blaOXA-244"),
    ("blaOXA-181", "blaOXA-244")
]

try:
    annotator = Annotator(ax, pairs, data=df_plot, x="Group", y="Value")
    annotator.configure(test="t-test_ind", text_format="star", loc="inside")
    annotator.apply_and_annotate()
except ImportError:
    print("未安装statannotations，启用简易t检验标注")
    add_statistical_annotations_simple(ax, df_plot)

# 坐标轴美化
ax.set_xlabel("")
ax.set_ylim(0, 150)
ax.set_yticks(np.arange(0, 151, 30))
ax.set_ylabel("Number of VGs", fontsize=20, fontweight="bold")
ax.grid(False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

plt.tight_layout()
plt.savefig(r"fig4_OXA_VGs_violin.svg", dpi=600, format="svg")
plt.show()
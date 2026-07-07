# 更专业的统计学标注版本
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statannotations.Annotator import Annotator

# -------------------------- 1. 读取原始数据 --------------------------
df = pd.read_csv(r"ARGs_result.csv")
# 目标OXA亚型
major_OXA_variants = ["blaOXA-48","blaOXA-181","blaOXA-244"]

# -------------------------- 分层抽样函数：单亚型内按国家抽样，每国最多20株 --------------------------
def sample_by_country_per_variant(sub_df, max_n=20, seed=42):
    """
    输入某一种OXA亚型的全部菌株，按country分层抽样
    每个国家随机最多取max_n个，不足则全保留
    """
    np.random.seed(seed)
    sample_list = []
    for _, cnt_group in sub_df.groupby("country"):
        if len(cnt_group) > max_n:
            samp = cnt_group.sample(n=max_n, random_state=seed)
        else:
            samp = cnt_group.copy()
        sample_list.append(samp)
    return pd.concat(sample_list, ignore_index=True)

# 分别对三种OXA亚型分层重采样
sampled_subsets = []
for var in major_OXA_variants:
    var_df = df[df["variants"] == var].copy()
    var_samp = sample_by_country_per_variant(var_df, max_n=20, seed=42)
    sampled_subsets.append(var_samp)
# 合并重采样后全部数据
df_sampled = pd.concat(sampled_subsets, ignore_index=True)
print(f"重采样后总样本量：{df_sampled.shape[0]}")

# -------------------------- ARGs矩阵提取 --------------------------
ARGs_df = df_sampled.loc[:,"aac(2')-IIa":"tet(X)"]
non_zero_ratio = (ARGs_df != 0).mean().sort_values(ascending=False)
print(non_zero_ratio.head(20))
non_zero_counts = (ARGs_df != 0).sum().sort_values(ascending=False)
print(non_zero_counts.head(20))

# 筛选目标OXA，提取耐药基因矩阵
major_ARGs_df = df_sampled[df_sampled["variants"].isin(major_OXA_variants)]
major_ARGs_df = major_ARGs_df.loc[:,"variants":"tet(X)"]
major_ARGs_df = major_ARGs_df.loc[:,major_ARGs_df.sum()!=0]
major_ARGs_df.set_index("variants",inplace=True)
major_ARGs_df.to_csv(r"major_ARGs_sampled.csv")
print("三个变异子共有耐药基因数目",len(major_ARGs_df.columns))

# 分亚型提取耐药矩阵并计算单株耐药基因总数
OXA48_ARGs_df = major_ARGs_df[major_ARGs_df.index == "blaOXA-48"]
OXA48_ARGs_df = OXA48_ARGs_df.loc[:,OXA48_ARGs_df.sum() != 0]
OXA48_ARGs_sum = OXA48_ARGs_df.sum(axis=1)
print("OXA48耐药基因组数目",len(OXA48_ARGs_df.columns))

OXA181_ARGs_df = major_ARGs_df[major_ARGs_df.index == "blaOXA-181"]
OXA181_ARGs_df = OXA181_ARGs_df.loc[:,OXA181_ARGs_df.sum() != 0]
OXA181_ARGs_sum = OXA181_ARGs_df.sum(axis=1)
print("OXA181耐药基因组数目",len(OXA181_ARGs_df.columns))

OXA244_ARGs_df = major_ARGs_df[major_ARGs_df.index == "blaOXA-244"]
OXA244_ARGs_df = OXA244_ARGs_df.loc[:,OXA244_ARGs_df.sum() != 0]
OXA244_ARGs_sum = OXA244_ARGs_df.sum(axis=1)
print("OXA244耐药基因组数目",len(OXA244_ARGs_df.columns))

# 配色字典
color_dict = {'blaOXA-48': (0.984313725490196, 0.7058823529411765, 0.6823529411764706),
              'blaOXA-181': (0.7019607843137254, 0.803921568627451, 0.8901960784313725),
              'blaOXA-244': (0.8, 0.9215686274509803, 0.7725490196078432),
              'blaOXA-484': (0.8705882352941177, 0.796078431372549, 0.8941176470588236),
              'blaOXA-232': (0.996078431372549, 0.8509803921568627, 0.6509803921568628),
              'blaOXA-204': (1.0, 1.0, 0.8),
              'blaOXA-162': (0.8980392156862745, 0.8470588235294118, 0.7411764705882353),
              'blaOXA-245': (0.9921568627450981, 0.8549019607843137, 0.9254901960784314),
              'blaOXA-163': (0.9490196078431372, 0.9490196078431372, 0.9490196078431372)}

# 全局字体配置
import matplotlib as mpl
mpl.rcParams['svg.hashsalt'] = 'optimize'
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["font.family"] = "Times New Roman"
mpl.rcParams["axes.labelsize"] = 24
mpl.rcParams["axes.labelweight"] = "bold"
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12
mpl.rcParams["figure.titlesize"] = 24
mpl.rcParams["axes.titlesize"] = 24

# 数据转长格式
def prepare_long_format_data(OXA48_ARGs_sum,OXA181_ARGs_sum,OXA244_ARGs_sum):
    data_list = []
    for value in OXA48_ARGs_sum.values:
        data_list.append({'Group': 'blaOXA-48', 'Value': value})
    for value in OXA181_ARGs_sum.values:
        data_list.append({'Group': 'blaOXA-181', 'Value': value})
    for value in OXA244_ARGs_sum.values:
        data_list.append({'Group': 'blaOXA-244', 'Value': value})
    return pd.DataFrame(data_list)

# 无statannotations备用显著性标注
def add_statistical_annotations_simple(ax, df):
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
            x1 = list(groups).index(group1)
            x2 = list(groups).index(group2)
            y_line = y_max + (idx + 1) * y_step
            ax.plot([x1, x2], [y_line, y_line], color='black', linewidth=2)
            ax.plot([x1, x1], [y_line - y_step * 0.3, y_line], color='black', linewidth=2)
            ax.plot([x2, x2], [y_line - y_step * 0.3, y_line], color='black', linewidth=2)
            if p_val < 0.001:
                stars = '***'
            elif p_val < 0.01:
                stars = '**'
            elif p_val < 0.05:
                stars = '*'
            else:
                stars = 'ns'
            ax.text((x1 + x2) / 2, y_line + y_step * 0.2, stars,
                    ha='center', va='bottom', fontsize=16, fontweight='bold')

# -------------------------- 小提琴图绘图 --------------------------
fig, ax = plt.subplots(figsize=(6, 5), dpi=600)
df_plot = prepare_long_format_data(OXA48_ARGs_sum,OXA181_ARGs_sum,OXA244_ARGs_sum)

sns.violinplot(
    data=df_plot,
    x='Group',
    y='Value',
    palette=[color_dict["blaOXA-48"], color_dict["blaOXA-181"], color_dict["blaOXA-244"]],
    ax=ax
)

pairs = [("blaOXA-48", "blaOXA-181"),
         ("blaOXA-48", "blaOXA-244"),
         ("blaOXA-181", "blaOXA-244")]

try:
    annotator = Annotator(ax, pairs, data=df_plot, x='Group', y='Value')
    annotator.configure(test='t-test_ind', text_format='star', loc='inside')
    annotator.apply_and_annotate()
except ImportError:
    print("statannotations库未安装，使用简易t检验标注")
    add_statistical_annotations_simple(ax, df_plot)

ax.set_xlabel("")
ax.set_ylim(0,70)
ax.set_yticks(np.arange(0, 71, 10))
ax.set_ylabel("Number of ARGs", fontsize=20, fontweight='bold')
ax.grid(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig(r"fig4_OXA_ARGs_violin.svg", dpi=600)
plt.show()
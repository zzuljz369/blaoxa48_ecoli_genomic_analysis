# coding:utf-8
# Author:Ljz
# date:2025/10/24 16:53
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib.ticker import PercentFormatter
matplotlib.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.labelsize"] = 11
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["xtick.labelsize"] = 11
plt.rcParams["ytick.labelsize"] = 11
plt.rcParams["figure.titlesize"] = 11
plt.rcParams["axes.titlesize"] = 11
df = pd.read_csv(r"../dataset_4886_ecoli_blaOXA_metadata.csv")
df = df[df["MLST"] != "Unknown"]
df = df[df["Year"] != "2025"]
df = df[df["Year"] != "Unknown"]
OXA_variants = ["blaOXA-48","blaOXA-181","blaOXA-244"]
STs = ["ST10","ST1284","ST131","ST167","ST354","ST361","ST38","ST405","ST410","ST69","ST940"]
df['MLST'] = "ST" + df["MLST"]

df['MLST'] = np.where(df['MLST'].isin(STs),
                         df['MLST'],
                         'Other STs')
df["Year"] = df["Year"].apply(lambda x: x if int(x) > 2013 else 2013)
color_dict = {'blaOXA-48': (0.984313725490196, 0.7058823529411765, 0.6823529411764706), 'blaOXA-181': (0.7019607843137254, 0.803921568627451, 0.8901960784313725), 'blaOXA-244': (0.8, 0.9215686274509803, 0.7725490196078432), 'blaOXA-484': (0.8705882352941177, 0.796078431372549, 0.8941176470588236), 'blaOXA-232': (0.996078431372549, 0.8509803921568627, 0.6509803921568628), 'blaOXA-204': (1.0, 1.0, 0.8), 'blaOXA-162': (0.8980392156862745, 0.8470588235294118, 0.7411764705882353), 'blaOXA-245': (0.9921568627450981, 0.8549019607843137, 0.9254901960784314), 'blaOXA-163': (0.9490196078431372, 0.9490196078431372, 0.9490196078431372), 'USA': '#E995E6', 'France': '#BAEE8C', 'Others': '#D2D84C', 'United Kingdom': '#5D05E8', 'Singapore': '#89CFB1', 'Unknown': '#DCDCDC', 'Germany': '#5F632D', 'Norway': '#8639BD', 'Netherlands': '#2B39D0', 'Spain': '#155325', 'Switzerland': '#9BE795', 'Canada': '#2A8095', 'Australia': '#562F8A', 'India': '#73F78E', 'Qatar': '#587984', 'Lebanon': '#3EE78A', 'Slovenia': '#B4CAD0', 'New Zealand': '#B6754C', 'Ireland': '#2FCC5E', 'Poland': '#C7AC5A', 'Other STs': '#E6E6FA', 'ST410': (0.12156862745098039, 0.47058823529411764, 0.7058823529411765, 1.0), 'ST361': (0.984313725490196, 0.6039215686274509, 0.6, 1.0), 'ST448': '#4CE60C', 'ST69': (0.6509803921568628, 0.807843137254902, 0.8901960784313725, 1.0), 'ST38': (0.2, 0.6274509803921569, 0.17254901960784313, 1.0), 'ST1722': (0.9921568627450981, 0.7490196078431373, 0.43529411764705883, 1.0), 'ST131': (0.792156862745098, 0.6980392156862745, 0.8392156862745098, 1.0), 'ST405': (0.6980392156862745, 0.8745098039215686, 0.5411764705882353, 1.0), 'ST617': '#D51F6E', 'ST167': (1.0, 0.4980392156862745, 0.0, 1.0), 'ST10': (0.41568627450980394, 0.23921568627450981, 0.6039215686274509, 1.0), 'ST648': '#1B4612', 'ST1284': '#5D7799', 'ST940': '#95AC90', 'ST354': (0.8901960784313725, 0.10196078431372549, 0.10980392156862745, 1.0), 'ST205': '#D09ED8', 'ST58': '#3ADEB3', 'ST127': '#4972BE'}

ST_order = STs + ["Other STs"]

fig = plt.figure(figsize=(15, 5))
plt.subplots_adjust(wspace=0.3, hspace=0.4)

i = 1
for OXA_variant in OXA_variants:
    OXA_df = df[df['variants'] == OXA_variant]
    result = (OXA_df.groupby('Year')['MLST']
              .value_counts(normalize=True)
              .mul(100)
              .round(2)
              .reset_index(name='percentage'))

    result = result.pivot(index='Year', columns='MLST', values='percentage')
    result.fillna(0, inplace=True)
    result = result.reindex(columns=ST_order, fill_value=0)
    result.to_csv(OXA_variant + '.csv')
    # 将年份索引转换为数值类型
    years = result.index.astype(int)  # 确保年份是数值类型
    x_positions = np.arange(len(years))  # 创建数值型x轴位置

    bottom_y = np.zeros(len(result.index))

    ax = fig.add_subplot(1, 3, i)

    for item in ST_order:
        ax.bar(x_positions, result[item].values, bottom=bottom_y,
               color=color_dict[item], label=item, edgecolor='white', linewidth=0.5)
        bottom_y += result[item].values

    # 设置x轴刻度和标签
    ax.set_xticks(x_positions)
    ax.set_xticklabels(years,rotation=60)  # 使用原始年份作为标签
    ax.set_title(f'{OXA_variant}', fontsize=12, pad=10)
    ax.set_xlabel('Year', fontsize=10)
    ax.set_ylabel('Percentage (%)', fontsize=10)
    # ax.grid(axis='y', linestyle='--', alpha=0.3)
    ax.grid(False)
    ax.set_ylim(0, 100)
    ax.tick_params("x",which="minor",length=0,width=0)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=100))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if i == 1:
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, title='MLST Types',
                   # bbox_to_anchor=(1.05, 0.5),
                   loc='lower center',
                   ncol=6,
                   fontsize=20, frameon=False)


    i += 1

plt.tight_layout()
plt.savefig("fig3_OXA_each_gene_st_trend.svg",dpi=600)
plt.show()




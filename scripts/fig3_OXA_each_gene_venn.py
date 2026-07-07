# coding:utf-8
# Author:Ljz
# Date:2025/7/30
import pandas
import pandas as pd
from venn import venn
import matplotlib.pyplot as plt
from proplot import rc
from matplotlib.colors import ListedColormap

rc["font.family"] = "Times New Roman"
rc["axes.labelsize"] = 24
rc["axes.labelweight"] = "bold"
rc["tick.labelsize"] = 24
rc["suptitle.size"] = 24
rc["title.size"] = 24
df = pd.read_excel(r"dataset_4886_ecoli_blaOXA_metadata.xlsx")
print(df.columns)
OXA_variants = ["blaOXA-48","blaOXA-181","blaOXA-244"]
df = df[df["variants"].isin(OXA_variants)]
df = df[df["MLST"] != "Unknown"]

venn_dict3 = dict(
set_OXA_48 = set(df["MLST"][df["variants"] == "blaOXA-48"].values),
set_OXA_181 = set(df["MLST"][df["variants"] == "blaOXA-181"].values),
set_OXA_244 = set(df["MLST"][df["variants"] == "blaOXA-244"].values)
)

fig,ax = plt.subplots(figsize=(5,4),dpi=600)
color_dict = {'blaOXA-48': (0.984313725490196, 0.7058823529411765, 0.6823529411764706), 'blaOXA-181': (0.7019607843137254, 0.803921568627451, 0.8901960784313725), 'blaOXA-244': (0.8, 0.9215686274509803, 0.7725490196078432), 'blaOXA-484': (0.8705882352941177, 0.796078431372549, 0.8941176470588236), 'blaOXA-232': (0.996078431372549, 0.8509803921568627, 0.6509803921568628), 'blaOXA-204': (1.0, 1.0, 0.8), 'blaOXA-162': (0.8980392156862745, 0.8470588235294118, 0.7411764705882353), 'blaOXA-245': (0.9921568627450981, 0.8549019607843137, 0.9254901960784314), 'blaOXA-163': (0.9490196078431372, 0.9490196078431372, 0.9490196078431372)}
colors = [color_dict[i] for i in OXA_variants]
cmap = ListedColormap(colors)
venn_figure = venn(venn_dict3,
               ax=ax,
               cmap=cmap,
               alpha=0.8)
plt.tight_layout()
plt.savefig(r"venn.svg",dpi=600)
plt.show()

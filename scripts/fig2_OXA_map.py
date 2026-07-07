# -*- coding:utf-8 -*-
"""
https://blog.csdn.net/weixin_40340586/article/details/120576295
��װ����ο�https://blog.csdn.net/Yichar/article/details/109072367  �� https://blog.csdn.net/weixin_41608080/article/details/114494953
http://www.syrr.cn/news/33791.html?action=onClick
"""
#世界地图数据下载地址 https://img.hcharts.cn/mapdata/
#重点参考https://github.com/DingWB/folium_documentation
import matplotlib.colors
import matplotlib.pyplot as plt
from matplotlib import patches
import pandas as pd
import matplotlib as mpl
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import copy
from matplotlib.collections import PatchCollection
import shapely
import seaborn as sns
import matplotlib
font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 12}
mpl.rc('font',**font)
import numpy as np
fig = plt.figure(figsize=(10,8),dpi=600)
ax = fig.add_subplot(111)
# 加载内置的世界地图数据
world = gpd.read_file("世界国家.shp")
# 查看数据的前几行
print(world.head())
world.crs = "epsg:4326"
world.to_crs(crs="epsg:3857")

#读取数据，并对名称进行修正
country_data = pd.read_csv(r"dataset_4886_ecoli_blaOXA_metadata.csv",encoding="gbk")
country_data = country_data.loc[country_data["Country"] != "Unknown",'Country']
country_counts = country_data.value_counts()
# world = pd.merge(world,species_map,on="name",how="left")
print(world.columns)
world["NAME"] = world["NAME"].apply(lambda x:x.lower().title() if not pd.isna(x) else np.nan)
world.NAME.replace({'Russian Federation':'Russia','Cote D¡¯Ivoire':"Cote D'Ivoire'",
                  "Korea,Democratic People'S Republic Of":'North Korea',
                   'Macedonia,The Former Yugoslav Republic Of':'North Macedonia',
                   'Syrian Arab Republic':'Syria','Korea, Republic Of':'South Korea',
                   'Congo,The Democratic Republic Of The':'Dem. Rep. Congo',
                    'United States':'USA'
                  },inplace=True)

country_counts = country_counts.reset_index(name="size")
print(country_counts)
differences = [ x for x in list(country_counts["Country"].values) if x not in list(world["NAME"].values)]
world = world.merge(country_counts,how="left",left_on=["NAME"],right_on=["Country"])
world = world[world["NAME"] != "Antarctica"]
# line_geoms = world.geometry.boundary[world.geometry.boundary.geom_type == 'LineString']
# ax = line_geoms.plot(edgecolor="black", linewidth=0.5)
ax = world.geometry.boundary.plot(edgecolor="black",
                                  linewidth=0.3)


new_colors = matplotlib.colormaps["Blues"]
new_colors = new_colors(np.linspace(0.3,1,256))
new_colors = colors.ListedColormap(new_colors)
# https://blog.csdn.net/weixin_46090057/article/details/119882021
# def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=624):
#     new_cmap = colors.LinearSegmentedColormap.from_list(
#         "trunc({n},{a:.2f},{b:.2f})".format(n=cmap.name, a=minval, b=maxval),
#         cmap(np.linspace(minval, maxval, n)),
#     )
#     return new_cmap
# cmap = truncate_colormap(cmap,0.3,1.0)

# https://zhuanlan.zhihu.com/p/536879741

ax = world.plot(linewidth=0.3,
                edgecolor="black",
                ax=ax,
                facecolor="white",
                column="size",
                cmap=new_colors,
                legend=False,

                )
# 获取颜色映射的范围
vmin = world['size'].min()  # 数据最小值
vmax = world['size'].max()  # 数据最大值

# 创建一个 ScalarMappable 对象，用于与 colorbar 关联
import matplotlib.cm as cm
sm = cm.ScalarMappable(cmap=new_colors, norm=plt.Normalize(vmin=vmin, vmax=1000))
sm.set_array([])  # 必须调用 set_array，即使没有数据传入

# 手动添加 colorbar
cbar = fig.colorbar(sm, ax=ax, orientation="vertical")
# 设置刻度值和刻度标签
tick_values = [0, 200, 400, 600, 800, 1000]  # 刻度值
tick_labels = ["0", "200", "400", "600", "800", "1000"]  # 刻度标签
cbar.set_ticks(tick_values)  # 设置刻度值
cbar.set_ticklabels(tick_labels)  # 设置刻度标签
# 设置 colorbar 标题
cbar.set_label("Size (0-1000)")

ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])
plt.tight_layout()
plt.savefig(r"OXA_all.svg",format="svg",dpi=600)
plt.show()

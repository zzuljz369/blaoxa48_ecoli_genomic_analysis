# -*- coding: utf-8 -*-
"""
Created on 2025/7/12
@author: ljzzz
"""
"""对每个基因绘制地图"""
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
#读取数据，并对名称进行修正
country_data = pd.read_csv(r"dataset_4886_ecoli_blaOXA_metadata",encoding="gbk")
def OXA_map(OXA_variant,country_data):
    fig = plt.figure(figsize=(10,8),dpi=600)
    ax = fig.add_subplot(111)
    # 加载内置的世界地图数据
    world = gpd.read_file("世界国家.shp")
    # 查看数据的前几行
    print(world.head())
    world.crs = "epsg:4326"
    world.to_crs(crs="epsg:3857")
    country_data = country_data.loc[country_data["variants"] == OXA_variant]
    country_data = country_data.loc[country_data["Country"] != "Unknown",'Country']
    country_counts = country_data.value_counts()
    world["NAME"] = world["NAME"].apply(lambda x:x.lower().title() if not pd.isna(x) else np.nan)
    world.NAME.replace({'Russian Federation':'Russia','Cote D¡¯Ivoire':"Cote D'Ivoire'",
                      "Korea,Democratic People'S Republic Of":'North Korea',
                       'Macedonia,The Former Yugoslav Republic Of':'North Macedonia',
                       'Syrian Arab Republic':'Syria','Korea, Republic Of':'South Korea',
                       'Congo,The Democratic Republic Of The':'Dem. Rep. Congo',
                        'United States':'USA'
                      },inplace=True)

    country_counts = country_counts.reset_index(name="size")
    differences = [ x for x in list(country_counts["index"].values) if x not in list(world["NAME"].values)]
    world = world.merge(country_counts,how="left",left_on=["NAME"],right_on=["index"])
    world = world[world["NAME"] != "Antarctica"]
    # line_geoms = world.geometry.boundary[world.geometry.boundary.geom_type == 'LineString']
    # ax = line_geoms.plot(edgecolor="black", linewidth=0.5)
    ax = world.geometry.boundary.plot(edgecolor="black",
                                      linewidth=0.3)
    size_max = world["size"].max()
    stop = 0.3+ 0.7/1000*size_max
    new_colors = matplotlib.colormaps["Blues"]
    new_colors = new_colors(np.linspace(0.3,stop,256))
    new_colors = colors.ListedColormap(new_colors)


    ax = world.plot(linewidth=0.3,
                    edgecolor="black",
                    ax=ax,
                    facecolor="white",
                    column="size",
                    cmap=new_colors,
                    legend=False,
                    )
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(OXA_variant)
    tiff_name = OXA_variant + ".tif"
    plt.tight_layout()
    plt.savefig(tiff_name,format="tif",dpi=600)
    plt.show()

item = "blaOXA-484"
OXA_map(item,country_data)


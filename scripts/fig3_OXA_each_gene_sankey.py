# -*- coding:utf-8 -*-
import matplotlib
import matplotlib.pyplot as plt
from pysankey2 import Sankey
import pandas as pd
import numpy as np
matplotlib.rcParams['font.family'] = 'Times New Roman '
A = pd.read_csv(r"sankey.csv")
df = A.sort_values(by='Layer1',ascending=False)
color_dict = {'blaOXA-48': (0.984313725490196, 0.7058823529411765, 0.6823529411764706), 'blaOXA-181': (0.7019607843137254, 0.803921568627451, 0.8901960784313725), 'blaOXA-244': (0.8, 0.9215686274509803, 0.7725490196078432), 'blaOXA-484': (0.8705882352941177, 0.796078431372549, 0.8941176470588236), 'blaOXA-232': (0.996078431372549, 0.8509803921568627, 0.6509803921568628), 'blaOXA-204': (1.0, 1.0, 0.8), 'blaOXA-162': (0.8980392156862745, 0.8470588235294118, 0.7411764705882353), 'blaOXA-245': (0.9921568627450981, 0.8549019607843137, 0.9254901960784314), 'blaOXA-163': (0.9490196078431372, 0.9490196078431372, 0.9490196078431372)}



#随机生成颜色
def random_color(number):
  color = []
  intnum = [str(x) for x in np.arange(10)]
  #Out[138]: ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
  alphabet = [chr(x) for x in (np.arange(6) + ord('A'))]
  #Out[139]: ['A', 'B', 'C', 'D', 'E', 'F']
  colorArr = np.hstack((intnum, alphabet))
  #Out[142]: array(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C','D', 'E', 'F'], dtype='<U1')
  print(colorArr)
  for j in range(number):
    color_single = '#'
    for i in range(6):
      index = np.random.randint(len(colorArr))
      color_single += colorArr[index]
    #Out[148]: '#EDAB33'
    color.append(color_single)
  return color

lst_all = list(df["Layer3"].unique())  + list(df["Layer2"].unique())
colors = random_color(len(lst_all))
for item,color in zip(lst_all,colors):
    color_dict[item] = color

color_dict = {'blaOXA-48': (0.984313725490196, 0.7058823529411765, 0.6823529411764706), 'blaOXA-181': (0.7019607843137254, 0.803921568627451, 0.8901960784313725), 'blaOXA-244': (0.8, 0.9215686274509803, 0.7725490196078432), 'blaOXA-484': (0.8705882352941177, 0.796078431372549, 0.8941176470588236), 'blaOXA-232': (0.996078431372549, 0.8509803921568627, 0.6509803921568628), 'blaOXA-204': (1.0, 1.0, 0.8), 'blaOXA-162': (0.8980392156862745, 0.8470588235294118, 0.7411764705882353), 'blaOXA-245': (0.9921568627450981, 0.8549019607843137, 0.9254901960784314), 'blaOXA-163': (0.9490196078431372, 0.9490196078431372, 0.9490196078431372), 'USA': '#E995E6', 'France': '#BAEE8C', 'Others': '#D2D84C', 'United Kingdom': '#5D05E8', 'Singapore': '#89CFB1', 'Unknown': '#956E59', 'Germany': '#5F632D', 'Norway': '#8639BD', 'Netherlands': '#2B39D0', 'Spain': '#155325', 'Switzerland': '#9BE795', 'Canada': '#2A8095', 'Australia': '#562F8A', 'India': '#73F78E', 'Qatar': '#587984', 'Lebanon': '#3EE78A', 'Slovenia': '#B4CAD0', 'New Zealand': '#B6754C', 'Ireland': '#2FCC5E', 'Poland': '#C7AC5A', 'Other STs': '#7E769D', 'ST410': '#B75D2C', 'ST361': '#920837', 'ST448': '#4CE60C', 'ST69': '#72CB38', 'ST38': '#E37DED', 'ST1722': '#822F7A', 'ST131': '#1542CE', 'ST405': '#AA6747', 'ST617': '#D51F6E', 'ST167': '#9901A6', 'ST10': '#C03C13', 'ST648': '#1B4612', 'ST1284': '#5D7799', 'ST940': '#95AC90', 'ST354': '#40D3FF', 'ST205': '#D09ED8', 'ST58': '#3ADEB3', 'ST127': '#4972BE'}
ColorDict = {'ST69': (0.6509803921568628, 0.807843137254902, 0.8901960784313725, 1.0), 'ST410': (0.12156862745098039, 0.47058823529411764, 0.7058823529411765, 1.0), 'ST405': (0.6980392156862745, 0.8745098039215686, 0.5411764705882353, 1.0), 'ST38': (0.2, 0.6274509803921569, 0.17254901960784313, 1.0), 'ST361': (0.984313725490196, 0.6039215686274509, 0.6, 1.0), 'ST354': (0.8901960784313725, 0.10196078431372549, 0.10980392156862745, 1.0), 'ST1722': (0.9921568627450981, 0.7490196078431373, 0.43529411764705883, 1.0), 'ST167': (1.0, 0.4980392156862745, 0.0, 1.0), 'ST131': (0.792156862745098, 0.6980392156862745, 0.8392156862745098, 1.0), 'ST10': (0.41568627450980394, 0.23921568627450981, 0.6039215686274509, 1.0), 'Other STs': '#E6E6FA', 'Unknown': '#DCDCDC'}
for key,values in ColorDict.items():
  if key in color_dict:
    color_dict[key] = ColorDict[key]
sky_auto_global_colors = Sankey(df,colorDict=color_dict,colorMode="global",stripColor="left")
print(sky_auto_global_colors.colorDict)
print(colors)

fig,ax = sky_auto_global_colors.plot(figSize=(12, 8),fontSize=10,boxWidth=.6,text_kws={"family":"Times New Roman","weight":"bold"})
plt.tight_layout()
plt.savefig(r"sankey.svg",format="svg",dpi=600,bbox_inches="tight")
plt.show()

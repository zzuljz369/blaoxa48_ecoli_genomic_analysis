# coding:utf-8
# Author:Ljz
# date:2025/11/2 23:55
# coding:utf-8
# Author:Ljz
# date:2025/11/2 23:35

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity
from scipy.spatial.distance import pdist, squareform
from scipy import stats
import matplotlib.patches as patches
import warnings
# from proplot import rc
import matplotlib as mpl
warnings.filterwarnings('ignore')

# 设置中文字体和图形样式
# mpl.rcParams["font.family"] = "Times New Roman"
# mpl.rcParams["axes.labelsize"] = 24
# mpl.rcParams["axes.labelweight"] = "bold"
# mpl.rcParams["tick.labelsize"] = 24
# mpl.rcParams["suptitle.size"] = 24
# mpl.rcParams["title.size"] = 24
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

class PCoAAnalyzer:
    def __init__(self):
        self.results = {}

    def load_data(self, file_path, group_column='Group', sep=','):
        """加载数据"""
        self.df = pd.read_csv(file_path, sep=sep)
        self.group_column = group_column

        # 分离数值数据和分组信息
        numeric_columns = [col for col in self.df.columns if col != group_column]
        self.numeric_data = self.df[numeric_columns]
        self.groups = self.df[group_column]

        print("数据基本信息:")
        print(f"数据维度: {self.df.shape}")
        print(f"数值数据维度: {self.numeric_data.shape}")
        print(f"分组信息: {self.groups.unique().tolist()}")
        print(f"各组样本数:\n{self.groups.value_counts()}")

        return self.df

    def perform_pcoa(self, distance_metric='bray'):
        """执行PCoA分析"""
        print(f"\n开始PCoA分析，距离度量: {distance_metric}")

        # 计算距离矩阵
        if distance_metric == 'braycurtis':
            from skbio.diversity import beta_diversity
            dist_matrix = beta_diversity(distance_metric, self.numeric_data.values,
                                         self.numeric_data.index.astype(str))
            dist_matrix = dist_matrix.data
        else:
            dist_matrix = pdist(self.numeric_data, metric=distance_metric)

        # 执行PCoA
        self.pcoa_result = pcoa(dist_matrix)

        # 计算方差解释比例
        eigenvalues = self.pcoa_result.eigvals
        variance_explained = (eigenvalues / eigenvalues.sum()) * 100

        # 创建结果数据框
        self.pcoa_df = pd.DataFrame({
            'Sample': self.numeric_data.index,
            'Group': self.groups.values,
            'PCo1': self.pcoa_result.samples.iloc[:, 0],
            'PCo2': self.pcoa_result.samples.iloc[:, 1],
            'PCo3': self.pcoa_result.samples.iloc[:, 2]
        })

        print("PCoA分析结果:")
        print(f"PCo1解释方差: {variance_explained[0]:.2f}%")
        print(f"PCo2解释方差: {variance_explained[1]:.2f}%")
        print(f"PCo3解释方差: {variance_explained[2]:.2f}%")

        self.results['variance_explained'] = variance_explained
        self.results['distance_matrix'] = dist_matrix

        return self.pcoa_df, variance_explained

    def calculate_confidence_ellipse(self, x, y, group, n_std=1.96):
        """计算置信椭圆坐标"""
        if len(x) <= 2:
            return None, None, None, None

        covariance = np.cov(x, y)
        pearson = covariance[0, 1] / np.sqrt(covariance[0, 0] * covariance[1, 1])

        # 椭圆半径
        ell_radius_x = np.sqrt(1 + pearson)
        ell_radius_y = np.sqrt(1 - pearson)

        # 椭圆参数
        theta = np.linspace(0, 2 * np.pi, 100)
        ellipse_coords = np.column_stack([
            ell_radius_x * np.cos(theta),
            ell_radius_y * np.sin(theta)
        ])

        # 椭圆旋转
        rotation = np.array([
            [np.cos(np.pi / 4), np.sin(np.pi / 4)],
            [-np.sin(np.pi / 4), np.cos(np.pi / 4)]
        ])

        ellipse_coords = ellipse_coords @ rotation

        # 缩放和平移
        scale = n_std * np.sqrt(np.diag(covariance))
        ellipse_coords[:, 0] = ellipse_coords[:, 0] * scale[0] + np.mean(x)
        ellipse_coords[:, 1] = ellipse_coords[:, 1] * scale[1] + np.mean(y)

        return ellipse_coords[:, 0], ellipse_coords[:, 1], np.mean(x), np.mean(y)

    def permanova_test(self):
        """PERMANOVA检验组间差异"""
        try:
            from skbio.stats.distance import permanova
            # 需要将距离矩阵转换为DistanceMatrix对象
            from skbio.stats.distance import DistanceMatrix
            dm = DistanceMatrix(squareform(self.results['distance_matrix']))
            result = permanova(dm, self.groups.reset_index(drop=True))
            return result
        except:
            print("PERMANOVA测试需要scikit-bio的完整功能")
            return None

    def plot_pcoa_with_ellipses(self, figsize=(12, 8), show_labels=True, save_path=None):
        """绘制带置信椭圆的PCoA图"""

        # 设置颜色和形状
        unique_groups = self.pcoa_df['Group'].unique()
        n_groups = len(unique_groups)

        # 颜色方案
        if n_groups == 3:
            colors = ['#E41A1C', '#377EB8', '#4DAF4A']  # 红，蓝，绿
            markers = ['o', 'o', 'o']  # 圆形，方形，三角形
        else:
            colors = sns.color_palette("Set1", n_groups)
            markers = ['o'] * n_groups

        # 创建图形
        fig, ax = plt.subplots(figsize=figsize)

        # 为每个组绘制置信椭圆和散点
        ellipses_data = []

        for i, group in enumerate(unique_groups):
            group_data = self.pcoa_df[self.pcoa_df['Group'] == group]
            x, y = group_data['PCo1'], group_data['PCo2']
            color = colors[i]
            marker = markers[i]

            # 绘制散点
            ax.scatter(x, y, c=color, marker=marker, s=80,
                       label=group, alpha=0.8, edgecolors='white', linewidth=1)

            # 计算并绘制置信椭圆
            if len(x) > 2:
                ellipse_x, ellipse_y, center_x, center_y = self.calculate_confidence_ellipse(
                    x, y, group, n_std=1.96
                )
                if ellipse_x is not None:
                    ax.fill(ellipse_x, ellipse_y, color=color, alpha=0.2)
                    ax.plot(ellipse_x, ellipse_y, color=color, alpha=0.6, linewidth=1)

                    # 标记中心点
                    ax.scatter(center_x, center_y, c=color, marker='*', s=150,
                               alpha=0.8, edgecolors='black')

            # 添加样本标签
            if show_labels and len(group_data) <= 20:
                texts = []
                for _, row in group_data.iterrows():
                    texts.append(ax.text(row['PCo1'], row['PCo2'], row['Sample'],
                                         fontsize=8, alpha=0.8))


        # 设置坐标轴标签
        variance_exp = self.results['variance_explained']
        ax.set_xlabel(f'PCo1 ({variance_exp[0]:.2f}%)', fontsize=24, fontweight='bold')
        ax.set_ylabel(f'PCo2 ({variance_exp[1]:.2f}%)', fontsize=24, fontweight='bold')

        # 添加标题和图例
        ax.set_title('PCoA Analysis with 95% Confidence Ellipses',
                     fontsize=16, fontweight='bold', pad=20)

        # 添加网格和图例
        ax.grid(True, alpha=0.3)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        # 设置坐标轴相等
        ax.set_aspect('equal')

        plt.tight_layout()

        # 保存图片
        if save_path:
            plt.savefig(save_path, dpi=600, bbox_inches='tight')
            print(f"图片已保存至: {save_path}")

        plt.show()

        return fig, ax

    def statistical_summary(self):
        """生成统计摘要"""
        print("\n" + "=" * 60)
        print("PCoA分析统计摘要")
        print("=" * 60)

        # 各组中心位置
        group_centers = self.pcoa_df.groupby('Group')[['PCo1', 'PCo2']].mean()
        print("各组中心坐标:")
        print(group_centers)

        # 方差解释
        variance_exp = self.results['variance_explained']
        print(f"\n方差解释比例:")
        print(f"前3个主坐标累计解释方差: {variance_exp[:3].sum():.2f}%")

        # 尝试PERMANOVA检验
        permanova_result = self.permanova_test()
        if permanova_result is not None:
            print(f"\nPERMANOVA检验结果:")
            print(f"F-value: {permanova_result[4]:.4f}")
            print(f"p-value: {permanova_result[5]:.4f}")


# 使用示例
def main():
    # 创建分析器实例
    analyzer = PCoAAnalyzer()

    # 加载数据（请修改为您的文件路径）
    try:
        df = analyzer.load_data('major_VFs.csv',
                                group_column='Group', sep=',')
    except Exception as e:
        print(f"数据加载错误: {e}")
        print("创建示例数据进行演示...")
        # 创建示例数据
        np.random.seed(42)
        n_samples = 30
        n_features = 20

        data = {
            'Sample': [f'S{i + 1}' for i in range(n_samples)],
            'Group': ['Group_A'] * 10 + ['Group_B'] * 10 + ['Group_C'] * 10
        }

        # 添加模拟的特征数据
        for i in range(n_features):
            data[f'Feature_{i + 1}'] = np.random.exponential(1, n_samples)

        df = pd.DataFrame(data)
        df.set_index('Sample', inplace=True)
        analyzer.df = df
        analyzer.numeric_data = df.iloc[:, 1:]
        analyzer.groups = df['Group']
        analyzer.group_column = 'Group'

    # 执行PCoA分析
    pcoa_df, variance_exp = analyzer.perform_pcoa(distance_metric='braycurtis')

    # 绘制图形
    analyzer.plot_pcoa_with_ellipses(
        figsize=(12, 8),
        show_labels=True,
        save_path='VGs_pcoa_plot.svg'
    )

    # 生成统计摘要
    analyzer.statistical_summary()

    # 保存坐标数据
    pcoa_df.to_csv('pcoa_coordinates.csv', index=False)
    print(f"\n坐标数据已保存至: pcoa_coordinates.csv")


if __name__ == "__main__":
    main()
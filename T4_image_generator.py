import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


wide = 0.3
num_all = 223
length_head = 2.86
length_body = 1.65


def draw_rotated_rectangles(domian, t):
    fig, ax = plt.subplots(figsize=(10, 10))

    # 遍历 domian 中的每个矩形
    for rect in domian:
        # 提取矩形的四个顶点坐标
        x1, y1 = rect[0]
        x2, y2 = rect[2]
        x3, y3 = rect[1]
        x4, y4 = rect[3]

        # 创建一个包含四个顶点的多边形 (不规则矩形)
        polygon = Polygon(
            [(x1, y1), (x2, y2), (x3, y3), (x4, y4)], edgecolor="blue", facecolor="none"
        )

        # 添加多边形到图形中
        ax.add_patch(polygon)

    # 设置坐标轴比例
    ax.set_aspect("equal", "box")

    # 反转 y 轴
    ax.invert_yaxis()

    # 设置图形显示范围
    ax.autoscale()

    # 添加标题和标签
    plt.title(f"t = {t-100}")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.savefig(f"result44_test_{t}.png")


def domain_cal(point_ahead_x, point_ahead_y, point_back_x, point_back_y, length):
    return [
        (
            point_ahead_x
            + (0.15 / length) * (point_back_y - point_ahead_y) / 2
            + (0.275 / length) * (point_ahead_x - point_back_x) / 2,
            point_ahead_y
            + (0.15 / length) * (point_ahead_x - point_back_x) / 2
            + (0.275 / length) * (point_ahead_y - point_back_y) / 2,
        ),
        (
            point_back_x
            + (0.15 / length) * (point_ahead_y - point_back_y) / 2
            + (0.275 / length) * (point_back_x - point_ahead_x) / 2,
            point_back_y
            + (0.15 / length) * (point_back_x - point_ahead_x) / 2
            + (0.275 / length) * (point_back_y - point_ahead_y) / 2,
        ),
        (
            point_ahead_x
            - (0.15 / length) * (point_back_y - point_ahead_y) / 2
            + (0.275 / length) * (point_ahead_x - point_back_x) / 2,
            point_ahead_y
            - (0.15 / length) * (point_ahead_x - point_back_x) / 2
            + (0.275 / length) * (point_ahead_y - point_back_y) / 2,
        ),
        (
            point_back_x
            - (0.15 / length) * (point_ahead_y - point_back_y) / 2
            + (0.275 / length) * (point_back_x - point_ahead_x) / 2,
            point_back_y
            - (0.15 / length) * (point_back_x - point_ahead_x) / 2
            + (0.275 / length) * (point_back_y - point_ahead_y) / 2,
        ),
    ]


import pandas as pd

# 读取Excel文件
data = pd.read_excel("result44_test.xlsx", sheet_name=0)

# 提取奇数行和偶数行
# 奇数行 (0, 2, 4,... => 1, 3, 5,...)
xx = data.iloc[1::2].to_numpy()  # 提取奇数行作为 x 坐标
# 偶数行 (1, 3, 5,... => 2, 4, 6,...)
yy = data.iloc[0::2].to_numpy()  # 提取偶数行作为 y 坐标

for t in range(1, 201):
    domain = [[(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)] for _ in range(num_all)]
    if t % 5 == 0:
        for i in range(0, num_all - 1):
            if i == 0:
                domain[i] = domain_cal(
                    xx[i][t], yy[i][t], xx[i + 1][t], yy[i + 1][t], length_head / 2
                )
            else:
                domain[i] = domain_cal(
                    xx[i][t], yy[i][t], xx[i + 1][t], yy[i + 1][t], length_body / 2
                )
        draw_rotated_rectangles(domain, t)

import math
import pprint
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

right_num = 1e-13

wide = 0.3
V_head = 1.0
num_all = 223
length_head = 2.86
length_body = 1.65
a = 1.7 / (2 * math.pi)
theta_start = 90 / 17 * math.pi


def S_calculation(theta):
    return (
        0.5
        * a
        * (
            theta * math.sqrt(1 + theta * theta)
            + math.log(theta + math.sqrt(1 + theta * theta))
        )
    )


def x_y_cal(theta):
    return a * theta * math.cos(theta), a * theta * math.sin(theta)


def theta_cal(S, theta_start):
    low = 0
    high = theta_start
    while high - low > right_num:  # 精度控制
        mid = (low + high) / 2
        if S_calculation(theta_start) - S_calculation(mid) >= S:
            low = mid
        else:
            high = mid
    return high


def theta_cal_mus(S, theta_start):
    low = theta_start
    high = theta_start * 10
    while high - low > right_num:  # 精度控制
        mid = (low + high) / 2
        if S_calculation(mid) - S_calculation(theta_start) >= S:
            high = mid
        else:
            low = mid
    return high


points = [[(0, 0) for _ in range(num_all + 1)] for _ in range(3000)]
Speed = [[0 for _ in range(num_all + 1)] for _ in range(3000 + 1)]
Theta = [[0 for _ in range(num_all + 1)] for _ in range(3000 + 1)]


def distance(theta1, theta2):
    r1 = a * theta1
    r2 = a * theta2
    return math.sqrt((r1 * r1) + (r2 * r2) - 2 * r1 * r2 * math.cos(theta1 - theta2))


def theta_i_cal(length, theta):
    low = theta
    high = theta + 10
    while high - low > right_num:  # 精度控制
        mid = (low + high) / 2
        if distance(mid, theta) >= length:
            high = mid
        else:
            low = mid
    return high


def dis(x1, y1, x2, y2):
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


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

    # 设置图形显示范围
    ax.autoscale()

    # 添加标题和标签
    plt.title(f"t = {t}")
    plt.xlabel("X")
    plt.ylabel("Y")

    # 显示图形
    plt.show()


def domain_cal(rect1, rect2, length):
    point_ahead_x, point_ahead_y = rect1
    point_back_x, point_back_y = rect2
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


def dis_point_line(point, point_1, point_2):
    return abs(
        (point_2[1] - point_1[1]) * point[0]
        - (point_2[0] - point_1[0]) * point[1]
        + point_2[0] * point_1[1]
        - point_2[1] * point_1[0]
    ) / math.sqrt((point_2[1] - point_1[1]) ** 2 + (point_2[0] - point_1[0]) ** 2)


def rectangles_overlap(rect1, rect2):
    for i in range(0, 4):
        if (
            dis_point_line(rect1[i], rect2[3], rect2[0]) <= wide
            and dis_point_line(rect1[i], rect2[1], rect2[2]) <= wide
            and dis_point_line(rect1[i], rect2[2], rect2[0]) <= length_body
            and dis_point_line(rect1[i], rect2[1], rect2[3]) <= length_body
        ):
            return False

    return True


theta_0 = theta_start


def draw_points(points_head_x, points_head_y):
    plt.figure(figsize=(20, 20))
    plt.scatter(points_head_x, points_head_y, color="blue", marker="o")

    # 添加每个点的标号
    for i, (x, y) in enumerate(zip(points_head_x, points_head_y)):
        plt.text(
            x, y, str(i), fontsize=12, ha="right", va="bottom"
        )  # 使用点的索引作为标号

    plt.title("Scatter Plot of Points")
    plt.xlabel("X Coordinates")
    plt.ylabel("Y Coordinates")

    plt.grid()
    plt.show()


def theta_circle_cal(x_0, y_0, x_1, y_1, r):
    if y_1 > y_0:
        return math.acos((x_1 - x_0) / r)
    else:
        return 2 * math.pi - math.acos((x_1 - x_0) / r)


def judge(t):
    domain = [[(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)] for _ in range(num_all)]
    for i in range(0, num_all):
        if i == 0:
            domain[i] = domain_cal(
                points[t][i],
                points[t][i + 1],
                length_head / 2,
            )
        else:
            domain[i] = domain_cal(
                points[t][i],
                points[t][i + 1],
                length_body / 2,
            )

    if t % 10 == 0:
        draw_rotated_rectangles(domain, t)

    for i in range(0, num_all - 3):
        for j in range(i + 2, num_all):
            if (
                rectangles_overlap(domain[j], domain[i]) == False
                or rectangles_overlap(domain[i], domain[j]) == False
            ):
                draw_rotated_rectangles(domain, t)
                return False

    for i in range(0, num_all):
        if Speed[t][i] > 2:
            return False

    return True


def x_y_t(
    mid, theta_0, theta_21, theta_20, x_2r, y_2r, x_r, y_r, theta_31, theta_30, theta
):
    if mid > 0 and mid <= 100:
        theta_this = theta_cal_mus((100 - mid) * V_head, theta_0)
        x_this, y_this = x_y_cal(theta_this)
    elif mid > 100 and mid <= 200:
        theta_2 = theta_21 - (mid - 100) / 100 * theta_20 * V_head
        x_this = x_2r + 2 * r * math.cos(theta_2)
        y_this = y_2r + 2 * r * math.sin(theta_2)
    elif mid > 200 and mid <= 300:
        theta_3 = theta_31 + theta_30 * (mid - 200) * V_head / 100
        x_this = x_r + r * math.cos(theta_3)
        y_this = y_r + r * math.sin(theta_3)
    else:
        x_this, y_this = x_y_cal(theta)
        x_this = -x_this
        y_this = -y_this
    return x_this, y_this


def x_y_find(
    length,
    x_last,
    y_last,
    t,
    theta_0,
    theta_21,
    theta_20,
    x_2r,
    y_2r,
    x_r,
    y_r,
    theta_31,
    theta_30,
    theta,
):
    low = 0
    hight = t
    while hight - low > right_num:
        mid = (hight + low) / 2
        if mid > 0 and mid <= 100:
            theta_this = theta_cal_mus(100 - mid, theta_0)
            x_this, y_this = x_y_cal(theta_this)
        elif mid > 100 and mid <= 200:
            theta_2 = theta_21 - (mid - 100) * V_head / 100 * theta_20
            x_this = x_2r + 2 * r * math.cos(theta_2)
            y_this = y_2r + 2 * r * math.sin(theta_2)
        elif mid > 200 and mid <= 300:
            theta_3 = theta_31 + theta_30 * (mid - 200) * V_head / 100
            x_this = x_r + r * math.cos(theta_3)
            y_this = y_r + r * math.sin(theta_3)
        else:
            x_this, y_this = x_y_cal(theta)
            x_this = -x_this
            y_this = -y_this

        if dis(x_this, y_this, x_last, y_last) >= length:
            low = mid
        else:
            hight = mid

    return x_y_t(
        hight,
        theta_0,
        theta_21,
        theta_20,
        x_2r,
        y_2r,
        x_r,
        y_r,
        theta_31,
        theta_30,
        theta,
    )


low = 0
high = 100


while high - low > right_num:
    V_head = (low + high) / 2

    x_head = [0 for _ in range(30000)]
    y_head = [0 for _ in range(30000)]

    points = [[(0, 0) for _ in range(num_all + 1)] for _ in range(3000)]
    Speed = [[0 for _ in range(num_all + 1)] for _ in range(3000 + 1)]
    Theta = [[0 for _ in range(num_all + 1)] for _ in range(3000 + 1)]

    for t in range(0, 101):

        theta_this = theta_cal_mus((100 - t) * V_head, theta_0)
        x_this, y_this = x_y_cal(theta_this)

        points[t][0] = (x_this, y_this)
        x_head[t] = x_this
        y_head[t] = y_this

        theta_i = theta_i_cal(length_head, theta_this)
        points[t][1] = x_y_cal(theta_i)
        Theta[t][1] = theta_i

        for i in range(1, num_all):
            theta_i = theta_i_cal(length_body, theta_i)
            points[t][i + 1] = x_y_cal(theta_i)
            Theta[t][i + 1] = theta_i

        if judge(t) == 0:
            flag = 0
            break

    k_1 = math.cos(theta_0) - math.sin(theta_0) * theta_0
    k_2 = math.sin(theta_0) + math.cos(theta_0) * theta_0
    k = k_2 / k_1
    d = abs(
        a * theta_0 * math.cos(theta_0) * k - a * theta_0 * math.sin(theta_0)
    ) / math.sqrt(k * k + 1)
    r = a * a * theta_0 * theta_0 / 3.0 / d

    x_2r = a * theta_0 * math.cos(theta_0) - 2 * r * k_2 / math.sqrt(
        k_1 * k_1 + k_2 * k_2
    )
    y_2r = a * theta_0 * math.sin(theta_0) + 2 * r * k_1 / math.sqrt(
        k_1 * k_1 + k_2 * k_2
    )

    x_r = -a * theta_0 * math.cos(theta_0) + r * k_2 / math.sqrt(k_1 * k_1 + k_2 * k_2)
    y_r = -a * theta_0 * math.sin(theta_0) - r * k_1 / math.sqrt(k_1 * k_1 + k_2 * k_2)

    x_in = a * theta_0 * math.cos(theta_0)
    y_in = a * theta_0 * math.sin(theta_0)

    x_qie = x_2r / 3 + 2 * x_r / 3
    y_qie = y_2r / 3 + 2 * y_r / 3

    theta_21 = theta_circle_cal(x_2r, y_2r, x_in, y_in, 2 * r)
    theta_22 = theta_circle_cal(x_2r, y_2r, x_qie, y_qie, 2 * r)

    if theta_21 > theta_22:
        theta_20 = theta_21 - theta_22
    else:
        theta_20 = theta_21 - theta_22 + 2 * math.pi

    t_1 = theta_20 * 2 * r

    theta_31 = theta_circle_cal(x_r, y_r, x_qie, y_qie, r)
    theta_32 = theta_circle_cal(x_r, y_r, -x_in, -y_qie, r)

    if theta_31 > theta_32:
        theta_30 = theta_32 - theta_31 + 2 * math.pi
    else:
        theta_30 = theta_32 - theta_31

    t_2 = t_1 + theta_30 * r

    for i in range(101, 401):
        if i >= 101 and i < 201:
            t = i - 100
            theta_2 = theta_21 - t / 100 * theta_20 * V_head
            points[i][0] = (
                x_2r + 2 * r * math.cos(theta_2),
                y_2r + 2 * r * math.sin(theta_2),
            )

        elif i >= 201 and i < 301:
            t = i - 200
            theta_3 = theta_31 + theta_30 * t / 100 * V_head
            points[i][0] = (
                x_r + r * math.cos(theta_3),
                y_r + r * math.sin(theta_3),
            )
        elif i >= 301 and i < 401:
            t = i - 300
            theta_this = theta_cal(1, theta_0)
            x_this, y_this = x_y_cal(theta_this)
            points[i][0] = (-x_this, -y_this)

        x_head[i] = points[i][0][0]
        y_head[i] = points[i][0][1]
        if i == 130:
            draw_points(x_head, y_head)

        points[i][1] = x_y_find(
            length_head,
            points[i][0][0],
            points[i][0][1],
            i,
            theta_0,
            theta_21,
            theta_20,
            x_2r,
            y_2r,
            x_r,
            y_r,
            theta_31,
            theta_30,
            theta_this,
        )

        for j in range(2, num_all):
            points[i][j] = x_y_find(
                length_head,
                points[i][j - 1][0],
                points[i][j - 1][1],
                i,
                theta_0,
                theta_21,
                theta_20,
                x_2r,
                y_2r,
                x_r,
                y_r,
                theta_31,
                theta_30,
                theta_this,
            )

        theta_this_t = theta_cal(0.000001, theta_this)
        Speed[t][0] = 1

        theta_i_t = theta_i_cal(length_head, theta_this_t)
        Speed[t][1] = (
            -(S_calculation(theta_i_t) - S_calculation(Theta[t][1])) / 0.000001
        )

        for j in range(1, num_all):
            theta_i_t = theta_i_cal(length_body, theta_i_t)
            Speed[i][j + 1] = (
                -(S_calculation(theta_i_t) - S_calculation(Theta[i][j + 1])) / 0.000001
            )

        if judge(0) == True:
            low = V_head
        else:
            high = V_head

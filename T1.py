import math
import pprint
import matplotlib.pyplot as plt

"""
初始化
"""
right_num = 1e-13

V_head = 1.0  # 龙头速度
num_all = 223 # 总数
T_total = 300  # 总时间
length_head = 2.86  # 龙头长度
length_body = 1.65  # 其余长度
a = 0.55 / (2 * math.pi)  # 螺线参数
theta_start = 32 * math.pi  # 龙头初始角度

"""
求弧长公式
"""
def S_calculation(theta):
    return 0.5 * a * (theta * math.sqrt(1 + theta * theta) + math.log(theta + math.sqrt(1 + theta * theta)))

"""
计算当前 S 值下的角度
"""
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

"""
绘制每个时间段龙头前端手把所在位置
"""
def draw_points(points):
    # 绘制极坐标图
    plt.figure(figsize=(8, 8))
    ax = plt.subplot(111, projection='polar')

    # 提取 x 和 y 坐标
    for t, (x, y) in enumerate(points):
        r = math.sqrt(x**2 + y**2)  # 计算极坐标中的 r
        theta = math.atan2(y, x)  # 计算极坐标中的 θ
        ax.plot(theta, r, marker='o')  # 绘制点
        ax.text(theta, r, str(t), fontsize=8)  # 标注点的时间 t

    # 设置图形标题
    ax.set_title("极坐标图", va='bottom')
    plt.show()

"""
计算两点之间距离
"""
def distance(theta1, theta2):
    r1 = a * theta1
    r2 = a * theta2
    return math.sqrt( (r1 * r1) + ( r2 * r2 ) - 2 * r1 * r2 * math.cos(theta1 - theta2) )

"""
计算当前i的角度
"""

def theta_i_cal( length ,theta):
    low = theta
    high = theta + 10
    while high - low > right_num:  # 精度控制
        mid = (low + high) / 2
        if distance(mid,theta) >= length:
            high = mid
        else:
            low = mid
    return high


S_start = S_calculation(theta_start)  # 初始位置

points_Head = [] # 存储每个 t 对应的 龙头(x, y) 坐标
points = [[(0,0) for _ in range(T_total+1)] for _ in range(T_total+1)] # 存储每个 t 对应的 i 的(x, y) 坐标
Speed = [[ 0 for _ in range(T_total+1)] for _ in range(T_total+1)] # 存储每个 t 对应的 i 的(x, y) 速度
Theta = [[ 0 for _ in range(T_total+1)] for _ in range(T_total+1)] # 存储每个 t 对应的 i 的角度

theta_this = 32 * math.pi

for t in range(0, T_total+1):
    if t == 0:
        theta_this = theta_start
    else :
        theta_this = theta_cal(1,theta_this)  # 当前角度
    x_this = a * theta_this * math.cos(theta_this) # 当前 x 坐标
    y_this = a * theta_this * math.sin(theta_this) # 当前 y 坐标

    points_Head.append((x_this, y_this)) # 存储点

    points[t][0] = (x_this, y_this) # 存储点

    theta_i = theta_i_cal(length_head , theta_this)
    points[t][1] = (a * theta_i * math.cos(theta_i), a * theta_i * math.sin(theta_i))
    Theta[t][1] = theta_i

    """
    求每个点的位置
    """

    for i in range(1 , num_all):
        theta_i = theta_i_cal(length_body , theta_i)
        points[t][i+1] = (a*theta_i*math.cos(theta_i), a*theta_i*math.sin(theta_i))
        Theta[t][i+1] = theta_i

    """
    求每个点的速度
    """

    theta_this_t = theta_cal(0.000001,theta_this)
    Speed[t][0] = 1

    theta_i_t = theta_i_cal(length_head , theta_this_t)
    Speed[t][1] = - ( S_calculation(theta_i_t)- S_calculation(Theta[t][1]) ) / 0.000001

    for i in range(1, num_all):
        theta_i_t = theta_i_cal(length_body, theta_i_t)
        if( theta_i_t > Theta[t][i+1] ):
            print("t:",t,"i:",i,"theta_i_t:",theta_i_t,"Theta[t][i+1]:",Theta[t][i+1])
        Speed[t][i+1] = - ( S_calculation(theta_i_t) - S_calculation(Theta[t][i+1]) ) / 0.000001


draw_points(points_Head) #绘制
for i in range(0,T_total):
    if i == 0 or i == 1 or i == 51 or i == 101 or i == 151 or i == 201:
        draw_points(points[i])


import pandas as pd

# 创建一个 Excel 文件
with pd.ExcelWriter('output.xlsx') as writer:
    # 将 Speed 导出为 Excel
    speed_df = pd.DataFrame(Speed)
    speed_df.to_excel(writer, sheet_name='Speed', index=False, header=[f'Point {i}' for i in range(T_total+1)])

    # 将 points 导出为 Excel
    points_df = pd.DataFrame(points[0])
    points_df.columns = ['X', 'Y']
    points_df.to_excel(writer, sheet_name='Points0', index=False)


    points_df = pd.DataFrame(points[60])
    points_df.columns = ['X', 'Y']
    points_df.to_excel(writer, sheet_name='Points60', index=False)

    points_df = pd.DataFrame(points[120])
    points_df.columns = ['X', 'Y']
    points_df.to_excel(writer, sheet_name='Points120', index=False)

    points_df = pd.DataFrame(points[180])
    points_df.columns = ['X', 'Y']
    points_df.to_excel(writer, sheet_name='Points180', index=False)

    points_df = pd.DataFrame(points[240])
    points_df.columns = ['X', 'Y']
    points_df.to_excel(writer, sheet_name='Points240', index=False)

    points_df = pd.DataFrame(points[300])
    points_df.columns = ['X', 'Y']
    points_df.to_excel(writer, sheet_name='Points300', index=False)

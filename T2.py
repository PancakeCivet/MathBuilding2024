import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

right_num = 1e-13

V_head = 1.0  # 龙头速度
num_all = 223  # 总数
length_head = 2.86  # 龙头长度
length_body = 1.65  # 其余长度
a = 0.55 / (2 * math.pi)  # 螺线参数
theta_start = 32 * math.pi  # 龙头初始角度
wide = 0.3

theta_this = theta_start

points = [[(0,0) for _ in range(num_all+1)] for _ in range(300000)]
Theta =  [[ 0 for _ in range(num_all+1)] for _ in range(300000)]

def S_calculation(theta):
    return 0.5 * a * (theta * math.sqrt(1 + theta * theta) + math.log(theta + math.sqrt(1 + theta * theta)))

S_start = S_calculation(theta_start)
theta_this = theta_start

cnt = 0

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

def draw_points(points,time):
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
    ax.set_title(f"t={time}", va='bottom')
    plt.show()

def distance(theta1, theta2):
    r1 = a * theta1
    r2 = a * theta2
    return math.sqrt( (r1 * r1) + ( r2 * r2 ) - 2 * r1 * r2 * math.cos(theta1 - theta2) )

def dis(x1, y1, x2, y2):
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

def theta_i_cal( length ,theta):
    low = theta
    high = theta + 15
    while high - low > right_num:  # 精度控制
        mid = (low + high) / 2
        if distance(mid,theta) >= length:
            high = mid
        else:
            low = mid
    return high

def draw_rotated_rectangles(domian,t):
    fig, ax = plt.subplots(figsize=(20, 20))

    # 遍历 domian 中的每个矩形
    for rect in domian:
        # 提取矩形的四个顶点坐标
        x1, y1 = rect[0]
        x2, y2 = rect[2]
        x3, y3 = rect[1]
        x4, y4 = rect[3]

        # 创建一个包含四个顶点的多边形 (不规则矩形)
        polygon = Polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4)], edgecolor='blue', facecolor='none')

        # 添加多边形到图形中
        ax.add_patch(polygon)

    # 设置坐标轴比例
    ax.set_aspect('equal', 'box')

    # 设置图形显示范围
    ax.autoscale()

    # 添加标题和标签
    plt.title(f"t = {t}")
    plt.xlabel("X")
    plt.ylabel("Y")

    # 显示图形
    plt.show()

def dis_point_line(point, point_1,point_2):
    return abs( (point_2[1] - point_1[1]) * point[0] - (point_2[0] - point_1[0]) * point[1] + point_2[0] * point_1[1] - point_2[1] * point_1[0] ) / math.sqrt( (point_2[1] - point_1[1])**2 + (point_2[0] - point_1[0])**2 )

def rectangles_overlap(rect1, rect2):
    for i in range(0,4):
        if( dis_point_line(rect1[i],rect2[3],rect2[0]) <= wide and dis_point_line(rect1[i],rect2[1],rect2[2]) <= wide and
                dis_point_line(rect1[i],rect2[2],rect2[0]) <= length_body and dis_point_line(rect1[i],rect2[1],rect2[3]) <= length_body ):
            return False

    return True

def domain_cal( rect1 , rect2 , length ):
    point_ahead_x , point_ahead_y = rect1
    point_back_x , point_back_y = rect2
    return [ ( point_ahead_x + ( 0.15 / length ) * ( point_back_y - point_ahead_y ) / 2 + ( 0.275 / length ) * ( point_ahead_x - point_back_x ) / 2
               , point_ahead_y + ( 0.15 / length ) * ( point_ahead_x - point_back_x ) / 2 + ( 0.275 / length ) * ( point_ahead_y - point_back_y ) / 2 )
        , ( point_back_x + ( 0.15 / length ) * ( point_ahead_y - point_back_y ) / 2 + ( 0.275 / length ) * ( point_back_x - point_ahead_x ) / 2
            , point_back_y + ( 0.15 / length ) * ( point_back_x - point_ahead_x ) / 2 + ( 0.275 / length ) * ( point_back_y - point_ahead_y ) / 2 )
        ,(  point_ahead_x - ( 0.15 / length ) * ( point_back_y - point_ahead_y ) / 2 + ( 0.275 / length ) * ( point_ahead_x - point_back_x ) / 2
            ,   point_ahead_y - ( 0.15 / length ) * ( point_ahead_x - point_back_x ) / 2 + ( 0.275 / length ) * ( point_ahead_y - point_back_y ) / 2 )
        ,(  point_back_x - ( 0.15 / length ) * ( point_ahead_y - point_back_y ) / 2 + ( 0.275 / length ) * ( point_back_x - point_ahead_x ) / 2
            , point_back_y - ( 0.15 / length ) * ( point_back_x - point_ahead_x ) / 2 + ( 0.275 / length ) * ( point_back_y - point_ahead_y ) / 2 ) ]

def time_cal(theta_this):
    S_this = S_calculation(theta_start)- S_calculation(theta_this)
    return S_this * 1.0

def judge(cnt):
    domain = [[(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)] for _ in range(num_all)]

    for i in range(0,num_all):
        if i == 0 :
            domain[i] = domain_cal( points[cnt][i] , points[cnt][i+1] , length_head / 2 )
        else :
            domain[i] = domain_cal( points[cnt][i] , points[cnt][i+1] , length_body / 2 )

    for i in range(0,num_all-3):
        for j in range(i+2,num_all):
            if rectangles_overlap(domain[j], domain[i]) == False or rectangles_overlap(domain[i], domain[j]) == False :
                draw_rotated_rectangles(domain,time_cal(Theta[cnt][0]))
                return False

    if round(distance(Theta[cnt][0],Theta[cnt][1]),2) < 2.86:
        draw_rotated_rectangles(domain,time_cal(Theta[cnt][0]))
        return False

    for i in range(0,num_all,1):
        if cnt != 0 and Theta[cnt-1][i] < Theta[cnt][i]:
            draw_rotated_rectangles(domain,time_cal(Theta[cnt][0]))
            return False

    return True

Speed = [[ 0 for _ in range(num_all+1)] for _ in range(300000)] # 存储每个 t 对应的 i 的(x, y) 速度

while True:

    print(time_cal(theta_this))

    x_this = a * theta_this * math.cos(theta_this)
    y_this = a * theta_this * math.sin(theta_this)
    points[cnt][0] = (x_this,y_this)

    theta_i = theta_i_cal(length_head,theta_this)
    points[cnt][1] = (a * theta_i * math.cos(theta_i),a * theta_i * math.sin(theta_i))

    Theta[cnt][0] = theta_this
    Theta[cnt][1] = theta_i

    for i in range(1 , num_all):
        theta_i = theta_i_cal(length_body , theta_i)
        points[cnt][i+1] = (a*theta_i*math.cos(theta_i), a*theta_i*math.sin(theta_i))
        Theta[cnt][i+1] = theta_i

    """
   求每个点的速度
   """

    theta_this_t = theta_cal(0.000001,theta_this)
    Speed[cnt][0] = 1

    theta_i_t = theta_i_cal(length_head , theta_this_t)
    Speed[cnt][1] = - ( S_calculation(theta_i_t)- S_calculation(Theta[cnt][1]) ) / 0.000001

    for i in range(1, num_all):
        theta_i_t = theta_i_cal(length_body, theta_i_t)
        Speed[cnt][i+1] = - ( S_calculation(theta_i_t) - S_calculation(Theta[cnt][i+1]) ) / 0.000001


    if judge(cnt) == 0 :
        break

    if time_cal(theta_this) > 412.473900 :
        theta_this = round( theta_this - 0.000001, 13 )

    elif time_cal(theta_this) > 412 :
        theta_this = round( theta_this - 0.0001, 13 )
    else:
        theta_this = round( theta_this - 0.1, 13 )


import pandas as pd

# 假设 Speed 和 points 是你已经定义好的数据

# 将 Speed 导出为一个独立的 Excel 文件
Speed_man = Speed[:224]

speed_df = pd.DataFrame(Speed_man)
header = [f'Point {i}' for i in range(speed_df.shape[1])]
speed_df.to_excel('speed_output.xlsx', sheet_name='Speed', index=False, header=header)

points_df = pd.DataFrame(points[cnt])
points_df.columns = ['X', 'Y']  # 为 points_df 设置列名
points_df.to_excel(f'points_output_{cnt}.xlsx', sheet_name=f'Points{cnt}', index=False)
#!/usr/bin/python 
# -*- coding: utf-8 -*-
"""
 @Author: Jane
 @FileName: 二维对流弥散模拟-显示上游加权.py
 @DateTime: 2022/7/21 21:09
 @SoftWare: PyCharm
"""
import numpy as np
import matplotlib.pyplot as plt

# 模型相关参数设置
u = 2             # 地下水水平渗流速度，m/d
w = 1.37e-4       # 地下水垂向渗流速度，m/d
n = 0.5           # 土壤孔隙度
pb = 1370         # 土壤干密度，kg/m3
aL = 1            # 纵向弥散度，m
aT = 0.3          # 横向弥散度，m
Kd = 3.63e-3      # 污染物吸附分配系数，m3/kg
k = 0.003         # 污染物反应速率常数，d-1
t = 30           # 模拟时长，d

# 计算参数
Dx = aL*u         # 纵向弥散系数，m2/d
Dz = aT*w         # 横向弥散系数，m2/d
R = 1+pb/n*Kd     # 阻滞因子/延迟因子

# 模拟区域情况
lenX = 30         # x方向长度，m
lenZ = 5          # z方向深度，m
qs = 1            # 源汇处单位体积含水层的流量，d-1
Cs = 0.25         # 源汇项的浓度，kg/m3
S = qs*Cs         # 源汇项的量，kg/m3/d

# 空间、时间离散
dx = 0.5           # x方向空间步长，m
dz = 0.1           # z方向空间步长，m
numx = int(lenX/dx)+1
numz = int(lenZ/dz)+1
x = np.linspace(0, lenX, numx)
z = np.linspace(0, lenZ, numz)
xx, zz = np.meshgrid(x, z)
dt = 0.1        # 时间步长，d
numt = int(t/dt)

# 初始条件
C = np.zeros((numx, numz))
C[int(numx/2), -int(2/dz)] = 0.25
print(C.shape)

# 初始条件下作图
# C = C.transpose()
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1.set_title('Initial Concertation')
# ax1.set_xlabel('x')
# ax1.set_ylabel('z')
# cf = ax1.contourf(xx, zz, C)
# fig.colorbar(cf)
# plt.show()

# 迭代计算
# todo 核心问题所在
for h in range(0, numt+1):
    C = C.copy()
    C[int(numx/2), -int(2/dz)] = 0.25
    for i in range(1, numx-1):
        for j in range(1, numz-1):
            C[i, j] = -u*dt/R/dx*(C[i, j]-C[i-1, j])-w*dt/R/dz*(C[i, j+1]-C[i, j]) + \
                      Dx*dt/R/dx/dx*(C[i+1, j]-2*C[i, j]+C[i-1, j])+Dz*dt/R/dz/dz*(C[i, j+1]-2*C[i, j]+C[i, j-1]) + \
                      -k*dt*C[i, j]+C[i, j]
# 边界条件
    C[:, -1] = 0  # 上边界
    C[:, 0] = 0  # 下边界
    C[-1, :] = 0  # 右边界
    C[0, -1] = 0  # 左边界

C = C.transpose()

# 面向对象的作图
fig = plt.figure()
ax = fig.add_subplot(111)
cf = ax.contourf(xx, zz, C)
ax.set_xlabel('x')
ax.set_ylabel('z')
ax.set_title('Concerntation')
fig.colorbar(cf)
plt.show()

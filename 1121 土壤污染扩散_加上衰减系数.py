# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 21:48:32 2021

@author: 
"""
import numpy as np
import pandas as pd
import matplotlib as plt
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime






def LJY_alc_con(x, z):#x,z分别表示场地的水平距离和深度(m)
    yorn=input('请确认是否要自主输入参数（是输入y，否输入n）\n')
    if yorn=='y':


        x0=float(input('请输入污染源的水平长度（以 m 计）\n'))
        z0=float(input('请输入污染源的深度（以 m 计）\n'))
        t=float(input('请输入需模拟的时间（以 d 计）\n'));
        # m=float(input('请输入释放的时间（以s计）\n'));
        dt=float(input('请输入微元时间间隔dt(d)\n'));
        c=float(input('请输入初始污染物浓度（以 mg/L 计）\n'));
        Dz=float(input('请输入垂向弥散系数（以 m2/d 计）\n'));
        Dx=float(input('请输入纵向弥散系数（以 m2/d 计）\n'));
        
        u=float(input('请输入纵向速度（以 m/d 计）\n'));
        w=float(input('请输入垂向速度（以 m/d 计）\n'));
        
        dxz=float(input('请输入网格计算步长（以 m 计）\n'));
        
        kk=float(input('请输入降解系数（以小数计，为负值）\n'));
        
    else :
        
        x0=20.0
        z0=2
        t=10
        # m=float(input('请输入释放的时间（以s计）\n'));
        dt=0.05
        c=50.0
        Dz=1.0
        Dx=1.0
        
        u=0.2
        w=0.1
        
        dxz=2
        
        kk=float(-0.3)

    X = int(x / dxz);
    Z = int(z / dxz);
    N = int(t / dt);
    
    # dist=int(dist);
    # dd=int(dist/10+1);
    x1 = np.arange(0, x, dxz)
    z1 = np.arange(0, z, dxz)
    xx, zz = np.meshgrid(x1,z1)
    
    y=np.zeros((X,Z));
    
    yk=np.zeros((X,Z));
    
    for i in range(int(x0/dxz)):
        
        for j in range(int(z0/dxz)):
            
            y[i,j]=c;
    
    yn=list();
    
    
    for k in range(N):
        
        for i in range(1,X-1):
            
            for j in range(1,Z-1):
                i=int(i);
                j=int(j);
                yk[i,j]=float((1+dt*kk)*y[i,j] + dt*Dz*(y[i,j+1]-2*y[i,j]+y[i,j-1])/dxz/dxz + dt*Dx*(y[i+1,j]-2*y[i,j]+y[i-1,j])/dxz/dxz-dt*u*(y[i+1,j]-y[i-1,j])/2/dxz-dt*w*(y[i,j+1]-y[i,j-1])/2/dxz);
        
        
        for i in range(X):
            yk[i,0]=yk[i,1];
            yk[i,Z-1]=yk[i,Z-2];
        
        for j in range(Z):
            yk[0,j]=yk[1,j];
            yk[X-1,j]=yk[X-2,j];
        
        y=yk;
        yn.append(yk);
        
    
    #%%
    
    
    
    
    yk=yk.transpose();
    
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
    fig, ax = plt.subplots(1, 1)
    
    C = plt.contour(xx,zz,yk,8, cmap='plasma')
    plt.contourf(xx,zz,yk,8, cmap='plasma')
    # 添加标记，标记处不显示轮廓线，颜色为黑色，保留两位小数
    plt.clabel(C, inline=True, colors='k', fmt='%1.2f')
    # 显示颜色条
    
    
    cb=plt.colorbar(C)
    # cb.ax.tick_params(labelsize=16)
    plt.title('等高线图颜色填充')
    plt.xlabel('x axis label')
    plt.ylabel('y axis label')
    ax.set_ylim((z-2, 0))
    plt.rcParams['savefig.dpi'] = 1000 # 图片像素
    plt.rcParams['figure.dpi'] = 1000 # 分辨率
    plt.show()
    
    
    # fig, ax = plt.subplots(1, 1)
    # C = plt.contour(xx,zz,yk,8, cmap='plasma')
    # # 添加标记，标记处不显示轮廓线，颜色为黑色，保留两位小数
    # plt.clabel(C, inline=True, colors='k', fmt='%1.2f')
    # # 显示颜色条
    # plt.colorbar(C)
    # plt.title('等高线图设置渐变色')
    # plt.xlabel('x axis label')
    # plt.ylabel('y axis label')
    # ax.set_ylim((z-2, 0))
    
    # plt.show()
    
    
    
    #%%

    
    # writer = pd.ExcelWriter('D:\\卤代烃迁移转化模拟\\' + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-") +'air_pollution_spread.xls')
    yk=pd.DataFrame(yk);
    
    yk.to_csv('H:\\08 土壤课题-个人推进情况\\03 Python模型代码\\' + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-") +'Soil_pollution_spread.csv');
    
    # writer.save()


LJY_alc_con(500, 20)

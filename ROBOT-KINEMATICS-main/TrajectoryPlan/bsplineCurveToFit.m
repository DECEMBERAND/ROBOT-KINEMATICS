

% 此文件为b样条曲线插值算法
%k = 0 Bi,0(t) = {1, t=[ti,ti+1]
%                 0, otherwise
%k > 0 Bi,k(t) = (t-ti)/(ti+k-ti) * Bi,k-1(t) + (ti+k+1-t)/(ti+k+1-ti+1) * Bi+1,k-1(t)
%
% 均匀 B 样条：节点均匀分布
% 准均匀 B 样条：在开始和结束处的节点可重复，中间节点均匀分布
% 非均匀 B 样条：节点非均匀分布
% 分段贝塞尔曲线
% PLUS：B样条无法描述圆锥曲线，为解决此问题，产生了非均匀有理B样条（non-uniform rational b-spline, NURBS）

clear;
clc;
 
% 定义控制点
x = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115];
y = [0, 0, 0, -1, -3, -4, -6, -8, -11, -12, -11, -11, -11, -11, -12, -13, -15, -18, -19, -21, -23, -22, -24, -24];
z = [-2, -1, -0, -1, -2, -3, -2.3, -1, 2, 1, 2.1, 1.5, 2, -1, -2.4, -2.6, -4, -3, -2, -1, 0, 1, 2, 3];
% 过末端点需要重复最后一个数据
x = [x, x(end)];
y = [y, y(end)];
z = [z, z(end)];
% P3d = randi(100,5,3);
% x = P3d(:,1)';
% y = P3d(:,2)';
% z = P3d(:,3)';

n=length(x)-1;% B样条曲线有n+1个控制点
k = 4;% 三次B样条，因此阶数=次数+1

% 用于储存生成的点
X=[];
Y=[];
Z=[];
 
Nik_u=zeros(1,n+1);% 用于储存基函数
nodevector=U_quasi_uniform(n,k);% 生成节点向量

for u=0:0.001:1-0.01% 循环次数决定组成B样条曲线点的总数量
    for j=1:n+1% 生成Bi,3这一系列基函数
        Nik_u(j)=BaseFunction(j,k,u,nodevector);
    end
    % 由于B样条曲线具有局部支撑性，因此只有u周围一定区间会参与到计算
    X=[X,x*Nik_u'];
    Y=[Y,y*Nik_u'];
    Z=[Z,z*Nik_u'];
end
figure
plot3(X,Y,Z,'r')
hold on
plot3(x,y,z, 'g.-', 'MarkerSize',40)
plot3(X(1), Y(1), Z(1), 'b*', MarkerSize= 40)
plot3(X(end), Y(end), Z(end), 'b*', MarkerSize= 40)
title(num2str(k-1),'次BSpline-3d')
axis manual
hold off
% 生成准均匀B样条曲线
function nodevector=U_quasi_uniform(n, k)
nodevector=zeros(1,n+k+1);% 节点数=控制点的个数+阶数
piecewise = n - k + 2;% B样条曲线的段数=控制点个数-次数
if piecewise==1% 只有一段曲线时，直接末尾重复度k
    nodevector(n+1+1:n+k+1)=1;
else
    for i=1:n-k+1% 中间段内节点均匀分布：两端共2k个节点，中间还剩(n+k+1-2k=n-k+1）个节点
        nodevector(k+i)=nodevector(k+i-1)+1/piecewise;
    end
    nodevector(n+1+1:n+k+1)=1;% 末尾重复度k
end
end
%定义第i个k阶B样条基函数
function Nik_u=BaseFunction(i, k, u, NodeVector)
if k==1% 定义Bi,0这一系列基函数
    if u>=NodeVector(i)&&u<NodeVector(i+1)
        Nik_u=1;
    else
        Nik_u=0;
    end
else
    % 公式中的两个分母
    denominator_1 = NodeVector(i + k - 1) - NodeVector(i);
    denominator_2 = NodeVector(i + k) - NodeVector( i + 1);
    % 如果遇到分母为0的情况，定义为1以便继续计算
    if denominator_1==0
        denominator_1=1;
    end
    if denominator_2==0
        denominator_2=1;
    end
    % 不断递归
    Nik_u=(u - NodeVector(i)) / denominator_1 * BaseFunction(i, k - 1, u, NodeVector) ...
    + (NodeVector(i + k) - u) / denominator_2 * BaseFunction(i + 1, k - 1, u, NodeVector);
end
end
% https://blog.csdn.net/qq_61517948/article/details/128631732?spm=1001.2101.3001.6650.3&utm_medium=distribute.pc_relevant.none-task-blog
% -2%7Edefault%7ECTRLIST%7ERate-3-128631732-blog-108495272.235%5Ev39%5Epc_relevant_anti_t3&depth_1-utm_source=distribute.pc_relevant.
% none-task-blog-2%7Edefault%7ECTRLIST%7ERate-3-128631732-blog-108495272.235%5Ev39%5Epc_relevant_anti_t3&utm_relevant_index=6

%一版为python代码
% import numpy as np
% import matplotlib.pyplot as plt
% from mpl_toolkits.mplot3d import Axes3D
% 
% # 计算在某一特定t下的 B_{i,k}
% def getBt(controlPoints, knots, t):
%     # calculate m,n,k
%     m = knots.shape[0]-1
%     n = controlPoints.shape[0]-1
%     k = m - n - 1
%     # initialize B by zeros 
%     B = np.zeros((k+1, m))
% 
%     # get t region
%     tStart = 0
%     for x in range(m+1):
%         if t==1:
%             tStart = m-1
%         if knots[x] > t:
%             tStart = x-1
%             break
%      
%     # calculate B(t)
%     for _k in range(k+1):
%         if _k == 0:
%             B[_k, tStart] = 1
%         else:
%             for i in range(m-_k):
%                 if knots[i+_k]-knots[i]== 0:
%                     w1 = 0
%                 else:
%                     w1 = (t-knots[i])/(knots[i+_k]-knots[i]) 
%                 if knots[i+_k+1]-knots[i+1] == 0:
%                     w2 = 0
%                 else:
%                     w2 = (knots[i+_k+1]-t)/(knots[i+_k+1]-knots[i+1])
%                 B[_k,i] = w1*B[_k-1, i] + w2*B[_k-1, i+1]
%     return B
% 
% # 绘制 B_{i,k}(t)函数
% def plotBt(Bt,num, i,k):
%     print(k,i)
%     Bt = np.array(Bt)
%     tt = np.linspace(0,1,num)
%     yy = [Bt[t,k,i] for t in range(num)]
%     plt.plot(tt, yy)
% 
% # 根据最后一列（最高阶次）的 B(t)，即权重，乘以控制点坐标，从而求出曲线上点坐标
% def getPt(Bt, controlPoints):
%     Bt = np.array(Bt)
%     ptArray = Bt.reshape(-1,1) * controlPoints
%     pt = ptArray.sum(axis = 0)
%     return pt
% 
% # 绘制出生成的样条曲线: useReg 表示是否使用曲线有效定义域[t_k, t_{m-k}]
% def main1(useReg = False):
%     controlPoints = np.array([[50,50,50], [100,300,40], [300,100,30], [380,200,20], [400,600,10]])
%     knots = np.array([0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9,1])
%     m = knots.shape[0]-1
%     n = controlPoints.shape[0]-1
%     k = m - n - 1
%     # print('knots.shape: ',knots.shape)
%     print('n:',n)
%     print('m:',m)
%     print('k:',k)
%     # print('t: ', np.linspace(0,1,100))
%     fig = plt.figure()
%     ax = fig.add_subplot(111, projection='3d')
%     for t in np.linspace(0,1,100):
%         if useReg and not(t >= knots[k] and t<= knots[n+1]):
%             continue
%         Bt = getBt(controlPoints, knots, t)
%         Pt = getPt(Bt[k, :n+1], controlPoints)
%         ax.scatter(Pt[0],Pt[1], Pt[2] ,color='b')
%         # if t==1:
%         #     print('Bt: ',Bt)
%     ax.scatter(controlPoints[:,0], controlPoints[:,1], controlPoints[:,2],color = 'r')
%     plt.show()
% 
% # 绘制 B_{i,k} 变化图:如果不给定{i,k}则显示所有B{i,k}(t)图像
% def main2(i=-1,k=-1):
%     controlPoints = np.array([[50,50], [100,300], [300,100], [380,200], [400,600]])
%     knots = np.array([0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9,1])
%     m = knots.shape[0]-1
%     n = controlPoints.shape[0]-1
%     k = m - n - 1
%     print('n:',n)
%     print('m:',m)
%     print('k:',k)
%     B = []
%     num = 100 # 离散点数目
%     for t in np.linspace(0,1,num):
%         Bt = getBt(controlPoints, knots, t)
%         B.append(list(Bt))
% 
%     figure1 = plt.figure('B_{i,k}')
%     if i==-1:
%         fig = []
%         for i in range(n+1):
%             for k in range(k+1):
%                 plotBt(B,num, i,k)
%                 fig.append('B_{%d,%d}'%(i,k))
%     else:
%         plotBt(B,num, i,k)
%         fig.append('B_{%d,%d}'%(i,k))
%     plt.legend(fig)
%     plt.show()   
%     
% if __name__ == '__main__':
%     main1()
%     # main2()




       



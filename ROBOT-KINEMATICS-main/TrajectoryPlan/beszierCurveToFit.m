% 此文件为贝塞尔曲线拟合，输出n个控制点，输出一条光滑连续的轨迹，
% 此轨迹只过起点和终点，不经过其他的点，改变其中一个控制点会影响整个轨迹

% 公式：  n       (n     i       n-i
%       Bi (t) =   i) * t * (1-t)         t=[0,1]
% n = 贝塞尔曲线的次数 = 控制点数 - 1
% i = [0,n]
close all
clc;
clear
% 定义控制点
x = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115];
y = [0, 0, 0, -1, -3, -4, -6, -8, -11, -12, -11, -11, -11, -11, -12, -13, -15, -18, -19, -21, -23, -22, -24, -24];
z = [-2, -1, -0, -1, -2, -3, -2.3, -1, 2, 1, 2.1, 1.5, 2, -1, -2.4, -2.6, -4, -3, -2, -1, 0, 1, 2, 3];
P2d = [x;y]';
P3d = [x;y;z]';
% controlPoibtNum = 5;
% P2d = randi(100, controlPoibtNum, 2); % 5*2
% P3d = randi(100, controlPoibtNum, 3); % 5*3

position2d = computeCubicBezierCurve(P2d, 1000);
position3d = computeCubicBezierCurve(P3d,1000);
figure

plot(P2d(:,1),P2d(:,2),'g.-','MarkerSize',40)
hold on
plot(position2d(:,1),position2d(:,2))
hold on
plot(position2d(1,1),position2d(1,2),'*','MarkerFaceColor','red')
hold off
axis manual
title('Bezier-2d')
% pic_num=1;
% for k = 2:1000
%     p.XData = position2d(k,1);
%     p.YData = position2d(k,2);
% 
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     if pic_num==1
%     imwrite(I,map,'test.gif','gif','Loopcount',inf,'DelayTime',0.2);
%     elseif mod(pic_num,3)==1
%     imwrite(I,map,'test.gif','gif','WriteMode','append','DelayTime',0.2);
%     end
%     pic_num = pic_num + 1;
%  
%     drawnow limitrate
% end

figure
plot3(P3d(:,1), P3d(:,2), P3d(:,3),'g.-','MarkerSize',40)
hold on
plot3(position3d(:,1), position3d(:,2), position3d(:,3))
hold on
plot3(position3d(1,1), position3d(1,2), position3d(1,3))
hold off
axis manual
title('Bezeir-3d')
% https://blog.csdn.net/m0_68926749/article/details/134837709
function points = computeCubicBezierCurve(controlPoints, numPoints)
    numControlPoints = size(controlPoints, 1);
    t = linspace(0, 1, numPoints); % 在0到1之间生成numPoints个等间距的参数值
    [r,c] = size(controlPoints);
    points = zeros(numPoints, c); % 存储曲线上的点坐标
    
    for i = 1:numPoints
        B = zeros(1, c);
        for j = 1:numControlPoints
            B = B + nchoosek(numControlPoints-1, j-1) * t(i)^(j-1) * (1 - t(i))^(numControlPoints-j) * controlPoints(j, :);
        end
        points(i, :) = B; % 存储点坐标
    end
end

%-----------------------------------------------------------------------------------
%N = 3: P = (1-t)^2P0 + 2(1-t)tP1 + t^2*P2
%N = 4: P = (1-t)^3P0 + 3(1-t)^2tP1 + 3(1-t)t^2P2 + t^3*P3
%N = 5: P = (1-t)^4P0 + 4(1-t)^3tP1 + 6(1-t)2*t2P2 + 4(1-t)t^3P3 + t^4*P4
%可将贝塞尔曲线一般参数公式中的表达式用如下方式表示：
%设有常数 a,b 和 c，则该表达式可统一表示为如下形式：
%a * (1 - t)^b * t^c * Pn;
%根据上面的分析就可以总结出 a,b,c 对应的取值规则：
%b: (N - 1) 递减到 0 (b 为 1-t 的幂)
%c: 0 递增到 (N - 1) (c 为 t 的幂)
% a: 在 N 分别为 1,2,3,4,5 时将其值用如下形式表示： 
% N=1:---------1
% N=2:--------1 1
% N=3:------1 2 1
% N=4:-----1 3 3 1
% N=5:---1 4 6 4 1
% a 值的改变规则为： 杨辉三角
function ta = YHTriangle(control_points)
    N=length(control_points);
    ta=zeros(N,N);%%对数组进行初始化
    %%杨辉三角左右两边的值赋1
     
    % 杨辉三角的数的规律
    % 1
    % 1 1
    % 1 2 1
    % 1 3 3 1
    % 1 4 6 4 1
    for i=1:N
        ta(i,1)=1;
        ta(i,i)=1;
    end
    %%从第二个数开始，也就是从第三行开始，等于前列的左边加上正上方的一个
    for row=2:N
        for col=2:row
            ta(row,col)=ta(row-1,col-1)+ta(row-1,col);
        end
    end
end

function points = computeCubicBezierCurve1(control_points, numPoints) %德卡斯特里奥算法 De Casteljau
    ta = YHTriangle(control_points);
    M = numPoints;N=length(control_points);
    points = zeros(numPoints, 2);
    for i=1:M
        t=i/M;%%确定每一个点的比例
        for k=0:N-1
            c=k;%分别确定a,b,c三个系数
            b=N-c-1;%分别确定a,b,c三个系数
            a=ta(N,k+1);%分别确定a,b,c三个系数
                 
            points(i,1)=points(i,1)+a*(1-t)^b*t^c*control_points(k+1,1);%确定点的x坐标
           
            points(i,2)=points(i,2)+a*(1-t)^b*t^c*control_points(k+1,2);%确定点的y坐标

            %points(i,3)=points(i,3)+a*(1-t)^b*t^c*control_points(k+1,3);%确定点的z坐标
       end
      
    end
end














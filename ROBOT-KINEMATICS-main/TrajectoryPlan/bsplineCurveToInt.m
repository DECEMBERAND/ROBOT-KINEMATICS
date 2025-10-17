% 给定一组数据点，即型值点，反求其控制顶点，并构造通过这些型值点的2D或3D nurbs曲线，首末型值点与控制顶点重合，称为nurbs曲线插值
% 具有n个型值点的k次nurbs插值曲线，将有n+2个未知控制顶点，首末节点取重复度r=k+1,从而具有n+2+k+1个节点
% 三次Nurbs曲线的首末节点取重复度r=k+1=4，u1=u2=u3=u4=0,Un+3=Un+4=Un+5=Un+6=1
% https://blog.csdn.net/shs1992shs/article/details/107748375?spm=1001.2014.3001.5506
clear
pt=load('data1.txt'); 
[row,column]=size(pt);  %每行代表一个型值点，采用列向量表示x,y,z坐标，便于线性方程组求解
n=length(pt);  %型值点数量
k=3;  %样条阶数
U=zeros(1,n+k+3);  %节点矢量
%第一步：计算节点矢量********************
if(column == 2)   % 2D curve
    x=pt(:,1);
    y=pt(:,2);
else              % 3D curve
    x=pt(:,1);
    y=pt(:,2);
    z=pt(:,3);
end

%n个数据点，利用规范积累弦长法进行型值点参数化处理
temp=zeros(1,n-1);
for i=1:n-1  %n个数据点，n-1段弦长
   if(column == 2)    % 2D curve
   temp(i)=sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
   else               % 3D curve
   temp(i)=sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2+(z(i+1)-z(i))^2);
   end
end
sumtemp=sum(temp);%顺序相邻两型值点之间距离的和,即总弦长
for i=1:k+1      %前k+1个节点为0
    U(i)=0;
end
for i=n+k:n+k+3  %后k+1个节点为1
    U(i)=1;
end
 %n个数据点，内节点为n-2个，U(k+1)作为初始值
for i=k+1:n+k-2
    U(i+1)=U(i)+temp(i-k)/sumtemp;  %n-1段弦长，U（k+1）到U（n+k）共n个节点
end

%第二步：反算n+2个控制点，采用切矢边界*****************************
%控制顶点的首末端点和给定型值点的首末端点重合
if(column == 2)
    dpt1=[0 1];%给定首数据点切矢
    dptn=[-1 0];%给定末数据点切矢
else
    dpt1=[0 0 1];%给定首数据点切矢
    dptn=[-1 0 0];%给定末数据点切矢
end
dU=zeros(1,n+k+2);  %节点增量，△U=Ui+1-Ui
for i=k+1:n+k-1
    dU(i)=U(i+1)-U(i);
end

%求解线性方程组获得控制顶点向量，A*D=E,A为系数矩阵，元素为B样条基函数的值；D是控制顶点列向量；E是列向量
A=zeros(n);
if(column == 2)
    E=zeros(n,2);  %2D curve
else
    E=zeros(n,3);  %3D curve
end
A(1,1)=1;  %切矢条件a1=1,b1=c1=0
A(n,n)=1;  %切矢条件an=bn=0,cn=1
E(1,:)=pt(1,:)+(dU(4)/3)*dpt1;    %首端点条件
E(n,:)=pt(n,:)-(dU(n+2)/3)*dptn;  %末端点条件
%计算系数矩阵A的元素a,b,c以及列向量E的元素e的值
for i=2:n-1  
    A(i,i-1)=dU(i+3)^2/(dU(i+1)+dU(i+2)+dU(i+3));                  %a的值
    A(i,i)=dU(i+3)*(dU(i+1)+dU(i+2))/(dU(i+1)+dU(i+2)+dU(i+3))+...  %b的值
        dU(i+2)*(dU(i+3)+dU(i+4))/(dU(i+2)+dU(i+3)+dU(i+4));
    A(i,i+1)=dU(i+2)^2/(dU(i+2)+dU(i+3)+dU(i+4));                  %c的值
    E(i,:)=(dU(i+2)+dU(i+3))*pt(i,:);                               %e的值
end
D=A\E;  %解方程组，获得去除首末端点的控制顶点向量
D=[pt(1,:);D;pt(n,:)]; %控制顶点向量，控制顶点比数据点多两个，加上首末端点
[s,t]=size(D);  %s：控制点数量

%第三步：求出Nurbs曲线*****************************
dt=0.001;  %Nurbs曲线插值密度，dt越小，越光滑
P=[];      %Nurbs曲线
syms dx;   %定义节点间的递增变量为变量dx
for i=k+1:s
    u=U;
    d=sym(D);
    %每次迭代后将d恢复初值
    for m=1:k
        for j=i-k:i-m
            alpha(j)=(dx-u(j+m))/(u(j+k+1)-u(j+m));
            d(j,:)=(1-alpha(j))*d(j,:)+alpha(j)*d(j+1,:);
        end
    end
    %带入变量dx求出各节点区间[ui,ui+1]内m=3的控制点d
    M=subs(d(i-k,:),dx,(u(i):dt:(u(i+1)-dt))');    %节点区间内的插值点替换dx
    P=[P;double(M)];
end
M=subs(d(s-k,:),dx,1);
P=[P;double(M)];%补上最后一个点

%第四步、绘制给定型值点、控制点及其Nurbs曲线*****************************
figure
if(column == 2)  %2D curve
    plot(pt(:,1),pt(:,2),'*r');  %绘制型值点
    hold on;
    plot(D(:,1),D(:,2),'b-o');   %绘制控制点
    hold on;
    plot(P(:,1),P(:,2),'r');     %绘制nurbs曲线
    hold on;  
else             %3D curve
    plot3(pt(:,1),pt(:,2),pt(:,3),'*r'); %绘制型值点
    hold on;
    plot3(D(:,1),D(:,2),D(:,3),'b-o');   %绘制控制点

    hold on;
    plot3(P(:,1),P(:,2),P(:,3),'r');     %绘制nurbs曲线
    hold on;
end
axis equal;
grid on;
title('切失边界条件')

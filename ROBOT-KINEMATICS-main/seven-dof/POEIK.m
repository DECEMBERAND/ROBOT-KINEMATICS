% 此构型冗余机械臂为拟人臂，以KUKA iiwa7 为例做正运动学分析
% https://zhuandan.zhihu.com/p/345258933
dbs = 340; dse = 400; dew = 400; dwf = 126;
wlist = [0 0 1;
         0 1 0;
         0 0 1;
         0 -1 0;
         0 0 1;
         0 1 0;
         0 0 1];  % 旋转轴在基坐标系下的表示
qlist = [0 0 0;
         0 0 dbs;
         0 0 0;
         0 0 dbs+dse;
         0 0 0;
         0 0 dbs+dse+dew;
         0 0 0];
vlist = - cross(wlist, qlist);
Slist = [wlist vlist]';
M = [1 0 0 0;
     0 1 0 0;
     0 0 1 dbs+dse+dew+dwf;
     0 0 0 1];

thetalist0 = [10 20 30 40 50 60 70]'*pi/180;
% 正运动学
T = FKinSpace(M, Slist, thetalist0);
% 逆运动学
[thetalist, success] = IKinSpace(Slist, M, T, thetalist0, 0.01, 0.001);
clear; 
clc;
%-------------ABB IRB2600 20(12)\1.65 POE----------------
a1 =150;
a2 = 700;
a3 = 115;
d1 = 445;
d4 = 795;
d6 = 85;
M = [0 0 1 a1+d4+d6;
     0 1 0 0;
    -1 0 0 d1+a2+a3;
     0 0 0 1];
wlist =[ 0 0 1;
         0 1 0;
         0 1 0;
         1 0 0;
         0 1 0;
         1 0 0];       
qlist = [0 0 0;
        a1 0 d1;
        a1 0 d1+a2;
        0 0 d1+a2+a3;
        a1+d4 0 d1+a2+a3;
        0 0 d1+a2+a3];
vlist = zeros (size(qlist));
for i = 1:6
    vlist(i,:)=-cross(wlist(i,:),qlist(i,:));
end
Slist = [wlist';vlist'];
thetalist1 = [0 0 0 0 0 0];
thetalist2 = deg2rad([10 20 20 20 20 20]');
thetalist0 = [0 0 0 0 0 0]; %验证逆解需要给定一个初值
eomg = 0.001;
ev = 0.0001;
T1 = FKinSpace(M, Slist, thetalist1) % 正解1
[thetalist, success] = IKinSpace(Slist, M, T1, thetalist0, eomg, ev);%逆解1
rad2deg(thetalist)
q1 = rotm2quat(T1(1:3,1:3));%旋转矩阵对应的四元数，因为ABB是四元数显示
T2 = FKinSpace(M, Slist, thetalist2)
thetalist0 = thetalist2 - 0.5;
[thetalist, success] = IKinSpace(Slist, M, T2, thetalist0, eomg, ev);
rad2deg(thetalist)
q2 = rotm2quat (T2(1:3,1:3));
 % 用逆解结果重新计算正解，看是否匹配原T1
T1_check = FKinSpace(M, Slist, deg2rad([10, 20, 20, -1420, 20, 1460]));
disp('误差:'); disp(T1 - T1_check);
% -------------------ABB IRB2600 MDH-------------------
d1 = 445;
a1 = 150;
a2 = -700;
a3 = -115;
d4 = 795;
d6 = 85;
%   alpha (i-1)  ai    di   offset
DH = [0          0     d1    0;
     -pi/2       a1     0   pi/2;
      0          a2     0    0;
     pi/2        a3    d4    0;
    -pi/2        0     0     0;
     pi/2        0     d6    0];
thetalist1 = [0 0 0 0 0 0];
thetalist2 = deg2rad ([10 20 20 20 20 20]');
T1 = FK(thetalist1, DH)
q1 = rotm2quat (T1 (1:3,1:3));
T2 = FK(thetalist2, DH)% 解析解算正解 T06 =  T61*T12*T23*T34*T45*T56；
q2 = rotm2quat (T2(1:3,1:3));
T1_check = FKinSpace(M, Slist, deg2rad([10, 20, 20, -1420, 20, 1460]));
disp('误差:'); disp(T1 - T1_check);

function T = FK(thetalist,DH) % 6*1
        dof = max(size(thetalist));
        thetalist = thetalist + DH(:,4);
        T = eye(4);
        for i = 1:dof
            T_tmp = [cos(thetalist(i)) -sin(thetalist(i)) 0 DH(i,2);
                     sin(thetalist(i))*cos(DH(i,1)) cos(thetalist(i))*cos(DH(i,1)) -sin(DH(i,1)) -DH(i,3)*sin(DH(i,1));
                     sin(thetalist(i))*sin(DH(i,1)) cos(thetalist(i))*sin(DH(i,1)) cos(DH(i,1)) DH(i,3)*cos(DH(i,1));
                     0 0 0 1];
            T = T*T_tmp;
        end

end
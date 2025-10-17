% 图片中给出的是MDH的建系，SDH要把图中的X1和X2互换，X3，X5同理
% SDH参数：
% i     a      adpha       d            theta
% 1     0       -90        dbs =340        0  
% 2     0       90         0               0  
% 3     0       90        dse = 400       0   
% 4     0       -90         0              0  
% 5     0       -90        dew = 400       0  
% 6     0       90         0               0  
% 7     0       0          dwf = 126       0  
% MDH参数：
% i     a      adpha       d            theta
% 1     0       0          dbs =340        0   x0->x1
% 2     0       -90         0               0  x1->x2
% 3     0       90        dse = 400        0   x2->x3
% 4     0       90         0              0    x3->x4
% 5     0       -90        dew = 400       0   x4->x5
% 6     0       -90         0               0   x5->x6
% 7     0       90        dwf = 126       0   x6->x7

%% DH法建立模型,关节转角，关节距离，连杆长度，连杆转角，关节类型（0转动，1移动），'standard'：建立标准型D-H参数
theta1 = 0.1745;   D1 = 340;   A1 = 0;     alpha1 = -pi/2;   offset1 = 0;
theta2 = 0.3491;   D2 = 0;     A2 = 0;     alpha2 = pi/2;    offset2 = 0;
theta3 = 0.5236;   D3 = 400;   A3 = 0;     alpha3 = pi/2;    offset3 = 0;
theta4 = 0.6981;   D4 = 0;     A4 = 0;     alpha4 = -pi/2;   offset4 = 0;
theta5 = 0.8727;   D5 = 400;   A5 = 0;     alpha5 = -pi/2;   offset5 = 0;
theta6 = 1.0472;   D6 = 0;     A6 = 0;     alpha6 = pi/2;    offset6 = 0;
theta7 = 1.2217;   D7 = 126;   A7 = 0;     alpha7 = 0;       offset7 = 0;
L1 = Link('d', D1, 'a', A1, 'alpha', alpha1, 'standard');
L2 = Link('d', D2, 'a', A2, 'alpha', alpha2, 'standard');
L3 = Link('d', D3, 'a', A3, 'alpha', alpha3, 'standard');
L4 = Link('d', D4, 'a', A4, 'alpha', alpha4, 'standard');
L5 = Link('d', D5, 'a', A5, 'alpha', alpha5, 'standard');
L6 = Link('d', D6, 'a', A6, 'alpha', alpha6, 'standard');
L7 = Link('d', D7, 'a', A7, 'alpha', alpha7, 'standard');

% theta1 = 0.1745;   D1 = 340;   A1 = 0;     alpha1 = 0;   offset1 = 0;
% theta2 = 0.3491;   D2 = 0;     A2 = 0;     alpha2 = -pi/2;    offset2 = 0;
% theta3 = 0.5236;   D3 = 400;   A3 = 0;     alpha3 = pi/2;    offset3 = 0;
% theta4 = 0.6981;   D4 = 0;     A4 = 0;     alpha4 = pi/2;   offset4 = 0;
% theta5 = 0.8727;   D5 = 400;   A5 = 0;     alpha5 = -pi/2;   offset5 = 0;
% theta6 = 1.0472;   D6 = 0;     A6 = 0;     alpha6 = -pi/2;    offset6 = 0;
% theta7 = 1.2217;   D7 = 126;   A7 = 0;     alpha7 = pi/2;       offset7 = 0;
% L1 = Link('d', D1, 'a', A1, 'alpha', alpha1, 'modified');
% L2 = Link('d', D2, 'a', A2, 'alpha', alpha2, 'modified');
% L3 = Link('d', D3, 'a', A3, 'alpha', alpha3, 'modified');
% L4 = Link('d', D4, 'a', A4, 'alpha', alpha4, 'modified');
% L5 = Link('d', D5, 'a', A5, 'alpha', alpha5, 'modified');
% L6 = Link('d', D6, 'a', A6, 'alpha', alpha6, 'modified');
% L7 = Link('d', D7, 'a', A7, 'alpha', alpha7, 'modified');

%加入teach指令，则可调整各个关节角度
robot = SerialLink([L1, L2, L3, L4, L5, L6, L7],'name','seven');
figure(1)
%robot.plot(theta);
robot.teach
title('七轴机械臂模型可调节2');

%% urdf验证
% https://www.mathworks.com/help/robotics/ref/rigidbodytree.gettransform.html
% 姿态和建模不一致的原因是杆长有些差别
robot = importrobot('iiwa7.urdf');
q = randomConfiguration(robot); %任意姿态
q(1).JointPosition = 0.1745;
q(2).JointPosition = 0.3491;
q(3).JointPosition = 0.5236;
q(4).JointPosition = 0.6981;
q(5).JointPosition = 0.8727;
q(6).JointPosition = 1.0472;
q(7).JointPosition = 1.2217;
figure
show(robot,q);
transform = getTransform(robot,q,"iiwa_link_0","iiwa_link_7");

%% 公式计算T
dhflag = 1; % 0--SDH, 1--MDH
theta = [theta1 theta2 theta3 theta4 theta5 theta6 theta7];
d = [D1 D2 D3 D4 D5 D6 D7];
a = [A1 A2 A3 A4 A5 A6 A7];
alpha = [alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7];
dh = [a', alpha', theta', d'];
[T,T10,T20,T30,T40,T50,T60,T70] = myfunFKTransMatrix(dh,dhflag);

function [T,T10,T20,T30,T40,T50,T60,T70]=myfunFKTransMatrix(dh,DHflg)
for k=1:7
    for i=1:k
        if DHflg==0
             T(:,:,k)=myfunSDHMatrix( dh(i,1),dh(i,2),dh(i,3),dh(i,4));  % SDH
        elseif DHflg==1
             T(:,:,k)=myfunMDHMatrix( dh(i,1),dh(i,2),dh(i,3),dh(i,4));  % MDH
        elseif (DHflg~=1 && DHflg~=1)
            disp('Error input');
            break;
        end
    end
end

%disp('display each transform matrix Tn:');
% transform matrix
T10=(T(:,:,1));T21=(T(:,:,2));T32=(T(:,:,3));T43=(T(:,:,4));T54=(T(:,:,5));T65=(T(:,:,6));T76=(T(:,:,7));
T20=T10*T21;T30=T20*T32;T40=T30*T43;T50=T40*T54;T60=T50*T65;T70=T60*T76;
end

function [T]=myfunMDHMatrix(a,alpha,theta,d)
T=[ 
                   cos(theta),          -sin(theta),          0,           a;
        sin(theta)*cos(alpha),cos(theta)*cos(alpha),-sin(alpha),-sin(alpha)*d;
        sin(theta)*sin(alpha),cos(theta)*sin(alpha), cos(alpha), cos(alpha)*d;
                            0,                    0,          0,            1
       ];
end

function [T]=myfunSDHMatrix(a,alpha,theta,d)
T=[cos(theta),-sin(theta)*cos(alpha), sin(theta)*sin(alpha),a*cos(theta);
    sin(theta), cos(theta)*cos(alpha),-cos(theta)*sin(alpha),a*sin(theta);
         0,        sin(alpha),        cos(alpha),       d;
         0,                 0,                 0,       1
   ];
end
%边界条件为抛物线边界
%参考：https://blog.csdn.net/qq_32059343/article/details/86408359?spm=1001.2014.3001.5506

clc;clear;
% originPoints = [1,1,0;2,2,0;3,5,0;4,2,0;4.5,1.3,0;5,1,0;7,0.5,0;8,-1.5,0];
originPoints = load('wayPoints1.txt');
[ctrlPoints, curvePoints] = generateCurvePoints(originPoints);
figure
plot3(originPoints(:,1), originPoints(:,2), originPoints(:,3),'r*','MarkerSize',10)
hold on
plot3(ctrlPoints(:,1), ctrlPoints(:,2),ctrlPoints(:,3),'b-o')
plot3(curvePoints(:,1), curvePoints(:,2), curvePoints(:,3), 'r')
axis equal;
grid on;
title('抛物线边界条件')

function result = Nu(i, u, k, U)
    if k == 0
        if u >= U(i) && u < U(i + 1)
            result = 1;
        else
            result = 0;
        end
    else
        a1 = (u - U(i));
        c1 = (U(i + k) - U(i));
        b1 = (U(i + k + 1) - u);
        d1 = (U(i + k + 1) - U(i + 1));
        if d1 < 0.0001
            b1 = 0;
            d1 = 1;
        end
        if c1 < 0.0001
            a1 = 0;
            c1 = 1;
        end
        nu = (a1 / c1) * Nu(i, u, k - 1, U) + (b1 / d1) * Nu(i + 1, u, k - 1, U);
        result = nu;
    end
end
function [ctrlPoints, curvePoints] = generateCurvePoints(originPoints)
    % 参数初始化
    m = length(originPoints) - 1;
    k = 3;  % 代码只支持k=4次
    [row,col] = size(originPoints);
    % 积累弦长
    L = 0;
    l = zeros(1, m);
    for i = 1:m
        l(i) = sqrt((originPoints(i + 1, 1) - originPoints(i, 1))^2 + (originPoints(i + 1, 2) - originPoints(i, 2))^2 + (originPoints(i + 1, 3) - originPoints(i, 3))^2);
        L = L + l(i)
    end

    % 节点矢量
    U = zeros(1, m + 1 + 2 * k);
    for i = 1:k + 1
        U(i) = 0;
        U(i + m + k) = 1;
    end
    for i = k + 1:m + k
        L2 = 0;
        for j = 1:i - k
            L2 = L2 + l(j);
        end
        U(i+1) = L2 / L;
    end
    disp('U=');
    disp(U);

    % Delta，公式中的三角形的东西
    delta = zeros(1, m + 2 * k);
    for i = 1:m + 2 * k
        delta(i) = U(i + 1) - U(i);
    end
    disp('Delta=');
    disp(delta);

    % 系数矩阵的参数
    a = zeros(1, m + 1);
    b = zeros(1, m + 1);
    c = zeros(1, m + 1);
    e = zeros(col, m + 1);
%     f = zeros(1, m + 1);

    % 抛物线条件
    a(1) = 1 - delta(4) * delta(5) / ((delta(4) + delta(5))^2);
    b(1) = delta(4) / (delta(4) + delta(5)) * (delta(5) / (delta(4) + delta(5)) - delta(4) / (delta(4) + delta(5) + delta(6)));
    c(1) = delta(4) * delta(4) / ((delta(4) + delta(5)) * (delta(4) + delta(5) + delta(6)));
    e(:,1) = 1 / 3 * (originPoints(1, :) + 2 * originPoints(2, :));
%     f(1) = 1 / 3 * (originPoints(1, 2) + 2 * originPoints(2, 2));
    a(m + 1) = -delta(m + 3) * delta(m + 3) / ((delta(m + 2) + delta(m + 3)) * (delta(m + 2) + delta(m + 2) + delta(m + 3)));
    b(m + 1) = delta(m + 3) / (delta(m + 2) + delta(m + 3)) * (delta(m + 3) / (delta(m + 2) + delta(m + 2) + delta(m + 3)) - delta(m + 2) / (delta(m + 2) + delta(m + 3)));
    c(m + 1) = delta(m + 2) * delta(m + 3) / ((delta(m + 2) + delta(m + 3))^2) - 1;
    e(:,m + 1) = -1 / 3 * (originPoints(m+1, :) + 2 * originPoints(m, :));
%     f(m + 1) = -1 / 3 * (originPoints(m+1, 2) + 2 * originPoints(m, 2));

    for i = 2:m
        a(i) = delta(i + 3)^2 / (delta(i + 1) + delta(i + 2) + delta(i + 3));
        b(i) = delta(i + 3) * (delta(i + 1) + delta(i + 2)) / (delta(i + 1) + delta(i + 2) + delta(i + 3)) + delta(i + 2) * (delta(i + 3) + delta(i + 4)) / (delta(i + 2) + delta(i + 3) + delta(i + 4));
        c(i) = delta(i + 2)^2 / (delta(i + 2) + delta(i + 3) + delta(i + 4));
        e(:,i) = (delta(i + 2) + delta(i + 3)) * originPoints(i, :);
%         f(i) = (delta(i + 2) + delta(i + 3)) * originPoints(i, 2);
    end
%     disp('a=');
%     disp(a);
%     disp('b=');
%     disp(b);
%     disp('c=');
%     disp(c);
%     disp('e=');
%     disp(e);
%     disp('f=');
%     disp(f);

    % 构造系数矩阵
    matrix = zeros(m + 1, m + 1);
    matrix(1, 1) = a(1);
    matrix(1, 2) = b(1);
    matrix(1, 3) = c(1);
    matrix(m + 1, m - 1) = a(m + 1);
    matrix(m + 1, m) = b(m + 1);
    matrix(m + 1, m + 1) = c(m + 1);
    for i = 2:m
        matrix(i, i - 1) = a(i);
        matrix(i, i) = b(i);
        matrix(i, i + 1) = c(i);
    end
    disp('matrix=');
    disp(matrix);

    % 求逆
    matrix_inv = inv(matrix);
%     disp('matrix_inv=');
%     disp(matrix_inv);
    % 验证求逆结果
    % disp('求逆结果');
    % disp(matrix * matrix_inv);

    % 求控制点
%     ctrl_xs = matrix_inv * e';
%     ctrl_ys = matrix_inv * f';
    ctrlPoints = matrix_inv * e';
    % 在最前面添加，第一个控制点和第一个型值点重合
%     ctrl_xs = [originPoints(1, 1); ctrl_xs];
%     ctrl_ys = [originPoints(1, 2); ctrl_ys];
    ctrlPoints = [originPoints(1,:); ctrlPoints];
    % 在最后面添加，最后一个控制点和最后一个型值点重合
%     ctrl_xs = [ctrl_xs; originPoints(m+1, 1)];
%     ctrl_ys = [ctrl_ys; originPoints(m+1, 2)];
%     ctrlPoints = [ctrl_xs, ctrl_ys];
    ctrlPoints = [ctrlPoints; originPoints(end,:)];
%     disp('control points=');
%     disp(ctrlPoints);

    % 求插值点
    curvePoints = [];
    i = 4;
    u = 0;
    while u < 1
        tmp = zeros(1,col);
        if u > U(i + 1)
            i = i + 1;
        end
        for k = 0:3
%             t = [ctrlPoints(i - k, 1), ctrlPoints(i - k, 2)];
%             t(1) = t(1) * Nu(i - k, u, 3, U);
%             t(2) = t(2) * Nu(i - k, u, 3, U);
%             tmp(1) = tmp(1) + t(1);
%             tmp(2) = tmp(2) + t(2);
            t = [ctrlPoints(i-k,:)];
%             nu = Nu(i-k, u, 3, U)
            t = t * Nu(i-k, u, 3, U);
            tmp = tmp + t;
        end
        curvePoints = [curvePoints; tmp];
        % u = u + 0.01; % u可以自定义取值大小
        u = u + 1 / m / 20;
    end
    curvePoints = [curvePoints; originPoints(m+1, :)];
    
    % 返回控制点和插值点
%     disp('curve points=');
%     format long
%     disp(curvePoints);
end

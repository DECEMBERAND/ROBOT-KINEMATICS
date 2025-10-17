% 此文件为贝塞尔曲线插值算法，常用的为二次、三次贝塞尔 
% 二次：B(t) = (1-t)^2*P0 + 2*t*(1-t)*P1 + t^2*P2 t∈[0,1]
% 三次：B(t) = (1-t)^3*P0 + 3*t*(1-t)^2*P1 + 3*t^2*(1-t)*P2 + t^2*P3 t∈[0,1]
% 对于贝塞尔曲线，有三个要点：
%         其是通过控制点来生成的，控制点不全在最终曲线上。
%         控制点首末的两点是最终曲线的端点（意味着首末控制点会在最终曲线上），且各自与相邻点的连线同最终曲线相切。
%         两个贝塞尔曲线如果平滑连接，则需要连接点与其左右相邻两端点共线。
%参考博客：https://blog.csdn.net/deepsprings/article/details/107881698

points = [1, 1; 3, 6; 6, 3; 8, 0; 11, 6; 12, 12];
controlP = getControlPointList(points, -1, 1)';
l = size(points, 1) - 1;

figure;
t = linspace(0, 1, 50);

for i = 1:l
    p = [points(i, :); controlP(2*i-1, :); controlP(2*i, :); points(i+1, :)];
    interPoints = getInterpolationPoints(p, t);
    x = interPoints(:, 1);
    y = interPoints(:, 2);
    plot(x, y);
    hold on;
end

scatter(points(:, 1), points(:, 2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g');
hold off;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
title('Bezier Curve Interpolation');

function interPoints = getInterpolationPoints(controlPoints, tList)
    n = length(controlPoints) - 1;
    interPoints = zeros(length(tList), 2);
    
    for ti = 1:length(tList)
        t = tList(ti);
        Bt = zeros(1, 2);
        
        for i = 1:length(controlPoints)
            Bt = Bt + nchoosek(n, i-1) * (1-t)^(n-i+1) * t^(i-1) * controlPoints(i, :);
        end
        
        interPoints(ti, :) = Bt;
    end
end

function controlPoints = getControlPointList(pointsArray, k1, k2)
    points = pointsArray;
    index = size(points, 1) - 2;
    res = zeros(0, 2);
    
    for i = 1:index
        tmp = points(i:i+2, :);
        p1 = tmp(1, :);
        p2 = tmp(2, :);
        p3 = tmp(3, :);
        
        if k1 == -1
            l1 = norm(p1 - p2);
            l2 = norm(p2 - p3);
            k1 = l1 / (l1 + l2);
            k2 = l2 / (l1 + l2);
        end
        
        p01 = k1*p1 + (1-k1)*p2;
        p02 = (1-k2)*p2 + k2*p3;
        p00 = k2*p01 + (1-k2)*p02;
        
        sub = p2 - p00;
        p12 = p01 + sub;
        p21 = p02 + sub;
        
        res = [res; p12; p21];
    end
    
    pFirst = points(1, :) + 0.1 * (res(1, :) - points(1, :));
    pEnd = points(end, :) + 0.1 * (res(end, :) - points(end, :));
    res = [pFirst; res; pEnd];
    
    controlPoints = reshape(res', [], 2)';
end
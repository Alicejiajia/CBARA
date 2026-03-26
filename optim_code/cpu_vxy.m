function [v_x,v_y] = cpu_vxy(start,ed,T,N)
%输入起点和终点坐标，计算沿x,y轴分量的速度
% start=起点坐标,end=终点坐标,T=采样时间,N=采样点数
%计算合成速度v
d = norm(start-ed);
v = d/(T*N);
x = ed(1,1) - start(1,1);
y = ed(1,2) - start(1,2);
sin_th = y/d;
cos_th = x/d;
%计算分量速度
v_x = v*cos_th;
v_y = v*sin_th;
end


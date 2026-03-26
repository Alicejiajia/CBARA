function [h] = cpu_h(Xi,base)
%通过状态计算观测值,Xi为4*1的矩阵,base为基站坐标1*2
% 
xb=base(1,1)
yb= base(1,2)
x = Xi(1)
y = Xi(2)
vx = Xi(3)
vy = Xi(4)
d = sqrt((x-xb)^2 + (y-yb)^2)
v = ((x-xb)*vx+(y-yb)*vy)/d
theta = atan((y-yb)/(x-xb))
h = [d theta v].'
end


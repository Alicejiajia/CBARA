function [Jacobih] = cpu_Jacobih(Xi,base)
%计算估计参数向量的雅克比矩阵
syms x y vx vy
sy = [x y vx vy]
d = sqrt((base(1)-x)^2+(base(2)-y)^2)
v = (vx*(base(1)-x) + vy*(base(2)-y))/d
theta = atan((base(2)-y)/(base(1)-x))
h = [d ; v ; theta]
J = jacobian(h,sy)
jx = subs(J,sy,Xi)
Jacobih = double(jx)
% [J,Jf,var,jx] = Jacobi(h,Xi)

end


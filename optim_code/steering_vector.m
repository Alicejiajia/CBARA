function [a] = steering_vector(Nt,d1,d2)
%用于形成传输和接收阵列矢量,d11为接收目标位置，d22为发射目标位置

a = zeros(1,Nt)
thta = atan((d1(1,2)-d2(1,2))/(d1(1,1)-d2(1,1)))
%thta=rad2deg(thta)
for i=1:Nt
    a(1,i) = (1/Nt)^(1/2)*exp(j*pi*sin(thta)*(i-1))
end
a = a.'
end


function [s_T1,s_T2,s_T3, s_T4, s_T5] = state_init(dist,N,T,BS)
%对三个目标的轨迹进行初始化，N为采样时隙个数
%   此处显示详细说明

%可行的解决方案N=30，****************************三个用户轨迹1************************************
%T1
% dist=100
% x1 = linspace(-dist/2,0,N).'
% y1 = linspace(dist/4*sqrt(3),-dist/4*sqrt(3)+15,N).'
% vx1 = linspace(3.3,3.3,N).'
% vy1 = linspace(5.8,5.8,N).'
% t_T1 = [x1 y1]
% s_T1 = [x1 y1 vx1 vy1]
% 
% %T2
% x2 = linspace(55,55,N).'
% y2 = linspace(10,70,N).'
% t_T2 = [x2 y2]
% vx2 = linspace(0,0,N).'
% vy2 = linspace(4,4,N).'
% s_T2 = [x2 y2 vx2 vy2]
% 
% x3 = linspace(-60,60,N).'
% y3 = linspace(0,0,N).'
% t_T3 = [x3 y3]
% vx3 = linspace(8,8,N).'
% vy3 = linspace(0,0,N).'
% s_T3 = [x3 y3 vx3 vy3]
% 
% plot(t_T1(:,1),t_T1(:,2),'ro','linewidth',1);
% text(t_T1(1,1)+3,t_T1(1,2)-3,'T1-sens')
% plot(t_T2(:,1),t_T2(:,2),'bo','linewidth',1);
% text(t_T2(1,1)+3,t_T2(1,2),'T2-sens')
% plot(t_T3(:,1),t_T3(:,2),'go','linewidth',1);
% text(t_T3(1,1)-3,t_T3(1,2)-3,'T3-isac')
%**************************************************扩大范围
% target 1
x1 = linspace(-dist/2,0,N).'
y1 = linspace(dist/4*sqrt(3),-dist/4*sqrt(3)+15,N).'
start1 = [-dist/2, dist/4*sqrt(3)] %计算速度分量
ed1 = [0, -dist/4*sqrt(3)+15]
[vx1,vy1] = cpu_vxy(start1,ed1,T,N)
vx1 = linspace(vx1,vx1,N).'
vy1 = linspace(vy1,vy1,N).'
t_T1 = [x1 y1]
s_T1 = [x1 y1 vx1 vy1]

%Target 2
x2 = linspace(50,-40,N).'
y2 = linspace(90,70,N).'
t_T2 = [x2 y2]
start2 = [50, 90] %计算速度分量
ed2 = [-40, 70]
[vx2,vy2] = cpu_vxy(start2,ed2,T,N)
vx2 = linspace(vx2,vx2,N).'
vy2 = linspace(vy2,vy2,N).'
s_T2 = [x2 y2 vx2 vy2]

%targrt 3
x3 = linspace(66,66,N).'
y3 = linspace(12,84,N).'
t_T3 = [x3 y3]
start3 = [66, 12] %计算速度分量
ed3 = [66, 84]
[vx3,vy3] = cpu_vxy(start3,ed3,T,N)
vx3 = linspace(vx3,vx3,N).'
vy3 = linspace(vy3,vy3,N).'
s_T3 = [x3 y3 vx3 vy3]

%ISAC 4
x4 = linspace(-72,72,N).'
y4 = linspace(0,0,N).'
t_T4 = [x4 y4]
start4 = [-72, 0] %计算速度分量
ed4 = [72, 0]
[vx4,vy4] = cpu_vxy(start4,ed4,T,N)
vx4 = linspace(vx4,vx4,N).'
vy4 = linspace(vy4,vy4,N).'
s_T4 = [x4 y4 vx4 vy4]

%ISAC 5
x5 = linspace(-50,80,N).'
y5 = linspace(-80,-20,N).'
t_T5 = [x5 y5]
start5 = [-50, -80] %计算速度分量
ed5 = [80, -20]
[vx5,vy5] = cpu_vxy(start5,ed5,T,N)
vx5 = linspace(vx5,vx5,N).'
vy5 = linspace(vy5,vy5,N).'
s_T5 = [x5 y5 vx5 vy5]

% plot(t_T1(1,1),t_T1(1,2),'.','MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerEdgeColor',[0.4940 0.1840 0.5560],'markersize',18);
p1 = plot(t_T1(:,1),t_T1(:,2),'--d','color',[0.4940 0.1840 0.5560],'linewidth',1.5,'MarkerIndices',[1 N],'markersize',4);
text(t_T1(1,1)-8,t_T1(1,2)+10,'Target 1','FontName','Times New Roman','FontSize',11);
% plot(t_T1(N,1),t_T1(N,2),'.','MarkerFaceColor',[0.4940 0.1840 0.5560]','MarkerEdgeColor',[0.4940 0.1840 0.5560],'markersize',18);
%坐标归一化
x11 = [(t_T1(1,1)+25+100)/200,(t_T1(10,1)+25+100)/200]
y11 = [(t_T1(1,2)+100-5)/200,(t_T1(10,2)+100-5)/200]
annotation('arrow',x11,y11,'color',[0.4940 0.1840 0.5560]) 

% plot(t_T2(1,1),t_T2(1,2),'.','MarkerFaceColor', [0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],'markersize',18);
p2 = plot(t_T2(:,1),t_T2(:,2),'--d','color',[0 0.4470 0.7410],'linewidth',1.5,'MarkerIndices',[1 N],'markersize',4);
text(t_T2(1,1)+5,t_T2(1,2)+2,'Target 2','FontName','Times New Roman','FontSize',11);
% plot(t_T2(N,1),t_T2(N,2),'.','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],'markersize',18);
x22 = [(t_T2(1,1)+100-5)/200,(t_T2(10,1)+100-5)/200]
y22 = [(t_T2(1,2)+100-20)/200,(t_T2(10,2)+100-20)/200]
annotation('arrow',x22,y22,'color',[0 0.4470 0.7410]) 

p3 = plot(t_T3(:,1),t_T3(:,2),'--d','color',[0.85 0.325 0.098],'linewidth',1.5,'MarkerIndices',[1 N],'markersize',4);
text(t_T3(1,1)+7,t_T3(1,2)+10,'Target 3','FontName','Times New Roman','FontSize',11);
% plot(t_T2(N,1),t_T2(N,2),'.','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],'markersize',18);
x33 = [(t_T3(1,1)+100-5.5)/200,(t_T3(10,1)+100-6)/200]
y33 = [(t_T3(1,2)+100+15)/200,(t_T3(10,2)+100+15)/200]
annotation('arrow',x33,y33,'color',[0.85 0.325 0.098]) 

% plot(t_T3(1,1),t_T3(1,2),'.','MarkerFaceColor', [0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880]','markersize',18);
p4 = plot(t_T4(:,1),t_T4(:,2),'--d','color',[0.4660 0.6740 0.1880],'linewidth',1.5,'MarkerIndices',[1 N],'markersize',4);
text(t_T4(1,1)-5,t_T4(1,2)-8,'ISAC 4','FontName','Times New Roman','FontSize',11);
% plot(t_T3(N,1),t_T3(N,2),'.','MarkerFaceColor', [0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880],'markersize',18);
x44 = [(t_T4(1,1)+100+70)/200,(t_T4(10,1)+100+70)/200]
y44 = [(t_T4(1,2)+100-5)/200,(t_T4(10,2)+100-5)/200]
annotation('arrow',x44,y44,'color',[0.4660 0.6740 0.1880]) 

p5 = plot(t_T5(:,1),t_T5(:,2),'--d','color',[0.301 0.745 0.933],'linewidth',1.5,'MarkerIndices',[1 N],'markersize',4);
text(t_T5(1,1)-15,t_T5(1,2)-8,'ISAC 5','FontName','Times New Roman','FontSize',11);
% plot(t_T3(N,1),t_T3(N,2),'.','MarkerFaceColor', [0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880],'markersize',18);
x55 = [(t_T5(1,1)+100+20)/200,(t_T5(10,1)+100+20)/200]
y55 = [(t_T5(1,2)+100+10)/200,(t_T5(10,2)+100+10)/200]
annotation('arrow',x55,y55,'color',[0.301 0.745 0.933]) 

% legend([BS,p1 p2 p3],{'Base Station','Target 1','Target 2','Target 3'},'FontName','Times New Roman','FontSize',7);
%************************************************4个用户轨迹2
% dist=100
% x1 = linspace(40,-40,N).'
% y1 = linspace(-100,80,N).'
% vx1 = linspace(3.3,3.3,N).'
% vy1 = linspace(5.8,5.8,N).'
% t_T1 = [x1 y1]
% s_T1 = [x1 y1 vx1 vy1]
% 
% %T2
% x2 = linspace(80,80,N).'
% y2 = linspace(-60,110,N).'
% t_T2 = [x2 y2]
% vx2 = linspace(0,0,N).'
% vy2 = linspace(4,4,N).'
% s_T2 = [x2 y2 vx2 vy2]
% 
% x3 = linspace(-60,60,N).'
% y3 = linspace(dist/4*sqrt(3),dist/4*sqrt(3),N).'
% t_T3 = [x3 y3]
% vx3 = linspace(8,8,N).'
% vy3 = linspace(0,0,N).'
% s_T3 = [x3 y3 vx3 vy3]
% 
% x4 = linspace(40,-60,N).'
% y4 = linspace(60,-80,N).'
% t_T4 = [x4 y4]
% vx4 = linspace(8,8,N).'
% vy4 = linspace(0,0,N).'
% s_T4 = [x4 y4 vx4 vy4]
% 
% plot(t_T1(:,1),t_T1(:,2),'ro','linewidth',1);
% text(t_T1(1,1)+3,t_T1(1,2)-3,'T1-sens')
% plot(t_T2(:,1),t_T2(:,2),'bo','linewidth',1);
% text(t_T2(1,1)+3,t_T2(1,2),'T2-sens')
% plot(t_T3(:,1),t_T3(:,2),'go','linewidth',1);
% text(t_T3(1,1)-3,t_T3(1,2)-3,'T3-isac')
% plot(t_T4(:,1),t_T4(:,2),'mo','linewidth',1);
% text(t_T4(1,1)-3,t_T3(1,2)-3,'T4-sens')

end


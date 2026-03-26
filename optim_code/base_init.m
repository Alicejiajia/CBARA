function [base,BS] = base_init(dist)
%*******************************三个基站**********************************************
base1=[0 dist/4*sqrt(3)];
base2=[dist/2 -dist/4*sqrt(3)];
base3=[-dist/2 -dist/4*sqrt(3)];
figure
hold on
set(gca, 'box', 'on','FontName','Times New Roman','FontSize',11,'FontWeight','bold');
BS = plot(base1(1),base1(2),'p','MarkerFaceColor', [0.9290 0.6940 0.1250],'MarkerEdgeColor',[0.9290 0.6940 0.1250],'markersize',15,'LineWidth',1.5);
text(base1(1)-5,base1(2)-10,'BS 1','FontName','Times New Roman','FontSize',11);
plot(base2(1),base2(2),'p','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor',[0.9290 0.6940 0.1250],'markersize',15','LineWidth',1.5);
text(base2(1)-5,base2(2)-10,'BS 2','FontName','Times New Roman','FontSize',11);
plot(base3(1),base3(2),'p','MarkerFaceColor', [0.9290 0.6940 0.1250],'MarkerEdgeColor',[0.9290 0.6940 0.1250],'markersize',15,'LineWidth',1.5);
text(base3(1)-5,base3(2)-10,'BS 3','FontName','Times New Roman','FontSize',11);
base = [base1;base2;base3];
ylim([-100 100])
xlim([-100 100])
xlabel('X axis (m)')
ylabel('Y axis (m)')
end


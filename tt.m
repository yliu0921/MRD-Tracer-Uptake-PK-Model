%计算cellnum和RK的趋势

clc;close all
clear
cellNum = 100000;             % 小圆球数量
d = 8e-6;            % 小圆球直径 (单位: 米)
r_inner = 8e-6;      % 内圆柱半径 (单位: 米)
fill_efficiency = 0.5; % 填充效率

Rk = cal_Rk(cellNum,8e-6,8e-6,0.5);
disp(['外圆半径: ', num2str(Rk), ' m']);
%%
Rk = [];
elem = 20;
cellNum_array = logspace(2,8, elem);
for i = 1:elem
    cellNum = cellNum_array(i);
    Rk(i) = cal_Rk(cellNum,8e-6,8e-6,0.5)*10^6;
    
end

figure(1)
plot(cellNum_array,Rk,'LineWidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
xlabel('Cell Number (个)')
ylabel('Rk (um)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Tumor Size')
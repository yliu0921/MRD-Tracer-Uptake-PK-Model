%tumor pk model
%绘制三维时间曲线
clc; clear; close all
warning('off')

%Description:
%INPUT PARAMETERS PAYLOAD
L = 0 ;
LIDcircV = 15E-8;
NR = 1E3;%Intrautmoral Receptor Density on Cells - [receptors/cell]

Vb = 2E-3;%L,2mL

%Receptor Rate Constants
kon_R = 1.00E5;                     % [1/M/s] - on rate for IL-2 and IL-2RB from NKTR paper
koff_R = 0.001;                      % [1/s] - off rate for IL-2 and IL-2RB NKTR paper
kendo_R = 0.3/(60);                % [1/s] - endocytic rate of IL-2R with ligand

% ODE solver options
options = odeset('RelTol',1e-14,'AbsTol',[1e-14]);
tspan = [0 1*48*60*60];
tplot = linspace(0,1*24,1000);

%% 绘制三维图像
% figure(1);
% for Group_Num = 1:elem
%
%     x = tplot; % 数据长度为100，将全部显示出来
%
%     y = koff_R_array(Group_Num) * ones(size(x));  % 共elem组数据(组别号：1:elem)
%
%     z = ID_L(1:1000,Group_Num); % 每组数据的曲线值
%
%     plot3(x,y,z);   % 绘图
%     hold on;        % 保持
% end
%
% grid on;
% title(['多组样本值变化曲线']);
% xlabel("时间");
%
% ylabel("tumor size");zlabel("曲线值");

%%
%%%% Affinity and time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elem =150;                           %number of elements in the array iteration
% MW_array = linspace(1, 150, elem);    %k iterations through for-loop
% cellNum_array = logspace(2,5 , elem);    %k iterations through for-loop

cellNum = 10000;
koff_R_array = logspace(-9,-2,elem);
Kd_array = koff_R_array/kon_R;   %g iterations through for-loop
figure(1)
MW = [0.5,7,15,150];
for i = 1:4
    for k = 1:elem
        koff_R = koff_R_array(k);
        [p, y0] = Inputs(MW(i),[],NR,kon_R,koff_R,kendo_R,cellNum,L);
        [t,y] = ode15s(@odefun,tspan,y0,options,p);%dt from tspa
        %         L_tumor = (y(:,1) + y(:,3))/LIDcircV;
        L_tumor = ((y(:,1) + y(:,3)))/LIDcircV;
        Ltumor = interp1(t/(60*60),L_tumor,tplot);
        ID_L(1:1000,k) = Ltumor;


    end


    subplot(2,2,i)
    [Time,KD]=meshgrid(tplot,Kd_array);
    surf(Time,KD,ID_L'.*100)
    shading flat
    shading flat
    colormap('parula')
    cmap = colorbar;
    set(gca,'XScale','linear')
    set(gca,'YScale','log')
    set(gca, 'fontsize', 18)
% zlim([0, 10])
    xlabel('Time (hr)')
    ylabel('Affinity (M)')
    zlabel('%ID/g')

    set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
    title(['Molecular Weight = ', num2str(MW(i)),' kDa'])

end


%%
%%%% CellNum and time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
elem =150;                           %number of elements in the array iteration
% MW_array = linspace(1, 150, elem);    %k iterations through for-loop
MW = 20;
cellNum_array = logspace(2,5 , elem);    %k iterations through for-loop

% cellNum = 2500;
% koff_R_array = logspace(-9,-2,elem);
% Kd_array = koff_R_array/kon_R;   %g iterations through for-loop
figure(2)
koff_R = [10E-8 10E-4 10E-2 10E0];
for i = 1:4
    for k = 1:elem
        cellNum = cellNum_array(k);
        [p, y0] = Inputs(MW,[],NR,kon_R,koff_R(i),kendo_R,cellNum,L);
        [t,y] = ode15s(@odefun,tspan,y0,options,p);%dt from tspa
        %         L_tumor = (y(:,1) + y(:,3))/LIDcircV;
        L_tumor = ((y(:,1) + y(:,3)))/LIDcircV;
        Ltumor = interp1(t/(60*60),L_tumor,tplot);
        ID_L(1:1000,k) = Ltumor;


    end

    subplot(2,2,i)
    [Time,CNum]=meshgrid(tplot,cellNum_array);
    surf(Time,CNum,ID_L'.*100)
    shading flat
    shading flat
    colormap('parula')
    cmap = colorbar;
    set(gca,'XScale','linear')
    set(gca,'YScale','log')
    set(gca, 'fontsize', 18)
% zlim([0, 10])
    xlabel('Time (hr)')
    ylabel('Cell Number (个)')
    zlabel('%ID/g')

    set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
    title(['MW = ', num2str(MW),' kDa','   Affinity = ', num2str(koff_R(i)/kon_R),' M'])

end
%%
%%%% MW and time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elem =150;                           %number of elements in the array iteration
MW_array = linspace(1, 150, elem);    %k iterations through for-loop
% cellNum_array = logspace(2,5 , elem);    %k iterations through for-loop

cellNum = 10000;
koff_R_array = logspace(-9,-2,elem);
Kd_array = koff_R_array/kon_R;   %g iterations through for-loop
figure(3)
koff_R = [10E-8 10E-4 10E-2 10E0];
for i = 1:4
    for k = 1:elem
        MW = MW_array(k);
        [p, y0] = Inputs(MW,[],NR,kon_R,koff_R(i),kendo_R,cellNum,L);
        [t,y] = ode15s(@odefun,tspan,y0,options,p);%dt from tspa
        %         L_tumor = (y(:,1) + y(:,3))/LIDcircV;
        L_tumor = ((y(:,1) + y(:,3)))/LIDcircV;
        Ltumor = interp1(t/(60*60),L_tumor,tplot);
        ID_L(1:1000,k) = Ltumor;


    end


    subplot(2,2,i)
    [Time,MW]=meshgrid(tplot,MW_array);
    surf(Time,MW,ID_L'.*100)
    shading flat
    shading flat
    colormap('parula')
    cmap = colorbar;
    set(gca,'XScale','linear')
    set(gca,'YScale','linear')
    set(gca, 'fontsize', 18)
% zlim([0, 10])
    xlabel('Time (hr)')
    ylabel('Molecular Weight (kDa)')
    zlabel('%ID/g')

    set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
    title(['Affinity = ', num2str(koff_R(i)/kon_R),' M'])

end

%%
%%%% CellNum and time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
elem =150;                           %number of elements in the array iteration
% MW_array = linspace(1, 150, elem);    %k iterations through for-loop

cellNum_array = logspace(2,5 , elem);    %k iterations through for-loop
koff_R = 10E-7;
% cellNum = 2500;
% koff_R_array = logspace(-9,-2,elem);
% Kd_array = koff_R_array/kon_R;   %g iterations through for-loop
figure(4)
MW = [0.5,7,15,150];
for i = 1:4
    for k = 1:elem
        cellNum = cellNum_array(k);

        [p, y0] = Inputs(MW(i),[],NR,kon_R,koff_R,kendo_R,cellNum,L);
        [t,y] = ode15s(@odefun,tspan,y0,options,p);%dt from tspa
        %         L_tumor = (y(:,1) + y(:,3))/LIDcircV;
        L_tumor = ((y(:,1) + y(:,3)))/LIDcircV;
        Ltumor = interp1(t/(60*60),L_tumor,tplot,'pchip');
        ID_L(1:1000,k) = Ltumor;


    end

    subplot(2,2,i)
    [Time,CNum]=meshgrid(tplot,cellNum_array);
    surf(Time,CNum,ID_L'.*100)
    shading flat
    shading flat
    colormap('parula')
    cmap = colorbar;
    set(gca,'XScale','linear')
    set(gca,'YScale','log')
    set(gca, 'fontsize', 18)
% zlim([0, 10])
    xlabel('Time (hr)')
    ylabel('Cell Number (个)')
    zlabel('%ID/g')

    set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
    title(['MW = ', num2str(MW(i)),' kDa','   Affinity = ', num2str(koff_R/kon_R),' M'])

end
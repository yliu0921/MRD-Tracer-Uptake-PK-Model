clc; clear; close all

%INPUT PARAMETERS PAYLOAD
L = 0 ;
LIDcircV = 15E-8;
%Intrautmoral Receptor+ Cells - [cells/mL]
% cells = 150E3;
%Intrautmoral Receptor Density on Cells - [receptors/cell]
NR = 1000;

%Receptor Rate Constants
kon_R = 1E2;                     % [1/M/s] - on rate for IL-2 and IL-2RB from NKTR paper
koff_R = 7E-4;                      % [1/s] - off rate for IL-2 and IL-2RB NKTR paper
kendo_R = 0.3/(60);                % [1/s] - endocytic rate of IL-2R with ligand

% ODE solver options
options = odeset('RelTol',1e-14,'AbsTol',[1e-14]);
tspan = [0 1*24*60*60];                        %1d，24hours

elem =9;                           %number of elements in the array iteration
% MW_array = linspace(1,150 , elem);    %k iterations through for-loop
MW_array = [0.5,1,7,15,27,55,80,110,150];    %k iterations through for-loop
% cellNum_array = logspace(2,5 , elem);
cellNum = 25000; %2500-75,um
% MW= 15;
disp(num2str(koff_R/kon_R))
figure(1)
subplot(2,1,1)
%tumor
for k = 1:elem
    MW = MW_array(k);
% cellNum= cellNum_array(k);
        [p, y0] = Inputs(MW,[],NR,kon_R,koff_R,kendo_R,cellNum,L);
        [t,y] = ode15s(@odefun_new,tspan,y0,options,p);%dt from tspan
        ID_L = (y(:,1) + y(:,3))/LIDcircV*100;%tumor

        plot(t/(60*60), ID_L,'LineWidth',2)
        legend(num2str(MW_array')+" kDa")
        title('Tumor TAC')
        xlabel('Time (hr)')
        ylabel('% Injected Dose')
        xlim([-1, 24])
%         ylim([0, 80])
        set(gca, 'fontsize', 18)
        set(gca,'LineWidth',1.5,'TickLength',[0.02 0.02]);
        hold on
end

% subplot(3,1,2)
% %tumor
% for k = 1:elem
%     MW = MW_array(k);
% % cellNum= cellNum_array(k);
%         [p, y0] = Inputs(MW,[],NR,kon_R,koff_R,kendo_R,cellNum,L);
%         [t,y] = ode15s(@odefun_new,tspan,y0,options,p);%dt from tspan
%         ID_L = (y(:,2) )/LIDcircV*100;%tumor
% 
%         plot(t/(60*60), ID_L,'LineWidth',2)
%         legend(num2str(MW_array')+" kDa")
%         title('Blood TAC')
%         xlabel('Time (hr)')
%         ylabel('% Injected Dose')
%         xlim([-1, 24])
% %         ylim([0, 80])
%         set(gca, 'fontsize', 18)
%         set(gca,'LineWidth',1.5,'TickLength',[0.02 0.02]);
%         hold on
% end

% figure(2)
subplot(2,1,2)
for k = 1:elem
    MW = MW_array(k);
% cellNum= cellNum_array(k);
        [p, y0] = Inputs(MW,[],NR,kon_R,koff_R,kendo_R,cellNum,L);
        [t,y] = ode15s(@odefun_new,tspan,y0,options,p);%dt from tspan
        L_circL = ((y(:,1)+y(:,3))-y(:,2))/LIDcircV*100;%blood
      % L_circL = ((y(:,1)+y(:,3))/y(:,2))/LIDcircV*100;%blood
        plot(t/(60*60), L_circL,'LineWidth',2)
        legend(num2str(MW_array')+" kDa")
%         legend(num2str(round(MW_array)')+" kDa")
        title('Tumor - Blood TAC')
        xlabel('Time (hr)')
        ylabel('% Injected Dose')
        xlim([-1, 24])
        ylim([-1, 20])
        set(gca, 'fontsize', 18)
        set(gca,'LineWidth',1.5,'TickLength',[0.02 0.02]);
        hold on
end

%tumor pk model

clc; clear; close all
warning('off')

%Description:
%INPUT PARAMETERS PAYLOAD
L = 0 ;
LIDcircV = 15E-8;
%INPUT PARAMETERS RECEPTOR (QSS)

cellNum = 2500;
%Intrautmoral Receptor Density on Cells - [receptors/cell]
NR = 1000;

%Receptor Rate Constants
kon_R = 1.00E5;                     % [1/M/s] - on rate 
koff_R = 0.001;                      % [1/s] - off rate
kendo_R = 0.3/(60);                % [1/s] - endocytic rate of IL-2R with ligand

%ODE solver options
options = odeset('RelTol',1e-14,'AbsTol',[1e-14]);
tspan = [0 1*24*60*60];                        %time span

%% %Molecular Weight, Affinity and Fractional Activity Durations Arrays
elem =17;                           %number of elements in the array iteration
MW_array = linspace(1,150,elem);    %k iterations through for-loop
koff_R_array = logspace(-9,-2,elem);
Kd_array = koff_R_array/kon_R;   %g iterations through for-loop

ID_L = ones(elem);
ID_Ltime = ones(elem);
Act_array = ones(elem);
Sys_Act_array = ones(elem);
Exp_array = ones(elem);

for k = 1:elem
   MW = MW_array(k);

   for g = 1:elem
       koff_R = koff_R_array(g);
       [p, y0] = Inputs(MW,[],NR,kon_R,koff_R,kendo_R,cellNum,L);
       [t,y] = ode15s(@ODEs,tspan,y0,options,p);%dt from tspan
       [maxval, maxloca] = max((y(:,1) + y(:,3)-y(:,2))/LIDcircV);
       ID_L(k,g) = maxval.*100;%tumor blood contrast max %ID
       ID_Ltime(k,g) = t(maxloca)/(60*60);%time to reach max contrast
       Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);%y1y3,[L,LC]
       Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);%y4[Lesc]
        
    end
end

figure(1)
subplot(211)%Figure4b
contourf(Kd_array,MW_array,Act_array)
colormap('jet')
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Molecular Weight (kDa)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Tumor Activity Durations (hr)')

% hold on
% % scatter([15.29E-12, 3.71E-12], [15, 70],100, 'ok','MarkerFaceColor', 'w')
% scatter([649.8E-12, 751.5E-12], [20, 70],100, 'ok','MarkerFaceColor', 'w')
% hold off

subplot(212)%Figure4c
contourf(Kd_array,MW_array,Sys_Act_array)
colormap('jet')
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Molecular Weight (kDa)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Blood Durations (hr)')


%% %Molecular Weight, Affinity and Peak %Uptake Arrays

figure(2)
subplot(211)%Figure4b
contourf(Kd_array,MW_array,ID_L)
colormap('jet')
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Molecular Weight (kDa)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Peak Tumor-blood Contrast (%ID)')


subplot(212)%Figure4c
contourf(Kd_array,MW_array,ID_Ltime)
colormap('jet')
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Molecular Weight (kDa)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('the Time Reach Peak Contrast (hr)')

%Activity duration(in hours) in Tumor Function
function cA = tumorA(t1,y1,kon_R, koff_R,tspan)
EC50tumor = koff_R/kon_R;
L_tumor = y1(:,1) + y1(:,3);%from bonded and free
Activity_L = 1./(1+EC50tumor./L_tumor);
cA = trapz(t1, Activity_L)/(60*60);%in hours
end

%Activity duration(in hours) in Circulation Function
function cA = sysA(t1,y1,kon_R,koff_R,tspan)
EC50sys = koff_R/kon_R;
L_sys = y1(:,2); %from Circulation
Activity_sys = 1./(1+EC50sys./L_sys);
cA = trapz(t1, Activity_sys)/(60*60);%in hours
end
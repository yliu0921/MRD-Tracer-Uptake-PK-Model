%tumor pk model

clc; clear; close all
warning('off')

%Description:
%INPUT PARAMETERS PAYLOAD
L = 0 ;
LIDcircV = 15E-8;
%INPUT PARAMETERS RECEPTOR (QSS)
%Intrautmoral Receptor+ Cells - [cells/mL]
cells = 150E3;
%Intrautmoral Receptor Density on Cells - [receptors/cell]
NR = 1000;
Vb = 2E-3;%L,2mL
MW = 1.04212; %kDa 
% Vt = 524E-12;%L,100uL

%Receptor Rate Constants
kon_R = 1E14;                     % [1/M/s] - on rate for IL-2 and IL-2RB from NKTR paper
koff_R = 7E4;                      % [1/s] - off rate for IL-2 and IL-2RB NKTR paper
kendo_R = 0.3/(60);                % [1/s] - endocytic rate of IL-2R with ligand

%ODE solver options
options = odeset('RelTol',1e-14,'AbsTol',[1e-14]);
tspan = [0 1*24*60*60];                        %3d，72hours

%%
%%%% Affinity and CellNum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% % % % %Molecular Weight, Affinity and Fractional Activity Arrays
elem =20;                           %number of elements in the array iteration
% MW_array = logspace(0, 3, elem);    %k iterations through for-loop
cellNum_array = logspace(2,5,elem);    %k iterations through for-loop
% Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop

koff_R_array = logspace(-9,-2,elem);
% Rk_array = linspace(20e-6,200e-6,elem);  %20um-200um

Kd_array = koff_R_array/kon_R;   %g iterations through for-loop

ID_L = ones(elem);
ID_Ltime = ones(elem);
Act_array = ones(elem);
Sys_Act_array = ones(elem);
Exp_array = ones(elem);
for k = 1:elem
    cellNum = cellNum_array(k);

    for g = 1:elem
        koff_R = koff_R_array(g);
        % Rk = Rk_array(g);
        %         disp(k)
        [p, y0] = Inputs(20,[],NR,kon_R,koff_R,kendo_R,cellNum,L);
        %[p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
        %p = [k1, k2, k3, k4, k5, k6, k12, k13, k14, k15, A, k16,k17];
        %y0 = [L,C,LC,Lesc,Cint, LCint, LR, LCR, LRint,Reng,R];

        [t,y] = ode15s(@odefun,tspan,y0,options,p);%dt from tspan

        [maxval, maxloca] = max((y(:,1) + y(:,3))/LIDcircV*100);%%tumor max ID
        ID_L(k,g) = maxval;
        ID_Ltime(k,g) = t(maxloca)/(60*60);
        Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);%y1y3y7y8,[L,LC, LR, LCR]
        Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);%y4[Lesc]
        Exp_array(k,g) = tumorExposure(t,y,kon_R,koff_R,tspan);
        
    end
end

figure(1)
subplot(211)
% s = mesh(Kd_array,MW_array,Act_array);
contourf(Kd_array,cellNum_array,Act_array)
colormap('jet')
% caxis([0, 1])
% shading flat
% s.FaceColor = 'flat';
% shading interp
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Cell Number (个)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Tumor Retention (hr)')

% hold on
% % scatter([15.29E-12, 3.71E-12], [15, 70],100, 'ok','MarkerFaceColor', 'w')
% scatter([649.8E-12, 751.5E-12], [20, 70],100, 'ok','MarkerFaceColor', 'w')
% hold off

subplot(212)
% sl = surfl(Kd_array,MW_array,Sys_Act_array);
contourf(Kd_array,cellNum_array,Sys_Act_array)
% sl.EdgeColor = 'none';
% shading flat
% shading interp
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Cell Number (个)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Blood Retention (hr)')
% hold on
% scatter([649.8E-12, 751.5E-12], [20, 70],100, 'ok','MarkerFaceColor', 'w')
% hold off

figure(2)

figure(2)
subplot(211)
% s=mesh(Kd_array,MW_array,ID_L);
contourf(Kd_array,cellNum_array,ID_L)
colormap('jet')
% shading flat
% s.FaceColor = 'flat';
% shading interp
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Cell Number (个)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Peak Tumor Uptake (%ID)')

subplot(212)
% s=mesh(Kd_array,MW_array,ID_Ltime);
contourf(Kd_array,cellNum_array,ID_Ltime)
% shading flat
% s.FaceColor = 'flat';
% shading interp
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Cell Number (个)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('the Time of Peaktumor Uptake (hr)')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Supplemental
%         [p, y0_LLM] = Inputs(14E-9,137,[],Th_C,0,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [tLLM,yLLM] = ode15s(@odefun,tspan,y0_LLM,options,p);
%
%         [p, y0_LLxM] = Inputs(130E-9,137,[],Th_C,0,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [tLLxM,yLLxM] = ode15s(@odefun,tspan,y0_LLxM,options,p);
%
%         [p, y0_LxLxM] = Inputs(5E-6,137,[],Th_C,0,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [tLxLxM,yLxLxM] = ode15s(@odefun,tspan,y0_LxLxM,options,p);
%
%         [p, y0_L] = Inputs(221E-9,32,[],Th_C,0,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [tL,yL] = ode15s(@odefun,tspan,y0_L,options,p);
%
%         [p, y0_Lx] = Inputs(5E-6,32,[],Th_C,0,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [tLx,yLx] = ode15s(@odefun,tspan,y0_L,options,p);
%
%         L_tumorLLM = yLLM(:,1) + yLLM(:,3) + yLLM(:,7) + yLLM(:,8);
%         L_tumorLLxM = yLLxM(:,1) + yLLxM(:,3) + yLLxM(:,7) + yLLxM(:,8);
%         L_tumorLxLxM = yLxLxM(:,1) + yLxLxM(:,3) + yLxLxM(:,7) + yLxLxM(:,8);
%         L_tumorL = yL(:,1) + yL(:,3) + yL(:,7) + yL(:,8);
%         L_tumorLx = yLx(:,1) + yLx(:,3) + yLx(:,7) + yLx(:,8);
%
%         L0 = 1E-6;
%         ID_LLM = L_tumorLLM/(L0)*100;
%         ID_LLxM = L_tumorLLxM/(L0)*100;
%         ID_LxLxM = L_tumorLxLxM/(L0)*100;
%         ID_L = L_tumorL/(L0)*100;
%         ID_Lx = L_tumorLx/(L0)*100;
%
%         figure
%         subplot(5,1,1), plot(tLLM/(60*60), ID_LLM)
%         ylim([0, 100])
%         title('Supplemental Figure 7a')
%         subplot(5,1,2), plot(tLLxM/(60*60), ID_LLxM)
%         ylim([0, 100])
%         subplot(5,1,3), plot(tLxLxM/(60*60), ID_LxLxM)
%         ylim([0, 100])
%         subplot(5,1,4), plot(tL/(60*60), ID_L)
%         ylim([0, 100])
%         subplot(5,1,5), plot(tLx/(60*60), ID_Lx)
%         ylim([0, 100])

%%
%%%% CellNUm and MW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% % % % %Molecular Weight, Affinity and Fractional Activity Arrays
elem =50;                           %number of elements in the array iteration
MW_array = logspace(0, 2, elem);    %k iterations through for-loop
cellNum_array = logspace(2,5 , elem);    %k iterations through for-loop

ID_L = ones(elem);
ID_Ltime = ones(elem);
Act_array = ones(elem);
Sys_Act_array = ones(elem);
Exp_array = ones(elem);

disp(['Affinity: ', num2str((koff_R/kon_R)*10^9),' nM']);

for k = 1:elem
    cellNum = cellNum_array(k);
    for g = 1:elem
        MW = MW_array(g);
        % Rk = Rk_array(g);
        %disp(k)
        [p, y0] = Inputs(MW,[],NR,kon_R,koff_R,kendo_R,cellNum,L);
        %[p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
        %p = [k1, k2, k3, k4, k5, k6, k12, k13, k14, k15, A, k16,k17];
        %y0 = [L,C,LC,Lesc,Cint, LCint, LR, LCR, LRint,Reng,R];

        [t,y] = ode15s(@odefun,tspan,y0,options,p);%dt from tspan

        [maxval, maxloca] = max((y(:,1) + y(:,3))/LIDcircV*100);%%tumor max ID
        ID_L(k,g) = maxval;
        ID_Ltime(k,g) = t(maxloca)/(60*60);
        Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);%y1y3y7y8,[L,LC, LR, LCR]
        Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);%y4[Lesc]
        Exp_array(k,g) = tumorExposure(t,y,kon_R,koff_R,tspan);
        
    end
end

figure(1)
subplot(211)
% s = mesh(Kd_array,MW_array,Act_array);
contourf(MW_array,cellNum_array,Exp_array)
colormap('jet')
% caxis([0, 1])
% shading flat
% s.FaceColor = 'flat';
% shading interp
colorbar('eastoutside')
set(gca,'XScale','linear')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
xlabel('Molecular Weight (kDa)')
ylabel('Cell Number (个)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Tumor Retention (hr)')

% hold on
% % scatter([15.29E-12, 3.71E-12], [15, 70],100, 'ok','MarkerFaceColor', 'w')
% scatter([649.8E-12, 751.5E-12], [20, 70],100, 'ok','MarkerFaceColor', 'w')
% hold off

subplot(212)
% sl = surfl(Kd_array,MW_array,Sys_Act_array);
contourf(MW_array,cellNum_array,Act_array)
% sl.EdgeColor = 'none';
% shading flat
% shading interp
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','linear')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
xlabel('Molecular Weight (kDa)')
ylabel('Cell Number (个)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Blood Retention (hr)')
% hold on
% scatter([649.8E-12, 751.5E-12], [20, 70],100, 'ok','MarkerFaceColor', 'w')
% hold off

figure(2)

figure(2)
subplot(211)
% s=mesh(Kd_array,MW_array,ID_L);
contourf(MW_array,cellNum_array,ID_L)
colormap('jet')
% shading flat
% s.FaceColor = 'flat';
% shading interp
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','linear')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
xlabel('Molecular Weight (kDa)')
ylabel('Cell Number (个)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Peak Tumor Uptake (%ID)')

subplot(212)
% s=mesh(Kd_array,MW_array,ID_Ltime);
contourf(MW_array,cellNum_array,ID_Ltime)
% shading flat
% s.FaceColor = 'flat';
% shading interp
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','linear')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
xlabel('Molecular Weight (kDa)')
ylabel('Cell Number (个)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('the Time of Peaktumor Uptake (hr)')

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% % % % %Molecular Weight, Affinity and Fractional Activity Arrays
elem =20;                           %number of elements in the array iteration
% MW_array = logspace(0, 3, elem);    %k iterations through for-loop
MW_array = linspace(1,150 , elem);    %k iterations through for-loop
% Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop
cellNum = 2500;
koff_R_array = logspace(-9,-2,elem);
% Rk_array = linspace(20e-6,200e-6,elem);%20um-200um

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
        % Rk = Rk_array(g);
        %         disp(k)
        [p, y0] = Inputs(MW,[],NR,kon_R,koff_R,kendo_R,cellNum,L);
        %[p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
        %p = [k1, k2, k3, k4, k5, k6, k12, k13, k14, k15, A, k16,k17];
        %y0 = [L,C,LC,Lesc,Cint, LCint, LR, LCR, LRint,Reng,R];

        [t,y] = ode15s(@odefun,tspan,y0,options,p);%dt from tspan
        [maxval, maxloca] = max((y(:,1) + y(:,3))/LIDcircV*100);%%tumor max ID
        ID_L(k,g) = maxval;
        ID_Ltime(k,g) = t(maxloca)/(60*60);
        Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);%y1y3y7y8,[L,LC, LR, LCR]
        Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);%y4[Lesc]
        Exp_array(k,g) = tumorExposure(t,y,kon_R,koff_R,tspan);
        
    end
end

figure(1)
subplot(211)
% s = mesh(Kd_array,MW_array,Act_array);
contourf(Kd_array,MW_array,Act_array)
colormap('jet')
% caxis([0, 1])
% shading flat
% s.FaceColor = 'flat';
% shading interp
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Molecular Weight (kDa)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Tumor Retention (hr)')

% hold on
% % scatter([15.29E-12, 3.71E-12], [15, 70],100, 'ok','MarkerFaceColor', 'w')
% scatter([649.8E-12, 751.5E-12], [20, 70],100, 'ok','MarkerFaceColor', 'w')
% hold off

subplot(212)
% sl = surfl(Kd_array,MW_array,Sys_Act_array);
contourf(Kd_array,MW_array,Sys_Act_array)
% sl.EdgeColor = 'none';
% shading flat
% shading interp
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Molecular Weight (kDa)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Blood Retention (hr)')
% hold on
% scatter([649.8E-12, 751.5E-12], [20, 70],100, 'ok','MarkerFaceColor', 'w')
% hold off

figure(2)

figure(2)
subplot(211)
% s=mesh(Kd_array,MW_array,ID_L);
contourf(Kd_array,MW_array,ID_L)
colormap('jet')
% shading flat
% s.FaceColor = 'flat';
% shading interp
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Molecular Weight (kDa)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Peak Tumor Uptake (%ID)')

subplot(212)
% s=mesh(Kd_array,MW_array,ID_Ltime);
contourf(Kd_array,MW_array,ID_Ltime)
% shading flat
% s.FaceColor = 'flat';
% shading interp
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
xlabel('Receptor Affinity (M)')
ylabel('Molecular Weight (kDa)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('the Time of Peaktumor Uptake (hr)')
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Supplementary
% ydot = [dL;dLesc*k5;dLR;dLRint;dR]
% koff_R = 300E-7;
cellNum = 2500;
[p, y0_L] = Inputs(150,[],NR,kon_R,koff_R,kendo_R,cellNum,L);

%Inputs(MW,Th_L,cells,NR,kon_R,koff_R,kendo_R,L)
[t,yL] = ode15s(@odefun,tspan,y0_L,options,p);

% LIDcircV = 5E-8; %injected total[M]
% 
% L_RInterL = yL(:,4)/LIDcircV*100;
ID_L = (yL(:,1) + yL(:,3))/LIDcircV*100;%tumor
L_circL = (yL(:,2))/LIDcircV*100;%blood

figure
subplot(2,1,1),plot(t/(60*60), L_circL,'LineWidth',2)
xlim([0, 24])
subtitle('Blood %ID')
xlabel('Time (hr)')
ylabel('% Injected Dose')
set(gca, 'fontsize', 18)
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
subplot(2,1,2),plot(t/(60*60), ID_L,'LineWidth',2)
xlim([0, 24])
subtitle('Tumor %ID')
xlabel('Time (hr)')
ylabel('% Injected Dose')
set(gca, 'fontsize', 18)
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
%%
%Affinity and Endo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
close all
% % % % %Molecular Weight, Affinity and Fractional Activity Arrays
elem =10;                           %number of elements in the array iteration

% MW_array = logspace(0, 3, elem);    %k iterations through for-loop
MW = 20;
% Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop
kon_R = 1.26E5;                     % [1/M/s] - on rate for IL-2 and IL-2RB from NKTR paper
koff_R = 0.301;                      % [1/s] - off rate for IL-2 and IL-2RB NKTR paper

koff_R_array = logspace(-10,-2,elem);
Kd_array = koff_R_array/kon_R;   %g iterations through for-loop
% [1/s] - endocytic rate iterations through for-loop
kendo_array = linspace(0.01,0.4, elem);   
kendo_array = kendo_array./60;

ID_L = ones(elem);
ID_Ltime = ones(elem);
Act_array = ones(elem);
Sys_Act_array = ones(elem);

for k = 1:elem
    koff_R = koff_R_array(k);
% disp(['KD=',num2str(koff_R/kon_R)])
for g = 1:elem
    kendo_R = kendo_array(g);
    [p, y0] = Inputs(MW,[],cells,NR,kon_R,koff_R,kendo_R,L);
    %[p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
    %p = [k1, k2, k3, k4, k5, k6, k12, k13, k14, k15, A, k16,k17];
    %y0 = [L,C,LC,Lesc,Cint, LCint, LR, LCR, LRint,Reng,R];

    [t,y] = ode15s(@odefun,tspan,y0,options,p);%dt from tspan
    [maxval maxloca] = max((y(:,1) + y(:,3))/LIDcircV*100);%%tumor max ID
    ID_L(k,g) = maxval;
    ID_Ltime(k,g) = t(maxloca)/(60*60);
    Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);%y1y3y7y8,[L,LC, LR, LCR]
    Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);%y4[Lesc]
end
end
figure(1)
subplot(211)
contourf(Kd_array,kendo_array.*60,Act_array)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
ylabel('Internalized Rate (hr^{-1})')
xlabel('Receptor Affinity (M)')
% set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Tumor Retention (hr)')

subplot(212)
contourf(Kd_array,kendo_array.*60,Sys_Act_array)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
ylabel('Internalized Rate (hr^{-1})')
xlabel('Receptor Affinity (M)')
% set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Blood Retention (hr)')
figure(2)
subplot(211)
contourf(Kd_array,kendo_array.*60,ID_L)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
ylabel('Internalized Rate (hr^{-1})')
xlabel('Receptor Affinity (M)')
% set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Peak Tumor Uptake (%ID)')

subplot(212)
contourf(Kd_array,kendo_array.*60,ID_Ltime)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca, 'fontsize', 18)
ylabel('Internalized Rate (hr^{-1})')
xlabel('Receptor Affinity (M)')
% set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('the Time of Peaktumor Uptake (hr)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Cells and NR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
close all
% % % % %Molecular Weight, Affinity and Fractional Activity Arrays
elem =10;                           %number of elements in the array iteration

% MW_array = logspace(0, 3, elem);    %k iterations through for-loop
MW = 20;
% Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop
kon_R = 1.26E11;                     % [1/M/s] - on rate for IL-2 and IL-2RB from NKTR paper
koff_R = 0.301;                      % [1/s] - off rate for IL-2 and IL-2RB NKTR paper
kendo_R = 0.3/(60);  
disp(['KD = ',num2str(koff_R/kon_R)])
%Intrautmoral Receptor+ Cells - [cells/mL]
cells_array = logspace(0, 6, elem);
%Intrautmoral Receptor Density on Cells - [receptors/cell]
NR_array = logspace(2.3, 8, elem);

ID_L = ones(elem);
ID_Ltime = ones(elem);
Act_array = ones(elem);
Sys_Act_array = ones(elem);

for k = 1:elem
    cells = cells_array(k);
% disp(['KD=',num2str(koff_R/kon_R)])
for g = 1:elem
    NR = NR_array(g);
    [p, y0] = Inputs(MW,[],cells,NR,kon_R,koff_R,kendo_R,L);
    %[p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
    %p = [k1, k2, k3, k4, k5, k6, k12, k13, k14, k15, A, k16,k17];
    %y0 = [L,C,LC,Lesc,Cint, LCint, LR, LCR, LRint,Reng,R];

    [t,y] = ode15s(@odefun,tspan,y0,options,p);%dt from tspan
    [maxval maxloca] = max((y(:,1) + y(:,3))/LIDcircV*100);%%tumor max ID
    ID_L(k,g) = maxval;
    ID_Ltime(k,g) = t(maxloca)/(60*60);
    Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);%y1y3y7y8,[L,LC, LR, LCR]
    Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);%y4[Lesc]
end
end
figure(1)
subplot(211)
contourf(cells_array,NR_array,Act_array)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
ylabel('Recpetor Density per Cell')
xlabel('#Cells')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Tumor Retention (hr)')
hold on
scatter(150E3,1000,100, 'ok','MarkerFaceColor', 'w')
hold off

subplot(212)
contourf(cells_array,NR_array,Sys_Act_array)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
ylabel('Recpetor Density per Cell')
xlabel('#Cells')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Blood Retention (hr)')
hold on
scatter(150E3,1000,100, 'ok','MarkerFaceColor', 'w')
hold off

figure(2)
subplot(211)
contourf(cells_array,NR_array,ID_L)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
ylabel('Recpetor Density per Cell')
xlabel('#Cells')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Peak Tumor Uptake (%ID)')
hold on
scatter(150E3,1000,100, 'ok','MarkerFaceColor', 'w')
hold off

subplot(212)
contourf(cells_array,NR_array,ID_Ltime)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
ylabel('Recpetor Density per Cell')
xlabel('#Cells')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('the Time of Peaktumor Uptake (hr)')
hold on
scatter(150E3,1000,100, 'ok','MarkerFaceColor', 'w')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%KD and NR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
close all
% % % % %Molecular Weight, Affinity and Fractional Activity Arrays
elem =10;                           %number of elements in the array iteration

% MW_array = logspace(0, 3, elem);    %k iterations through for-loop
MW = 20;
% Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop
kon_R = 1.26E5;                     % [1/M/s] - on rate for IL-2 and IL-2RB from NKTR paper
koff_R_array = logspace(-10,-2,elem);

Kd_array = koff_R_array/kon_R;   %g iterations through for-loop
kendo_R = 0.3/(60);  
% disp(['KD = ',num2str(koff_R/kon_R)])
%Intrautmoral Receptor Density on Cells - [receptors/cell]
NR_array = logspace(2, 5, elem);

ID_L = ones(elem);
ID_Ltime = ones(elem);
Act_array = ones(elem);
Sys_Act_array = ones(elem);

for k = 1:elem
    koff_R = koff_R_array(k);
% disp(['KD=',num2str(koff_R/kon_R)])
for g = 1:elem
    NR = NR_array(g);
    [p, y0] = Inputs(MW,[],cells,NR,kon_R,koff_R,kendo_R,L);
    %[p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
    %p = [k1, k2, k3, k4, k5, k6, k12, k13, k14, k15, A, k16,k17];
    %y0 = [L,C,LC,Lesc,Cint, LCint, LR, LCR, LRint,Reng,R];

    [t,y] = ode15s(@odefun,tspan,y0,options,p);%dt from tspan
    [maxval maxloca] = max((y(:,1) + y(:,3))/LIDcircV*100);%%tumor max ID
    ID_L(k,g) = maxval;
    ID_Ltime(k,g) = t(maxloca)/(60*60);
    Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);%y1y3y7y8,[L,LC, LR, LCR]
    Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);%y4[Lesc]
end
end
figure(1)
subplot(211)
contourf(Kd_array,NR_array,Act_array)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
ylabel('Recpetor Density per Cell')
xlabel('Receptor Affinity (M)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Tumor Retention (hr)')
% hold on
% scatter(150E3,1000,100, 'ok','MarkerFaceColor', 'w')
% hold off

subplot(212)
contourf(Kd_array,NR_array,Sys_Act_array)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
ylabel('Recpetor Density per Cell')
xlabel('Receptor Affinity (M)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Blood Retention (hr)')
% hold on
% scatter(150E3,1000,100, 'ok','MarkerFaceColor', 'w')
% hold off

figure(2)
subplot(211)
contourf(Kd_array,NR_array,ID_L)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
ylabel('Recpetor Density per Cell')
xlabel('Receptor Affinity (M)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('Peak Tumor Uptake (%ID)')
% hold on
% scatter(150E3,1000,100, 'ok','MarkerFaceColor', 'w')
% hold off

subplot(212)
contourf(Kd_array,NR_array,ID_Ltime)
colormap('jet')
% caxis([0, 1])
colorbar('eastoutside')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca, 'fontsize', 18)
ylabel('Recpetor Density per Cell')
xlabel('Receptor Affinity (M)')
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
title('the Time of Peaktumor Uptake (hr)')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Supplementary
%         ticker = 0;
%         ticker = ticker + 1;
%         A = 0
%         % % % % %Molecular Weight, Affinity and Fractional Activity Arrays
%         elem =10; %40                          %number of elements in the array iteration
%         MW_array = logspace(0, 3, elem);    %k iterations through for-loop
%         Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop
%         Act_array = ones(elem);
%
%         for k = 1:elem
%                 MW = MW_array(k);
%                 for g = 1:elem
%                     Kd = Kd_array(g);
%                     [p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,1,1,0,kon_R,koff_R,0,L,C);
%                     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%                     Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);
%                 end
%         end
%
%         figure
%         contourf(Kd_array,MW_array,Act_array)
%         colormap('jet')
%         %caxis([0,2.5])
%         colorbar('eastoutside')
%         set(gca,'XScale','log')
%         set(gca,'YScale','log')
%         set(gca, 'fontsize', 18)
%         xlabel('Collagen Affinity (M)')
%         ylabel('Molecular Weight (kDa)')
%         set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
%         title('Supplementary Figure 8b')
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Exposure Instead
%         %Counter
%         ticker = 0;
%         ticker = ticker + 1;
%         A = 0;
%         % % % % %Molecular Weight, Affinity and Fractional Activity Arrays
%         elem =100;                           %number of elements in the array iteration
%         MW_array = logspace(0, 3, elem);    %k iterations through for-loop
%         Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop
%         Exp_array = ones(elem);
%         %Sys_Act_array = ones(elem);
%
%         for k = 1:elem
%                 MW = MW_array(k);
%                 for g = 1:elem
%                     Kd = Kd_array(g);
%                     [p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%                     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%                     Exp_array(k,g) = tumorExposure(t,y,kon_R,koff_R,tspan);
%          %           Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);
%                 end
%         end
% %
% %         [Xkd,Ymw] = meshgrid(1:elem,1:elem)
% %         [Xkd2,Ymw2] = meshgrid(1:0.01:elem, 1:0.1:elem);
% %         outData = interp2(Xkd, Ymw, Act_array, Xkd2, Ymw2, 'linear');
% %
% %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Additional Functions Used in this Code
%Fractional Activity in Tumor Function
function cA = tumorA(t1,y1,kon_R, koff_R,tspan)
EC50tumor = koff_R/kon_R;
% L_tumor = y1(:,1) + y1(:,3) + y1(:,7) + y1(:,8);%y1y3y7y8,[L,LC, LR, LCR]
L_tumor = y1(:,1) + y1(:,3);                    %y1y3    ,[L,    LR     ]
Activity_L = 100./(1+EC50tumor./L_tumor);
cA = trapz(t1, Activity_L)/(tspan(end)*100)*tspan(end)/(60*60);%in hours
end

function thalf = tumorExposure(t1,y1,kon_R, koff_R,tspan)
%Half-life
EC50tumor = koff_R/kon_R;
% L_tumor = (y1(:,1) + y1(:,3) + y1(:,7) + y1(:,8))/(1E-6);
L_tumor = y1(:,1) + y1(:,3)/(1E-6);%ID%
f = fit(t1/(60*60),L_tumor,'exp2');
thalf = (log(2)/-f.b)+(log(2)/-f.d);%%sec->hour
end

%Seconds to Days Function
function day = d(t)   %convert s to days
day = t/(24*60*60);
end

%Fractional Activity in Circulation Function
function cA = sysA(t1,y1,kon_R,koff_R,tspan)
EC50sys = koff_R/kon_R;
L_sys = y1(:,2); %from Lesc
Activity_sys = 100./(1+EC50sys./L_sys);
cA = trapz(t1, Activity_sys)/(tspan(end)*100)*tspan(end)/(60*60);%in hours
end





function ydot = odefun(~,y,p)
% Collect param values in cell array and redefine params with names
paramsCell = num2cell(p);
%[k1, k2, k3, k4, k5, k6, k12, k13, k14, k15, A, k16,k17]=paramsCell{:};
[         k3, k4, k5,     k12, k13, k14        , k16, k17]=paramsCell{:};
%k16=Vt
% Collect y-vals in cell array and redefine y-vals with names
yCell = num2cell(y);
% [L,C,LC,Lesc,Cint, LCint, LR, LCR, LRint,Reng, R] = yCell{:};
[L,     Lesc,             LR,        LRint,      R] = yCell{:};

%%dL intumor
% dL = k2*LC - k1*C*L/k5 - k12*L*(R)/k5 + k13*LR + k3*k5*(Lesc-L)-k17*L;
% dL = k3*k5*(Lesc-L)- k12*L*(R)/k5 + k13*LR -k17*L;

dL = k3*k5*(Lesc-L)- k12*L*(R)/k5 + k13*LR;%tracer不考虑k17.也就是自己解体
% dCint = k6*C;

%%dLC
% dLC = k1*C*L/k5 -k2*LC - k6*A*LC - k12*LC*(R)/k5 + k13*LCR;
% dLCint = k6*A*LC;

%%Lblood
%Vtumor=100E-6，Vblood=2E-3；k17=degradation of cytokine in serum measured[1/s]
%k4 = func(Payload serum half-life[days],血浆清除)
% dLesc = k3*k5*(L/k5-Lesc)*(100E-6*k5/2E-3)-k4*Lesc- k17*Lesc;
dLesc = k3*k5*(L/k5-Lesc)*(k16*k5/3)-k4*Lesc;%tracer不考虑k17.也就是自己解体

%dLR = k12*L*(R)/k5  - k13*LR - k1*LR*C/k5 + k2*LCR - k14*(LR);
dLR = k12*L*(R)/k5  - k13*LR - k14*(LR);

%%dLCR
% dLCR = k12*LC*(R)/k5 - k13*LCR - k2*LCR + k1*LR*C/k5;

%%dC
% dC = k2*LC - k1*C*L/k5  + dCint + dLCint - k6*C + k2*LCR - k1*LR*C/k5;

dLRint = k14*(LR);
%dReng= LR+LCR;
%%dR
%dR = dLRint + k16*(dReng+R)*(1-(dReng+R)/(k15*10000)) - k12*L*R/k5 + k13*LR -k12*LC*R/k5 + k13*LCR;
dR  = dLRint - k12*L*R/k5 + k13*LR;

% ydot = [dL; dC; dLC; dLesc*k5; dCint; dLCint; dLR; dLCR; dLRint; dReng; dR].*(1/k5);
ydot = [dL;          dLesc*k5;                dLR;       dLRint;        dR].*(1/k5);
end


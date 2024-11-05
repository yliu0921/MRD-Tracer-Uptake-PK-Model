%supp fig8
clc;clear;close all
%radius vs kclear
% matrix=[
% 0.5	    7.13
% 0.62	4.14
% 0.99	5.58
% 1.19	3.87
% 15.6	2.19
% 2.28	2.35
% 2.74	2.88
% 2.74	2.35
% 2.74	1.49
% 3.47	0.87
% 3.47	0.40
% 3.47	0.32
% 3.47	1.29
% 3.93	0.26
% 3.93	0.20
% 4.37	0.37
% 4.37	0.17
% 4.85	0.02
% 4.85	0.02
% 4.85	0.02
% 4.62	0.36
% 4.62	0.20
% 6.11	0.07
% 14.20	0.03
% 14.20	0.03
% 43.60	0.03];
%[MW,Radius,k_clear]
matrix=[0.5	0.73	7.13
0.62	0.78	4.14
0.99	0.91	5.58
1.19	0.97	3.87
7	1.74	2.48
15.6	2.28	2.19
27	2.74	1.18
27	2.74	2.38
33	2.74	2.85
55	3.47	1.49
55	3.47	0.87
55	3.47	0.4
55	3.47	1.29
80	3.93	0.32
80	3.93	0.26
80	3.93	0.2
110	4.37	0.2
110	4.37	0.37
130	4.62	0.17
130	4.62	0.36
150	4.85	0.2
300	6.11	0.07
150	4.85	0.02
150	4.85	0.02
150	4.85	0.02];

MW = 0.5:0.1:300;
Mol_R = 0.912*(MW).^0.333;

k4 = [];
for i = 1:length(Mol_R)
    k4(i) = cal_kcl(Mol_R(i));
end
figure(1)
subplot(2,1,1)
plot(Mol_R,k4,LineWidth=2);hold on
scatter(matrix(:,2),matrix(:,3),'x',LineWidth=2)
set(gca, 'LineWidth', 2,'FontSize',16); % 将当前轴的线宽设置为2
xlabel('Molecular Radius (nm)')
ylabel('kclear (h-1)')

subplot(2,1,2)
plot(MW,k4,LineWidth=2);hold on
scatter(matrix(:,1),matrix(:,3),'x',LineWidth=2)
set(gca, 'LineWidth', 2,'FontSize',16); % 将当前轴的线宽设置为2
xlabel('Molecular Weight (kDa)')
ylabel('kclear (h-1)')
%%

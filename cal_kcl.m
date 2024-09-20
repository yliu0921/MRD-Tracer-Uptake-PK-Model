% clc;close all;clear
function k4 = cal_kcl(mol_R)
%%
matrix=[
0.73	7.13
0.78	4.14
0.91	5.58
0.97	3.87
1.94	2.19
2.28	2.35
2.74	2.88
2.74	2.35
2.74	1.49
3.47	0.87
3.47	0.40
3.47	0.32
3.47	1.29
3.93	0.26
3.93	0.20
4.37	0.37
4.37	0.17
4.85	0.02
4.85	0.02
4.85	0.02
4.62	0.36
4.62	0.20
6.11	0.07
14.20	0.03
14.20	0.03
43.60	0.03];
% figure(1)
% scatter(data(:,1),data(:,2))
% f = fit(data(:,1),data(:,2),'exp2');
% plot(f,data(:,1),data(:,2))


% elem = 100;
% MW = linspace(1,150,elem);
% for i = 1:elem
%     Mol_R(i) = 0.912*(MW(i))^0.333;
%     CL1(i) = (exp(-3.3+4.9/(1+exp((log(Mol_R(i))-1.4)/0.25))));
%     CL2(i) = f.a*exp(f.b*Mol_R(i)) + f.c*exp(f.d*Mol_R(i));
% %     CL2(i) = 18/(0.3-exp(Mol_R(i)));
% 
% end
% figure(2)
% scatter(data(:,1),data(:,2))
% plot(MW,CL1)
% hold on
% plot(MW,CL2)
% 
% hold on
% plot(MW,f(Mol_R))
% legend('ta','nihe','f')
% 假设matrix是你的矩阵，goal是你想要找到与之最接近的数
 
% 计算每个元素与目标数的差的绝对值
differences = abs(matrix(:,1) - mol_R);
 
% 找到最小的差值
[minVal, minIdx] = min(differences(:));
 
% 找出矩阵中最接近目标数的元素
[row, col] = ind2sub(size(matrix), minIdx);
closestNumber = matrix(row, 2);
 
% 显示结果
% fprintf(['最接近',num2str(mol_R),'的数是', num2str(matrix(row, 1)),',对应的kclear是',num2str(matrix(row, 2)),'1/hr','\n']);


k4 = closestNumber;
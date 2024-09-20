function ydot = odefun_new(~,y,p)
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
dLesc = k3*k5*(L/k5-Lesc)*(100E-6*k5/2E-3)-k4*Lesc;%tracer不考虑k17.也就是自己解体 k16

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

function ydot = ODEs(~,y,p)
% Collect param values in cell array and redefine params with names
paramsCell = num2cell(p);

[k1, k2, k3,k4, k5, k6, k7, k8]=paramsCell{:};
% Collect y-vals in cell array and redefine y-vals with names
yCell = num2cell(y);
[L,Lesc,LR,LRint,R] = yCell{:};

%%dL intumor
dL = k1*k3*(Lesc-L)- k4*L*(R)/k3 + k5*LR;%tracer
%%Lblood
dLesc = k1*k3*(L/k3-Lesc)*(k7*k3/3)-k2*Lesc;%tracer
dLR = k4*L*(R)/k3  - k5*LR - k6*(LR);
dLRint = k6*(LR);
dR  = dLRint - k4*L*R/k3 + k5*LR;
ydot = [dL;dLesc*k3;dLR;dLRint;dR].*(1/k3);
end

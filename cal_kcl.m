function k4 = cal_kcl(mol_R)
%nonlin fit
B0=7.494;
B1=-3.069;
B2=0.3872;
B3=-0.01415;
k4=B0 + B1*mol_R +B2*mol_R^2 + B3*mol_R^3 ;

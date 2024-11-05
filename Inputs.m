% function [p,y0] = Inputs(Kd,MW,Th_L,Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C)
function [p,y0] = Inputs(MW,Th_L,NR,kon_R,koff_R,kendo_R,cellNum,L)
%INTERNAL CALL FUNCTIONS
%   %   SchmidtPerm
%   %   SchmidtVoid

%INPUT PARAMETERS
%   %   Kd          Antigen/Payload Affinity - [M]
%   %   MW          Payload Molecular Weight - [kDa]
%   %   Th_L        Payload serum half-life  - [days]
%   %   cells       Intrautmoral Receptor+ Cells Density- [cells/mL]
%   %   NR          Intrautmoral Receptor Density on Cells - [receptors/cell]

%   %   kon_R       Payload/Receptor on rate - [1/M/s]
%   %   koff_R      Payload/Receptor off rate - [1/s]
%   %   kendo_R     Payload/Receptor endocytic rate - [1/s]
%   %   L           Payload initial concentration - [M]
%   %   R           Antigen initial concentration - [M]

%OUTPUT
%   %   p           p matrix for ODE solver
%   %   y0          y0 matrix for ODE solver 

% Initial concentrations in [M]
Lesc = 15E-8; %blood    
LR = 0;
LRint = 0;


%ASSEMBLING VARIABLES
k4 = kon_R;                                    % [1/M/s] - on rate payload and receptor 
k5 = koff_R;                                   % [1/s] - off rate payload and receptor
k6 = kendo_R;                                  % [1/s] - endocytic rate of payload/receptor
k8 = 0.01/(60*60);                             % [1/s] - degradation of cytokine in serum measured

%Escape Mediated by Diffusion 
Mol_R = 0.912*(MW)^0.333;            %[nm] molecular radius in [nm]
P = SchmidtPerm(Mol_R)*10^-2;        %[m/s] permierbility
e = SchmidtVoid(Mol_R);              %void fraction
Rcap = 8E-6;                         %[m] capillary radius
Rcell = 8E-6;                        %[m] tumor cell radius
Vi = 0.5;
Rk = cal_Rk(cellNum,Rcell,Rcap,Vi);  %[m] kroghs cylinder radius
k1 = 2*P*Rcap/(Rk^2)*(1/e);         %out the tumor
k3 = e;
Vt = pi*Rk^3;%m^3
k7 = Vt*10^3;%L

cells = cellNum/(Vt*10^6);%Ccell 
R = cells*NR*10^3/(6.022E23);   

if isempty(Th_L)
    k2 = cal_kcl(Mol_R)/3600;
else
    k2 = log(2)/(Th_L*24*60*60);                                    % [1/s] kclearance, unit conversion
end

p = [k1, k2, k3,k4, k5, k6,k7, k8];
y0 = [L,Lesc,LR,LRint,R];

end
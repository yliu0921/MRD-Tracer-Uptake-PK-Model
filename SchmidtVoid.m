
%%   References
%   (1) Schmidt and Wittrup, Mol Cancer Ther. 2009
%   (2) Nugent and Jain et al 1983 
%   (3) Paine and Scherr et al 1975
%%  Input and Outputs
%   Input: Mol_R - molecular radius in [nm]
%   Output: P - tumor capillary permeability  [cm/s]
%   Assumptions: Fit parameters (ref 1) for two-pore model of the capillary wall
%%
function e = SchmidtVoid(Mol_R)


%Two-Pore Model of Capillary Wall Fit Parameters (ref 1)
Rsmall = 13.8;      % small capillary radius [nm]
fasmall = 9/10;        % ratio of small pores
Rlarge = 1000;      % large capillary radius [nm]
falarge = 1/10;        % ratio of large pores
Vi = 0.5;           % interstitial fluid volume fraction

%For the Small Capillary
L = Mol_R/Rsmall;   %ratio of molecular radius to pore radius

    %Conditional Dpore/Dfree (DpDf) and partition coefficient (sig)
    if L<0.6 
        sigsmall = (1-L)^2;
    elseif 0.6<=L && L<= 0.9
%         [minD, minI] = min(abs(L-PS(:,1)));
        sigsmall = (1-L)^2;
    elseif 0.9<L
        sigsmall = 0;
    end
    
%For the Large Capillary
L = Mol_R/Rlarge;   %ratio of molecular radius to pore radius

    %Conditional Dpore/Dfree (DpDf) and partition coefficient (sig)
    if L < 0.6 
%         DpDflarge = (1 - 2.105*L + 2.0865*L^3 - 1.7068*L^5 + 0.72603*L^6)/(1-0.75857*L^5);
        siglarge = (1-L)^2;
    elseif 0.6<=L && L<= 0.9
%         [minD, minI] = min(abs(L-PS(:,1)));
%         DpDflarge = (1-L)^2/PS(minI,2);
        siglarge = (1-L)^2;
    elseif 0.9<L
%         DpDflarge = 0;
        siglarge = 0;
    end

%Pore Permeability [cm^2/s]
Plarge = Vi*siglarge*falarge;
Psmall = Vi*sigsmall*fasmall;

%Total Permeability 
e = Psmall + Plarge;

end
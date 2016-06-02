function [rho] = BJ08_corrNew(T1, T2)

% Created by Jack Baker, Oct 19, 2015
% Compute the correlation of epsilons for a number of IMs, using the Baker
% and Jayaram (2008) model. Updated to handle vector inputs.
%
% Documentation is provided in the following document:
% 
% Baker, J. W., and Jayaram, N. (2008). ?Correlation of spectral 
% acceleration values from NGA ground motion models.? Earthquake Spectra, 
% 24(1), 299?317.
%
% INPUT
%
%   T1, T2      = Vectors of the two sets of Sa periods of interest. 
%
% OUTPUT
%
%   rho         = The predicted correlation matrix. If length(T1)=n and
%                   length(T2)=m, then size(rho)= [m n]


for i = 1:length(T1)
    for j = 1:length(T2)
        
        T_min = min(T1(i), T2(j));
        T_max = max(T1(i), T2(j));
        
        C1 = (1-cos(pi/2 - log(T_max/max(T_min, 0.109)) * 0.366 ));
        if T_max < 0.2
            C2 = 1 - 0.105*(1 - 1./(1+exp(100*T_max-5)))*(T_max-T_min)/(T_max-0.0099);
        end
        if T_max < 0.109
            C3 = C2;
        else
            C3 = C1;
        end
        C4 = C1 + 0.5 * (sqrt(C3) - C3) * (1 + cos(pi*(T_min)/(0.109)));
        
        if T_max <= 0.109
            rho(i,j) = C2;
        elseif T_min > 0.109
            rho(i,j) = C1;
        elseif T_max < 0.2
            rho(i,j) = min(C2, C4);
        else
            rho(i,j) = C4;
        end
    end
end


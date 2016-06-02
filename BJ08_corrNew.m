function [rho] = BJ08_corrNew(IM1, IM2)

% Created by Jack Baker, Oct 19, 2015
% Compute the correlation of epsilons for a number of IMs, using the Baker
% and Jayaram (2008) model. Updated to handle vector inputs.
%
% Documentation is provided in the following document:
%
% INPUT
%
%   IM1, IM2      = Vectors of the two IMs of interest. IMs are defined as
%                   follows:
%                     IM1(i) > 0   --> Sa(IM1(i) s )
%                     elseif IM1(i) == -1, then the IM is Ia
%                     ...  
%
% OUTPUT
%
%   rho         = The predicted correlation matrix. If length(IM1)=n and
%                   length(IM2)=m, then size(rho)= [m n]

 
%%%% NOTE: at present, this function just replicates Baker Jayaram (2008)
%%%% and doesn't work for non-SA IMs

for i = 1:length(IM1)
    for j = 1:length(IM2)
        
        T_min = min(IM1(i), IM2(j));
        T_max = max(IM1(i), IM2(j));
        
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


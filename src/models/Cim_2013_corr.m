function [ rho ] = Cim_2013_corr( T1, T2 )
% Horizontal spectral correlation coefficients from the following
% reference
% 
% Cimellaro, G. P. (2013). "Correlation in spectral accelerations for earthquakes in Europe." Earthquake Engineering & Structural Dynamics, 42(4), 623?633.

A0 = -0.0798;
A1 = 0.1147;
A2 = 0.0431; 
A3 = -0.3502;
A4 = 0.0153;

Tmin = min([T1 T2]);
Tmax = max([T1 T2]);


if (Tmin>=0.05 && Tmax <=2.5) 

    numer = A0 + A2*log10(Tmin) + A4*(log10(Tmax))^2;
    denom = 1  + A1*log10(Tmax) + A3*(log10(Tmin))^2;
    
    rho = 1 - (numer/denom) * log(Tmin/Tmax);
    
else % invalid periods
    %     fprintf('Periods outside of allowable range (0.001 to 4s) \n')
    rho = nan;
end

end


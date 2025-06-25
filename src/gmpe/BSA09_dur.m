function [ median, sigma, tau, phi ] = BSA09_dur( Def, M, R_rup, Vs30, Z_tor )
%
% Created by Jack Baker, October 19, 2015
%
% Duration prediction from the following model
%
% Bommer, J. J., Stafford, P. J., and Alarcón, J. E. (2009). ?Empirical 
% Equations for the Prediction of the Significant, Bracketed, and Uniform 
% Duration of Earthquake Ground Motion.? Bulletin of the Seismological 
% Society of America, 99(6), 3217?3233.
%
% INPUTS
%
% Def   = 1 for Ds5-75, =2 for DS5-95
% M     = magnitude
% R_rup = closest distance to rupture (km)
% Vs30  = site Vs30
% Z_tor = depth to top of rupture (km)
%
% OUTPUTS



% check for valid "Def" input
if (Def ~= 1) & (Def ~=2)
    fprintf('Error--invalide value for input paramter ''Def'' \n');
    return
end

% check for other valid inputs
if min([M, R_rup, Vs30, Z_tor]) < 0 % any inputs are '-999'
    median = nan;
    sigma = nan;
    tau = nan;
    phi = nan;
    return
end


% coefficients (from Table 2)
c0 = [-5.6298 -2.2393];
m1 = [1.2619   0.9368];
r1 = [2.0063   1.5686];
r2 = [-0.252  -0.1953];
h1 = [-2.3316  2.5];
v1 = [-0.29   -0.3478];
z1 = [-0.0522 -0.0365];
tauCoeff  = [0.3527 0.3252];
phiCoeff  = [0.4304 0.3460]; % (named 'sigma' in Table 2, but this is intra-event std and so is renamed to phi here) 
sigma_c   = [0.1729 0.1114];
sigma_Tgm = [0.5289 0.4616];

% log-median (equation 5)
i=Def; % rename variable for convenience below
lnDur = c0(i) + m1(i)*M + (r1(i) + r2(i)*M)*log(sqrt(R_rup.^2+ h1(i).^2)) + v1(i)*log(Vs30) + z1(i)*Z_tor;
median = exp(lnDur);

% return proper values of standard deviation terms
sigma = sigma_Tgm(i);
tau = tauCoeff(i);
phi = phiCoeff(i);




end


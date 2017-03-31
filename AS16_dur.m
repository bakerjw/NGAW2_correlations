function [ median, sigma, tau, phi ] = AS16_dur( Def, M, R, Vs30, Mech, Z1, CJ )
%
% Created by Jack Baker, August 9, 2016
% Based on code from Kioumars Afshari
%
% Duration prediction from the following model
%
% Afshari, K., and Stewart, J. P. (2016). "Physically Parameterized
% Prediction Equations for Significant Duration in Active Crustal Regions."
% Earthquake Spectra, Vol. 32, No. 4, pp. 2057-2081.  
% doi: http://dx.doi.org/10.1193/063015EQS106M 
%
% INPUTS
%
% Def   = 1 for Ds5-75,
%       = 2 for DS5-95
%       = 3 for DS20-80
% M     = magnitude
% R     = closest distance to rupture (km)
% Vs30  = site Vs30
% Mech  = Rupture mechanism (0=unknown, 1=Normal, 2=Reverse, 3=Strike-slip)
% Z1    = Depth to shear wave velovity of 1 km/s isosurface (m)
%         Enter -999 if unknown
% CJ    = Region (0 = California,  1 = Japan, -999 = other)
%
% OUTPUTS
%
% median = median predicted duration
% sigma  = log standard deviation of predicted duration
% tau    = within-event log standard deviation
% phi    = between-event log standard deviation


% check "Def" input
if (Def ~= 1) && (Def ~=2) && (Def ~=3)
    fprintf('Error--invalide value for input paramter ''Def'' \n');
    return
end

% estimate median basin depth from Vs30
if CJ==0 % California (eq 11)
    mu_z1=exp(-7.15/4*log((Vs30^4+570.94^4)/(1360^4+570.94^4))-log(1000));
else % other regions (eq 12)
    mu_z1=exp(-5.23/4*log((Vs30^2+412.39^2)/(1360^2+412.39^2))-log(1000));
end

% differential basin depth (eq 10)
dz1 = Z1 - mu_z1;
if Z1 == -999 || CJ == -999
    dz1 = 0;
end

% get coefficients
[M1, M2, b0, b1, b2, b3, Mstar, c1, RR1, RR2, c2, c3, c4, Vref, V1, c5, dz1ref] = get_coeffs(Def, Mech);

% Source term (eq 3)
if M<M1
    F_E = b0; % constant duration at small M
else
    % Stress index parameter (eq 6)
    if M<M2
        deltaSigma=exp(b1+b2*(M-Mstar));
    else
        deltaSigma=exp(b1+b2*(M2-Mstar)+b3*(M-M2));
    end
    
    M_0 = 10^(1.5*M+16.05); % seismic moment (eq 5)
    f_0 = 4.9E6 * 3.2 * (deltaSigma / M_0)^(1/3); % corner frequency (eq 4)
    F_E = 1/f_0;
end


% Path term (eq 7)
if R<RR1
    F_P=c1*R;
elseif R<RR2
    F_P=c1*RR1+c2*(R-RR1);
else
    F_P=c1*RR1+c2*(RR2-RR1)+c3*(R-RR2);
end

% F_deltaz term (eq 9)
if dz1<= dz1ref
    F_deltaz = c5*dz1;
else
    F_deltaz = c5*dz1ref;
end

% Site term (eq 8)
if Vs30<V1
    F_S=c4*log(Vs30/Vref) + F_deltaz;
else
    F_S=c4*log(V1/Vref) + F_deltaz;
end


% median duration (eq 2)
ln_dur = log(F_E + F_P) + F_S;
median = exp(ln_dur);

%  standard deviation terms
[phi, tau] = get_standard_dev(Def, M);
sigma = sqrt(phi^2 + tau^2); % total standard deviation (eq 13)

end


function [M1, M2, b0, b1, b2, b3, Mstar, c1, RR1, RR2, c2, c3, c4, Vref, V1, c5, dz1ref] = get_coeffs(Def, Mech)
% local function to get definition-specific coefficients

% ordering of all coefficients below is [5-75, 5-95, 20-80]

% source coefficients
M1 = [5.35, 5.2, 5.2];
M2 = [7.15, 7.4, 7.4];
b00 = [1.28, 2.182, 0.8822];
b01 = [1.555, 2.541, 1.409];
b02 = [0.7806, 1.612, 0.7729];
b03 = [1.279, 2.302, 0.8804];
b10 = [5.576, 3.628, 6.182];
b11 = [4.992, 3.17, 4.778];
b12 = [7.061, 4.536, 6.579];
b13 = [5.578, 3.467, 6.188];
b2 = [0.9011, 0.9443, 0.7414];
b3 = [-1.684, -3.911, -3.164];
Mstar = [6, 6, 6];

% path coefficients
c1 = [0.1159, 0.3165, 0.0646];
RR1 = [10, 10, 10];
RR2 = [50, 50, 50];
c2 = [0.1065, 0.2539, 0.0865];
c3 = [0.0682, 0.0932, 0.0373];

% site coefficients
c4 = [-0.2246, -0.3183, -0.4237];
Vref = [368.2, 369.9, 369.6];
V1 = [600, 600, 600];
c5 = [0.0006, 0.0006, 0.0005];
dz1ref = [200, 200, 200];

% compute definition-specific coefficients
M1		= M1(Def);
M2		= M2(Def);
b00		= b00(Def);
b01		= b01(Def);
b02		= b02(Def);
b03		= b03(Def);
b10		= b10(Def);
b11		= b11(Def);
b12		= b12(Def);
b13		= b13(Def);
b2		= b2(Def);
b3		= b3(Def);
Mstar	= Mstar(Def);
c1		= c1(Def);
RR1		= RR1(Def);
RR2		= RR2(Def);
c2		= c2(Def);
c3		= c3(Def);
c4		= c4(Def);
Vref	= Vref(Def);
V1		= V1(Def);
c5		= c5(Def);
dz1ref	= dz1ref(Def);

% mechanism-based coefficients
if Mech==0
    b1=b10;
    b0=b00;
elseif Mech==1
    b1=b11;
    b0=b01;
elseif Mech==2
    b1=b12;
    b0=b02;
elseif Mech==3
    b1=b13;
    b0=b03;
end

end

function [phi, tau] = get_standard_dev(Def, M)

% standard deviation coefficients
phi1	= [	0.54	0.43	0.56];
phi2	= [	0.41	0.35	0.45];
tau1	= [	0.28	0.25	0.3];
tau2	= [	0.25	0.19	0.19];

% compute phi (eq 15)
if M<5.5
    phi = phi1(Def);
elseif M<5.75
    phi = phi1(Def) + (phi2(Def) - phi1(Def)) * (M-5.5)/(5.75-5.5);
else
    phi = phi2(Def);
end

% compute tau (eq 14)
if M<6.5
    tau = tau1(Def);
elseif M<7
    tau = tau1(Def) + (tau2(Def) - tau1(Def)) * (M-6.5)/(7-6.5);
else
    tau = tau2(Def);
end

end



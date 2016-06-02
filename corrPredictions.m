% compute reference correlation predictions for the considered IMs
% Jack Baker
% last modified 2 June 2016


clear; close all; clc;

load allIMsResids Periods SaIDX nonSaIDX

rhoPredAll = zeros(109,109);

%% Sa correlations 
% from Baker, J. W., and Jayaram, N. (2008). ?Correlation of spectral acceleration values from NGA ground motion models.? Earthquake Spectra, 24(1), 299?317.

rhoPredAll(1:105,1:105) = BJ08_corrNew(Periods, Periods);


%% d5-75 correlations
imIDX = 106;

% Bradley, B. A. (2011a). ?Correlation of Significant Duration with Amplitude and Cumulative Intensity Measures and Its Use in Ground Motion Selection.? Journal of Earthquake Engineering, 15(6), 809?832.

%D5-75 coefficients
a=[0    -0.45 -0.39 -0.39 -0.06 0.16];
b=[0    -0.39 -0.39 -0.06  0.16  0.0];
e=[0.01  0.09  0.30  1.40   6.5 10.0];
for i=1:length(Periods)
    for j=2:length(a)
        if Periods(i)<=e(j)
            rhoPredAll(i,imIDX)=a(j)   + (b(j)  -a(j))  /log(e(j)  /e(j-1))  *log(Periods(i)/e(j-1));
            break
        end
    end
end

%% d5-95 correlations
imIDX = 107;

% D5-95 coefficients
a95=[0    -0.41 -0.41 -0.38 -0.35 -0.02 0.23];
b95=[0    -0.41 -0.38 -0.35 -0.02  0.23 0.02];
e95=[0.01  0.04  0.08  0.26  1.40   6.0 10.0];
for i=1:length(Periods)
    for j=2:length(a95)
        if Periods(i)<=e95(j)
            %rho_Ds595_SA_fit(i,1)=a95(j) + (b95(j)-a95(j))/log(e95(j)/e95(j-1))*log(Periods(i)/e95(j-1));
            rhoPredAll(i,imIDX)=a95(j) + (b95(j)-a95(j))/log(e95(j)/e95(j-1))*log(Periods(i)/e95(j-1));
            break
        end
    end
end

%% PGA correlations
imIDX = 108;

a_3=[1 0.97];
b_3=[0.895 0.25];
c_3=[0.06 0.8];
d_3=[1.6 0.8];
e_3=[0.2 10];
for i=1:length(Periods)
    for j=1:2
        if Periods(i)<=e_3(j)
            %rho_PGA_SA_fit(i,1)=(a_3(j)+b_3(j))/2-(a_3(j)-b_3(j))/2*tanh(d_3(j)*log((Periods(i)/c_3(j))));
            rhoPredAll(i,imIDX)=(a_3(j)+b_3(j))/2-(a_3(j)-b_3(j))/2*tanh(d_3(j)*log((Periods(i)/c_3(j))));
            break
        end
    end
end


%% PGV correlations
imIDX = 109;

a=[0.73 0.54 0.80 0.76];
b=[0.54 0.81 0.76 0.7];
c=[0.045 0.28 1.1 5];
d=[1.8 1.5 3.0 3.2];
e=[0.1 0.75 2.5 10];
for i=1:length(Periods)
    for j=1:4
        if Periods(i)<=e(j)
            %rho_PGV_SA_fit(i,1)=(a(j)+b(j))/2-(a(j)-b(j))/2*tanh(d(j)*log((Periods(i)/c(j))));
            rhoPredAll(i,imIDX)=(a(j)+b(j))/2-(a(j)-b(j))/2*tanh(d(j)*log((Periods(i)/c(j))));
            break
        end
    end
end

%% add lower-left terms to matrix

rhoPredAll(nonSaIDX,SaIDX) = rhoPredAll(SaIDX,nonSaIDX)';


%% lower right terms

% correlations compiled by JWB from Bradley papers for a number of IMs
%        1      2       3       4       5       6       7       8       9
%        PGA    AI      PGV     ASI     SI      D575    D595    CAV     DSI
rho22 = [1      0.83	0.733	0.928	0.599	-0.442	-0.405	0.7     0.395; ...
        0.83	1       0.73	0.81	0.68	-0.19	-0.2	0.89	0.51; ...
        0.733	0.73	1       0.729	0.89	-0.259	-0.211	0.691	0.8; ...
        0.928	0.81	0.729	1       0.641	-0.411	-0.37	0.703	0.376; ...
        0.599	0.68	0.89	0.641	1       -0.131	-0.079	0.681	0.782; ...
        -0.442	-0.19	-0.259	-0.411	-0.131	1       0.843	0.077	0.074; ...
        -0.405	-0.2	-0.211	-0.37	-0.079	0.843	1       0.122	0.163; ...
        0.7     0.89	0.691	0.703	0.681	0.077	0.122	1       0.565; ...
        0.395	0.51	0.8     0.376	0.782	0.074	0.163	0.565	1];

consideredIDX = [6 7 1 3]; % these are the only IMs used in this work    
rhoPredAll(nonSaIDX,nonSaIDX) = rho22(consideredIDX,consideredIDX);

save rhoPredAll rhoPredAll

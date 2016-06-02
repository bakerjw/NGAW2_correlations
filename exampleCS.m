function [ ] = exampleCS()

% Compute example conditional spectra to illustrate the effect of
% correlation models
% Jack Baker
% 2 June 2016


clear
load allIMsResids 

%% specify an analysis scenario
M = 7; Dist_jb = 20; Vs30 = 500; Fault_Type = 0; z1 = -999; region = 0; % unkonwn fault and z1, global model 

% conditioning Sa
tStar = 1; % conditioning period
tIdx = find(Periods == tStar);
saStar = 0.75; % target Sa at conditioning period


%% load means and covariances
% get correlations  
load rhoDataAll rhoData
rhoData = rhoData(tIdx,1:105); % Sa correlations for just the conditioning period
rhoModel = BJ08_corrNew(tStar, Periods); % model correlations

% predicted medians and standard deviations for Chiou and Youngs
[medianCY, sigmaCY] = CY_2014_nga_mod(M, Periods, Dist_jb, Dist_jb, Dist_jb, 1, 90, 90, 1, Vs30, 0, 1, 1);

% predicted medians and standard deviations for Boore et al.
[medianBSSA, sigmaBSSA] = BSSA_2014_nga_mod(M, Periods, Dist_jb, Fault_Type, region, z1, Vs30);



%% compute conditional spectra

[condMean{1,1}, condSigma{1,1}] = fn_CMS(saStar, tIdx, medianCY, sigmaCY, rhoModel);
[condMean{2,1}, condSigma{2,1}] = fn_CMS(saStar, tIdx, medianCY, sigmaCY, rhoData);

[condMean{1,2}, condSigma{1,2}] = fn_CMS(saStar, tIdx, medianBSSA, sigmaBSSA, rhoModel);
[condMean{2,2}, condSigma{2,2}] = fn_CMS(saStar, tIdx, medianBSSA, sigmaBSSA, rhoData);



%% plot results

i = 1; % CY GMPE
j = 2; % data rho

figure
h1 = plot(Periods, exp(condMean{j,i}), '-k');
hold on
plot(Periods, exp(condMean{j,i} + 2*condSigma{j,i}), '--k')
plot(Periods, exp(condMean{j,i} - 2*condSigma{j,i}), '--k')

i = 2; % BSSA GMPE
h2 = loglog(Periods, exp(condMean{j,i}), '-b');
plot(Periods, exp(condMean{j,i} + 2*condSigma{j,i}), '--b')
plot(Periods, exp(condMean{j,i} - 2*condSigma{j,i}), '--b')
legend([h1 h2], 'Chiou and Youngs, data corr', 'BSSA, data corr')
title('Two GMPEs')
xlabel('Period (s)')
ylabel('Spectral acceleration (g)')
FormatFigure


i = 1; % CY GMPE
j = 1; % model rho

figure
h1 = plot(Periods, exp(condMean{j,i}), '-k');
hold on
plot(Periods, exp(condMean{j,i} + 2*condSigma{j,i}), '--k')
plot(Periods, exp(condMean{j,i} - 2*condSigma{j,i}), '--k')

j = 2; % data rho
h2 = loglog(Periods, exp(condMean{j,i}), '-b');
plot(Periods, exp(condMean{j,i} + 2*condSigma{j,i}), '--b')
plot(Periods, exp(condMean{j,i} - 2*condSigma{j,i}), '--b')
legend([h1 h2], 'Chiou and Youngs, BJ corr', 'Chiou and Youngs, data corr')
title('Two Rho models')
xlabel('Period (s)')
ylabel('Spectral acceleration (g)')
FormatFigure


end

function [condMean, condSigma] = fn_CMS(saStar, tIdx, medianSa, sigmaSa, rho)
% helper function to compute conditional means and standard deviations
    
eps_targ = (log(saStar) - log(medianSa(tIdx)))/sigmaSa(tIdx) % target epsilon
condMean = log(medianSa) + sigmaSa.*eps_targ.*rho;
condSigma = sigmaSa .* sqrt(1 - rho.^2);

   

end


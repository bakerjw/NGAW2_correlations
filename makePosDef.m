% export correlation matrices and ensure positive definiteness
% Jack Baker
% Last modified 2 June 2016

clear


load allIMsResids % load empirical correlation data
load rhoPredAll % load reference model correlations

% screen records and compute correlations
maxDistance = 100; % maximum distance to consider
mSplit = 5; % minimum magnitude to consider
notAllowed = find(magnitude < mSplit | Rjb > maxDistance);
notAllowedEvent = find(eventMagLong < mSplit); 
rhoData = fnGetRho(notAllowed, notAllowedEvent, sigma, tau, phi, residWithin, residBetweenLong );


% compute positive definite matrices

b =ones(109,1); % initialize parameters for function call
tau = 1.0e-5; % initialize parameters for function call
tol = 1.0e-6; % initialize parameters for function call
rhoDataPD = nearestCorrelationMatrix(rhoData, b, tau, tol);
rhoPredPD = nearestCorrelationMatrix(rhoPredAll, b, tau, tol);

% compute minimum eigenvalues
fprintf('Eigenvalues for original data and model: \n')
[min(eig(rhoData)) min(eig(rhoPredAll))]

fprintf('Eigenvalues for updated data and model: \n')
[ min(eig(rhoDataPD)) min(eig(rhoPredPD))]

% quick check of how much the matrices changed
fprintf('Max absolute changes in data and model matrices: \n')
[max(max(abs(rhoData - rhoDataPD))) max(max(abs(rhoPredAll - rhoPredPD)))]

save rhoDataAll rhoData rhoDataPD rhoPredAll rhoPredPD


%% write CSV files of outputs
csvwrite('eSupplement/rhoData.csv',rhoData)
csvwrite('eSupplement/rhoDataPD.csv',rhoDataPD)
csvwrite('eSupplement/rhoPredicted.csv',rhoPredAll)
csvwrite('eSupplement/rhoPredictedPD.csv',rhoPredPD)




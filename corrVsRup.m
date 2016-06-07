% Compute correlations from residuals, plotted versus rupture and site
% parameters
% Jack Baker
% 17 May 2016


clear

% load data
load allIMsResids
load rhoDataAll rhoPredAll % predicted correlations


%% make a simple plot of correlations from two magnitude ranges

% screen records 
maxDistance = 100; % maximum distance to consider
mHi = 6; % minimum magnitude to consider
notAllowed = find(magnitude <= mHi | Rjb > maxDistance);
notAllowedEvent = find(eventMagLong < mHi); 
rhoLargeTotal = fnGetRho(notAllowed, notAllowedEvent, sigma, tau, phi, residWithin, residBetweenLong );

% small magnitude data
mLow = 4.5;
notAllowed = find(magnitude >= mLow | Rjb > maxDistance);
notAllowedEvent = find(eventMagLong >= mLow); 
rhoSmallTotal = fnGetRho(notAllowed, notAllowedEvent, sigma, tau, phi, residWithin, residBetweenLong );

% plot results
tStar = [ 0.1 1];
for i=1:length(tStar)
    tIdx(i) = find(Periods == tStar(i));
end
figure
h1 = semilogx(Periods, rhoLargeTotal(tIdx, SaIDX), '-k', 'linewidth', 2);
hold on
h2 = semilogx(Periods, rhoSmallTotal(tIdx, SaIDX), '-g', 'linewidth', 2);
h3 = plot(Periods, rhoPredAll(tIdx, SaIDX), '--b');
set(gca, 'ylim', [0 1])
set(gca,'xticklabel', [0.01 0.1 1 10])
hx = xlabel('T1 (s)');
hy = ylabel('\rho');
legend([h1(1) h2(1) h3(1)], ['M > ' num2str(mHi) ' NGA-West2 data'], ['M < ' num2str(mLow) ' NGA-West2 data'], 'Baker and Jayaram (2008)', 'location', 'southwest');
FormatFigure
print('-dpdf', ['Figures/sa_correlations_hiLowM.pdf']); % save the figure to a file



%% compute windowed correlations

% screen records 
maxDistance = 100; % maximum distance to consider
minMag = 5;

% windowed Vs30
vs30s = 200:50:750;
dVs30 = 100;
for i = 1:length(vs30s)
    notAllowed = find(soil_Vs30 < vs30s(i)-dVs30 | soil_Vs30 > vs30s(i)+dVs30 | Rjb > maxDistance | magnitude < minMag);
    notAllowedEvent = find(eventMagLong < minMag);
    rhoVs30{i} = fnGetRho(notAllowed, notAllowedEvent, sigma, tau, phi, residWithin, residBetweenLong );
end

% windowed magnitudes
mags = 3:0.25:7.5;
dMag = 0.5;
for i = 1:length(mags)
    notAllowed = find(magnitude < mags(i)-dMag | magnitude > mags(i)+dMag | Rjb > maxDistance);
    notAllowedEvent = find(eventMagLong < mags(i)-dMag | eventMagLong > mags(i)+dMag);
    rhoMag{i} = fnGetRho(notAllowed, notAllowedEvent, sigma, tau, phi, residWithin, residBetweenLong );
end

% windowed distances
notAllowedEvent = find(eventMagLong < minMag); 
dists = 00:5:100;
dDist = 10;
for i = 1:length(dists)
    notAllowed = find(Rjb < dists(i)-dDist | Rjb > dists(i)+dDist | magnitude < minMag);
    rhoDist{i} = fnGetRho(notAllowed, notAllowedEvent, sigma, tau, phi, residWithin, residBetweenLong );
end
save rhoVsRupData vs30s rhoVs30 mags rhoMag dists rhoDist
    


%% pair of conditioning periods
tPairs =      [ 0.3     0.5; ...
                0.1     0.3; ...
                2     0.2; ...
                0.1     3];
         
xLabelText = 'Vs30 (m/s)';
fn_windowed_sa_corr( tPairs, Periods, rhoVs30, rhoPredAll, vs30s, xLabelText )
print('-dpdf', ['Figures/sa_correlations_vsVs30.pdf']); % save the figure to a file


xLabelText = 'Magnitude';
fn_windowed_sa_corr( tPairs, Periods, rhoMag, rhoPredAll, mags, xLabelText )
print('-dpdf', ['Figures/sa_correlations_vsMag.pdf']); % save the figure to a file

xLabelText = 'Distance (km)';
fn_windowed_sa_corr( tPairs, Periods, rhoDist, rhoPredAll, dists, xLabelText )
print('-dpdf', ['Figures/sa_correlations_vsDist.pdf']); % save the figure to a file



%% alternate IM 
imIDX = 106;
tVals =      [ 0.02; 0.1; 1; 5]; % Sa periods to condition against

xLabelText = 'Vs30 (m/s)';
fn_windowed_IM_corr( imIDX, tVals, Periods, rhoVs30, rhoPredAll, vs30s, xLabelText, IMLabel)
print('-dpdf', ['Figures/dur_correlations_vsVs30.pdf']); % save the figure to a file

xLabelText = 'Magnitude';
fn_windowed_IM_corr( imIDX, tVals, Periods, rhoMag, rhoPredAll, mags, xLabelText , IMLabel)
set(gca, 'xlim', [5.5 7.5])
print('-dpdf', ['Figures/dur_correlations_vsMag.pdf']); % save the figure to a file

xLabelText = 'Distance (km)';
fn_windowed_IM_corr( imIDX, tVals, Periods, rhoDist, rhoPredAll, dists, xLabelText, IMLabel )
print('-dpdf', ['Figures/dur_correlations_vsDist.pdf']); % save the figure to a file




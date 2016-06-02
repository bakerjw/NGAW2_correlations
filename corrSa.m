% Compute correlations from residuals
% Jack Baker
% 24 August 2015


clear; close all; clc;


% load data
load allIMsResids
load mixedEffectsResids resid_RotD100Within resid_RotD100BetweenLong
load rhoPredAll % correlation predictions

SaIDX = 1:105; %indices of Sa IMs

% screen records 
maxDistance = 100; % maximum distance to consider
mSplit = 5; % minimum magnitude to consider
notAllowed = find(magnitude < mSplit | Rjb > maxDistance);
notAllowedEvent = find(eventMagLong < mSplit); 
[rhoD50, rhoWI, rhoBTW] = fnGetRho(notAllowed, notAllowedEvent, sigma(SaIDX), tau(SaIDX), phi(SaIDX), residWithin(:,SaIDX), residBetweenLong(:,SaIDX) );
rhoD100 = fnGetRho(notAllowed, notAllowedEvent, sigma(SaIDX), tau(SaIDX), phi(SaIDX), resid_RotD100Within, resid_RotD100BetweenLong );

residWithin(notAllowed,SaIDX) = nan; % flag all not allowed



%% Plots --  compare at a given conditioning period
tStar = [ 0.1 0.3 1 3];
for i=1:length(tStar)
    tIdx(i) = find(Periods == tStar(i));
end
figure
h1 = semilogx(Periods, rhoD50(tIdx, :), '-k');
hold on
h2 = plot(Periods, rhoD100(tIdx, :), '-r');
h3 = plot(Periods, rhoPredAll(tIdx, SaIDX), '--b');
set(gca, 'ylim', [0 1])
hx = xlabel('T1 (s)');
hy = ylabel('\rho');
for i=1:length(tStar)
    text(tStar(i)*.75, 1.05, ['T2=' num2str(tStar(i)) 's'], 'Fontsize', 12);
end
axis([0.01 10 0 1])
legh = legend([h1(1) h2(1) h3(1)], 'Sa_{RotD50} correlation', 'Sa_{RotD100} correlation', 'Baker Jayaram (2008)', 'location', 'southwest');
FormatFigure
print('-dpdf', ['Figures/sa_correlations.pdf']); % save the figure to a file
print('-dpng', ['Figures/sa_correlations.png']); % save the figure to a file

%% compare at a given conditioning period, including within and between event correlations

tStar = [ 0.1 0.3 1 3];
for i=1:length(tStar)
    tIdx(i) = find(Periods == tStar(i));
end
figure
h1 = semilogx(Periods, rhoWI(tIdx, SaIDX), '-k');
hold on
h2 = plot(Periods, rhoBTW(tIdx, SaIDX), '-g');
h3 = plot(Periods, rhoPredAll(tIdx, SaIDX), '--b');
set(gca, 'ylim', [-0.2 1])
hx = xlabel('Period (s)');
hy = ylabel('\rho');
legh = legend([h1(1) h2(1) h3(1)], 'Sa_{RotD50,WI}', 'Sa_{RotD50,BTW}', 'Prediction', 'location', 'southeast');
axis([0.01 10 0 1])
FormatFigure



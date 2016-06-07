% Compute correlations from residuals
% Jack Baker
% 17 May 2016


clear

% load data
load allIMsResids
load rhoPredAll % correlation predictions

numPeriods = length(Periods);
numIMs = length(tau);

% screen records
maxDistance = 100; % maximum distance to consider
minMagnitude = 5; % minimum magnitude to consider
notAllowed = find(magnitude < minMagnitude | Rjb > maxDistance);
notAllowedEventLong = find(eventMagLong < minMagnitude);


% compute correlations
[rhoTotal, rhoWI, rhoBTW] = fnGetRho(notAllowed, notAllowedEventLong, sigma, tau, phi, residWithin, residBetweenLong );
residWithin(notAllowed,:) = nan; % flag all not allowed
residBetweenLong(notAllowedEventLong,:) = nan; % flag all not allowed



%% observations versus predictions

imIdx = 109;

figure
h1 = semilogx(Periods, rhoTotal(imIdx, SaIDX), 'linewidth', 2);
hold on
h2 = semilogx(Periods, rhoWI(imIdx, SaIDX), 'linewidth', 1);
h3 = semilogx(Periods, rhoBTW(imIdx, SaIDX), 'linewidth', 1);
h4 = semilogx(Periods, rhoPredAll(imIdx, SaIDX), '--k');
plot([0.01 10], [0 0], '-k')

set(gca, 'ylim', [-0.5 1.5])
hx = xlabel('Period (s)');
hy = ylabel('\rho');
legend(IMLabel{imIdx}, 'Within-event', 'Between-event', 'Bradley prediction', 'location', 'northwest');
FormatFigure

%% figures for paper

imIdx = 106:107;

figure
h1 = semilogx(Periods, rhoTotal(imIdx(1), SaIDX), '-b', 'linewidth', 2);
hold on
h2 = semilogx(Periods, rhoPredAll(imIdx(1), SaIDX), '--b');
h3 = semilogx(Periods, rhoTotal(imIdx(2), SaIDX), '-r', 'linewidth', 2);
h4 = semilogx(Periods, rhoPredAll(imIdx(2), SaIDX), '--r');
plot([0.01 10], [0 0], '-k')

set(gca, 'ylim', [-0.5 0.5])
set(gca, 'ytick', [-0.5:0.1:0.5])
set(gca,'xticklabel', [0.01 0.1 1 10])
hx = xlabel('T (s)');
hy = ylabel('\rho');
legend(IMLabel{imIdx(1)}, 'Bradley (2011a)', IMLabel{imIdx(2)}, 'Bradley (2011a)', 'location', 'northwest');
FormatFigure
print('-dpdf', ['Figures/dur_correlations.pdf']); % save the figure to a file



imIdx = 108:109;

figure
h1 = semilogx(Periods, rhoTotal(imIdx(1), SaIDX), '-b', 'linewidth', 2);
hold on
h2 = semilogx(Periods, rhoPredAll(imIdx(1), SaIDX), '--b');
h3 = semilogx(Periods, rhoTotal(imIdx(2), SaIDX), '-r', 'linewidth', 2);
h4 = semilogx(Periods, rhoPredAll(imIdx(2), SaIDX), '--r');
plot([0.01 10], [0 0], '-k')

set(gca, 'ylim', [-0 1])
set(gca,'xticklabel', [0.01 0.1 1 10])
hx = xlabel('T (s)');
hy = ylabel('\rho');
legend(IMLabel{imIdx(1)}, 'Bradley (2011b)', IMLabel{imIdx(2)}, 'Bradley (2012)', 'location', 'southwest');
FormatFigure
print('-dpdf', ['Figures/pga_correlations.pdf']); % save the figure to a file



%% scatter plots
indices = [95 imIdx];


% don't consider notAllowed records (by setting others to nan)
residWithin(notAllowed,:)       = nan; % flag all not allowed
residTotal(notAllowed,:)       = nan; % flag all not allowed
residBetweenLong(notAllowedEventLong,:) = nan; % flag all not allowed
nEvents = sum(~isnan(sum(residBetweenLong(:,indices),2)));
nRecs = sum(~isnan(sum(residWithin(:,indices),2)));



SIGMA = corr(residBetweenLong(:,indices), 'rows', 'pairwise');
figure
plot(residBetweenLong(:,indices(1)), residBetweenLong(:,indices(2)), 'o')
xlabel(['Between event residual, ' IMLabel{indices(1)}])
ylabel(['Between event residual, ' IMLabel{indices(2)}])
axis([-4 4 -4 4])
title(['\rho = ' num2str(SIGMA(1,2),2), ', ' num2str(nEvents) ' observations']);
FormatFigure

SIGMA = corr(residWithin(:,indices), 'rows', 'pairwise');
figure
plot(residWithin(:,indices(1)), residWithin(:,indices(2)), 'o')
xlabel(['Within event residual, ' IMLabel{indices(1)}])
ylabel(['Within event residual, ' IMLabel{indices(2)}])
axis([-4 4 -4 4])
title(['\rho = ' num2str(SIGMA(1,2),2), ', ' num2str(nRecs) ' observations']);
FormatFigure



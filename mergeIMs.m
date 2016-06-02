% Build matrices of all IM values
% Jack Baker
% 2/18/2016, updated 4/11/2016

% master IMs spec:
% residXX(:,1:105) -   SaRotD50 (BSSA 14)
% IMLabel{106} = 'Ds575';
% IMLabel{107} = 'Ds595';
% IMLabel{108} = 'PGA';
% IMLabel{109} = 'PGV';

clear
load mixedEffectsResids

numIMs = 109;
numRecs = length(magnitude);
numEvents = max(eqid);
SaIDX = 1:105;
nonSaIDX = 106:109;


%% remove event terms when there is only one ground motion
resid_RotD50BetweenLong(resid_BetweenNumRecs <= 1) = nan;
residOtherBetweenLong(residOtherBetweenNumRecs <= 1) = nan;


% make big matrices of correlations
residWithin = nan*ones(numRecs, numIMs); % initialize
residWithin(:,1:105)   = resid_RotD50Within;
residWithin(:,106:109) = resid_OtherWithin;

residBetweenLong = nan*ones(numEvents, numIMs);  % initialize
residBetweenLong(:,1:105)   = resid_RotD50BetweenLong;
residBetweenLong(:,106:109) = residOtherBetweenLong;

residTotal = nan*ones(numRecs, numIMs); % initialize
residTotal(:,1:105)   = resid_RotD50Total;
residTotal(:,106:109) = resid_OtherTotal;


% store label information
IMLabel = [];
for i = 1:length(Periods)
    IMLabel{i} = ['Sa(' num2str(Periods(i)) 's)'];
end
IMLabel{106} = 'Ds575';
IMLabel{107} = 'Ds595';
IMLabel{108} = 'PGA';
IMLabel{109} = 'PGV';


% get phi, tau and sigma values for a "typical" event 
M = 7; Dist_jb = 20; Vs30 = 500; Fault_Type = 0; z1 = -999; region = 0; % unkonwn fault and z1, global model 

[~, ~, phi, tau] = CY_2014_nga_mod(M, Periods, Dist_jb, Dist_jb, Dist_jb, 1, 90, 90, 1, Vs30, 0, 1, 1);
[~, ~, tau(106), phi(106) ] = BSA09_dur( 1, M, Dist_jb, Vs30, 5 );
[~, ~, tau(107), phi(107) ] = BSA09_dur( 2, M, Dist_jb, Vs30, 5 );
[~, ~, phi(108), tau(108)] = CY_2014_nga_mod(M, 0, Dist_jb, Dist_jb, Dist_jb, 1, 90, 90, 1, Vs30, 0, 1, 1);
[~, ~, phi(109), tau(109)] = CY_2014_nga_mod(M, -1, Dist_jb, Dist_jb, Dist_jb, 1, 90, 90, 1, Vs30, 0, 1, 1);

% recompute sigmas for internal consistency
sigma = sqrt(tau.^2 + phi.^2);


%% save everything
save allIMsResids residWithin residBetweenLong  residTotal sigma phi tau eventMagLong IMLabel Periods magnitude soil_Vs30 eqid Rjb SaIDX nonSaIDX

%% some extra plots

% histograms of residuals
imIdx = 109;

figure
hist(residWithin(:,imIdx),100)
xlabel(['Within-event residuals, ' IMLabel{imIdx} ])
ylabel('Number of observations')
FormatFigure

figure
hist(residBetweenLong(:,imIdx))
xlabel(['Between-event residuals, ' IMLabel{imIdx} ])
ylabel('Number of observations')
FormatFigure

% event terms vs magnitude
figure
plot(eventMagLong, residBetweenLong(:,95), '.b')
hold on
plot([3 8], [0 0], '-k')
xlabel('Magnitude')
ylabel([IMLabel{imIdx} ' within-event residuals'])
legend('RotD50 Within-event')
set(gca, 'ylim', [-5 5])
FormatFigure





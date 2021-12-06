% Compute within- and between-event residuals
% Jack Baker
% Last updated 2 June 2016
% Updated 6 December 2021 to utilize corrected Z1_CVMH values

clear

%% load raw data
load NGA_W2_corr_meta_data % needed data from the NGA-West2 flatfile
load Z1_Z25_updated % load corrected values of the Z1 variables
Z1_CVMH(Z1_CVMH>0) = Z1_CVMH(Z1_CVMH>0)/1000; % the flatfile values are in units of m, the GMM below expects units of km (leave the -999 values unchanged)

load durObs % 5-75% (column 1) and 5-95% (column 2) durations for the first 8163 ground motions in the database


%% initialize matrices

% trim long periods so periods match the GMPE range
tIndex = find(Periods<=10); 
Periods = Periods(tIndex);
Sa_RotD50 = Sa_RotD50(:,tIndex);
Sa_RotD100 = Sa_RotD100(:,tIndex);

% get sizes of matrices
numRecs = length(eqid);
numT = length(Periods);

[eventEQID,eventIdx] = unique(eqid); % get indices of unique events
eventEQID(1) = []; % throw out the first case, which is eqid = -999
eventIdx(1) = []; % throw out the first case, which is eqid = -999
for i = 1:length(eventEQID)
   idx = find(eqid == eventEQID(i),1);
   eventMag(i,1) = magnitude(idx); % get event magnitudes 
end

% initialize residual matrices
resid_RotD100Within = nan*ones(numRecs, numT);
resid_RotD50Within  = nan*ones(numRecs, numT);
resid_RotD50Total   = nan*ones(numRecs, numT);
resid_OtherWithin      = nan*ones(numRecs, 4);
resid_OtherTotal      = nan*ones(numRecs, 4);

resid_RotD100BetweenLong = nan*ones(max(eqid), numT); % leave a row for each eqid
resid_RotD50BetweenLong  = nan*ones(max(eqid), numT); % leave a row for each eqid
residOtherBetweenLong      = nan*ones(max(eqid), 4);
residOtherBetweenNumRecs   = nan*ones(max(eqid), 4);

%% get predictions

% get fw/hw terms for CY 2014 GMPE
FwHw = zeros(numRecs,1);
for i = 1:numRecs
    if strcmp(Fhw{i},'hw')
        FwHw(i) = 1;
    end
end

% Predicted RotD100/RotD50 ratios
[ rotD100Factor, sigmaRotD100Factor , phiRotD100Factor, tauRotD100Factor] = SB_2014_ratios( Periods );

% get max usable period for each ground motion
maxUsableT = 1./lowest_usable_freq;

%% Sa residuals
for i = 1:numT % period index
    i
    allowableFilter = (maxUsableT > Periods(i));
    idxTemp = find(allowableFilter & chiouYoungsUsed); % find records used by Chiou and Youngs, and with usable data at this period
    [~, eqIdx] = sort(eqid(idxTemp)); % now adjust index so that it puts the data in order by eqid
    idx = idxTemp(eqIdx);

    numUsableRecs(i) = length(idx);
    
    %%%%%%%%%%%%%%%%%%%%%% RotD50 %%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 1:length(idx)
        [median_GMPE(k,1), sigma_GMPE(k,1), phi_GMPE(k,1), tau_GMPE(k,1)] = CY_2014_nga_mod(magnitude(idx(k)), Periods(i), closest_D(idx(k)), ...
                                                    Rjb(idx(k)), Rx(idx(k)), Z_tor(idx(k)), dip(idx(k)), rakeAngle(idx(k)), ...
                                                    Z1_CVMH(idx(k)), soil_Vs30(idx(k)), FwHw(idx(k)), 1, Region_BSSA(idx(k)));
    end
    
    % build data structure
    input.imObservations = Sa_RotD50(idx,i);
    input.imMedian   = median_GMPE;
    input.totalSigma = sigma_GMPE;
    input.interSigma = tau_GMPE;
    input.intraSigma = phi_GMPE;
    input.eventIds   = eqid(idx);
    Out=GetInterIntraEventResiduals(input);

    
    % store results
    resid_RotD50Total(idx,i)  = (log(Sa_RotD50(idx,i)) - log(median_GMPE))./sigma_GMPE;
    resid_RotD50Within(idx,i)  = Out.intraEventResidualsNormalized;

    idxNgms = find(Out.eventData.eventNumGms <= 1);
    Out.eventData.interEventResidualNormalized(idxNgms) = nan;

    resid_RotD50BetweenLong(Out.eventData.eventId ,i) = Out.eventData.interEventResidualNormalized; % alternate set of event residuals (for debugging)
    resid_BetweenNumRecs(Out.eventData.eventId ,i) = Out.eventData.eventNumGms;

    %%%%%%%%%%%%%%%%%%%%%% RotD100 %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % build data structure
    input.imObservations = Sa_RotD100(idx,i);
    input.imMedian   = median_GMPE .* rotD100Factor(i);
    input.totalSigma = sqrt(sigma_GMPE.^2 +  sigmaRotD100Factor(i).^2 );
    input.interSigma = sqrt(phi_GMPE.^2 +    phiRotD100Factor(i).^2   );
    input.intraSigma = sqrt(tau_GMPE.^2 +    tauRotD100Factor(i).^2   );
    input.eventIds   = eqid(idx);
    Out=GetInterIntraEventResiduals(input);

    % store results
    resid_RotD100Within(idx,i)  = Out.intraEventResidualsNormalized;
    idxNgms = find(Out.eventData.eventNumGms <= 1);
    Out.eventData.interEventResidualNormalized(idxNgms) = nan;

    resid_RotD100BetweenLong(Out.eventData.eventId ,i) = Out.eventData.interEventResidualNormalized; % alternate set of event residuals (for debugging)

    clear median_GMPE sigma_GMPE phi_GMPE tau_GMPE % clean up variables before moving to the next period
end



%% Duration residuals

allowableFilter = (maxUsableT > 1); % put some limit on filter frequencies
idxTemp = find(NGA_num <= 8163 & allowableFilter & chiouYoungsUsed);
[~, eqIdx] = sort(eqid(idxTemp)); % now adjust index so that it puts the data in order by eqid
idx = idxTemp(eqIdx);
numUsableRecsDur = length(idx);

% get mechanism codes for AS16
MechAS16 = zeros(size(mechanism));
MechAS16(mechanism==-999) = 0; % unknown
MechAS16(mechanism==1) = 1; % Normal
MechAS16(mechanism==2) = 2; % Reverse
MechAS16(mechanism==0) = 3; % Strike-slip


for i = 1:2 % two duration definitions
    for k = 1:length(idx)
        %[med_DurBSA09(k,1), sigmaDur(k,1), tauDur(k,1), phiDur(k,1) ] = BSA09_dur(i, magnitude(idx(k)), closest_D(idx(k)), soil_Vs30(idx(k)), Z_tor(idx(k)) );
        [med_Dur(k,1), sigmaDur(k,1), tauDur(k,1), phiDur(k,1) ] = AS16_dur( i, magnitude(idx(k)), closest_D(idx(k)), soil_Vs30(idx(k)), MechAS16(idx(k)), Z1_CVMH(idx(k)), Region_BSSA(idx(k)) );
    end
    
    % build data structure
    input.imObservations = durObs(idx,i);
    input.imMedian   = med_Dur;
    input.totalSigma = sigmaDur;
    input.interSigma = tauDur;
    input.intraSigma = phiDur;
    input.eventIds   = eqid(idx);
    Out=GetInterIntraEventResiduals(input);

    resid_OtherTotal(idx,i)  = (log(durObs(idx,i)) - log(med_Dur))./sigmaDur;
    resid_OtherWithin(idx,i)  = Out.intraEventResidualsNormalized;
    idxNgms = find(Out.eventData.eventNumGms <= 1);
    Out.eventData.interEventResidualNormalized(idxNgms) = nan;
    residOtherBetweenLong(Out.eventData.eventId ,i) = Out.eventData.interEventResidualNormalized;
    residOtherBetweenNumRecs(Out.eventData.eventId ,i) = Out.eventData.eventNumGms;

    clear med_Dur sigmaDur tauDur phiDur % clean up variables
end

%% PGA residuals

i=3; % column index in nonSa data matrix
Period = 0;
allowableFilter = (maxUsableT > 0.1); % put some limit on filter frequencies

idxTemp = find(allowableFilter & chiouYoungsUsed);
[~, eqIdx] = sort(eqid(idxTemp)); % now adjust index so that it puts the data in order by eqid
idx = idxTemp(eqIdx);
numUsableRecsPGA = length(idx);

for k = 1:length(idx)
    [median_GMPE(k,1), sigma_GMPE(k,1), phi_GMPE(k,1), tau_GMPE(k,1)] = CY_2014_nga_mod(magnitude(idx(k)), Period, closest_D(idx(k)), ...
        Rjb(idx(k)), Rx(idx(k)), Z_tor(idx(k)), dip(idx(k)), rakeAngle(idx(k)), ...
        Z1_CVMH(idx(k)), soil_Vs30(idx(k)), FwHw(idx(k)), 1, Region_BSSA(idx(k)));
end

% build data structure
input.imObservations = PGA_RotD50(idx);
input.imMedian   = median_GMPE;
input.totalSigma = sigma_GMPE;
input.interSigma = tau_GMPE;
input.intraSigma = phi_GMPE;
input.eventIds   = eqid(idx);
Out=GetInterIntraEventResiduals(input);

resid_OtherTotal(idx,i)  = (log(PGA_RotD50(idx)) - log(median_GMPE))./sigma_GMPE;
resid_OtherWithin(idx,i)  = Out.intraEventResidualsNormalized;
idxNgms = find(Out.eventData.eventNumGms <= 1);
Out.eventData.interEventResidualNormalized(idxNgms) = nan;
residOtherBetweenLong(Out.eventData.eventId ,i) = Out.eventData.interEventResidualNormalized;
residOtherBetweenNumRecs(Out.eventData.eventId ,i) = Out.eventData.eventNumGms;

clear median_GMPE sigma_GMPE phi_GMPE tau_GMPE % clean up variables

%% PGV residuals

i=4; % column index in nonSa data matrix
Period = -1;
allowableFilter = (maxUsableT > 1); % put some limit on filter frequencies

idxTemp = find(allowableFilter & chiouYoungsUsed);
[~, eqIdx] = sort(eqid(idxTemp)); % now adjust index so that it puts the data in order by eqid
idx = idxTemp(eqIdx);
numUsableRecsPGA = length(idx);

for k = 1:length(idx)
    [median_GMPE(k,1), sigma_GMPE(k,1), phi_GMPE(k,1), tau_GMPE(k,1)] = CY_2014_nga_mod(magnitude(idx(k)), Period, closest_D(idx(k)), ...
        Rjb(idx(k)), Rx(idx(k)), Z_tor(idx(k)), dip(idx(k)), rakeAngle(idx(k)), ...
        Z1_CVMH(idx(k)), soil_Vs30(idx(k)), FwHw(idx(k)), 1, Region_BSSA(idx(k)));
end

% build data structure
input.imObservations = PGV_RotD50(idx);
input.imMedian   = median_GMPE;
input.totalSigma = sigma_GMPE;
input.interSigma = tau_GMPE;
input.intraSigma = phi_GMPE;
input.eventIds   = eqid(idx);
Out=GetInterIntraEventResiduals(input);

resid_OtherTotal(idx,i)  = (log(PGV_RotD50(idx)) - log(median_GMPE))./sigma_GMPE;
resid_OtherWithin(idx,i)  = Out.intraEventResidualsNormalized;
idxNgms = find(Out.eventData.eventNumGms <= 1);
Out.eventData.interEventResidualNormalized(idxNgms) = nan;
residOtherBetweenLong(Out.eventData.eventId ,i) = Out.eventData.interEventResidualNormalized;
residOtherBetweenNumRecs(Out.eventData.eventId ,i) = Out.eventData.eventNumGms;

clear median_GMPE sigma_GMPE phi_GMPE tau_GMPE % clean up variables



%% get event magnitudes for the between event residuals
eventMagLong = nan*ones(max(eqid),1);
for i = 1:length(eventEQID)
    eventMagLong(eventEQID(i)) = eventMag(i);
end


%% save data

save('mixedEffectsResids.mat', 'resid_RotD50Within',  'resid_RotD100Within',  ...
                               'resid_RotD50BetweenLong', 'resid_RotD100BetweenLong', ...
                               'resid_OtherWithin', 'residOtherBetweenLong', ...
                               'resid_BetweenNumRecs', 'residOtherBetweenNumRecs', ...
                               'resid_OtherTotal', 'resid_RotD50Total', ...
                               'Periods', 'eventMagLong', 'eqid', 'magnitude', 'soil_Vs30', 'Rjb', 'numUsableRecsDur', 'numUsableRecs')

%% plots
% number of records per period
figure
semilogx([Periods ], [numUsableRecs ], '-b')
hold on
semilogx([ 11], [ numUsableRecsDur], 'ob')
axis([0.01 11 0 15000])
xlabel('Period (s)')
ylabel('Number of ground motions')
FormatFigure




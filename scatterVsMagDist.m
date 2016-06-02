% Compute correlations from residuals
% Jack Baker
% 17 May 2016


clear; close all; clc

% load data
load allIMsResids

% screen records 
maxDistance = 150; % maximum distance to consider
imIndices = [95 106];
% imIndices = [68 51 ];
minM = 5;

% mRange = [5 6]; % minimum magnitude to consider
% figTitle = [num2str(mRange(1)) ' < M < ' num2str(mRange(2))];
% notAllowed = find(magnitude < mRange(1) | magnitude > mRange(2) | Rjb > maxDistance);
% notAllowedEvent = find(eventMagLong < mRange(1) | eventMagLong > mRange(2)); 
% fnScatterPlot(notAllowed, notAllowedEvent, residWithin, residBetweenLong, imIndices, IMLabel, figTitle );


% % see what the underlying data looks like in our figure
% mags = 3:0.25:7.5;
% dMag = 0.5;
% for i = 1:2:length(mags)
%     notAllowed = find(magnitude < mags(i)-dMag | magnitude > mags(i)+dMag | Rjb > maxDistance);
%     notAllowedEvent = find(eventMagLong < mags(i)-dMag | eventMagLong > mags(i)+dMag); 
%     figTitle = [num2str(mags(i)-dMag) ' < M < ' num2str(mags(i)+dMag)];
%     fnScatterPlot(notAllowed, notAllowedEvent, residWithin, residBetweenLong, imIndices, IMLabel, figTitle );
% end



%% windowed Vs30 
vs30s = 200:100:900;
dVs30 = 100;
for i = 1:length(vs30s)
    notAllowed = find(soil_Vs30 < vs30s(i)-dVs30 | soil_Vs30 > vs30s(i)+dVs30 | Rjb > maxDistance | magnitude < minM);
    notAllowedEvent = find(eventMagLong < minM); 
    
    fnScatterPlot(notAllowed, notAllowedEvent, residWithin, residBetweenLong, imIndices, IMLabel, num2str(vs30s(i)) );

    rhoVs30{i} = fnGetRho(notAllowed, notAllowedEvent, sigma, tau, phi, residWithin, residBetweenLong );
end


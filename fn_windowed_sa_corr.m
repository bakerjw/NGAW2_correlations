function [ output_args ] = fn_windowed_sa_corr( tPairs, Periods, rhoDataAll, rhoPrediction, paramVals, xLabelText )
% Helper function to create a plot of Sa correlations, versus some rupture
% or site parameter that is being varied
%
% Jack Baker
% June 2, 2016

nT = length(tPairs);
for i=1:nT
    for j = 1:2
        imIDX(i,j) = find(Periods == tPairs(i,j));
    end
end

% plot versus soil_Vs30
for i=1:length(rhoDataAll)
    for k = 1:nT
        rhoData{k}(i) = rhoDataAll{i}(imIDX(k,1),imIDX(k,2));
    end
end

figure
plot(paramVals, rhoData{1}, '-k', 'linewidth', 2)
hold on
plot(paramVals, rhoData{2}, '-b', 'linewidth', 2)
plot(paramVals, rhoData{3}, '-g', 'linewidth', 2)
plot(paramVals, rhoData{4}, '-r', 'linewidth', 2)

plot(paramVals, ones(size(paramVals))*rhoPrediction(imIDX(1,1),imIDX(1,2)), '--k', 'linewidth', 1)
plot(paramVals, ones(size(paramVals))*rhoPrediction(imIDX(2,1),imIDX(2,2)), '--b', 'linewidth', 1)
plot(paramVals, ones(size(paramVals))*rhoPrediction(imIDX(3,1),imIDX(3,2)), '--g', 'linewidth', 1)
plot(paramVals, ones(size(paramVals))*rhoPrediction(imIDX(4,1),imIDX(4,2)), '--r', 'linewidth', 1)

axis([min(paramVals) max(paramVals) 0 1])
set(gca, 'ytick', 0:0.2:1)
legend([num2str(tPairs(1,1),2) 's, ' num2str(tPairs(1,2),2) 's'], ...
       [num2str(tPairs(2,1),2) 's, ' num2str(tPairs(2,2),2) 's'], ...
       [num2str(tPairs(3,1),2) 's, ' num2str(tPairs(3,2),2) 's'], ...
       [num2str(tPairs(4,1),2) 's, ' num2str(tPairs(4,2),2) 's'] ...
   );
hx = xlabel(xLabelText);
hy = ylabel('\rho');
FormatFigure


end


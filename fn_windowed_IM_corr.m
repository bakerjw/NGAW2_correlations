function [ output_args ] = fn_windowed_IM_corr( imIDX, tVals, Periods, rhoDataAll, rhoPrediction, paramVals, xLabelText, IMLabel )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here




nT = length(tVals);
for i=1:nT
    tIDX(i) = find(Periods == tVals(i));
end

% plot versus predictor parameter
for i=1:length(rhoDataAll)
    for k = 1:nT
        rhoData{k}(i) = rhoDataAll{i}(imIDX,tIDX(k));
    end
end

figure
plot(paramVals, rhoData{1}, '-k', 'linewidth', 2)
hold on
plot(paramVals, rhoData{2}, '-b', 'linewidth', 2)
plot(paramVals, rhoData{3}, '-g', 'linewidth', 2)
plot(paramVals, rhoData{4}, '-r', 'linewidth', 2)

plot(paramVals, ones(size(paramVals))*rhoPrediction(imIDX,tIDX(1)), '--k', 'linewidth', 1)
plot(paramVals, ones(size(paramVals))*rhoPrediction(imIDX,tIDX(2)), '--b', 'linewidth', 1)
plot(paramVals, ones(size(paramVals))*rhoPrediction(imIDX,tIDX(3)), '--g', 'linewidth', 1)
plot(paramVals, ones(size(paramVals))*rhoPrediction(imIDX,tIDX(4)), '--r', 'linewidth', 1)

axis([min(paramVals) max(paramVals) -0.7 0.7])
% set(gca, 'ytick', 0:0.2:1)
legend([num2str(tVals(1),2) 's vs ' IMLabel{imIDX}], ...
       [num2str(tVals(2),2) 's vs ' IMLabel{imIDX}], ...
       [num2str(tVals(3),2) 's vs ' IMLabel{imIDX}], ...
       [num2str(tVals(4),2) 's vs ' IMLabel{imIDX}] ...
   );
hx = xlabel(xLabelText);
hy = ylabel('\rho');
FormatFigure


end


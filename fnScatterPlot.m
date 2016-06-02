function [  ] = fnScatterPlot(notAllowed, notAllowedEvent, residWithin, residBetween, imIndices, IMLabel, figTitle )

% scatter plots of mixed-effects residuals


% remove un-allowed indices
residWithin(notAllowed,:)       = nan; % flag all not allowed
residBetween(notAllowedEvent,:) = nan; % flag all not allowed

% scatter plot of within-event residuals
figure
plot(residWithin(:,imIndices(1)), residWithin(:,imIndices(2)), '.b')
axis square
axis([-4 4 -4 4])
xlabel([ IMLabel(imIndices(1)) ' within-event residual'])
ylabel([ IMLabel(imIndices(2)) ' within-event residual'])
title(figTitle)
FormatFigure

% scatter plot of between-event residuals
% figure
% plot(residBetween(:,imIndices(1)), residBetween(:,imIndices(2)), '.b')
% axis square
% axis([-4 4 -4 4])
% xlabel([ IMLabel(imIndices(1)) ' between-event residual'])
% ylabel([ IMLabel(imIndices(2)) ' between-event residual'])
% title(figTitle)
% FormatFigure


end


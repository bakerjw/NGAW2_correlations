function [ rho, rhoWithin, rhoBetween ] = fnGetRho(notAllowed, notAllowedEvent, sigma, tau, phi, residWithin, residBetween )

% Get correlation matrix from mixed-effects residuals


% remove un-allowed indices
residWithin(notAllowed,:)       = nan; % flag all not allowed
residBetween(notAllowedEvent,:) = nan; % flag all not allowed

% correlations
rhoWithin =  corr(residWithin,  'rows', 'pairwise');
rhoBetween = corr(residBetween, 'rows', 'pairwise');


% aggregate
numIMs = length(sigma);
for i = 1:numIMs
    for j = 1:numIMs
        rho(i,j) = tau(i)/sigma(i)*tau(j)/sigma(j) * rhoBetween(i,j) + ...
                   phi(i)/sigma(i)*phi(j)/sigma(j) * rhoWithin(i,j);
    end
end



end


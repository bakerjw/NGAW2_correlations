% Plots of correlation data and models from other sources
% Jack Baker
% last modified 2 June 2016

clear; close all; clc;

load ../rhoDataAll rhoData % load NGA-W2 observations
load ../allIMsResids Periods SaIDX
rhoData = rhoData(SaIDX, SaIDX); % keep only spectral correlations

%% compute targets

rho_BJ = BJ08_corrNew( Periods, Periods);
rho_GA = ga_09_correlation( Periods, Periods);

for i =1:length(Periods)
    for j = 1:length(Periods)
        rho_ASA(i,j) = ASA_2014_corr( Periods(i), Periods(j));
        rho_Cim(i,j) = Cim_2013_corr( Periods(i), Periods(j));
        rho_AlAtik(i,j) = AlAtik_corr( Periods(i), Periods(j));
    end
end

%% plots for a conditioning T2

tStar = [ 0.1 0.3 1 3]; % periods to plot
legendtext{1} = 'NGA-West2 data'; % legend data
legendtext{2} = 'Al Atik (2011)';
legendtext{3} = 'Akkar et al. (2014)';
legendtext{4} = 'Baker and Jayaram (2008)';
legendtext{5} = 'Cimellaro (2013)';
legendtext{6} = 'Goda and Atkinson (2009)';

figure

for i=1:length(tStar)
    tIdx = find(Periods == tStar(i));
    
    figure
    h1 = semilogx(Periods, rhoData(tIdx, :), '-b', 'linewidth', 2);
    hold on
    plot(Periods, rho_AlAtik(tIdx, :));
    plot(Periods, rho_ASA(tIdx, :), '--');
    plot(Periods, rho_BJ(tIdx, :), ':k', 'linewidth', 1);
    plot(Periods, rho_Cim(tIdx, :));
    plot(Periods, rho_GA(tIdx, :), '-.');
    set(gca, 'ylim', [0 1])
    axis([0.01 10 0 1])
    set(gca, 'xtick', 10.^[-2:1])
    if(i>2) % add x tick labels
        set(gca,'xticklabel', [0.01 0.1 1 10])
        xlabel('T1 (s)');
    else
        set(gca,'xticklabel', [])
    end
    
    if mod(i,2) % add y axis labels
        ylabel('\rho');
    else
        set(gca,'yticklabel', [])
    end
   
    if i==1
        legend(legendtext, 'location', 'southwest');
    end
    FormatFigure
    print('-dpdf', ['../Figures/saCorrCompare' num2str(i) '.pdf']); % save the figure to a file



end



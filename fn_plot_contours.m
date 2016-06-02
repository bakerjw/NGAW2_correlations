function plot_contours(R, T)

% Created by Jack Baker, 2/28/07
% basic function to create a contour plot of correlation coefficients
% versus the period pairs of interest

% R(n x n) = correlation matrix
% T (1 x n) = vector of period values



% v = [0:0.1:0.9]; % user-defined contour levels
v = [0.1:0.1:0.8 0.89]; % I used this because there is a funny jump between the 0.89 level 0.9 
[C,h] = contour(T,T,R,v, 'k');
hold on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
% set(gca, 'XDir', 'reverse')
% set(gca, 'YDir', 'reverse')
set(gca,'XTick',[0.01;0.1;1;10])
set(gca,'YTick',[0.01;0.1;1;10])
hx = xlabel('T_1');
set(hx, 'fontsize', 18)
hy = ylabel('T_2');
set(hy, 'fontsize', 18)
set(gca,'zlim', [0 1])
set(gca, 'fontsize', 18)
% h = clabel(C,h,'manual'); %uncomment this line if you want to add labels to the contours
% set(h, 'fontsize', 18)
set(gca,'PlotBoxAspectRatioMode', 'manual')

% clabel(C,h);

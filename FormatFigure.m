%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file sets the defaults for the figures
% 6-6-05
% CBH of Stanford University, haselton@stanford.edu
% modified by Jack Baker, most recently 12/10/13
%
% This assumes that the figure is open and you have 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% use fixed figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3.25 3.25]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 3.25 3.25]);


% Alter the plot
label_size = 10;
axis_size = 10;

% axis labels
set(gca, 'FontSize', axis_size);

% set labels' font size
axLabels = get(gca,{'XLabel', 'YLabel', 'ZLabel'});
set([axLabels{:}], 'FontSize', label_size);

% set legend font size
legH = findobj(gca, 'tag', 'legend');
if ~isempty(legH)
    set(legH, 'FontSize', label_size);
end

% set title font size
titleH = get(gca, 'title');
if ~isempty(titleH)
    if iscell(titleH)
        for i = 1:length(titleH)
            set(titleH{i}, 'FontSize', label_size);
        end
    else
        set(titleH, 'FontSize', label_size);
    end
end







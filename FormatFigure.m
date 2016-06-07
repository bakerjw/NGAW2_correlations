%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reformat figures for better display
% last modified by Jack Baker, October 6, 2015
%
% This assumes that the figure is open 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% target font sizes
label_size = 16;
axis_size = 14;

% use fixed figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 4.5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.5]);


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




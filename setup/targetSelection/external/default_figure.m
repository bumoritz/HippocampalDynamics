function HF = default_figure(pos_inches)
% lateral start, start from below, width, height
% 20, 1.5 is bottom of new screen
% stolen and modified from https://github.com/MouseLand/stringer-pachitariu-et-al-2018b/blob/master/utils/default_figure.m

if nargin < 1
    pos_inches = [25,3.5,7,5];
end

set(0, 'DefaultFigureColor', 'White', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 10, ...
    'DefaultAxesFontSizeMode', 'manual', ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultAxesFontWeight', 'normal', ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesTickDir', 'Out', ...
    'DefaultAxesTickDirMode', 'manual',...
    'DefaultAxesBox', 'Off', ...
    'DefaultLineColor', 'Black', ...
    'DefaultLineLineWidth', 1, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 15, ...
    'DefaultTextFontSizeMode', 'manual', ...
    'DefaultTextVerticalAlignment', 'top', ...
    'DefaultTextFontWeight', 'normal', ...
    'DefaultTextUnits', 'normalized',...
    'DefaultTextFontName', 'Arial', ...
    'DefaultLegendFontSize', 8,...
    'DefaultLegendFontSizeMode','manual',...
    'DefaultFigureUnits','inches',...
    'DefaultFigurePosition',[5,3.5,7,5],...
    'DefaultFigurePaperPositionMode','auto',...
'DefaultFigurePaperSize',[11 11]);

HF=figure('Position', pos_inches);
pos = get(HF,'Position');
set(HF,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+.15, pos(4)+.4])
HF.Renderer = 'Painters';
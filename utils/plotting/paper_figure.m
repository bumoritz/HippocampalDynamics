function HF = paper_figure(pos_inches)
% lateral start, start from below, width, height in inches
% 20, 1.5 is bottom of new screen
% taken and modified from https://github.com/MouseLand/stringer-pachitariu-et-al-2018b/blob/master/utils/default_figure.m

if nargin < 1
    pos_inches = [0,0.5,20,9.9]; % ~full second screen: [20,0.5,20,9.9] % smaller: [25,3.5,7,5]
end

set(0, 'DefaultFigureColor', 'White', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 6, ...
    'DefaultAxesFontSizeMode', 'manual', ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultAxesFontWeight', 'normal', ...
    'DefaultAxesTitleFontSize', 6, ...
    'DefaultAxesTitleFontWeight', 'normal', ...
    'DefaultAxesLabelFontSize', 1,...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesTickDir', 'Out', ...
    'DefaultAxesTickDirMode', 'manual',...
    'DefaultAxesBox', 'Off', ...
    'DefaultLineColor', 'Black', ...
    'DefaultLineLineWidth', 1, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 6, ...
    'DefaultTextFontSizeMode', 'manual', ...
    'DefaultTextVerticalAlignment', 'top', ...
    'DefaultTextFontWeight', 'normal', ...
    'DefaultTextUnits', 'normalized',...
    'DefaultTextFontName', 'Arial', ...
    'DefaultLegendFontSize', 6,...
    'DefaultLegendFontSizeMode','manual',...
    'DefaultFigureUnits','inches',...
    'DefaultFigurePosition',[5,3.5,7,5],...
    'DefaultFigurePaperPositionMode','auto',...
    'DefaultFigurePaperSize',[11 11],...
    'DefaultFigureRendererMode','manual');

HF = figure('Position', pos_inches);
pos = get(HF,'Position');
set(HF,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
HF.Renderer = 'Painters';



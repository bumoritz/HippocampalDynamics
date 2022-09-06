function F = myHistogram(data,p,plt)

%% Set defaults

if ~exist('plt','var')
    plt = struct();
end
if ~isfield(plt,'illustrator')
    plt.illustrator = false;
end
if ~isfield(plt,'figure')
    plt.figure = false;
end
if ~isfield(plt,'norm')
    plt.norm = 'probability';
end
if ~isfield(plt,'plotEdges')
    plt.edgeColor = 'none';
elseif plt.plotEdges
    plt.edgeColor = p.col.black;
end
if ~isfield(plt,'ylabel')
    plt.ylabel = 'Proportion';
end
if ~isfield(plt,'xlim_method')
    if isfield(plt,'xlim')
        plt.xlim_method = 'constant';
    else
        plt.xlim_method = 'minmax';
    end
end
if ~isfield(plt,'ylim_method')
    if isfield(plt,'ylim')
        plt.ylim_method = 'constant';
    else
        plt.ylim_method = 'none';
    end
end
if ~isfield(plt,'color')
    plt.color = p.col.gray;
end


%% Preparations

if plt.figure
    F = default_figure;
else
    F = [];
end

if ~isfield(plt,'xlim')
    if strcmp(plt.xlim_method,'minmax')
        plt.xlim = [nanmin(data),nanmax(data)];
    elseif strcmp(plt.xlim_method,'zeromax')
        plt.xlim = [0,nanmax(data)];
    end
end


%% Core

hold on
if isfield(plt,'binWidth')
    h=histogram(data,'BinWidth',plt.binWidth,'Normalization',plt.norm,'EdgeColor',plt.edgeColor,'FaceColor',plt.color);
elseif isfield(plt,'edges')
    h=histogram(data,plt.edges,'Normalization',plt.norm,'EdgeColor',plt.edgeColor,'FaceColor',plt.color);
elseif isfield(plt,'numBins')
    h=histogram(data,plt.numBins,'Normalization',plt.norm,'EdgeColor',plt.edgeColor,'FaceColor',plt.color);
else
    h=histogram(data,'Normalization',plt.norm,'EdgeColor',plt.edgeColor,'FaceColor',plt.color);
end
if isfield(plt,'criterion')
    xline(plt.criterion,'Color',p.col.photostim,'LineStyle','-','LineWidth',3);
end
if isfield(plt,'dataPoint') && (~isnan(plt.dataPoint))
    xline(plt.dataPoint,'Color',p.col.black,'LineStyle','-','LineWidth',3);
end
hold off

if ~isfield(plt,'ylim')
    if strcmp(plt.ylim_method,'minmax')
        plt.ylim = [nanmin(h.Values),nanmax(h.Values)];
    elseif strcmp(plt.ylim_method,'zeromax')
        plt.ylim = [0,nanmax(h.Values)];
    elseif strcmp(plt.ylim_method,'none')
        plt.ylim = [];
    end
end

xlim(plt.xlim)
if ~isempty(plt.ylim)
    ylim(plt.ylim)
end
if isfield(plt,'xticks')
    xticks(plt.xticks)
end
if isfield(plt,'yticks')
    yticks(plt.yticks)
end
if ~plt.illustrator
    xlabel(plt.xlabel)
    ylabel(plt.ylabel)
end

end
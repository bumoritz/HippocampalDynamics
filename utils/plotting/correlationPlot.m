function F = correlationPlot(square_x,square_y,p,plt)

%% Set defaults

if ~exist('plt','var')
    plt = struct();
end
if ~isfield(plt,'illustrator')
    plt.illustrator = false;
end
if ~isfield(plt,'figure')
    plt.figure = true;
end
if ~isfield(plt,'type')
    plt.type = 'subplots'; % subplots or combined or subplot1 or subplot2
end
if ~isfield(plt,'fit')
    plt.fit = 'linear';
end
if ~isfield(plt,'ylabel')
    plt.ydatatype = 'perc';
elseif any([strcmp(plt.ylabel,'Performance'),strcmp(plt.ylabel,'Hit rate'),strcmp(plt.ylabel,'Correct rejection rate')])
    plt.ydatatype = 'perc';
else
    plt.ydatatype = 'normal';
end
if ~isfield(plt,'ylabel')
    plt.ylabel = 'Performance';
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
        if strcmp(plt.ydatatype,'perc')
            plt.ylim_method = 'constant';
            plt.ylim = [0,100];
        else
            plt.ylim_method = 'minmax';
        end
    end
end


%% Core

if plt.figure
    F = default_figure;
else
    F = [];
end

if strcmp(plt.type,'subplots')
    subplot(1,2,1)
    correlationPlot_core(square_x,square_y,p,plt,1);
    subplot(1,2,2)
    correlationPlot_core(square_x,square_y,p,plt,2);
elseif strcmp(plt.type,'subplot1')
    correlationPlot_core(square_x,square_y,p,plt,1);
elseif strcmp(plt.type,'subplot2')
    correlationPlot_core(square_x,square_y,p,plt,2);
elseif strcmp(plt.type,'combined')
    correlationPlot_core(square_x,square_y,p,plt,0);
end


end
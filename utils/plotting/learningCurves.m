function F = learningCurves(dat,plt,p,dat_col)

%% Set defaults

if ~exist('plt','var')
    plt = struct();
end
if ~isfield(plt,'illustrator')
    plt.illustrator = false;
end
if ~isfield(plt,'type')
    plt.type = 'subplots';
end
if ~isfield(plt,'traces')
    plt.traces = 'mean+-sem';
end
if ~isfield(plt,'col')
    plt.col = {};
    for k=1:size(dat,1)
        for j=1:size(dat,2)
            plt.col{k,j} = [0,0,0];
        end
    end
end
if ~isfield(plt,'ylim')
    plt.ylim = [0,100];
end
if ~isfield(plt,'xlabel')
    plt.xlabel = 'Trial block';
end
if ~isfield(plt,'ylabel')
    plt.ylabel = 'Performance';
end

if exist('dat_col','var')
    plt.colormap = 'copper';
    plt.colormap_flip = true;
    plt.col_nan = [1,1,1]*2/3;
    plt.clim = [0,120];
    plt.cbar_numTicks = 7;
    plt.clabel = 'Average running speed (cm/s)';
    if size(dat,1)>1
        error('Only one condition (k) can be plotted at the same time when using dat_col.')
    end
    if ~strcmp(plt.traces,'individuals')
        error('plt.traces has to be set to individuals when using dat_col.')
    end
end

%% Pre-processing

for k=1:size(dat,1)
    for j=1:size(dat,2)
        dat{k,j} = smoothdata(dat{k,j},2,'gaussian',p.lcs.smoothing_sd*5);
    end
end


%% Core

F = default_figure([0,3.5,20,5]);

if strcmp(plt.type,'subplots')

    for j=1:size(dat,2)
        subplot(1,size(dat,2),j)
        yline(50,'Color',[0,0,0],'LineStyle',':');
        hold on
            
        for k=1:size(dat,1)            
            if strcmp(plt.traces,'mean+-sem')
                shadedErrorBar(1:size(dat{k,j},2),nanmean(dat{k,j},1),nansem(dat{k,j},1),'lineProps',plt.col{k,j});             
            elseif strcmp(plt.traces,'individuals')
                if exist('dat_col','var')
                    for i=1:size(dat{1,j},1)
                        if ~isnan(dat_col{1,j}(i))
                            if plt.colormap_flip
                                temp = flipud(eval(plt.colormap));
                            else
                                temp = eval(plt.colormap);
                            end
                            [~,temp2] = discretize(plt.clim,size(temp,1)); % returns NaN leading to error when e.g. speed is below 0
                            plot(dat{k,j}(i,:)','Color',temp(discretize(dat_col{1,j}(i),temp2),:))
                        else
                            plot(dat{k,j}(i,:)','Color',plt.col_nan)
                        end
                    end
                else
                    for i=1:size(dat{k,j},1)
                        plot(dat{k,j}(i,:)','Color',plt.col{k,j})
                    end
                end
            end
        end
            
        xlim([0,size(dat{k,j},2)])
        ylim(plt.ylim)
        ytickformat('percentage')
        if ~plt.illustrator
            xlabel(plt.xlabel)
            ylabel(plt.ylabel)
            title(['Day ',num2str(j)])
        end
    end
	if ~plt.illustrator && exist('dat_col','var')
        colormap(plt.colormap); 
        if plt.colormap_flip
            temp=colorbar; set(gca,'CLim',plt.clim); temp.Ticks = [plt.clim(1):plt.clim(2)/(plt.cbar_numTicks-1):plt.clim(2)]; temp.TickLabels = flip(temp.TickLabels); temp.Direction='reverse';
        else
            temp=colorbar; set(gca,'CLim',plt.clim); temp.Ticks = [plt.clim(1):plt.clim(2)/(plt.cbar_numTicks-1):plt.clim(2)];
        end
        temp.Label.String = plt.clabel;
    end
end


end
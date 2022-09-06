function correlationPlot_core(square_x,square_y,p,plt,subplot_number)

if ~isfield(plt,'xlim')
    if strcmp(plt.xlim_method,'minmax')
        plt.xlim = [nanmin([square_x.seq_switch;square_x.seq_postswitch;square_x.ctrl_switch;square_x.ctrl_postswitch]),nanmax([square_x.seq_switch;square_x.seq_postswitch;square_x.ctrl_switch;square_x.ctrl_postswitch])];
    elseif strcmp(plt.xlim_method,'zeromax')
        plt.xlim = [0,nanmax([square_x.seq_switch;square_x.seq_postswitch;square_x.ctrl_switch;square_x.ctrl_postswitch])];
    end
end
if ~isfield(plt,'ylim')
    if strcmp(plt.ylim_method,'minmax')
        plt.ylim = [nanmin([square_y.seq_switch;square_y.seq_postswitch;square_y.ctrl_switch;square_y.ctrl_postswitch]),nanmax([square_y.seq_switch;square_y.seq_postswitch;square_y.ctrl_switch;square_y.ctrl_postswitch])];
    elseif strcmp(plt.ylim_method,'zeromax')
        plt.ylim = [0,nanmax([square_y.seq_switch;square_y.seq_postswitch;square_y.ctrl_switch;square_y.ctrl_postswitch])];
    end
    if strcmp(plt.ydatatype,'perc')
        plt.ylim = plt.ylim*100;
    end
end

if subplot_number==0
    data_x_seq = rmmissing([square_x.seq_switch;square_x.seq_postswitch]);
    data_x_ctrl = rmmissing([square_x.ctrl_switch;square_x.ctrl_postswitch]);
    data_y_seq = rmmissing([square_y.seq_switch;square_y.seq_postswitch]);
    data_y_ctrl = rmmissing([square_y.ctrl_switch;square_y.ctrl_postswitch]);
    data_x_seq_switch = rmmissing(square_x.seq_switch);
    data_x_ctrl_switch = rmmissing(square_x.ctrl_switch);
    data_y_seq_switch = rmmissing(square_y.seq_switch);
    data_y_ctrl_switch = rmmissing(square_y.ctrl_switch);
    data_x_seq_postswitch = rmmissing(square_x.seq_postswitch);
    data_x_ctrl_postswitch = rmmissing(square_x.ctrl_postswitch);    
    data_y_seq_postswitch = rmmissing(square_y.seq_postswitch);
    data_y_ctrl_postswitch = rmmissing(square_y.ctrl_postswitch);
elseif subplot_number==1
    data_x_seq = rmmissing(square_x.seq_switch);
    data_x_ctrl = rmmissing(square_x.ctrl_switch);
    data_y_seq = rmmissing(square_y.seq_switch);
    data_y_ctrl = rmmissing(square_y.ctrl_switch);
elseif subplot_number==2
    data_x_seq = rmmissing(square_x.seq_postswitch);
    data_x_ctrl = rmmissing(square_x.ctrl_postswitch);    
    data_y_seq = rmmissing(square_y.seq_postswitch);
    data_y_ctrl = rmmissing(square_y.ctrl_postswitch);
end

if strcmp(plt.ylabel,'Performance')
    yline(50,'Color',p.col.black,'LineStyle',':');
end
if strcmp(plt.ydatatype,'perc')
    data_y_seq = data_y_seq*100;
    data_y_ctrl = data_y_ctrl*100;
    if subplot_number==0
        data_y_seq_switch = data_y_seq_switch*100;
        data_y_ctrl_switch = data_y_ctrl_switch*100;
        data_y_seq_postswitch = data_y_seq_postswitch*100;
        data_y_ctrl_postswitch = data_y_ctrl_postswitch*100;
    end
end
hold on

[corr_r_seq,corr_p_seq] = fitLine(data_x_seq,data_y_seq,p.col.seq);
[corr_r_ctrl,corr_p_ctrl] = fitLine(data_x_ctrl,data_y_ctrl,p.col.ctrl);

if subplot_number==0
    scatter(data_x_seq_switch,data_y_seq_switch,'d','MarkerEdgeColor',p.col.black,'MarkerFaceColor',p.col.seq)
    scatter(data_x_ctrl_switch,data_y_ctrl_switch,'d','MarkerEdgeColor',p.col.black,'MarkerFaceColor',p.col.ctrl)
    scatter(data_x_seq_postswitch,data_y_seq_postswitch,'s','MarkerEdgeColor',p.col.black,'MarkerFaceColor',p.col.seq)
    scatter(data_x_ctrl_postswitch,data_y_ctrl_postswitch,'s','MarkerEdgeColor',p.col.black,'MarkerFaceColor',p.col.ctrl)
elseif subplot_number==1
    scatter(data_x_seq,data_y_seq,'d','MarkerEdgeColor',p.col.black,'MarkerFaceColor',p.col.seq)
    scatter(data_x_ctrl,data_y_ctrl,'d','MarkerEdgeColor',p.col.black,'MarkerFaceColor',p.col.ctrl)
elseif subplot_number==2
    scatter(data_x_seq,data_y_seq,'s','MarkerEdgeColor',p.col.black,'MarkerFaceColor',p.col.seq)
    scatter(data_x_ctrl,data_y_ctrl,'s','MarkerEdgeColor',p.col.black,'MarkerFaceColor',p.col.ctrl)
end

xlim(plt.xlim)
ylim(plt.ylim)
if isfield(plt,'xticks')
    xticks(plt.xticks)
end
if isfield(plt,'yticks')
    yticks(plt.yticks)
end
if strcmp(plt.ydatatype,'perc')
    ytickformat('percentage')
end
if ~plt.illustrator
    xlabel(plt.xlabel)
    ylabel(plt.ylabel)
    if subplot_number==1
        title('Switch days')
    elseif subplot_number==2
        title('Post-switch days')
    end
    
    %
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['seq: n = ',num2str(length(data_x_seq)),', \rho = ',num2str(corr_r_seq,2),', p = ',num2str(corr_p_seq,2),newline,...
        'ctrl: n = ',num2str(length(data_x_ctrl)),', \rho = ',num2str(corr_r_ctrl,2),', p = ',num2str(corr_p_ctrl,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
end
        
end

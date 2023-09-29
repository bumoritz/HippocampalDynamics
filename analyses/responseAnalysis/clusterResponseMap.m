function F = clusterResponseMap(s2p_meta,iscell,trg,resp,clusterID,issubplot,firstCluster,plt,paper,exclusionZone,info)

if nargin < 10
    exclusionZone = false;
end
if nargin < 9
    paper = false;
end
if nargin < 7
    firstCluster = false;
end

c_abs = 2;%3;
c_label = '\Delta_{resp-base} z-scored \DeltaF/F'; %'\Delta_{resp-base} \DeltaF/F';
titleText = ['Cluster ',num2str(clusterID)];

meds = zeros(length(s2p_meta.stat),2);
for i = 1:length(s2p_meta.stat)
    meds(i,:) = [double(s2p_meta.stat{i}.med(2)),double(s2p_meta.stat{i}.med(1))];
end

if paper
	F = paper_figure([0,0.5,mm2inch(3*34),mm2inch(3*34)]); hold on;
elseif ~issubplot
    F = default_figure([-20,0.5,20,9.9]);
else
    F = NaN;
end

if paper
    for n=1:length(iscell)
        if iscell(n)==1
            if ~isnan(resp.avgAct_net(n,clusterID))
                x = s2p_meta.stat{n}.xpix';
                y = s2p_meta.stat{n}.ypix';
                temp = boundary(double(x),double(y));
                temp0 = eval('redblue');
                [~,temp2] = discretize([-c_abs,c_abs],size(temp0,1)); % returns NaN leading to error when e.g. speed is below 0
                if resp.avgAct_net(n,clusterID) > c_abs
                    temp3 = c_abs;
                elseif resp.avgAct_net(n,clusterID) < -c_abs
                    temp3 = -c_abs;
                else
                    temp3 = resp.avgAct_net(n,clusterID);
                end
                this_col = temp0(discretize(temp3,temp2),:);
                h=patch(x(temp),y(temp),this_col,'FaceAlpha',1,'EdgeColor',[0.8,0.8,0.8]);
            end
        end
    end
else
    if firstCluster
        img = zeros(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
        for n=1:length(iscell)
            if iscell(n)==1
                for i=1:length(s2p_meta.stat{n}.xpix)
                    if ~isnan(resp.avgAct_net(n,clusterID))
                        img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = resp.firstCluster.avgAct_net(n,clusterID);
                    end
                end
            end
        end
    else
        img = zeros(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
        for n=1:length(iscell)
            if iscell(n)==1
                for i=1:length(s2p_meta.stat{n}.xpix)
                    if ~isnan(resp.avgAct_net(n,clusterID))
                        img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = resp.avgAct_net(n,clusterID);
                    end
                end
            end
        end
    end
end

if ~paper
    imagesc(img);
    set(gca,'CLim',[-c_abs,c_abs]);
    colormap(redblue);
    temp=colorbar; 
    temp.Label.String = c_label;
end
daspect([1,1,1]);
hold on

if issubplot
%    scatter(meds(find(resp.stats_sig_pos(:,clusterID)==1),1),meds(find(resp.stats_sig_pos(:,clusterID)==1),2),'g*','LineWidth',0.5,'SizeData',25)
else
    if paper
        
    elseif firstCluster
        scatter(meds(find(resp.firstCluster.stats_sig_neg(:,clusterID)==1),1),meds(find(resp.firstCluster.stats_sig_neg(:,clusterID)==1),2),'y.','LineWidth',1,'SizeData',50)
        scatter(meds(find(resp.firstCluster.stats_sig_pos(:,clusterID)==1),1),meds(find(resp.firstCluster.stats_sig_pos(:,clusterID)==1),2),'g.','LineWidth',1,'SizeData',50)
        scatter(meds(find(resp.firstCluster.responders_main.respAll(:,clusterID)==1),1),meds(find(resp.firstCluster.responders_main.respAll(:,clusterID)==1),2),'g*','LineWidth',1,'SizeData',70)
    else
        scatter(meds(find(resp.stats_sig_neg(:,clusterID)==1),1),meds(find(resp.stats_sig_neg(:,clusterID)==1),2),'y.','LineWidth',1,'SizeData',50)
        scatter(meds(find(resp.stats_sig_pos(:,clusterID)==1),1),meds(find(resp.stats_sig_pos(:,clusterID)==1),2),'g.','LineWidth',1,'SizeData',50)
        scatter(meds(find(resp.responders_main.respAll(:,clusterID)==1),1),meds(find(resp.responders_main.respAll(:,clusterID)==1),2),'g*','LineWidth',1,'SizeData',70)
    end
end
temp = resp.responders_main.targetedCells(:,clusterID);
temp = temp(~isnan(temp));
if paper
	%scatter(meds(temp,1),meds(temp,2),'k+','LineWidth',0.5,'SizeData',5)
    scatter(trg.laser_y(:,clusterID),trg.laser_x(:,clusterID),'k+','LineWidth',0.5,'SizeData',5)
elseif issubplot
    scatter(meds(temp,1),meds(temp,2),'k+','LineWidth',1,'SizeData',15)
    legend('targeted cells')
else
    scatter(trg.laser_y(:,clusterID),trg.laser_x(:,clusterID),'co','LineWidth',1,'SizeData',50)    
    scatter(meds(temp,1),meds(temp,2),'k+','LineWidth',2,'SizeData',30)
    legend('passed stats (-)','passed stats (+)','passed stats and amp (+)','laser','targeted')
end

% plotting stuff for the lab meeting
if false
    scatter(meds(find(resp.responders_main.respAll(:,clusterID)==1),1),meds(find(resp.responders_main.respAll(:,clusterID)==1),2),'g*','LineWidth',1,'SizeData',70)
    scatter(meds(find(flw.stats_sig_pos(:,clusterID)==1),1),meds(find(flw.stats_sig_pos(:,clusterID)==1),2),'r.','LineWidth',1,'SizeData',50)
    scatter(meds(find(flw.stats_sig_neg(:,clusterID)==1),1),meds(find(flw.stats_sig_neg(:,clusterID)==1),2),'b.','LineWidth',1,'SizeData',50)
    temp = resp.responders_main.targetedCells(:,clusterID);
    temp = temp(~isnan(temp));
    scatter(meds(temp,1),meds(temp,2),'co','LineWidth',1,'SizeData',50)
    scatter(trg.laser_y(:,clusterID),trg.laser_x(:,clusterID),'k+','LineWidth',2,'SizeData',50)

end

if exclusionZone
    for k=1:6
        plot((50*info.scope.fovSize_pix/info.scope.fovSize_um)*cos(0:pi/50:2*pi)+trg.laser_y(k,clusterID),(50*info.scope.fovSize_pix/info.scope.fovSize_um)*sin(0:pi/50:2*pi)+trg.laser_x(k,clusterID),'k-','LineWidth',0.5)
    end
end

xlim([0,s2p_meta.ops.Lx])
ylim([0,s2p_meta.ops.Ly])
if ~paper
    title(titleText)
end

box on;
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);

if ~issubplot
    drawnow;
end
end

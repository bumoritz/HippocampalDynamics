function F = clusterResponseMap(s2p_meta,iscell,trg,resp,clusterID,issubplot,plt)

img = zeros(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
for n=1:length(iscell)
    if iscell(n)
        for i=1:length(s2p_meta.stat{n}.xpix)
            img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = resp.avgAct_net(n,clusterID);
        end
    end
end

meds = zeros(length(s2p_meta.stat),2);
for i = 1:length(s2p_meta.stat)
    meds(i,:) = [double(s2p_meta.stat{i}.med(2)),double(s2p_meta.stat{i}.med(1))];
end

if nargin < 7
    c_abs = 1;%3;
    c_label = '\Delta_{resp-base} z-scored \DeltaF/F'; %'\Delta_{resp-base} \DeltaF/F';
    titleText = ['Cluster ',num2str(clusterID)];
else
    c_abs = plt.c_abs;
    c_label = plt.c_label;
    titleText = plt.titleText;
end

if ~issubplot
    F = default_figure([20,0.5,20,9.9]);
else
    F = NaN;
end

imagesc(img);
set(gca,'CLim',[-c_abs,c_abs]);
colormap(redblue);
temp=colorbar; 
temp.Label.String = c_label;
daspect([1,1,1]);
hold on

if issubplot
%    scatter(meds(find(resp.stats_sig_pos(:,clusterID)==1),1),meds(find(resp.stats_sig_pos(:,clusterID)==1),2),'g*','LineWidth',0.5,'SizeData',25)
else
    scatter(meds(find(resp.stats_sig_neg(:,clusterID)==1),1),meds(find(resp.stats_sig_neg(:,clusterID)==1),2),'y.','LineWidth',1,'SizeData',50)
    scatter(meds(find(resp.stats_sig_pos(:,clusterID)==1),1),meds(find(resp.stats_sig_pos(:,clusterID)==1),2),'g.','LineWidth',1,'SizeData',50)
    scatter(meds(find(resp.responders_main.respAll(:,clusterID)==1),1),meds(find(resp.responders_main.respAll(:,clusterID)==1),2),'g*','LineWidth',1,'SizeData',70)
end
temp = resp.responders_main.targetedCells(:,clusterID);
temp = temp(~isnan(temp));
if issubplot
    scatter(meds(temp,1),meds(temp,2),'k+','LineWidth',1,'SizeData',15)
    legend('targeted cells')
else
    scatter(trg.laser_y(:,clusterID),trg.laser_x(:,clusterID),'co','LineWidth',1,'SizeData',50)    
    scatter(meds(temp,1),meds(temp,2),'k+','LineWidth',2,'SizeData',30)
    legend('passed stats (-)','passed stats (+)','passed stats and amp (+)','laser','targeted')
end

xlim([0,s2p_meta.ops.Lx])
ylim([0,s2p_meta.ops.Ly])
title(titleText)

set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);

if ~issubplot
    drawnow;
end
end

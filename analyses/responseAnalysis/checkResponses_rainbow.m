%function F = checkResponses_rainbow(trg,resp,p,plt)
stimType = 'seq';

if strcmp(stimType,'seq')
    sta = resp.avgSta_seq;
    cols = p.col.seq_rainbow;
    clusters = trg.grouping(:,1);
elseif strcmp(stimType,'ctrl')
    sta = resp.avgSta_ctrl;
    cols = p.col.ctrl_rainbow;
    clusters = trg.grouping(:,2);
end

F = default_figure([20,0.5,20,9.9]);

for i=1:length(clusters)
    temp = sta(rmmissing(trg.idcs_targetedCells(:,clusters(i))),:);
    for j=1:size(temp,1)
        plot(resp.sta_t,temp(j,:)+length(clusters)+1-i,'Color',cols(i,:))
        hold on
    end
end
hold off
xlim([floor(resp.sta_t(1)),ceil(resp.sta_t(end))])
ylim([0,length(clusters)+3])
title('STA sequence')
xlabel('time (s)')
ylabel('z-scored \DeltaF/F (stacked by stim cluster)')











%%

img = zeros(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
for n=1:length(s2p_meta.iscell(:,1))
    if s2p_meta.iscell(n,1)
        for i=1:length(s2p_meta.stat{n}.xpix)
            img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = resp.avgAct_net_byCluster(n,clusterID);
        end
    end
end

meds = zeros(length(s2p_meta.stat),2);
for i = 1:length(s2p_meta.stat)
    meds(i,:) = [double(s2p_meta.stat{i}.med(2)),double(s2p_meta.stat{i}.med(1))];
end

if nargin < 5
    c_abs = 3;
    c_label = '\Delta_{resp-base} z-scored \DeltaF/F';
    suptitleText = ['Cluster ',num2str(clusterID)];
else
    c_abs = plt.c_abs;
    c_label = plt.c_label;
    suptitleText = plt.suptitleText;
end

F = default_figure([20,0.5,20,9.9]);

imagesc(img);
set(gca,'CLim',[-c_abs,c_abs]);
colormap(redblue);
temp=colorbar; 
temp.Label.String = c_label;
daspect([1,1,1]);
hold on

scatter(meds(find(resp.stats_sig(:,clusterID)==1),1),meds(find(resp.stats_sig(:,clusterID)==1),2),'g*','LineWidth',1,'SizeData',50)
temp = trg.idcs_targetedCells(:,clusterID);
temp = temp(~isnan(temp));
scatter(meds(temp,1),meds(temp,2),'k+','LineWidth',2,'SizeData',50)

xlim([0,s2p_meta.ops.Lx])
ylim([0,s2p_meta.ops.Ly])
legend('significant','targeted')
suptitle(suptitleText)

drawnow;
%end

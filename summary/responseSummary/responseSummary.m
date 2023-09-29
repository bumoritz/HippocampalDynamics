%% Get data

this_struct = 'flw'; % 'resp'

responders = nan(d_info.numAnimals,d_info.numDays);
responders_seq = nan(d_info.numAnimals,d_info.numDays);
responders_ctrl = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            responders(i,j) = d{i,j}.(this_struct).responders_main.numRespTargeted / d{i,j}.(this_struct).responders_main.numIdentifiedTargeted; %d{i,j}.(this_struct).responders_main.numRespTargeted;
            if d_info.group(i)==7 && (j==2 || j==3)
                responders_seq(i,j) = responders(i,j);
            elseif d_info.group(i)==8 && (j==4 || j==5)
                responders_seq(i,j) = responders(i,j);
            elseif d_info.group(i)==7 && (j==4 || j==5)
                responders_ctrl(i,j) = responders(i,j);
            elseif d_info.group(i)==8 && (j==2 || j==3)
            	responders_ctrl(i,j) = responders(i,j);
            end
        catch
        end
    end
end
responders_seq_cmpr = nanmean(responders_seq(:,[2,4]),2); %nanmean(responders_seq(:,2:5),2);
responders_ctrl_cmpr = nanmean(responders_ctrl(:,[2,4]),2); %nanmean(responders_ctrl(:,2:5),2);
responders_pooled_cmpr = nanmean([responders_seq_cmpr,responders_ctrl_cmpr],2);

offTargets = nan(d_info.numAnimals,d_info.numDays);
offTargets_seq = nan(d_info.numAnimals,d_info.numDays);
offTargets_ctrl = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            offTargets(i,j) = d{i,j}.(this_struct).responders_main.proportion_RespNonTargeted_RespTargeted;
            if d_info.group(i)==7 && (j==2 || j==3)
                offTargets_seq(i,j) = offTargets(i,j);
            elseif d_info.group(i)==8 && (j==4 || j==5)
                offTargets_seq(i,j) = offTargets(i,j);
            elseif d_info.group(i)==7 && (j==4 || j==5)
                offTargets_ctrl(i,j) = offTargets(i,j);
            elseif d_info.group(i)==8 && (j==2 || j==3)
            	offTargets_ctrl(i,j) = offTargets(i,j);
            end
        catch
        end
    end
end
offTargets_seq_cmpr = nanmean(offTargets_seq(:,[2,4]),2); %nanmean(offTargets_seq(:,2:5),2);
offTargets_ctrl_cmpr = nanmean(offTargets_ctrl(:,[2,4]),2); %nanmean(offTargets_ctrl(:,2:5),2);
offTargets_pooled_cmpr = nanmean([offTargets_seq_cmpr,offTargets_ctrl_cmpr],2);

stability = nan(d_info.numAnimals,d_info.numDays,40);
stability_seq = nan(d_info.numAnimals,d_info.numDays,40);
stability_ctrl = nan(d_info.numAnimals,d_info.numDays,40);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            temp = nanmean(d{i,j}.(this_struct).trgRespAmps_cells_bw,1);
            stability(i,j,1:length(temp)) = temp;
            if d_info.group(i)==7 && (j==2 || j==3)
                stability_seq(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==8 && (j==4 || j==5)
                stability_seq(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==7 && (j==4 || j==5)
                stability_ctrl(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==8 && (j==2 || j==3)
            	stability_ctrl(i,j,1:length(temp)) = temp;
            end
        catch
        end
    end
end
stability_seq_cmpr = squeeze(nanmean(stability_seq(:,[2,4],:),2)); %squeeze(nanmean(stability_seq(:,2:5,:),2));
stability_ctrl_cmpr = squeeze(nanmean(stability_ctrl(:,[2,4],:),2)); %squeeze(nanmean(stability_ctrl(:,2:5,:),2));
stability_pooled_cmpr = nanmean(cat(3,stability_seq_cmpr,stability_ctrl_cmpr),3);

resolution = nan(d_info.numAnimals,d_info.numDays,20);
resolution_seq = nan(d_info.numAnimals,d_info.numDays,20);
resolution_ctrl = nan(d_info.numAnimals,d_info.numDays,20);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            temp1 = d{i,j}.(this_struct).responders_ext.respAll(:); %d{i,j}.(this_struct).responders_ext.respAll(:); %d{i,j}.(this_struct).responders_main.respAll(:);
            temp2 = discretize(d{i,j}.resp.dist_closestLaser(:),[0:2.5:50]);
            temp = nan(nanmax(temp2),1);
            for k=1:nanmax(temp2)
                temp(k) = nanmean(temp1(find(temp2==k)));
            end
            resolution(i,j,1:length(temp)) = temp;
            if d_info.group(i)==7 && (j==2 || j==3)
                resolution_seq(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==8 && (j==4 || j==5)
                resolution_seq(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==7 && (j==4 || j==5)
                resolution_ctrl(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==8 && (j==2 || j==3)
            	resolution_ctrl(i,j,1:length(temp)) = temp;
            end
        catch
        end
    end
end
resolution_seq_cmpr = squeeze(nanmean(resolution_seq(:,[2,4],:),2)); %squeeze(nanmean(resolution_seq(:,2:5,:),2));
resolution_ctrl_cmpr = squeeze(nanmean(resolution_ctrl(:,[2,4],:),2)); %squeeze(nanmean(resolution_ctrl(:,2:5,:),2));
resolution_pooled_cmpr = nanmean(cat(3,resolution_seq_cmpr,resolution_ctrl_cmpr),3);

inhibition = nan(d_info.numAnimals,d_info.numDays,40);
inhibition_seq = nan(d_info.numAnimals,d_info.numDays,40);
inhibition_ctrl = nan(d_info.numAnimals,d_info.numDays,40);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            temp1 = d{i,j}.(this_struct).stats_sig_neg(:);
            temp2 = discretize(d{i,j}.resp.dist_closestLaser(:),[0:10:400]);
            temp = nan(nanmax(temp2),1);
            for k=1:nanmax(temp2)
                temp(k) = nanmean(temp1(find(temp2==k)));
            end
            inhibition(i,j,1:length(temp)) = temp;
            if d_info.group(i)==7 && (j==2 || j==3)
                inhibition_seq(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==8 && (j==4 || j==5)
                inhibition_seq(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==7 && (j==4 || j==5)
                inhibition_ctrl(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==8 && (j==2 || j==3)
            	inhibition_ctrl(i,j,1:length(temp)) = temp;
            end
        catch
        end
    end
end
inhibition_seq_cmpr = squeeze(nanmean(inhibition_seq(:,2:5,:),2));
inhibition_ctrl_cmpr = squeeze(nanmean(inhibition_ctrl(:,2:5,:),2));
inhibition_pooled_cmpr = nanmean(cat(3,inhibition_seq_cmpr,inhibition_ctrl_cmpr),3);

excitation = nan(d_info.numAnimals,d_info.numDays,40);
excitation_seq = nan(d_info.numAnimals,d_info.numDays,40);
excitation_ctrl = nan(d_info.numAnimals,d_info.numDays,40);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            temp1 = d{i,j}.(this_struct).stats_sig_pos(:);
            temp2 = discretize(d{i,j}.resp.dist_closestLaser(:),[0:10:400]);
            temp = nan(nanmax(temp2),1);
            for k=1:nanmax(temp2)
                temp(k) = nanmean(temp1(find(temp2==k)));
            end
            excitation(i,j,1:length(temp)) = temp;
            if d_info.group(i)==7 && (j==2 || j==3)
                excitation_seq(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==8 && (j==4 || j==5)
                excitation_seq(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==7 && (j==4 || j==5)
                excitation_ctrl(i,j,1:length(temp)) = temp;
            elseif d_info.group(i)==8 && (j==2 || j==3)
            	excitation_ctrl(i,j,1:length(temp)) = temp;
            end
        catch
        end
    end
end
excitation_seq_cmpr = squeeze(nanmean(excitation_seq(:,2:5,:),2));
excitation_ctrl_cmpr = squeeze(nanmean(excitation_ctrl(:,2:5,:),2));
excitation_pooled_cmpr = nanmean(cat(3,excitation_seq_cmpr,excitation_ctrl_cmpr),3);


%% Measures that are only in flw struct

avgPosAmp_by_avgTrgAmpBin = nan(d_info.numAnimals,d_info.numDays,length(p.resp.ampBins_x));
avgPosAmp_by_avgTrgAmpBin_seq = nan(d_info.numAnimals,d_info.numDays,length(p.resp.ampBins_x));
avgPosAmp_by_avgTrgAmpBin_ctrl = nan(d_info.numAnimals,d_info.numDays,length(p.resp.ampBins_x));
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            avgPosAmp_by_avgTrgAmpBin(i,j,:) = d{i,j}.flw.ampByResp.avgPosAmp_by_avgTrgAmpBin;
            if d_info.group(i)==7 && (j==2 || j==3)
                avgPosAmp_by_avgTrgAmpBin_seq(i,j,:) = temp;
            elseif d_info.group(i)==8 && (j==4 || j==5)
                avgPosAmp_by_avgTrgAmpBin_seq(i,j,:) = temp;
            elseif d_info.group(i)==7 && (j==4 || j==5)
                avgPosAmp_by_avgTrgAmpBin_ctrl(i,j,:) = temp;
            elseif d_info.group(i)==8 && (j==2 || j==3)
            	avgPosAmp_by_avgTrgAmpBin_ctrl(i,j,:) = temp;
            end
        catch
        end
    end
end
avgPosAmp_by_avgTrgAmpBin_seq_cmpr = squeeze(nanmean(avgPosAmp_by_avgTrgAmpBin_seq(:,2:5,:),2));
avgPosAmp_by_avgTrgAmpBin_ctrl_cmpr = squeeze(nanmean(avgPosAmp_by_avgTrgAmpBin_ctrl(:,2:5,:),2));
avgPosAmp_by_avgTrgAmpBin_pooled_cmpr = nanmean(cat(3,avgPosAmp_by_avgTrgAmpBin_seq_cmpr,avgPosAmp_by_avgTrgAmpBin_ctrl_cmpr),3);


%% Response summary figure - seq vs ctrl

these_animals = ~any(d_info.presponsive==0,2);
these_labels = categorical({'seq','ctrl'});
these_labels = reordercats(these_labels,{'seq','ctrl'});

nrows = 1; ncols = 4; m=0;
F = default_figure([-20,0.5,20,3.5]);

% a) Target responses
m = m+1; subplot(nrows,ncols,m); hold on;
this_data_seq = responders_seq_cmpr(these_animals);
this_data_ctrl = responders_ctrl_cmpr(these_animals);
v = bar(these_labels,nanmean([this_data_seq,this_data_ctrl],1));
for i=1:length(this_data_seq)
    plot(these_labels,[this_data_seq,this_data_ctrl],'ok')
end
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.ctrl;
% yline(240,'k:');
% ylim([0,240])
% yticks([0,40,80,120,160,200,240])
% ylabel('Number of responsive targets')
ylabel('Proportion of responsive targets')
title('Responsive targets')

% b) Off-target responses
m = m+1; subplot(nrows,ncols,m); hold on;
this_data_seq = offTargets_seq_cmpr(these_animals);
this_data_ctrl = offTargets_ctrl_cmpr(these_animals);
v = bar(these_labels,nanmean([this_data_seq,this_data_ctrl],1));
for i=1:length(this_data_seq)
    plot(these_labels,[this_data_seq,this_data_ctrl],'ok')
end
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.ctrl;
ylabel('Number of resp. off-targets per resp. target')
title('Responsive off-targets')

% c) Resolution
m = m+1; subplot(nrows,ncols,m); hold on;
this_data_ctrl = resolution_ctrl_cmpr(these_animals,:);
this_data_seq = resolution_seq_cmpr(these_animals,:);
shadedErrorBar(1:size(this_data_ctrl,2),nanmean(this_data_ctrl,1)*100,nansem(this_data_ctrl,1)*100,'lineProps',p.col.ctrl);
shadedErrorBar(1:size(this_data_seq,2),nanmean(this_data_seq,1)*100,nansem(this_data_seq,1)*100,'lineProps',p.col.seq);
xlim([1,size(this_data_seq,2)+1])
xticks([0:4:length([0:2.5:50])-1]+1)
xticklabels({'0','10','20','30','40','50'})
ylim([0,100])
ytickformat('percentage')
xlabel('Distance from closest laser beamlet (\mum)')
ylabel('Response probability')
title('Photostimulation resolution')

% d) Stability
m = m+1; subplot(nrows,ncols,m); hold on;
this_data_ctrl = stability_ctrl_cmpr(these_animals,:);
this_data_seq = stability_seq_cmpr(these_animals,:);
shadedErrorBar(1:size(this_data_ctrl,2),nanmean(this_data_ctrl,1),nansem(this_data_ctrl,1),'lineProps',p.col.ctrl);
shadedErrorBar(1:size(this_data_seq,2),nanmean(this_data_seq,1),nansem(this_data_seq,1),'lineProps',p.col.seq);
xlim([0,size(this_data_seq,2)])
ylim([0,inf])
xlabel('Trial block')
ylabel(['Mean response amplitude of targeted cells',newline,'(z-scored \DeltaF/F)'])
title('Photostimulation stability')

savefig(F,[path.root_summary,'plots\responseSummary\responseSummary_daysPooled_seqvsctrl.fig']);
saveas(F,[path.root_summary,'plots\responseSummary\responseSummary_daysPooled_seqvsctrl.png']);


%% Response summary figure

these_animals = ~any(d_info.presponsive==0,2);
these_labels = categorical({'stim'});
these_labels = reordercats(these_labels,{'stim'});

nrows = 1; ncols = 4; m=0;
F = default_figure([-20,0.5,20,3.5]);

% a) Target responses
m = m+1; subplot(nrows,ncols,m); hold on;
this_data = responders_pooled_cmpr(these_animals);
v = bar(these_labels,nanmean(this_data,1));
for i=1:length(this_data)
    plot(these_labels,this_data,'ok')
end
v.FaceColor = 'flat'; v.CData(1,:) = p.col.photostim;
% yline(240,'k:');
% ylim([0,240])
% yticks([0,40,80,120,160,200,240])
ylabel('Number of responsive targets')
title('Responsive targets')

% b) Off-target responses
m = m+1; subplot(nrows,ncols,m); hold on;
this_data = offTargets_pooled_cmpr(these_animals);
v = bar(these_labels,nanmean(this_data,1));
for i=1:length(this_data)
    plot(these_labels,this_data,'ok')
end
v.FaceColor = 'flat'; v.CData(1,:) = p.col.photostim;
ylabel('Number of resp. off-targets per resp. target')
title('Responsive off-targets')

% c) Resolution
m = m+1; subplot(nrows,ncols,m); hold on;
% this_data = resolution_pooled_cmpr(these_animals,:);
% shadedErrorBar(1:size(this_data,2),nanmean(this_data,1)*100,nansem(this_data,1)*100,'lineProps',p.col.darkGray);
this_data = resolution_pooled_cmpr(these_animals,:);
shadedErrorBar(1:size(this_data,2),nanmean(this_data,1)*100,nansem(this_data,1)*100,'lineProps',p.col.photostim);
xlim([1,size(this_data,2)+1])
xticks([0:4:length([0:2.5:50])-1]+1)
xticklabels({'0','10','20','30','40','50'})
ylim([0,100])
ytickformat('percentage')
xlabel('Distance from closest laser beamlet (\mum)')
ylabel('Response probability')
title('Photostimulation resolution')

% d) Stability
m = m+1; subplot(nrows,ncols,m); hold on;
this_data = stability_pooled_cmpr(these_animals,:);
shadedErrorBar(1:size(this_data,2),nanmean(this_data,1),nansem(this_data,1),'lineProps',p.col.photostim);
xlim([0,size(this_data,2)])
ylim([0,inf])
xlabel('Trial block')
ylabel(['Mean response amplitude of targeted cells',newline,'(z-scored \DeltaF/F)'])
title('Photostimulation stability')

savefig(F,[path.root_summary,'plots\responseSummary\responseSummary_daysPooled_pooled.fig']);
saveas(F,[path.root_summary,'plots\responseSummary\responseSummary_daysPooled_pooled.png']);


%% Response summary figure - seq vs ctrl

these_animals = ~any(d_info.presponsive==0,2);
these_labels = categorical({'seq','ctrl'});
these_labels = reordercats(these_labels,{'seq','ctrl'});

nrows = 1; ncols = 4; m=0;
F = default_figure([-20,0.5,20,3.5]);

% a) Negative followers
m = m+1; subplot(nrows,ncols,m); hold on;
this_data_ctrl = inhibition_ctrl_cmpr(these_animals,:);
this_data_seq = inhibition_seq_cmpr(these_animals,:);
plot(this_data_ctrl'*100,'Color',p.col.ctrl)
plot(this_data_seq'*100,'Color',p.col.seq)
xticks([0:(length([0:10:400]))/8:length([0:10:400])]+1)
xticklabels({'0','50','100','150','200','250','300','350','400'})
xlim([1,length([0:10:400])+1])
% ylim([0,100])
ytickformat('percentage')
xlabel('Distance from closest laser beamlet (\mum)')
ylabel('Response probability')
title('Negative followers')

% b) Negative followers
m = m+1; subplot(nrows,ncols,m); hold on;
this_data_ctrl = inhibition_ctrl_cmpr(these_animals,:);
this_data_seq = inhibition_seq_cmpr(these_animals,:);
shadedErrorBar(1:size(this_data_ctrl,2),nanmean(this_data_ctrl,1)*100,nansem(this_data_ctrl,1)*100,'lineProps',p.col.ctrl);
shadedErrorBar(1:size(this_data_seq,2),nanmean(this_data_seq,1)*100,nansem(this_data_seq,1)*100,'lineProps',p.col.seq);
xticks([0:(length([0:10:400]))/8:length([0:10:400])]+1)
xticklabels({'0','50','100','150','200','250','300','350','400'})
xlim([1,length([0:10:400])+1])
% ylim([0,100])
ytickformat('percentage')
xlabel('Distance from closest laser beamlet (\mum)')
ylabel('Response probability')
title('Negative followers')

% c) Positive followers
m = m+1; subplot(nrows,ncols,m); hold on;
this_data_ctrl = excitation_ctrl_cmpr(these_animals,:);
this_data_seq = excitation_seq_cmpr(these_animals,:);
plot(this_data_ctrl'*100,'Color',p.col.ctrl)
plot(this_data_seq'*100,'Color',p.col.seq)
xticks([0:(length([0:10:400]))/8:length([0:10:400])]+1)
xticklabels({'0','50','100','150','200','250','300','350','400'})
xlim([1,length([0:10:400])+1])
% ylim([0,100])
ytickformat('percentage')
xlabel('Distance from closest laser beamlet (\mum)')
ylabel('Response probability')
title('Positive followers')

% d) Positive followers
m = m+1; subplot(nrows,ncols,m); hold on;
this_data_ctrl = excitation_ctrl_cmpr(these_animals,:);
this_data_seq = excitation_seq_cmpr(these_animals,:);
shadedErrorBar(1:size(this_data_ctrl,2),nanmean(this_data_ctrl,1)*100,nansem(this_data_ctrl,1)*100,'lineProps',p.col.ctrl);
shadedErrorBar(1:size(this_data_seq,2),nanmean(this_data_seq,1)*100,nansem(this_data_seq,1)*100,'lineProps',p.col.seq);
xticks([0:(length([0:10:400]))/8:length([0:10:400])]+1)
xticklabels({'0','50','100','150','200','250','300','350','400'})
xlim([1,length([0:10:400])+1])
% ylim([0,100])
ytickformat('percentage')
xlabel('Distance from closest laser beamlet (\mum)')
ylabel('Response probability')
title('Positive followers')




%%  (ctrl)
% 
% these_animals = ~any(d_info.presponsive==0,2);
% these_labels = categorical({'seq','ctrl'});
% these_labels = reordercats(these_labels,{'seq','ctrl'});
% 
% nrows = 1; ncols = 4; m=0;
% F = default_figure([-20,0.5,20,3.5]);
% 
% % a) Negative followers
% m = m+1; subplot(nrows,ncols,m); hold on;
% this_data_ctrl = inhibition_ctrl_cmpr(these_animals,:);
% this_data_seq = inhibition_seq_cmpr(these_animals,:);
% plot(this_data_ctrl'*100,'Color',p.col.ctrl)
% plot(this_data_seq'*100,'Color',p.col.seq)
% xticks([0:(length([0:10:400]))/8:length([0:10:400])]+1)
% xticklabels({'0','50','100','150','200','250','300','350','400'})
% xlim([1,length([0:10:400])+1])
% % ylim([0,100])
% ytickformat('percentage')
% xlabel('Distance from closest laser beamlet (\mum)')
% ylabel('Response probability')
% title('Negative followers')


%% FIGURE DUMP

%% Fig5_ResponseAnalysis_StabilityAcrossDays

F = paper_figure([0,0.5,mm2inch(50),mm2inch(50)]); hold on;

this_data_seq_1 = amplitude_seq_paired_sd12(:,1,1);
this_data_seq_2 = amplitude_seq_paired_sd12(:,1,2);
this_data_ctrl_1 = amplitude_ctrl_paired_sd12(:,1,1);
this_data_ctrl_2 = amplitude_ctrl_paired_sd12(:,1,2);
v = bar(these_labels2,nanmean([this_data_seq_1,this_data_seq_2,this_data_ctrl_1,this_data_ctrl_2],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq_1)
    plot(these_labels2,[this_data_seq_1(i),this_data_seq_2(i)],'-k','LineWidth',1)
    plot(these_labels2,[this_data_ctrl_1(i),this_data_ctrl_2(i)],'-k','LineWidth',1)
end
ylim([0,3])
yticks([0,1,2,3])
ylabel('Mean target response amplitude (S)')

savefig(F,[save_root_fig,'\Fig5_ResponseAnalysis_StabilityAcrossDays.fig']);
saveas(F,[save_root_png,'\Fig5_ResponseAnalysis_StabilityAcrossDays.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ResponseAnalysis_StabilityAcrossDays.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_ResponseAnalysis_StabilityAcrossDays.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_StabilityAcrossDays\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nsignrank p=',num2str(signrank(this_data_seq,this_data_ctrl),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_StabilityAcrossDays\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nsignrank p=',num2str(signrank(this_data_seq,this_data_ctrl),4),'\n']);


%% Fig5_ResponseAnalysis_StabilityAcrossDays

F = paper_figure([0,0.5,mm2inch(50),mm2inch(50)]); hold on;

this_data_seq = amplitude_seq_paired_sd12(:,1,2) ./ amplitude_seq_paired_sd12(:,1,1);
this_data_ctrl = amplitude_ctrl_paired_sd12(:,1,2) ./ amplitude_ctrl_paired_sd12(:,1,1);
v = bar(these_labels,nanmean([this_data_seq,this_data_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.ctrl; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq)
    plot(these_labels,[this_data_seq(i),this_data_ctrl(i)],'-k','LineWidth',1)
end
yline(100,'k:','LineWidth',1);
ylim([0,150])
yticks([0,50,100,150])
ytickformat('percentage')
ylabel('Proportion of initial')



ylim([0,3])
yticks([0,1,2,3])
ylabel('Mean target response amplitude (S)')

savefig(F,[save_root_fig,'\Fig5_ResponseAnalysis_StabilityAcrossDays.fig']);
saveas(F,[save_root_png,'\Fig5_ResponseAnalysis_StabilityAcrossDays.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ResponseAnalysis_StabilityAcrossDays.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_ResponseAnalysis_StabilityAcrossDays.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_StabilityAcrossDays\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nsignrank p=',num2str(signrank(this_data_seq,this_data_ctrl),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_StabilityAcrossDays\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nsignrank p=',num2str(signrank(this_data_seq,this_data_ctrl),4),'\n']);


F = paper_figure([0,0.5,mm2inch(50),mm2inch(50)]); hold on;

this_data_seq_1 = amplitude_seq_paired_sd12(:,1,1);
this_data_seq_2 = amplitude_seq_paired_sd12(:,1,2);
this_data_ctrl_1 = amplitude_ctrl_paired_sd12(:,1,1);
this_data_ctrl_2 = amplitude_ctrl_paired_sd12(:,1,2);
v = bar(these_labels2,nanmean([this_data_seq_1,this_data_seq_2,this_data_ctrl_1,this_data_ctrl_2],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq_1)
    plot(these_labels2,[this_data_seq_1(i),this_data_seq_2(i),this_data_ctrl_1(i),this_data_ctrl_2(i)],'-k','LineWidth',1)
end
ylim([0,3])
yticks([0,1,2,3])
ylabel('Mean target response amplitude (S)')

savefig(F,[save_root_fig,'\Fig5_ResponseAnalysis_StabilityAcrossDays.fig']);
saveas(F,[save_root_png,'\Fig5_ResponseAnalysis_StabilityAcrossDays.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ResponseAnalysis_StabilityAcrossDays.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_ResponseAnalysis_StabilityAcrossDays.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_StabilityAcrossDays\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nsignrank p=',num2str(signrank(this_data_seq,this_data_ctrl),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_StabilityAcrossDays\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nsignrank p=',num2str(signrank(this_data_seq,this_data_ctrl),4),'\n']);

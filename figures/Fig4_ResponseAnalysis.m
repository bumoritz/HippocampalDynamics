%% Fig4_ResponseAnalysis

% import data using Summary_Master with ops.do_responseSummary = true;
% d_info = selectDataset(d_info,'-g78-d2345-e01-r01-p01-l01-i00',sheet,path,ops);

save_root_fig = [path.root_summary,'figures\Fig4_fig\'];
save_root_png = [path.root_summary,'figures\Fig4_png\'];
save_root_pdf = [path.root_summary,'figures\Fig4_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig4_txt\'];


%% Get data

this_struct = 'resp';
these_animals_paired_firstDay = ~any(d_info.presponsive(:,[2,4])==0,2);
these_animals_paired_bothDays = ~any(d_info.presponsive(:,2:5)==0,2);

running_seq = nan(d_info.numAnimals,d_info.numDays);
running_ctrl = nan(d_info.numAnimals,d_info.numDays);
running_unpaired = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            running_unpaired(i,j) = d_info.running(i,j);
            try
                if d_info.group(i)==7 && (j==2 || j==3)
                    running_seq(i,j) = d_info.running(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    running_seq(i,j) = d_info.running(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    running_ctrl(i,j) = d_info.running(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    running_ctrl(i,j) = d_info.running(i,j);
                end
            catch
            end
        end
    end
end
running_seq_unpaired_sd1 = nanmean(running_seq(:,[2,4]),2);
running_seq_paired_sd1 = nanmean(running_seq(:,[2,4]),2);
running_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
running_ctrl_unpaired_sd1 = nanmean(running_ctrl(:,[2,4]),2);
running_ctrl_paired_sd1 = nanmean(running_ctrl(:,[2,4]),2);
running_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;

engagement_seq = nan(d_info.numAnimals,d_info.numDays);
engagement_ctrl = nan(d_info.numAnimals,d_info.numDays);
engagement_unpaired = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            engagement_unpaired(i,j) = d_info.engagement(i,j);
            try
                if d_info.group(i)==7 && (j==2 || j==3)
                    engagement_seq(i,j) = d_info.engagement(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    engagement_seq(i,j) = d_info.engagement(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    engagement_ctrl(i,j) = d_info.engagement(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    engagement_ctrl(i,j) = d_info.engagement(i,j);
                end
            catch
            end
        end
    end
end
engagement_seq_unpaired_sd1 = nanmean(engagement_seq(:,[2,4]),2);
engagement_seq_paired_sd1 = nanmean(engagement_seq(:,[2,4]),2);
engagement_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
engagement_ctrl_unpaired_sd1 = nanmean(engagement_ctrl(:,[2,4]),2);
engagement_ctrl_paired_sd1 = nanmean(engagement_ctrl(:,[2,4]),2);
engagement_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;

responders = nan(d_info.numAnimals,d_info.numDays);
responders_seq = nan(d_info.numAnimals,d_info.numDays);
responders_ctrl = nan(d_info.numAnimals,d_info.numDays);
responders_unpaired = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                responders(i,j) = d{i,j}.(this_struct).responders_main.numRespAll/40;
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
            responders_unpaired(i,j) = responders(i,j);
        end
    end
    if (isnan(responders(i,2)) || isnan(responders(i,4)))
        these_animals_paired_firstDay(i) = 0;
    end
    if any([isnan(responders(i,2)),isnan(responders(i,3)),isnan(responders(i,4)),isnan(responders(i,5))])
        these_animals_paired_bothDays(i) = 0;
    end
end
responders_seq_unpaired_sd1 = nanmean(responders_seq(:,[2,4]),2);
responders_seq_paired_sd1 = nanmean(responders_seq(:,[2,4]),2);
responders_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
responders_ctrl_unpaired_sd1 = nanmean(responders_ctrl(:,[2,4]),2);
responders_ctrl_paired_sd1 = nanmean(responders_ctrl(:,[2,4]),2);
responders_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
responders_seq_paired_sd12 = [nanmean(responders_seq(:,[2,4]),2),nanmean(responders_seq(:,[3,5]),2)];
responders_seq_paired_sd12(~these_animals_paired_bothDays) = NaN;
responders_ctrl_paired_sd12 = [nanmean(responders_ctrl(:,[2,4]),2),nanmean(responders_ctrl(:,[3,5]),2)];
responders_ctrl_paired_sd12(~these_animals_paired_bothDays) = NaN;

targets = nan(d_info.numAnimals,d_info.numDays);
targets_seq = nan(d_info.numAnimals,d_info.numDays);
targets_ctrl = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                targets(i,j) = d{i,j}.(this_struct).responders_main.numRespTargeted/40; % / d{i,j}.(this_struct).responders_main.numIdentifiedTargeted;
                if d_info.group(i)==7 && (j==2 || j==3)
                    targets_seq(i,j) = targets(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    targets_seq(i,j) = targets(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    targets_ctrl(i,j) = targets(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    targets_ctrl(i,j) = targets(i,j);
                end
            catch
            end
        end
    end
end
targets_seq_paired_sd1 = nanmean(targets_seq(:,[2,4]),2);
targets_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
targets_ctrl_paired_sd1 = nanmean(targets_ctrl(:,[2,4]),2);
targets_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;

offTargets = nan(d_info.numAnimals,d_info.numDays);
offTargets_seq = nan(d_info.numAnimals,d_info.numDays);
offTargets_ctrl = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                offTargets(i,j) = d{i,j}.(this_struct).responders_main.numRespNonTargeted/40;
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
end
offTargets_seq_paired_sd1 = nanmean(offTargets_seq(:,[2,4]),2);
offTargets_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
offTargets_ctrl_paired_sd1 = nanmean(offTargets_ctrl(:,[2,4]),2);
offTargets_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;

offTargetsPerTarget = nan(d_info.numAnimals,d_info.numDays);
offTargetsPerTarget_seq = nan(d_info.numAnimals,d_info.numDays);
offTargetsPerTarget_ctrl = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                offTargetsPerTarget(i,j) = d{i,j}.(this_struct).responders_main.proportion_RespNonTargeted_RespTargeted;
                if d_info.group(i)==7 && (j==2 || j==3)
                    offTargetsPerTarget_seq(i,j) = offTargetsPerTarget(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    offTargetsPerTarget_seq(i,j) = offTargetsPerTarget(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    offTargetsPerTarget_ctrl(i,j) = offTargetsPerTarget(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    offTargetsPerTarget_ctrl(i,j) = offTargetsPerTarget(i,j);
                end
            catch
            end
        end
    end
end
offTargetsPerTarget_seq_paired_sd1 = nanmean(offTargetsPerTarget_seq(:,[2,4]),2);
offTargetsPerTarget_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
offTargetsPerTarget_ctrl_paired_sd1 = nanmean(offTargetsPerTarget_ctrl(:,[2,4]),2);
offTargetsPerTarget_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;

amplitude = nan(d_info.numAnimals,d_info.numDays,8);
amplitude_unpaired = nan(d_info.numAnimals,d_info.numDays,8);
amplitude_seq = nan(d_info.numAnimals,d_info.numDays,8);
amplitude_ctrl = nan(d_info.numAnimals,d_info.numDays,8);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp0 = nanmean(d{i,j}.(this_struct).respAmps_cells_bw,1); %nanmean(d{i,j}.(this_struct).trgRespAmps_cells_bw,1);
                temp = nanmean(reshape(temp0,5,length(temp0)/5),1);
                amplitude(i,j,1:length(temp)) = temp;
                if (d_info.group(i)==7 && (j==2 || j==3))
                    amplitude_seq(i,j,1:length(temp)) = temp;
                elseif (d_info.group(i)==8 && (j==4 || j==5))
                    amplitude_seq(i,j,1:length(temp)) = temp;
                elseif (d_info.group(i)==7 && (j==4 || j==5))
                    amplitude_ctrl(i,j,1:length(temp)) = temp;
                elseif (d_info.group(i)==8 && (j==2 || j==3))
                    amplitude_ctrl(i,j,1:length(temp)) = temp;
                end
            catch
            end
            amplitude_unpaired(i,j,:) = amplitude(i,j,:);
        end
    end
end
amplitude_seq_unpaired_sd1 = squeeze(nanmean(amplitude_seq(:,[2,4],:),2));
amplitude_seq_paired_sd1 = squeeze(nanmean(amplitude_seq(:,[2,4],:),2));
amplitude_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
amplitude_ctrl_unpaired_sd1 = squeeze(nanmean(amplitude_ctrl(:,[2,4],:),2));
amplitude_ctrl_paired_sd1 = squeeze(nanmean(amplitude_ctrl(:,[2,4],:),2));
amplitude_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
amplitude_seq_paired_sd12 = cat(3,squeeze(nanmean(amplitude_seq(:,[2,4],:),2)),squeeze(nanmean(amplitude_seq(:,[3,5],:),2)));
amplitude_seq_paired_sd12(~these_animals_paired_bothDays,:,:) = NaN;
amplitude_ctrl_paired_sd12 = cat(3,squeeze(nanmean(amplitude_ctrl(:,[2,4],:),2)),squeeze(nanmean(amplitude_ctrl(:,[3,5],:),2)));
amplitude_ctrl_paired_sd12(~these_animals_paired_bothDays,:,:) = NaN;

resolution = nan(d_info.numAnimals,d_info.numDays,20);
resolution_seq = nan(d_info.numAnimals,d_info.numDays,20);
resolution_ctrl = nan(d_info.numAnimals,d_info.numDays,20);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp1 = d{i,j}.(this_struct).responders_main.respAll(:);
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
end
resolution_seq_paired_sd1 = squeeze(nanmean(resolution_seq(:,[2,4],:),2));
resolution_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
resolution_ctrl_paired_sd1 = squeeze(nanmean(resolution_ctrl(:,[2,4],:),2));
resolution_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;

amplitudeResolution = nan(d_info.numAnimals,d_info.numDays,20);
amplitudeResolution_seq = nan(d_info.numAnimals,d_info.numDays,20);
amplitudeResolution_ctrl = nan(d_info.numAnimals,d_info.numDays,20);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp1 = d{i,j}.resp.avgAct_resp(:);
                temp2 = discretize(d{i,j}.resp.dist_closestLaser(:),[0:2.5:50]);
                temp = nan(nanmax(temp2),1);
                for k=1:nanmax(temp2)
                    temp(k) = nanmean(temp1(find(temp2==k)));
                end
                resolution(i,j,1:length(temp)) = temp;
                if d_info.group(i)==7 && (j==2 || j==3)
                    amplitudeResolution_seq(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    amplitudeResolution_seq(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    amplitudeResolution_ctrl(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    amplitudeResolution_ctrl(i,j,1:length(temp)) = temp;
                end
            catch
            end
        end
    end
end
amplitudeResolution_seq_paired_sd1 = squeeze(nanmean(amplitudeResolution_seq(:,[2,4],:),2));
amplitudeResolution_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
amplitudeResolution_ctrl_paired_sd1 = squeeze(nanmean(amplitudeResolution_ctrl(:,[2,4],:),2));
amplitudeResolution_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;


%% --- Main figure ---

%% Fig4_ResponseAnalysis_Responders

these_labels = categorical({'Consistent','Shuffled'});
these_labels = reordercats(these_labels,{'Consistent','Shuffled'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq = responders_seq_paired_sd1;
this_data_ctrl = responders_ctrl_paired_sd1;
v = bar(these_labels,nanmean([this_data_seq,this_data_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.ctrl; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq)
    plot(these_labels,[this_data_seq(i),this_data_ctrl(i)],'-k','LineWidth',1)
end
ylim([0,30])
yticks([0,10,20,30])
ylabel({'Photoactivated cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_Responders.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_Responders.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_Responders.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_Responders.txt'],'wt');
[temp1,~,temp2] = signrank(this_data_seq,this_data_ctrl);
fprintf(fid,['\nResponseAnalysis_Responders\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nmean(seq)=',num2str(nanmean(this_data_seq),4),'\nmean(ctrl)=',num2str(nanmean(this_data_ctrl),4),...
    '\nstd(seq)=',num2str(nanstd(this_data_seq),4),'\nstd(ctrl)=',num2str(nanstd(this_data_ctrl),4),...
    '\nsem(seq)=',num2str(nansem(this_data_seq,1),4),'\nsem(ctrl)=',num2str(nansem(this_data_ctrl,1),4),...
    '\nsignrank p=',num2str(temp1,4),...
    '\nsignrank T=',num2str(temp2.signedrank,4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_Responders\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nmean(seq)=',num2str(nanmean(this_data_seq),4),'\nmean(ctrl)=',num2str(nanmean(this_data_ctrl),4),...
    '\nstd(seq)=',num2str(nanstd(this_data_seq),4),'\nstd(ctrl)=',num2str(nanstd(this_data_ctrl),4),...
    '\nsem(seq)=',num2str(nansem(this_data_seq,1),4),'\nsem(ctrl)=',num2str(nansem(this_data_ctrl,1),4),...
    '\nsignrank p=',num2str(temp1,4),...
    '\nsignrank T=',num2str(temp2.signedrank,4),'\n']);


%% Fig4_ResponseAnalysis_Amplitude

these_labels = categorical({'Consistent','Shuffled'});
these_labels = reordercats(these_labels,{'Consistent','Shuffled'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq = nanmean(amplitude_seq_paired_sd1,2); %amplitude_seq_paired_sd1(:,1);
this_data_ctrl = nanmean(amplitude_ctrl_paired_sd1,2); %amplitude_ctrl_paired_sd1(:,1);
v = bar(these_labels,nanmean([this_data_seq,this_data_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.ctrl; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq)
    plot(these_labels,[this_data_seq(i),this_data_ctrl(i)],'-k','LineWidth',1)
end
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Photoactivation amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_Amplitude.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_Amplitude.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_Amplitude.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_Amplitude.txt'],'wt');
[temp1,~,temp2] = signrank(this_data_seq,this_data_ctrl);
fprintf(fid,['\nResponseAnalysis_Amplitude\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nmean(seq)=',num2str(nanmean(this_data_seq),4),'\nmean(ctrl)=',num2str(nanmean(this_data_ctrl),4),...
    '\nstd(seq)=',num2str(nanstd(this_data_seq),4),'\nstd(ctrl)=',num2str(nanstd(this_data_ctrl),4),...
    '\nsem(seq)=',num2str(nansem(this_data_seq,1),4),'\nsem(ctrl)=',num2str(nansem(this_data_ctrl,1),4),...
    '\nsignrank p=',num2str(temp1,4),...
    '\nsignrank T=',num2str(temp2.signedrank,4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_Amplitude\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nmean(seq)=',num2str(nanmean(this_data_seq),4),'\nmean(ctrl)=',num2str(nanmean(this_data_ctrl),4),...
    '\nstd(seq)=',num2str(nanstd(this_data_seq),4),'\nstd(ctrl)=',num2str(nanstd(this_data_ctrl),4),...
    '\nsem(seq)=',num2str(nansem(this_data_seq,1),4),'\nsem(ctrl)=',num2str(nansem(this_data_ctrl,1),4),...
    '\nsignrank p=',num2str(temp1,4),...
    '\nsignrank T=',num2str(temp2.signedrank,4),'\n']);


%% --- Supplementary figure ---

%% Fig4_ResponseAnalysis_Resolution_ind

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq = resolution_seq_paired_sd1;
this_data_ctrl = resolution_ctrl_paired_sd1;
plot(this_data_ctrl'*100,'Color',mean([p.col.ctrl;p.col.white]),'LineWidth',0.5)
plot(this_data_seq'*100,'Color',mean([p.col.seq;p.col.white]),'LineWidth',0.5)
plot(nanmean(this_data_ctrl',2)'*100,'Color',p.col.ctrl,'LineWidth',1.5)
plot(nanmean(this_data_seq',2)*100,'Color',p.col.seq,'LineWidth',1.5)
xlim([1,size(this_data_seq,2)+1]) % xlim([1,size(this_data_seq,2)+1]-0.5)
xticks([0:4:length([0:2.5:50])]+1) % xticks([0:4:length([0:2.5:50])]+0.5)
xticklabels({'0','10','20','30','40','50'})
yticks([0,50,100])
ylim([0,100])
ytickformat('percentage')
xlabel({'Distance from closest';'laser beamlet (Um)'})
ylabel('Photoactivation probability')

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_Resolution_ind.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_Resolution_ind.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_Resolution_ind.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_Resolution_ind.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_Resolution\nn(seq)=',num2str(size(rmmissing(this_data_seq),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_ctrl),1)),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_Resolution\nn(seq)=',num2str(size(rmmissing(this_data_seq),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_ctrl),1)),'\n']);

% HWHM
% this_data = this_data_ctrl;
% this_max = this_data(:,1);
% this_hwhm = nan(size(this_max));
% for i=1:length(this_max)
%     if ~isnan(this_max(i))
%         [~,this_hwhm(i)] = nanmin(abs(this_data(i,:)-0.5*this_max(i)));
%     end
% end
% nanmean(this_hwhm)
% nanstd(this_hwhm)

%% Fig4_ResponseAnalysis_Resolution

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq = resolution_seq_paired_sd1;
this_data_ctrl = resolution_ctrl_paired_sd1;
shadedErrorBar(1:size(this_data_ctrl,2),nanmean(this_data_ctrl'*100,2),nansem(this_data_ctrl'*100,2),'lineProps',p.col.ctrl);
shadedErrorBar(1:size(this_data_seq,2),nanmean(this_data_seq'*100,2),nansem(this_data_seq'*100,2),'lineProps',p.col.seq);
xlim([1,size(this_data_seq,2)+1]) % xlim([1,size(this_data_seq,2)+1]-0.5)
xticks([0:4:length([0:2.5:50])]+1) % xticks([0:4:length([0:2.5:50])]+0.5)
xticklabels({'0','10','20','30','40','50'})
yticks([0,50,100])
ylim([0,100])
ytickformat('percentage')
xlabel({'Distance from closest';'laser beamlet (Um)'})
ylabel('Photoactivation probability')

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_Resolution.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_Resolution.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_Resolution.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_Resolution.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_Resolution\nn(seq)=',num2str(size(rmmissing(this_data_seq),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_ctrl),1)),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_Resolution\nn(seq)=',num2str(size(rmmissing(this_data_seq),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_ctrl),1)),'\n']);


%% Fig4_ResponseAnalysis_Targets

these_labels = categorical({'Consistent','Shuffled'});
these_labels = reordercats(these_labels,{'Consistent','Shuffled'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq = targets_seq_paired_sd1;
this_data_ctrl = targets_ctrl_paired_sd1;
v = bar(these_labels,nanmean([this_data_seq,this_data_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.ctrl; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq)
    plot(these_labels,[this_data_seq(i),this_data_ctrl(i)],'-k','LineWidth',1)
end
yline(6,'k:','LineWidth',1);
ylim([0,6])
yticks([0,2,4,6])
ylabel({'Photoactivated';'targeted cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_Targets.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_Targets.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_Targets.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_Targets.txt'],'wt');
[temp1,~,temp2] = signrank(this_data_seq,this_data_ctrl);
fprintf(fid,['\nResponseAnalysis_Targets\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nmean(seq)=',num2str(nanmean(this_data_seq),4),'\nmean(ctrl)=',num2str(nanmean(this_data_ctrl),4),...
    '\nstd(seq)=',num2str(nanstd(this_data_seq),4),'\nstd(ctrl)=',num2str(nanstd(this_data_ctrl),4),...
    '\nsem(seq)=',num2str(nansem(this_data_seq,1),4),'\nsem(ctrl)=',num2str(nansem(this_data_ctrl,1),4),...
    '\nsignrank p=',num2str(temp1,4),...
    '\nsignrank T=',num2str(temp2.signedrank,4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_Targets\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nmean(seq)=',num2str(nanmean(this_data_seq),4),'\nmean(ctrl)=',num2str(nanmean(this_data_ctrl),4),...
    '\nstd(seq)=',num2str(nanstd(this_data_seq),4),'\nstd(ctrl)=',num2str(nanstd(this_data_ctrl),4),...
    '\nsem(seq)=',num2str(nansem(this_data_seq,1),4),'\nsem(ctrl)=',num2str(nansem(this_data_ctrl,1),4),...
    '\nsignrank p=',num2str(temp1,4),...
    '\nsignrank T=',num2str(temp2.signedrank,4),'\n']);


%% Fig4_ResponseAnalysis_OffTargetsPerTarget

these_labels = categorical({'Consistent','Shuffled'});
these_labels = reordercats(these_labels,{'Consistent','Shuffled'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq = offTargetsPerTarget_seq_paired_sd1;
this_data_ctrl = offTargetsPerTarget_ctrl_paired_sd1;
v = bar(these_labels,nanmean([this_data_seq,this_data_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.ctrl; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq)
    plot(these_labels,[this_data_seq(i),this_data_ctrl(i)],'-k','LineWidth',1)
end
ylim([0,6])
yticks([0,2,4,6])
ylabel({'Photoactivated non-targeted cells';'per photoactivated target cell'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_OffTargetsPerTarget.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_OffTargetsPerTarget.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_OffTargetsPerTarget.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_OffTargetsPerTarget.txt'],'wt');
[temp1,~,temp2] = signrank(this_data_seq,this_data_ctrl);
fprintf(fid,['\nResponseAnalysis_OffTargetsPerTarget\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nmean(seq)=',num2str(nanmean(this_data_seq),4),'\nmean(ctrl)=',num2str(nanmean(this_data_ctrl),4),...
    '\nstd(seq)=',num2str(nanstd(this_data_seq),4),'\nstd(ctrl)=',num2str(nanstd(this_data_ctrl),4),...
    '\nsem(seq)=',num2str(nansem(this_data_seq,1),4),'\nsem(ctrl)=',num2str(nansem(this_data_ctrl,1),4),...
    '\nsignrank p=',num2str(temp1,4),...
    '\nsignrank T=',num2str(temp2.signedrank,4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_OffTargetsPerTarget\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nmean(seq)=',num2str(nanmean(this_data_seq),4),'\nmean(ctrl)=',num2str(nanmean(this_data_ctrl),4),...
    '\nstd(seq)=',num2str(nanstd(this_data_seq),4),'\nstd(ctrl)=',num2str(nanstd(this_data_ctrl),4),...
    '\nsem(seq)=',num2str(nansem(this_data_seq,1),4),'\nsem(ctrl)=',num2str(nansem(this_data_ctrl,1),4),...
    '\nsignrank p=',num2str(temp1,4),...
    '\nsignrank T=',num2str(temp2.signedrank,4),'\n']);


%% Fig4_ResponseAnalysis_Targets_2_OffTargets

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq_x = targets_seq_paired_sd1;
this_data_ctrl_x = targets_ctrl_paired_sd1;
this_data_seq_y = offTargets_seq_paired_sd1;
this_data_ctrl_y = offTargets_ctrl_paired_sd1;
scatter(this_data_ctrl_x,this_data_ctrl_y,'.','SizeData',100,'MarkerEdgeColor',p.col.ctrl)
scatter(this_data_seq_x,this_data_seq_y,'.','SizeData',100,'MarkerEdgeColor',p.col.seq)
[this_corr_r,this_corr_p] = fitLine([this_data_seq_x;this_data_ctrl_x],[this_data_seq_y;this_data_ctrl_y],p.col.black);
xlim([3,6])
xticks([3,4,5,6])
xlabel({'Photoactivated';'targeted cells per cluster'})
ylim([0,20])
yticks([0,10,20])
ylabel({'Photoactivated';'non-targeted cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_Targets_2_OffTargets.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_Targets_2_OffTargets.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_Targets_2_OffTargets.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_Targets_2_OffTargets.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_Targets_2_OffTargets\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_Targets_2_OffTargets\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_ResponseAnalysis_Responders_2_Amplitude

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq_x = responders_seq_paired_sd1;
this_data_ctrl_x = responders_ctrl_paired_sd1;
this_data_seq_y = nanmean(amplitude_seq_paired_sd1,2); %amplitude_seq_paired_sd1(:,1);
this_data_ctrl_y = nanmean(amplitude_ctrl_paired_sd1,2); %amplitude_ctrl_paired_sd1(:,1);
scatter(this_data_ctrl_x,this_data_ctrl_y,'.','SizeData',100,'MarkerEdgeColor',p.col.ctrl)
scatter(this_data_seq_x,this_data_seq_y,'.','SizeData',100,'MarkerEdgeColor',p.col.seq)
[this_corr_r,this_corr_p] = fitLine([this_data_seq_x;this_data_ctrl_x],[this_data_seq_y;this_data_ctrl_y],p.col.black);
xlim([0,30])
xticks([0,10,20,30])
xlabel({'Photoactivated cells per cluster'})
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Photoactivation amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_Responders_2_Amplitude.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_Responders_2_Amplitude.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_Responders_2_Amplitude.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_Responders_2_Amplitude.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_Responders_2_Amplitude\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_Responders_2_Amplitude\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_ResponseAnalysis_StabilityWithinDay

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq = amplitude_seq_paired_sd1;
this_data_ctrl = amplitude_ctrl_paired_sd1;
shadedErrorBar(1:size(this_data_ctrl,2),nanmean(this_data_ctrl',2),nansem(this_data_ctrl',2),'lineProps',p.col.ctrl);
shadedErrorBar(1:size(this_data_seq,2),nanmean(this_data_seq',2),nansem(this_data_seq',2),'lineProps',p.col.seq);
xlim([1,size(this_data_seq,2)])
xticks([2:2:8])
xticklabels({'2','4','6','8'})
yticks([0,1,2,3])
ylim([0,3])
xlabel('Trial block')
ylabel({'Photoactivation amplitude (S)'})

temp1 = [this_data_seq;this_data_ctrl];
temp2 = repmat(1:8,size(temp1,1),1);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete');
savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_StabilityWithinDay.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_StabilityWithinDay.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_StabilityWithinDay.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_StabilityWithinDay.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_StabilityWithinDay\nn(seq)=',num2str(size(rmmissing(this_data_seq),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_ctrl),1)),...
	'\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_StabilityWithinDay\nn(seq)=',num2str(size(rmmissing(this_data_seq),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_ctrl),1)),...
	'\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_ResponseAnalysis_StabilityWithinDay_ind

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq = amplitude_seq_paired_sd1;
this_data_ctrl = amplitude_ctrl_paired_sd1;
plot(this_data_ctrl','Color',mean([p.col.ctrl;p.col.white]),'LineWidth',0.5)
plot(this_data_seq','Color',mean([p.col.seq;p.col.white]),'LineWidth',0.5)
plot(nanmean(this_data_ctrl',2)','Color',p.col.ctrl,'LineWidth',1.5)
plot(nanmean(this_data_seq',2),'Color',p.col.seq,'LineWidth',1.5)
xlim([1,size(this_data_seq,2)])
xticks([2:2:8])
xticklabels({'2','4','6','8'})
yticks([0,1,2,3])
ylim([0,3])
xlabel('Trial block')
ylabel({'Photoactivation amplitude (S)'})

% [temp1n,temp3n,temp2n] = kruskalwallis(this_data_seq); close;
% [c,m,h] = multcompare(temp2n,'CType','dunn-sidak');
% temp = rmmissing(rmmissing(this_data_seq,1,'MinNumMissing',4),2,'MinNumMissing',4);
% temp0 = repmat(1:size(temp,2),size(temp,1),1)
% statsp=kwtest([temp(:),temp0(:)]);
% % dunn(temp(:)',temp0(:)');
% 
% [temp1n,temp3n,temp2n] = kruskalwallis(this_data_ctrl); close;
% [c,m,h] = multcompare(temp2n,'CType','dunn-sidak');
% temp = rmmissing(rmmissing(this_data_ctrl,1,'MinNumMissing',4),2,'MinNumMissing',4);
% temp0 = repmat(1:size(temp,2),size(temp,1),1)
% statsp=kwtest([temp(:),temp0(:)]);
% % dunn(temp(:)',temp0(:)');

% this_data_seq_ = this_data_seq(:,8)./this_data_seq(:,1);
% this_data_ctrl_ = this_data_ctrl(:,8)./this_data_ctrl(:,1);
% [temp1,~,temp2] = signrank(this_data_seq_,this_data_ctrl_);

temp1 = [this_data_seq;this_data_ctrl];
temp2 = repmat(1:8,size(temp1,1),1);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete');
savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_StabilityWithinDay_ind.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_StabilityWithinDay_ind.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_StabilityWithinDay_ind.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_StabilityWithinDay_ind.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_StabilityWithinDay\nn(seq)=',num2str(size(rmmissing(this_data_seq),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_ctrl),1)),...
	'\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_StabilityWithinDay\nn(seq)=',num2str(size(rmmissing(this_data_seq),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_ctrl),1)),...
	'\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_ResponseAnalysis_RespondersAcrossDays_unpaired

these_labels = categorical({'Day 2','Day 3','Day 4','Day 5'});
these_labels = reordercats(these_labels,{'Day 2','Day 3','Day 4','Day 5'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = responders_unpaired;
v = bar(these_labels,nanmean([this_data(:,2),this_data(:,3),this_data(:,4),this_data(:,5)],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.darkGray; v.CData(3,:) = p.col.darkGray; v.CData(4,:) = p.col.darkGray; v.EdgeColor = 'none'; 
plot(these_labels,this_data(:,2:5),'-','Color',p.col.black,'LineWidth',1)
for i=1:d_info.numAnimals
    for j=2:5
        if ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
            plot(these_labels(j-1),this_data(i,j),'.','Color',p.col.seq)
        elseif ((d_info.group(i)==8 && (j==2 || j==3)) || (d_info.group(i)==7 && (j==4 || j==5)))
            plot(these_labels(j-1),this_data(i,j),'.','Color',p.col.ctrl)
        end
    end
end
ylim([0,30])
yticks([0,10,20,30])
ylabel('Photoactivated cells per cluster')

[temp1n,temp3n,temp2n] = kruskalwallis(this_data); close;
[c,m,h] = multcompare(temp2n,'CType','dunn-sidak');
temp = rmmissing(rmmissing(this_data,1,'MinNumMissing',4),2,'MinNumMissing',4);
temp0 = repmat(1:size(temp,2),size(temp,1),1)
statsp=kwtest([temp(:),temp0(:)]);
dunn(temp(:)',temp0(:)');

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_RespondersAcrossDays_unpaired.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_RespondersAcrossDays_unpaired.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_RespondersAcrossDays_unpaired.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_RespondersAcrossDays_unpaired.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_RespondersAcrossDays_unpaired\nn(day2)=',num2str(length(rmmissing(this_data(:,2)))),'\nn(dat3)=',num2str(length(rmmissing(this_data(:,3)))),'\nn(day4)=',num2str(length(rmmissing(this_data(:,4)))),'\nn(day5)=',num2str(length(rmmissing(this_data(:,5)))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_RespondersAcrossDays_unpaired\nn(day2)=',num2str(length(rmmissing(this_data(:,2)))),'\nn(dat3)=',num2str(length(rmmissing(this_data(:,3)))),'\nn(day4)=',num2str(length(rmmissing(this_data(:,4)))),'\nn(day5)=',num2str(length(rmmissing(this_data(:,5)))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig4_ResponseAnalysis_AmplitudeAcrossDays_unpaired

these_labels = categorical({'Day 2','Day 3','Day 4','Day 5'});
these_labels = reordercats(these_labels,{'Day 2','Day 3','Day 4','Day 5'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nanmean(amplitude_unpaired,3);
v = bar(these_labels,nanmean([this_data(:,2),this_data(:,3),this_data(:,4),this_data(:,5)],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.darkGray; v.CData(3,:) = p.col.darkGray; v.CData(4,:) = p.col.darkGray; v.EdgeColor = 'none'; 
plot(these_labels,this_data(:,2:5),'-','Color',p.col.black,'LineWidth',1)
for i=1:d_info.numAnimals
    for j=2:5
        if ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
            plot(these_labels(j-1),this_data(i,j),'.','Color',p.col.seq)
        elseif ((d_info.group(i)==8 && (j==2 || j==3)) || (d_info.group(i)==7 && (j==4 || j==5)))
            plot(these_labels(j-1),this_data(i,j),'.','Color',p.col.ctrl)
        end
    end
end
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Photoactivation amplitude (S)'})

[temp1n,temp3n,temp2n] = kruskalwallis(this_data); close;
[c,m,h] = multcompare(temp2n,'CType','dunn-sidak');
temp = rmmissing(rmmissing(this_data,1,'MinNumMissing',4),2,'MinNumMissing',4);
temp0 = repmat(1:size(temp,2),size(temp,1),1)
statsp=kwtest([temp(:),temp0(:)]);
% dunn(temp(:)',temp0(:)');

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_AmplitudeAcrossDays_unpaired.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_AmplitudeAcrossDays_unpaired.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_AmplitudeAcrossDays_unpaired.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_AmplitudeAcrossDays_unpaired.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_AmplitudeAcrossDays_unpaired\nn(day2)=',num2str(length(rmmissing(this_data(:,2)))),'\nn(dat3)=',num2str(length(rmmissing(this_data(:,3)))),'\nn(day4)=',num2str(length(rmmissing(this_data(:,4)))),'\nn(day5)=',num2str(length(rmmissing(this_data(:,5)))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_AmplitudeAcrossDays_unpaired\nn(day2)=',num2str(length(rmmissing(this_data(:,2)))),'\nn(dat3)=',num2str(length(rmmissing(this_data(:,3)))),'\nn(day4)=',num2str(length(rmmissing(this_data(:,4)))),'\nn(day5)=',num2str(length(rmmissing(this_data(:,5)))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% --- Reserve ---

%% Fig4_ResponseAnalysis_OffTargets

these_labels = categorical({'Consistent','Shuffled'});
these_labels = reordercats(these_labels,{'Consistent','Shuffled'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq = offTargets_seq_paired_sd1;
this_data_ctrl = offTargets_ctrl_paired_sd1;
v = bar(these_labels,nanmean([this_data_seq,this_data_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.ctrl; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq)
    plot(these_labels,[this_data_seq(i),this_data_ctrl(i)],'-k','LineWidth',1)
end
ylim([0,20])
yticks([0,10,20])
ylabel({'Photoactivated';'non-targeted cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_OffTargets.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_OffTargets.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_OffTargets.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_OffTargets.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_OffTargets\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nmean(seq)=',num2str(nanmean(this_data_seq),4),'\nmean(ctrl)=',num2str(nanmean(this_data_ctrl),4),...
    '\nstd(seq)=',num2str(nanstd(this_data_seq),4),'\nstd(ctrl)=',num2str(nanstd(this_data_ctrl),4),...
    '\nsem(seq)=',num2str(nansem(this_data_seq,1),4),'\nsem(ctrl)=',num2str(nansem(this_data_ctrl,1),4),...
    '\nsignrank p=',num2str(signrank(this_data_seq,this_data_ctrl),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_OffTargets\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nmean(seq)=',num2str(nanmean(this_data_seq),4),'\nmean(ctrl)=',num2str(nanmean(this_data_ctrl),4),...
    '\nstd(seq)=',num2str(nanstd(this_data_seq),4),'\nstd(ctrl)=',num2str(nanstd(this_data_ctrl),4),...
    '\nsem(seq)=',num2str(nansem(this_data_seq,1),4),'\nsem(ctrl)=',num2str(nansem(this_data_ctrl,1),4),...
    '\nsignrank p=',num2str(signrank(this_data_seq,this_data_ctrl),4),'\n']);


%% Fig4_ResponseAnalysis_Targets_2_OffTargetsPerTarget

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq_x = targets_seq_paired_sd1;
this_data_ctrl_x = targets_ctrl_paired_sd1;
this_data_seq_y = offTargetsPerTarget_seq_paired_sd1;
this_data_ctrl_y = offTargetsPerTarget_ctrl_paired_sd1;
scatter(this_data_ctrl_x,this_data_ctrl_y,'.','SizeData',100,'MarkerEdgeColor',p.col.ctrl)
scatter(this_data_seq_x,this_data_seq_y,'.','SizeData',100,'MarkerEdgeColor',p.col.seq)
[this_corr_r,this_corr_p] = fitLine([this_data_seq_x;this_data_ctrl_x],[this_data_seq_y;this_data_ctrl_y],p.col.black);
xlim([3,6])
xticks([3,4,5,6])
xlabel({'Photoactivated';'targeted cells per cluster'})
ylim([0,4])
yticks([0,2,4])
ylabel({'Photoactivated non-targeted cells';'per photoactivated target cell'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_Targets_2_OffTargetsPerTarget.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_Targets_2_OffTargetsPerTarget.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_Targets_2_OffTargetsPerTarget.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_Targets_2_OffTargetsPerTarget.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_Targets_2_OffTargetsPerTarget\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_Targets_2_OffTargetsPerTarget\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_ResponseAnalysis_RespondersAcrossDays

these_labels = categorical({'1st','2nd','1s+','2n+'});
these_labels = reordercats(these_labels,{'1st','2nd','1s+','2n+'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq_1 = responders_seq_paired_sd12(:,1);
this_data_seq_2 = responders_seq_paired_sd12(:,2);
this_data_ctrl_1 = responders_ctrl_paired_sd12(:,1);
this_data_ctrl_2 = responders_ctrl_paired_sd12(:,2);
v = bar(these_labels,nanmean([this_data_seq_1,this_data_seq_2,this_data_ctrl_1,this_data_ctrl_2],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq_1)
    plot(these_labels,[this_data_seq_1(i),this_data_seq_2(i),this_data_ctrl_1(i),this_data_ctrl_2(i)],'-k','LineWidth',1)
end
ylim([0,30])
yticks([0,10,20,30])
ylabel({'Photoactivated cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_RespondersAcrossDays.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_RespondersAcrossDays.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_RespondersAcrossDays.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_RespondersAcrossDays.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_RespondersAcrossDays\nn(seq1)=',num2str(length(rmmissing(this_data_seq_1))),'\nn(seq2)=',num2str(length(rmmissing(this_data_seq_2))),'\nn(ctrl1)=',num2str(length(rmmissing(this_data_ctrl_1))),'\nn(ctrl2)=',num2str(length(rmmissing(this_data_ctrl_2))),...
    '\nsignrank seq p=',num2str(signrank(this_data_seq_1,this_data_seq_2),4),'\nsignrank ctrl p=',num2str(signrank(this_data_ctrl_1,this_data_ctrl_2),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_RespondersAcrossDays\nn(seq1)=',num2str(length(rmmissing(this_data_seq_1))),'\nn(seq2)=',num2str(length(rmmissing(this_data_seq_2))),'\nn(ctrl1)=',num2str(length(rmmissing(this_data_ctrl_1))),'\nn(ctrl2)=',num2str(length(rmmissing(this_data_ctrl_2))),...
    '\nsignrank seq p=',num2str(signrank(this_data_seq_1,this_data_seq_2),4),'\nsignrank ctrl p=',num2str(signrank(this_data_ctrl_1,this_data_ctrl_2),4),'\n']);


%% Fig4_ResponseAnalysis_StabilityAcrossDays

these_labels = categorical({'1st','2nd','1s+','2n+'});
these_labels = reordercats(these_labels,{'1st','2nd','1s+','2n+'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq_1 = nanmean(amplitude_seq_paired_sd12(:,:,1),2); %amplitude_seq_paired_sd12(:,1,1);
this_data_seq_2 = nanmean(amplitude_seq_paired_sd12(:,:,2),2); %amplitude_seq_paired_sd12(:,1,2);
this_data_ctrl_1 = nanmean(amplitude_ctrl_paired_sd12(:,:,1),2); %amplitude_ctrl_paired_sd12(:,1,1);
this_data_ctrl_2 = nanmean(amplitude_ctrl_paired_sd12(:,:,2),2); %amplitude_ctrl_paired_sd12(:,1,2);
v = bar(these_labels,nanmean([this_data_seq_1,this_data_seq_2,this_data_ctrl_1,this_data_ctrl_2],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq_1)
    plot(these_labels,[this_data_seq_1(i),this_data_seq_2(i),this_data_ctrl_1(i),this_data_ctrl_2(i)],'-k','LineWidth',1)
end
ylim([0,3])
yticks([0,1,2,3])
ylabel('Photoactivation amplitude (S)')

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_StabilityAcrossDays.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_StabilityAcrossDays.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_StabilityAcrossDays.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_StabilityAcrossDays.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_StabilityAcrossDays\nn(seq1)=',num2str(length(rmmissing(this_data_seq_1))),'\nn(seq2)=',num2str(length(rmmissing(this_data_seq_2))),'\nn(ctrl1)=',num2str(length(rmmissing(this_data_ctrl_1))),'\nn(ctrl2)=',num2str(length(rmmissing(this_data_ctrl_2))),...
    '\nsignrank seq p=',num2str(signrank(this_data_seq_1,this_data_seq_2),4),'\nsignrank ctrl p=',num2str(signrank(this_data_ctrl_1,this_data_ctrl_2),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_StabilityAcrossDays\nn(seq1)=',num2str(length(rmmissing(this_data_seq_1))),'\nn(seq2)=',num2str(length(rmmissing(this_data_seq_2))),'\nn(ctrl1)=',num2str(length(rmmissing(this_data_ctrl_1))),'\nn(ctrl2)=',num2str(length(rmmissing(this_data_ctrl_2))),...
    '\nsignrank seq p=',num2str(signrank(this_data_seq_1,this_data_seq_2),4),'\nsignrank ctrl p=',num2str(signrank(this_data_ctrl_1,this_data_ctrl_2),4),'\n']);


%% Fig4_ResponseAnalysis_RespondersByRunning

these_labels = categorical({'NR','R','nr','r'});
these_labels = reordercats(these_labels,{'NR','R','nr','r'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq_nonrunner = responders_seq_unpaired_sd1(running_seq_unpaired_sd1==0);
this_data_seq_runner = responders_seq_unpaired_sd1(running_seq_unpaired_sd1==1);
this_data_ctrl_nonrunner = responders_ctrl_unpaired_sd1(running_ctrl_unpaired_sd1==0);
this_data_ctrl_runner = responders_ctrl_unpaired_sd1(running_ctrl_unpaired_sd1==1);
this_data = nan(10,4);
this_data(1:length(this_data_seq_nonrunner),1) = this_data_seq_nonrunner;
this_data(1:length(this_data_seq_runner),2) = this_data_seq_runner;
this_data(1:length(this_data_ctrl_nonrunner),3) = this_data_ctrl_nonrunner;
this_data(1:length(this_data_ctrl_runner),4) = this_data_ctrl_runner;
v = bar(these_labels,[nanmean(this_data_seq_nonrunner),nanmean(this_data_seq_runner),nanmean(this_data_ctrl_nonrunner),nanmean(this_data_ctrl_runner)]);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.nonrunner; v.CData(2,:) = p.col.runner; v.CData(3,:) = p.col.nonrunner; v.CData(4,:) = p.col.runner; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq_nonrunner)
    plot(these_labels(1),this_data_seq_nonrunner,'.','Color',p.col.seq)
end
for i=1:length(this_data_seq_runner)
    plot(these_labels(2),this_data_seq_runner,'.','Color',p.col.seq)
end
for i=1:length(this_data_ctrl_nonrunner)
    plot(these_labels(3),this_data_ctrl_nonrunner,'.','Color',p.col.ctrl)
end
for i=1:length(this_data_ctrl_runner)
    plot(these_labels(4),this_data_ctrl_runner,'.','Color',p.col.ctrl)
end
ylim([0,30])
yticks([0,10,20,30])
ylabel({'Photoactivated cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_RespondersByRunning.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_RespondersByRunning.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_RespondersByRunning.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_RespondersByRunning.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_RespondersByRunning\nn(seq_nonrunner)=',num2str(length(rmmissing(this_data_seq_nonrunner))),'\nn(seq_runner)=',num2str(length(rmmissing(this_data_seq_runner))),'\nn(ctrl_nonrunner)=',num2str(length(rmmissing(this_data_ctrl_nonrunner))),'\nn(ctrl_runner)=',num2str(length(rmmissing(this_data_ctrl_runner))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_RespondersByRunning\nn(seq_nonrunner)=',num2str(length(rmmissing(this_data_seq_nonrunner))),'\nn(seq_runner)=',num2str(length(rmmissing(this_data_seq_runner))),'\nn(ctrl_nonrunner)=',num2str(length(rmmissing(this_data_ctrl_nonrunner))),'\nn(ctrl_runner)=',num2str(length(rmmissing(this_data_ctrl_runner))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig4_ResponseAnalysis_AmplitudeByRunning

these_labels = categorical({'NR','R','nr','r'});
these_labels = reordercats(these_labels,{'NR','R','nr','r'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq_nonrunner = nanmean(amplitude_seq_unpaired_sd1(running_seq_unpaired_sd1==0,:),2);
this_data_seq_runner = nanmean(amplitude_seq_unpaired_sd1(running_seq_unpaired_sd1==1,:),2);
this_data_ctrl_nonrunner = nanmean(amplitude_ctrl_unpaired_sd1(running_ctrl_unpaired_sd1==0,:),2);
this_data_ctrl_runner = nanmean(amplitude_ctrl_unpaired_sd1(running_ctrl_unpaired_sd1==1,:),2);
this_data = nan(10,4);
this_data(1:length(this_data_seq_nonrunner),1) = this_data_seq_nonrunner;
this_data(1:length(this_data_seq_runner),2) = this_data_seq_runner;
this_data(1:length(this_data_ctrl_nonrunner),3) = this_data_ctrl_nonrunner;
this_data(1:length(this_data_ctrl_runner),4) = this_data_ctrl_runner;
v = bar(these_labels,[nanmean(this_data_seq_nonrunner),nanmean(this_data_seq_runner),nanmean(this_data_ctrl_nonrunner),nanmean(this_data_ctrl_runner)]);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.nonrunner; v.CData(2,:) = p.col.runner; v.CData(3,:) = p.col.nonrunner; v.CData(4,:) = p.col.runner; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq_nonrunner)
    plot(these_labels(1),this_data_seq_nonrunner,'.','Color',p.col.seq)
end
for i=1:length(this_data_seq_runner)
    plot(these_labels(2),this_data_seq_runner,'.','Color',p.col.seq)
end
for i=1:length(this_data_ctrl_nonrunner)
    plot(these_labels(3),this_data_ctrl_nonrunner,'.','Color',p.col.ctrl)
end
for i=1:length(this_data_ctrl_runner)
    plot(these_labels(4),this_data_ctrl_runner,'.','Color',p.col.ctrl)
end
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Photoactivation amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_AmplitudeByRunning.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_AmplitudeByRunning.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_AmplitudeByRunning.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_AmplitudeByRunning.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_AmplitudeByRunning\nn(seq_nonrunner)=',num2str(length(rmmissing(this_data_seq_nonrunner))),'\nn(seq_runner)=',num2str(length(rmmissing(this_data_seq_runner))),'\nn(ctrl_nonrunner)=',num2str(length(rmmissing(this_data_ctrl_nonrunner))),'\nn(ctrl_runner)=',num2str(length(rmmissing(this_data_ctrl_runner))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_AmplitudeByRunning\nn(seq_nonrunner)=',num2str(length(rmmissing(this_data_seq_nonrunner))),'\nn(seq_runner)=',num2str(length(rmmissing(this_data_seq_runner))),'\nn(ctrl_nonrunner)=',num2str(length(rmmissing(this_data_ctrl_nonrunner))),'\nn(ctrl_runner)=',num2str(length(rmmissing(this_data_ctrl_runner))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig4_ResponseAnalysis_RespondersByEngagement

these_labels = categorical({'NE','E','ne','e'});
these_labels = reordercats(these_labels,{'NE','E','ne','e'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq_nonengaged = responders_seq_unpaired_sd1(engagement_seq_unpaired_sd1==0);
this_data_seq_engaged = responders_seq_unpaired_sd1(engagement_seq_unpaired_sd1==1);
this_data_ctrl_nonengaged = responders_ctrl_unpaired_sd1(engagement_ctrl_unpaired_sd1==0);
this_data_ctrl_engaged = responders_ctrl_unpaired_sd1(engagement_ctrl_unpaired_sd1==1);
this_data = nan(10,4);
this_data(1:length(this_data_seq_nonengaged),1) = this_data_seq_nonengaged;
this_data(1:length(this_data_seq_engaged),2) = this_data_seq_engaged;
this_data(1:length(this_data_ctrl_nonengaged),3) = this_data_ctrl_nonengaged;
this_data(1:length(this_data_ctrl_engaged),4) = this_data_ctrl_engaged;
v = bar(these_labels,[nanmean(this_data_seq_nonengaged),nanmean(this_data_seq_engaged),nanmean(this_data_ctrl_nonengaged),nanmean(this_data_ctrl_engaged)]);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.gray; v.CData(3,:) = p.col.darkGray; v.CData(4,:) = p.col.gray; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq_nonengaged)
    plot(these_labels(1),this_data_seq_nonengaged,'.','Color',p.col.seq)
end
for i=1:length(this_data_seq_engaged)
    plot(these_labels(2),this_data_seq_engaged,'.','Color',p.col.seq)
end
for i=1:length(this_data_ctrl_nonengaged)
    plot(these_labels(3),this_data_ctrl_nonengaged,'.','Color',p.col.ctrl)
end
for i=1:length(this_data_ctrl_runner)
    plot(these_labels(4),this_data_ctrl_engaged,'.','Color',p.col.ctrl)
end
ylim([0,30])
yticks([0,10,20,30])
ylabel({'Photoactivated cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_RespondersByEngagement.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_RespondersByEngagement.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_RespondersByEngagement.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_RespondersByEngagement.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_RespondersByEngagement\nn(seq_nonengaged)=',num2str(length(rmmissing(this_data_seq_nonengaged))),'\nn(seq_engaged)=',num2str(length(rmmissing(this_data_seq_engaged))),'\nn(ctrl_nonengaged)=',num2str(length(rmmissing(this_data_ctrl_nonengaged))),'\nn(ctrl_engaged)=',num2str(length(rmmissing(this_data_ctrl_engaged))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_RespondersByEngagement\nn(seq_nonengaged)=',num2str(length(rmmissing(this_data_seq_nonengaged))),'\nn(seq_engaged)=',num2str(length(rmmissing(this_data_seq_engaged))),'\nn(ctrl_nonengaged)=',num2str(length(rmmissing(this_data_ctrl_nonengaged))),'\nn(ctrl_engaged)=',num2str(length(rmmissing(this_data_ctrl_engaged))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig4_ResponseAnalysis_AmplitudeByEngagement

these_labels = categorical({'NE','E','ne','e'});
these_labels = reordercats(these_labels,{'NE','E','ne','e'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq_nonengaged = nanmean(amplitude_seq_unpaired_sd1(engagement_seq_unpaired_sd1==0,:),2);
this_data_seq_engaged = nanmean(amplitude_seq_unpaired_sd1(engagement_seq_unpaired_sd1==1,:),2);
this_data_ctrl_nonengaged = nanmean(amplitude_ctrl_unpaired_sd1(engagement_ctrl_unpaired_sd1==0,:),2);
this_data_ctrl_engaged = nanmean(amplitude_ctrl_unpaired_sd1(engagement_ctrl_unpaired_sd1==1,:),2);
this_data = nan(10,4);
this_data(1:length(this_data_seq_nonengaged),1) = this_data_seq_nonengaged;
this_data(1:length(this_data_seq_engaged),2) = this_data_seq_engaged;
this_data(1:length(this_data_ctrl_nonengaged),3) = this_data_ctrl_nonengaged;
this_data(1:length(this_data_ctrl_engaged),4) = this_data_ctrl_engaged;
v = bar(these_labels,[nanmean(this_data_seq_nonengaged),nanmean(this_data_seq_engaged),nanmean(this_data_ctrl_nonengaged),nanmean(this_data_ctrl_engaged)]);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.gray; v.CData(3,:) = p.col.darkGray; v.CData(4,:) = p.col.gray; v.EdgeColor = 'none'; 
for i=1:length(this_data_seq_nonengaged)
    plot(these_labels(1),this_data_seq_nonengaged,'.','Color',p.col.seq)
end
for i=1:length(this_data_seq_engaged)
    plot(these_labels(2),this_data_seq_engaged,'.','Color',p.col.seq)
end
for i=1:length(this_data_ctrl_nonengaged)
    plot(these_labels(3),this_data_ctrl_nonengaged,'.','Color',p.col.ctrl)
end
for i=1:length(this_data_ctrl_engaged)
    plot(these_labels(4),this_data_ctrl_engaged,'.','Color',p.col.ctrl)
end
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Photoactivation amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_AmplitudeByEngagement.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_AmplitudeByEngagement.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_AmplitudeByEngagement.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_AmplitudeByEngagement.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_AmplitudeByEngagement\nn(seq_nonengaged)=',num2str(length(rmmissing(this_data_seq_nonengaged))),'\nn(seq_engaged)=',num2str(length(rmmissing(this_data_seq_engaged))),'\nn(ctrl_nonengaged)=',num2str(length(rmmissing(this_data_ctrl_nonengaged))),'\nn(ctrl_engaged)=',num2str(length(rmmissing(this_data_ctrl_engaged))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_AmplitudeByEngagement\nn(seq_nonengaged)=',num2str(length(rmmissing(this_data_seq_nonengaged))),'\nn(seq_engaged)=',num2str(length(rmmissing(this_data_seq_engaged))),'\nn(ctrl_nonengaged)=',num2str(length(rmmissing(this_data_ctrl_nonengaged))),'\nn(ctrl_engaged)=',num2str(length(rmmissing(this_data_ctrl_engaged))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig4_ResponseAnalysis_AmplitudeResolution

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_seq = amplitudeResolution_seq_paired_sd1;
this_data_ctrl = amplitudeResolution_ctrl_paired_sd1;
plot(this_data_ctrl','Color',mean([p.col.ctrl;p.col.white]),'LineWidth',0.5)
plot(this_data_seq','Color',mean([p.col.seq;p.col.white]),'LineWidth',0.5)
plot(nanmean(this_data_ctrl',2)','Color',p.col.ctrl,'LineWidth',1.5)
plot(nanmean(this_data_seq',2),'Color',p.col.seq,'LineWidth',1.5)
xlim([1,size(this_data_seq,2)+1]) % xlim([1,size(this_data_seq,2)+1]-0.5)
xticks([0:4:length([0:2.5:50])]+1) % xticks([0:4:length([0:2.5:50])]+0.5)
xticklabels({'0','10','20','30','40','50'})
% yticks([0,50,100])
% ylim([0,100])
% ytickformat('percentage')
xlabel({'Distance from closest';'laser beamlet (Um)'})
ylabel('Response amplitude (all cells, S)')

savefig(F,[save_root_fig,'\Fig4_ResponseAnalysis_AmplitudeResolution.fig']);
saveas(F,[save_root_png,'\Fig4_ResponseAnalysis_AmplitudeResolution.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ResponseAnalysis_AmplitudeResolution.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_ResponseAnalysis_AmplitudeResolution.txt'],'wt');
fprintf(fid,['\nResponseAnalysis_AmplitudeResolution\nn(seq)=',num2str(size(rmmissing(this_data_seq),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_ctrl),1)),'\n']); fclose(fid);
fprintf(['\nResponseAnalysis_AmplitudeResolution\nn(seq)=',num2str(size(rmmissing(this_data_seq),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_ctrl),1)),'\n']);


%% Mixed-effects model

% response variables
temp = responders_unpaired(:,2:5);
tbl_in.resp = temp(:);
temp = amplitude_unpaired(:,2:5);
tbl_in.amp = temp(:);

% fixed effects
temp = nan(d_info.numAnimals,4);
for i=1:d_info.numAnimals
    if d_info.group(i)==7
        temp(i,:) = [1,1,0,0];
    elseif d_info.group(i)==8
        temp(i,:) = [0,0,1,1];
    end
end
tbl_in.seq = nominal(temp(:));
temp = running_unpaired(:,2:5);
tbl_in.running = nominal(temp(:));
temp = engagement_unpaired(:,2:5);
tbl_in.engagement = nominal(temp(:));

% random effects
temp = repmat((1:d_info.numAnimals)',1,4);
tbl_in.animal = nominal(temp(:));
temp = repmat([1,1,2,2],d_info.numAnimals,1);
tbl_in.switch = nominal(temp(:));
temp = repmat([1,2,1,2],d_info.numAnimals,1);
tbl_in.switchday = nominal(temp(:));

tbl = table(tbl_in.resp,tbl_in.amp,tbl_in.seq,tbl_in.running,tbl_in.engagement,tbl_in.animal,tbl_in.switch,tbl_in.switchday,...
    'VariableNames',{'resp','amp','seq','running','engagement','animal','switch','switchday'});

lme = fitlme(tbl,'resp ~ seq + running + engagement + (1|animal) + (1|switch) + (1|switchday)')
%     Name                    Estimate    SE         tStat      DF    pValue        Lower      Upper 
%     {'(Intercept)' }          9.5626     1.9121      5.001    29    2.5295e-05     5.6518    13.473
%     {'seq_1'       }          1.8635    0.57547     3.2383    29     0.0030086    0.68657    3.0405
%     {'running_1'   }        -0.94603     1.4692    -0.6439    29        0.5247    -3.9509    2.0589
%     {'engagement_1'}          3.8225     1.2856     2.9732    29     0.0058782     1.1931    6.4519

lme = fitlme(tbl,'amp ~ seq + running + engagement + (1|animal) + (1|switch) + (1|switchday)')
%     Name                    Estimate     SE          tStat       DF    pValue        Lower        Upper   
%     {'(Intercept)' }           1.4962    0.099255      15.074    29    2.9571e-15       1.2932      1.6992
%     {'seq_1'       }         0.012092    0.032906     0.36749    29       0.71592    -0.055207    0.079392
%     {'running_1'   }        -0.048352    0.084342    -0.57329    29       0.57087     -0.22085     0.12415
%     {'engagement_1'}        0.0081458    0.073698     0.11053    29       0.91275     -0.14258     0.15887


%% --- In progress ---






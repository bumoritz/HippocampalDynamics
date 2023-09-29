%% Fig5_BehaviouralEffectByRunning

% import data using Summary_Master with ops.do_learningCurveSummary = true;
% d_info = selectDataset(d_info,'-g78-d2345-e01-r01-p01-l01-i00',sheet,path,ops);

save_root_fig = [path.root_summary,'figures\Fig5_fig\'];
save_root_png = [path.root_summary,'figures\Fig5_png\'];
save_root_pdf = [path.root_summary,'figures\Fig5_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig5_txt\'];


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

perfStim = nan(d_info.numAnimals,d_info.numDays,8);
perfStim_unpaired = nan(d_info.numAnimals,d_info.numDays,8);
perfStim_seq = nan(d_info.numAnimals,d_info.numDays,8);
perfStim_ctrl = nan(d_info.numAnimals,d_info.numDays,8);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp = d{i,j}.perf.blocks_general.correct_stim;
                temp1 = nanmean(reshape(temp,5,length(temp)/5),1);
                perfStim(i,j,1:length(temp1)) = temp1;
                if (d_info.group(i)==7 && (j==2 || j==3))
                    perfStim_seq(i,j,1:length(temp)) = temp;
                elseif (d_info.group(i)==8 && (j==4 || j==5))
                    perfStim_seq(i,j,1:length(temp)) = temp;
                elseif (d_info.group(i)==7 && (j==4 || j==5))
                    perfStim_ctrl(i,j,1:length(temp)) = temp;
                elseif (d_info.group(i)==8 && (j==2 || j==3))
                    perfStim_ctrl(i,j,1:length(temp)) = temp;
                end
            catch
            end
            perfStim_unpaired(i,j,:) = perfStim(i,j,:);
        end
    end
end
perfStim_seq_unpaired_sd1 = squeeze(nanmean(perfStim_seq(:,[2,4],:),2));
perfStim_seq_paired_sd1 = squeeze(nanmean(perfStim_seq(:,[2,4],:),2));
perfStim_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
perfStim_ctrl_unpaired_sd1 = squeeze(nanmean(perfStim_ctrl(:,[2,4],:),2));
perfStim_ctrl_paired_sd1 = squeeze(nanmean(perfStim_ctrl(:,[2,4],:),2));
perfStim_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;

perfStim_seq_paired_sd12 = cat(3,squeeze(nanmean(perfStim_seq(:,[2,4],:),2)),squeeze(nanmean(perfStim_seq(:,[3,5],:),2)));
perfStim_seq_paired_sd12(~these_animals_paired_bothDays,:,:) = NaN;
perfStim_ctrl_paired_sd12 = cat(3,squeeze(nanmean(perfStim_ctrl(:,[2,4],:),2)),squeeze(nanmean(perfStim_ctrl(:,[3,5],:),2)));
perfStim_ctrl_paired_sd12(~these_animals_paired_bothDays,:,:) = NaN;












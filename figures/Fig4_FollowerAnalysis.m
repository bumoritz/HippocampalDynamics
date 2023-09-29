%% Fig4_FollowerAnalysis

% import data using Summary_Master with ops.do_responseSummary = true;
% d_info = selectDataset(d_info,'-g78-d2345-e01-r01-p01-l01-i00',sheet,path,ops);

save_root_fig = [path.root_summary,'figures\Fig4_fig\'];
save_root_png = [path.root_summary,'figures\Fig4_png\'];
save_root_pdf = [path.root_summary,'figures\Fig4_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig4_txt\'];


%% Get data

followerBinSize = 25;
maxDistanceFromClosestTarget = 300;
posFollowers_skipFirstXbins = 3;
negFollowers_skipFirstXbins = 0;
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
                responders(i,j) = d{i,j}.resp.responders_main.numRespAll/40;
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

amplitude = nan(d_info.numAnimals,d_info.numDays,8);
amplitude_unpaired = nan(d_info.numAnimals,d_info.numDays,8);
amplitude_seq = nan(d_info.numAnimals,d_info.numDays,8);
amplitude_ctrl = nan(d_info.numAnimals,d_info.numDays,8);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp0 = nanmean(d{i,j}.resp.respAmps_cells_bw,1);
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

negFollowerProfile = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
negFollowerProfile_seq = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
negFollowerProfile_ctrl = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp1 = d{i,j}.flw.stats_sig_neg(:);
                temp2 = discretize(d{i,j}.resp.dist_closestLaser(:),[0:followerBinSize:maxDistanceFromClosestTarget]);
                temp = nan(nanmax(temp2),1);
                for k=negFollowers_skipFirstXbins+1:nanmax(temp2)
                    temp(k) = nanmean(temp1(find(temp2==k)));
                end
                negFollowerProfile(i,j,1:length(temp)) = temp;
                if d_info.group(i)==7 && (j==2 || j==3)
                    negFollowerProfile_seq(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    negFollowerProfile_seq(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    negFollowerProfile_ctrl(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    negFollowerProfile_ctrl(i,j,1:length(temp)) = temp;
                end
            catch
            end
        end
    end
end
negFollowerProfile_seq_paired_sd1 = squeeze(nanmean(negFollowerProfile_seq(:,[2,4],:),2));
negFollowerProfile_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
negFollowerProfile_ctrl_paired_sd1 = squeeze(nanmean(negFollowerProfile_ctrl(:,[2,4],:),2));
negFollowerProfile_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;

posFollowerProfile = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
posFollowerProfile_seq = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
posFollowerProfile_ctrl = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp1 = d{i,j}.flw.stats_sig_pos(:);
                temp2 = discretize(d{i,j}.resp.dist_closestLaser(:),[0:followerBinSize:maxDistanceFromClosestTarget]);
                temp = nan(nanmax(temp2),1);
                for k=posFollowers_skipFirstXbins+1:nanmax(temp2)
                    temp(k) = nanmean(temp1(find(temp2==k)));
                end
                posFollowerProfile(i,j,1:length(temp)) = temp;
                if d_info.group(i)==7 && (j==2 || j==3)
                    posFollowerProfile_seq(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    posFollowerProfile_seq(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    posFollowerProfile_ctrl(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    posFollowerProfile_ctrl(i,j,1:length(temp)) = temp;
                end
            catch
            end
        end
    end
end
posFollowerProfile_seq_paired_sd1 = squeeze(nanmean(posFollowerProfile_seq(:,[2,4],:),2));
posFollowerProfile_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
posFollowerProfile_ctrl_paired_sd1 = squeeze(nanmean(posFollowerProfile_ctrl(:,[2,4],:),2));
posFollowerProfile_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;

followerHeatmap = nan(d_info.numAnimals,d_info.numDays,2*maxDistanceFromClosestTarget/followerBinSize,2*maxDistanceFromClosestTarget/followerBinSize);
followerHeatmap_seq = nan(d_info.numAnimals,d_info.numDays,2*maxDistanceFromClosestTarget/followerBinSize,2*maxDistanceFromClosestTarget/followerBinSize);
followerHeatmap_ctrl = nan(d_info.numAnimals,d_info.numDays,2*maxDistanceFromClosestTarget/followerBinSize,2*maxDistanceFromClosestTarget/followerBinSize);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp1 = d{i,j}.flw.stats_sig_pos(:) - d{i,j}.flw.stats_sig_neg(:);
                
                temp = zeros(2*maxDistanceFromClosestTarget/followerBinSize);
                tempX = zeros(2*maxDistanceFromClosestTarget/followerBinSize);
                for x=1:2*maxDistanceFromClosestTarget/followerBinSize
                    for y=1:2*maxDistanceFromClosestTarget/followerBinSize
                
                        temp2 = discretize(d{i,j}.resp.pos_closestLaser_x(:),[-maxDistanceFromClosestTarget:followerBinSize:maxDistanceFromClosestTarget]);
                        temp3 = discretize(d{i,j}.resp.pos_closestLaser_y(:),[-maxDistanceFromClosestTarget:followerBinSize:maxDistanceFromClosestTarget]);
                        temp4 = intersect(find(temp2==x),find(temp3==y));
                        
                        % plot values
                        tempX(y,x) = length(rmmissing(temp1(temp4)));
                        if tempX(y,x) >= 30 % 10%30%20
                        	temp(y,x) = nanmean(temp1(temp4));
                        end

                        % draw circle of nan for exclusion zone (before smoothing)
                        if any([(x==10 && ismember(y,[12:13])),(x==11 && ismember(y,[11:14])),(x==12 && ismember(y,[10:15])),...
                                (x==15 && ismember(y,[12:13])),(x==14 && ismember(y,[11:14])),(x==13 && ismember(y,[10:15]))])
                            temp(y,x) = NaN;
                        end
                    end
                end
                temp = nanconv(temp,fspecial('gaussian',7,2*30/50),'nanout');
                %temp(isnan(temp))=0;
                
                followerHeatmap(i,j,:,:) = temp;
                if d_info.group(i)==7 && (j==2 || j==3)
                    followerHeatmap_seq(i,j,:,:) = temp;
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    followerHeatmap_seq(i,j,:,:) = temp;
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    followerHeatmap_ctrl(i,j,:,:) = temp;
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    followerHeatmap_ctrl(i,j,:,:) = temp;
                end
            catch
            end
        end
    end
end
followerHeatmap_ctrl_paired_sd1 = squeeze(nanmean(followerHeatmap_ctrl(:,[2,4],:,:),2));
followerHeatmap_ctrl_paired_sd1(~these_animals_paired_firstDay,:,:) = NaN;
followerHeatmap_ctrl_paired_sd1_avg = zeros(size(followerHeatmap_ctrl_paired_sd1,2));
for x=1:size(followerHeatmap_ctrl_paired_sd1,2)
    for y=1:size(followerHeatmap_ctrl_paired_sd1,2)
        if sum(followerHeatmap_ctrl_paired_sd1(:,y,x)~=0 & ~isnan(followerHeatmap_ctrl_paired_sd1(:,y,x))) >= 3
            followerHeatmap_ctrl_paired_sd1_avg(y,x) = nanmean(followerHeatmap_ctrl_paired_sd1(:,y,x));
        end
    end
end

avgSTAneg = nan(d_info.numAnimals,d_info.numDays,402);
avgSTAneg_seq = nan(d_info.numAnimals,d_info.numDays,402);
avgSTAneg_ctrl = nan(d_info.numAnimals,d_info.numDays,402);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                avgSTAneg(i,j,:) = nanmean(d{i,j}.flwsta.neg_sta,1);
                if d_info.group(i)==7 && (j==2 || j==3)
                    avgSTAneg_seq(i,j,:) = avgSTAneg(i,j,:);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    avgSTAneg_seq(i,j,:) = avgSTAneg(i,j,:);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    avgSTAneg_ctrl(i,j,:) = avgSTAneg(i,j,:);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    avgSTAneg_ctrl(i,j,:) = avgSTAneg(i,j,:);
                end
            catch
            end
        end
    end
end
avgSTAneg_seq_paired_sd1 = squeeze(nanmean(avgSTAneg_seq(:,[2,4],:),2));
avgSTAneg_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
avgSTAneg_ctrl_paired_sd1 = squeeze(nanmean(avgSTAneg_ctrl(:,[2,4],:),2));
avgSTAneg_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
avgSTAneg_ctrl_paired_sd1 = smoothdata(avgSTAneg_ctrl_paired_sd1,2,'gaussian',10);

avgSTApos = nan(d_info.numAnimals,d_info.numDays,402);
avgSTApos_seq = nan(d_info.numAnimals,d_info.numDays,402);
avgSTApos_ctrl = nan(d_info.numAnimals,d_info.numDays,402);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                avgSTApos(i,j,:) = nanmean(d{i,j}.flwsta.pos_sta,1);
                if d_info.group(i)==7 && (j==2 || j==3)
                    avgSTApos_seq(i,j,:) = avgSTApos(i,j,:);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    avgSTApos_seq(i,j,:) = avgSTApos(i,j,:);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    avgSTApos_ctrl(i,j,:) = avgSTApos(i,j,:);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    avgSTApos_ctrl(i,j,:) = avgSTApos(i,j,:);
                end
            catch
            end
        end
    end
end
avgSTApos_seq_paired_sd1 = squeeze(nanmean(avgSTApos_seq(:,[2,4],:),2));
avgSTApos_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
avgSTApos_ctrl_paired_sd1 = squeeze(nanmean(avgSTApos_ctrl(:,[2,4],:),2));
avgSTApos_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
avgSTApos_ctrl_paired_sd1 = smoothdata(avgSTApos_ctrl_paired_sd1,2,'gaussian',10);

avgSTAresp = nan(d_info.numAnimals,d_info.numDays,402);
avgSTAresp_seq = nan(d_info.numAnimals,d_info.numDays,402);
avgSTAresp_ctrl = nan(d_info.numAnimals,d_info.numDays,402);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                avgSTAresp(i,j,:) = nanmean(d{i,j}.flwsta.resp_sta,1);
                if d_info.group(i)==7 && (j==2 || j==3)
                    avgSTAresp_seq(i,j,:) = avgSTAresp(i,j,:);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    avgSTAresp_seq(i,j,:) = avgSTAresp(i,j,:);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    avgSTAresp_ctrl(i,j,:) = avgSTAresp(i,j,:);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    avgSTAresp_ctrl(i,j,:) = avgSTAresp(i,j,:);
                end
            catch
            end
        end
    end
end
avgSTAresp_seq_paired_sd1 = squeeze(nanmean(avgSTAresp_seq(:,[2,4],:),2));
avgSTAresp_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
avgSTAresp_ctrl_paired_sd1 = squeeze(nanmean(avgSTAresp_ctrl(:,[2,4],:),2));
avgSTAresp_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
avgSTAresp_ctrl_paired_sd1 = smoothdata(avgSTAresp_ctrl_paired_sd1,2,'gaussian',5);

% negFollowerAmplitudeProfile = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
% negFollowerAmplitudeProfile_seq = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
% negFollowerAmplitudeProfile_ctrl = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
% for i=1:d_info.numAnimals
%     for j=1:d_info.numDays
%         if d_info.presponsive(i,j)==1
%             try
%                 temp1 = d{i,j}.flw.stats_sig_neg(:);
%                 temp2 = discretize(d{i,j}.resp.dist_closestLaser(:),[0:followerBinSize:maxDistanceFromClosestTarget]);
%                 temp = nan(nanmax(temp2),1);
%                 for k=negFollowers_skipFirstXbins+1:nanmax(temp2)
%                     temp(k) = nanmean(temp1(find(temp2==k)));
%                 end
%                 negFollowerAmplitudeProfile(i,j,1:length(temp)) = temp;
%                 if d_info.group(i)==7 && (j==2 || j==3)
%                     negFollowerAmplitudeProfile_seq(i,j,1:length(temp)) = temp;
%                 elseif d_info.group(i)==8 && (j==4 || j==5)
%                     negFollowerAmplitudeProfile_seq(i,j,1:length(temp)) = temp;
%                 elseif d_info.group(i)==7 && (j==4 || j==5)
%                     negFollowerAmplitudeProfile_ctrl(i,j,1:length(temp)) = temp;
%                 elseif d_info.group(i)==8 && (j==2 || j==3)
%                     negFollowerAmplitudeProfile_ctrl(i,j,1:length(temp)) = temp;
%                 end
%             catch
%             end
%         end
%     end
% end
% negFollowerAmplitudeProfile_seq_paired_sd1 = squeeze(nanmean(negFollowerAmplitudeProfile_seq(:,[2,4],:),2));
% negFollowerAmplitudeProfile_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
% negFollowerAmplitudeProfile_ctrl_paired_sd1 = squeeze(nanmean(negFollowerAmplitudeProfile_ctrl(:,[2,4],:),2));
% negFollowerAmplitudeProfile_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
% 
% posFollowerAmplitudeProfile = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
% posFollowerAmplitudeProfile_seq = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
% posFollowerAmplitudeProfile_ctrl = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
% for i=1:d_info.numAnimals
%     for j=1:d_info.numDays
%         if d_info.presponsive(i,j)==1
%             try
%                 temp1 = d{i,j}.flw.stats_sig_pos(:);
%                 temp2 = discretize(d{i,j}.resp.dist_closestLaser(:),[0:followerBinSize:maxDistanceFromClosestTarget]);
%                 temp = nan(nanmax(temp2),1);
%                 for k=posFollowers_skipFirstXbins+1:nanmax(temp2)
%                     temp(k) = nanmean(temp1(find(temp2==k)));
%                 end
%                 posFollowerAmplitudeProfile(i,j,1:length(temp)) = temp;
%                 if d_info.group(i)==7 && (j==2 || j==3)
%                     posFollowerAmplitudeProfile_seq(i,j,1:length(temp)) = temp;
%                 elseif d_info.group(i)==8 && (j==4 || j==5)
%                     posFollowerAmplitudeProfile_seq(i,j,1:length(temp)) = temp;
%                 elseif d_info.group(i)==7 && (j==4 || j==5)
%                     posFollowerAmplitudeProfile_ctrl(i,j,1:length(temp)) = temp;
%                 elseif d_info.group(i)==8 && (j==2 || j==3)
%                     posFollowerAmplitudeProfile_ctrl(i,j,1:length(temp)) = temp;
%                 end
%             catch
%             end
%         end
%     end
% end
% posFollowerAmplitudeProfile_seq_paired_sd1 = squeeze(nanmean(posFollowerAmplitudeProfile_seq(:,[2,4],:),2));
% posFollowerAmplitudeProfile_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
% posFollowerAmplitudeProfile_ctrl_paired_sd1 = squeeze(nanmean(posFollowerAmplitudeProfile_ctrl(:,[2,4],:),2));
% posFollowerAmplitudeProfile_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;

followerAmplitudeProfile = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
followerAmplitudeProfile_seq = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
followerAmplitudeProfile_ctrl = nan(d_info.numAnimals,d_info.numDays,maxDistanceFromClosestTarget/followerBinSize);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp1 = d{i,j}.resp.avgAct_resp(:);
                temp2 = discretize(d{i,j}.resp.dist_closestLaser(:),[0:followerBinSize:maxDistanceFromClosestTarget]);
                temp = nan(nanmax(temp2),1);
                for k=posFollowers_skipFirstXbins+1:nanmax(temp2)
                    temp(k) = nanmean(temp1(find(temp2==k)));
                end
                followerAmplitudeProfile(i,j,1:length(temp)) = temp;
                if d_info.group(i)==7 && (j==2 || j==3)
                    followerAmplitudeProfile_seq(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    followerAmplitudeProfile_seq(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    followerAmplitudeProfile_ctrl(i,j,1:length(temp)) = temp;
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    followerAmplitudeProfile_ctrl(i,j,1:length(temp)) = temp;
                end
            catch
            end
        end
    end
end
followerAmplitudeProfile_seq_paired_sd1 = squeeze(nanmean(followerAmplitudeProfile_seq(:,[2,4],:),2));
followerAmplitudeProfile_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
followerAmplitudeProfile_ctrl_paired_sd1 = squeeze(nanmean(followerAmplitudeProfile_ctrl(:,[2,4],:),2));
followerAmplitudeProfile_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;

negFollowers = nan(d_info.numAnimals,d_info.numDays);
negFollowers_seq = nan(d_info.numAnimals,d_info.numDays);
negFollowers_ctrl = nan(d_info.numAnimals,d_info.numDays);
negFollowers_unpaired = nan(d_info.numAnimals,d_info.numDays);
negFollowers_cw = nan(d_info.numAnimals,d_info.numDays,40);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                %negFollowers(i,j) = nansum(d{i,j}.flw.stats_sig_neg(:)==1)/40;
                for k=1:40
                    temp1 = d{i,j}.flw.stats_sig_neg(:,k);
                    temp2 = discretize(d{i,j}.resp.dist_closestLaser(:,k),[0:followerBinSize:maxDistanceFromClosestTarget]);
                    temp = nan(length(temp1),1);
                    temp(temp2>negFollowers_skipFirstXbins) = temp1(temp2>negFollowers_skipFirstXbins);
                    negFollowers_cw(i,j,k) = nansum(temp);
                end
                negFollowers(i,j) = nanmean(negFollowers_cw(i,j,:),3);
                if d_info.group(i)==7 && (j==2 || j==3)
                    negFollowers_seq(i,j) = negFollowers(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    negFollowers_seq(i,j) = negFollowers(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    negFollowers_ctrl(i,j) = negFollowers(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    negFollowers_ctrl(i,j) = negFollowers(i,j);
                end
            catch
            end
            negFollowers_unpaired(i,j) = negFollowers(i,j);
        end
    end
    if (isnan(negFollowers(i,2)) || isnan(negFollowers(i,4)))
        these_animals_paired_firstDay(i) = 0;
    end
    if any([isnan(negFollowers(i,2)),isnan(negFollowers(i,3)),isnan(negFollowers(i,4)),isnan(negFollowers(i,5))])
        these_animals_paired_bothDays(i) = 0;
    end
end
negFollowers_seq_unpaired_sd1 = nanmean(negFollowers_seq(:,[2,4]),2);
negFollowers_seq_paired_sd1 = nanmean(negFollowers_seq(:,[2,4]),2);
negFollowers_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
negFollowers_ctrl_unpaired_sd1 = nanmean(negFollowers_ctrl(:,[2,4]),2);
negFollowers_ctrl_paired_sd1 = nanmean(negFollowers_ctrl(:,[2,4]),2);
negFollowers_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
negFollowers_seq_paired_sd12 = [nanmean(negFollowers_seq(:,[2,4]),2),nanmean(negFollowers_seq(:,[3,5]),2)];
negFollowers_seq_paired_sd12(~these_animals_paired_bothDays) = NaN;
negFollowers_ctrl_paired_sd12 = [nanmean(negFollowers_ctrl(:,[2,4]),2),nanmean(negFollowers_ctrl(:,[3,5]),2)];
negFollowers_ctrl_paired_sd12(~these_animals_paired_bothDays) = NaN;

posFollowers = nan(d_info.numAnimals,d_info.numDays);
posFollowers_seq = nan(d_info.numAnimals,d_info.numDays);
posFollowers_ctrl = nan(d_info.numAnimals,d_info.numDays);
posFollowers_unpaired = nan(d_info.numAnimals,d_info.numDays);
posFollowers_cw = nan(d_info.numAnimals,d_info.numDays,40);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                %posFollowers(i,j) = nansum(d{i,j}.flw.stats_sig_pos(:)==1)/40;
                for k=1:40
                    temp1 = d{i,j}.flw.stats_sig_pos(:,k);
                    temp2 = discretize(d{i,j}.resp.dist_closestLaser(:,k),[0:followerBinSize:maxDistanceFromClosestTarget]);
                    temp = nan(length(temp1),1);
                    temp(temp2>posFollowers_skipFirstXbins) = temp1(temp2>posFollowers_skipFirstXbins);
                    posFollowers_cw(i,j,k) = nansum(temp);
                end
                posFollowers(i,j) = nanmean(posFollowers_cw(i,j,:),3);
                if d_info.group(i)==7 && (j==2 || j==3)
                    posFollowers_seq(i,j) = posFollowers(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    posFollowers_seq(i,j) = posFollowers(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    posFollowers_ctrl(i,j) = posFollowers(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    posFollowers_ctrl(i,j) = posFollowers(i,j);
                end
            catch
            end
            posFollowers_unpaired(i,j) = posFollowers(i,j);
        end
    end
    if (isnan(posFollowers(i,2)) || isnan(posFollowers(i,4)))
        these_animals_paired_firstDay(i) = 0;
    end
    if any([isnan(posFollowers(i,2)),isnan(posFollowers(i,3)),isnan(posFollowers(i,4)),isnan(posFollowers(i,5))])
        these_animals_paired_bothDays(i) = 0;
    end
end
posFollowers_seq_unpaired_sd1 = nanmean(posFollowers_seq(:,[2,4]),2);
posFollowers_seq_paired_sd1 = nanmean(posFollowers_seq(:,[2,4]),2);
posFollowers_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
posFollowers_ctrl_unpaired_sd1 = nanmean(posFollowers_ctrl(:,[2,4]),2);
posFollowers_ctrl_paired_sd1 = nanmean(posFollowers_ctrl(:,[2,4]),2);
posFollowers_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
posFollowers_seq_paired_sd12 = [nanmean(posFollowers_seq(:,[2,4]),2),nanmean(posFollowers_seq(:,[3,5]),2)];
posFollowers_seq_paired_sd12(~these_animals_paired_bothDays) = NaN;
posFollowers_ctrl_paired_sd12 = [nanmean(posFollowers_ctrl(:,[2,4]),2),nanmean(posFollowers_ctrl(:,[3,5]),2)];
posFollowers_ctrl_paired_sd12(~these_animals_paired_bothDays) = NaN;

negFollowerAmplitude = nan(d_info.numAnimals,d_info.numDays);
negFollowerAmplitude_seq = nan(d_info.numAnimals,d_info.numDays);
negFollowerAmplitude_ctrl = nan(d_info.numAnimals,d_info.numDays);
negFollowerAmplitude_unpaired = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp0 = [];
                for k=1:40
                    temp1 = d{i,j}.flw.stats_sig_neg(:,k);
                    temp2 = discretize(d{i,j}.resp.dist_closestLaser(:,k),[0:followerBinSize:maxDistanceFromClosestTarget]);
                    these_idcs = intersect(find(temp1==1),find(temp2>negFollowers_skipFirstXbins));
                    if ~isempty(these_idcs)
                        temp0 = [temp0; d{i,j}.flw.avgAct_net(these_idcs,k)];
                    end
                end
                negFollowerAmplitude(i,j) = nanmean(temp0);
                if d_info.group(i)==7 && (j==2 || j==3)
                    negFollowerAmplitude_seq(i,j) = negFollowerAmplitude(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    negFollowerAmplitude_seq(i,j) = negFollowerAmplitude(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    negFollowerAmplitude_ctrl(i,j) = negFollowerAmplitude(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    negFollowerAmplitude_ctrl(i,j) = negFollowerAmplitude(i,j);
                end
            catch
            end
            negFollowerAmplitude_unpaired(i,j) = negFollowerAmplitude(i,j);
        end
    end
    if (isnan(negFollowerAmplitude(i,2)) || isnan(negFollowerAmplitude(i,4)))
        these_animals_paired_firstDay(i) = 0;
    end
    if any([isnan(negFollowerAmplitude(i,2)),isnan(negFollowerAmplitude(i,3)),isnan(negFollowerAmplitude(i,4)),isnan(negFollowerAmplitude(i,5))])
        these_animals_paired_bothDays(i) = 0;
    end
end
negFollowerAmplitude_seq_unpaired_sd1 = nanmean(negFollowerAmplitude_seq(:,[2,4]),2);
negFollowerAmplitude_seq_paired_sd1 = nanmean(negFollowerAmplitude_seq(:,[2,4]),2);
negFollowerAmplitude_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
negFollowerAmplitude_ctrl_unpaired_sd1 = nanmean(negFollowerAmplitude_ctrl(:,[2,4]),2);
negFollowerAmplitude_ctrl_paired_sd1 = nanmean(negFollowerAmplitude_ctrl(:,[2,4]),2);
negFollowerAmplitude_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
negFollowerAmplitude_seq_paired_sd12 = [nanmean(negFollowerAmplitude_seq(:,[2,4]),2),nanmean(negFollowerAmplitude_seq(:,[3,5]),2)];
negFollowerAmplitude_seq_paired_sd12(~these_animals_paired_bothDays) = NaN;
negFollowerAmplitude_ctrl_paired_sd12 = [nanmean(negFollowerAmplitude_ctrl(:,[2,4]),2),nanmean(negFollowerAmplitude_ctrl(:,[3,5]),2)];
negFollowerAmplitude_ctrl_paired_sd12(~these_animals_paired_bothDays) = NaN;

posFollowerAmplitude = nan(d_info.numAnimals,d_info.numDays);
posFollowerAmplitude_seq = nan(d_info.numAnimals,d_info.numDays);
posFollowerAmplitude_ctrl = nan(d_info.numAnimals,d_info.numDays);
posFollowerAmplitude_unpaired = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp0 = [];
                for k=1:40
                    temp1 = d{i,j}.flw.stats_sig_pos(:,k);
                    temp2 = discretize(d{i,j}.resp.dist_closestLaser(:,k),[0:followerBinSize:maxDistanceFromClosestTarget]);
                    these_idcs = intersect(find(temp1==1),find(temp2>posFollowers_skipFirstXbins));
                    if ~isempty(these_idcs)
                        temp0 = [temp0; d{i,j}.flw.avgAct_net(these_idcs,k)];
                    end
                end
                posFollowerAmplitude(i,j) = nanmean(temp0);
                if d_info.group(i)==7 && (j==2 || j==3)
                    posFollowerAmplitude_seq(i,j) = posFollowerAmplitude(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    posFollowerAmplitude_seq(i,j) = posFollowerAmplitude(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    posFollowerAmplitude_ctrl(i,j) = posFollowerAmplitude(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    posFollowerAmplitude_ctrl(i,j) = posFollowerAmplitude(i,j);
                end
            catch
            end
            posFollowerAmplitude_unpaired(i,j) = posFollowerAmplitude(i,j);
        end
    end
    if (isnan(posFollowerAmplitude(i,2)) || isnan(posFollowerAmplitude(i,4)))
        these_animals_paired_firstDay(i) = 0;
    end
    if any([isnan(posFollowerAmplitude(i,2)),isnan(posFollowerAmplitude(i,3)),isnan(posFollowerAmplitude(i,4)),isnan(posFollowerAmplitude(i,5))])
        these_animals_paired_bothDays(i) = 0;
    end
end
posFollowerAmplitude_seq_unpaired_sd1 = nanmean(posFollowerAmplitude_seq(:,[2,4]),2);
posFollowerAmplitude_seq_paired_sd1 = nanmean(posFollowerAmplitude_seq(:,[2,4]),2);
posFollowerAmplitude_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
posFollowerAmplitude_ctrl_unpaired_sd1 = nanmean(posFollowerAmplitude_ctrl(:,[2,4]),2);
posFollowerAmplitude_ctrl_paired_sd1 = nanmean(posFollowerAmplitude_ctrl(:,[2,4]),2);
posFollowerAmplitude_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
posFollowerAmplitude_seq_paired_sd12 = [nanmean(posFollowerAmplitude_seq(:,[2,4]),2),nanmean(posFollowerAmplitude_seq(:,[3,5]),2)];
posFollowerAmplitude_seq_paired_sd12(~these_animals_paired_bothDays) = NaN;
posFollowerAmplitude_ctrl_paired_sd12 = [nanmean(posFollowerAmplitude_ctrl(:,[2,4]),2),nanmean(posFollowerAmplitude_ctrl(:,[3,5]),2)];
posFollowerAmplitude_ctrl_paired_sd12(~these_animals_paired_bothDays) = NaN;

negFollowerAmplitudeStability = nan(d_info.numAnimals,d_info.numDays,8);
negFollowerAmplitudeStability_seq = nan(d_info.numAnimals,d_info.numDays,8);
negFollowerAmplitudeStability_ctrl = nan(d_info.numAnimals,d_info.numDays,8);
negFollowerAmplitudeStability_unpaired = nan(d_info.numAnimals,d_info.numDays,8);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp0 = [];
                for k=1:40
                    temp1 = d{i,j}.flw.stats_sig_neg(:,k);
                    temp2 = discretize(d{i,j}.resp.dist_closestLaser(:,k),[0:followerBinSize:maxDistanceFromClosestTarget]);
                    these_idcs = intersect(find(temp1==1),find(temp2>negFollowers_skipFirstXbins));
                    if ~isempty(these_idcs)
                        temp0 = [temp0; ...
                            nanmean(d{i,j}.flw.act_net(these_idcs,k,1:40),3),nanmean(d{i,j}.flw.act_net(these_idcs,k,41:80),3),...
                            nanmean(d{i,j}.flw.act_net(these_idcs,k,81:120),3),nanmean(d{i,j}.flw.act_net(these_idcs,k,121:160),3),...
                            nanmean(d{i,j}.flw.act_net(these_idcs,k,161:200),3),nanmean(d{i,j}.flw.act_net(these_idcs,k,201:240),3),...
                            nanmean(d{i,j}.flw.act_net(these_idcs,k,241:280),3),nanmean(d{i,j}.flw.act_net(these_idcs,k,281:320),3)];
                    end
                end
                negFollowerAmplitudeStability(i,j,:) = nanmean(temp0,1);
                if d_info.group(i)==7 && (j==2 || j==3)
                    negFollowerAmplitudeStability_seq(i,j,:) = negFollowerAmplitudeStability(i,j,:);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    negFollowerAmplitudeStability_seq(i,j,:) = negFollowerAmplitudeStability(i,j,:);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    negFollowerAmplitudeStability_ctrl(i,j,:) = negFollowerAmplitudeStability(i,j,:);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    negFollowerAmplitudeStability_ctrl(i,j,:) = negFollowerAmplitudeStability(i,j,:);
                end
            catch
            end
            negFollowerAmplitudeStability_unpaired(i,j,:) = negFollowerAmplitudeStability(i,j,:);
        end
    end
    if (isnan(negFollowerAmplitudeStability(i,2)) || isnan(negFollowerAmplitudeStability(i,4)))
        these_animals_paired_firstDay(i) = 0;
    end
    if any([isnan(negFollowerAmplitudeStability(i,2)),isnan(negFollowerAmplitudeStability(i,3)),isnan(negFollowerAmplitudeStability(i,4)),isnan(negFollowerAmplitudeStability(i,5))])
        these_animals_paired_bothDays(i) = 0;
    end
end
negFollowerAmplitudeStability_seq_unpaired_sd1 = squeeze(nanmean(negFollowerAmplitudeStability_seq(:,[2,4],:),2));
negFollowerAmplitudeStability_seq_paired_sd1 = squeeze(nanmean(negFollowerAmplitudeStability_seq(:,[2,4],:),2));
negFollowerAmplitudeStability_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
negFollowerAmplitudeStability_ctrl_unpaired_sd1 = squeeze(nanmean(negFollowerAmplitudeStability_ctrl(:,[2,4],:),2));
negFollowerAmplitudeStability_ctrl_paired_sd1 = squeeze(nanmean(negFollowerAmplitudeStability_ctrl(:,[2,4],:),2));
negFollowerAmplitudeStability_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;




posFollowerAmplitudeStability = nan(d_info.numAnimals,d_info.numDays,8);
posFollowerAmplitudeStability_seq = nan(d_info.numAnimals,d_info.numDays,8);
posFollowerAmplitudeStability_ctrl = nan(d_info.numAnimals,d_info.numDays,8);
posFollowerAmplitudeStability_unpaired = nan(d_info.numAnimals,d_info.numDays,8);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp0 = [];
                for k=1:40
                    temp1 = d{i,j}.flw.stats_sig_pos(:,k);
                    temp2 = discretize(d{i,j}.resp.dist_closestLaser(:,k),[0:followerBinSize:maxDistanceFromClosestTarget]);
                    these_idcs = intersect(find(temp1==1),find(temp2>posFollowers_skipFirstXbins));
                    if ~isempty(these_idcs)
                        temp0 = [temp0; ...
                            nanmean(d{i,j}.flw.act_net(these_idcs,k,1:40),3),nanmean(d{i,j}.flw.act_net(these_idcs,k,41:80),3),...
                            nanmean(d{i,j}.flw.act_net(these_idcs,k,81:120),3),nanmean(d{i,j}.flw.act_net(these_idcs,k,121:160),3),...
                            nanmean(d{i,j}.flw.act_net(these_idcs,k,161:200),3),nanmean(d{i,j}.flw.act_net(these_idcs,k,201:240),3),...
                            nanmean(d{i,j}.flw.act_net(these_idcs,k,241:280),3),nanmean(d{i,j}.flw.act_net(these_idcs,k,281:320),3)];
                    end
                end
                posFollowerAmplitudeStability(i,j,:) = nanmean(temp0,1);
                if d_info.group(i)==7 && (j==2 || j==3)
                    posFollowerAmplitudeStability_seq(i,j,:) = posFollowerAmplitudeStability(i,j,:);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    posFollowerAmplitudeStability_seq(i,j,:) = posFollowerAmplitudeStability(i,j,:);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    posFollowerAmplitudeStability_ctrl(i,j,:) = posFollowerAmplitudeStability(i,j,:);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    posFollowerAmplitudeStability_ctrl(i,j,:) = posFollowerAmplitudeStability(i,j,:);
                end
            catch
            end
            posFollowerAmplitudeStability_unpaired(i,j,:) = posFollowerAmplitudeStability(i,j,:);
        end
    end
    if (isnan(posFollowerAmplitudeStability(i,2)) || isnan(posFollowerAmplitudeStability(i,4)))
        these_animals_paired_firstDay(i) = 0;
    end
    if any([isnan(posFollowerAmplitudeStability(i,2)),isnan(posFollowerAmplitudeStability(i,3)),isnan(posFollowerAmplitudeStability(i,4)),isnan(posFollowerAmplitudeStability(i,5))])
        these_animals_paired_bothDays(i) = 0;
    end
end
posFollowerAmplitudeStability_seq_unpaired_sd1 = squeeze(nanmean(posFollowerAmplitudeStability_seq(:,[2,4],:),2));
posFollowerAmplitudeStability_seq_paired_sd1 = squeeze(nanmean(posFollowerAmplitudeStability_seq(:,[2,4],:),2));
posFollowerAmplitudeStability_seq_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
posFollowerAmplitudeStability_ctrl_unpaired_sd1 = squeeze(nanmean(posFollowerAmplitudeStability_ctrl(:,[2,4],:),2));
posFollowerAmplitudeStability_ctrl_paired_sd1 = squeeze(nanmean(posFollowerAmplitudeStability_ctrl(:,[2,4],:),2));
posFollowerAmplitudeStability_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;

negFollowerAmplitudeRespCorr = nan(d_info.numAnimals,d_info.numDays);
negFollowerAmplitudeRespCorr_seq = nan(d_info.numAnimals,d_info.numDays);
negFollowerAmplitudeRespCorr_ctrl = nan(d_info.numAnimals,d_info.numDays);
negFollowerAmplitudeRespCorr_unpaired = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp0 = [];
                for k=1:40
                    temp1 = d{i,j}.flw.stats_sig_neg(:,k);
                    temp2 = discretize(d{i,j}.resp.dist_closestLaser(:,k),[0:followerBinSize:maxDistanceFromClosestTarget]);
                    these_idcs = intersect(find(temp1==1),find(temp2>negFollowers_skipFirstXbins));
                    if ~isempty(these_idcs)
                        for m=1:length(these_idcs)
                            temp = corr(squeeze(nanmean(d{i,j}.flw.act_net(d{i,j}.resp.responders_main.idcs_resp{k},k,:),1)),squeeze(d{i,j}.flw.act_net(these_idcs(m),k,:)),'Type','Pearson','Rows','Complete');
                            temp0 = [temp0; temp];
                        end
                    end
                end
                negFollowerAmplitudeRespCorr(i,j) = nanmean(temp0);
                if d_info.group(i)==7 && (j==2 || j==3)
                    negFollowerAmplitudeRespCorr_seq(i,j) = negFollowerAmplitudeRespCorr(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    negFollowerAmplitudeRespCorr_seq(i,j) = negFollowerAmplitudeRespCorr(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    negFollowerAmplitudeRespCorr_ctrl(i,j) = negFollowerAmplitudeRespCorr(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    negFollowerAmplitudeRespCorr_ctrl(i,j) = negFollowerAmplitudeRespCorr(i,j);
                end
            catch
            end
            negFollowerAmplitudeRespCorr_unpaired(i,j) = negFollowerAmplitudeRespCorr(i,j);
        end
    end
    if (isnan(negFollowerAmplitudeRespCorr(i,2)) || isnan(negFollowerAmplitudeRespCorr(i,4)))
        these_animals_paired_firstDay(i) = 0;
    end
    if any([isnan(negFollowerAmplitudeRespCorr(i,2)),isnan(negFollowerAmplitudeRespCorr(i,3)),isnan(negFollowerAmplitudeRespCorr(i,4)),isnan(negFollowerAmplitudeRespCorr(i,5))])
        these_animals_paired_bothDays(i) = 0;
    end
end
negFollowerAmplitudeRespCorr_seq_unpaired_sd1 = nanmean(negFollowerAmplitudeRespCorr_seq(:,[2,4]),2);
negFollowerAmplitudeRespCorr_seq_paired_sd1 = nanmean(negFollowerAmplitudeRespCorr_seq(:,[2,4]),2);
negFollowerAmplitudeRespCorr_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
negFollowerAmplitudeRespCorr_ctrl_unpaired_sd1 = nanmean(negFollowerAmplitudeRespCorr_ctrl(:,[2,4]),2);
negFollowerAmplitudeRespCorr_ctrl_paired_sd1 = nanmean(negFollowerAmplitudeRespCorr_ctrl(:,[2,4]),2);
negFollowerAmplitudeRespCorr_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
negFollowerAmplitudeRespCorr_seq_paired_sd12 = [nanmean(negFollowerAmplitudeRespCorr_seq(:,[2,4]),2),nanmean(negFollowerAmplitudeRespCorr_seq(:,[3,5]),2)];
negFollowerAmplitudeRespCorr_seq_paired_sd12(~these_animals_paired_bothDays) = NaN;
negFollowerAmplitudeRespCorr_ctrl_paired_sd12 = [nanmean(negFollowerAmplitudeRespCorr_ctrl(:,[2,4]),2),nanmean(negFollowerAmplitudeRespCorr_ctrl(:,[3,5]),2)];
negFollowerAmplitudeRespCorr_ctrl_paired_sd12(~these_animals_paired_bothDays) = NaN;

posFollowerAmplitudeRespCorr = nan(d_info.numAnimals,d_info.numDays);
posFollowerAmplitudeRespCorr_seq = nan(d_info.numAnimals,d_info.numDays);
posFollowerAmplitudeRespCorr_ctrl = nan(d_info.numAnimals,d_info.numDays);
posFollowerAmplitudeRespCorr_unpaired = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp0 = [];
                for k=1:40
                    temp1 = d{i,j}.flw.stats_sig_pos(:,k);
                    temp2 = discretize(d{i,j}.resp.dist_closestLaser(:,k),[0:followerBinSize:maxDistanceFromClosestTarget]);
                    these_idcs = intersect(find(temp1==1),find(temp2>posFollowers_skipFirstXbins));
                    if ~isempty(these_idcs)
                        for m=1:length(these_idcs)
                            temp = corr(squeeze(nanmean(d{i,j}.flw.act_net(d{i,j}.resp.responders_main.idcs_resp{k},k,:),1)),squeeze(d{i,j}.flw.act_net(these_idcs(m),k,:)),'Type','Pearson','Rows','Complete');
                            temp0 = [temp0; temp];
                        end
                    end
                end
                posFollowerAmplitudeRespCorr(i,j) = nanmean(temp0);
                if d_info.group(i)==7 && (j==2 || j==3)
                    posFollowerAmplitudeRespCorr_seq(i,j) = posFollowerAmplitudeRespCorr(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    posFollowerAmplitudeRespCorr_seq(i,j) = posFollowerAmplitudeRespCorr(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    posFollowerAmplitudeRespCorr_ctrl(i,j) = posFollowerAmplitudeRespCorr(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    posFollowerAmplitudeRespCorr_ctrl(i,j) = posFollowerAmplitudeRespCorr(i,j);
                end
            catch
            end
            posFollowerAmplitudeRespCorr_unpaired(i,j) = posFollowerAmplitudeRespCorr(i,j);
        end
    end
    if (isnan(posFollowerAmplitudeRespCorr(i,2)) || isnan(posFollowerAmplitudeRespCorr(i,4)))
        these_animals_paired_firstDay(i) = 0;
    end
    if any([isnan(posFollowerAmplitudeRespCorr(i,2)),isnan(posFollowerAmplitudeRespCorr(i,3)),isnan(posFollowerAmplitudeRespCorr(i,4)),isnan(posFollowerAmplitudeRespCorr(i,5))])
        these_animals_paired_bothDays(i) = 0;
    end
end
posFollowerAmplitudeRespCorr_seq_unpaired_sd1 = nanmean(posFollowerAmplitudeRespCorr_seq(:,[2,4]),2);
posFollowerAmplitudeRespCorr_seq_paired_sd1 = nanmean(posFollowerAmplitudeRespCorr_seq(:,[2,4]),2);
posFollowerAmplitudeRespCorr_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
posFollowerAmplitudeRespCorr_ctrl_unpaired_sd1 = nanmean(posFollowerAmplitudeRespCorr_ctrl(:,[2,4]),2);
posFollowerAmplitudeRespCorr_ctrl_paired_sd1 = nanmean(posFollowerAmplitudeRespCorr_ctrl(:,[2,4]),2);
posFollowerAmplitudeRespCorr_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
posFollowerAmplitudeRespCorr_seq_paired_sd12 = [nanmean(posFollowerAmplitudeRespCorr_seq(:,[2,4]),2),nanmean(posFollowerAmplitudeRespCorr_seq(:,[3,5]),2)];
posFollowerAmplitudeRespCorr_seq_paired_sd12(~these_animals_paired_bothDays) = NaN;
posFollowerAmplitudeRespCorr_ctrl_paired_sd12 = [nanmean(posFollowerAmplitudeRespCorr_ctrl(:,[2,4]),2),nanmean(posFollowerAmplitudeRespCorr_ctrl(:,[3,5]),2)];
posFollowerAmplitudeRespCorr_ctrl_paired_sd12(~these_animals_paired_bothDays) = NaN;


%% --- Main figure ---

%% Fig4_FollowerAnalysis_Followers

these_labels = categorical({'Positive','Negative'});
these_labels = reordercats(these_labels,{'Positive','Negative'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos = posFollowers_ctrl_unpaired_sd1;
this_data_neg = negFollowers_ctrl_unpaired_sd1;
v = bar(these_labels,nanmean([this_data_pos,this_data_neg],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.pos; v.CData(2,:) = p.col.neg; v.EdgeColor = 'none'; 
for i=1:length(this_data_pos)
    plot(these_labels,[this_data_pos(i),this_data_neg(i)],'-k','LineWidth',1)
end
yticks([0,1,2,3])
ylim([0,3])
ylabel({'Follower cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Followers.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Followers.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Followers.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Followers.txt'],'wt');
[temp1,~,temp2] = signrank(this_data_pos,this_data_neg);
fprintf(fid,['\nFollowerAnalysis_Followers\nn(pos)=',num2str(length(rmmissing(this_data_pos))),'\nn(neg)=',num2str(length(rmmissing(this_data_neg))),...
    '\nmean(pos)=',num2str(nanmean(this_data_pos),4),'\nmean(neg)=',num2str(nanmean(this_data_neg),4),...
    '\nstd(pos)=',num2str(nanstd(this_data_pos),4),'\nstd(neg)=',num2str(nanstd(this_data_neg),4),...
    '\nsem(pos)=',num2str(nansem(this_data_pos,1),4),'\nsem(neg)=',num2str(nansem(this_data_neg,1),4),...
    '\nsignrank p=',num2str(temp1,4),...
    '\nsignrank T=',num2str(temp2.signedrank,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Followers\nn(pos)=',num2str(length(rmmissing(this_data_pos))),'\nn(neg)=',num2str(length(rmmissing(this_data_neg))),...
    '\nmean(pos)=',num2str(nanmean(this_data_pos),4),'\nmean(neg)=',num2str(nanmean(this_data_neg),4),...
    '\nstd(pos)=',num2str(nanstd(this_data_pos),4),'\nstd(neg)=',num2str(nanstd(this_data_neg),4),...
    '\nsem(pos)=',num2str(nansem(this_data_pos,1),4),'\nsem(neg)=',num2str(nansem(this_data_neg,1),4),...
    '\nsignrank p=',num2str(temp1,4),...
    '\nsignrank T=',num2str(temp2.signedrank,4),'\n']);


%% Fig4_FollowerAnalysis_FollowerProfile

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos = squeeze(nanmean(posFollowerProfile_ctrl(:,[2,4],:),2));
this_data_neg = -squeeze(nanmean(negFollowerProfile_ctrl(:,[2,4],:),2));
yline(0,'k-');
shadedErrorBar(1:size(this_data_pos,2),nanmean(this_data_pos'*100,2),nansem(this_data_pos'*100,2),'lineProps',p.col.pos);
shadedErrorBar(1:size(this_data_neg,2),nanmean(this_data_neg'*100,2),nansem(this_data_neg'*100,2),'lineProps',p.col.neg);
xlim([1,size(this_data_neg,2)+1]) % xlim([1,size(this_data_neg,2)+1]-0.5)
% xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+1) % xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+0.5)
% xticklabels({'0','50','100','150','200','250','300','350'})
xticks([0:4:length([0:followerBinSize:maxDistanceFromClosestTarget])]+1) % xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+0.5)
xticklabels({'0','100','200','300'})
yticks([-0.3,0,0.3])
ylim([-0.3,0.3])
ytickformat('percentage')
yticklabels({'0.3%','0%','0.3%'})
xlabel({'Distance from closest laser beamlet (Um)'})
%ylabel({'Follower recruitment p'})
ylabel({'Follower recruitment probability'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_FollowerProfile.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_FollowerProfile.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_FollowerProfile.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_FollowerProfile.txt'],'wt');
[temp1p,temp3p,temp2p] = kruskalwallis(this_data_pos); close;
[c,m,h] = multcompare(temp2p,'CType','dunn-sidak');
temp = rmmissing(rmmissing(this_data_pos,1,'MinNumMissing',4),2,'MinNumMissing',4);
temp0 = repmat(1:size(temp,2),size(temp,1),1)
statsp=kwtest([temp(:),temp0(:)]);
[temp1n,temp3n,temp2n] = kruskalwallis(this_data_neg); close;
[c,m,h] = multcompare(temp2n,'CType','dunn-sidak');
temp = rmmissing(rmmissing(this_data_neg,1,'MinNumMissing',4),2,'MinNumMissing',4);
temp0 = repmat(1:size(temp,2),size(temp,1),1)
statsp=kwtest([temp(:),temp0(:)]);
dunn(temp(:)',temp0(:)');
temp = [this_data_neg(:,5);this_data_neg(:,12)];
[~,pval,~,zstat] = ztest(temp,nanmean(temp),nanstd(temp));
fprintf(fid,['\nFollowerAnalysis_FollowerProfile\nn=',num2str(size(rmmissing(this_data_pos,'MinNumMissing',10),1)),...
    '\nkruskal-wallis p=',num2str(temp1p,4),...
    '\nkruskal-wallis H=',num2str(temp2p.sumt,4),...
    '\nkruskal-wallis p=',num2str(temp1n,4),...
    '\nkruskal-wallis H=',num2str(temp2n.sumt,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_FollowerProfile\nn=',num2str(size(rmmissing(this_data_pos,'MinNumMissing',10),1)),...
    '\nkruskal-wallis p=',num2str(temp1p,4),...
    '\nkruskal-wallis p=',num2str(temp1n,4),'\n']);


%% Fig4_FollowerAnalysis_FollowerProfile_ind

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos = squeeze(nanmean(posFollowerProfile_ctrl(:,[2,4],:),2));
this_data_neg = -squeeze(nanmean(negFollowerProfile_ctrl(:,[2,4],:),2));
yline(0,'k-');
plot(1:size(this_data_pos,2),this_data_pos'*100,'Color',mean([p.col.pos;p.col.white]),'LineWidth',0.5)
plot(1:size(this_data_neg,2),this_data_neg'*100,'Color',mean([p.col.neg;p.col.white]),'LineWidth',0.5)
plot(1:size(this_data_pos,2),nanmean(this_data_pos',2)*100,'Color',p.col.pos,'LineWidth',1.5)
plot(1:size(this_data_neg,2),nanmean(this_data_neg',2)*100,'Color',p.col.neg,'LineWidth',1.5)
xlim([1,size(this_data_neg,2)+1]) % xlim([1,size(this_data_neg,2)+1]-0.5)
% xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+1) % xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+0.5)
% xticklabels({'0','50','100','150','200','250','300','350'})
xticks([0:4:length([0:followerBinSize:maxDistanceFromClosestTarget])]+1) % xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+0.5)
xticklabels({'0','100','200','300'})
yticks([-0.6,-0.3,0,0.3,0.6])
ylim([-0.6,0.6])
ytickformat('percentage')
yticklabels({'0.6%','0.3%','0%','0.3%','0.6%'})
xlabel({'Distance from closest laser beamlet (Um)'})
%ylabel({'Follower recruitment p'})
ylabel({'Follower recruitment probability'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_FollowerProfile_ind.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_FollowerProfile_ind.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_FollowerProfile_ind.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_FollowerProfile_ind.txt'],'wt');
[temp1p,temp3p,temp2p] = kruskalwallis(this_data_pos); close;
[c,m,h] = multcompare(temp2p,'CType','dunn-sidak');
temp = rmmissing(rmmissing(this_data_pos,1,'MinNumMissing',4),2,'MinNumMissing',4);
temp0 = repmat(1:size(temp,2),size(temp,1),1)
statsp=kwtest([temp(:),temp0(:)]);
[temp1n,temp3n,temp2n] = kruskalwallis(this_data_neg); close;
[c,m,h] = multcompare(temp2n,'CType','dunn-sidak');
temp = rmmissing(rmmissing(this_data_neg,1,'MinNumMissing',4),2,'MinNumMissing',4);
temp0 = repmat(1:size(temp,2),size(temp,1),1)
statsp=kwtest([temp(:),temp0(:)]);
dunn(temp(:)',temp0(:)');
temp = [this_data_neg(:,5);this_data_neg(:,12)];
[~,pval,~,zstat] = ztest(temp,nanmean(temp),nanstd(temp));
fprintf(fid,['\nFollowerAnalysis_FollowerProfile\nn=',num2str(size(rmmissing(this_data_pos,'MinNumMissing',10),1)),...
    '\nkruskal-wallis p=',num2str(temp1p,4),...
    '\nkruskal-wallis H=',num2str(temp2p.sumt,4),...
    '\nkruskal-wallis p=',num2str(temp1n,4),...
    '\nkruskal-wallis H=',num2str(temp2n.sumt,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_FollowerProfile\nn=',num2str(size(rmmissing(this_data_pos,'MinNumMissing',10),1)),...
    '\nkruskal-wallis p=',num2str(temp1p,4),...
    '\nkruskal-wallis p=',num2str(temp1n,4),'\n']);


%% Fig4_FollowerAnalysis_Responders_2_Followers

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos_resp = responders_ctrl_unpaired_sd1;
this_data_pos_flw = posFollowers_ctrl_unpaired_sd1;
this_data_neg_resp = responders_ctrl_unpaired_sd1;
this_data_neg_flw = -negFollowers_ctrl_unpaired_sd1;
yline(0,'k-');
[this_corr_r_pos,this_corr_p_pos] = fitLine(this_data_pos_resp,this_data_pos_flw,p.col.black);
scatter(this_data_pos_resp,this_data_pos_flw,'.','SizeData',100,'MarkerEdgeColor',p.col.pos)
[this_corr_r_neg,this_corr_p_neg] = fitLine(this_data_neg_resp,this_data_neg_flw,p.col.black);
scatter(this_data_neg_resp,this_data_neg_flw,'.','SizeData',100,'MarkerEdgeColor',p.col.neg)
xlim([0,30])
xticks([0,10,20,30])
xlabel({'Photoactivated cells per cluster'})
ylim([-4,4])
yticks([-4,-2,0,2,4])
yticklabels({'4','2','0','2','4'})
ylabel({'Follower cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Responders_2_Followers.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Responders_2_Followers.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Responders_2_Followers.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Responders_2_Followers.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_Responders_2_Followers\nn(pos_resp)=',num2str(length(rmmissing(this_data_pos_resp))),'\nn(pos_flw)=',num2str(length(rmmissing(this_data_pos_flw))),'\nn(neg_resp)=',num2str(length(rmmissing(this_data_neg_resp))),'\nn(neg_flw)=',num2str(length(rmmissing(this_data_neg_flw))),...
    '\nPearson rho=',num2str(this_corr_r_pos,4),'\nPearson p=',num2str(this_corr_p_pos,4),'\nPearson rho=',num2str(this_corr_r_neg,4),'\nPearson p=',num2str(this_corr_p_neg,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Responders_2_Followers\nn(pos_resp)=',num2str(length(rmmissing(this_data_pos_resp))),'\nn(pos_flw)=',num2str(length(rmmissing(this_data_pos_flw))),'\nn(neg_resp)=',num2str(length(rmmissing(this_data_neg_resp))),'\nn(neg_flw)=',num2str(length(rmmissing(this_data_neg_flw))),...
    '\nPearson rho=',num2str(this_corr_r_pos,4),'\nPearson p=',num2str(this_corr_p_pos,4),'\nPearson rho=',num2str(this_corr_r_neg,4),'\nPearson p=',num2str(this_corr_p_neg,4),'\n']);


%% Fig4_FollowerAnalysis_AvgSTA

F = paper_figure([0,0.5,mm2inch(0.5*34),mm2inch(34)]); hold on;

info = get_info;
distanceScaling = 1; ymax = 5; numTraces = 3;
this_t = p.general.t_unbinned;
patch([interp1(this_t,1:length(this_t),0)-5,interp1(this_t,1:length(this_t),0)-5,interp1(this_t,1:length(this_t),0)-2,interp1(this_t,1:length(this_t),0)-2],...
    [0,ymax,ymax,0],mean([p.col.gray;p.col.white]),'EdgeColor','none');
patch([interp1(this_t,1:length(this_t),0)+4,interp1(this_t,1:length(this_t),0)+4,interp1(this_t,1:length(this_t),0)+7,interp1(this_t,1:length(this_t),0)+7],...
    [0,ymax,ymax,0],mean([p.col.gray;p.col.white]),'EdgeColor','none');
patch([interp1(this_t,1:length(this_t),0)-1,interp1(this_t,1:length(this_t),0)-1,interp1(this_t,1:length(this_t),0)+3,interp1(this_t,1:length(this_t),0)+3],...
    [0,ymax,ymax,0],mean([p.col.photostim;p.col.white]),'EdgeColor','none');

i=1; shadedErrorBar(1:length(p.general.t_unbinned),nanmean(avgSTAresp_ctrl_paired_sd1,1)+(numTraces+1-i)*distanceScaling,nansem(avgSTAresp_ctrl_paired_sd1,1),'lineProps',p.col.black);
i=2; shadedErrorBar(1:length(p.general.t_unbinned),nanmean(avgSTApos_ctrl_paired_sd1,1)+(numTraces+1-i)*distanceScaling,nansem(avgSTApos_ctrl_paired_sd1,1),'lineProps',p.col.pos);
i=3; shadedErrorBar(1:length(p.general.t_unbinned),nanmean(avgSTAneg_ctrl_paired_sd1,1)+(numTraces+1-i)*distanceScaling,nansem(avgSTAneg_ctrl_paired_sd1,1),'lineProps',p.col.neg);

%xlim([86,120])
%xlim([60,150])
xlim([80,120])
ylim([0.65,ymax-0.1])
box off; axis off;

nanmax(nanmean(avgSTAresp_ctrl_paired_sd1,1)) - nanmin(nanmean(avgSTAresp_ctrl_paired_sd1,1));

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_AvgSTA.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_AvgSTA.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_AvgSTA.pdf']); set(gcf,'Color',[1,1,1])


%% Fig4_FollowerHeatmap

this_data = squeeze(nanmean(followerHeatmap_ctrl_paired_sd1,1)); %followerHeatmap_ctrl_paired_sd1_avg; %squeeze(nanmean(followerHeatmap_ctrl_paired_sd1,1)); %temp;

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

h=imagesc(flipud(this_data),[-0.002,0.002]);
set(h,'AlphaData',~isnan(flipud(this_data)));
set(gca,'color',p.col.gray);
colormap(redblue)

daspect([1,1,1]);
xlim([0.5,24.5])
ylim([0.5,24.5])
box on;
set(gca,'xtick',[0.5,12.5,24.5]);
set(gca,'xticklabel',{'-300','0','300'});
set(gca,'ytick',[0.5,12.5,24.5]);
set(gca,'yticklabel',{'-300','0','300'});
scatter(12.5,12.5,'k+','LineWidth',0.5,'SizeData',5)
% xlabel({'Distance (x) from closest';'laser beamlet (Um)'})
% ylabel({'Distance (y) from closest';'laser beamlet (Um)'})
xlabel({'Distance (x) from closest'})
ylabel({'Distance (y) from closest'})

savefig(F,[save_root_fig,'\Fig4_FollowerHeatmap.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerHeatmap.png']);
saveas(F,[save_root_pdf,'\Fig4_FollowerHeatmap.pdf']);


%% --- Supplementary figure ---

%% Fig4_FollowerAnalysis_PosFollowers_2_NegFollowers

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos = posFollowers_ctrl_unpaired_sd1;
this_data_neg = negFollowers_ctrl_unpaired_sd1;
scatter(this_data_pos,this_data_neg,'.','SizeData',100,'MarkerEdgeColor',p.col.darkGray)
[this_corr_r,this_corr_p] = fitLine(this_data_pos,this_data_neg,p.col.black);
xlim([0,1])
xticks([0,0.5,1])
xlabel({'Positive follower';'cells per cluster'})
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Negative follower';'cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_PosFollowers_2_NegFollowers.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_PosFollowers_2_NegFollowers.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_PosFollowers_2_NegFollowers.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_PosFollowers_2_NegFollowers.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_PosFollowers_2_NegFollowers\nn(pos)=',num2str(length(rmmissing(this_data_pos))),'\nn(neg)=',num2str(length(rmmissing(this_data_neg))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_PosFollowers_2_NegFollowers\nn(pos)=',num2str(length(rmmissing(this_data_pos))),'\nn(neg)=',num2str(length(rmmissing(this_data_neg))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_FollowerAnalysis_Amplitude

these_labels = categorical({'Positive','Negative'});
these_labels = reordercats(these_labels,{'Positive','Negative'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos = posFollowerAmplitude_ctrl_unpaired_sd1;
this_data_neg = negFollowerAmplitude_ctrl_unpaired_sd1;
v = bar(these_labels,nanmean([this_data_pos,this_data_neg],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.pos; v.CData(2,:) = p.col.neg; v.EdgeColor = 'none'; 
for i=1:length(this_data_pos)
    plot(these_labels,[this_data_pos(i),this_data_neg(i)],'-k','LineWidth',1)
end
yticks([-0.5,-0.25,0,0.25,0.5])
ylim([-0.5,0.5])
ylabel({'Follower amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_FollowerAmplitude.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_FollowerAmplitude.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_FollowerAmplitude.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_FollowerAmplitude.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_FollowerAmplitude\nn(pos)=',num2str(length(rmmissing(this_data_pos))),'\nn(neg)=',num2str(length(rmmissing(this_data_neg))),...
    '\nmean(pos)=',num2str(nanmean(this_data_pos),4),'\nmean(neg)=',num2str(nanmean(this_data_neg),4),...
    '\nstd(pos)=',num2str(nanstd(this_data_pos),4),'\nstd(neg)=',num2str(nanstd(this_data_neg),4),...
    '\nsem(pos)=',num2str(nansem(this_data_pos,1),4),'\nsem(neg)=',num2str(nansem(this_data_neg,1),4),...
    '\nsignrank p=',num2str(signrank(this_data_pos,this_data_neg),4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_FollowerAmplitude\nn(pos)=',num2str(length(rmmissing(this_data_pos))),'\nn(neg)=',num2str(length(rmmissing(this_data_neg))),...
    '\nmean(pos)=',num2str(nanmean(this_data_pos),4),'\nmean(neg)=',num2str(nanmean(this_data_neg),4),...
    '\nstd(pos)=',num2str(nanstd(this_data_pos),4),'\nstd(neg)=',num2str(nanstd(this_data_neg),4),...
    '\nsem(pos)=',num2str(nansem(this_data_pos,1),4),'\nsem(neg)=',num2str(nansem(this_data_neg,1),4),...
    '\nsignrank p=',num2str(signrank(this_data_pos,this_data_neg),4),'\n']);


%% --- Reserve ---

%% Fig4_FollowerAnalysis_NegFollowerProfile

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_neg = squeeze(nanmean(negFollowerProfile_ctrl(:,[2,4],:),2));
plot(1:size(this_data_neg,2),this_data_neg'*100,'Color',mean([p.col.neg;p.col.white]),'LineWidth',0.5)
plot(1:size(this_data_neg,2),nanmean(this_data_neg',2)*100,'Color',p.col.neg,'LineWidth',1.5)
xlim([1,size(this_data_neg,2)+1]) % xlim([1,size(this_data_neg,2)+1]-0.5)
% xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+1) % xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+0.5)
% xticklabels({'0','50','100','150','200','250','300','350'})
xticks([0:4:length([0:followerBinSize:maxDistanceFromClosestTarget])]+1) % xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+0.5)
xticklabels({'0','100','200','300'})
yticks([0,0.3,0.6])
ylim([0,0.6])
ytickformat('percentage')
xlabel({'Distance from closest';'laser beamlet (Um)'})
ylabel({'Negative follower';'recruitment probability'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_NegFollowerProfile.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_NegFollowerProfile.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_NegFollowerProfile.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_NegFollowerProfile.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_NegFollowerProfile\nn=',num2str(size(rmmissing(this_data_neg,'MinNumMissing',10),1)),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_NegFollowerProfile\nn=',num2str(size(rmmissing(this_data_neg,'MinNumMissing',10),1)),'\n']);


%% Fig4_FollowerAnalysis_PosFollowerProfile

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos = squeeze(nanmean(posFollowerProfile_ctrl(:,[2,4],:),2));
plot(1:size(this_data_pos,2),this_data_pos'*100,'Color',mean([p.col.pos;p.col.white]),'LineWidth',0.5)
plot(1:size(this_data_pos,2),nanmean(this_data_pos',2)*100,'Color',p.col.pos,'LineWidth',1.5)
xlim([1,size(this_data_pos,2)+1]) % xlim([1,size(this_data_neg,2)+1]-0.5)
% xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+1) % xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+0.5)
% xticklabels({'0','50','100','150','200','250','300','350'})
xticks([0:4:length([0:followerBinSize:maxDistanceFromClosestTarget])]+1) % xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+0.5)
xticklabels({'0','100','200','300'})
yticks([0,0.3,0.6])
ylim([0,0.6])
ytickformat('percentage')
xlabel({'Distance from closest';'laser beamlet (Um)'})
ylabel({'Positive follower';'recruitment probability'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_PosFollowerProfile.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_PosFollowerProfile.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_PosFollowerProfile.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_PosFollowerProfile.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_PosFollowerProfile\nn=',num2str(size(rmmissing(this_data_pos,'MinNumMissing',10),1)),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_PosFollowerProfile\nn=',num2str(size(rmmissing(this_data_pos,'MinNumMissing',10),1)),'\n']);


%% Fig4_FollowerAnalysis_Responders_2_PosFollowers

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_resp = responders_ctrl_unpaired_sd1;
this_data_flw = posFollowers_ctrl_unpaired_sd1;
scatter(this_data_resp,this_data_flw,'.','SizeData',100,'MarkerEdgeColor',p.col.pos)
[this_corr_r,this_corr_p] = fitLine(this_data_resp,this_data_flw,p.col.black);
xlim([0,30])
xticks([0,10,20,30])
xlabel({'Photoactivated cells per cluster'})
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Positive follower';'cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Responders_2_PosFollowers.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Responders_2_PosFollowers.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Responders_2_PosFollowers.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Responders_2_PosFollowers.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_Responders_2_PosFollowers\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Responders_2_PosFollowers\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_FollowerAnalysis_Responders_2_NegFollowers

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_resp = responders_ctrl_unpaired_sd1;
this_data_flw = negFollowers_ctrl_unpaired_sd1;
scatter(this_data_resp,this_data_flw,'.','SizeData',100,'MarkerEdgeColor',p.col.neg)
[this_corr_r,this_corr_p] = fitLine(this_data_resp,this_data_flw,p.col.black);
xlim([0,30])
xticks([0,10,20,30])
xlabel({'Photoactivated cells per cluster'})
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Negative follower';'cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Responders_2_NegFollowers.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Responders_2_NegFollowers.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Responders_2_NegFollowers.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Responders_2_NegFollowers.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_Responders_2_NegFollowers\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Responders_2_NegFollowers\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_FollowerAnalysis_Responders_2_PosFollowerAmplitude

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_resp = responders_ctrl_unpaired_sd1;
this_data_flw = posFollowerAmplitude_ctrl_unpaired_sd1;
scatter(this_data_resp,this_data_flw,'.','SizeData',100,'MarkerEdgeColor',p.col.pos)
[this_corr_r,this_corr_p] = fitLine(this_data_resp,this_data_flw,p.col.black);
xlim([0,30])
xticks([0,10,20,30])
xlabel({'Photoactivated cells per cluster'})
yticks([0,0.25,0.5])
ylim([0,0.5])
ylabel({'Positive follower amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Responders_2_PosFollowerAmplitude.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Responders_2_PosFollowerAmplitude.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Responders_2_PosFollowerAmplitude.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Responders_2_PosFollowerAmplitude.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_Responders_2_PosFollowerAmplitude\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Responders_2_PosFollowerAmplitude\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_FollowerAnalysis_Responders_2_NegFollowerAmplitude

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_resp = responders_ctrl_unpaired_sd1;
this_data_flw = negFollowerAmplitude_ctrl_unpaired_sd1;
scatter(this_data_resp,this_data_flw,'.','SizeData',100,'MarkerEdgeColor',p.col.neg)
[this_corr_r,this_corr_p] = fitLine(this_data_resp,this_data_flw,p.col.black);
xlim([0,30])
xticks([0,10,20,30])
xlabel({'Photoactivated cells per cluster'})
yticks([-0.5,-0.25,0])
ylim([-0.5,0])
ylabel({'Negative follower amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Responders_2_NegFollowerAmplitude.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Responders_2_NegFollowerAmplitude.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Responders_2_NegFollowerAmplitude.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Responders_2_NegFollowerAmplitude.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_Responders_2_NegFollowerAmplitude\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Responders_2_NegFollowerAmplitude\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_FollowerAnalysis_Amplitude_2_PosFollowerAmplitude

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_resp = amplitude_ctrl_unpaired_sd1(:,1);
this_data_flw = posFollowerAmplitude_ctrl_unpaired_sd1;
scatter(this_data_resp,this_data_flw,'.','SizeData',100,'MarkerEdgeColor',p.col.pos)
[this_corr_r,this_corr_p] = fitLine(this_data_resp,this_data_flw,p.col.black);
xlim([0,3])
xticks([0,1,2,3])
xlabel({'Photoactivation amplitude (S)'})
yticks([0,0.25,0.5])
ylim([0,0.5])
ylabel({'Positive follower amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Amplitude_2_PosFollowerAmplitude.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Amplitude_2_PosFollowerAmplitude.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Amplitude_2_PosFollowerAmplitude.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Amplitude_2_PosFollowerAmplitude.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_Amplitude_2_PosFollowerAmplitude\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Amplitude_2_PosFollowerAmplitude\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_FollowerAnalysis_Amplitude_2_NegFollowerAmplitude

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_resp = amplitude_ctrl_unpaired_sd1(:,1);
this_data_flw = negFollowerAmplitude_ctrl_unpaired_sd1;
scatter(this_data_resp,this_data_flw,'.','SizeData',100,'MarkerEdgeColor',p.col.neg)
[this_corr_r,this_corr_p] = fitLine(this_data_resp,this_data_flw,p.col.black);
xlim([0,3])
xticks([0,1,2,3])
xlabel({'Photoactivation amplitude (S)'})
yticks([-0.5,-0.25,0])
ylim([-0.5,0])
ylabel({'Negative follower amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Amplitude_2_NegFollowerAmplitude.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Amplitude_2_NegFollowerAmplitude.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Amplitude_2_NegFollowerAmplitude.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Amplitude_2_NegFollowerAmplitude.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_Amplitude_2_NegFollowerAmplitude\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Amplitude_2_NegFollowerAmplitude\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_FollowerAnalysis_Amplitude_2_PosFollowers

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_resp = amplitude_ctrl_unpaired_sd1(:,1);
this_data_flw = posFollowers_ctrl_unpaired_sd1;
scatter(this_data_resp,this_data_flw,'.','SizeData',100,'MarkerEdgeColor',p.col.pos)
[this_corr_r,this_corr_p] = fitLine(this_data_resp,this_data_flw,p.col.black);
xlim([0,3])
xticks([0,1,2,3])
xlabel({'Photoactivation amplitude (S)'})
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Positive follower';'cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Amplitude_2_PosFollowers.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Amplitude_2_PosFollowers.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Amplitude_2_PosFollowers.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Amplitude_2_PosFollowers.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_Amplitude_2_PosFollowers\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Amplitude_2_PosFollowers\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_FollowerAnalysis_Amplitude_2_NegFollowers

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_resp = amplitude_ctrl_unpaired_sd1(:,1);
this_data_flw = negFollowers_ctrl_unpaired_sd1;
scatter(this_data_resp,this_data_flw,'.','SizeData',100,'MarkerEdgeColor',p.col.neg)
[this_corr_r,this_corr_p] = fitLine(this_data_resp,this_data_flw,p.col.black);
xlim([0,3])
xticks([0,1,2,3])
xlabel({'Photoactivation amplitude (S)'})
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Negative follower';'cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Amplitude_2_NegFollowers.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Amplitude_2_NegFollowers.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Amplitude_2_NegFollowers.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Amplitude_2_NegFollowers.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_Amplitude_2_NegFollowers\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Amplitude_2_NegFollowers\nn(resp)=',num2str(length(rmmissing(this_data_resp))),'\nn(flw)=',num2str(length(rmmissing(this_data_flw))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_FollowerAnalysis_Amplitude_2_FollowerAmplitude_trialwise

these_labels = categorical({'Positive','Negative'});
these_labels = reordercats(these_labels,{'Positive','Negative'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos = posFollowerAmplitudeRespCorr_ctrl_unpaired_sd1;
this_data_neg = negFollowerAmplitudeRespCorr_ctrl_unpaired_sd1;
v = bar(these_labels,nanmean([this_data_pos,this_data_neg],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.pos; v.CData(2,:) = p.col.neg; v.EdgeColor = 'none'; 
for i=1:length(this_data_pos)
    plot(these_labels,[this_data_pos(i),this_data_neg(i)],'-k','LineWidth',1)
end
yticks([-0.1,-0.05,0,0.05,0.1])
ylim([-0.1,0.1])
ylabel({'Correlation between photoactivation';'and follower amplitude'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_Amplitude_2_FollowerAmplitude_trialwise.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_Amplitude_2_FollowerAmplitude_trialwise.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_Amplitude_2_FollowerAmplitude_trialwise.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_Amplitude_2_FollowerAmplitude_trialwise.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_Amplitude_2_FollowerAmplitude_trialwise\nn(pos)=',num2str(length(rmmissing(this_data_pos))),'\nn(neg)=',num2str(length(rmmissing(this_data_neg))),...
    '\nsignrank p=',num2str(signrank(this_data_pos,this_data_neg),4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_Amplitude_2_FollowerAmplitude_trialwise\nn(pos)=',num2str(length(rmmissing(this_data_pos))),'\nn(neg)=',num2str(length(rmmissing(this_data_neg))),...
    '\nsignrank p=',num2str(signrank(this_data_pos,this_data_neg),4),'\n']);


%% Fig4_FollowerAnalysis_FollowersByRunning

these_labels = categorical({'NR','R','nr','r'});
these_labels = reordercats(these_labels,{'NR','R','nr','r'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos_nonrunner = posFollowers_ctrl_unpaired_sd1(running_seq_unpaired_sd1==0);
this_data_pos_runner = posFollowers_ctrl_unpaired_sd1(running_seq_unpaired_sd1==1);
this_data_neg_nonrunner = negFollowers_ctrl_unpaired_sd1(running_ctrl_unpaired_sd1==0);
this_data_neg_runner = negFollowers_ctrl_unpaired_sd1(running_ctrl_unpaired_sd1==1);
this_data = nan(20,4);
this_data(1:length(this_data_pos_nonrunner),1) = this_data_pos_nonrunner;
this_data(1:length(this_data_pos_runner),2) = this_data_pos_runner;
this_data(1:length(this_data_neg_nonrunner),3) = this_data_neg_nonrunner;
this_data(1:length(this_data_neg_runner),4) = this_data_neg_runner;
v = bar(these_labels,[nanmean(this_data_pos_nonrunner),nanmean(this_data_pos_runner),nanmean(this_data_neg_nonrunner),nanmean(this_data_neg_runner)]);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.nonrunner; v.CData(2,:) = p.col.runner; v.CData(3,:) = p.col.nonrunner; v.CData(4,:) = p.col.runner; v.EdgeColor = 'none'; 
for i=1:length(this_data_pos_nonrunner)
    plot(these_labels(1),this_data_pos_nonrunner,'.','Color',p.col.pos)
end
for i=1:length(this_data_pos_runner)
    plot(these_labels(2),this_data_pos_runner,'.','Color',p.col.pos)
end
for i=1:length(this_data_neg_nonrunner)
    plot(these_labels(3),this_data_neg_nonrunner,'.','Color',p.col.neg)
end
for i=1:length(this_data_neg_runner)
    plot(these_labels(4),this_data_neg_runner,'.','Color',p.col.neg)
end
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Follower cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_FollowersByRunning.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_FollowersByRunning.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_FollowersByRunning.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_FollowersByRunning.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_FollowersByRunning\nn(pos_nonrunner)=',num2str(length(rmmissing(this_data_pos_nonrunner))),'\nn(pos_runner)=',num2str(length(rmmissing(this_data_pos_runner))),'\nn(neg_nonrunner)=',num2str(length(rmmissing(this_data_neg_nonrunner))),'\nn(neg_runner)=',num2str(length(rmmissing(this_data_neg_runner))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_FollowersByRunning\nn(pos_nonrunner)=',num2str(length(rmmissing(this_data_pos_nonrunner))),'\nn(pos_runner)=',num2str(length(rmmissing(this_data_pos_runner))),'\nn(neg_nonrunner)=',num2str(length(rmmissing(this_data_neg_nonrunner))),'\nn(neg_runner)=',num2str(length(rmmissing(this_data_neg_runner))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig4_FollowerAnalysis_FollowerAmplitudeByRunning

these_labels = categorical({'NR','R','nr','r'});
these_labels = reordercats(these_labels,{'NR','R','nr','r'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos_nonrunner = posFollowerAmplitude_ctrl_unpaired_sd1(running_seq_unpaired_sd1==0);
this_data_pos_runner = posFollowerAmplitude_ctrl_unpaired_sd1(running_seq_unpaired_sd1==1);
this_data_neg_nonrunner = negFollowerAmplitude_ctrl_unpaired_sd1(running_ctrl_unpaired_sd1==0);
this_data_neg_runner = negFollowerAmplitude_ctrl_unpaired_sd1(running_ctrl_unpaired_sd1==1);
this_data = nan(20,4);
this_data(1:length(this_data_pos_nonrunner),1) = this_data_pos_nonrunner;
this_data(1:length(this_data_pos_runner),2) = this_data_pos_runner;
this_data(1:length(this_data_neg_nonrunner),3) = this_data_neg_nonrunner;
this_data(1:length(this_data_neg_runner),4) = this_data_neg_runner;
v = bar(these_labels,[nanmean(this_data_pos_nonrunner),nanmean(this_data_pos_runner),nanmean(this_data_neg_nonrunner),nanmean(this_data_neg_runner)]);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.nonrunner; v.CData(2,:) = p.col.runner; v.CData(3,:) = p.col.nonrunner; v.CData(4,:) = p.col.runner; v.EdgeColor = 'none'; 
for i=1:length(this_data_pos_nonrunner)
    plot(these_labels(1),this_data_pos_nonrunner,'.','Color',p.col.pos)
end
for i=1:length(this_data_pos_runner)
    plot(these_labels(2),this_data_pos_runner,'.','Color',p.col.pos)
end
for i=1:length(this_data_neg_nonrunner)
    plot(these_labels(3),this_data_neg_nonrunner,'.','Color',p.col.neg)
end
for i=1:length(this_data_neg_runner)
    plot(these_labels(4),this_data_neg_runner,'.','Color',p.col.neg)
end
yticks([-0.5,-0.25,0,0.25,0.5])
ylim([-0.5,0.5])
ylabel({'Follower amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_FollowerAmplitudeByRunning.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_FollowerAmplitudeByRunning.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_FollowerAmplitudeByRunning.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_FollowerAmplitudeByRunning.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_FollowerAmplitudeByRunning\nn(pos_nonrunner)=',num2str(length(rmmissing(this_data_pos_nonrunner))),'\nn(pos_runner)=',num2str(length(rmmissing(this_data_pos_runner))),'\nn(neg_nonrunner)=',num2str(length(rmmissing(this_data_neg_nonrunner))),'\nn(neg_runner)=',num2str(length(rmmissing(this_data_neg_runner))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_FollowerAmplitudeByRunning\nn(pos_nonrunner)=',num2str(length(rmmissing(this_data_pos_nonrunner))),'\nn(pos_runner)=',num2str(length(rmmissing(this_data_pos_runner))),'\nn(neg_nonrunner)=',num2str(length(rmmissing(this_data_neg_nonrunner))),'\nn(neg_runner)=',num2str(length(rmmissing(this_data_neg_runner))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig4_FollowerAnalysis_FollowersByEngagement

these_labels = categorical({'NE','E','ne','e'});
these_labels = reordercats(these_labels,{'NE','E','ne','e'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos_nonengaged = posFollowers_ctrl_unpaired_sd1(engagement_seq_unpaired_sd1==0);
this_data_pos_engaged = posFollowers_ctrl_unpaired_sd1(engagement_seq_unpaired_sd1==1);
this_data_neg_nonengaged = negFollowers_ctrl_unpaired_sd1(engagement_ctrl_unpaired_sd1==0);
this_data_neg_engaged = negFollowers_ctrl_unpaired_sd1(engagement_ctrl_unpaired_sd1==1);
this_data = nan(20,4);
this_data(1:length(this_data_pos_nonengaged),1) = this_data_pos_nonengaged;
this_data(1:length(this_data_pos_engaged),2) = this_data_pos_engaged;
this_data(1:length(this_data_neg_nonengaged),3) = this_data_neg_nonengaged;
this_data(1:length(this_data_neg_engaged),4) = this_data_neg_engaged;
v = bar(these_labels,[nanmean(this_data_pos_nonengaged),nanmean(this_data_pos_engaged),nanmean(this_data_neg_nonengaged),nanmean(this_data_neg_engaged)]);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.gray; v.CData(3,:) = p.col.darkGray; v.CData(4,:) = p.col.gray; v.EdgeColor = 'none'; 
for i=1:length(this_data_pos_nonengaged)
    plot(these_labels(1),this_data_pos_nonengaged,'.','Color',p.col.pos)
end
for i=1:length(this_data_pos_engaged)
    plot(these_labels(2),this_data_pos_engaged,'.','Color',p.col.pos)
end
for i=1:length(this_data_neg_nonengaged)
    plot(these_labels(3),this_data_neg_nonengaged,'.','Color',p.col.neg)
end
for i=1:length(this_data_neg_engaged)
    plot(these_labels(4),this_data_neg_engaged,'.','Color',p.col.neg)
end
ylim([0,3])
yticks([0,1,2,3])
ylabel({'Follower cells per cluster'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_FollowersByEngagement.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_FollowersByEngagement.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_FollowersByEngagement.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_FollowersByEngagement.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_FollowersByEngagement\nn(pos_nonengaged)=',num2str(length(rmmissing(this_data_pos_nonengaged))),'\nn(pos_engaged)=',num2str(length(rmmissing(this_data_pos_engaged))),'\nn(neg_nonengaged)=',num2str(length(rmmissing(this_data_neg_nonengaged))),'\nn(neg_engaged)=',num2str(length(rmmissing(this_data_neg_engaged))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_FollowersByEngagement\nn(pos_nonengaged)=',num2str(length(rmmissing(this_data_pos_nonengaged))),'\nn(pos_engaged)=',num2str(length(rmmissing(this_data_pos_engaged))),'\nn(neg_nonengaged)=',num2str(length(rmmissing(this_data_neg_nonengaged))),'\nn(neg_engaged)=',num2str(length(rmmissing(this_data_neg_engaged))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig4_FollowerAnalysis_FollowerAmplitudeByEngagement

these_labels = categorical({'NE','E','ne','e'});
these_labels = reordercats(these_labels,{'NE','E','ne','e'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos_nonengaged = posFollowerAmplitude_ctrl_unpaired_sd1(engagement_seq_unpaired_sd1==0);
this_data_pos_engaged = posFollowerAmplitude_ctrl_unpaired_sd1(engagement_seq_unpaired_sd1==1);
this_data_neg_nonengaged = negFollowerAmplitude_ctrl_unpaired_sd1(engagement_ctrl_unpaired_sd1==0);
this_data_neg_engaged = negFollowerAmplitude_ctrl_unpaired_sd1(engagement_ctrl_unpaired_sd1==1);
this_data = nan(20,4);
this_data(1:length(this_data_pos_nonengaged),1) = this_data_pos_nonengaged;
this_data(1:length(this_data_pos_engaged),2) = this_data_pos_engaged;
this_data(1:length(this_data_neg_nonengaged),3) = this_data_neg_nonengaged;
this_data(1:length(this_data_neg_engaged),4) = this_data_neg_engaged;
v = bar(these_labels,[nanmean(this_data_pos_nonengaged),nanmean(this_data_pos_engaged),nanmean(this_data_neg_nonengaged),nanmean(this_data_neg_engaged)]);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.gray; v.CData(3,:) = p.col.darkGray; v.CData(4,:) = p.col.gray; v.EdgeColor = 'none'; 
for i=1:length(this_data_pos_nonengaged)
    plot(these_labels(1),this_data_pos_nonengaged,'.','Color',p.col.pos)
end
for i=1:length(this_data_pos_engaged)
    plot(these_labels(2),this_data_pos_engaged,'.','Color',p.col.pos)
end
for i=1:length(this_data_neg_nonengaged)
    plot(these_labels(3),this_data_neg_nonengaged,'.','Color',p.col.neg)
end
for i=1:length(this_data_neg_engaged)
    plot(these_labels(4),this_data_neg_engaged,'.','Color',p.col.neg)
end
yticks([-0.5,-0.25,0,0.25,0.5])
ylim([-0.5,0.5])
ylabel({'Follower amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_FollowerAmplitudeByEngagement.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_FollowerAmplitudeByEngagement.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_FollowerAmplitudeByEngagement.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_FollowerAmplitudeByEngagement.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_FollowerAmplitudeByEngagement\nn(pos_nonengaged)=',num2str(length(rmmissing(this_data_pos_nonengaged))),'\nn(pos_engaged)=',num2str(length(rmmissing(this_data_pos_engaged))),'\nn(neg_nonengaged)=',num2str(length(rmmissing(this_data_neg_nonengaged))),'\nn(neg_engaged)=',num2str(length(rmmissing(this_data_neg_engaged))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_FollowerAmplitudeByEngagement\nn(pos_nonengaged)=',num2str(length(rmmissing(this_data_pos_nonengaged))),'\nn(pos_engaged)=',num2str(length(rmmissing(this_data_pos_engaged))),'\nn(neg_nonengaged)=',num2str(length(rmmissing(this_data_neg_nonengaged))),'\nn(neg_engaged)=',num2str(length(rmmissing(this_data_neg_engaged))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig4_FollowerAnalysis_PosFollowerAmplitude_2_NegFollowerAmplitude

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos = posFollowerAmplitude_ctrl_unpaired_sd1;
this_data_neg = negFollowerAmplitude_ctrl_unpaired_sd1;
scatter(this_data_pos,this_data_neg,'.','SizeData',100,'MarkerEdgeColor',p.col.darkGray)
[this_corr_r,this_corr_p] = fitLine(this_data_pos,this_data_neg,p.col.black);
xticks([0,0.25,0.5])
xlim([0,0.5])
xlabel({'Positive follower amplitude (S)'})
yticks([-0.5,-0.25,0])
ylim([-0.5,0])
ylabel({'Negative follower amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_PosFollowerAmplitude_2_NegFollowerAmplitude.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_PosFollowerAmplitude_2_NegFollowerAmplitude.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_PosFollowerAmplitude_2_NegFollowerAmplitude.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_PosFollowerAmplitude_2_NegFollowerAmplitude.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_PosFollowerAmplitude_2_NegFollowerAmplitude\nn(pos)=',num2str(length(rmmissing(this_data_pos))),'\nn(neg)=',num2str(length(rmmissing(this_data_neg))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_PosFollowerAmplitude_2_NegFollowerAmplitude\nn(pos)=',num2str(length(rmmissing(this_data_pos))),'\nn(neg)=',num2str(length(rmmissing(this_data_neg))),...
    '\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);


%% Fig4_FollowerAnalysis_FollowerAmplitudeProfile

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos = squeeze(nanmean(followerAmplitudeProfile_ctrl(:,[2,4],:),2));
plot(1:size(this_data_pos,2),this_data_pos','Color',mean([p.col.pos;p.col.white]),'LineWidth',0.5)
plot(1:size(this_data_pos,2),nanmean(this_data_pos',2),'Color',p.col.pos,'LineWidth',1.5)
xlim([1,size(this_data_pos,2)+1]) % xlim([1,size(this_data_neg,2)+1]-0.5)
% xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+1) % xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+0.5)
% xticklabels({'0','50','100','150','200','250','300','350'})
xticks([0:4:length([0:followerBinSize:maxDistanceFromClosestTarget])]+1) % xticks([0:2:length([0:followerBinSize:maxDistanceFromClosestTarget])]+0.5)
xticklabels({'0','100','200','300'})
% yticks([0,0.3,0.6])
% ylim([0,0.6])
% ytickformat('percentage')
xlabel({'Distance from closest';'laser beamlet (Um)'})
ylabel({'Response amplitude (S)'})

savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_FollowerAmplitudeProfile.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_FollowerAmplitudeProfile.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_FollowerAmplitudeProfile.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_FollowerAmplitudeProfile.txt'],'wt');
fprintf(fid,['\nFollowerAnalysis_FollowerAmplitudeProfile\nn=',num2str(size(rmmissing(this_data_pos,'MinNumMissing',10),1)),'\n']); fclose(fid);
fprintf(['\nFollowerAnalysis_FollowerAmplitudeProfile\nn=',num2str(size(rmmissing(this_data_pos,'MinNumMissing',10),1)),'\n']);


%% --- In progress ---









%% Fig4_FollowerAnalysis_StabilityWithinDay

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_pos = posFollowerAmplitudeStability_ctrl_unpaired_sd1;
this_data_neg = negFollowerAmplitudeStability_ctrl_unpaired_sd1;
shadedErrorBar(1:size(this_data_neg,2),nanmean(this_data_neg',2),nansem(this_data_neg',2),'lineProps',p.col.neg);
shadedErrorBar(1:size(this_data_pos,2),nanmean(this_data_pos',2),nansem(this_data_pos',2),'lineProps',p.col.pos);
yline(0,'k-');
xlim([1,size(this_data_pos,2)])
xticks([2:2:8])
xticklabels({'2','4','6','8'})
yticks([-1,0,1])
ylim([-1,1])
xlabel('Trial block')
ylabel({'Follower amplitude (sd)'})

temp1 = [this_data_pos;this_data_neg];
temp2 = repmat(1:8,size(temp1,1),1);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete');
savefig(F,[save_root_fig,'\Fig4_FollowerAnalysis_StabilityWithinDay.fig']);
saveas(F,[save_root_png,'\Fig4_FollowerAnalysis_StabilityWithinDay.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_FollowerAnalysis_StabilityWithinDay.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig4_FollowerAnalysis_StabilityWithinDay.txt'],'wt');
fprintf(fid,['\Fig4_FollowerAnalysis_StabilityWithinDay\nn(seq)=',num2str(size(rmmissing(this_data_pos),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_neg),1)),...
	'\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']); fclose(fid);
fprintf(['\Fig4_FollowerAnalysis_StabilityWithinDay\nn(seq)=',num2str(size(rmmissing(this_data_pos),1)),'\nn(ctrl)=',num2str(size(rmmissing(this_data_neg),1)),...
	'\nPearson rho=',num2str(this_corr_r,4),'\nPearson p=',num2str(this_corr_p,4),'\n']);








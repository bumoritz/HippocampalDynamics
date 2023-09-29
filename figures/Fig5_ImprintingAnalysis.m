%% Fig5_ImprintingAnalysis

% import data using Summary_Master with ops.do_imprSummary = true;
% d_info = selectDataset(d_info,'-g78-d2345-e01-r01-p1-l01-i00',sheet,path,ops); + one somehow needs imaging data as well

save_root_fig = [path.root_summary,'figures\Fig5_fig\'];
save_root_png = [path.root_summary,'figures\Fig5_png\'];
save_root_pdf = [path.root_summary,'figures\Fig5_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig5_txt\'];

resp_struct = 'resp';
tng_struct = 'tng_all_stimVersion';
tng_100t_struct = 'tng_100t_stimVersion';

followerBinSize = 25;
maxDistanceFromClosestTarget = 300;
posFollowers_skipFirstXbins = 3;
negFollowers_skipFirstXbins = 0;
cellTypeEdges = [0,1.3,5.3];
t_stim = 0.3:0.25:5.05;
these_animals_paired_firstDay = ~any(d_info.presponsive(:,[2,4])==0,2);
these_animals_paired_secondDay = ~any(d_info.presponsive(:,[3,5])==0,2);
these_animals_paired_bothDays = ~any(d_info.presponsive(:,2:5)==0,2);


%% Get context data

running_seq = nan(d_info.numAnimals,d_info.numDays);
running_ctrl = nan(d_info.numAnimals,d_info.numDays);
running_img = nan(d_info.numAnimals,d_info.numDays);
running_unpaired = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        running_unpaired(i,j) = d_info.running(i,j);
        if d_info.group(i)==2
            running_img(i,j) = d_info.running(i,j);
        end
        if d_info.presponsive(i,j)==1
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
running_img_unpaired_sd1 = nanmean(running_img(:,[2,4]),2);
running_img_unpaired_sd2 = nanmean(running_img(:,[3,5]),2);
running_seq_unpaired_sd1 = nanmean(running_seq(:,[2,4]),2);
running_seq_unpaired_sd2 = nanmean(running_seq(:,[3,5]),2);
running_seq_paired_sd1 = nanmean(running_seq(:,[2,4]),2);
running_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
running_seq_paired_sd2 = nanmean(running_seq(:,[3,5]),2);
running_seq_paired_sd2(~these_animals_paired_secondDay) = NaN;
running_ctrl_unpaired_sd1 = nanmean(running_ctrl(:,[2,4]),2);
running_ctrl_unpaired_sd2 = nanmean(running_ctrl(:,[3,5]),2);
running_ctrl_paired_sd1 = nanmean(running_ctrl(:,[2,4]),2);
running_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
running_ctrl_paired_sd2 = nanmean(running_ctrl(:,[3,5]),2);
running_ctrl_paired_sd2(~these_animals_paired_secondDay) = NaN;


%% Performance

correct = nan(d_info.numAnimals,d_info.numDays);
correct_img = nan(d_info.numAnimals,d_info.numDays);
correct_seq = nan(d_info.numAnimals,d_info.numDays);
correct_ctrl = nan(d_info.numAnimals,d_info.numDays);
correct_stim = nan(d_info.numAnimals,d_info.numDays);
correct_stim_img = nan(d_info.numAnimals,d_info.numDays);
correct_stim_seq = nan(d_info.numAnimals,d_info.numDays);
correct_stim_ctrl = nan(d_info.numAnimals,d_info.numDays);
correct_catch = nan(d_info.numAnimals,d_info.numDays);
correct_catch_img = nan(d_info.numAnimals,d_info.numDays);
correct_catch_seq = nan(d_info.numAnimals,d_info.numDays);
correct_catch_ctrl = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            correct(i,j) = d{i,j}.perf.general.correct; %nanmean(d{i,j}.perf.blocks_general.correct(1:5)); %d{i,j}.perf.general.correct;
            correct_stim(i,j) = d{i,j}.perf.general.correct_stim; %nanmean(d{i,j}.perf.blocks_general.correct_stim(1:5)); %d{i,j}.perf.general.correct_stim;
            correct_catch(i,j) = d{i,j}.perf.general.correct_catch; %nanmean(d{i,j}.perf.blocks_general.correct_catch(1:5)); %d{i,j}.perf.general.correct_catch;
            if d_info.group(i)==2
                correct_img(i,j) = correct(i,j);
                correct_stim_img(i,j) = correct_stim(i,j);
                correct_catch_img(i,j) = correct_catch(i,j);
            elseif d_info.group(i)==7 && (j==2 || j==3)
                correct_seq(i,j) = correct(i,j);
                correct_stim_seq(i,j) = correct_stim(i,j);
                correct_catch_seq(i,j) = correct_catch(i,j);
            elseif d_info.group(i)==8 && (j==4 || j==5)
                correct_seq(i,j) = correct(i,j);
                correct_stim_seq(i,j) = correct_stim(i,j);
                correct_catch_seq(i,j) = correct_catch(i,j);
            elseif d_info.group(i)==7 && (j==4 || j==5)
                correct_ctrl(i,j) = correct(i,j);
                correct_stim_ctrl(i,j) = correct_stim(i,j);
                correct_catch_ctrl(i,j) = correct_catch(i,j);
            elseif d_info.group(i)==8 && (j==2 || j==3)
                correct_ctrl(i,j) = correct(i,j);
                correct_stim_ctrl(i,j) = correct_stim(i,j);
                correct_catch_ctrl(i,j) = correct_catch(i,j);
            end
        catch
        end
    end
end


%% Get response and follower analysis data

responders = nan(d_info.numAnimals,d_info.numDays);
responders_seq = nan(d_info.numAnimals,d_info.numDays);
responders_ctrl = nan(d_info.numAnimals,d_info.numDays);
responders_unpaired = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                responders(i,j) = d{i,j}.(resp_struct).responders_main.numRespAll/40;
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
end
responders_unpaired_sd1 = nanmean(responders_unpaired(:,[2,4]),2); 
responders_seq_unpaired_sd1 = nanmean(responders_seq(:,[2,4]),2);
responders_seq_paired_sd1 = nanmean(responders_seq(:,[2,4]),2);
responders_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
responders_seq_paired_sd2 = nanmean(responders_seq(:,[3,5]),2);
responders_seq_paired_sd2(~these_animals_paired_secondDay) = NaN;
responders_ctrl_unpaired_sd1 = nanmean(responders_ctrl(:,[2,4]),2);
responders_ctrl_paired_sd1 = nanmean(responders_ctrl(:,[2,4]),2);
responders_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
responders_ctrl_paired_sd2 = nanmean(responders_ctrl(:,[3,5]),2);
responders_ctrl_paired_sd2(~these_animals_paired_secondDay) = NaN;
responders_seq_paired_sd12 = [nanmean(responders_seq(:,[2,4]),2),nanmean(responders_seq(:,[3,5]),2)];
responders_seq_paired_sd12(~these_animals_paired_bothDays) = NaN;
responders_ctrl_paired_sd12 = [nanmean(responders_ctrl(:,[2,4]),2),nanmean(responders_ctrl(:,[3,5]),2)];
responders_ctrl_paired_sd12(~these_animals_paired_bothDays) = NaN;

targets_unpaired = nan(d_info.numAnimals,d_info.numDays);
targets_seq = nan(d_info.numAnimals,d_info.numDays);
targets_ctrl = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                targets_unpaired(i,j) = d{i,j}.(resp_struct).responders_main.numRespTargeted/40; % / d{i,j}.(this_struct).responders_main.numIdentifiedTargeted;
                if d_info.group(i)==7 && (j==2 || j==3)
                    targets_seq(i,j) = targets_unpaired(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    targets_seq(i,j) = targets_unpaired(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    targets_ctrl(i,j) = targets_unpaired(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    targets_ctrl(i,j) = targets_unpaired(i,j);
                end
            catch
            end
        end
    end
end
targets_unpaired_sd1 = nanmean(targets_unpaired(:,[2,4]),2); 
targets_seq_unpaired_sd1 = nanmean(targets_seq(:,[2,4]),2);
targets_ctrl_unpaired_sd1 = nanmean(targets_ctrl(:,[2,4]),2);
targets_seq_paired_sd1 = nanmean(targets_seq(:,[2,4]),2);
targets_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
targets_seq_paired_sd2 = nanmean(targets_seq(:,[3,5]),2);
targets_seq_paired_sd2(~these_animals_paired_secondDay) = NaN;
targets_ctrl_paired_sd1 = nanmean(targets_ctrl(:,[2,4]),2);
targets_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
targets_ctrl_paired_sd2 = nanmean(targets_ctrl(:,[3,5]),2);
targets_ctrl_paired_sd2(~these_animals_paired_secondDay) = NaN;

amplitude = nan(d_info.numAnimals,d_info.numDays,8);
amplitude_unpaired = nan(d_info.numAnimals,d_info.numDays,8);
amplitude_seq = nan(d_info.numAnimals,d_info.numDays,8);
amplitude_ctrl = nan(d_info.numAnimals,d_info.numDays,8);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                temp0 = nanmean(d{i,j}.(resp_struct).respAmps_cells_bw,1); %nanmean(d{i,j}.(this_struct).trgRespAmps_cells_bw,1);
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
amplitude_seq_paired_sd2 = squeeze(nanmean(amplitude_seq(:,[3,5],:),2));
amplitude_seq_paired_sd2(~these_animals_paired_secondDay,:) = NaN;
amplitude_ctrl_unpaired_sd1 = squeeze(nanmean(amplitude_ctrl(:,[2,4],:),2));
amplitude_ctrl_paired_sd1 = squeeze(nanmean(amplitude_ctrl(:,[2,4],:),2));
amplitude_ctrl_paired_sd1(~these_animals_paired_firstDay,:) = NaN;
amplitude_ctrl_paired_sd2 = squeeze(nanmean(amplitude_ctrl(:,[3,5],:),2));
amplitude_ctrl_paired_sd2(~these_animals_paired_secondDay,:) = NaN;
amplitude_seq_paired_sd12 = cat(3,squeeze(nanmean(amplitude_seq(:,[2,4],:),2)),squeeze(nanmean(amplitude_seq(:,[3,5],:),2)));
amplitude_seq_paired_sd12(~these_animals_paired_bothDays,:,:) = NaN;
amplitude_ctrl_paired_sd12 = cat(3,squeeze(nanmean(amplitude_ctrl(:,[2,4],:),2)),squeeze(nanmean(amplitude_ctrl(:,[3,5],:),2)));
amplitude_ctrl_paired_sd12(~these_animals_paired_bothDays,:,:) = NaN;

negFollowers = nan(d_info.numAnimals,d_info.numDays);
negFollowers_seq = nan(d_info.numAnimals,d_info.numDays);
negFollowers_ctrl = nan(d_info.numAnimals,d_info.numDays);
negFollowers_unpaired = nan(d_info.numAnimals,d_info.numDays);
negFollowers_cw = nan(d_info.numAnimals,d_info.numDays,40);
idcs.followers.negFollowers = cell(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                idcs.followers.negFollowers{i,j} = [];
                for k=1:40
                    temp1 = d{i,j}.flw.stats_sig_neg(:,k);
                    temp2 = discretize(d{i,j}.resp.dist_closestLaser(:,k),[0:followerBinSize:maxDistanceFromClosestTarget]);
                    temp = nan(length(temp1),1);
                    temp(temp2>negFollowers_skipFirstXbins) = temp1(temp2>negFollowers_skipFirstXbins);
                    negFollowers_cw(i,j,k) = nansum(temp);
                    idcs.followers.negFollowers{i,j} = [idcs.followers.negFollowers{i,j}; find(temp==1)];
                end
                idcs.followers.negFollowers{i,j} = sort(unique(idcs.followers.negFollowers{i,j}));
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


%% Get cell structs

% idcs.iscells
temp = extractVariable(d,[tng_struct,'.prop.iscell'],'cell','all');
idcs.iscells = cell(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            idcs.iscells{i,j} = find(temp{i,j}==1);
        catch
        end
    end
end

% idcs.passed
these_conditions = {'Acatchonly','Xcatchonly','pref','none'};
for k=1:length(these_conditions)-2
    temp = extractVariable(d,[tng_struct,'.passed_catch.AW.',these_conditions{k}],'cell','all');
    idcs.passed.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            try
                idcs.passed.(these_conditions{k}){i,j} = find(temp{i,j}==1);
            catch
            end
        end
    end
end
idcs.passed.pref = cell(d_info.numAnimals,d_info.numDays);
idcs.passed.none = cell(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            idcs.passed.pref{i,j} = [idcs.passed.('Acatchonly'){i,j};idcs.passed.('Xcatchonly'){i,j}];
            idcs.passed.none{i,j} = setdiff(idcs.iscells{i,j},idcs.passed.pref{i,j});
        catch
        end
    end
end

% idcs.resp
these_conditions = {'A','X','O'};
for k=1:length(these_conditions)
    temp_rigid = extractVariable(d,['trg_rigid.grouping'],'cell','all');
    temp_nonrigid = extractVariable(d,['trg_nonrigid.grouping'],'cell','all');        
    temp = extractVariable(d,[resp_struct,'.responders_main.idcs_resp'],'cell','all');

    idcs.resp.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            try
                if strcmp(these_conditions{k},'A')
                    if ~isempty(temp_rigid{i,j})
                        temp0 = temp{i,j}(temp_rigid{i,j}(:,1));
                        idcs.resp.(these_conditions{k}){i,j} = unique(sort(vertcat(temp0{:})));
                    elseif ~isempty(temp_nonrigid{i,j})
                        temp0 = temp{i,j}(temp_nonrigid{i,j}(:,1));
                        idcs.resp.(these_conditions{k}){i,j} = unique(sort(vertcat(temp0{:})));
                    end
                elseif strcmp(these_conditions{k},'X')
                    if ~isempty(temp_rigid{i,j})
                        temp0 = temp{i,j}(temp_rigid{i,j}(:,2));
                        idcs.resp.(these_conditions{k}){i,j} = unique(sort(vertcat(temp0{:})));
                    elseif ~isempty(temp_nonrigid{i,j})
                        temp0 = temp{i,j}(temp_nonrigid{i,j}(:,2));
                        idcs.resp.(these_conditions{k}){i,j} = unique(sort(vertcat(temp0{:})));
                    end
                elseif strcmp(these_conditions{k},'O')
                    temp0 = temp{i,j};
                    idcs.resp.(these_conditions{k}){i,j} = setdiff(idcs.iscells{i,j},unique(sort(vertcat(temp0{:}))));
                end
            catch
            end
        end
    end
end

% idcs.respTrg
these_conditions = {'A','X','O'};
for k=1:length(these_conditions)
    temp_rigid = extractVariable(d,['trg_rigid.grouping'],'cell','all');
    temp_nonrigid = extractVariable(d,['trg_nonrigid.grouping'],'cell','all');        
    temp = extractVariable(d,[resp_struct,'.responders_main.idcs_respTargeted'],'cell','all');

    idcs.respTrg.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            try
                if strcmp(these_conditions{k},'A')
                    if ~isempty(temp_rigid{i,j})
                        temp0 = temp{i,j}(temp_rigid{i,j}(:,1));
                        idcs.respTrg.(these_conditions{k}){i,j} = unique(sort(vertcat(temp0{:})));
                    elseif ~isempty(temp_nonrigid{i,j})
                        temp0 = temp{i,j}(temp_nonrigid{i,j}(:,1));
                        idcs.respTrg.(these_conditions{k}){i,j} = unique(sort(vertcat(temp0{:})));
                    end
                elseif strcmp(these_conditions{k},'X')
                    if ~isempty(temp_rigid{i,j})
                        temp0 = temp{i,j}(temp_rigid{i,j}(:,2));
                        idcs.respTrg.(these_conditions{k}){i,j} = unique(sort(vertcat(temp0{:})));
                    elseif ~isempty(temp_nonrigid{i,j})
                        temp0 = temp{i,j}(temp_nonrigid{i,j}(:,2));
                        idcs.respTrg.(these_conditions{k}){i,j} = unique(sort(vertcat(temp0{:})));
                    end
                elseif strcmp(these_conditions{k},'O')
                    temp0 = temp{i,j};
                    idcs.respTrg.(these_conditions{k}){i,j} = setdiff(idcs.iscells{i,j},unique(sort(vertcat(temp0{:}))));
                end
            catch
            end
        end
    end
end

% idcs.recruited
these_conditions_1 = {'Acatchonly','Xcatchonly','pref','none'};
these_conditions_2 = {'A','X','O'};
for k=1:length(these_conditions_1)
    for m=1:length(these_conditions_2)
        idcs.recruited.([these_conditions_1{k},'_',these_conditions_2{m}]) = cell(d_info.numAnimals,d_info.numDays); 
        for i=1:d_info.numAnimals
            for j=1:d_info.numDays
                try
                    idcs.recruited.([these_conditions_1{k},'_',these_conditions_2{m}]){i,j} = intersect(idcs.passed.(these_conditions_1{k}){i,j},idcs.resp.(these_conditions_2{m}){i,j});
                catch
                end
            end
        end
    end
end

% idcs.recruited_trg
these_conditions_1 = {'Acatchonly','Xcatchonly','pref','none'};
these_conditions_2 = {'A','X','O'};
for k=1:length(these_conditions_1)
    for m=1:length(these_conditions_2)
        idcs.recruited_trg.([these_conditions_1{k},'_',these_conditions_2{m}]) = cell(d_info.numAnimals,d_info.numDays); 
        for i=1:d_info.numAnimals
            for j=1:d_info.numDays
                try
                    idcs.recruited_trg.([these_conditions_1{k},'_',these_conditions_2{m}]){i,j} = intersect(idcs.passed.(these_conditions_1{k}){i,j},idcs.respTrg.(these_conditions_2{m}){i,j});
                catch
                end
            end
        end
    end
end


%% Number of sequence cells

% numCells
numCells_img = nan(d_info.numAnimals,d_info.numDays);
numCells_seq = nan(d_info.numAnimals,d_info.numDays);
numCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
numCells_unpaired = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if ~isempty(idcs.iscells{i,j})
            numCells_unpaired(i,j) = length(idcs.iscells{i,j});
            try
                if d_info.group(i)==2
                    numCells_img(i,j) = numCells_unpaired(i,j);
                elseif d_info.group(i)==7 && (j==2 || j==3)
                    numCells_seq(i,j) = numCells_unpaired(i,j);
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    numCells_seq(i,j) = numCells_unpaired(i,j);
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    numCells_ctrl(i,j) = numCells_unpaired(i,j);
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    numCells_ctrl(i,j) = numCells_unpaired(i,j);
                end
            catch
            end
        end
    end
end
numCells_img_unpaired_sd1 = nanmean(numCells_img(:,[2,4]),2);
numCells_img_unpaired_sd2 = nanmean(numCells_img(:,[3,5]),2);
numCells_seq_unpaired_sd1 = nanmean(numCells_seq(:,[2,4]),2);
numCells_seq_unpaired_sd2 = nanmean(numCells_seq(:,[3,5]),2);
numCells_ctrl_unpaired_sd1 = nanmean(numCells_ctrl(:,[2,4]),2);
numCells_ctrl_unpaired_sd2 = nanmean(numCells_ctrl(:,[3,5]),2);
numCells_img_paired_sd1 = nanmean(numCells_img(:,[2,4]),2);
numCells_img_paired_sd1(~these_animals_paired_firstDay) = NaN;
numCells_img_paired_sd2 = nanmean(numCells_img(:,[3,5]),2);
numCells_img_paired_sd2(~these_animals_paired_secondDay) = NaN;
numCells_seq_paired_sd1 = nanmean(numCells_seq(:,[2,4]),2);
numCells_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
numCells_seq_paired_sd2 = nanmean(numCells_seq(:,[3,5]),2);
numCells_seq_paired_sd2(~these_animals_paired_secondDay) = NaN;
numCells_ctrl_paired_sd1 = nanmean(numCells_ctrl(:,[2,4]),2);
numCells_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
numCells_ctrl_paired_sd2 = nanmean(numCells_ctrl(:,[3,5]),2);
numCells_ctrl_paired_sd2(~these_animals_paired_secondDay) = NaN;

% numSeqCells, propSeqCells
these_fields = fields(idcs.passed);
for k=1:length(these_fields)
    numSeqCells_img.(these_fields{k}) = nan(d_info.numAnimals,d_info.numDays);
    numSeqCells_seq.(these_fields{k}) = nan(d_info.numAnimals,d_info.numDays);
    numSeqCells_ctrl.(these_fields{k}) = nan(d_info.numAnimals,d_info.numDays);
    numSeqCells_unpaired.(these_fields{k}) = nan(d_info.numAnimals,d_info.numDays);
    propSeqCells_img.(these_fields{k}) = nan(d_info.numAnimals,d_info.numDays);
    propSeqCells_seq.(these_fields{k}) = nan(d_info.numAnimals,d_info.numDays);
    propSeqCells_ctrl.(these_fields{k}) = nan(d_info.numAnimals,d_info.numDays);
    propSeqCells_unpaired.(these_fields{k}) = nan(d_info.numAnimals,d_info.numDays);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if ~isempty(idcs.iscells{i,j})
                numSeqCells_unpaired.(these_fields{k})(i,j) = length(idcs.passed.(these_fields{k}){i,j});
                propSeqCells_unpaired.(these_fields{k})(i,j) = numSeqCells_unpaired.(these_fields{k})(i,j) / numCells_unpaired(i,j);
                try
                    if d_info.group(i)==2
                        numSeqCells_img.(these_fields{k})(i,j) = numSeqCells_unpaired.(these_fields{k})(i,j);
                        propSeqCells_img.(these_fields{k})(i,j) = propSeqCells_unpaired.(these_fields{k})(i,j);
                    elseif d_info.group(i)==7 && (j==2 || j==3)
                        numSeqCells_seq.(these_fields{k})(i,j) = numSeqCells_unpaired.(these_fields{k})(i,j);
                        propSeqCells_seq.(these_fields{k})(i,j) = propSeqCells_unpaired.(these_fields{k})(i,j);
                    elseif d_info.group(i)==8 && (j==4 || j==5)
                        numSeqCells_seq.(these_fields{k})(i,j) = numSeqCells_unpaired.(these_fields{k})(i,j);
                        propSeqCells_seq.(these_fields{k})(i,j) = propSeqCells_unpaired.(these_fields{k})(i,j);
                    elseif d_info.group(i)==7 && (j==4 || j==5)
                        numSeqCells_ctrl.(these_fields{k})(i,j) = numSeqCells_unpaired.(these_fields{k})(i,j);
                        propSeqCells_ctrl.(these_fields{k})(i,j) = propSeqCells_unpaired.(these_fields{k})(i,j);
                    elseif d_info.group(i)==8 && (j==2 || j==3)
                        numSeqCells_ctrl.(these_fields{k})(i,j) = numSeqCells_unpaired.(these_fields{k})(i,j);
                        propSeqCells_ctrl.(these_fields{k})(i,j) = propSeqCells_unpaired.(these_fields{k})(i,j);
                    end
                catch
                end
            end
        end
    end
    numSeqCells_img_unpaired_sd1.(these_fields{k}) = nanmean(numSeqCells_img.(these_fields{k})(:,[2,4]),2);
    numSeqCells_img_unpaired_sd2.(these_fields{k}) = nanmean(numSeqCells_img.(these_fields{k})(:,[3,5]),2);
    numSeqCells_seq_unpaired_sd1.(these_fields{k}) = nanmean(numSeqCells_seq.(these_fields{k})(:,[2,4]),2);
    numSeqCells_seq_unpaired_sd2.(these_fields{k}) = nanmean(numSeqCells_seq.(these_fields{k})(:,[3,5]),2);
    numSeqCells_ctrl_unpaired_sd1.(these_fields{k}) = nanmean(numSeqCells_ctrl.(these_fields{k})(:,[2,4]),2);
    numSeqCells_ctrl_unpaired_sd2.(these_fields{k}) = nanmean(numSeqCells_ctrl.(these_fields{k})(:,[3,5]),2);
    propSeqCells_img_unpaired_sd1.(these_fields{k}) = nanmean(propSeqCells_img.(these_fields{k})(:,[2,4]),2);
    propSeqCells_img_unpaired_sd2.(these_fields{k}) = nanmean(propSeqCells_img.(these_fields{k})(:,[3,5]),2);
    propSeqCells_seq_unpaired_sd1.(these_fields{k}) = nanmean(propSeqCells_seq.(these_fields{k})(:,[2,4]),2);
    propSeqCells_seq_unpaired_sd2.(these_fields{k}) = nanmean(propSeqCells_seq.(these_fields{k})(:,[3,5]),2);
    propSeqCells_ctrl_unpaired_sd1.(these_fields{k}) = nanmean(propSeqCells_ctrl.(these_fields{k})(:,[2,4]),2);
    propSeqCells_ctrl_unpaired_sd2.(these_fields{k}) = nanmean(propSeqCells_ctrl.(these_fields{k})(:,[3,5]),2);
end
k=3;
propSeqCells_img_paired_sd1 = nanmean(propSeqCells_img.(these_fields{k})(:,[2,4]),2);
propSeqCells_img_paired_sd1(~these_animals_paired_firstDay) = NaN;
propSeqCells_img_paired_sd2 = nanmean(propSeqCells_img.(these_fields{k})(:,[3,5]),2);
propSeqCells_img_paired_sd2(~these_animals_paired_secondDay) = NaN;
propSeqCells_seq_paired_sd1 = nanmean(propSeqCells_seq.(these_fields{k})(:,[2,4]),2);
propSeqCells_seq_paired_sd1(~these_animals_paired_firstDay) = NaN;
propSeqCells_seq_paired_sd2 = nanmean(propSeqCells_seq.(these_fields{k})(:,[3,5]),2);
propSeqCells_seq_paired_sd2(~these_animals_paired_secondDay) = NaN;
propSeqCells_ctrl_paired_sd1 = nanmean(propSeqCells_ctrl.(these_fields{k})(:,[2,4]),2);
propSeqCells_ctrl_paired_sd1(~these_animals_paired_firstDay) = NaN;
propSeqCells_ctrl_paired_sd2 = nanmean(propSeqCells_ctrl.(these_fields{k})(:,[3,5]),2);
propSeqCells_ctrl_paired_sd2(~these_animals_paired_secondDay) = NaN;


%% Get firing field times

% stimTime_seq
temp_rigid = extractVariable(d,['trg_rigid.sequenceClusters'],'cell','all');
temp_nonrigid = extractVariable(d,['trg_nonrigid.sequenceClusters'],'cell','all');
these_conditions = fields(idcs.recruited);
for k=1:length(these_conditions)
    stimTime_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    %stimTime_all_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
%     stimTime_earlyIdcs_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
%     stimTime_lateIdcs_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    stimTime_cleanEarlyIdcs_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    stimTime_cleanLateIdcs_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    stimTime_cleanEarlyTimes_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    stimTime_cleanLateTimes_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if ~isempty(idcs.iscells{i,j})
                if ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
                    if isfield(d{i,j},'trg_rigid')
                        if strcmp(these_conditions{k},'Acatchonly_A')
                            temp = idcs.recruited.Acatchonly_A{i,j};
                            stimTime_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [];
                            for m=1:length(temp)
                                temp2 = find(d{i,j}.resp.responders_main.respAll(temp(m),:)==1);
                                if length(temp2)==1
                                    try
                                        stimTime_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,1)==temp2));
                                        if stimTime_seq.(these_conditions{k}){i,j}(m) < cellTypeEdges(2)
                                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        elseif stimTime_seq.(these_conditions{k}){i,j}(m) > cellTypeEdges(2)
                                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        end
                                    catch
                                    end
                                end
                            end
                        elseif strcmp(these_conditions{k},'Xcatchonly_X')
                            temp = idcs.recruited.Xcatchonly_X{i,j};
                            stimTime_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [];
                            for m=1:length(temp)
                                temp2 = find(d{i,j}.resp.responders_main.respAll(temp(m),:)==1);
                                if length(temp2)==1
                                    try
                                        stimTime_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,2)==temp2));
                                        if stimTime_seq.(these_conditions{k}){i,j}(m) < cellTypeEdges(2)
                                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        elseif stimTime_seq.(these_conditions{k}){i,j}(m) > cellTypeEdges(2)
                                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        end
                                    catch
                                    end
                                end
                            end
                        elseif strcmp(these_conditions{k},'Acatchonly_X')
                            temp = idcs.recruited.Acatchonly_X{i,j};
                            stimTime_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [];
                            for m=1:length(temp)
                                temp2 = find(d{i,j}.resp.responders_main.respAll(temp(m),:)==1);
                                if length(temp2)==1
                                    try
                                        stimTime_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,2)==temp2));
                                        if stimTime_seq.(these_conditions{k}){i,j}(m) < cellTypeEdges(2)
                                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        elseif stimTime_seq.(these_conditions{k}){i,j}(m) > cellTypeEdges(2)
                                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        end
                                    catch
                                    end
                                end
                            end
                        elseif strcmp(these_conditions{k},'Xcatchonly_A')
                            temp = idcs.recruited.Xcatchonly_A{i,j};
                            stimTime_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [];
                            for m=1:length(temp)
                                temp2 = find(d{i,j}.resp.responders_main.respAll(temp(m),:)==1);
                                if length(temp2)==1
                                    try
                                        stimTime_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,1)==temp2));
                                        if stimTime_seq.(these_conditions{k}){i,j}(m) < cellTypeEdges(2)
                                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        elseif stimTime_seq.(these_conditions{k}){i,j}(m) > cellTypeEdges(2)
                                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        end
                                    catch
                                    end
                                end
                            end
                        end
                    elseif isfield(d{i,j},'trg_nonrigid')
                        if strcmp(these_conditions{k},'Acatchonly_A')
                            temp = idcs.recruited.Acatchonly_A{i,j};
                            stimTime_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [];
                            for m=1:length(temp)
                                temp2 = find(d{i,j}.resp.responders_main.respAll(temp(m),:)==1);
                                if length(temp2)==1
                                    try
                                        stimTime_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,1)==temp2));
                                        if stimTime_seq.(these_conditions{k}){i,j}(m) < cellTypeEdges(2)
                                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        elseif stimTime_seq.(these_conditions{k}){i,j}(m) > cellTypeEdges(2)
                                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        end
                                    catch
                                    end
                                end
                            end
                        elseif strcmp(these_conditions{k},'Xcatchonly_X')
                            temp = idcs.recruited.Xcatchonly_X{i,j};
                            stimTime_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [];
                            for m=1:length(temp)
                                temp2 = find(d{i,j}.resp.responders_main.respAll(temp(m),:)==1);
                                if length(temp2)==1
                                    try
                                        stimTime_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,2)==temp2));
                                        if stimTime_seq.(these_conditions{k}){i,j}(m) < cellTypeEdges(2)
                                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        elseif stimTime_seq.(these_conditions{k}){i,j}(m) > cellTypeEdges(2)
                                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        end
                                    catch
                                    end
                                end
                            end
                        elseif strcmp(these_conditions{k},'Acatchonly_X')
                            temp = idcs.recruited.Acatchonly_X{i,j};
                            stimTime_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [];
                            for m=1:length(temp)
                                temp2 = find(d{i,j}.resp.responders_main.respAll(temp(m),:)==1);
                                if length(temp2)==1
                                    try
                                        stimTime_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,2)==temp2));
                                        if stimTime_seq.(these_conditions{k}){i,j}(m) < cellTypeEdges(2)
                                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        elseif stimTime_seq.(these_conditions{k}){i,j}(m) > cellTypeEdges(2)
                                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        end
                                    catch
                                    end
                                end
                            end
                        elseif strcmp(these_conditions{k},'Xcatchonly_A')
                            temp = idcs.recruited.Xcatchonly_A{i,j};
                            stimTime_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [];
                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [];
                            for m=1:length(temp)
                                temp2 = find(d{i,j}.resp.responders_main.respAll(temp(m),:)==1);
                                if length(temp2)==1
                                    try
                                        stimTime_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,1)==temp2));
                                        if stimTime_seq.(these_conditions{k}){i,j}(m) < cellTypeEdges(2)
                                            stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanEarlyTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        elseif stimTime_seq.(these_conditions{k}){i,j}(m) > cellTypeEdges(2)
                                            stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateIdcs_seq.(these_conditions{k}){i,j}; temp(m)];
                                            stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j} = [stimTime_cleanLateTimes_seq.(these_conditions{k}){i,j}; stimTime_seq.(these_conditions{k}){i,j}(m)];
                                        end
                                    catch
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% stimTime_trg_seq
temp_rigid = extractVariable(d,['trg_rigid.sequenceClusters'],'cell','all');
temp_nonrigid = extractVariable(d,['trg_nonrigid.sequenceClusters'],'cell','all');
these_conditions = fields(idcs.recruited_trg);
t_stim = 0.3:0.25:5.05;
for k=1:length(these_conditions)
    stimTime_trg_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if ~isempty(idcs.iscells{i,j})
                if ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
                    if isfield(d{i,j},'trg_rigid')
                        if strcmp(these_conditions{k},'Acatchonly_A')
                            temp = idcs.recruited_trg.Acatchonly_A{i,j};
                            stimTime_trg_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            for m=1:length(temp)
                                try
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,1)==min(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                catch
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,1)==max(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                end
                            end
                        elseif strcmp(these_conditions{k},'Xcatchonly_X')
                            temp = idcs.recruited_trg.Xcatchonly_X{i,j};
                            stimTime_trg_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            for m=1:length(temp)
                                try
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,2)==min(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                catch
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,2)==max(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                end
                            end
                        elseif strcmp(these_conditions{k},'Acatchonly_X')
                            temp = idcs.recruited_trg.Acatchonly_X{i,j};
                            stimTime_trg_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            for m=1:length(temp)
                                try
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,2)==min(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                catch
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,2)==max(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                end
                            end
                        elseif strcmp(these_conditions{k},'Xcatchonly_A')
                            temp = idcs.recruited_trg.Xcatchonly_A{i,j};
                            stimTime_trg_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            for m=1:length(temp)
                                try
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,1)==min(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                catch
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_rigid.sequenceClusters(:,1)==max(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                end
                            end
                        end
                    elseif isfield(d{i,j},'trg_nonrigid')
                        if strcmp(these_conditions{k},'Acatchonly_A')
                            temp = idcs.recruited_trg.Acatchonly_A{i,j};
                            stimTime_trg_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            for m=1:length(temp)
                                % problem i 26, j 2, m 4: cell should be
                                % part of cluster 8, but cluster 8 seems to
                                % be from X session
                                try
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,1)==min(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                catch
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,1)==max(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                end
                            end
                        elseif strcmp(these_conditions{k},'Xcatchonly_X')
                            temp = idcs.recruited_trg.Xcatchonly_X{i,j};
                            stimTime_trg_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            for m=1:length(temp)
                                try
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,2)==min(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                    % stimTime_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,2)==find(nanmean(d{i,j}.trg_nonrigid.idcs_targetedCells==temp(m),1))));
                                catch
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,2)==max(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                end
                            end
                        elseif strcmp(these_conditions{k},'Acatchonly_X')
                            temp = idcs.recruited_trg.Acatchonly_X{i,j};
                            stimTime_trg_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            for m=1:length(temp)
                                try
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,2)==min(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                catch
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,2)==max(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                end
                            end
                        elseif strcmp(these_conditions{k},'Xcatchonly_A')
                            temp = idcs.recruited_trg.Xcatchonly_A{i,j};
                            stimTime_trg_seq.(these_conditions{k}){i,j} = nan(length(temp),1);
                            for m=1:length(temp)
                                try
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,1)==min(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                catch
                                    stimTime_trg_seq.(these_conditions{k}){i,j}(m) = t_stim(find(d{i,j}.trg_nonrigid.sequenceClusters(:,1)==max(find(nanmean(d{i,j}.resp.responders_main.targetedCells==temp(m),1)))));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% peakTime
temp_A = extractVariable(d,[tng_struct,'.firingField.Acatch_AW.peakLocation_s'],'cell','all');
temp_X = extractVariable(d,[tng_struct,'.firingField.Xcatch_AW.peakLocation_s'],'cell','all');
these_conditions = [fields(idcs.recruited);'Acatchonly';'Xcatchonly';'AcatchonlyResp';'XcatchonlyResp'];
for k=1:length(these_conditions)
    peakTime.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    peakTime_img.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    peakTime_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    peakTime_ctrl.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if ~isempty(idcs.iscells{i,j})
                if strcmp(these_conditions{k},'Acatchonly')
                    peakTime.(these_conditions{k}){i,j} = temp_A{i,j}(idcs.passed.Acatchonly{i,j});
                elseif strcmp(these_conditions{k},'Xcatchonly')
                    peakTime.(these_conditions{k}){i,j} = temp_X{i,j}(idcs.passed.Xcatchonly{i,j});
                elseif strcmp(these_conditions{k},'AcatchonlyResp')
                    peakTime.(these_conditions{k}){i,j} = temp_A{i,j}(intersect(idcs.passed.Acatchonly{i,j},unique([idcs.resp.A{i,j};idcs.resp.X{i,j}])));
                elseif strcmp(these_conditions{k},'XcatchonlyResp')
                    peakTime.(these_conditions{k}){i,j} = temp_X{i,j}(intersect(idcs.passed.Xcatchonly{i,j},unique([idcs.resp.A{i,j};idcs.resp.X{i,j}])));
                elseif strcmp(these_conditions{k},'Acatchonly_A')
                    peakTime.(these_conditions{k}){i,j} = temp_A{i,j}(idcs.recruited.Acatchonly_A{i,j});
                elseif strcmp(these_conditions{k},'Xcatchonly_X')
                    peakTime.(these_conditions{k}){i,j} = temp_X{i,j}(idcs.recruited.Xcatchonly_X{i,j});
                elseif strcmp(these_conditions{k},'Acatchonly_X')
                    peakTime.(these_conditions{k}){i,j} = temp_A{i,j}(idcs.recruited.Acatchonly_X{i,j});
                elseif strcmp(these_conditions{k},'Xcatchonly_A')
                    peakTime.(these_conditions{k}){i,j} = temp_X{i,j}(idcs.recruited.Xcatchonly_A{i,j});
                end
                if d_info.group(i)==2
                    peakTime_img.(these_conditions{k}){i,j} = peakTime.(these_conditions{k}){i,j};
                elseif d_info.group(i)==7 && (j==2 || j==3)
                    peakTime_seq.(these_conditions{k}){i,j} = peakTime.(these_conditions{k}){i,j};
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    peakTime_seq.(these_conditions{k}){i,j} = peakTime.(these_conditions{k}){i,j};
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    peakTime_ctrl.(these_conditions{k}){i,j} = peakTime.(these_conditions{k}){i,j};
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    peakTime_ctrl.(these_conditions{k}){i,j} = peakTime.(these_conditions{k}){i,j};
                end
            end
        end
    end
end

% peakTime_trg
temp_A = extractVariable(d,[tng_struct,'.firingField.Acatch_AW.peakLocation_s'],'cell','all');
temp_X = extractVariable(d,[tng_struct,'.firingField.Xcatch_AW.peakLocation_s'],'cell','all');
these_conditions = [fields(idcs.recruited);'Acatchonly';'Xcatchonly';'AcatchonlyResp';'XcatchonlyResp'];
for k=1:length(these_conditions)
    peakTime_trg.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    peakTime_trg_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    peakTime_trg_ctrl.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if ~isempty(idcs.iscells{i,j})
                if strcmp(these_conditions{k},'Acatchonly')
                    peakTime.(these_conditions{k}){i,j} = temp_A{i,j}(idcs.passed.Acatchonly{i,j});
                elseif strcmp(these_conditions{k},'Xcatchonly')
                    peakTime.(these_conditions{k}){i,j} = temp_X{i,j}(idcs.passed.Xcatchonly{i,j});
                elseif strcmp(these_conditions{k},'AcatchonlyResp')
                    peakTime.(these_conditions{k}){i,j} = temp_A{i,j}(intersect(idcs.passed.Acatchonly{i,j},unique([idcs.respTrg.A{i,j};idcs.respTrg.X{i,j}])));
                elseif strcmp(these_conditions{k},'XcatchonlyResp')
                    peakTime.(these_conditions{k}){i,j} = temp_X{i,j}(intersect(idcs.passed.Xcatchonly{i,j},unique([idcs.respTrg.A{i,j};idcs.respTrg.X{i,j}])));
                elseif strcmp(these_conditions{k},'Acatchonly_A')
                    peakTime_trg.(these_conditions{k}){i,j} = temp_A{i,j}(idcs.recruited_trg.Acatchonly_A{i,j});
                elseif strcmp(these_conditions{k},'Xcatchonly_X')
                    peakTime_trg.(these_conditions{k}){i,j} = temp_X{i,j}(idcs.recruited_trg.Xcatchonly_X{i,j});
                elseif strcmp(these_conditions{k},'Acatchonly_X')
                    peakTime_trg.(these_conditions{k}){i,j} = temp_A{i,j}(idcs.recruited_trg.Acatchonly_X{i,j});
                elseif strcmp(these_conditions{k},'Xcatchonly_A')
                    peakTime_trg.(these_conditions{k}){i,j} = temp_X{i,j}(idcs.recruited_trg.Xcatchonly_A{i,j});
                end
                if d_info.group(i)==7 && (j==2 || j==3)
                    peakTime_trg_seq.(these_conditions{k}){i,j} = peakTime_trg.(these_conditions{k}){i,j};
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    peakTime_trg_seq.(these_conditions{k}){i,j} = peakTime_trg.(these_conditions{k}){i,j};
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    peakTime_trg_ctrl.(these_conditions{k}){i,j} = peakTime_trg.(these_conditions{k}){i,j};
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    peakTime_trg_ctrl.(these_conditions{k}){i,j} = peakTime_trg.(these_conditions{k}){i,j};
                end
            end
        end
    end
end

% numByCellType
these_fields = fields(peakTime);
for k=1:length(these_fields)
    numByCellType.(these_fields{k}) = nan(d_info.numAnimals,d_info.numDays,length(cellTypeEdges)-1);
    idcsByCellType.(these_fields{k}) = cell(d_info.numAnimals,d_info.numDays,length(cellTypeEdges)-1);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if ~isempty(idcs.iscells{i,j})
                temp = discretize(peakTime.(these_fields{k}){i,j},cellTypeEdges);
                for n=1:length(cellTypeEdges)-1
                    numByCellType.(these_fields{k})(i,j,n) = nansum(temp==n);
                end
                
                % added to get idcsByCellType
                if (strcmp(these_conditions{k},'Acatchonly') || strcmp(these_conditions{k},'Xcatchonly'))                   
                    idcsByCellType.(these_fields{k}){i,j,1} = idcs.passed.(these_fields{k}){i,j}(find(peakTime.(these_fields{k}){i,j} <= cellTypeEdges(2)));
                    idcsByCellType.(these_fields{k}){i,j,2} = idcs.passed.(these_fields{k}){i,j}(find(peakTime.(these_fields{k}){i,j} > cellTypeEdges(2)));
                end
            end
        end
    end
end

% numByCellType_trg
these_fields = fields(peakTime_trg);
for k=1:length(these_fields)
    numByCellType_trg.(these_fields{k}) = nan(d_info.numAnimals,d_info.numDays,length(cellTypeEdges)-1);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if ~isempty(idcs.iscells{i,j})
                temp = discretize(peakTime_trg.(these_fields{k}){i,j},cellTypeEdges);
                for n=1:length(cellTypeEdges)-1
                    numByCellType_trg.(these_fields{k})(i,j,n) = nansum(temp==n);
                end
            end
        end
    end
end

% comTime
temp_A = extractVariable(d,[tng_struct,'.firingField.Acatch_AW.comLocation_s'],'cell','all');
temp_X = extractVariable(d,[tng_struct,'.firingField.Xcatch_AW.comLocation_s'],'cell','all');
these_conditions = fields(idcs.recruited);
for k=1:length(these_conditions)
    comTime.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    comTime_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    comTime_ctrl.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if ~isempty(idcs.iscells{i,j})
                if strcmp(these_conditions{k},'Acatchonly_A')
                    comTime.(these_conditions{k}){i,j} = temp_A{i,j}(idcs.recruited.Acatchonly_A{i,j});
                elseif strcmp(these_conditions{k},'Xcatchonly_X')
                    comTime.(these_conditions{k}){i,j} = temp_X{i,j}(idcs.recruited.Xcatchonly_X{i,j});
                elseif strcmp(these_conditions{k},'Acatchonly_X')
                    comTime.(these_conditions{k}){i,j} = temp_A{i,j}(idcs.recruited.Acatchonly_X{i,j});
                elseif strcmp(these_conditions{k},'Xcatchonly_A')
                    comTime.(these_conditions{k}){i,j} = temp_X{i,j}(idcs.recruited.Xcatchonly_A{i,j});
                end
                if d_info.group(i)==7 && (j==2 || j==3)
                    comTime_seq.(these_conditions{k}){i,j} = comTime.(these_conditions{k}){i,j};
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    comTime_seq.(these_conditions{k}){i,j} = comTime.(these_conditions{k}){i,j};
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    comTime_ctrl.(these_conditions{k}){i,j} = comTime.(these_conditions{k}){i,j};
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    comTime_ctrl.(these_conditions{k}){i,j} = comTime.(these_conditions{k}){i,j};
                end
            end
        end
    end
end

% comTime_trg
temp_A = extractVariable(d,[tng_struct,'.firingField.Acatch_AW.comLocation_s'],'cell','all');
temp_X = extractVariable(d,[tng_struct,'.firingField.Xcatch_AW.comLocation_s'],'cell','all');
these_conditions = fields(idcs.recruited);
for k=1:length(these_conditions)
    comTime_trg.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    comTime_trg_seq.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    comTime_trg_ctrl.(these_conditions{k}) = cell(d_info.numAnimals,d_info.numDays);
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if ~isempty(idcs.iscells{i,j})
                if strcmp(these_conditions{k},'Acatchonly_A')
                    comTime_trg.(these_conditions{k}){i,j} = temp_A{i,j}(idcs.recruited_trg.Acatchonly_A{i,j});
                elseif strcmp(these_conditions{k},'Xcatchonly_X')
                    comTime_trg.(these_conditions{k}){i,j} = temp_X{i,j}(idcs.recruited_trg.Xcatchonly_X{i,j});
                elseif strcmp(these_conditions{k},'Acatchonly_X')
                    comTime_trg.(these_conditions{k}){i,j} = temp_A{i,j}(idcs.recruited_trg.Acatchonly_X{i,j});
                elseif strcmp(these_conditions{k},'Xcatchonly_A')
                    comTime_trg.(these_conditions{k}){i,j} = temp_X{i,j}(idcs.recruited_trg.Xcatchonly_A{i,j});
                end
                if d_info.group(i)==7 && (j==2 || j==3)
                    comTime_trg_seq.(these_conditions{k}){i,j} = comTime_trg.(these_conditions{k}){i,j};
                elseif d_info.group(i)==8 && (j==4 || j==5)
                    comTime_trg_seq.(these_conditions{k}){i,j} = comTime_trg.(these_conditions{k}){i,j};
                elseif d_info.group(i)==7 && (j==4 || j==5)
                    comTime_trg_ctrl.(these_conditions{k}){i,j} = comTime_trg.(these_conditions{k}){i,j};
                elseif d_info.group(i)==8 && (j==2 || j==3)
                    comTime_trg_ctrl.(these_conditions{k}){i,j} = comTime_trg.(these_conditions{k}){i,j};
                end
            end
        end
    end
end


%% AvgTraces

% avgTraces    
avgTraces.A_var0 = cell(d_info.numAnimals,d_info.numDays);
avgTraces.A_var1 = cell(d_info.numAnimals,d_info.numDays);
avgTraces.X_var0 = cell(d_info.numAnimals,d_info.numDays);
avgTraces.X_var1 = cell(d_info.numAnimals,d_info.numDays);
avgTraces.nft_bl_avg_binned = cell(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            if ~isempty(idcs.iscells{i,j})
                avgTraces.A_var0{i,j} = d{i,j}.str.avgTraces_all_stimVersion.A_var0;
                avgTraces.A_var1{i,j} = d{i,j}.str.avgTraces_all_stimVersion.A_var1;
                avgTraces.X_var0{i,j} = d{i,j}.str.avgTraces_all_stimVersion.X_var0;
                avgTraces.X_var1{i,j} = d{i,j}.str.avgTraces_all_stimVersion.X_var1;
                avgTraces.nft_bl_avg_binned{i,j} = d{i,j}.str.nft_bl_avg_binned;
            end
        catch
        end
    end
end

% avgTraces_recruited
avgTraces_recruited.A_var0__Acatchonly_A = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.A_var0__Acatchonly_X = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.A_var0__Xcatchonly_X = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.A_var0__Xcatchonly_A = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.A_var1__Acatchonly_A = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.A_var1__Acatchonly_X = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.A_var1__Xcatchonly_X = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.A_var1__Xcatchonly_A = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.X_var0__Acatchonly_A = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.X_var0__Acatchonly_X = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.X_var0__Xcatchonly_X = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.X_var0__Xcatchonly_A = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.X_var1__Acatchonly_A = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.X_var1__Acatchonly_X = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.X_var1__Xcatchonly_X = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.X_var1__Xcatchonly_A = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.nft_bl_avg_binned__Acatchonly_A = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.nft_bl_avg_binned__Acatchonly_X = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.nft_bl_avg_binned__Xcatchonly_X = cell(d_info.numAnimals,d_info.numDays);
avgTraces_recruited.nft_bl_avg_binned__Xcatchonly_A = cell(d_info.numAnimals,d_info.numDays);
these_fields = fields(avgTraces_recruited);
for k=1:length(these_fields)
    avgTraces_recruited_seq.(these_fields{k}) = cell(d_info.numAnimals,d_info.numDays);
end
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            avgTraces_recruited.A_var0__Acatchonly_A{i,j} = avgTraces.A_var0{i,j}(idcs.recruited.Acatchonly_A{i,j},:);
            avgTraces_recruited.A_var0__Acatchonly_X{i,j} = avgTraces.A_var0{i,j}(idcs.recruited.Acatchonly_X{i,j},:);
            avgTraces_recruited.A_var0__Xcatchonly_X{i,j} = avgTraces.A_var0{i,j}(idcs.recruited.Xcatchonly_X{i,j},:);
            avgTraces_recruited.A_var0__Xcatchonly_A{i,j} = avgTraces.A_var0{i,j}(idcs.recruited.Xcatchonly_A{i,j},:);
            avgTraces_recruited.A_var1__Acatchonly_A{i,j} = avgTraces.A_var1{i,j}(idcs.recruited.Acatchonly_A{i,j},:);
            avgTraces_recruited.A_var1__Acatchonly_X{i,j} = avgTraces.A_var1{i,j}(idcs.recruited.Acatchonly_X{i,j},:);
            avgTraces_recruited.A_var1__Xcatchonly_X{i,j} = avgTraces.A_var1{i,j}(idcs.recruited.Xcatchonly_X{i,j},:);
            avgTraces_recruited.A_var1__Xcatchonly_A{i,j} = avgTraces.A_var1{i,j}(idcs.recruited.Xcatchonly_A{i,j},:);
            avgTraces_recruited.X_var0__Acatchonly_A{i,j} = avgTraces.X_var0{i,j}(idcs.recruited.Acatchonly_A{i,j},:);
            avgTraces_recruited.X_var0__Acatchonly_X{i,j} = avgTraces.X_var0{i,j}(idcs.recruited.Acatchonly_X{i,j},:);
            avgTraces_recruited.X_var0__Xcatchonly_X{i,j} = avgTraces.X_var0{i,j}(idcs.recruited.Xcatchonly_X{i,j},:);
            avgTraces_recruited.X_var0__Xcatchonly_A{i,j} = avgTraces.X_var0{i,j}(idcs.recruited.Xcatchonly_A{i,j},:);
            avgTraces_recruited.X_var1__Acatchonly_A{i,j} = avgTraces.X_var1{i,j}(idcs.recruited.Acatchonly_A{i,j},:);
            avgTraces_recruited.X_var1__Acatchonly_X{i,j} = avgTraces.X_var1{i,j}(idcs.recruited.Acatchonly_X{i,j},:);
            avgTraces_recruited.X_var1__Xcatchonly_X{i,j} = avgTraces.X_var1{i,j}(idcs.recruited.Xcatchonly_X{i,j},:);
            avgTraces_recruited.X_var1__Xcatchonly_A{i,j} = avgTraces.X_var1{i,j}(idcs.recruited.Xcatchonly_A{i,j},:);
            avgTraces_recruited.nft_bl_avg_binned__Acatchonly_A{i,j} = avgTraces.nft_bl_avg_binned{i,j}(idcs.recruited.Acatchonly_A{i,j},:);
            avgTraces_recruited.nft_bl_avg_binned__Acatchonly_X{i,j} = avgTraces.nft_bl_avg_binned{i,j}(idcs.recruited.Acatchonly_X{i,j},:);
            avgTraces_recruited.nft_bl_avg_binned__Xcatchonly_X{i,j} = avgTraces.nft_bl_avg_binned{i,j}(idcs.recruited.Xcatchonly_X{i,j},:);
            avgTraces_recruited.nft_bl_avg_binned__Xcatchonly_A{i,j} = avgTraces.nft_bl_avg_binned{i,j}(idcs.recruited.Xcatchonly_A{i,j},:);
            
            these_fields = fields(avgTraces_recruited);
            if ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
                for k=1:length(these_fields)
                    avgTraces_recruited_seq.(these_fields{k}){i,j} = avgTraces_recruited.(these_fields{k}){i,j};
                end
            end
        catch
        end
    end
end


%% Properties

% activationProbability
activationProbabilityPhotoactivated_img = nan(d_info.numAnimals,d_info.numDays);
activationProbabilityPhotoactivated_seq = nan(d_info.numAnimals,d_info.numDays);
activationProbabilityPhotoactivated_ctrl = nan(d_info.numAnimals,d_info.numDays);
activationProbabilityPhotoactivated = nan(d_info.numAnimals,d_info.numDays);
activationProbabilityOther_img = nan(d_info.numAnimals,d_info.numDays);
activationProbabilityOther_seq = nan(d_info.numAnimals,d_info.numDays);
activationProbabilityOther_ctrl = nan(d_info.numAnimals,d_info.numDays);
activationProbabilityOther = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            % calculate activation probability
            tempA = []; tempX = [];
            if ~isempty(idcs.recruited.Acatchonly_A{i,j})
                tempA = d{i,j}.(tng_struct).firingField.Acatch_AW.activationProbability_ipsi(idcs.recruited.Acatchonly_A{i,j});
            end
            if ~isempty(idcs.recruited.Xcatchonly_X{i,j})
                tempX = d{i,j}.(tng_struct).firingField.Xcatch_AW.activationProbability_ipsi(idcs.recruited.Xcatchonly_X{i,j});
            end
            activationProbabilityPhotoactivated(i,j) = nanmean([tempA;tempX]);
            tempA = []; tempX = [];
            if ~isempty(idcs.passed.Acatchonly{i,j})
                tempA = d{i,j}.(tng_struct).firingField.Acatch_AW.activationProbability_ipsi(setdiff(idcs.passed.Acatchonly{i,j},idcs.recruited.Acatchonly_A{i,j}));
            end
            if ~isempty(idcs.passed.Xcatchonly{i,j})
                tempX = d{i,j}.(tng_struct).firingField.Xcatch_AW.activationProbability_ipsi(setdiff(idcs.passed.Xcatchonly{i,j},idcs.recruited.Xcatchonly_X{i,j}));
            end
            activationProbabilityOther(i,j)= nanmean([tempA;tempX]);
            
            % split by seq vs ctrl
            if d_info.group(i)==2
                activationProbabilityPhotoactivated_img(i,j) = activationProbabilityPhotoactivated(i,j);
                activationProbabilityOther_img(i,j) = activationProbabilityOther(i,j);
            elseif((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
                activationProbabilityPhotoactivated_seq(i,j) = activationProbabilityPhotoactivated(i,j);
                activationProbabilityOther_seq(i,j) = activationProbabilityOther(i,j);
            elseif ((d_info.group(i)==7 && (j==4 || j==5)) || (d_info.group(i)==8 && (j==2 || j==3)))
                activationProbabilityPhotoactivated_ctrl(i,j) = activationProbabilityPhotoactivated(i,j);
                activationProbabilityOther_ctrl(i,j) = activationProbabilityOther(i,j);
            end
        catch
        end
    end
end

% selectivity
selectivityPhotoactivated_img = nan(d_info.numAnimals,d_info.numDays);
selectivityPhotoactivated_seq = nan(d_info.numAnimals,d_info.numDays);
selectivityPhotoactivated_ctrl = nan(d_info.numAnimals,d_info.numDays);
selectivityPhotoactivated = nan(d_info.numAnimals,d_info.numDays);
selectivityOther_img = nan(d_info.numAnimals,d_info.numDays);
selectivityOther_seq = nan(d_info.numAnimals,d_info.numDays);
selectivityOther_ctrl = nan(d_info.numAnimals,d_info.numDays);
selectivityOther = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            % calculate selectivity
            tempA = []; tempX = [];
            if ~isempty(idcs.recruited.Acatchonly_A{i,j})
                tempA = d{i,j}.(tng_struct).firingField.Acatch_AW.selectivity_ipsi(idcs.recruited.Acatchonly_A{i,j});
            end
            if ~isempty(idcs.recruited.Xcatchonly_X{i,j})
                tempX = d{i,j}.(tng_struct).firingField.Xcatch_AW.selectivity_ipsi(idcs.recruited.Xcatchonly_X{i,j});
            end
            selectivityPhotoactivated(i,j) = nanmean([tempA;tempX]);
            tempA = []; tempX = [];
            if ~isempty(idcs.passed.Acatchonly{i,j})
                tempA = d{i,j}.(tng_struct).firingField.Acatch_AW.selectivity_ipsi(setdiff(idcs.passed.Acatchonly{i,j},idcs.recruited.Acatchonly_A{i,j}));
            end
            if ~isempty(idcs.passed.Xcatchonly{i,j})
                tempX = d{i,j}.(tng_struct).firingField.Xcatch_AW.selectivity_ipsi(setdiff(idcs.passed.Xcatchonly{i,j},idcs.recruited.Xcatchonly_X{i,j}));
            end
            selectivityOther(i,j)= nanmean([tempA;tempX]);
            
            % split by seq vs ctrl
            if d_info.group(i)==2
                selectivityPhotoactivated_img(i,j) = selectivityPhotoactivated(i,j);
                selectivityOther_img(i,j) = selectivityOther(i,j);
            elseif((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
                selectivityPhotoactivated_seq(i,j) = selectivityPhotoactivated(i,j);
                selectivityOther_seq(i,j) = selectivityOther(i,j);
            elseif ((d_info.group(i)==7 && (j==4 || j==5)) || (d_info.group(i)==8 && (j==2 || j==3)))
                selectivityPhotoactivated_ctrl(i,j) = selectivityPhotoactivated(i,j);
                selectivityOther_ctrl(i,j) = selectivityOther(i,j);
            end
        catch
        end
    end
end

% meanAmplitudeActive
meanAmplitudeActivePhotoactivated_img = nan(d_info.numAnimals,d_info.numDays);
meanAmplitudeActivePhotoactivated_seq = nan(d_info.numAnimals,d_info.numDays);
meanAmplitudeActivePhotoactivated_ctrl = nan(d_info.numAnimals,d_info.numDays);
meanAmplitudeActivePhotoactivated = nan(d_info.numAnimals,d_info.numDays);
meanAmplitudeActiveOther_img = nan(d_info.numAnimals,d_info.numDays);
meanAmplitudeActiveOther_seq = nan(d_info.numAnimals,d_info.numDays);
meanAmplitudeActiveOther_ctrl = nan(d_info.numAnimals,d_info.numDays);
meanAmplitudeActiveOther = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            % calculate meanAmplitudeActive
            tempA = []; tempX = [];
            if ~isempty(idcs.recruited.Acatchonly_A{i,j})
                tempA = d{i,j}.(tng_struct).firingField.Acatch_AW.meanAmplitude_blSub_active_ipsi(idcs.recruited.Acatchonly_A{i,j});
            end
            if ~isempty(idcs.recruited.Xcatchonly_X{i,j})
                tempX = d{i,j}.(tng_struct).firingField.Xcatch_AW.meanAmplitude_blSub_active_ipsi(idcs.recruited.Xcatchonly_X{i,j});
            end
            meanAmplitudeActivePhotoactivated(i,j) = nanmean([tempA;tempX]);
            tempA = []; tempX = [];
            if ~isempty(idcs.passed.Acatchonly{i,j})
                tempA = d{i,j}.(tng_struct).firingField.Acatch_AW.meanAmplitude_blSub_active_ipsi(setdiff(idcs.passed.Acatchonly{i,j},idcs.recruited.Acatchonly_A{i,j}));
            end
            if ~isempty(idcs.passed.Xcatchonly{i,j})
                tempX = d{i,j}.(tng_struct).firingField.Xcatch_AW.meanAmplitude_blSub_active_ipsi(setdiff(idcs.passed.Xcatchonly{i,j},idcs.recruited.Xcatchonly_X{i,j}));
            end
            meanAmplitudeActiveOther(i,j)= nanmean([tempA;tempX]);
            
            % split by seq vs ctrl
            if d_info.group(i)==2
                meanAmplitudeActivePhotoactivated_img(i,j) = meanAmplitudeActivePhotoactivated(i,j);
                meanAmplitudeActiveOther_img(i,j) = meanAmplitudeActiveOther(i,j);
            elseif((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
                meanAmplitudeActivePhotoactivated_seq(i,j) = meanAmplitudeActivePhotoactivated(i,j);
                meanAmplitudeActiveOther_seq(i,j) = meanAmplitudeActiveOther(i,j);
            elseif ((d_info.group(i)==7 && (j==4 || j==5)) || (d_info.group(i)==8 && (j==2 || j==3)))
                meanAmplitudeActivePhotoactivated_ctrl(i,j) = meanAmplitudeActivePhotoactivated(i,j);
                meanAmplitudeActiveOther_ctrl(i,j) = meanAmplitudeActiveOther(i,j);
            end
        catch
        end
    end
end


%% Assemble timing matrix

%matrixBinEdges = 0.3:0.5:5.3; % 0.299:0.25:5.299;%
matrixBinEdges = 0:1:5.3;

% scatter plot (n=191)
% this_data_stim = [vertcat(stimTime_trg_seq.Acatchonly_A{:,[2,4]});vertcat(stimTime_trg_seq.Xcatchonly_X{:,[2,4]});vertcat(stimTime_trg_seq.Acatchonly_X{:,[2,4]});vertcat(stimTime_trg_seq.Xcatchonly_A{:,[2,4]})];
% this_data_img = [vertcat(peakTime_trg_seq.Acatchonly_A{:,[2,4]});vertcat(peakTime_trg_seq.Xcatchonly_X{:,[2,4]});vertcat(peakTime_trg_seq.Acatchonly_X{:,[2,4]});vertcat(peakTime_trg_seq.Xcatchonly_A{:,[2,4]})];

% heatmap, responding targets (n=180)
% this_data_idcs_photoactivationTimes = idcs.recruited_trg; % contains Acatchonly_A, Xcatchonly_X, Acatchonly_X, Xcatchonly_A
% this_data_photoactivationTimes = stimTime_trg_seq; % NaN for all idcs with more than one responsive time; contains Acatchonly_A, Xcatchonly_X, Acatchonly_X, Xcatchonly_A
% this_data_firingFieldPeaks = peakTime_trg; % contains Acatchonly_A, Xcatchonly_X, Acatchonly_X, Xcatchonly_A

% heatmap, responders (n=201)
this_data_idcs_photoactivationTimes = idcs.recruited; % contains Acatchonly_A, Xcatchonly_X, Acatchonly_X, Xcatchonly_A
this_data_photoactivationTimes = stimTime_seq; % NaN for all idcs with more than one responsive time; contains Acatchonly_A, Xcatchonly_X, Acatchonly_X, Xcatchonly_A
this_data_firingFieldPeaks = peakTime_seq; % contains Acatchonly_A, Xcatchonly_X, Acatchonly_X, Xcatchonly_A

these_conditions = {'ipsi','contra','all'};
for k=1:length(these_conditions)
    
    % initialise output structs
    timingMatrix.data.(these_conditions{k}) = cell(length(matrixBinEdges)-1);
    timingMatrix.metadata.animal.(these_conditions{k}) = cell(length(matrixBinEdges)-1);
    timingMatrix.metadata.switch.(these_conditions{k}) = cell(length(matrixBinEdges)-1);
    timingMatrix.metadata.switchday.(these_conditions{k}) = cell(length(matrixBinEdges)-1);
    timingMatrix.metadata.trialType.(these_conditions{k}) = cell(length(matrixBinEdges)-1);
    timingMatrix.metadata.running.(these_conditions{k}) = cell(length(matrixBinEdges)-1);
    timingMatrix.count.(these_conditions{k}) = nan(length(matrixBinEdges)-1);
    timingMatrix.norm_withinRow.(these_conditions{k}) = nan(length(matrixBinEdges)-1);
    timingMatrix.norm_catchSequence.(these_conditions{k}) = nan(length(matrixBinEdges)-1);
    timingMatrix.norm_imgSequence.(these_conditions{k}) = nan(length(matrixBinEdges)-1);
    timingMatrix.p_img.(these_conditions{k}) = nan(length(matrixBinEdges)-1);
    timingMatrix.p_stim.(these_conditions{k}) = nan(length(matrixBinEdges)-1);
    timingMatrix.delta_imgSequence.(these_conditions{k}) = nan(length(matrixBinEdges)-1);
    timingMatrix.percdelta_imgSequence.(these_conditions{k}) = nan(length(matrixBinEdges)-1);
    
    % pool data from all sessions
    these_idcs_photoactivationTimes = [];
    these_photoactivationTimes = [];
    these_firingFieldPeaks = [];
    these_firingFieldPeaks_catchSequence = [];
    these_firingFieldPeaks_imgSequence = [];
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if d_info.group(i)==2 && (j==2 || j==4) % first switch days of imaging animals
                these_firingFieldPeaks_imgSequence = [these_firingFieldPeaks_imgSequence; [peakTime.Acatchonly{i,j};peakTime.Xcatchonly{i,j}]];  
            elseif ((d_info.group(i)==7 && j==2) || (d_info.group(i)==8 && j==4)) % first seq stim day of stim animals
                these_firingFieldPeaks_catchSequence = [these_firingFieldPeaks_catchSequence; [peakTime.Acatchonly{i,j};peakTime.Xcatchonly{i,j}]];
                if strcmp(these_conditions{k},'ipsi')
                    these_idcs_photoactivationTimes = [these_idcs_photoactivationTimes; this_data_idcs_photoactivationTimes.Acatchonly_A{i,j};this_data_idcs_photoactivationTimes.Xcatchonly_X{i,j}];
                    these_photoactivationTimes = [these_photoactivationTimes; this_data_photoactivationTimes.Acatchonly_A{i,j};this_data_photoactivationTimes.Xcatchonly_X{i,j}];
                    these_firingFieldPeaks = [these_firingFieldPeaks; this_data_firingFieldPeaks.Acatchonly_A{i,j};this_data_firingFieldPeaks.Xcatchonly_X{i,j}];
                elseif strcmp(these_conditions{k},'contra')
                    these_idcs_photoactivationTimes = [these_idcs_photoactivationTimes; this_data_idcs_photoactivationTimes.Acatchonly_X{i,j};this_data_idcs_photoactivationTimes.Xcatchonly_A{i,j}];
                    these_photoactivationTimes = [these_photoactivationTimes; this_data_photoactivationTimes.Acatchonly_X{i,j};this_data_photoactivationTimes.Xcatchonly_A{i,j}];
                    these_firingFieldPeaks = [these_firingFieldPeaks; this_data_firingFieldPeaks.Acatchonly_X{i,j};this_data_firingFieldPeaks.Xcatchonly_A{i,j}];
                elseif strcmp(these_conditions{k},'all')
                    these_idcs_photoactivationTimes = [these_idcs_photoactivationTimes; this_data_idcs_photoactivationTimes.Acatchonly_A{i,j};this_data_idcs_photoactivationTimes.Xcatchonly_X{i,j};this_data_idcs_photoactivationTimes.Acatchonly_X{i,j};this_data_idcs_photoactivationTimes.Xcatchonly_A{i,j}];
                    these_photoactivationTimes = [these_photoactivationTimes; this_data_photoactivationTimes.Acatchonly_A{i,j};this_data_photoactivationTimes.Xcatchonly_X{i,j};this_data_photoactivationTimes.Acatchonly_X{i,j};this_data_photoactivationTimes.Xcatchonly_A{i,j}];
                    these_firingFieldPeaks = [these_firingFieldPeaks; this_data_firingFieldPeaks.Acatchonly_A{i,j};this_data_firingFieldPeaks.Xcatchonly_X{i,j};this_data_firingFieldPeaks.Acatchonly_X{i,j};this_data_firingFieldPeaks.Xcatchonly_A{i,j}];
%                     [a,b,c] = unique(these_idcs_photoactivationTimes)
%                     if length(these_idcs_photoactivationTimes) ~= length(unique(these_idcs_photoactivationTimes))
%                         warning([num2str(k),',',num2str(i),',',num2str(j)])
%                     end
                end
            end
        end
    end
    these_idcs_photoactivationTimes = these_idcs_photoactivationTimes(~isnan(these_photoactivationTimes));
    these_firingFieldPeaks = these_firingFieldPeaks(~isnan(these_photoactivationTimes));
    these_photoactivationTimes = these_photoactivationTimes(~isnan(these_photoactivationTimes));
    
    % catch sequence without recruited neurons
    these_firingFieldPeaks_catchSequenceWoResp = these_firingFieldPeaks_catchSequence;
    for i=1:length(these_firingFieldPeaks)
        temp = randsample(find(these_firingFieldPeaks_catchSequenceWoResp==these_firingFieldPeaks(i)),1);
        these_firingFieldPeaks_catchSequenceWoResp(temp) = NaN;
    end
    these_firingFieldPeaks_catchSequenceWoResp = rmmissing(these_firingFieldPeaks_catchSequenceWoResp);
    
    % divide into matrix bins
	temp_img = discretize(these_firingFieldPeaks,matrixBinEdges);
    temp_img_catchSequence = discretize(these_firingFieldPeaks_catchSequence,matrixBinEdges);
    temp_img_catchSequenceWoResp = discretize(these_firingFieldPeaks_catchSequenceWoResp,matrixBinEdges);
    temp_img_imgSequence = discretize(these_firingFieldPeaks_imgSequence,matrixBinEdges);
    temp_stim = discretize(these_photoactivationTimes,matrixBinEdges);
    for m=1:length(matrixBinEdges)-1 % firing field peaks
        for n=1:length(matrixBinEdges)-1 % photoactivation times
%             timingMatrix.data.(these_conditions{k}){m,n} = 
            timingMatrix.count.(these_conditions{k})(m,n) = nansum(((temp_img==m) & (temp_stim==n)));
        end
        timingMatrix.norm_withinRow.(these_conditions{k})(m,:) = timingMatrix.count.(these_conditions{k})(m,:) ./ nansum(timingMatrix.count.(these_conditions{k})(m,:));
        timingMatrix.norm_catchSequence.(these_conditions{k})(m,:) = timingMatrix.count.(these_conditions{k})(m,:) ./ nansum(temp_img_catchSequence==m); %./ nansum(temp_img_catchSequence==m);
        timingMatrix.norm_imgSequence.(these_conditions{k})(m,:) = timingMatrix.count.(these_conditions{k})(m,:) ./ nansum(temp_img_imgSequence==m);
    end
    
    % change in probability of firing field formation
    for m=1:length(matrixBinEdges)-1 % firing field peaks
        for n=1:length(matrixBinEdges)-1 % photoactivation times
            temp1 = timingMatrix.count.(these_conditions{k})(m,n) / nansum(timingMatrix.count.(these_conditions{k})(:,n));
            
            % _imgSequence
            temp2 = nansum(temp_img_imgSequence==m) / length(temp_img_imgSequence);
            timingMatrix.p_stim.(these_conditions{k})(m,n) = temp1;
            timingMatrix.p_img.(these_conditions{k})(m,n) = temp2;
            timingMatrix.delta_imgSequence.(these_conditions{k})(m,n) = temp1 - temp2;
            timingMatrix.percdelta_imgSequence.(these_conditions{k})(m,n) = (temp1 - temp2) ./ temp2;
            
            % _catchSequenceWoResp
            temp2 = nansum(temp_img_catchSequenceWoResp==m) / length(temp_img_catchSequenceWoResp);
            timingMatrix.p_stim.(these_conditions{k})(m,n) = temp1;
            timingMatrix.p_catchSequenceWoResp.(these_conditions{k})(m,n) = temp2;
            timingMatrix.delta_catchSequenceWoResp.(these_conditions{k})(m,n) = temp1 - temp2;
            timingMatrix.percdelta_catchSequenceWoResp.(these_conditions{k})(m,n) = (temp1 - temp2) ./ temp2;
        end
    end
    
    % delta change in probability of firing field formation histogram data
    temp = spdiags(timingMatrix.norm_catchSequence.(these_conditions{k}));
    timingMatrix.norm_catchSequence_hist.(these_conditions{k}) = nan(1,size(temp,2));
    for m=1:size(temp,2)
        timingMatrix.norm_catchSequence_hist.(these_conditions{k})(m) = nanmean(nonzeros(temp(:,m)));
    end
    temp = spdiags(timingMatrix.norm_imgSequence.(these_conditions{k}));
    timingMatrix.norm_imgSequence_hist.(these_conditions{k}) = nan(1,size(temp,2));
    for m=1:size(temp,2)
        timingMatrix.norm_imgSequence_hist.(these_conditions{k})(m) = nanmean(nonzeros(temp(:,m)));
    end
    temp = spdiags(timingMatrix.delta_imgSequence.(these_conditions{k}));
    timingMatrix.delta_imgSequence_hist.(these_conditions{k}) = nan(1,size(temp,2));
    for m=1:size(temp,2)
        timingMatrix.delta_imgSequence_hist.(these_conditions{k})(m) = nanmean(nonzeros(temp(:,m)));
    end
    temp = spdiags(timingMatrix.percdelta_imgSequence.(these_conditions{k}));
    timingMatrix.percdelta_imgSequence_hist.(these_conditions{k}) = nan(1,size(temp,2));
    for m=1:size(temp,2)
        timingMatrix.percdelta_imgSequence_hist.(these_conditions{k})(m) = nanmean(nonzeros(temp(:,m)));
    end
    temp = [[0;0;0;0;0],[0;0;0;0;0],spdiags(timingMatrix.count.(these_conditions{k}))];
    timingMatrix.count_hist.(these_conditions{k}) = nan(1,size(temp,2));
    for m=1:size(temp,2)
        timingMatrix.count_hist.(these_conditions{k})(m) = nansum(nonzeros(temp(:,m)));
    end
    % _catchSequenceWoResp
%     temp = spdiags(timingMatrix.norm_imgSequence.(these_conditions{k}));
%     timingMatrix.norm_imgSequence_hist.(these_conditions{k}) = nan(1,size(temp,2));
%     for m=1:size(temp,2)
%         timingMatrix.norm_imgSequence_hist.(these_conditions{k})(m) = nanmean(nonzeros(temp(:,m)));
%     end
    temp = spdiags(timingMatrix.delta_catchSequenceWoResp.(these_conditions{k}));
    timingMatrix.delta_catchSequenceWoResp_hist.(these_conditions{k}) = nan(1,size(temp,2));
    for m=1:size(temp,2)
        timingMatrix.delta_catchSequenceWoResp_hist.(these_conditions{k})(m) = nanmean(nonzeros(temp(:,m)));
    end
    temp = spdiags(timingMatrix.percdelta_catchSequenceWoResp.(these_conditions{k}));
    timingMatrix.percdelta_catchSequenceWoResp_hist.(these_conditions{k}) = nan(1,size(temp,2));
    for m=1:size(temp,2)
        timingMatrix.percdelta_catchSequenceWoResp_hist.(these_conditions{k})(m) = nanmean(nonzeros(temp(:,m)));
    end
end


%% Shuffling photoactivation time labels and recomputing timing matrix

numShuffles = 1000;

these_conditions = {'all'};
for k=1:length(these_conditions)
    
    % initialise output structs
    timingMatrix.count_shuffled.(these_conditions{k}) = nan(length(matrixBinEdges)-1,length(matrixBinEdges)-1,numShuffles);
    timingMatrix.norm_catchSequence_shuffled.(these_conditions{k}) = nan(length(matrixBinEdges)-1,length(matrixBinEdges)-1,numShuffles);
    timingMatrix.delta_catchSequenceWoResp_shuffled.(these_conditions{k}) = nan(length(matrixBinEdges)-1,length(matrixBinEdges)-1,numShuffles);
    timingMatrix.percdelta_catchSequenceWoResp_shuffled.(these_conditions{k}) = nan(length(matrixBinEdges)-1,length(matrixBinEdges)-1,numShuffles);
    timingMatrix.delta_imgSequence_shuffled.(these_conditions{k}) = nan(length(matrixBinEdges)-1,length(matrixBinEdges)-1,numShuffles);
    timingMatrix.percdelta_imgSequence_shuffled.(these_conditions{k}) = nan(length(matrixBinEdges)-1,length(matrixBinEdges)-1,numShuffles);

    % pool data from all sessions
    these_idcs_photoactivationTimes = [];
    these_photoactivationTimes = [];
    these_firingFieldPeaks = [];
    these_firingFieldPeaks_catchSequence = [];
    these_firingFieldPeaks_imgSequence = [];
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if d_info.group(i)==2 && (j==2 || j==4) % first switch days of imaging animals
                these_firingFieldPeaks_imgSequence = [these_firingFieldPeaks_imgSequence; [peakTime.Acatchonly{i,j};peakTime.Xcatchonly{i,j}]];  
            elseif ((d_info.group(i)==7 && j==2) || (d_info.group(i)==8 && j==4)) % first seq stim day of stim animals
                these_firingFieldPeaks_catchSequence = [these_firingFieldPeaks_catchSequence; [peakTime.Acatchonly{i,j};peakTime.Xcatchonly{i,j}]];
                if strcmp(these_conditions{k},'ipsi')
                    these_idcs_photoactivationTimes = [these_idcs_photoactivationTimes; this_data_idcs_photoactivationTimes.Acatchonly_A{i,j};this_data_idcs_photoactivationTimes.Xcatchonly_X{i,j}];
                    these_photoactivationTimes = [these_photoactivationTimes; this_data_photoactivationTimes.Acatchonly_A{i,j};this_data_photoactivationTimes.Xcatchonly_X{i,j}];
                    these_firingFieldPeaks = [these_firingFieldPeaks; this_data_firingFieldPeaks.Acatchonly_A{i,j};this_data_firingFieldPeaks.Xcatchonly_X{i,j}];
                elseif strcmp(these_conditions{k},'contra')
                    these_idcs_photoactivationTimes = [these_idcs_photoactivationTimes; this_data_idcs_photoactivationTimes.Acatchonly_X{i,j};this_data_idcs_photoactivationTimes.Xcatchonly_A{i,j}];
                    these_photoactivationTimes = [these_photoactivationTimes; this_data_photoactivationTimes.Acatchonly_X{i,j};this_data_photoactivationTimes.Xcatchonly_A{i,j}];
                    these_firingFieldPeaks = [these_firingFieldPeaks; this_data_firingFieldPeaks.Acatchonly_X{i,j};this_data_firingFieldPeaks.Xcatchonly_A{i,j}];
                elseif strcmp(these_conditions{k},'all')
                    these_idcs_photoactivationTimes = [these_idcs_photoactivationTimes; this_data_idcs_photoactivationTimes.Acatchonly_A{i,j};this_data_idcs_photoactivationTimes.Xcatchonly_X{i,j};this_data_idcs_photoactivationTimes.Acatchonly_X{i,j};this_data_idcs_photoactivationTimes.Xcatchonly_A{i,j}];
                    these_photoactivationTimes = [these_photoactivationTimes; this_data_photoactivationTimes.Acatchonly_A{i,j};this_data_photoactivationTimes.Xcatchonly_X{i,j};this_data_photoactivationTimes.Acatchonly_X{i,j};this_data_photoactivationTimes.Xcatchonly_A{i,j}];
                    these_firingFieldPeaks = [these_firingFieldPeaks; this_data_firingFieldPeaks.Acatchonly_A{i,j};this_data_firingFieldPeaks.Xcatchonly_X{i,j};this_data_firingFieldPeaks.Acatchonly_X{i,j};this_data_firingFieldPeaks.Xcatchonly_A{i,j}];
%                     [a,b,c] = unique(these_idcs_photoactivationTimes)
%                     if length(these_idcs_photoactivationTimes) ~= length(unique(these_idcs_photoactivationTimes))
%                         warning([num2str(k),',',num2str(i),',',num2str(j)])
%                     end
                end
            end
        end
    end
    these_idcs_photoactivationTimes = these_idcs_photoactivationTimes(~isnan(these_photoactivationTimes));
    these_firingFieldPeaks = these_firingFieldPeaks(~isnan(these_photoactivationTimes));
    these_photoactivationTimes = these_photoactivationTimes(~isnan(these_photoactivationTimes));
    
    % catch sequence without recruited neurons
    these_firingFieldPeaks_catchSequenceWoResp = these_firingFieldPeaks_catchSequence;
    for i=1:length(these_firingFieldPeaks)
        temp = randsample(find(these_firingFieldPeaks_catchSequenceWoResp==these_firingFieldPeaks(i)),1);
        these_firingFieldPeaks_catchSequenceWoResp(temp) = NaN;
    end
    these_firingFieldPeaks_catchSequenceWoResp = rmmissing(these_firingFieldPeaks_catchSequenceWoResp);
    
    for l=1:numShuffles
        
        % divide into matrix bins
        temp_img = discretize(these_firingFieldPeaks,matrixBinEdges);
        temp_img_catchSequence = discretize(these_firingFieldPeaks_catchSequence,matrixBinEdges);
        temp_img_catchSequenceWoResp = discretize(these_firingFieldPeaks_catchSequenceWoResp,matrixBinEdges);
        temp_stim = discretize(these_photoactivationTimes(randperm(length(these_photoactivationTimes))),matrixBinEdges);
        for m=1:length(matrixBinEdges)-1 % firing field peaks
            for n=1:length(matrixBinEdges)-1 % photoactivation times
                timingMatrix.count_shuffled.(these_conditions{k})(m,n,l) = nansum(((temp_img==m) & (temp_stim==n)));
            end
            timingMatrix.norm_catchSequence_shuffled.(these_conditions{k})(m,:,l) = timingMatrix.count_shuffled.(these_conditions{k})(m,:,l) ./ nansum(temp_img_catchSequence==m);
        end
        
        % change in probability of firing field formation
        for m=1:length(matrixBinEdges)-1 % firing field peaks
            for n=1:length(matrixBinEdges)-1 % photoactivation times
                temp1 = timingMatrix.count_shuffled.(these_conditions{k})(m,n,l) / nansum(timingMatrix.count_shuffled.(these_conditions{k})(:,n,l));
                
                % _imgSequence
                temp2 = nansum(temp_img_imgSequence==m) / length(temp_img_imgSequence);
                timingMatrix.delta_imgSequence_shuffled.(these_conditions{k})(m,n,l) = temp1 - temp2;
                timingMatrix.percdelta_imgSequence_shuffled.(these_conditions{k})(m,n,l) = (temp1 - temp2) ./ temp2;
                
                % _catchSequenceWoResp
                temp2 = nansum(temp_img_catchSequenceWoResp==m) / length(temp_img_catchSequenceWoResp);
                timingMatrix.delta_catchSequenceWoResp_shuffled.(these_conditions{k})(m,n,l) = temp1 - temp2;
                timingMatrix.percdelta_catchSequenceWoResp_shuffled.(these_conditions{k})(m,n,l) = (temp1 - temp2) ./ temp2;
            end
        end
    end
    
    % delta change in probability of firing field formation histogram data
    temp = spdiags(timingMatrix.delta_imgSequence_shuffled.(these_conditions{k})(:,:,l));
    timingMatrix.delta_imgSequence_shuffled_hist.(these_conditions{k}) = nan(numShuffles,size(temp,2));
    for l=1:numShuffles
        temp = spdiags(timingMatrix.delta_imgSequence_shuffled.(these_conditions{k})(:,:,l));
        for m=1:size(temp,2)
            timingMatrix.delta_imgSequence_shuffled_hist.(these_conditions{k})(l,m) = nanmean(nonzeros(temp(:,m)));
        end
    end
    
    % delta change in probability of firing field formation histogram data
    temp = spdiags(timingMatrix.delta_catchSequenceWoResp_shuffled.(these_conditions{k})(:,:,l));
    timingMatrix.delta_catchSequenceWoResp_shuffled_hist.(these_conditions{k}) = nan(numShuffles,size(temp,2));
    for l=1:numShuffles
        temp = spdiags(timingMatrix.delta_catchSequenceWoResp_shuffled.(these_conditions{k})(:,:,l));
        for m=1:size(temp,2)
            timingMatrix.delta_catchSequenceWoResp_shuffled_hist.(these_conditions{k})(l,m) = nanmean(nonzeros(temp(:,m)));
        end
    end
end

% COMs for data matrix: 1.3, 2.1, 3.4, 3.5, 3.6 -> 0.3, 0.1, 0.4, -0.5, -1.4 -> avg = -0.23 (230 ms)
this_matrix = timingMatrix.norm_catchSequence.all;
timingMatrix.shifts.norm_catchSequence = nan(size(this_matrix,1),1);
for c=1:size(this_matrix,1)
    temp = this_matrix(:,c);
    timingMatrix.shifts.norm_catchSequence(c) = nansum((1:length(temp))'.*temp)/nansum(temp) - c;
end

% COMs for shuffled matrices: 
this_matrix = timingMatrix.norm_catchSequence_shuffled.all;
timingMatrix.shifts.norm_catchSequence_shuffled = nan(size(this_matrix,1),numShuffles);
for l=1:numShuffles
    for c=1:size(this_matrix,1)
        temp = this_matrix(:,c,l);
        timingMatrix.shifts.norm_catchSequence_shuffled(c,l) = nansum((1:length(temp))'.*temp)/nansum(temp) - c;
    end
end

% above vs below for data matrix
this_matrix = timingMatrix.norm_catchSequence.all;
timingMatrix.triangle.norm_catchSequence = nan;
temp1 = tril(this_matrix,-1); % upper left
temp2 = triu(this_matrix,1); % lower right
timingMatrix.triangle.norm_catchSequence = nansum(temp1(:)) - nansum(temp2(:));

% above vs below for shuffled matrices: 
this_matrix = timingMatrix.norm_catchSequence_shuffled.all;
timingMatrix.triangle.norm_catchSequence_shuffled = nan(numShuffles,1);
for l=1:numShuffles
    temp1 = tril(this_matrix(:,:,l),-1); % upper left
    temp2 = triu(this_matrix(:,:,l),1); % lower right
    timingMatrix.triangle.norm_catchSequence_shuffled(l) = nansum(temp1(:)) - nansum(temp2(:));
end

% delta: COMs for data matrix
this_matrix = timingMatrix.delta_imgSequence.all;
timingMatrix.shifts.delta_imgSequence = nan(size(this_matrix,1),1);
for c=1:size(this_matrix,1)
    temp = this_matrix(:,c);
    timingMatrix.shifts.delta_imgSequence(c) = nansum((1:length(temp))'.*temp)/nansum(temp) - c; % deal with negative numbers
end

% delta: COMs for shuffled matrices
this_matrix = timingMatrix.delta_imgSequence_shuffled.all;
timingMatrix.shifts.delta_imgSequence_shuffled = nan(size(this_matrix,1),numShuffles);
for l=1:numShuffles
    for c=1:size(this_matrix,1)
        temp = this_matrix(:,c,l);
        timingMatrix.shifts.delta_imgSequence_shuffled(c,l) = nansum((1:length(temp))'.*temp)/nansum(temp) - c;
    end
end

% percdelta: COMs for data matrix
this_matrix = timingMatrix.percdelta_imgSequence.all;
timingMatrix.shifts.percdelta_imgSequence = nan(size(this_matrix,1),1);
for c=1:size(this_matrix,1)
    temp = this_matrix(:,c);
    timingMatrix.shifts.percdelta_imgSequence(c) = nansum((1:length(temp))'.*temp)/nansum(temp) - c;
end

% percdelta: COMs for shuffled matrices
this_matrix = timingMatrix.percdelta_imgSequence_shuffled.all;
timingMatrix.shifts.percdelta_imgSequence_shuffled = nan(size(this_matrix,1),numShuffles);
for l=1:numShuffles
    for c=1:size(this_matrix,1)
        temp = this_matrix(:,c,l);
        timingMatrix.shifts.percdelta_imgSequence_shuffled(c,l) = nansum((1:length(temp))'.*temp)/nansum(temp) - c;
    end
end

% --- _catchSequenceWoResp

% delta: COMs for data matrix
this_matrix = timingMatrix.delta_catchSequenceWoResp.all;
timingMatrix.shifts.delta_catchSequenceWoResp = nan(size(this_matrix,1),1);
for c=1:size(this_matrix,1)
    temp = this_matrix(:,c);
    timingMatrix.shifts.delta_catchSequenceWoResp(c) = nansum((1:length(temp))'.*temp)/nansum(temp) - c; % deal with negative numbers
end

% delta: COMs for shuffled matrices
this_matrix = timingMatrix.delta_catchSequenceWoResp_shuffled.all;
timingMatrix.shifts.delta_catchSequenceWoResp_shuffled = nan(size(this_matrix,1),numShuffles);
for l=1:numShuffles
    for c=1:size(this_matrix,1)
        temp = this_matrix(:,c,l);
        timingMatrix.shifts.delta_catchSequenceWoResp_shuffled(c,l) = nansum((1:length(temp))'.*temp)/nansum(temp) - c;
    end
end

% percdelta: COMs for data matrix
this_matrix = timingMatrix.percdelta_catchSequenceWoResp.all;
timingMatrix.shifts.percdelta_catchSequenceWoResp = nan(size(this_matrix,1),1);
for c=1:size(this_matrix,1)
    temp = this_matrix(:,c);
    timingMatrix.shifts.percdelta_catchSequenceWoResp(c) = nansum((1:length(temp))'.*temp)/nansum(temp) - c;
end

% percdelta: COMs for shuffled matrices
this_matrix = timingMatrix.percdelta_catchSequenceWoResp_shuffled.all;
timingMatrix.shifts.percdelta_catchSequenceWoResp_shuffled = nan(size(this_matrix,1),numShuffles);
for l=1:numShuffles
    for c=1:size(this_matrix,1)
        temp = this_matrix(:,c,l);
        timingMatrix.shifts.percdelta_catchSequenceWoResp_shuffled(c,l) = nansum((1:length(temp))'.*temp)/nansum(temp) - c;
    end
end


%% Imprinting analysis

% num
num_cells_seq = nan(d_info.numAnimals,d_info.numDays);
num_cells_ctrl = nan(d_info.numAnimals,d_info.numDays);
num_cells = nan(d_info.numAnimals,d_info.numDays);
num_seqCells_seq = nan(d_info.numAnimals,d_info.numDays);
num_seqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
num_seqCells = nan(d_info.numAnimals,d_info.numDays);
num_respSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
num_respSeqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
num_respSeqCells = nan(d_info.numAnimals,d_info.numDays);
num_respIpsiSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
num_respIpsiSeqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
num_respIpsiSeqCells = nan(d_info.numAnimals,d_info.numDays);

% proportion of respCells
prop_respIpsiSeqCells_respCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_respIpsiSeqCells_respCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_respIpsiSeqCells_respCells = nan(d_info.numAnimals,d_info.numDays);
prop_respContraSeqCells_respCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_respContraSeqCells_respCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_respContraSeqCells_respCells = nan(d_info.numAnimals,d_info.numDays);
prop_respSeqCells_respCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_respSeqCells_respCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_respSeqCells_respCells = nan(d_info.numAnimals,d_info.numDays);
prop_respNonSeqCells_respCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_respNonSeqCells_respCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_respNonSeqCells_respCells = nan(d_info.numAnimals,d_info.numDays);

% proportion of respSeqCells
prop_respIpsiSeqCells_respSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_respIpsiSeqCells_respSeqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_respIpsiSeqCells_respSeqCells = nan(d_info.numAnimals,d_info.numDays);
prop_respContraSeqCells_respSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_respContraSeqCells_respSeqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_respContraSeqCells_respSeqCells = nan(d_info.numAnimals,d_info.numDays);

% proportion of nonRespCells
prop_nonRespSeqCells_nonRespCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespSeqCells_nonRespCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespSeqCells_nonRespCells = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespNonSeqCells_nonRespCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespNonSeqCells_nonRespCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespNonSeqCells_nonRespCells = nan(d_info.numAnimals,d_info.numDays);

% proportion of seqCells
prop_respSeqCells_seqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_respSeqCells_seqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_respSeqCells_seqCells = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespSeqCells_seqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespSeqCells_seqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespSeqCells_seqCells = nan(d_info.numAnimals,d_info.numDays);

% proportion of nonSeqCells
prop_respNonSeqCells_nonSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_respNonSeqCells_nonSeqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_respNonSeqCells_nonSeqCells = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespNonSeqCells_nonSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespNonSeqCells_nonSeqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_nonRespNonSeqCells_nonSeqCells = nan(d_info.numAnimals,d_info.numDays);

% calculate proportions of earlySeqCells or lateSeqCells
prop_respEarlySeqCells_earlySeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_respEarlySeqCells_earlySeqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_respEarlySeqCells_earlySeqCells = nan(d_info.numAnimals,d_info.numDays);
prop_respLateSeqCells_lateSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_respLateSeqCells_lateSeqCells_ctrl = nan(d_info.numAnimals,d_info.numDays);
prop_respLateSeqCells_lateSeqCells = nan(d_info.numAnimals,d_info.numDays);
prop_earlyRespEarlySeqCells_earlySeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_lateRespEarlySeqCells_earlySeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_earlyRespLateSeqCells_lateSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_lateRespLateSeqCells_lateSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_sameCellType_respSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
prop_oppositeCellType_respSeqCells_seq = nan(d_info.numAnimals,d_info.numDays);
% prop_earlyRespEarlySeqCells_earlyRespCells_seq = nan(d_info.numAnimals,d_info.numDays);
% prop_lateRespEarlySeqCells_lateRespCells_seq = nan(d_info.numAnimals,d_info.numDays);
% prop_earlyRespLateSeqCells_earlyRespCells_seq = nan(d_info.numAnimals,d_info.numDays);
% prop_lateRespLateSeqCells_lateRespCells_seq = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try

            % get broad metrics
            temp_allCells = length(idcs.iscells{i,j});
            temp_seqCells = length(idcs.passed.pref{i,j});
            temp_nonSeqCells = length(idcs.passed.none{i,j});
            temp_respCells = length(unique([idcs.resp.A{i,j};idcs.resp.X{i,j}]));
            temp_nonRespCells = length(idcs.resp.O{i,j});
            
            % get finer metrics
            temp_respIpsiSeqCells = length(unique([intersect(idcs.passed.Acatchonly{i,j},idcs.resp.A{i,j});intersect(idcs.passed.Xcatchonly{i,j},idcs.resp.X{i,j})]));
            temp_respContraSeqCells = length(unique([intersect(idcs.passed.Acatchonly{i,j},idcs.resp.X{i,j});intersect(idcs.passed.Xcatchonly{i,j},idcs.resp.A{i,j})]));
            temp_respSeqCells = length(unique([intersect(idcs.passed.Acatchonly{i,j},idcs.resp.A{i,j});intersect(idcs.passed.Xcatchonly{i,j},idcs.resp.X{i,j});intersect(idcs.passed.Acatchonly{i,j},idcs.resp.X{i,j});intersect(idcs.passed.Xcatchonly{i,j},idcs.resp.A{i,j})]));
            temp_respNonSeqCells = length(unique([intersect(idcs.passed.none{i,j},idcs.resp.A{i,j});intersect(idcs.passed.none{i,j},idcs.resp.X{i,j})]));
            temp_nonRespSeqCells = length(unique([intersect(idcs.passed.Acatchonly{i,j},idcs.resp.O{i,j});intersect(idcs.passed.Xcatchonly{i,j},idcs.resp.O{i,j})]));
            temp_nonRespNonSeqCells = length(intersect(idcs.passed.none{i,j},idcs.resp.O{i,j}));
            
            % get early vs late metrics
            temp_earlySeqCells = numByCellType.Acatchonly(i,j,1) + numByCellType.Xcatchonly(i,j,1);
            temp_lateSeqCells = numByCellType.Acatchonly(i,j,2) + numByCellType.Xcatchonly(i,j,2);
            temp_respEarlySeqCells = numByCellType.AcatchonlyResp(i,j,1) + numByCellType.XcatchonlyResp(i,j,1);
            temp_respLateSeqCells = numByCellType.AcatchonlyResp(i,j,2) + numByCellType.XcatchonlyResp(i,j,2);
            temp = unique([stimTime_cleanEarlyIdcs_seq.Acatchonly_A{i,j};stimTime_cleanEarlyIdcs_seq.Xcatchonly_X{i,j};stimTime_cleanEarlyIdcs_seq.Acatchonly_X{i,j};stimTime_cleanEarlyIdcs_seq.Xcatchonly_A{i,j}]);
            temp_earlyRespEarlySeqCells = length(intersect(temp,[idcsByCellType.Acatchonly{i,j,1};idcsByCellType.Xcatchonly{i,j,1}]));
            temp_earlyRespLateSeqCells = length(intersect(temp,[idcsByCellType.Acatchonly{i,j,2};idcsByCellType.Xcatchonly{i,j,2}]));
            temp = unique([stimTime_cleanLateIdcs_seq.Acatchonly_A{i,j};stimTime_cleanLateIdcs_seq.Xcatchonly_X{i,j};stimTime_cleanLateIdcs_seq.Acatchonly_X{i,j};stimTime_cleanLateIdcs_seq.Xcatchonly_A{i,j}]);
            temp_lateRespEarlySeqCells = length(intersect(temp,[idcsByCellType.Acatchonly{i,j,1};idcsByCellType.Xcatchonly{i,j,1}]));
            temp_lateRespLateSeqCells = length(intersect(temp,[idcsByCellType.Acatchonly{i,j,2};idcsByCellType.Xcatchonly{i,j,2}]));
%             temp = unique([stimTime_cleanEarlyIdcs_seq.Acatchonly_A{i,j};stimTime_cleanEarlyIdcs_seq.Xcatchonly_X{i,j};stimTime_cleanEarlyIdcs_seq.Acatchonly_X{i,j};stimTime_cleanEarlyIdcs_seq.Xcatchonly_A{i,j}]);
%             temp_earlyRespCells = length(temp);
%             temp = unique([stimTime_cleanLateIdcs_seq.Acatchonly_A{i,j};stimTime_cleanLateIdcs_seq.Xcatchonly_X{i,j};stimTime_cleanLateIdcs_seq.Acatchonly_X{i,j};stimTime_cleanLateIdcs_seq.Xcatchonly_A{i,j}]);
%             temp_lateRespCells = length(temp);
            
            % get num
            num_cells(i,j) = temp_allCells;
            num_seqCells(i,j) = temp_seqCells;
            num_respSeqCells(i,j) = temp_respSeqCells;
            num_respIpsiSeqCells(i,j) = temp_respIpsiSeqCells;
            
            % calculate proportions of respCells
            prop_respIpsiSeqCells_respCells(i,j) = temp_respIpsiSeqCells / temp_respCells;
            prop_respContraSeqCells_respCells(i,j) = temp_respContraSeqCells / temp_respCells;
            prop_respSeqCells_respCells(i,j) = temp_respSeqCells / temp_respCells;
            prop_respNonSeqCells_respCells(i,j) = temp_respNonSeqCells / temp_respCells;
            
            % calculate proportions of respSeqCells
            prop_respIpsiSeqCells_respSeqCells(i,j) = temp_respIpsiSeqCells / temp_respSeqCells;
            prop_respContraSeqCells_respSeqCells(i,j) = temp_respContraSeqCells / temp_respSeqCells;
            
            % calculate proportions of nonRespCells
            prop_nonRespSeqCells_nonRespCells(i,j) = temp_nonRespSeqCells / temp_nonRespCells;
            prop_nonRespNonSeqCells_nonRespCells(i,j) = temp_nonRespNonSeqCells / temp_nonRespCells;
                        
            % calculate proportions of seqCells
            prop_respSeqCells_seqCells(i,j) = temp_respSeqCells / temp_seqCells;
            prop_nonRespSeqCells_seqCells(i,j) = temp_nonRespSeqCells / temp_seqCells;
            
            % calculate proportions of nonSeqCells
            prop_respNonSeqCells_nonSeqCells(i,j) = temp_respNonSeqCells / temp_nonSeqCells;
            prop_nonRespNonSeqCells_nonSeqCells(i,j) = temp_nonRespNonSeqCells / temp_nonSeqCells;
            
            % calculate proportions of earlySeqCells or lateSeqCells
            prop_respEarlySeqCells_earlySeqCells(i,j) = temp_respEarlySeqCells / temp_earlySeqCells;
            prop_respLateSeqCells_lateSeqCells(i,j) = temp_respLateSeqCells / temp_lateSeqCells;
            
            if ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
                prop_respIpsiSeqCells_respCells_seq(i,j) = prop_respIpsiSeqCells_respCells(i,j);
                prop_respContraSeqCells_respCells_seq(i,j) = prop_respContraSeqCells_respCells(i,j);
                prop_respIpsiSeqCells_respSeqCells_seq(i,j) = prop_respIpsiSeqCells_respSeqCells(i,j);
                prop_respContraSeqCells_respSeqCells_seq(i,j) = prop_respContraSeqCells_respSeqCells(i,j);
                prop_respSeqCells_respCells_seq(i,j) = prop_respSeqCells_respCells(i,j);
                prop_respNonSeqCells_respCells_seq(i,j) = prop_respNonSeqCells_respCells(i,j);
                prop_nonRespSeqCells_nonRespCells_seq(i,j) = prop_nonRespSeqCells_nonRespCells(i,j);
                prop_nonRespNonSeqCells_nonRespCells_seq(i,j) = prop_nonRespNonSeqCells_nonRespCells(i,j);
                prop_respSeqCells_seqCells_seq(i,j) = prop_respSeqCells_seqCells(i,j);
                prop_nonRespSeqCells_seqCells_seq(i,j) = prop_nonRespSeqCells_seqCells(i,j);
                prop_respNonSeqCells_nonSeqCells_seq(i,j) = prop_respNonSeqCells_nonSeqCells(i,j);
                prop_nonRespNonSeqCells_nonSeqCells_seq(i,j) = prop_nonRespNonSeqCells_nonSeqCells(i,j);
                prop_respEarlySeqCells_earlySeqCells_seq(i,j) = prop_respEarlySeqCells_earlySeqCells(i,j);
                prop_respLateSeqCells_lateSeqCells_seq(i,j) = prop_respLateSeqCells_lateSeqCells(i,j);
                num_cells_seq(i,j) = num_cells(i,j);
                num_seqCells_seq(i,j) = num_seqCells(i,j);
                num_respSeqCells_seq(i,j) = num_respSeqCells(i,j);
                num_respIpsiSeqCells_seq(i,j) = num_respIpsiSeqCells(i,j);
                
                % calculate seq-specific proportions
                prop_earlyRespEarlySeqCells_earlySeqCells_seq(i,j) =  temp_earlyRespEarlySeqCells / temp_earlySeqCells;%(temp_earlyRespEarlySeqCells+temp_lateRespEarlySeqCells); %temp_respEarlySeqCells; %temp_earlySeqCells; %(temp_earlyRespEarlySeqCells+temp_lateRespEarlySeqCells);
                prop_lateRespEarlySeqCells_earlySeqCells_seq(i,j) =  temp_lateRespEarlySeqCells / temp_earlySeqCells;%(temp_earlyRespEarlySeqCells+temp_lateRespEarlySeqCells); %temp_respEarlySeqCells; %temp_earlySeqCells; %(temp_earlyRespEarlySeqCells+temp_lateRespEarlySeqCells);
                prop_earlyRespLateSeqCells_lateSeqCells_seq(i,j) =  temp_earlyRespLateSeqCells / temp_lateSeqCells;%(temp_earlyRespLateSeqCells+temp_lateRespLateSeqCells); %temp_respLateSeqCells; %temp_lateSeqCells; %(temp_earlyRespLateSeqCells+temp_lateRespLateSeqCells);
                prop_lateRespLateSeqCells_lateSeqCells_seq(i,j) =  temp_lateRespLateSeqCells / temp_lateSeqCells;%(temp_earlyRespLateSeqCells+temp_lateRespLateSeqCells); %temp_respLateSeqCells; %temp_lateSeqCells; %(temp_earlyRespLateSeqCells+temp_lateRespLateSeqCells);
                prop_sameCellType_respSeqCells_seq(i,j) =  prop_earlyRespEarlySeqCells_earlySeqCells_seq(i,j) + prop_lateRespLateSeqCells_lateSeqCells_seq(i,j);
                prop_oppositeCellType_respSeqCells_seq(i,j) =  prop_lateRespEarlySeqCells_earlySeqCells_seq(i,j) + prop_earlyRespLateSeqCells_lateSeqCells_seq(i,j);                
               
%                 prop_earlyRespEarlySeqCells_earlyRespCells_seq(i,j) =  temp_earlyRespEarlySeqCells / temp_respCells;%(temp_earlyRespEarlySeqCells+temp_lateRespEarlySeqCells); %temp_respEarlySeqCells; %temp_earlySeqCells; %(temp_earlyRespEarlySeqCells+temp_lateRespEarlySeqCells);
%                 prop_lateRespEarlySeqCells_lateRespCells_seq(i,j) =  temp_lateRespEarlySeqCells / temp_respCells;%(temp_earlyRespEarlySeqCells+temp_lateRespEarlySeqCells); %temp_respEarlySeqCells; %temp_earlySeqCells; %(temp_earlyRespEarlySeqCells+temp_lateRespEarlySeqCells);
%                 prop_earlyRespLateSeqCells_earlyRespCells_seq(i,j) =  temp_earlyRespLateSeqCells / temp_respCells;%(temp_earlyRespLateSeqCells+temp_lateRespLateSeqCells); %temp_respLateSeqCells; %temp_lateSeqCells; %(temp_earlyRespLateSeqCells+temp_lateRespLateSeqCells);
%                 prop_lateRespLateSeqCells_lateRespCells_seq(i,j) =  temp_lateRespLateSeqCells / temp_respCells;%(temp_earlyRespLateSeqCells+temp_lateRespLateSeqCells); %temp_respLateSeqCells; %temp_lateSeqCells; %(temp_earlyRespLateSeqCells+temp_lateRespLateSeqCells);                
%                 prop_sameCellType_respSeqCells_seq(i,j) =  prop_earlyRespEarlySeqCells_earlyRespCells_seq(i,j) + prop_lateRespLateSeqCells_lateRespCells_seq(i,j);
%                 prop_oppositeCellType_respSeqCells_seq(i,j) =  prop_lateRespEarlySeqCells_lateRespCells_seq(i,j) + prop_earlyRespLateSeqCells_earlyRespCells_seq(i,j);                

            elseif ((d_info.group(i)==7 && (j==4 || j==5)) || (d_info.group(i)==8 && (j==2 || j==3)))
                prop_respIpsiSeqCells_respCells_ctrl(i,j) = prop_respIpsiSeqCells_respCells(i,j);
                prop_respContraSeqCells_respCells_ctrl(i,j) = prop_respContraSeqCells_respCells(i,j);
                prop_respIpsiSeqCells_respSeqCells_ctrl(i,j) = prop_respIpsiSeqCells_respSeqCells(i,j);
                prop_respContraSeqCells_respSeqCells_ctrl(i,j) = prop_respContraSeqCells_respSeqCells(i,j);      
                prop_respSeqCells_respCells_ctrl(i,j) = prop_respSeqCells_respCells(i,j);
                prop_respNonSeqCells_respCells_ctrl(i,j) = prop_respNonSeqCells_respCells(i,j);
                prop_nonRespSeqCells_nonRespCells_ctrl(i,j) = prop_nonRespSeqCells_nonRespCells(i,j);
                prop_nonRespNonSeqCells_nonRespCells_ctrl(i,j) = prop_nonRespNonSeqCells_nonRespCells(i,j);
                prop_respSeqCells_seqCells_ctrl(i,j) = prop_respSeqCells_seqCells(i,j);
                prop_nonRespSeqCells_seqCells_ctrl(i,j) = prop_nonRespSeqCells_seqCells(i,j);
                prop_respNonSeqCells_nonSeqCells_ctrl(i,j) = prop_respNonSeqCells_nonSeqCells(i,j);
                prop_nonRespNonSeqCells_nonSeqCells_ctrl(i,j) = prop_nonRespNonSeqCells_nonSeqCells(i,j);
                prop_respEarlySeqCells_earlySeqCells_ctrl(i,j) = prop_respEarlySeqCells_earlySeqCells(i,j);
                prop_respLateSeqCells_lateSeqCells_ctrl(i,j) = prop_respLateSeqCells_lateSeqCells(i,j);
                num_cells_ctrl(i,j) = num_cells(i,j);
                num_seqCells_ctrl(i,j) = num_seqCells(i,j);
                num_respSeqCells_ctrl(i,j) = num_respSeqCells(i,j);
                num_respIpsiSeqCells_ctrl(i,j) = num_respIpsiSeqCells(i,j);
            end
        catch
        end
    end
end


%% Reactivation - Correlation of responsive cells

% ipsi
pwcorr.resp_ipsi_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_ipsi_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_ipsi_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.resp_ipsi_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_ipsi_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_ipsi_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.resp_ipsi_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_ipsi_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_ipsi_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.resp_ipsi_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_ipsi_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_ipsi_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                % resp_ipsi_pre_withinZ
                temp1 = d{i,j}.ppa.withinZ_pre(idcs.resp.A{i,j},idcs.resp.A{i,j});
                temp2 = d{i,j}.ppa.withinZ_pre(idcs.resp.X{i,j},idcs.resp.X{i,j});
                pwcorr.resp_ipsi_pre_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);
                
                % resp_ipsi_post_withinZ
                temp1 = d{i,j}.ppa.withinZ_post(idcs.resp.A{i,j},idcs.resp.A{i,j});
                temp2 = d{i,j}.ppa.withinZ_post(idcs.resp.X{i,j},idcs.resp.X{i,j});
                pwcorr.resp_ipsi_post_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);

                % resp_ipsi_pre_acrossZ
                temp1 = d{i,j}.ppa.acrossZ_pre(idcs.resp.A{i,j},idcs.resp.A{i,j});
                temp2 = d{i,j}.ppa.acrossZ_pre(idcs.resp.X{i,j},idcs.resp.X{i,j});
                pwcorr.resp_ipsi_pre_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);
                
                % resp_ipsi_post_acrossZ
                temp1 = d{i,j}.ppa.acrossZ_post(idcs.resp.A{i,j},idcs.resp.A{i,j});
                temp2 = d{i,j}.ppa.acrossZ_post(idcs.resp.X{i,j},idcs.resp.X{i,j});
                pwcorr.resp_ipsi_post_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);
                
                if ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
                    pwcorr_seq.resp_ipsi_pre_withinZ(i,j) = pwcorr.resp_ipsi_pre_withinZ(i,j);
                    pwcorr_seq.resp_ipsi_post_withinZ(i,j) = pwcorr.resp_ipsi_post_withinZ(i,j);
                    pwcorr_seq.resp_ipsi_pre_acrossZ(i,j) = pwcorr.resp_ipsi_pre_acrossZ(i,j);
                    pwcorr_seq.resp_ipsi_post_acrossZ(i,j) = pwcorr.resp_ipsi_post_acrossZ(i,j);
                elseif ((d_info.group(i)==7 && (j==4 || j==5)) || (d_info.group(i)==8 && (j==2 || j==3)))
                    pwcorr_ctrl.resp_ipsi_pre_withinZ(i,j) = pwcorr.resp_ipsi_pre_withinZ(i,j);
                    pwcorr_ctrl.resp_ipsi_post_withinZ(i,j) = pwcorr.resp_ipsi_post_withinZ(i,j);
                    pwcorr_ctrl.resp_ipsi_pre_acrossZ(i,j) = pwcorr.resp_ipsi_pre_acrossZ(i,j);
                    pwcorr_ctrl.resp_ipsi_post_acrossZ(i,j) = pwcorr.resp_ipsi_post_acrossZ(i,j);
                end
            catch
            end
        end
    end
end

% contra
pwcorr.resp_contra_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_contra_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_contra_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.resp_contra_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_contra_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_contra_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.resp_contra_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_contra_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_contra_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.resp_contra_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_contra_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_contra_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                % resp_contra_pre_withinZ
                temp1 = d{i,j}.ppa.withinZ_pre(idcs.resp.A{i,j},idcs.resp.X{i,j});
                temp2 = d{i,j}.ppa.withinZ_pre(idcs.resp.X{i,j},idcs.resp.A{i,j});
                pwcorr.resp_contra_pre_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);
                
                % resp_contra_post_withinZ
                temp1 = d{i,j}.ppa.withinZ_post(idcs.resp.A{i,j},idcs.resp.X{i,j});
                temp2 = d{i,j}.ppa.withinZ_post(idcs.resp.X{i,j},idcs.resp.A{i,j});
                pwcorr.resp_contra_post_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);

                % resp_contra_pre_acrossZ
                temp1 = d{i,j}.ppa.acrossZ_pre(idcs.resp.A{i,j},idcs.resp.X{i,j});
                temp2 = d{i,j}.ppa.acrossZ_pre(idcs.resp.X{i,j},idcs.resp.A{i,j});
                pwcorr.resp_contra_pre_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);
                
                % resp_contra_post_acrossZ
                temp1 = d{i,j}.ppa.acrossZ_post(idcs.resp.A{i,j},idcs.resp.X{i,j});
                temp2 = d{i,j}.ppa.acrossZ_post(idcs.resp.X{i,j},idcs.resp.A{i,j});
                pwcorr.resp_contra_post_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);
                
                if ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
                    pwcorr_seq.resp_contra_pre_withinZ(i,j) = pwcorr.resp_contra_pre_withinZ(i,j);
                    pwcorr_seq.resp_contra_post_withinZ(i,j) = pwcorr.resp_contra_post_withinZ(i,j);
                    pwcorr_seq.resp_contra_pre_acrossZ(i,j) = pwcorr.resp_contra_pre_acrossZ(i,j);
                    pwcorr_seq.resp_contra_post_acrossZ(i,j) = pwcorr.resp_contra_post_acrossZ(i,j);
                elseif ((d_info.group(i)==7 && (j==4 || j==5)) || (d_info.group(i)==8 && (j==2 || j==3)))
                    pwcorr_ctrl.resp_contra_pre_withinZ(i,j) = pwcorr.resp_contra_pre_withinZ(i,j);
                    pwcorr_ctrl.resp_contra_post_withinZ(i,j) = pwcorr.resp_contra_post_withinZ(i,j);
                    pwcorr_ctrl.resp_contra_pre_acrossZ(i,j) = pwcorr.resp_contra_pre_acrossZ(i,j);
                    pwcorr_ctrl.resp_contra_post_acrossZ(i,j) = pwcorr.resp_contra_post_acrossZ(i,j);
                end
            catch
            end
        end
    end
end

% all
pwcorr.resp_all_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_all_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_all_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.resp_all_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_all_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_all_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.resp_all_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_all_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_all_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.resp_all_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.resp_all_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.resp_all_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if d_info.presponsive(i,j)==1
            try
                % resp_all_pre_withinZ
                temp1 = d{i,j}.ppa.withinZ_pre(idcs.resp.A{i,j},idcs.resp.A{i,j});
                temp2 = d{i,j}.ppa.withinZ_pre(idcs.resp.X{i,j},idcs.resp.X{i,j});
                temp3 = d{i,j}.ppa.withinZ_pre(idcs.resp.A{i,j},idcs.resp.X{i,j});
                temp4 = d{i,j}.ppa.withinZ_pre(idcs.resp.X{i,j},idcs.resp.A{i,j});
                pwcorr.resp_all_pre_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1));nonzeros(triu(temp3,1));nonzeros(triu(temp4,1))]);
                
                % resp_all_post_withinZ
                temp1 = d{i,j}.ppa.withinZ_post(idcs.resp.A{i,j},idcs.resp.A{i,j});
                temp2 = d{i,j}.ppa.withinZ_post(idcs.resp.X{i,j},idcs.resp.X{i,j});
                temp3 = d{i,j}.ppa.withinZ_post(idcs.resp.A{i,j},idcs.resp.X{i,j});
                temp4 = d{i,j}.ppa.withinZ_post(idcs.resp.X{i,j},idcs.resp.A{i,j});
                pwcorr.resp_all_post_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1));nonzeros(triu(temp3,1));nonzeros(triu(temp4,1))]);

                % resp_all_pre_acrossZ
                temp1 = d{i,j}.ppa.acrossZ_pre(idcs.resp.A{i,j},idcs.resp.A{i,j});
                temp2 = d{i,j}.ppa.acrossZ_pre(idcs.resp.X{i,j},idcs.resp.X{i,j});
                temp3 = d{i,j}.ppa.acrossZ_pre(idcs.resp.A{i,j},idcs.resp.X{i,j});
                temp4 = d{i,j}.ppa.acrossZ_pre(idcs.resp.X{i,j},idcs.resp.A{i,j});
                pwcorr.resp_all_pre_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1));nonzeros(triu(temp3,1));nonzeros(triu(temp4,1))]);
                
                % resp_all_post_acrossZ
                temp1 = d{i,j}.ppa.acrossZ_post(idcs.resp.A{i,j},idcs.resp.A{i,j});
                temp2 = d{i,j}.ppa.acrossZ_post(idcs.resp.X{i,j},idcs.resp.X{i,j});
                temp3 = d{i,j}.ppa.acrossZ_post(idcs.resp.A{i,j},idcs.resp.X{i,j});
                temp4 = d{i,j}.ppa.acrossZ_post(idcs.resp.X{i,j},idcs.resp.A{i,j});
                pwcorr.resp_all_post_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1));nonzeros(triu(temp3,1));nonzeros(triu(temp4,1))]);
                
                if ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5)))
                    pwcorr_seq.resp_all_pre_withinZ(i,j) = pwcorr.resp_all_pre_withinZ(i,j);
                    pwcorr_seq.resp_all_post_withinZ(i,j) = pwcorr.resp_all_post_withinZ(i,j);
                    pwcorr_seq.resp_all_pre_acrossZ(i,j) = pwcorr.resp_all_pre_acrossZ(i,j);
                    pwcorr_seq.resp_all_post_acrossZ(i,j) = pwcorr.resp_all_post_acrossZ(i,j);
                elseif ((d_info.group(i)==7 && (j==4 || j==5)) || (d_info.group(i)==8 && (j==2 || j==3)))
                    pwcorr_ctrl.resp_all_pre_withinZ(i,j) = pwcorr.resp_all_pre_withinZ(i,j);
                    pwcorr_ctrl.resp_all_post_withinZ(i,j) = pwcorr.resp_all_post_withinZ(i,j);
                    pwcorr_ctrl.resp_all_pre_acrossZ(i,j) = pwcorr.resp_all_pre_acrossZ(i,j);
                    pwcorr_ctrl.resp_all_post_acrossZ(i,j) = pwcorr.resp_all_post_acrossZ(i,j);
                end
            catch
            end
        end
    end
end


%% Reactivation - Correlation of sequence cells

% ipsi
pwcorr.passed_ipsi_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_ipsi_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_ipsi_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_ipsi_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.passed_ipsi_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_ipsi_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_ipsi_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_ipsi_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.passed_ipsi_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_ipsi_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_ipsi_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_ipsi_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.passed_ipsi_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_ipsi_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_ipsi_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_ipsi_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            % passed_ipsi_pre_withinZ
            temp1 = d{i,j}.ppa.withinZ_pre(idcs.passed.Acatchonly{i,j},idcs.passed.Acatchonly{i,j});
            temp2 = d{i,j}.ppa.withinZ_pre(idcs.passed.Xcatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            pwcorr.passed_ipsi_pre_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);

            % passed_ipsi_post_withinZ
            temp1 = d{i,j}.ppa.withinZ_post(idcs.passed.Acatchonly{i,j},idcs.passed.Acatchonly{i,j});
            temp2 = d{i,j}.ppa.withinZ_post(idcs.passed.Xcatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            pwcorr.passed_ipsi_post_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);

            % passed_ipsi_pre_acrossZ
            temp1 = d{i,j}.ppa.acrossZ_pre(idcs.passed.Acatchonly{i,j},idcs.passed.Acatchonly{i,j});
            temp2 = d{i,j}.ppa.acrossZ_pre(idcs.passed.Xcatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            pwcorr.passed_ipsi_pre_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);

            % passed_ipsi_post_acrossZ
            temp1 = d{i,j}.ppa.acrossZ_post(idcs.passed.Acatchonly{i,j},idcs.passed.Acatchonly{i,j});
            temp2 = d{i,j}.ppa.acrossZ_post(idcs.passed.Xcatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            pwcorr.passed_ipsi_post_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);

%                 if any([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]>0.5)
%                     disp('Identical neurons')
%                 end

            if (d_info.group(i)==2 && (j>1))
                pwcorr_img.passed_ipsi_pre_withinZ(i,j) = pwcorr.passed_ipsi_pre_withinZ(i,j);
                pwcorr_img.passed_ipsi_post_withinZ(i,j) = pwcorr.passed_ipsi_post_withinZ(i,j);
                pwcorr_img.passed_ipsi_pre_acrossZ(i,j) = pwcorr.passed_ipsi_pre_acrossZ(i,j);
                pwcorr_img.passed_ipsi_post_acrossZ(i,j) = pwcorr.passed_ipsi_post_acrossZ(i,j);
            elseif ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5))) && d_info.presponsive(i,j)==1
                pwcorr_seq.passed_ipsi_pre_withinZ(i,j) = pwcorr.passed_ipsi_pre_withinZ(i,j);
                pwcorr_seq.passed_ipsi_post_withinZ(i,j) = pwcorr.passed_ipsi_post_withinZ(i,j);
                pwcorr_seq.passed_ipsi_pre_acrossZ(i,j) = pwcorr.passed_ipsi_pre_acrossZ(i,j);
                pwcorr_seq.passed_ipsi_post_acrossZ(i,j) = pwcorr.passed_ipsi_post_acrossZ(i,j);
            elseif ((d_info.group(i)==7 && (j==4 || j==5)) || (d_info.group(i)==8 && (j==2 || j==3))) && d_info.presponsive(i,j)==1
                pwcorr_ctrl.passed_ipsi_pre_withinZ(i,j) = pwcorr.passed_ipsi_pre_withinZ(i,j);
                pwcorr_ctrl.passed_ipsi_post_withinZ(i,j) = pwcorr.passed_ipsi_post_withinZ(i,j);
                pwcorr_ctrl.passed_ipsi_pre_acrossZ(i,j) = pwcorr.passed_ipsi_pre_acrossZ(i,j);
                pwcorr_ctrl.passed_ipsi_post_acrossZ(i,j) = pwcorr.passed_ipsi_post_acrossZ(i,j);
            end
        catch
        end
    end
end

% contra
pwcorr.passed_contra_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_contra_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_contra_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_contra_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.passed_contra_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_contra_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_contra_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_contra_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.passed_contra_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_contra_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_contra_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_contra_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.passed_contra_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_contra_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_contra_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_contra_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            % passed_contra_pre_withinZ
            temp1 = d{i,j}.ppa.withinZ_pre(idcs.passed.Acatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp2 = d{i,j}.ppa.withinZ_pre(idcs.passed.Xcatchonly{i,j},idcs.passed.Acatchonly{i,j});
            pwcorr.passed_contra_pre_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);

            % passed_contra_post_withinZ
            temp1 = d{i,j}.ppa.withinZ_post(idcs.passed.Acatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp2 = d{i,j}.ppa.withinZ_post(idcs.passed.Xcatchonly{i,j},idcs.passed.Acatchonly{i,j});
            pwcorr.passed_contra_post_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);

            % passed_contra_pre_acrossZ
            temp1 = d{i,j}.ppa.acrossZ_pre(idcs.passed.Acatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp2 = d{i,j}.ppa.acrossZ_pre(idcs.passed.Xcatchonly{i,j},idcs.passed.Acatchonly{i,j});
            pwcorr.passed_contra_pre_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);

            % passed_contra_post_acrossZ
            temp1 = d{i,j}.ppa.acrossZ_post(idcs.passed.Acatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp2 = d{i,j}.ppa.acrossZ_post(idcs.passed.Xcatchonly{i,j},idcs.passed.Acatchonly{i,j});
            pwcorr.passed_contra_post_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1))]);

            if (d_info.group(i)==2 && (j>1))
                pwcorr_img.passed_contra_pre_withinZ(i,j) = pwcorr.passed_contra_pre_withinZ(i,j);
                pwcorr_img.passed_contra_post_withinZ(i,j) = pwcorr.passed_contra_post_withinZ(i,j);
                pwcorr_img.passed_contra_pre_acrossZ(i,j) = pwcorr.passed_contra_pre_acrossZ(i,j);
                pwcorr_img.passed_contra_post_acrossZ(i,j) = pwcorr.passed_contra_post_acrossZ(i,j);
            elseif ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5))) && d_info.presponsive(i,j)==1
                pwcorr_seq.passed_contra_pre_withinZ(i,j) = pwcorr.passed_contra_pre_withinZ(i,j);
                pwcorr_seq.passed_contra_post_withinZ(i,j) = pwcorr.passed_contra_post_withinZ(i,j);
                pwcorr_seq.passed_contra_pre_acrossZ(i,j) = pwcorr.passed_contra_pre_acrossZ(i,j);
                pwcorr_seq.passed_contra_post_acrossZ(i,j) = pwcorr.passed_contra_post_acrossZ(i,j);
            elseif ((d_info.group(i)==7 && (j==4 || j==5)) || (d_info.group(i)==8 && (j==2 || j==3))) && d_info.presponsive(i,j)==1
                pwcorr_ctrl.passed_contra_pre_withinZ(i,j) = pwcorr.passed_contra_pre_withinZ(i,j);
                pwcorr_ctrl.passed_contra_post_withinZ(i,j) = pwcorr.passed_contra_post_withinZ(i,j);
                pwcorr_ctrl.passed_contra_pre_acrossZ(i,j) = pwcorr.passed_contra_pre_acrossZ(i,j);
                pwcorr_ctrl.passed_contra_post_acrossZ(i,j) = pwcorr.passed_contra_post_acrossZ(i,j);
            end
        catch
        end
    end
end

% all
pwcorr.passed_all_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_all_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_all_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_all_pre_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.passed_all_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_all_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_all_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_all_post_withinZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.passed_all_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_all_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_all_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_all_pre_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr.passed_all_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_img.passed_all_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_seq.passed_all_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
pwcorr_ctrl.passed_all_post_acrossZ = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        try
            % passed_all_pre_withinZ
            temp1 = d{i,j}.ppa.withinZ_pre(idcs.passed.Acatchonly{i,j},idcs.passed.Acatchonly{i,j});
            temp2 = d{i,j}.ppa.withinZ_pre(idcs.passed.Xcatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp3 = d{i,j}.ppa.withinZ_pre(idcs.passed.Acatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp4 = d{i,j}.ppa.withinZ_pre(idcs.passed.Xcatchonly{i,j},idcs.passed.Acatchonly{i,j});
            pwcorr.passed_all_pre_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1));nonzeros(triu(temp3,1));nonzeros(triu(temp4,1))]);

            % passed_all_post_withinZ
            temp1 = d{i,j}.ppa.withinZ_post(idcs.passed.Acatchonly{i,j},idcs.passed.Acatchonly{i,j});
            temp2 = d{i,j}.ppa.withinZ_post(idcs.passed.Xcatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp3 = d{i,j}.ppa.withinZ_post(idcs.passed.Acatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp4 = d{i,j}.ppa.withinZ_post(idcs.passed.Xcatchonly{i,j},idcs.passed.Acatchonly{i,j});
            pwcorr.passed_all_post_withinZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1));nonzeros(triu(temp3,1));nonzeros(triu(temp4,1))]);

            % passed_all_pre_acrossZ
            temp1 = d{i,j}.ppa.acrossZ_pre(idcs.passed.Acatchonly{i,j},idcs.passed.Acatchonly{i,j});
            temp2 = d{i,j}.ppa.acrossZ_pre(idcs.passed.Xcatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp3 = d{i,j}.ppa.acrossZ_pre(idcs.passed.Acatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp4 = d{i,j}.ppa.acrossZ_pre(idcs.passed.Xcatchonly{i,j},idcs.passed.Acatchonly{i,j});
            pwcorr.passed_all_pre_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1));nonzeros(triu(temp3,1));nonzeros(triu(temp4,1))]);

            % passed_all_post_acrossZ
            temp1 = d{i,j}.ppa.acrossZ_post(idcs.passed.Acatchonly{i,j},idcs.passed.Acatchonly{i,j});
            temp2 = d{i,j}.ppa.acrossZ_post(idcs.passed.Xcatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp3 = d{i,j}.ppa.acrossZ_post(idcs.passed.Acatchonly{i,j},idcs.passed.Xcatchonly{i,j});
            temp4 = d{i,j}.ppa.acrossZ_post(idcs.passed.Xcatchonly{i,j},idcs.passed.Acatchonly{i,j});
            pwcorr.passed_all_post_acrossZ(i,j) = nanmedian([nonzeros(triu(temp1,1));nonzeros(triu(temp2,1));nonzeros(triu(temp3,1));nonzeros(triu(temp4,1))]);

            if (d_info.group(i)==2 && (j>1))
                pwcorr_img.passed_all_pre_withinZ(i,j) = pwcorr.passed_all_pre_withinZ(i,j);
                pwcorr_img.passed_all_post_withinZ(i,j) = pwcorr.passed_all_post_withinZ(i,j);
                pwcorr_img.passed_all_pre_acrossZ(i,j) = pwcorr.passed_all_pre_acrossZ(i,j);
                pwcorr_img.passed_all_post_acrossZ(i,j) = pwcorr.passed_all_post_acrossZ(i,j);
            elseif ((d_info.group(i)==7 && (j==2 || j==3)) || (d_info.group(i)==8 && (j==4 || j==5))) && d_info.presponsive(i,j)==1
                pwcorr_seq.passed_all_pre_withinZ(i,j) = pwcorr.passed_all_pre_withinZ(i,j);
                pwcorr_seq.passed_all_post_withinZ(i,j) = pwcorr.passed_all_post_withinZ(i,j);
                pwcorr_seq.passed_all_pre_acrossZ(i,j) = pwcorr.passed_all_pre_acrossZ(i,j);
                pwcorr_seq.passed_all_post_acrossZ(i,j) = pwcorr.passed_all_post_acrossZ(i,j);
            elseif ((d_info.group(i)==7 && (j==4 || j==5)) || (d_info.group(i)==8 && (j==2 || j==3))) && d_info.presponsive(i,j)==1
                pwcorr_ctrl.passed_all_pre_withinZ(i,j) = pwcorr.passed_all_pre_withinZ(i,j);
                pwcorr_ctrl.passed_all_post_withinZ(i,j) = pwcorr.passed_all_post_withinZ(i,j);
                pwcorr_ctrl.passed_all_pre_acrossZ(i,j) = pwcorr.passed_all_pre_acrossZ(i,j);
                pwcorr_ctrl.passed_all_post_acrossZ(i,j) = pwcorr.passed_all_post_acrossZ(i,j);
            end
        catch
        end
    end
end


%% --- Main figure ---

%% Fig5_ImprintingAnalysis_sci_d1
% each dot is a mouse -> look how it would look session-wise (only different for img)

labels_sci = categorical({'Consistent','Shuffled','No-stim'});
labels_sci = reordercats(labels_sci,{'Consistent','Shuffled','No-stim'});

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,3);
this_data_seq = propSeqCells_seq_unpaired_sd1.pref;
this_data_ctrl = propSeqCells_ctrl_unpaired_sd1.pref;
this_data_img = propSeqCells_img_unpaired_sd1.pref;
this_data(1:length(this_data_seq),1) = this_data_seq;
this_data(1:length(this_data_ctrl),2) = this_data_ctrl;
this_data(1:length(this_data_img),3) = this_data_img;
            
v = violinplot(this_data*100,labels_sci,'bandwidth',4); % to find out what the default bandwidth for each violin is [density, value, bandwidth] = ksdensity(this_data(:,2)*100);
v(1).ViolinColor = p.col.seq; v(2).ViolinColor = p.col.ctrl; v(3).ViolinColor = p.col.imaging; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(3).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none'; v(3).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(3).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.seq; v(2).ScatterPlot.MarkerEdgeColor = p.col.ctrl; v(3).ScatterPlot.MarkerEdgeColor = p.col.imaging; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100; v(3).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(3).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none'; v(3).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15; v(3).MedianPlot.SizeData = 15;

xlim([0.5,3.5])
ylim([0,30])
yticks([0,15,30])
ytickformat('percentage')
ylabel(['Sequence participation probability']);

[temp1,~,temp2] = signrank(this_data_seq,this_data_ctrl);

[temp1n,temp3n,temp2n] = kruskalwallis(this_data); close;
[c,m,h] = multcompare(temp2n,'CType','dunn-sidak');
temp0 = repmat(1:size(this_data,2),size(this_data,1),1);
statsp=kwtest([this_data(:),temp0(:)]);

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_sci_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_sci_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_sci_d1.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_ImprintingAnalysis_sci_d1.txt'],'wt');
fprintf(fid,['\nImprintingAnalysis_sci_d1\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),'\nn(img)=',num2str(length(rmmissing(this_data_img))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nImprintingAnalysis_sci_d1\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),'\nn(img)=',num2str(length(rmmissing(this_data_img))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig5_ImprintingAnalysis_rn_stim_d1

these_labels = categorical({'Photoactivated','Other'});
these_labels = reordercats(these_labels,{'Photoactivated','Other'});
F = paper_figure([0,0.5,mm2inch(1.2*34),mm2inch(34)]); hold on;

this_data_resp_seq = nanmean(prop_respSeqCells_respCells_seq(:,[2,4]),2);
this_data_resp_ctrl = nanmean(prop_respSeqCells_respCells_ctrl(:,[2,4]),2);
this_data_resp = [this_data_resp_seq;this_data_resp_ctrl];
this_data_nonresp_seq = nanmean(prop_nonRespSeqCells_nonRespCells_seq(:,[2,4]),2);
this_data_nonresp_ctrl = nanmean(prop_nonRespSeqCells_nonRespCells_ctrl(:,[2,4]),2);
this_data_nonresp = [this_data_nonresp_seq;this_data_nonresp_ctrl];
v = bar(these_labels,nanmean([this_data_resp,this_data_nonresp],1)*100);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.darkGray; v.EdgeColor = 'none';
temp = propSeqCells_img.pref(:,[2,4]);
% yline(nanmean(temp(:))*100,':','Color',p.col.imaging,'LineWidth',1);
plot(these_labels,[this_data_resp_ctrl,this_data_nonresp_ctrl]*100,'-','Color',p.col.ctrl,'LineWidth',1)
plot(these_labels,[this_data_resp_seq,this_data_nonresp_seq]*100,'-','Color',p.col.seq,'LineWidth',1)
ylim([0,30])
yticks([0,15,30])
ytickformat('percentage')
ylabel({'Sequence participation probability'})

[temp1,~,temp2] = signrank(this_data_resp,this_data_nonresp);
ranksum(this_data_resp,propSeqCells_img_unpaired_sd1.pref);
ranksum(this_data_nonresp,propSeqCells_img_unpaired_sd1.pref);
temp = propSeqCells_img.pref(:,[2,4]);
ranksum(this_data_nonresp,temp(:));

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_rn_stim_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_rn_stim_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_rn_stim_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_ic_stim_d1

these_labels = categorical({'Same','Opposite'});
these_labels = reordercats(these_labels,{'Same','Opposite'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_ipsi_seq = nanmean(prop_respIpsiSeqCells_respCells_seq(:,[2,4]),2);
this_data_ipsi_ctrl = nanmean(prop_respIpsiSeqCells_respCells_ctrl(:,[2,4]),2);
this_data_ipsi = [this_data_ipsi_seq;this_data_ipsi_ctrl];
this_data_contra_seq = nanmean(prop_respContraSeqCells_respCells_seq(:,[2,4]),2);
this_data_contra_ctrl = nanmean(prop_respContraSeqCells_respCells_ctrl(:,[2,4]),2);
this_data_contra = [this_data_contra_seq;this_data_contra_ctrl];
v = bar(these_labels,nanmean([this_data_ipsi,this_data_contra],1)*100);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.darkGray; v.EdgeColor = 'none'; 
plot(these_labels,[this_data_ipsi_ctrl,this_data_contra_ctrl]*100,'-','Color',p.col.ctrl,'LineWidth',1)
plot(these_labels,[this_data_ipsi_seq,this_data_contra_seq]*100,'-','Color',p.col.seq,'LineWidth',1)
ylim([0,30])
yticks([0,15,30])
ytickformat('percentage')
ylabel({'Sequence participation probability'})

[temp1,~,temp2] = signrank(this_data_ipsi,this_data_contra);

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_ic_stim_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_ic_stim_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_ic_stim_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

% xline(1); xline(2); xline(3); xline(4); xline(5);

% this_data_stim = [vertcat(stimTime_trg_seq.Acatchonly_A{:,[2,4]});vertcat(stimTime_trg_seq.Xcatchonly_X{:,[2,4]});vertcat(stimTime_trg_seq.Acatchonly_X{:,[2,4]});vertcat(stimTime_trg_seq.Xcatchonly_A{:,[2,4]})];
% this_data_recr = [vertcat(peakTime_trg_seq.Acatchonly_A{:,[2,4]});vertcat(peakTime_trg_seq.Xcatchonly_X{:,[2,4]});vertcat(peakTime_trg_seq.Acatchonly_X{:,[2,4]});vertcat(peakTime_trg_seq.Xcatchonly_A{:,[2,4]})];
this_data_stim = [vertcat(stimTime_seq.Acatchonly_A{:,[2,4]});vertcat(stimTime_seq.Xcatchonly_X{:,[2,4]});vertcat(stimTime_seq.Acatchonly_X{:,[2,4]});vertcat(stimTime_seq.Xcatchonly_A{:,[2,4]})];
this_data_recr = [vertcat(peakTime_seq.Acatchonly_A{:,[2,4]});vertcat(peakTime_seq.Xcatchonly_X{:,[2,4]});vertcat(peakTime_seq.Acatchonly_X{:,[2,4]});vertcat(peakTime_seq.Xcatchonly_A{:,[2,4]})];
% n = length(rmmissing(this_data_stim))

these_xvals = rmmissing(unique(this_data_stim));
these_yvals = rmmissing(unique(this_data_recr));

these_counts = zeros(length(these_xvals),length(these_yvals));
for i=1:length(this_data_stim)
    %these_counts(find(these_xvals==this_data_stim(i)),find(these_yvals==this_data_recr(i))) = these_counts(find(these_xvals==this_data_stim(i)),find(these_yvals==this_data_recr(i))) + 1;
    these_counts(find(abs(these_xvals-this_data_stim(i))<0.01),find(abs(these_yvals-this_data_recr(i))<0.01)) = these_counts(find(abs(these_xvals-this_data_stim(i))<0.01),find(abs(these_yvals-this_data_recr(i))<0.01)) + 1;
end
this_max = nanmax(these_counts(:))

for x=1:length(these_xvals)
    for y=1:length(these_yvals)
        if these_counts(x,y)==7
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+0/9*(1-p.col.seq))
        elseif these_counts(x,y)==6
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+1/9*(1-p.col.seq))
        elseif these_counts(x,y)==5
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+2/9*(1-p.col.seq))
        elseif these_counts(x,y)==4
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+3/9*(1-p.col.seq))
        elseif these_counts(x,y)==3
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+4/9*(1-p.col.seq))
        elseif these_counts(x,y)==2
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+5/9*(1-p.col.seq))
        elseif these_counts(x,y)==1
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+6/9*(1-p.col.seq))
        end
    end
end
[this_rho,this_pval] = fitLine(this_data_stim,this_data_recr,p.col.black);

xlim([0,5.3])
xticks([0:2.5:5])
xlabel({'Photoactivation time (s)'})
ylim([0,5.3])
yticks([0:2.5:5])
ylabel({'Firing field peak (s)'})
%title(['d1, ipsi/contra, rho=',num2str(this_rho,2),', p=',num2str(this_pval,2)])

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1.pdf']); set(gcf,'Color',[1,1,1])

% F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
% cmap = jet(7);
% cmap(1,:) = p.col.seq+6/9*(1-p.col.seq);
% cmap(2,:) = p.col.seq+5/9*(1-p.col.seq);
% cmap(3,:) = p.col.seq+4/9*(1-p.col.seq);
% cmap(4,:) = p.col.seq+3/9*(1-p.col.seq);
% cmap(5,:) = p.col.seq+2/9*(1-p.col.seq);
% cmap(6,:) = p.col.seq+1/9*(1-p.col.seq);
% cmap(7,:) = p.col.seq+0/9*(1-p.col.seq);
% colormap(cmap);
% h=colorbar;
% 
% savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1_colorbar.fig']);
% saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1_colorbar.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1_colorbar.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_matrix = timingMatrix.delta_catchSequenceWoResp.all; %timingMatrix.delta_imgSequence_shuffled.all(:,:,10)     %timingMatrix.delta_imgSequence.all; %timingMatrix.norm_catchSequence.all; %this_matrix = timingMatrix.delta_imgSequence.all; %timingMatrix.norm_catchSequence.all;
%this_matrix = imgaussfilt(this_matrix,1);

imagesc(this_matrix,[-0.2,0.2]) %imagesc(this_matrix,[0,0.12]) %imagesc(this_matrix) %,[-1,1]) %imagesc(this_matrix,[0,0.12])
plot([0,6],[0,6],'Color',p.col.black,'LineWidth',1,'LineStyle',':');
% plot([0,6]-nanmean(timingMatrix.shifts.norm_catchSequence),[0,6],'Color',p.col.black,'LineWidth',1,'LineStyle','-');
colormap('redblue'); %colormap('hot') %colormap('redblue')%colormap('hot')
daspect([1,1,1])
xlim([0.5,5.5]) %xlim([0.5,10.5])
xticks([0.5,3,5.5]) %xticks([0.5,5.5,10.5])
xticklabels({'0','2.5','5'})
xlabel({'Photoactivation time (s)'})
ylim([0.5,5.5]) %ylim([0.5,10.5])
yticks([0.5,3,5.5]) %yticks([0.5,5.5,10.5])
yticklabels({'0','2.5','5'})
ylabel({'Firing field peak (s)'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1.pdf']); set(gcf,'Color',[1,1,1])

% F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
% imagesc(this_matrix,[0,0.04])
% colormap('hot')
% h=colorbar;
% 
% savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1_colorbar.fig']);
% saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1_colorbar.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1_colorbar.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_shiftHistogram_d1 (original)

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

% histogram(nanmean(timingMatrix.shifts.norm_catchSequence_shuffled,1)*1000,'Normalization','probability','FaceColor',p.col.darkGray,'FaceAlpha',1,'EdgeColor','none');
% xline(nanmean(timingMatrix.shifts.norm_catchSequence,1)*1000,'Color',p.col.seq,'LineWidth',1,'Alpha',1);
% pval = sum(nanmean(timingMatrix.shifts.norm_catchSequence_shuffled,1) < nanmean(timingMatrix.shifts.norm_catchSequence,1)) / length(nanmean(timingMatrix.shifts.norm_catchSequence_shuffled,1))

histogram(nanmean(timingMatrix.shifts.norm_catchSequenceWoResp_shuffled,1)*1000,'Normalization','probability','FaceColor',p.col.darkGray,'FaceAlpha',1,'EdgeColor','none');
xline(nanmean(timingMatrix.shifts.norm_catchSequenceWoResp,1)*1000,'Color',p.col.seq,'LineWidth',1,'Alpha',1);
sum(nanmean(timingMatrix.shifts.norm_catchSequenceWoResp_shuffled,1) < nanmean(timingMatrix.shifts.norm_catchSequenceWoResp,1)) / length(nanmean(timingMatrix.shifts.norm_catchSequenceWoResp_shuffled,1))

xlim([-450,450])
xticks([-400:200:400])
xlabel('Imprinting shift (ms)')
ytickformat('percentage')
ylim([0,0.15])
yticks([0:0.05:0.15])
yticklabels({'0%','5%','10%','15%'})
ylabel('Proportion')

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocationHistXimg_ipsicontra_stim_d1

histBinEdges = 0:0.2:5.3;

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); subplot(1,4,1:3)

h = histogram(these_firingFieldPeaks_catchSequenceWoResp,histBinEdges,'Normalization','probability','FaceColor',p.col.seq,'FaceAlpha',0.3,'EdgeColor','none'); % [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5.3]
%h = histogram(these_firingFieldPeaks_imgSequence,histBinEdges,'Normalization','probability','FaceColor',p.col.imaging,'FaceAlpha',0.3,'EdgeColor','none'); % [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5.3]

xlim([0,5.3])
xticks([0,2.5,5])
xticklabels({'','',''}) 
ylim([0,0.25])
yticks([0,0.25])
yticklabels({'0%','25%'})
ylabel({'Proportion'})
box off
set(gca,'view',[90 -90])

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationHistXimg_ipsicontra_stim_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationHistXimg_ipsicontra_stim_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationHistXimg_ipsicontra_stim_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_shiftDeltaPlot_d1

F = paper_figure([0,0.5,mm2inch(1.2*34),mm2inch(34)]);

% v = bar(timingMatrix.delta_imgSequence_hist.all*100);
% v.FaceColor = 'flat'; v.CData = p.col.seq; v.EdgeColor = 'none';

v = bar(timingMatrix.delta_catchSequenceWoResp_hist.all*100);
v.FaceColor = 'flat'; v.CData = p.col.seq; v.EdgeColor = 'none';

xticks([1:4:9])
xticklabels({'-4','0','4'})
xlabel({'Field formation shift (s)'}) % {'Photoactivation time -','firing field peak (s)'}
ytickformat('percentage')
ylim([-15,15]) % ylim([-8,8])
yticks([-15:7.5:15]) % yticks([-8:4:8])
ylabel({'Stimulation-related change in','field formation probability'})
box off

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_shiftDeltaPlot_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_shiftDeltaPlot_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_shiftDeltaPlot_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_shiftDeltaPlotWithChance_d1

F = paper_figure([0,0.5,mm2inch(1.2*34),mm2inch(34)]); hold on;

yline(0);
% errorbar(1:9,nanmean(timingMatrix.delta_catchSequenceWoResp_shuffled_hist.all,1)*100,nanstd(timingMatrix.delta_catchSequenceWoResp_shuffled_hist.all,1)*100,...
%     'Color',p.col.darkGray,'Marker','.');
v = violinplot(timingMatrix.delta_catchSequenceWoResp_shuffled_hist.all*100);
for i=1:9
    v(i).ViolinColor = p.col.darkGray; v(i).BoxColor = 'k'; v(i).ViolinPlot.EdgeColor = 'none';
    v(i).ScatterPlot.Marker = 'none'; v(i).ScatterPlot.MarkerEdgeColor = p.col.gray; v(i).ScatterPlot.SizeData = 100;
    v(i).BoxPlot.LineWidth = 1; v(i).BoxPlot.FaceColor = p.col.darkGray; v(i).BoxPlot.EdgeColor = p.col.darkGray; v(i).WhiskerPlot.Color = 'none'; 
    v(i).MedianPlot.SizeData = 3; v(i).MedianPlot.MarkerEdgeColor = p.col.darkGray;
    
%     pval = sum(timingMatrix.triangle.norm_catchSequence_shuffled < timingMatrix.triangle.norm_catchSequence) / length(timingMatrix.triangle.norm_catchSequence_shuffled)
end
a=plot(1:9,timingMatrix.delta_catchSequenceWoResp_hist.all*100,'.-','MarkerSize',8,'Color',p.col.seq); % 'MarkerSize',10

% lower = prctile(timingMatrix.delta_catchSequenceWoResp_shuffled_hist.all*100,97.5);
% upper = prctile(timingMatrix.delta_catchSequenceWoResp_shuffled_hist.all*100,2.5);

this_mean = nanmean(timingMatrix.delta_catchSequenceWoResp_shuffled_hist.all*100);
this_realDist2Mean = abs(timingMatrix.delta_catchSequenceWoResp_hist.all*100 - this_mean);
these_shuffledDist2Mean = abs(timingMatrix.delta_catchSequenceWoResp_shuffled_hist.all*100 - this_mean);
pval = sum((these_shuffledDist2Mean >= this_realDist2Mean))/numShuffles;


xlim([0.5,9.5])
xticks([1:4:9])
xticklabels({'-4','0','4'})
xlabel({'Field formation shift (s)'}) % {'Photoactivation time -','firing field peak (s)'}
ytickformat('percentage')
ylim([-20,20]) % ylim([-8,8])
yticks([-20:10:20]) % yticks([-8:4:8])
ylabel({'Photoactivation-related change in','field formation probability'})
box off

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_shiftDeltaPlotWithChance_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_shiftDeltaPlotWithChance_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_shiftDeltaPlotWithChance_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_shiftHistogram_d1 (delta)

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

% histogram(nanmean( timingMatrix.delta_imgSequence_shuffled_hist.all .* (-4:4) ,2)*1000,'Normalization','probability','FaceColor',p.col.darkGray,'FaceAlpha',1,'EdgeColor','none');
% xline(nanmean( timingMatrix.delta_imgSequence_hist.all .* (-4:4) )*1000,'Color',p.col.seq,'LineWidth',1,'Alpha',1);
histogram(nanmean( timingMatrix.delta_catchSequenceWoResp_shuffled_hist.all .* (-4:4) ,2)*1000,'Normalization','probability','FaceColor',p.col.darkGray,'FaceAlpha',1,'EdgeColor','none');
xline(nanmean( timingMatrix.delta_catchSequenceWoResp_hist.all .* (-4:4) )*1000,'Color',p.col.seq,'LineWidth',1,'Alpha',1);

% p value for shift comparison
% pval = sum(nanmean(timingMatrix.shifts.delta_imgSequence_shuffled,1) < nanmean(timingMatrix.shifts.delta_imgSequence,1)) / length(nanmean(timingMatrix.shifts.delta_imgSequence_shuffled,1))
% p value for above vs below diagonal
%pval = sum(timingMatrix.triangle.norm_catchSequence_shuffled < timingMatrix.triangle.norm_catchSequence) / length(timingMatrix.triangle.norm_catchSequence_shuffled)

% xlim([-450,450])
% xticks([-400:200:400])
xlabel('Avg. field formation shift (ms)')
% ytickformat('percentage')
% ylim([0,0.15])
% yticks([0:0.05:0.15])
% yticklabels({'0%','5%','10%','15%'})
ylabel('Proportion')

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.pdf']); set(gcf,'Color',[1,1,1])


%% --- Supplementary figure ---

%% Fig5_ImprintingAnalysis_sci_d2
% each dot is a mouse -> look how it would look session-wise (only different for img)

labels_sci = categorical({'Consistent','Shuffled','No-stim'});
labels_sci = reordercats(labels_sci,{'Consistent','Shuffled','No-stim'});

%F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;
F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,3);
this_data_seq = propSeqCells_seq_unpaired_sd2.pref;
this_data_ctrl = propSeqCells_ctrl_unpaired_sd2.pref;
this_data_img = propSeqCells_img_unpaired_sd2.pref;
this_data(1:length(this_data_seq),1) = this_data_seq;
this_data(1:length(this_data_ctrl),2) = this_data_ctrl;
this_data(1:length(this_data_img),3) = this_data_img;
            
v = violinplot(this_data*100,labels_sci,'bandwidth',4); % to find out what the default bandwidth for each violin is [density, value, bandwidth] = ksdensity(this_data(:,2)*100);
v(1).ViolinColor = p.col.seq; v(2).ViolinColor = p.col.ctrl; v(3).ViolinColor = p.col.imaging; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(3).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none'; v(3).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(3).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.seq; v(2).ScatterPlot.MarkerEdgeColor = p.col.ctrl; v(3).ScatterPlot.MarkerEdgeColor = p.col.imaging; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100; v(3).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(3).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none'; v(3).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15; v(3).MedianPlot.SizeData = 15;

xlim([0.5,3.5])
ylim([0,20])
yticks([0,10,20])
ytickformat('percentage')
ylabel(['Sequence participation probability']);

[temp1,~,temp2] = signrank(this_data_seq,this_data_ctrl);

% [temp1n,temp3n,temp2n] = kruskalwallis(this_data); close;
% [c,m,h] = multcompare(temp2n,'CType','dunn-sidak');
% temp0 = repmat(1:size(this_data,2),size(this_data,1),1);
% statsp=kwtest([this_data(:),temp0(:)]);

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_sci_d2.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_sci_d2.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_sci_d2.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_ImprintingAnalysis_sci_d2.txt'],'wt');
fprintf(fid,['\nImprintingAnalysis_sci_d2\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),'\nn(img)=',num2str(length(rmmissing(this_data_img))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']); fclose(fid);
fprintf(['\nImprintingAnalysis_sci_d2\nn(seq)=',num2str(length(rmmissing(this_data_seq))),'\nn(ctrl)=',num2str(length(rmmissing(this_data_ctrl))),'\nn(img)=',num2str(length(rmmissing(this_data_img))),...
    '\nkruskalwallis p=',num2str(kruskalwallis(this_data),4),'\n']);
close; close; close; close;


%% Fig5_ImprintingAnalysis_rn_stim_d2

these_labels = categorical({'Photoactivated','Other'});
these_labels = reordercats(these_labels,{'Photoactivated','Other'});
F = paper_figure([0,0.5,mm2inch(1.2*34),mm2inch(34)]); hold on;

this_data_resp_seq = nanmean(prop_respSeqCells_respCells_seq(:,[3,5]),2);
this_data_resp_ctrl = nanmean(prop_respSeqCells_respCells_ctrl(:,[3,5]),2);
this_data_resp = [this_data_resp_seq;this_data_resp_ctrl];
this_data_nonresp_seq = nanmean(prop_nonRespSeqCells_nonRespCells_seq(:,[3,5]),2);
this_data_nonresp_ctrl = nanmean(prop_nonRespSeqCells_nonRespCells_ctrl(:,[3,5]),2);
this_data_nonresp = [this_data_nonresp_seq;this_data_nonresp_ctrl];
v = bar(these_labels,nanmean([this_data_resp,this_data_nonresp],1)*100);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.darkGray; v.EdgeColor = 'none';
temp = propSeqCells_img.pref(:,[3,5]);
yline(nanmean(temp(:))*100,':','Color',p.col.imaging,'LineWidth',1);
plot(these_labels,[this_data_resp_ctrl,this_data_nonresp_ctrl]*100,'-','Color',p.col.ctrl,'LineWidth',1)
plot(these_labels,[this_data_resp_seq,this_data_nonresp_seq]*100,'-','Color',p.col.seq,'LineWidth',1)
ylim([0,20])
yticks([0,10,20])
ytickformat('percentage')
ylabel({'Sequence participation probability'})

[temp1,~,temp2] = signrank(this_data_resp,this_data_nonresp);
ranksum(this_data_resp,propSeqCells_img_unpaired_sd1.pref);
ranksum(this_data_nonresp,propSeqCells_img_unpaired_sd1.pref);
temp = propSeqCells_img.pref(:,[2,4]);
ranksum(this_data_nonresp,temp(:));

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_rn_stim_d2.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_rn_stim_d2.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_rn_stim_d2.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_ic_stim_d2

these_labels = categorical({'Same','Opposite'});
these_labels = reordercats(these_labels,{'Same','Opposite'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_ipsi_seq = nanmean(prop_respIpsiSeqCells_respCells_seq(:,[3,5]),2);
this_data_ipsi_ctrl = nanmean(prop_respIpsiSeqCells_respCells_ctrl(:,[3,5]),2);
this_data_ipsi = [this_data_ipsi_seq;this_data_ipsi_ctrl];
this_data_contra_seq = nanmean(prop_respContraSeqCells_respCells_seq(:,[3,5]),2);
this_data_contra_ctrl = nanmean(prop_respContraSeqCells_respCells_ctrl(:,[3,5]),2);
this_data_contra = [this_data_contra_seq;this_data_contra_ctrl];
v = bar(these_labels,nanmean([this_data_ipsi,this_data_contra],1)*100);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.darkGray; v.EdgeColor = 'none'; 
plot(these_labels,[this_data_ipsi_ctrl,this_data_contra_ctrl]*100,'-','Color',p.col.ctrl,'LineWidth',1)
plot(these_labels,[this_data_ipsi_seq,this_data_contra_seq]*100,'-','Color',p.col.seq,'LineWidth',1)
ylim([0,20])
yticks([0,10,20])
ytickformat('percentage')
ylabel({'Sequence participation probability'})

[temp1,~,temp2] = signrank(this_data_ipsi,this_data_contra);

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_ic_stim_d2.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_ic_stim_d2.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_ic_stim_d2.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d2

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

% this_data_stim = [vertcat(stimTime_trg_seq.Acatchonly_A{:,[2,4]});vertcat(stimTime_trg_seq.Xcatchonly_X{:,[2,4]});vertcat(stimTime_trg_seq.Acatchonly_X{:,[2,4]});vertcat(stimTime_trg_seq.Xcatchonly_A{:,[2,4]})];
% this_data_recr = [vertcat(peakTime_trg_seq.Acatchonly_A{:,[2,4]});vertcat(peakTime_trg_seq.Xcatchonly_X{:,[2,4]});vertcat(peakTime_trg_seq.Acatchonly_X{:,[2,4]});vertcat(peakTime_trg_seq.Xcatchonly_A{:,[2,4]})];
this_data_stim = [vertcat(stimTime_seq.Acatchonly_A{:,[3,5]});vertcat(stimTime_seq.Xcatchonly_X{:,[3,5]});vertcat(stimTime_seq.Acatchonly_X{:,[3,5]});vertcat(stimTime_seq.Xcatchonly_A{:,[3,5]})];
this_data_recr = [vertcat(peakTime_seq.Acatchonly_A{:,[3,5]});vertcat(peakTime_seq.Xcatchonly_X{:,[3,5]});vertcat(peakTime_seq.Acatchonly_X{:,[3,5]});vertcat(peakTime_seq.Xcatchonly_A{:,[3,5]})];

these_xvals = rmmissing(unique(this_data_stim));
these_yvals = rmmissing(unique(this_data_recr));

these_counts = zeros(length(these_xvals),length(these_yvals));
for i=1:length(this_data_stim)
    %these_counts(find(these_xvals==this_data_stim(i)),find(these_yvals==this_data_recr(i))) = these_counts(find(these_xvals==this_data_stim(i)),find(these_yvals==this_data_recr(i))) + 1;
    these_counts(find(abs(these_xvals-this_data_stim(i))<0.01),find(abs(these_yvals-this_data_recr(i))<0.01)) = these_counts(find(abs(these_xvals-this_data_stim(i))<0.01),find(abs(these_yvals-this_data_recr(i))<0.01)) + 1;
end
this_max = nanmax(these_counts(:))

for x=1:length(these_xvals)
    for y=1:length(these_yvals)
        if these_counts(x,y)==7
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+0/9*(1-p.col.seq))
        elseif these_counts(x,y)==6
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+1/9*(1-p.col.seq))
        elseif these_counts(x,y)==5
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+2/9*(1-p.col.seq))
        elseif these_counts(x,y)==4
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+3/9*(1-p.col.seq))
        elseif these_counts(x,y)==3
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+4/9*(1-p.col.seq))
        elseif these_counts(x,y)==2
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+5/9*(1-p.col.seq))
        elseif these_counts(x,y)==1
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+6/9*(1-p.col.seq))
        end
    end
end
[this_rho,this_pval] = fitLine(this_data_stim,this_data_recr,p.col.black);

xlim([0,5.3])
xticks([0:2.5:5])
xlabel({'Photoactivation time (s)'})
ylim([0,5.3])
yticks([0:2.5:5])
ylabel({'Firing field peak (s)'})
%title(['d1, ipsi/contra, rho=',num2str(this_rho,2),', p=',num2str(this_pval,2)])

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d2.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d2.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d2.pdf']); set(gcf,'Color',[1,1,1])

% F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
% cmap = jet(7);
% cmap(1,:) = p.col.seq+6/9*(1-p.col.seq);
% cmap(2,:) = p.col.seq+5/9*(1-p.col.seq);
% cmap(3,:) = p.col.seq+4/9*(1-p.col.seq);
% cmap(4,:) = p.col.seq+3/9*(1-p.col.seq);
% cmap(5,:) = p.col.seq+2/9*(1-p.col.seq);
% cmap(6,:) = p.col.seq+1/9*(1-p.col.seq);
% cmap(7,:) = p.col.seq+0/9*(1-p.col.seq);
% colormap(cmap);
% h=colorbar;
% 
% savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1_colorbar.fig']);
% saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1_colorbar.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1_colorbar.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocationHistX_ipsicontra_stim_d2

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]);

subplot(1,2,1)

this_data_recr = [vertcat(peakTime_seq.Acatchonly_A{:,[3,5]});vertcat(peakTime_seq.Xcatchonly_X{:,[3,5]});vertcat(peakTime_seq.Acatchonly_X{:,[3,5]});vertcat(peakTime_seq.Xcatchonly_A{:,[3,5]})];
h = histogram(this_data_recr,[0,1,2,3,4,5.3],'FaceColor',p.col.seq,'FaceAlpha',0.3,'EdgeColor','none'); % [0,1,2,3,4,5.3], [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5.3]

xlim([0,5.3])
xticks([0,2.5,5])
xticklabels({'','',''})
ylim([0,50])
yticks([0,50])
ylabel({'Count'})
box off
set(gca,'view',[90 -90])

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationHistX_ipsicontra_stim_d2.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationHistX_ipsicontra_stim_d2.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationHistX_ipsicontra_stim_d2.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocationHistY_ipsicontra_stim_d2

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]);

subplot(2,1,2)

this_data_stim = [vertcat(stimTime_seq.Acatchonly_A{:,[3,5]});vertcat(stimTime_seq.Xcatchonly_X{:,[3,5]});vertcat(stimTime_seq.Acatchonly_X{:,[3,5]});vertcat(stimTime_seq.Xcatchonly_A{:,[3,5]})];
h = histogram(this_data_stim,[0,1,2,3,4,5.3],'FaceColor',p.col.seq,'FaceAlpha',0.3,'EdgeColor','none'); % [0,1,2,3,4,5.3], [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5.3]

xlim([0,5.3])
xticks([0,2.5,5])
xticklabels({'','',''})
ylim([0,20])
yticks([0,20])
ylabel({'Count'})
box off

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationHistY_ipsicontra_stim_d2.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationHistY_ipsicontra_stim_d2.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationHistY_ipsicontra_stim_d2.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_matrix = repelem(timingMatrix.norm_catchSequence.all,1,1);
%this_matrix = imgaussfilt(this_matrix,1);

imagesc(this_matrix,[0,0.12]) %imagesc(this_matrix,[0,0.05]) % imagesc(this_matrix,[0,7])
plot([0,6],[0,6],'Color',p.col.black,'LineWidth',1,'LineStyle',':');
% plot([0,6]-nanmean(timingMatrix.shifts.norm_catchSequence),[0,6],'Color',p.col.black,'LineWidth',1,'LineStyle','-');
colormap('hot')
daspect([1,1,1])
xlim([0.5,5.5]) %xlim([0.5,10.5])
xticks([0.5,3,5.5]) %xticks([0.5,5.5,10.5])
xticklabels({'0','2.5','5'})
xlabel({'Photoactivation time (s)'})
ylim([0.5,5.5]) %ylim([0.5,10.5])
yticks([0.5,3,5.5]) %yticks([0.5,5.5,10.5])
yticklabels({'0','2.5','5'})
ylabel({'Firing field peak (s)'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1.pdf']); set(gcf,'Color',[1,1,1])

% F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
% imagesc(this_matrix,[0,0.04])
% colormap('hot')
% h=colorbar;
% 
% savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1_colorbar.fig']);
% saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1_colorbar.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationMatrix_ipsicontra_stim_d1_colorbar.pdf']); set(gcf,'Color',[1,1,1])


%% --- Reserve ---

%% RESEVE - Fig5_ImprintingAnalysis_recruitedLocation_ipsi_stim_d1

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

% this_data_stim = [vertcat(stimTime_trg_seq.Acatchonly_A{:,[2,4]});vertcat(stimTime_trg_seq.Xcatchonly_X{:,[2,4]});vertcat(stimTime_trg_seq.Acatchonly_X{:,[2,4]});vertcat(stimTime_trg_seq.Xcatchonly_A{:,[2,4]})];
% this_data_recr = [vertcat(peakTime_trg_seq.Acatchonly_A{:,[2,4]});vertcat(peakTime_trg_seq.Xcatchonly_X{:,[2,4]});vertcat(peakTime_trg_seq.Acatchonly_X{:,[2,4]});vertcat(peakTime_trg_seq.Xcatchonly_A{:,[2,4]})];
this_data_stim = [vertcat(stimTime_seq.Acatchonly_A{:,[2,4]});vertcat(stimTime_seq.Xcatchonly_X{:,[2,4]})];
this_data_recr = [vertcat(peakTime_seq.Acatchonly_A{:,[2,4]});vertcat(peakTime_seq.Xcatchonly_X{:,[2,4]})];

these_xvals = rmmissing(unique(this_data_stim));
these_yvals = rmmissing(unique(this_data_recr));

these_counts = zeros(length(these_xvals),length(these_yvals));
for i=1:length(this_data_stim)
    %these_counts(find(these_xvals==this_data_stim(i)),find(these_yvals==this_data_recr(i))) = these_counts(find(these_xvals==this_data_stim(i)),find(these_yvals==this_data_recr(i))) + 1;
    these_counts(find(abs(these_xvals-this_data_stim(i))<0.01),find(abs(these_yvals-this_data_recr(i))<0.01)) = these_counts(find(abs(these_xvals-this_data_stim(i))<0.01),find(abs(these_yvals-this_data_recr(i))<0.01)) + 1;
end
this_max = nanmax(these_counts(:))

for x=1:length(these_xvals)
    for y=1:length(these_yvals)
        if these_counts(x,y)==7
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+0/9*(1-p.col.seq))
        elseif these_counts(x,y)==6
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+1/9*(1-p.col.seq))
        elseif these_counts(x,y)==5
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+2/9*(1-p.col.seq))
        elseif these_counts(x,y)==4
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+3/9*(1-p.col.seq))
        elseif these_counts(x,y)==3
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+4/9*(1-p.col.seq))
        elseif these_counts(x,y)==2
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+5/9*(1-p.col.seq))
        elseif these_counts(x,y)==1
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+6/9*(1-p.col.seq))
        end
    end
end
[this_rho,this_pval] = fitLine(this_data_stim,this_data_recr,p.col.black)

xlim([0,5.3])
xticks([0:2.5:5])
xlabel({'Photoactivation time (s)'})
ylim([0,5.3])
yticks([0:2.5:5])
ylabel({'Firing field peak (s)'})
%title(['d1, ipsi/contra, rho=',num2str(this_rho,2),', p=',num2str(this_pval,2)])

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsi_stim_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsi_stim_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsi_stim_d1.pdf']); set(gcf,'Color',[1,1,1])

% F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
% cmap = jet(7);
% cmap(1,:) = p.col.seq+6/9*(1-p.col.seq);
% cmap(2,:) = p.col.seq+5/9*(1-p.col.seq);
% cmap(3,:) = p.col.seq+4/9*(1-p.col.seq);
% cmap(4,:) = p.col.seq+3/9*(1-p.col.seq);
% cmap(5,:) = p.col.seq+2/9*(1-p.col.seq);
% cmap(6,:) = p.col.seq+1/9*(1-p.col.seq);
% cmap(7,:) = p.col.seq+0/9*(1-p.col.seq);
% colormap(cmap);
% h=colorbar;
% 
% savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsi_stim_d1_colorbar.fig']);
% saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsi_stim_d1_colorbar.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsi_stim_d1_colorbar.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocationComp_ipsicontra_stim_nostim_d1

histBinEdges = 0:0.2:5.3;

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_ctrl = [vertcat(peakTime_ctrl.Acatchonly{:,[2,4]});vertcat(peakTime_ctrl.Xcatchonly{:,[2,4]})]; % default: Acatchonly and Xcatchonly
this_data_seq = [vertcat(peakTime_seq.Acatchonly{:,[2,4]});vertcat(peakTime_seq.Xcatchonly{:,[2,4]})];

this_data_curve_ctrl = histcounts(this_data_ctrl,histBinEdges,'Normalization','pdf'); % 'pdf'
plot(histBinEdges(1:end-1),this_data_curve_ctrl,'-','Color',p.col.ctrl)
this_data_curve_seq = histcounts(this_data_seq,histBinEdges,'Normalization','pdf'); % 'pdf'
plot(histBinEdges(1:end-1),this_data_curve_seq,'-','Color',p.col.seq)

this_data_ctrl = [vertcat(peakTime_ctrl.Acatchonly_A{:,[2,4]});vertcat(peakTime_ctrl.Xcatchonly_X{:,[2,4]})]; % default: Acatchonly and Xcatchonly
this_data_seq = [vertcat(peakTime_seq.Acatchonly_A{:,[2,4]});vertcat(peakTime_seq.Xcatchonly_X{:,[2,4]})];

this_data_curve_ctrl = histcounts(this_data_ctrl,histBinEdges,'Normalization','pdf'); % 'pdf'
plot(histBinEdges(1:end-1),this_data_curve_ctrl,':','Color',p.col.ctrl)
this_data_curve_seq = histcounts(this_data_seq,histBinEdges,'Normalization','pdf'); % 'pdf'
plot(histBinEdges(1:end-1),this_data_curve_seq,':','Color',p.col.seq)

this_data_img = [vertcat(peakTime_img.Acatchonly{:,[2,4]});vertcat(peakTime_img.Xcatchonly{:,[2,4]})]; % default: Acatchonly and Xcatchonly
this_data_curve_img = histcounts(this_data_img,histBinEdges,'Normalization','pdf'); % 'pdf'
plot(histBinEdges(1:end-1),this_data_curve_img,'-','Color',p.col.imaging)

xlim([0,5.3])
xlabel({'Firing field peak (s)'})
%ylim([0,5.3])
ylabel({'pdf'})
%title(['d1, ipsi/contra, rho=',num2str(this_rho,2),', p=',num2str(this_pval,2)])

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationComp_ipsicontra_stim_nostim_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationComp_ipsicontra_stim_nostim_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationComp_ipsicontra_stim_nostim_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocationComp_ipsicontra_stim_nostim_d2

histBinEdges = 0:0.2:5.3;

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_ctrl = [vertcat(peakTime_ctrl.Acatchonly{:,[3,5]});vertcat(peakTime_ctrl.Xcatchonly{:,[3,5]})]; % default: Acatchonly and Xcatchonly
this_data_seq = [vertcat(peakTime_seq.Acatchonly{:,[3,5]});vertcat(peakTime_seq.Xcatchonly{:,[3,5]})];

this_data_curve_ctrl = histcounts(this_data_ctrl,histBinEdges,'Normalization','pdf'); % 'pdf'
plot(histBinEdges(1:end-1),this_data_curve_ctrl,'-','Color',p.col.ctrl)
this_data_curve_seq = histcounts(this_data_seq,histBinEdges,'Normalization','pdf'); % 'pdf'
plot(histBinEdges(1:end-1),this_data_curve_seq,'-','Color',p.col.seq)

this_data_ctrl = [vertcat(peakTime_ctrl.Acatchonly_A{:,[3,5]});vertcat(peakTime_ctrl.Xcatchonly_X{:,[3,5]})]; % default: Acatchonly and Xcatchonly
this_data_seq = [vertcat(peakTime_seq.Acatchonly_A{:,[3,5]});vertcat(peakTime_seq.Xcatchonly_X{:,[3,5]})];

this_data_curve_ctrl = histcounts(this_data_ctrl,histBinEdges,'Normalization','pdf'); % 'pdf'
plot(histBinEdges(1:end-1),this_data_curve_ctrl,':','Color',p.col.ctrl)
this_data_curve_seq = histcounts(this_data_seq,histBinEdges,'Normalization','pdf'); % 'pdf'
plot(histBinEdges(1:end-1),this_data_curve_seq,':','Color',p.col.seq)

this_data_img = [vertcat(peakTime_img.Acatchonly{:,[3,5]});vertcat(peakTime_img.Xcatchonly{:,[3,5]})]; % default: Acatchonly and Xcatchonly
this_data_curve_img = histcounts(this_data_img,histBinEdges,'Normalization','pdf'); % 'pdf'
plot(histBinEdges(1:end-1),this_data_curve_img,'-','Color',p.col.imaging)

xlim([0,5.3])
xlabel({'Firing field peak (s)'})
%ylim([0,5.3])
ylabel({'pdf'})
%title(['d1, ipsi/contra, rho=',num2str(this_rho,2),', p=',num2str(this_pval,2)])

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationComp_ipsicontra_stim_nostim_d2.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationComp_ipsicontra_stim_nostim_d2.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationComp_ipsicontra_stim_nostim_d2.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocationHistX_ipsicontra_stim_d1

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]);

subplot(1,2,1)

this_data_recr = [vertcat(peakTime_seq.Acatchonly_A{:,[2,4]});vertcat(peakTime_seq.Xcatchonly_X{:,[2,4]});vertcat(peakTime_seq.Acatchonly_X{:,[2,4]});vertcat(peakTime_seq.Xcatchonly_A{:,[2,4]})];
h = histogram(this_data_recr,[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5.3],'FaceColor',p.col.seq,'FaceAlpha',0.3,'EdgeColor','none'); % [0,1,2,3,4,5.3], [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5.3]

xlim([0,5.3])
xticks([0,2.5,5])
xticklabels({'','',''})
ylim([0,250])
yticks([0,250])
ylabel({'Count'})
box off
set(gca,'view',[90 -90])

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationHistX_ipsicontra_stim_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationHistX_ipsicontra_stim_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationHistX_ipsicontra_stim_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_recruitedLocationHistY_ipsicontra_stim_d1

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]);

subplot(2,1,2)

this_data_stim = [vertcat(stimTime_seq.Acatchonly_A{:,[2,4]});vertcat(stimTime_seq.Xcatchonly_X{:,[2,4]});vertcat(stimTime_seq.Acatchonly_X{:,[2,4]});vertcat(stimTime_seq.Xcatchonly_A{:,[2,4]})];
h = histogram(this_data_stim,[0,1,2,3,4,5.3],'FaceColor',p.col.seq,'FaceAlpha',0.3,'EdgeColor','none'); % [0,1,2,3,4,5.3], [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5.3]

xlim([0,5.3])
xticks([0,2.5,5])
xticklabels({'','',''})
ylim([0,100])
yticks([0,50,100])
ylabel({'Count'})
box off

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocationHistY_ipsicontra_stim_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocationHistY_ipsicontra_stim_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocationHistY_ipsicontra_stim_d1.pdf']); set(gcf,'Color',[1,1,1])


%% --- In progress ---





%% Fig5_ImprintingAnalysis_reactivation_resp_ipsi_withinZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.resp_ipsi_pre_withinZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.resp_ipsi_post_withinZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.resp_ipsi_pre_withinZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.resp_ipsi_post_withinZ(:,[2,4]),2);

v = bar(these_labels,nanmean([this_data_pre_seq,this_data_post_seq,this_data_pre_ctrl,this_data_post_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.FaceColor = 'flat'; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none';
plot(these_labels(1:2),[this_data_pre_seq,this_data_post_seq],'-','Color',p.col.black,'LineWidth',1)
plot(these_labels(3:4),[this_data_pre_ctrl,this_data_post_ctrl],'-','Color',p.col.black,'LineWidth',1)

ylabel({'Corr between resp cells','ipsi, withinZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_resp_ipsi_withinZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_resp_ipsi_withinZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_resp_ipsi_withinZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_reactivation_resp_contra_withinZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.resp_contra_pre_withinZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.resp_contra_post_withinZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.resp_contra_pre_withinZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.resp_contra_post_withinZ(:,[2,4]),2);

v = bar(these_labels,nanmean([this_data_pre_seq,this_data_post_seq,this_data_pre_ctrl,this_data_post_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.FaceColor = 'flat'; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none';
plot(these_labels(1:2),[this_data_pre_seq,this_data_post_seq],'-','Color',p.col.black,'LineWidth',1)
plot(these_labels(3:4),[this_data_pre_ctrl,this_data_post_ctrl],'-','Color',p.col.black,'LineWidth',1)

ylabel({'Corr between resp cells','contra, withinZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_resp_contra_withinZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_resp_contra_withinZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_resp_contra_withinZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_reactivation_resp_ipsi_acrossZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.resp_ipsi_pre_acrossZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.resp_ipsi_post_acrossZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.resp_ipsi_pre_acrossZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.resp_ipsi_post_acrossZ(:,[2,4]),2);

v = bar(these_labels,nanmean([this_data_pre_seq,this_data_post_seq,this_data_pre_ctrl,this_data_post_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.FaceColor = 'flat'; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none';
plot(these_labels(1:2),[this_data_pre_seq,this_data_post_seq],'-','Color',p.col.black,'LineWidth',1)
plot(these_labels(3:4),[this_data_pre_ctrl,this_data_post_ctrl],'-','Color',p.col.black,'LineWidth',1)

ylabel({'Corr between resp cells','ipsi, acrossZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_resp_ipsi_acrossZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_resp_ipsi_acrossZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_resp_ipsi_acrossZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_reactivation_resp_contra_acrossZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.resp_contra_pre_acrossZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.resp_contra_post_acrossZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.resp_contra_pre_acrossZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.resp_contra_post_acrossZ(:,[2,4]),2);

v = bar(these_labels,nanmean([this_data_pre_seq,this_data_post_seq,this_data_pre_ctrl,this_data_post_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.FaceColor = 'flat'; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none';
plot(these_labels(1:2),[this_data_pre_seq,this_data_post_seq],'-','Color',p.col.black,'LineWidth',1)
plot(these_labels(3:4),[this_data_pre_ctrl,this_data_post_ctrl],'-','Color',p.col.black,'LineWidth',1)

ylabel({'Corr between resp cells','contra, acrossZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_resp_contra_acrossZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_resp_contra_acrossZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_resp_contra_acrossZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_reactivation_passed_ipsi_withinZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.passed_ipsi_pre_withinZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.passed_ipsi_post_withinZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.passed_ipsi_pre_withinZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.passed_ipsi_post_withinZ(:,[2,4]),2);

v = bar(these_labels,nanmean([this_data_pre_seq,this_data_post_seq,this_data_pre_ctrl,this_data_post_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.FaceColor = 'flat'; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none';
plot(these_labels(1:2),[this_data_pre_seq,this_data_post_seq],'-','Color',p.col.black,'LineWidth',1)
plot(these_labels(3:4),[this_data_pre_ctrl,this_data_post_ctrl],'-','Color',p.col.black,'LineWidth',1)

ylabel({'Corr between sequence cells','ipsi, withinZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_passed_ipsi_withinZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_passed_ipsi_withinZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_passed_ipsi_withinZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_reactivation_passed_contra_withinZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.passed_contra_pre_withinZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.passed_contra_post_withinZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.passed_contra_pre_withinZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.passed_contra_post_withinZ(:,[2,4]),2);

v = bar(these_labels,nanmean([this_data_pre_seq,this_data_post_seq,this_data_pre_ctrl,this_data_post_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.FaceColor = 'flat'; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none';
plot(these_labels(1:2),[this_data_pre_seq,this_data_post_seq],'-','Color',p.col.black,'LineWidth',1)
plot(these_labels(3:4),[this_data_pre_ctrl,this_data_post_ctrl],'-','Color',p.col.black,'LineWidth',1)

ylabel({'Corr between sequence cells','contra, withinZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_passed_contra_withinZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_passed_contra_withinZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_passed_contra_withinZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_reactivation_passed_ipsi_acrossZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.passed_ipsi_pre_acrossZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.passed_ipsi_post_acrossZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.passed_ipsi_pre_acrossZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.passed_ipsi_post_acrossZ(:,[2,4]),2);

v = bar(these_labels,nanmean([this_data_pre_seq,this_data_post_seq,this_data_pre_ctrl,this_data_post_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.FaceColor = 'flat'; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none';
plot(these_labels(1:2),[this_data_pre_seq,this_data_post_seq],'-','Color',p.col.black,'LineWidth',1)
plot(these_labels(3:4),[this_data_pre_ctrl,this_data_post_ctrl],'-','Color',p.col.black,'LineWidth',1)

ylabel({'Corr between sequence cells','ipsi, acrossZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_passed_ipsi_acrossZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_passed_ipsi_acrossZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_passed_ipsi_acrossZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_reactivation_passed_contra_acrossZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.passed_contra_pre_acrossZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.passed_contra_post_acrossZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.passed_contra_pre_acrossZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.passed_contra_post_acrossZ(:,[2,4]),2);

v = bar(these_labels,nanmean([this_data_pre_seq,this_data_post_seq,this_data_pre_ctrl,this_data_post_ctrl],1));
v.FaceColor = 'flat'; v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.FaceColor = 'flat'; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl; v.EdgeColor = 'none';
plot(these_labels(1:2),[this_data_pre_seq,this_data_post_seq],'-','Color',p.col.black,'LineWidth',1)
plot(these_labels(3:4),[this_data_pre_ctrl,this_data_post_ctrl],'-','Color',p.col.black,'LineWidth',1)

ylabel({'Corr between sequence cells','contra, acrossZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_passed_contra_acrossZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_passed_contra_acrossZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_passed_contra_acrossZ.pdf']); set(gcf,'Color',[1,1,1])




%% Fig5_ImprintingAnalysis_reactivation_v_passed_ipsi_withinZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post','no-pre','no-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post','no-pre','no-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.passed_ipsi_pre_withinZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.passed_ipsi_post_withinZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.passed_ipsi_pre_withinZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.passed_ipsi_post_withinZ(:,[2,4]),2);
this_data_pre_img = nanmean(pwcorr_img.passed_ipsi_pre_withinZ(:,[2,4]),2);
this_data_post_img = nanmean(pwcorr_img.passed_ipsi_post_withinZ(:,[2,4]),2);

this_data = nan(d_info.numAnimals,6);
this_data(1:length(this_data_pre_seq),1) = this_data_pre_seq;
this_data(1:length(this_data_post_seq),2) = this_data_post_seq;
this_data(1:length(this_data_pre_ctrl),3) = this_data_pre_ctrl;
this_data(1:length(this_data_post_ctrl),4) = this_data_post_ctrl;
this_data(1:length(this_data_pre_img),5) = this_data_pre_img;
this_data(1:length(this_data_post_img),6) = this_data_post_img;

v = violinplot(this_data,these_labels,'bandwidth',0.01); %,'bandwidth',4); % to find out what the default bandwidth for each violin is [density, value, bandwidth] = ksdensity(this_data(:,2)*100);
v(1).ViolinColor = p.col.seq; v(2).ViolinColor = p.col.seq; v(3).ViolinColor = p.col.ctrl; v(4).ViolinColor = p.col.ctrl; v(5).ViolinColor = p.col.imaging; v(6).ViolinColor = p.col.imaging;
v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(3).BoxColor = 'k'; v(4).BoxColor = 'k'; v(5).BoxColor = 'k'; v(6).BoxColor = 'k';
v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none'; v(3).ViolinPlot.EdgeColor = 'none'; v(4).ViolinPlot.EdgeColor = 'none'; v(5).ViolinPlot.EdgeColor = 'none'; v(6).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(3).ScatterPlot.Marker = '.'; v(4).ScatterPlot.Marker = '.'; v(5).ScatterPlot.Marker = '.'; v(6).ScatterPlot.Marker = '.';
v(1).ScatterPlot.MarkerEdgeColor =  p.col.seq; v(2).ScatterPlot.MarkerEdgeColor = p.col.seq; v(3).ScatterPlot.MarkerEdgeColor = p.col.ctrl; v(4).ScatterPlot.MarkerEdgeColor =  p.col.ctrl; v(5).ScatterPlot.MarkerEdgeColor = p.col.imaging; v(6).ScatterPlot.MarkerEdgeColor = p.col.imaging; 
v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100; v(3).ScatterPlot.SizeData = 100; v(4).ScatterPlot.SizeData = 100; v(5).ScatterPlot.SizeData = 100; v(6).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(3).BoxPlot.LineWidth = 1; v(4).BoxPlot.LineWidth = 1; v(5).BoxPlot.LineWidth = 1; v(6).BoxPlot.LineWidth = 1;
v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none'; v(3).WhiskerPlot.Color = 'none'; v(4).WhiskerPlot.Color = 'none'; v(5).WhiskerPlot.Color = 'none'; v(6).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15; v(3).MedianPlot.SizeData = 15; v(4).MedianPlot.SizeData = 15; v(5).MedianPlot.SizeData = 15; v(6).MedianPlot.SizeData = 15;

ylabel({'Corr between sequence cells','ipsi, withinZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_v_passed_ipsi_withinZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_v_passed_ipsi_withinZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_v_passed_ipsi_withinZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_reactivation_v_passed_contra_withinZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post','no-pre','no-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post','no-pre','no-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.passed_contra_pre_withinZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.passed_contra_post_withinZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.passed_contra_pre_withinZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.passed_contra_post_withinZ(:,[2,4]),2);
this_data_pre_img = nanmean(pwcorr_img.passed_contra_pre_withinZ(:,[2,4]),2);
this_data_post_img = nanmean(pwcorr_img.passed_contra_post_withinZ(:,[2,4]),2);

this_data = nan(d_info.numAnimals,6);
this_data(1:length(this_data_pre_seq),1) = this_data_pre_seq;
this_data(1:length(this_data_post_seq),2) = this_data_post_seq;
this_data(1:length(this_data_pre_ctrl),3) = this_data_pre_ctrl;
this_data(1:length(this_data_post_ctrl),4) = this_data_post_ctrl;
this_data(1:length(this_data_pre_img),5) = this_data_pre_img;
this_data(1:length(this_data_post_img),6) = this_data_post_img;

v = violinplot(this_data,these_labels,'bandwidth',0.01); %,'bandwidth',4); % to find out what the default bandwidth for each violin is [density, value, bandwidth] = ksdensity(this_data(:,2)*100);
v(1).ViolinColor = p.col.seq; v(2).ViolinColor = p.col.seq; v(3).ViolinColor = p.col.ctrl; v(4).ViolinColor = p.col.ctrl; v(5).ViolinColor = p.col.imaging; v(6).ViolinColor = p.col.imaging;
v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(3).BoxColor = 'k'; v(4).BoxColor = 'k'; v(5).BoxColor = 'k'; v(6).BoxColor = 'k';
v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none'; v(3).ViolinPlot.EdgeColor = 'none'; v(4).ViolinPlot.EdgeColor = 'none'; v(5).ViolinPlot.EdgeColor = 'none'; v(6).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(3).ScatterPlot.Marker = '.'; v(4).ScatterPlot.Marker = '.'; v(5).ScatterPlot.Marker = '.'; v(6).ScatterPlot.Marker = '.';
v(1).ScatterPlot.MarkerEdgeColor =  p.col.seq; v(2).ScatterPlot.MarkerEdgeColor = p.col.seq; v(3).ScatterPlot.MarkerEdgeColor = p.col.ctrl; v(4).ScatterPlot.MarkerEdgeColor =  p.col.ctrl; v(5).ScatterPlot.MarkerEdgeColor = p.col.imaging; v(6).ScatterPlot.MarkerEdgeColor = p.col.imaging; 
v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100; v(3).ScatterPlot.SizeData = 100; v(4).ScatterPlot.SizeData = 100; v(5).ScatterPlot.SizeData = 100; v(6).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(3).BoxPlot.LineWidth = 1; v(4).BoxPlot.LineWidth = 1; v(5).BoxPlot.LineWidth = 1; v(6).BoxPlot.LineWidth = 1;
v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none'; v(3).WhiskerPlot.Color = 'none'; v(4).WhiskerPlot.Color = 'none'; v(5).WhiskerPlot.Color = 'none'; v(6).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15; v(3).MedianPlot.SizeData = 15; v(4).MedianPlot.SizeData = 15; v(5).MedianPlot.SizeData = 15; v(6).MedianPlot.SizeData = 15;

ylabel({'Corr between sequence cells','contra, withinZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_v_passed_contra_withinZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_v_passed_contra_withinZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_v_passed_contra_withinZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_reactivation_v_passed_ipsi_acrossZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post','no-pre','no-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post','no-pre','no-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.passed_ipsi_pre_acrossZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.passed_ipsi_post_acrossZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.passed_ipsi_pre_acrossZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.passed_ipsi_post_acrossZ(:,[2,4]),2);
this_data_pre_img = nanmean(pwcorr_img.passed_ipsi_pre_acrossZ(:,[2,4]),2);
this_data_post_img = nanmean(pwcorr_img.passed_ipsi_post_acrossZ(:,[2,4]),2);

this_data = nan(d_info.numAnimals,6);
this_data(1:length(this_data_pre_seq),1) = this_data_pre_seq;
this_data(1:length(this_data_post_seq),2) = this_data_post_seq;
this_data(1:length(this_data_pre_ctrl),3) = this_data_pre_ctrl;
this_data(1:length(this_data_post_ctrl),4) = this_data_post_ctrl;
this_data(1:length(this_data_pre_img),5) = this_data_pre_img;
this_data(1:length(this_data_post_img),6) = this_data_post_img;

v = violinplot(this_data,these_labels,'bandwidth',0.01); %,'bandwidth',4); % to find out what the default bandwidth for each violin is [density, value, bandwidth] = ksdensity(this_data(:,2)*100);
v(1).ViolinColor = p.col.seq; v(2).ViolinColor = p.col.seq; v(3).ViolinColor = p.col.ctrl; v(4).ViolinColor = p.col.ctrl; v(5).ViolinColor = p.col.imaging; v(6).ViolinColor = p.col.imaging;
v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(3).BoxColor = 'k'; v(4).BoxColor = 'k'; v(5).BoxColor = 'k'; v(6).BoxColor = 'k';
v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none'; v(3).ViolinPlot.EdgeColor = 'none'; v(4).ViolinPlot.EdgeColor = 'none'; v(5).ViolinPlot.EdgeColor = 'none'; v(6).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(3).ScatterPlot.Marker = '.'; v(4).ScatterPlot.Marker = '.'; v(5).ScatterPlot.Marker = '.'; v(6).ScatterPlot.Marker = '.';
v(1).ScatterPlot.MarkerEdgeColor =  p.col.seq; v(2).ScatterPlot.MarkerEdgeColor = p.col.seq; v(3).ScatterPlot.MarkerEdgeColor = p.col.ctrl; v(4).ScatterPlot.MarkerEdgeColor =  p.col.ctrl; v(5).ScatterPlot.MarkerEdgeColor = p.col.imaging; v(6).ScatterPlot.MarkerEdgeColor = p.col.imaging; 
v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100; v(3).ScatterPlot.SizeData = 100; v(4).ScatterPlot.SizeData = 100; v(5).ScatterPlot.SizeData = 100; v(6).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(3).BoxPlot.LineWidth = 1; v(4).BoxPlot.LineWidth = 1; v(5).BoxPlot.LineWidth = 1; v(6).BoxPlot.LineWidth = 1;
v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none'; v(3).WhiskerPlot.Color = 'none'; v(4).WhiskerPlot.Color = 'none'; v(5).WhiskerPlot.Color = 'none'; v(6).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15; v(3).MedianPlot.SizeData = 15; v(4).MedianPlot.SizeData = 15; v(5).MedianPlot.SizeData = 15; v(6).MedianPlot.SizeData = 15;

ylabel({'Corr between sequence cells','ipsi, acrossZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_v_passed_ipsi_acrossZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_v_passed_ipsi_acrossZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_v_passed_ipsi_acrossZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_reactivation_v_passed_contra_acrossZ

these_labels = categorical({'seq-pre','seq-post','ctrl-pre','ctrl-post','no-pre','no-post'});
these_labels = reordercats(these_labels,{'seq-pre','seq-post','ctrl-pre','ctrl-post','no-pre','no-post'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_pre_seq = nanmean(pwcorr_seq.passed_contra_pre_acrossZ(:,[2,4]),2);
this_data_post_seq = nanmean(pwcorr_seq.passed_contra_post_acrossZ(:,[2,4]),2);
this_data_pre_ctrl = nanmean(pwcorr_ctrl.passed_contra_pre_acrossZ(:,[2,4]),2);
this_data_post_ctrl = nanmean(pwcorr_ctrl.passed_contra_post_acrossZ(:,[2,4]),2);
this_data_pre_img = nanmean(pwcorr_img.passed_contra_pre_acrossZ(:,[2,4]),2);
this_data_post_img = nanmean(pwcorr_img.passed_contra_post_acrossZ(:,[2,4]),2);

this_data = nan(d_info.numAnimals,6);
this_data(1:length(this_data_pre_seq),1) = this_data_pre_seq;
this_data(1:length(this_data_post_seq),2) = this_data_post_seq;
this_data(1:length(this_data_pre_ctrl),3) = this_data_pre_ctrl;
this_data(1:length(this_data_post_ctrl),4) = this_data_post_ctrl;
this_data(1:length(this_data_pre_img),5) = this_data_pre_img;
this_data(1:length(this_data_post_img),6) = this_data_post_img;

v = violinplot(this_data,these_labels,'bandwidth',0.01); %,'bandwidth',4); % to find out what the default bandwidth for each violin is [density, value, bandwidth] = ksdensity(this_data(:,2)*100);
v(1).ViolinColor = p.col.seq; v(2).ViolinColor = p.col.seq; v(3).ViolinColor = p.col.ctrl; v(4).ViolinColor = p.col.ctrl; v(5).ViolinColor = p.col.imaging; v(6).ViolinColor = p.col.imaging;
v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(3).BoxColor = 'k'; v(4).BoxColor = 'k'; v(5).BoxColor = 'k'; v(6).BoxColor = 'k';
v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none'; v(3).ViolinPlot.EdgeColor = 'none'; v(4).ViolinPlot.EdgeColor = 'none'; v(5).ViolinPlot.EdgeColor = 'none'; v(6).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(3).ScatterPlot.Marker = '.'; v(4).ScatterPlot.Marker = '.'; v(5).ScatterPlot.Marker = '.'; v(6).ScatterPlot.Marker = '.';
v(1).ScatterPlot.MarkerEdgeColor =  p.col.seq; v(2).ScatterPlot.MarkerEdgeColor = p.col.seq; v(3).ScatterPlot.MarkerEdgeColor = p.col.ctrl; v(4).ScatterPlot.MarkerEdgeColor =  p.col.ctrl; v(5).ScatterPlot.MarkerEdgeColor = p.col.imaging; v(6).ScatterPlot.MarkerEdgeColor = p.col.imaging; 
v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100; v(3).ScatterPlot.SizeData = 100; v(4).ScatterPlot.SizeData = 100; v(5).ScatterPlot.SizeData = 100; v(6).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(3).BoxPlot.LineWidth = 1; v(4).BoxPlot.LineWidth = 1; v(5).BoxPlot.LineWidth = 1; v(6).BoxPlot.LineWidth = 1;
v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none'; v(3).WhiskerPlot.Color = 'none'; v(4).WhiskerPlot.Color = 'none'; v(5).WhiskerPlot.Color = 'none'; v(6).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15; v(3).MedianPlot.SizeData = 15; v(4).MedianPlot.SizeData = 15; v(5).MedianPlot.SizeData = 15; v(6).MedianPlot.SizeData = 15;

ylabel({'Corr between sequence cells','contra, acrossZ'})

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_reactivation_v_passed_contra_acrossZ.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_reactivation_v_passed_contra_acrossZ.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_reactivation_v_passed_contra_acrossZ.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_shiftHistogram_d1 (delta)
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

histogram(nanmean(timingMatrix.shifts.delta_imgSequence_shuffled,1)*1000,'Normalization','probability','FaceColor',p.col.darkGray,'FaceAlpha',1,'EdgeColor','none');
xline(nanmean(timingMatrix.shifts.delta_imgSequence,1)*1000,'Color',p.col.seq,'LineWidth',1,'Alpha',1);
%pval = (sum(nanmean(timingMatrix.shifts.norm_catchSequence_shuffled,1) < nanmean(timingMatrix.shifts.norm_catchSequence,1)) + sum(nanmean(timingMatrix.shifts.norm_catchSequence_shuffled,1) > -nanmean(timingMatrix.shifts.norm_catchSequence,1))) / length(nanmean(timingMatrix.shifts.norm_catchSequence_shuffled,1))
% p value for shift comparison
pval = sum(nanmean(timingMatrix.shifts.delta_imgSequence_shuffled,1) < nanmean(timingMatrix.shifts.delta_imgSequence,1)) / length(nanmean(timingMatrix.shifts.delta_imgSequence_shuffled,1))
% p value for above vs below diagonal
%pval = sum(timingMatrix.triangle.norm_catchSequence_shuffled < timingMatrix.triangle.norm_catchSequence) / length(timingMatrix.triangle.norm_catchSequence_shuffled)

% xlim([-450,450])
% xticks([-400:200:400])
xlabel('Imprinting shift (ms)')
ytickformat('percentage')
ylim([0,0.15])
yticks([0:0.05:0.15])
yticklabels({'0%','5%','10%','15%'})
ylabel('Proportion')

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_ImprintingAnalysis_shiftHistogram_d1 (percdelta)
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

histogram(nanmean(timingMatrix.shifts.percdelta_imgSequence_shuffled,1)*1000,'Normalization','probability','FaceColor',p.col.darkGray,'FaceAlpha',1,'EdgeColor','none');
xline(nanmean(timingMatrix.shifts.percdelta_imgSequence,1)*1000,'Color',p.col.seq,'LineWidth',1,'Alpha',1);
%pval = (sum(nanmean(timingMatrix.shifts.norm_catchSequence_shuffled,1) < nanmean(timingMatrix.shifts.norm_catchSequence,1)) + sum(nanmean(timingMatrix.shifts.norm_catchSequence_shuffled,1) > -nanmean(timingMatrix.shifts.norm_catchSequence,1))) / length(nanmean(timingMatrix.shifts.norm_catchSequence_shuffled,1))
% p value for shift comparison
pval = sum(nanmean(timingMatrix.shifts.percdelta_imgSequence_shuffled,1) < nanmean(timingMatrix.shifts.percdelta_imgSequence,1)) / length(nanmean(timingMatrix.shifts.percdelta_imgSequence_shuffled,1))
% p value for above vs below diagonal
%pval = sum(timingMatrix.triangle.norm_catchSequence_shuffled < timingMatrix.triangle.norm_catchSequence) / length(timingMatrix.triangle.norm_catchSequence_shuffled)

% xlim([-450,450])
% xticks([-400:200:400])
xlabel('Imprinting shift (ms)')
ytickformat('percentage')
ylim([0,0.15])
yticks([0:0.05:0.15])
yticklabels({'0%','5%','10%','15%'})
ylabel('Proportion')

savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.fig']);
saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_shiftHistogram_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Sequence visualisation





%% Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
this_data_stim = [vertcat(stimTime_seq.Acatchonly_A{:,[2,4]});vertcat(stimTime_seq.Xcatchonly_X{:,[2,4]});vertcat(stimTime_seq.Acatchonly_X{:,[2,4]});vertcat(stimTime_seq.Xcatchonly_A{:,[2,4]})];
this_data_recr = [vertcat(peakTime_seq.Acatchonly_A{:,[2,4]});vertcat(peakTime_seq.Xcatchonly_X{:,[2,4]});vertcat(peakTime_seq.Acatchonly_X{:,[2,4]});vertcat(peakTime_seq.Xcatchonly_A{:,[2,4]})];

these_xvals = rmmissing(unique(this_data_stim));
these_yvals = rmmissing(unique(this_data_recr));

these_counts = zeros(length(these_xvals),length(these_yvals));
for i=1:length(this_data_stim)
    these_counts(find(abs(these_xvals-this_data_stim(i))<0.01),find(abs(these_yvals-this_data_recr(i))<0.01)) = these_counts(find(abs(these_xvals-this_data_stim(i))<0.01),find(abs(these_yvals-this_data_recr(i))<0.01)) + 1;
end
this_max = nanmax(these_counts(:))

for x=1:length(these_xvals)
    for y=1:length(these_yvals)
        if these_counts(x,y)==7
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+0/9*(1-p.col.seq))
        elseif these_counts(x,y)==6
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+1/9*(1-p.col.seq))
        elseif these_counts(x,y)==5
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+2/9*(1-p.col.seq))
        elseif these_counts(x,y)==4
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+3/9*(1-p.col.seq))
        elseif these_counts(x,y)==3
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+4/9*(1-p.col.seq))
        elseif these_counts(x,y)==2
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+5/9*(1-p.col.seq))
        elseif these_counts(x,y)==1
            scatter(these_xvals(x),these_yvals(y),'.','SizeData',100,'MarkerEdgeColor',p.col.seq+6/9*(1-p.col.seq))
        end
    end
end

xlim([0,5.3])
xticks([0:2.5:5])
xlabel({'Photoactivation time (s)'})
ylim([0,5.3])
yticks([0:2.5:5])
ylabel({'Firing field peak (s)'})

% savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1.fig']);
% saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_recruitedLocation_ipsicontra_stim_d1.pdf']); set(gcf,'Color',[1,1,1])



%% Internal revision comments

%% Fig5_ImprintingAnalysis_ic_stim_d1

these_labels = categorical({'Same','Opposite'});
these_labels = reordercats(these_labels,{'Same','Opposite'});
F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_ipsi_seq = nanmean(prop_respIpsiSeqCells_respSeqCells_seq(:,[2,4]),2); %nanmean(prop_respIpsiSeqCells_respCells_seq(:,[2,4]),2);
this_data_ipsi_ctrl = nanmean(prop_respIpsiSeqCells_respSeqCells_ctrl(:,[2,4]),2); %nanmean(prop_respIpsiSeqCells_respCells_ctrl(:,[2,4]),2);
this_data_contra_seq = nanmean(prop_respContraSeqCells_respSeqCells_seq(:,[2,4]),2); %nanmean(prop_respContraSeqCells_respCells_seq(:,[2,4]),2);
this_data_contra_ctrl = nanmean(prop_respContraSeqCells_respSeqCells_ctrl(:,[2,4]),2); %nanmean(prop_respContraSeqCells_respCells_ctrl(:,[2,4]),2);

this_data_ipsi = [this_data_ipsi_seq;this_data_ipsi_ctrl];
this_data_contra = [this_data_contra_seq;this_data_contra_ctrl];
v = bar(these_labels,nanmean([this_data_ipsi,this_data_contra],1)*100);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.darkGray; v.EdgeColor = 'none'; 
plot(these_labels,[this_data_ipsi_ctrl,this_data_contra_ctrl]*100,'-','Color',p.col.ctrl,'LineWidth',1)
plot(these_labels,[this_data_ipsi_seq,this_data_contra_seq]*100,'-','Color',p.col.seq,'LineWidth',1)
ylim([0,100])
yticks([0,50,100])
ytickformat('percentage')
ylabel({'Sequence participation probability ???'})

[temp1,~,temp2] = signrank(this_data_ipsi,this_data_contra);

% savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_ic_stim_d1.fig']);
% saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_ic_stim_d1.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_ic_stim_d1.pdf']); set(gcf,'Color',[1,1,1])




























function [resp] = responseAnalysis(info,iscell,ops,p,path,trg,task,s2p_meta,thor_beh,act)
% act = dFFn_beh; p = get_p;

%% Preparations

resp.stimType = info.stimType;

% pre-process activity measure
if p.resp.zscore
% 	act_smoothed = smoothdata(act,2,'gaussian',3*5);
%     this_mean = nanmean(act_smoothed,2);
%     this_std = nanstd(act_smoothed,[],2);
%     act = (act - this_mean) ./ this_std;
    act = nanzscore(act,[],2);
end

% extract important parameters
numRois = size(act,1);
numClusters = size(trg.clustering,2);
numTrialsPerGroup = trg.p.numTrialsPerGroup;
numPatternsPerSequence = size(thor_beh.stimPatternOnset,1);
numStimSequences = size(thor_beh.stimPatternOnset,2);

% store used curation and parameters
resp.iscell_used = iscell;
resp.p = p;
resp.p.maxTargetCentroidDistance_pix = p.resp.maxTargetCentroidDistance_um*info.scope.fovSize_pix/info.scope.fovSize_um;
resp.dist_closestLaser = trg.dist_closestLaser;


%% Get activity during analysis windows

disp('--- Extracting activity in analysis windows.')

% identify baseline and response window frames
resp.frames_base = nan(numPatternsPerSequence,numStimSequences,p.resp.numFrames);
resp.frames_resp = nan(numPatternsPerSequence,numStimSequences,p.resp.numFrames);
for i=1:numPatternsPerSequence
    for j=1:numStimSequences
        resp.frames_base(i,j,:) = thor_beh.stimPatternOnset(i,j)-p.resp.numFrames:thor_beh.stimPatternOnset(i,j)-1;
        resp.frames_resp(i,j,:) = thor_beh.stimPatternComplete(i,j)+1:thor_beh.stimPatternComplete(i,j)+p.resp.numFrames;
    end
end

% identify baseline and response window frames for no-stim trials
% noStimPatternOnset = repmat(thor_beh.sync(find(trg.seq(2,:)==0)),numPatternsPerSequence,1) + round(mean(thor_beh.stimPatternOnset - thor_beh.sync(find(trg.seq(2,:))),2));
% noStimPatternComplete = repmat(thor_beh.sync(find(trg.seq(2,:)==0)),numPatternsPerSequence,1) + round(mean(thor_beh.stimPatternComplete - thor_beh.sync(find(trg.seq(2,:))),2));
% resp.noStim.frames_base = nan(numPatternsPerSequence,size(noStimPatternOnset,2),p.resp.numFrames);
% resp.noStim.frames_resp = nan(numPatternsPerSequence,size(noStimPatternOnset,2),p.resp.numFrames);
% for i=1:numPatternsPerSequence
%     for j=1:size(noStimPatternOnset,2)
%         resp.noStim.frames_base(i,j,:) = noStimPatternOnset(i,j)-p.resp.numFrames:noStimPatternOnset(i,j)-1;
%         resp.noStim.frames_resp(i,j,:) = noStimPatternComplete(i,j)+1:noStimPatternComplete(i,j)+p.resp.numFrames;
%     end
% end

% get mean activity in analysis windows
act_base = nan(numRois,numPatternsPerSequence,numStimSequences);
act_resp = nan(numRois,numPatternsPerSequence,numStimSequences);
for i=1:numPatternsPerSequence
    for j=1:numStimSequences
        try
            if p.resp.subtractBlockwiseTuning
                this_block = floor(j/trg.p.numStimTrialsPerBlock);
                this_seq =  trg.seq(:,(this_block-1)*trg.info.trialsPerBlock+1:this_block*trg.info.trialsPerBlock);
                temp = find(trg.seq(2,:)==1);
                if (trg.seq(1,temp(j))==1 | trg.seq(1,temp(j))==3)
                    [~,temp,~] = intersect(find(this_seq(2,:)==0),intersect(find(this_seq(2,:)==0), find(this_seq(1,:)==1 | this_seq(1,:)==3)));
                elseif (trg.seq(1,temp(j))==2 | trg.seq(1,temp(j))==4)
                    [~,temp,~] = intersect(find(this_seq(2,:)==0),intersect(find(this_seq(2,:)==0), find(this_seq(1,:)==2 | this_seq(1,:)==4)));
                end
                these_noStimTrials =  temp + (this_block-1)*(trg.info.trialsPerBlock-trg.p.numStimTrialsPerBlock);
                temp1 = nanmean([nanmean(act(:,squeeze(resp.noStim.frames_base(i,these_noStimTrials(1),:))),2),nanmean(act(:,squeeze(resp.noStim.frames_base(i,these_noStimTrials(2),:))),2)],2);
                temp2 = nanmean([nanmean(act(:,squeeze(resp.noStim.frames_resp(i,these_noStimTrials(1),:))),2),nanmean(act(:,squeeze(resp.noStim.frames_resp(i,these_noStimTrials(2),:))),2)],2);
                act_base(:,i,j) = nanmean(act(:,squeeze(resp.frames_base(i,j,:))),2) - temp1;
                act_resp(:,i,j) = nanmean(act(:,squeeze(resp.frames_resp(i,j,:))),2) - temp2;
            elseif p.resp.subtractOverallTuning
                temp = find(trg.seq(2,:)==1);
                if (trg.seq(1,temp(j))==1 | trg.seq(1,temp(j))==3)
                    temp = find(trg.seq(2,:)==0 & (trg.seq(1,:)==1 | trg.seq(1,:)==3));
                elseif (trg.seq(1,temp(j))==2 | trg.seq(1,temp(j))==4)
                    temp = find(trg.seq(2,:)==0 & (trg.seq(1,:)==2 | trg.seq(1,:)==4));
                end
                [~,these_noStimTrials] = intersect(find(trg.seq(2,:)==0),temp);
                temp1 = nanmean(act(:,squeeze(resp.noStim.frames_base(i,these_noStimTrials,:))),2);
                temp2 = nanmean(act(:,squeeze(resp.noStim.frames_resp(i,these_noStimTrials,:))),2);
                act_base(:,i,j) = nanmean(act(:,squeeze(resp.frames_base(i,j,:))),2) - temp1;
                act_resp(:,i,j) = nanmean(act(:,squeeze(resp.frames_resp(i,j,:))),2) - temp2;
            else
                act_base(:,i,j) = nanmean(act(:,squeeze(resp.frames_base(i,j,:))),2);
                act_resp(:,i,j) = nanmean(act(:,squeeze(resp.frames_resp(i,j,:))),2);
            end
        catch
            act_base(:,i,j) = NaN;
            act_resp(:,i,j) = NaN;      
        end
    end
end

% get mean activity in analysis windows for no-stim trials
% act_base_noStim = nan(numRois,numPatternsPerSequence,size(noStimPatternOnset,2));
% act_resp_noStim = nan(numRois,numPatternsPerSequence,size(noStimPatternOnset,2));
% for i=1:numPatternsPerSequence
%     for j=1:size(noStimPatternOnset,2)
%         try
%             act_base_noStim(:,i,j) = nanmean(act(:,squeeze(resp.noStim.frames_base(i,j,:))),2);
%             act_resp_noStim(:,i,j) = nanmean(act(:,squeeze(resp.noStim.frames_resp(i,j,:))),2);
%         catch
%             act_base_noStim(:,i,j) = NaN;
%             act_resp_noStim(:,i,j) = NaN;      
%         end
%     end
% end

% separate and sort pattern responses by cluster
resp.act_base = nan(numRois,numClusters,numTrialsPerGroup);
resp.act_resp = nan(numRois,numClusters,numTrialsPerGroup);
resp.act_net = nan(numRois,numClusters,numTrialsPerGroup);
resp.firstCluster.clusterIdx = nan(numTrialsPerGroup,2);
resp.stimTrial2typeAtrial = nan(numStimSequences,1);
resp.stimTrial2typeXtrial = nan(numStimSequences,1);
resp.typeAtrial2stimTrial = nan(numTrialsPerGroup,1);
resp.typeXtrial2stimTrial = nan(numTrialsPerGroup,1);
stimTrials = nan(1,2*numTrialsPerGroup);
stimTrials(1:length(task.type(task.var==1))) = task.type(task.var==1);
this_stim1TrialCount = 0;
this_stim2TrialCount = 0;
for j=1:numStimSequences
    if ~isnan(stimTrials(j))

        this_sequence = trg.sequenceOrder(j);
        if (stimTrials(j)==1 | stimTrials(j)==3)
            this_stim1TrialCount = this_stim1TrialCount+1;
            this_trialIdx = this_stim1TrialCount;
            resp.firstCluster.clusterIdx(this_trialIdx,1) = trg.sequenceClusters(1,this_sequence);
            resp.typeAtrial2stimTrial(this_trialIdx) = j;
            resp.stimTrial2typeAtrial(j) = this_stim1TrialCount;
        elseif (stimTrials(j)==2 | stimTrials(j)==4)
            this_stim2TrialCount = this_stim2TrialCount+1;
            this_trialIdx = this_stim2TrialCount;
            resp.firstCluster.clusterIdx(this_trialIdx,2) = trg.sequenceClusters(1,this_sequence);
            resp.typeXtrial2stimTrial(this_trialIdx) = j;
            resp.stimTrial2typeXtrial(j) = this_stim2TrialCount;
        end

        these_clusters = trg.sequenceClusters(:,this_sequence);
        for i=1:numPatternsPerSequence
            resp.act_base(:,these_clusters(i),this_trialIdx) = act_base(:,i,j);
            resp.act_resp(:,these_clusters(i),this_trialIdx) = act_resp(:,i,j);
            resp.act_net(:,these_clusters(i),this_trialIdx) = resp.act_resp(:,these_clusters(i),this_trialIdx) - resp.act_base(:,these_clusters(i),this_trialIdx);
        end
    end
end

% separate and sort pattern responses by cluster for no-stim trials
% resp.noStim.act_base = nan(numRois,numClusters,size(noStimPatternOnset,2)/2);
% resp.noStim.act_resp = nan(numRois,numClusters,size(noStimPatternOnset,2)/2);
% resp.noStim.act_net = nan(numRois,numClusters,size(noStimPatternOnset,2)/2);
% noStimTrials = nan(1,size(noStimPatternOnset,2));
% noStimTrials(1:length(task.type(task.var==0))) = task.type(task.var==0);
% this_nostim1TrialCount = 0;
% this_nostim2TrialCount = 0;
% for j=1:size(noStimPatternOnset,2)
%     if ~isnan(noStimTrials(j))
%         if (noStimTrials(j)==1 | noStimTrials(j)==3)
%             this_nostim1TrialCount = this_nostim1TrialCount+1;
%             this_trialIdx = this_nostim1TrialCount;
%             if strcmp(resp.stimType,'seq')
%                 these_clusters = trg.sequenceClusters(:,1);
%             elseif strcmp(resp.stimType,'ctrl')
%                 temp = trg.sequenceClusters(:,1);
%                 these_clusters = temp(randperm(length(temp)));
%             end
%         elseif (noStimTrials(j)==2 | noStimTrials(j)==4)
%             this_nostim2TrialCount = this_nostim2TrialCount+1;
%             this_trialIdx = this_nostim2TrialCount;
%             if strcmp(resp.stimType,'seq')
%                 these_clusters = trg.sequenceClusters(:,2);
%             elseif strcmp(resp.stimType,'ctrl')
%                 temp = trg.sequenceClusters(:,2);
%                 these_clusters = temp(randperm(length(temp)));
%             end
%         end
%         
%         for i=1:numPatternsPerSequence
%             resp.noStim.act_base(:,these_clusters(i),this_trialIdx) = act_base_noStim(:,i,j);
%             resp.noStim.act_resp(:,these_clusters(i),this_trialIdx) = act_resp_noStim(:,i,j);
%             resp.noStim.act_net(:,these_clusters(i),this_trialIdx) = resp.noStim.act_resp(:,these_clusters(i),this_trialIdx) - resp.noStim.act_base(:,these_clusters(i),this_trialIdx);
%         end
%     end
% end

% calculate trial-averaged activities
resp.avgAct_base = nanmean(resp.act_base,3);
resp.avgAct_resp = nanmean(resp.act_resp,3);
resp.avgAct_net = nanmean(resp.act_net,3);

% calculate trial-averaged activities for no-stim trials
% resp.noStim.avgAct_base = nanmean(resp.noStim.act_base,3);
% resp.noStim.avgAct_resp = nanmean(resp.noStim.act_resp,3);
% resp.noStim.avgAct_net = nanmean(resp.noStim.act_net,3);

% % find trials were a given cluster was first
% resp.firstCluster.ipsiStimTrials = cell(numClusters,1);
% for i=1:numClusters
%     temp = resp.firstCluster.clusterIdx==i;
%     try
%         resp.firstCluster.ipsiStimTrials{i} = find(resp.firstCluster.clusterIdx(:,find(sum(temp,1)))==i);
%     catch
%         resp.firstCluster.ipsiStimTrials{i} = {};
%     end
% end
% 
% % isolate pattern responses for first clusters
% resp.firstCluster.act_base = nan(size(resp.act_base));
% resp.firstCluster.act_resp = nan(size(resp.act_resp));
% resp.firstCluster.act_net = nan(size(resp.act_net));
% for i=1:numClusters
%     try
%         resp.firstCluster.act_base(:,i,resp.firstCluster.ipsiStimTrials{i}) = resp.act_base(:,i,resp.firstCluster.ipsiStimTrials{i});
%         resp.firstCluster.act_resp(:,i,resp.firstCluster.ipsiStimTrials{i}) = resp.act_resp(:,i,resp.firstCluster.ipsiStimTrials{i});
%         resp.firstCluster.act_net(:,i,resp.firstCluster.ipsiStimTrials{i}) = resp.act_net(:,i,resp.firstCluster.ipsiStimTrials{i});
%     catch
%     end
% end
% resp.firstCluster.avgAct_base = nanmean(resp.firstCluster.act_base,3);
% resp.firstCluster.avgAct_resp = nanmean(resp.firstCluster.act_resp,3);
% resp.firstCluster.avgAct_net = nanmean(resp.firstCluster.act_net,3);


%% Statstical tests

disp('--- Identifying responders.')

% calculate p-values
resp.stats_p = nan(numRois,numClusters);
for n=1:numRois
    if iscell(n)
        for i=1:numClusters
            this_act_base = squeeze(resp.act_base(n,i,:));
            this_act_resp = squeeze(resp.act_resp(n,i,:));
            if sum(~isnan(this_act_base)) && sum(~isnan(this_act_resp))
                resp.stats_p(n,i) = signrank(this_act_base,this_act_resp);
            end
        end
    end
end
if ~isempty(find(resp.stats_p==0))
    warning('Something went wrong with p-value calculation.')
end

% calculate fdr and q-values
resp.stats_fdr = nan(numRois,numClusters);
resp.stats_q = nan(numRois,numClusters);
resp.stats_priori = nan(1,numClusters);
for i=1:numClusters
    [resp.stats_fdr(:,i),resp.stats_q(:,i),resp.stats_priori(i)] = mafdr(resp.stats_p(:,i)); % double-checked that it works correctly with NaNs
end
if any(resp.stats_priori>=1)
    warning('The priori for at least one cluster is greater than or equal to 1.')
end

% filter out 
if excludeOutsideExclusionZone
    for i=1:numClusters
        resp.stats_q(trg.dist_closestLaser(:,i)< p.resp.exclusionRadius,i) = NaN; % exclude neurons within 50 um radius
        %flw.stats_q(resp.responders_main.idcs_resp{i},i) = NaN; % exclude responders
        resp.stats_q(resp.responders_main.idcs_targeted{i},i) = NaN; % exclude targeted neurons
    end
end

% find significant followers/responders
resp.stats_sig = double(resp.stats_q<p.resp.qthres);
resp.stats_sig(isnan(resp.stats_q))=NaN;
resp.stats_sig_pos = double(double(resp.avgAct_net>0)+resp.stats_sig==2);
resp.stats_sig_pos(isnan(resp.stats_q))=NaN;
resp.stats_sig_neg = double(double(resp.avgAct_net<0)+resp.stats_sig==2);
resp.stats_sig_neg(isnan(resp.stats_q))=NaN;

resp.idcs_sig_neu = cell(1,numClusters);
resp.idcs_sig_pos = cell(1,numClusters);
resp.idcs_sig_neg = cell(1,numClusters);
for i=1:numClusters
    resp.idcs_sig_neu{i} = find(resp.stats_sig(:,i)==0);
    resp.idcs_sig_pos{i} = find(resp.stats_sig_pos(:,i)==1);
    resp.idcs_sig_neg{i} = find(resp.stats_sig_neg(:,i)==1);
end

% --- same for first clusters

% resp.firstCluster.stats_p = nan(numRois,numClusters);
% for n=1:numRois
%     if iscell(n)
%         for i=1:numClusters
%             this_act_base = rmmissing(squeeze(resp.firstCluster.act_base(n,i,:)));
%             this_act_resp = rmmissing(squeeze(resp.firstCluster.act_resp(n,i,:))); % HERE I NEED TO SUBSAMPLE TO MIN NUMBER OF TRIALS FOR SEQ CTRL COMPARISON!!!
%             if sum(~isnan(this_act_base)) && sum(~isnan(this_act_resp))
%                 resp.firstCluster.stats_p(n,i) = signrank(this_act_base,this_act_resp);
%             end
%         end
%     end
% end
% if ~isempty(find(resp.firstCluster.stats_p==0))
%     warning('Something went wrong with p-value calculation.')
% end
% 
% resp.firstCluster.stats_fdr = nan(numRois,numClusters);
% resp.firstCluster.stats_q = nan(numRois,numClusters);
% resp.firstCluster.stats_priori = nan(1,numClusters);
% for i=1:numClusters
%     if ~all(isnan(resp.firstCluster.stats_p(:,i)))
%         [resp.firstCluster.stats_fdr(:,i),resp.firstCluster.stats_q(:,i),resp.firstCluster.stats_priori(i)] = mafdr(resp.firstCluster.stats_p(:,i));
%     end
% end
% if any(resp.firstCluster.stats_priori>=1)
%     warning('The priori for at least one first cluster is greater than or equal to 1.')
% end
% 
% resp.firstCluster.stats_sig = double(resp.firstCluster.stats_q<p.resp.qthres);
% resp.firstCluster.stats_sig(isnan(resp.firstCluster.stats_q))=NaN;
% resp.firstCluster.stats_sig_pos = double(double(resp.firstCluster.avgAct_net>0)+resp.firstCluster.stats_sig==2);
% resp.firstCluster.stats_sig_pos(isnan(resp.firstCluster.stats_q))=NaN;
% resp.firstCluster.stats_sig_neg = double(double(resp.firstCluster.avgAct_net<0)+resp.firstCluster.stats_sig==2);
% resp.firstCluster.stats_sig_neg(isnan(resp.firstCluster.stats_q))=NaN;

% --- same for no-stim

% resp.noStim.stats_p = nan(numRois,numClusters);
% for n=1:numRois
%     if iscell(n)
%         for i=1:numClusters
%             this_act_base = rmmissing(squeeze(resp.noStim.act_base(n,i,:)));
%             this_act_resp = rmmissing(squeeze(resp.noStim.act_resp(n,i,:)));
%             if sum(~isnan(this_act_base)) && sum(~isnan(this_act_resp))
%                 resp.noStim.stats_p(n,i) = signrank(this_act_base,this_act_resp);
%             end
%         end
%     end
% end
% if ~isempty(find(resp.noStim.stats_p==0))
%     warning('Something went wrong with p-value calculation.')
% end
% 
% resp.noStim.stats_fdr = nan(numRois,numClusters);
% resp.noStim.stats_q = nan(numRois,numClusters);
% resp.noStim.stats_priori = nan(1,numClusters);
% for i=1:numClusters
%     if ~all(isnan(resp.noStim.stats_p(:,i)))
%         [resp.noStim.stats_fdr(:,i),resp.noStim.stats_q(:,i),resp.noStim.stats_priori(i)] = mafdr(resp.noStim.stats_p(:,i));
%     end
% end
% if any(resp.noStim.stats_priori>=1)
%     warning('The priori for at least one cluster is greater than or equal to 1.')
% end
% 
% resp.noStim.stats_sig = double(resp.noStim.stats_q<p.resp.qthres);
% resp.noStim.stats_sig(isnan(resp.noStim.stats_q))=NaN;
% resp.noStim.stats_sig_pos = double(double(resp.noStim.avgAct_net>0)+resp.noStim.stats_sig==2);
% resp.noStim.stats_sig_pos(isnan(resp.noStim.stats_q))=NaN;
% resp.noStim.stats_sig_neg = double(double(resp.noStim.avgAct_net<0)+resp.noStim.stats_sig==2);
% resp.noStim.stats_sig_neg(isnan(resp.noStim.stats_q))=NaN;
% 
% resp.noStim.idcs_sig_neu = cell(1,numClusters);
% resp.noStim.idcs_sig_pos = cell(1,numClusters);
% resp.noStim.idcs_sig_neg = cell(1,numClusters);
% for i=1:numClusters
%     resp.noStim.idcs_sig_neu{i} = find(resp.noStim.stats_sig(:,i)==0);
%     resp.noStim.idcs_sig_pos{i} = find(resp.noStim.stats_sig_pos(:,i)==1);
%     resp.noStim.idcs_sig_neg{i} = find(resp.noStim.stats_sig_neg(:,i)==1);
% end


%% Apply criteria

for i=1:length(p.resp.amplitudeCritiera)
    resp = applyResponderCriteria(trg,resp,s2p_meta,p.resp.amplitudeCritiera(i),0);
    resp = applyResponderCriteria(trg,resp,s2p_meta,p.resp.amplitudeCritiera(i),1);
end
resp.responders_main = resp.(p.resp.mainResponderCondition);

% for i=1:length(p.resp.amplitudeCritiera)
%     resp = applyResponderCriteria_firstCluster(trg,resp,s2p_meta,p.resp.amplitudeCritiera(i),0);
%     resp = applyResponderCriteria_firstCluster(trg,resp,s2p_meta,p.resp.amplitudeCritiera(i),1);
% end
% resp.firstCluster.responders_main = resp.firstCluster.(p.resp.mainResponderCondition);

% for i=1:length(p.resp.amplitudeCritiera)
%     resp = applyResponderCriteria_noStim(trg,resp,s2p_meta,p.resp.amplitudeCritiera(i),0);
%     resp = applyResponderCriteria_noStim(trg,resp,s2p_meta,p.resp.amplitudeCritiera(i),1);
% end
% resp.noStim.responders_main = resp.noStim.(p.resp.mainResponderCondition);


%% Photostimulation stability

n = 0;
resp.trgRespAmps_cells = nan(size(resp.responders_main.targetedCells,1)*size(resp.responders_main.targetedCells,2),size(resp.act_net,3));
resp.trgRespAmps_clusters = nan(size(resp.responders_main.targetedCells,2),size(resp.act_net,3));
for j=1:size(resp.responders_main.targetedCells,2)
    for i=1:size(resp.responders_main.targetedCells,1)
        n=n+1;
        if ~isnan(resp.responders_main.targetedCells(i,j))
            resp.trgRespAmps_cells(n,:) = squeeze(resp.act_net(resp.responders_main.targetedCells(i,j),j,:));
        end
    end
    resp.trgRespAmps_clusters(j,:) = nanmean(resp.trgRespAmps_cells(n-size(resp.responders_main.targetedCells,1)+1:n,:),1);
end
temp = reshape(resp.trgRespAmps_cells,[size(resp.trgRespAmps_cells,1),((info.task.trialsPerBlock/2)*(info.task.numStimTrials/info.task.numTrials)),info.task.numBlocks]);
resp.trgRespAmps_cells_bw = squeeze(nanmean(temp,2));
temp = reshape(resp.trgRespAmps_clusters,[size(resp.trgRespAmps_clusters,1),((info.task.trialsPerBlock/2)*(info.task.numStimTrials/info.task.numTrials)),info.task.numBlocks]);
resp.trgRespAmps_clusters_bw = squeeze(nanmean(temp,2));


%% Save

resp = orderfields(resp);
save([path.filepart_out,'resp.mat'],'resp','-v7.3');
disp(['--- Saved resp file as ',[path.filepart_out,'resp.mat'],'.'])


%% Make figures

% cluster responses maps
if ~ops.resp.skip_clusterResponseMaps
    if ~exist([path.filepart_outX,'plots\ClusterResponseMaps'],'dir')
        mkdir([path.filepart_outX,'plots\ClusterResponseMaps']);
    end
    for i=1:size(resp.responders_main.targetedCells,2)
        F = clusterResponseMap(s2p_meta,iscell,trg,resp,i,0); 
        savefig(F,[path.filepart_outX,'plots\ClusterResponseMaps\',info.animal,'_',info.date,'_','cluster_',num2str(i),'.fig']);
        saveas(F,[path.filepart_outX,'plots\ClusterResponseMaps\',info.animal,'_',info.date,'_','cluster_',num2str(i),'.png']);
    end
end

% response analysis figure
F = responseAnalysisFigure(s2p_meta,iscell,trg,resp,p,info);
% F = responseAnalysisFigure_reduced(s2p_meta,iscell,trg,resp,p,info);
savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','responseAnalysis.fig']);
saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','responseAnalysis.png']);
disp(['--- Saved response analysis figures to ',path.filepart_outX,'plots.'])

% response analysis figure for no-stim trials
% F = noStimResponseFigure(s2p_meta,iscell,trg,resp,p,info);
% savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','noStimResponse.fig']);
% saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','noStimResponse.png']);
% disp(['--- Saved no-stim response figures to ',path.filepart_outX,'plots.'])

% first cluster response figure
% F = firstClusterResponseFigure(s2p_meta,iscell,trg,resp,p,info);
% savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','firstClusterResponse.fig']);
% saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','firstClusterResponse.png']);
% disp(['--- Saved first cluster response figures to ',path.filepart_outX,'plots.'])
% if ~ops.resp.skip_clusterResponseMaps
%     if ~exist([path.filepart_outX,'plots\ClusterResponseMaps'],'dir')
%         mkdir([path.filepart_outX,'plots\ClusterResponseMaps']);
%     end
%     if strcmp(infio.stimType,'seq')
%         i = resp.firstCluster.clusterIdx(1,1);
%         F = clusterResponseMap(s2p_meta,iscell,trg,resp,i,0,1); 
%         savefig(F,[path.filepart_outX,'plots\ClusterResponseMaps\',info.animal,'_',info.date,'_','cluster_',num2str(i),'_firstClusterA.fig']);
%         saveas(F,[path.filepart_outX,'plots\ClusterResponseMaps\',info.animal,'_',info.date,'_','cluster_',num2str(i),'_firstClusterA.png']);
%         i = resp.firstCluster.clusterIdx(1,2);
%         F = clusterResponseMap(s2p_meta,iscell,trg,resp,i,0,1); 
%         savefig(F,[path.filepart_outX,'plots\ClusterResponseMaps\',info.animal,'_',info.date,'_','cluster_',num2str(i),'_firstClusterX.fig']);
%         saveas(F,[path.filepart_outX,'plots\ClusterResponseMaps\',info.animal,'_',info.date,'_','cluster_',num2str(i),'_firstClusterX.png']);
%     end
% end


%% Follower analysis

disp('--- Removing exclusion zones.')

% % store used curation and parameters
flw.iscell_used = iscell;
flw.p = p;
flw.p.maxTargetCentroidDistance_pix = p.resp.maxTargetCentroidDistance_um*info.scope.fovSize_pix/info.scope.fovSize_um;

% get activities (with exclusion zone from this and previous clusters removed)
flw.inExclusionZone = nan(numRois,numClusters,numTrialsPerGroup);
flw.act_base = nan(numRois,numClusters,numTrialsPerGroup);
flw.act_resp = nan(numRois,numClusters,numTrialsPerGroup);
flw.act_net = nan(numRois,numClusters,numTrialsPerGroup);
flw.numInclTrials = nan(numRois,numClusters);
flw.numInclTrials_afterThresholding = nan(numRois,numClusters);
for n=1:numRois
    if iscell(n)
        for i=1:numClusters
            
            % Arnd method (subtracting average cluster response)
            if false % subtracting average baseline and responses
                flw.act_base(n,i,:) = resp.act_base(n,i,:) - resp.avgAct_base(n,i);
                flw.act_resp(n,i,:) = resp.act_resp(n,i,:) - resp.avgAct_resp(n,i);
                flw.act_net(n,i,:) = resp.act_net(n,i,:) - resp.avgAct_net(n,i);
            elseif true % subtracting average baseline from baseline and response
                flw.act_base(n,i,:) = resp.act_base(n,i,:) - resp.avgAct_base(n,i);
                flw.act_resp(n,i,:) = resp.act_resp(n,i,:) - resp.avgAct_base(n,i);
                flw.act_net(n,i,:) = resp.act_net(n,i,:) - resp.avgAct_base(n,i);
            else
                flw.act_base(n,i,:) = resp.act_base(n,i,:);
                flw.act_resp(n,i,:) = resp.act_resp(n,i,:);
                flw.act_net(n,i,:) = resp.act_net(n,i,:);
            end
            
            % trial exclusions
            if false
                if find(sum(trg.grouping==i,1))==1
                    this_clusterType = "A";
                elseif find(sum(trg.grouping==i,1))==2
                    this_clusterType = "X";
                else
                    error('Something went wrong')
                end

                % exclude trial for neuron if it is in respective exclusion zone
                for k=1:numTrialsPerGroup
                    if this_clusterType=="A"
                        this_stimTrial = resp.typeAtrial2stimTrial(k);
                    elseif this_clusterType=="X"
                        this_stimTrial = resp.typeXtrial2stimTrial(k);
                    end
                    this_sequence = trg.sequenceOrder(this_stimTrial);
                    these_sequenceClusters = trg.sequenceClusters(:,this_sequence);
                    this_clusterRank = find(these_sequenceClusters==i);
                    if p.resp.excludeInSelfCluster
                        temp = this_clusterRank-p.resp.previousClustersForExclusionZone:this_clusterRank;
                    else
                        temp = this_clusterRank-p.resp.previousClustersForExclusionZone:this_clusterRank-1;
                    end
                    temp = temp(temp>0);
                    these_clustersForExclusion = these_sequenceClusters(temp);
                    flw.inExclusionZone(n,i,k) = any(trg.dist_closestLaser(n,these_clustersForExclusion) < p.resp.exclusionRadius);
                    if flw.inExclusionZone(n,i,k)
                        flw.act_base(n,i,k) = NaN;
                        flw.act_resp(n,i,k) = NaN;
                        flw.act_net(n,i,k) = NaN;
                    end
                end
                flw.numInclTrials(n,i) = length(rmmissing(squeeze(flw.act_net(n,i,:))));
                flw.numInclTrials_afterThresholding(n,i) = flw.numInclTrials(n,i);
                if flw.numInclTrials(n,i) < p.resp.minNumTrials
                    flw.numInclTrials_afterThresholding(n,i) = NaN;
                    flw.act_base(n,i,:) = NaN;
                    flw.act_resp(n,i,:) = NaN;
                    flw.act_net(n,i,:) = NaN;
                end
            end
        end
    end
end
flw.avgAct_base = nanmean(flw.act_base,3);
flw.avgAct_resp = nanmean(flw.act_resp,3);
flw.avgAct_net = nanmean(flw.act_net,3);


%% Statstical tests

disp('--- Identifying followers.')

% calculate p-values
flw.stats_p = nan(numRois,numClusters);
for n=1:numRois
    if iscell(n)
        for i=1:numClusters
            this_act_base = squeeze(flw.act_base(n,i,:));
            this_act_resp = squeeze(flw.act_resp(n,i,:));
            if (sum(~isnan(this_act_base)) && sum(~isnan(this_act_resp)))
                flw.stats_p(n,i) = signrank(this_act_base,this_act_resp);
            end
        end
    end
end

% calculate fdr and q-values
flw.stats_fdr = nan(numRois,numClusters);
flw.stats_q = nan(numRois,numClusters);
flw.stats_priori = nan(1,numClusters);
for i=1:numClusters
    [flw.stats_fdr(:,i),flw.stats_q_raw(:,i),flw.stats_priori(i)] = mafdr(flw.stats_p(:,i)); % double-checked that it works correctly with NaNs
end
if any(flw.stats_priori>=1)
    warning('The priori for at least one cluster is greater than or equal to 1.')
end

% filter out 
flw.stats_q = flw.stats_q_raw;
for i=1:numClusters
    flw.stats_q(trg.dist_closestLaser(:,i)< p.resp.exclusionRadius,i) = NaN; % exclude neurons within 50 um radius
    %flw.stats_q(resp.responders_main.idcs_resp{i},i) = NaN; % exclude responders
    flw.stats_q(resp.responders_main.idcs_targeted{i},i) = NaN; % exclude targeted neurons
end

% find significant followers/responders
flw.stats_sig = double(flw.stats_q<p.resp.qthres);
flw.stats_sig(isnan(flw.stats_q))=NaN;
flw.stats_sig_pos = double(double(flw.avgAct_net>0)+flw.stats_sig==2);
flw.stats_sig_pos(isnan(flw.stats_q))=NaN;
flw.stats_sig_neg = double(double(flw.avgAct_net<0)+flw.stats_sig==2);
flw.stats_sig_neg(isnan(flw.stats_q))=NaN;

flw.idcs_sig_neu = cell(1,numClusters);
flw.idcs_sig_pos = cell(1,numClusters);
flw.idcs_sig_neg = cell(1,numClusters);
for i=1:numClusters
    flw.idcs_sig_neu{i} = find(flw.stats_sig(:,i)==0);
    flw.idcs_sig_pos{i} = find(flw.stats_sig_pos(:,i)==1);
    flw.idcs_sig_neg{i} = find(flw.stats_sig_neg(:,i)==1);
end


%% Apply criteria

for i=1:length(p.resp.amplitudeCritiera)
    flw = applyResponderCriteria(trg,flw,s2p_meta,p.resp.amplitudeCritiera(i),0);
    flw = applyResponderCriteria(trg,flw,s2p_meta,p.resp.amplitudeCritiera(i),1);
end
flw.responders_main = flw.(p.resp.mainResponderCondition);


%% Photostimulation stability

n = 0;
flw.trgRespAmps_cells = nan(size(flw.responders_main.targetedCells,1)*size(flw.responders_main.targetedCells,2),size(flw.act_net,3));
flw.trgRespAmps_clusters = nan(size(flw.responders_main.targetedCells,2),size(flw.act_net,3));
for j=1:size(flw.responders_main.targetedCells,2)
    for i=1:size(flw.responders_main.targetedCells,1)
        n=n+1;
        if ~isnan(flw.responders_main.targetedCells(i,j))
            flw.trgRespAmps_cells(n,:) = squeeze(flw.act_net(flw.responders_main.targetedCells(i,j),j,:));
        end
    end
    flw.trgRespAmps_clusters(j,:) = nanmean(flw.trgRespAmps_cells(n-size(flw.responders_main.targetedCells,1)+1:n,:),1);
end
temp = reshape(flw.trgRespAmps_cells,[size(flw.trgRespAmps_cells,1),((info.task.trialsPerBlock/2)*(info.task.numStimTrials/info.task.numTrials)),info.task.numBlocks]);
flw.trgRespAmps_cells_bw = squeeze(nanmean(temp,2));
temp = reshape(flw.trgRespAmps_clusters,[size(flw.trgRespAmps_clusters,1),((info.task.trialsPerBlock/2)*(info.task.numStimTrials/info.task.numTrials)),info.task.numBlocks]);
flw.trgRespAmps_clusters_bw = squeeze(nanmean(temp,2));


%% Save

flw = orderfields(flw);
save([path.filepart_out,'flw.mat'],'flw','-v7.3');
disp(['--- Saved flw file as ',[path.filepart_out,'flw.mat'],'.'])


%% Make figures

% cluster responses maps
if ~ops.resp.skip_clusterResponseMaps
    if ~exist([path.filepart_outX,'plots\ClusterResponseMaps'],'dir')
        mkdir([path.filepart_outX,'plots\ClusterResponseMaps']);
    end
    for i=1:size(flw.responders_main.targetedCells,2)
        F = clusterResponseMap(s2p_meta,iscell,trg,flw,i,0); 
        savefig(F,[path.filepart_outX,'plots\ClusterResponseMaps\',info.animal,'_',info.date,'_','cluster_',num2str(i),'_flw.fig']);
        saveas(F,[path.filepart_outX,'plots\ClusterResponseMaps\',info.animal,'_',info.date,'_','cluster_',num2str(i),'_flw.png']);
    end
end

% response analysis figure
F = responseAnalysisFigure(s2p_meta,iscell,trg,flw,p,info);
savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','responseAnalysis_flw.fig']);
saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','responseAnalysis_flw.png']);
disp(['--- Saved response analysis (flw) figures to ',path.filepart_outX,'plots.'])


%% Follower amplitude by photostimulation response

% based on flw.act_net (not resp.act_net)
act_net = flw.act_net;

% get photostimulation response amplitude trial-by-trial
flw.ampByResp.avgTrgAmp = nan(numTrialsPerGroup,numClusters);
flw.ampByResp.avgRespAmp = nan(numTrialsPerGroup,numClusters);
flw.ampByResp.avgRespTrgAmp = nan(numTrialsPerGroup,numClusters);
for i=1:numClusters
    flw.ampByResp.avgTrgAmp(:,i) = squeeze(nanmean(act_net(resp.responders_main.idcs_targeted{i},i,:),1));
    flw.ampByResp.avgRespAmp(:,i) = squeeze(nanmean(act_net(resp.responders_main.idcs_resp{i},i,:),1));
    flw.ampByResp.avgRespTrgAmp(:,i) = squeeze(nanmean(act_net(resp.responders_main.idcs_respTargeted{i},i,:),1));
end

% get follower amplitude trial-by-trial
flw.ampByResp.avgNeuAmp = nan(numTrialsPerGroup,numClusters);
flw.ampByResp.avgPosAmp = nan(numTrialsPerGroup,numClusters);
flw.ampByResp.avgNegAmp = nan(numTrialsPerGroup,numClusters);
for i=1:numClusters
    flw.ampByResp.avgNeuAmp(:,i) = squeeze(nanmean(act_net(flw.idcs_sig_neu{i},i,:),1));
    flw.ampByResp.avgPosAmp(:,i) = squeeze(nanmean(act_net(flw.idcs_sig_pos{i},i,:),1));
    flw.ampByResp.avgNegAmp(:,i) = squeeze(nanmean(act_net(flw.idcs_sig_neg{i},i,:),1));
end


%% Correlation between follower amplitude and photostimulation response

temp1 = flw.ampByResp.avgPosAmp(:);
temp2 = discretize(flw.ampByResp.avgTrgAmp(:),p.resp.ampBins_edges);
flw.ampByResp.avgPosAmp_by_avgTrgAmpBin = nan(length(p.resp.ampBins_x),1);
for j=1:length(p.resp.ampBins_x)
    if length(find(temp2==j)) >= p.resp.ampBins_minNumDataPoints
        flw.ampByResp.avgPosAmp_by_avgTrgAmpBin(j) = nanmean(temp1(find(temp2==j)));
    end
end

temp1 = flw.ampByResp.avgNeuAmp(:);
temp2 = discretize(flw.ampByResp.avgTrgAmp(:),p.resp.ampBins_edges);
flw.ampByResp.avgNeuAmp_by_avgTrgAmpBin = nan(length(p.resp.ampBins_x),1);
for j=1:length(p.resp.ampBins_x)
    if length(find(temp2==j)) >= p.resp.ampBins_minNumDataPoints
        flw.ampByResp.avgNeuAmp_by_avgTrgAmpBin(j) = nanmean(temp1(find(temp2==j)));
    end
end

temp1 = flw.ampByResp.avgNegAmp(:);
temp2 = discretize(flw.ampByResp.avgTrgAmp(:),p.resp.ampBins_edges);
flw.ampByResp.avgNegAmp_by_avgTrgAmpBin = nan(length(p.resp.ampBins_x),1);
for j=1:length(p.resp.ampBins_x)
    if length(find(temp2==j)) >= p.resp.ampBins_minNumDataPoints
        flw.ampByResp.avgNegAmp_by_avgTrgAmpBin(j) = nanmean(temp1(find(temp2==j)));
    end
end

temp1 = flw.ampByResp.avgPosAmp(:);
temp2 = discretize(flw.ampByResp.avgRespAmp(:),p.resp.ampBins_edges);
flw.ampByResp.avgPosAmp_by_avgRespAmpBin = nan(length(p.resp.ampBins_x),1);
for j=1:length(p.resp.ampBins_x)
    if length(find(temp2==j)) >= p.resp.ampBins_minNumDataPoints
        flw.ampByResp.avgPosAmp_by_avgRespAmpBin(j) = nanmean(temp1(find(temp2==j)));
    end
end

temp1 = flw.ampByResp.avgNeuAmp(:);
temp2 = discretize(flw.ampByResp.avgRespAmp(:),p.resp.ampBins_edges);
flw.ampByResp.avgNeuAmp_by_avgRespAmpBin = nan(length(p.resp.ampBins_x),1);
for j=1:length(p.resp.ampBins_x)
    if length(find(temp2==j)) >= p.resp.ampBins_minNumDataPoints
        flw.ampByResp.avgNeuAmp_by_avgRespAmpBin(j) = nanmean(temp1(find(temp2==j)));
    end
end

temp1 = flw.ampByResp.avgNegAmp(:);
temp2 = discretize(flw.ampByResp.avgRespAmp(:),p.resp.ampBins_edges);
flw.ampByResp.avgNegAmp_by_avgRespAmpBin = nan(length(p.resp.ampBins_x),1);
for j=1:length(p.resp.ampBins_x)
    if length(find(temp2==j)) >= p.resp.ampBins_minNumDataPoints
        flw.ampByResp.avgNegAmp_by_avgRespAmpBin(j) = nanmean(temp1(find(temp2==j)));
    end
end

temp1 = flw.ampByResp.avgPosAmp(:);
temp2 = discretize(flw.ampByResp.avgRespTrgAmp(:),p.resp.ampBins_edges);
flw.ampByResp.avgPosAmp_by_avgRespTrgAmpBin = nan(length(p.resp.ampBins_x),1);
for j=1:length(p.resp.ampBins_x)
    if length(find(temp2==j)) >= p.resp.ampBins_minNumDataPoints
        flw.ampByResp.avgPosAmp_by_avgRespTrgAmpBin(j) = nanmean(temp1(find(temp2==j)));
    end
end

temp1 = flw.ampByResp.avgNeuAmp(:);
temp2 = discretize(flw.ampByResp.avgRespTrgAmp(:),p.resp.ampBins_edges);
flw.ampByResp.avgNeuAmp_by_avgRespTrgAmpBin = nan(length(p.resp.ampBins_x),1);
for j=1:length(p.resp.ampBins_x)
    if length(find(temp2==j)) >= p.resp.ampBins_minNumDataPoints
        flw.ampByResp.avgNeuAmp_by_avgRespTrgAmpBin(j) = nanmean(temp1(find(temp2==j)));
    end
end

temp1 = flw.ampByResp.avgNegAmp(:);
temp2 = discretize(flw.ampByResp.avgRespTrgAmp(:),p.resp.ampBins_edges);
flw.ampByResp.avgNegAmp_by_avgRespTrgAmpBin = nan(length(p.resp.ampBins_x),1);
for j=1:length(p.resp.ampBins_x)
    if length(find(temp2==j)) >= p.resp.ampBins_minNumDataPoints
        flw.ampByResp.avgNegAmp_by_avgRespTrgAmpBin(j) = nanmean(temp1(find(temp2==j)));
    end
end

% direct correlations
[flw.corr_ampByResp.avgPosAmp_by_avgTrgAmp_rho,flw.corr_ampByResp.avgPosAmp_by_avgTrgAmp_p] = corr(flw.ampByResp.avgPosAmp(:),flw.ampByResp.avgTrgAmp(:),'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNeuAmp_by_avgTrgAmp_rho,flw.corr_ampByResp.avgNeuAmp_by_avgTrgAmp_p] = corr(flw.ampByResp.avgNeuAmp(:),flw.ampByResp.avgTrgAmp(:),'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNegAmp_by_avgTrgAmp_rho,flw.corr_ampByResp.avgNegAmp_by_avgTrgAmp_p] = corr(flw.ampByResp.avgNegAmp(:),flw.ampByResp.avgTrgAmp(:),'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgPosAmp_by_avgRespAmp_rho,flw.corr_ampByResp.avgPosAmp_by_avgRespAmp_p] = corr(flw.ampByResp.avgPosAmp(:),flw.ampByResp.avgRespAmp(:),'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNeuAmp_by_avgRespAmp_rho,flw.corr_ampByResp.avgNeuAmp_by_avgRespAmp_p] = corr(flw.ampByResp.avgNeuAmp(:),flw.ampByResp.avgRespAmp(:),'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNegAmp_by_avgRespAmp_rho,flw.corr_ampByResp.avgNegAmp_by_avgRespAmp_p] = corr(flw.ampByResp.avgNegAmp(:),flw.ampByResp.avgRespAmp(:),'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgPosAmp_by_avgRespTrgAmp_rho,flw.corr_ampByResp.avgPosAmp_by_avgRespTrgAmp_p] = corr(flw.ampByResp.avgPosAmp(:),flw.ampByResp.avgRespTrgAmp(:),'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNeuAmp_by_avgRespTrgAmp_rho,flw.corr_ampByResp.avgNeuAmp_by_avgRespTrgAmp_p] = corr(flw.ampByResp.avgNeuAmp(:),flw.ampByResp.avgRespTrgAmp(:),'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNegAmp_by_avgRespTrgAmp_rho,flw.corr_ampByResp.avgNegAmp_by_avgRespTrgAmp_p] = corr(flw.ampByResp.avgNegAmp(:),flw.ampByResp.avgRespTrgAmp(:),'Type','Pearson','Rows','Complete');

% correlations with binned responsiveness metric
[flw.corr_ampByResp.avgPosAmp_by_avgTrgAmpBin_rho,flw.corr_ampByResp.avgPosAmp_by_avgTrgAmpBin_p] = corr(p.resp.ampBins_x,flw.ampByResp.avgPosAmp_by_avgTrgAmpBin,'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNeuAmp_by_avgTrgAmpBin_rho,flw.corr_ampByResp.avgNeuAmp_by_avgTrgAmpBin_p] = corr(p.resp.ampBins_x,flw.ampByResp.avgNeuAmp_by_avgTrgAmpBin,'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNegAmp_by_avgTrgAmpBin_rho,flw.corr_ampByResp.avgNegAmp_by_avgTrgAmpBin_p] = corr(p.resp.ampBins_x,flw.ampByResp.avgNegAmp_by_avgTrgAmpBin,'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgPosAmp_by_avgRespAmpBin_rho,flw.corr_ampByResp.avgPosAmp_by_avgRespAmpBin_p] = corr(p.resp.ampBins_x,flw.ampByResp.avgPosAmp_by_avgRespAmpBin,'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNeuAmp_by_avgRespAmpBin_rho,flw.corr_ampByResp.avgNeuAmp_by_avgRespAmpBin_p] = corr(p.resp.ampBins_x,flw.ampByResp.avgNeuAmp_by_avgRespAmpBin,'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNegAmp_by_avgRespAmpBin_rho,flw.corr_ampByResp.avgNegAmp_by_avgRespAmpBin_p] = corr(p.resp.ampBins_x,flw.ampByResp.avgNegAmp_by_avgRespAmpBin,'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgPosAmp_by_avgRespTrgAmpBin_rho,flw.corr_ampByResp.avgPosAmp_by_avgRespTrgAmpBin_p] = corr(p.resp.ampBins_x,flw.ampByResp.avgPosAmp_by_avgRespTrgAmpBin,'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNeuAmp_by_avgRespTrgAmpBin_rho,flw.corr_ampByResp.avgNeuAmp_by_avgRespTrgAmpBin_p] = corr(p.resp.ampBins_x,flw.ampByResp.avgNeuAmp_by_avgRespTrgAmpBin,'Type','Pearson','Rows','Complete');
[flw.corr_ampByResp.avgNegAmp_by_avgRespTrgAmpBin_rho,flw.corr_ampByResp.avgNegAmp_by_avgRespTrgAmpBin_p] = corr(p.resp.ampBins_x,flw.ampByResp.avgNegAmp_by_avgRespTrgAmpBin,'Type','Pearson','Rows','Complete');


%% Make figure

% follower analysis figure
F = followerAnalysisFigure(s2p_meta,iscell,trg,resp,flw,p,info);
savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','followerAnalysis.fig']);
saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','followerAnalysis.png']);
disp(['--- Saved follower analysis figures to ',path.filepart_outX,'plots.'])


%% Make custom STA
% 
% % n = 8; % neuron
% % i = 6; % cluster
% 
% act_woArtefact = act;
% %act_woArtefact(:,thor_beh.artefact) = NaN;
% 
% for n=1:numRois
%     for i=1:numClusters
%         try
%             if flw.stats_sig_neg(n,i)
% 
%                 frames_pre = 2*30;
%                 frames_post = 2*30;
% 
%                 % identify trials
%                 temp = trg.sequenceClusters(:,trg.sequenceOrder)==i;
%                 these_stimTrials = find(sum(temp));
%                 these_includedIpsiStimTrials = find(~squeeze(flw.inExclusionZone(n,i,:)));
%                 these_trials = these_stimTrials(these_includedIpsiStimTrials);
% 
%                 % identify cluster positions
%                 this_seqPosition = nan(length(these_trials),1);
%                 for k=1:length(these_trials)
%                     this_seqPosition(k) = find(temp(:,these_trials(k)));
%                 end
% 
%                 % make sta
%                 frames = -frames_pre+1:frames_post;
%                 sta = nan(length(these_trials),length(frames));
%                 for k=1:length(these_trials)
%                     this_anchor = thor_beh.stimPatternOnset(this_seqPosition(k),these_trials(k));
%                     sta(k,:) = act_woArtefact(n,this_anchor-frames_pre+1:this_anchor+frames_post);
%                 end
% 
%                 default_figure(); hold on;
%                 for k=1:length(these_trials)
%                     plot(frames,sta(k,:),'Color',p.col.gray)
%                 end
%                 shadedErrorBar(frames,nanmean(sta,1),nansem(sta,1),'lineProps',p.col.photostim)
%                 xline(0,'r-');
%                 xline(-1,'b-');
%                 xline(-4,'b-');
%                 xline(4,'b-');
%                 xline(7,'b-');
%                 %ylim([-1,2])
%             end
%         catch
%         end
%     end
% end



%% --- Old follower analysis (same as response analysis, but with absolute instead of relative baseline) ---

% disp('--- Follower analysis.')
% 
% % store used curation and parameters
% flw.iscell_used = iscell;
% flw.p = p;
% flw.p.maxTargetCentroidDistance_pix = p.resp.maxTargetCentroidDistance_um*info.scope.fovSize_pix/info.scope.fovSize_um;
% 
% 
% %% Get activity during analysis windows
% 
% disp('--- Extracting activity in analysis windows.')
% 
% % identify baseline and response window frames
% temp = find(task.var==1);
% flw.frames_base = nan(numPatternsPerSequence,numStimSequences,p.resp.numFrames);
% flw.frames_resp = nan(numPatternsPerSequence,numStimSequences,p.resp.numFrames);
% for i=1:numPatternsPerSequence
%     for j=1:numStimSequences
%         flw.frames_base(i,j,:) = thor_beh.sync(temp(j))-p.resp.numFrames:thor_beh.sync(temp(j))-1;
%         flw.frames_resp(i,j,:) = thor_beh.stimPatternComplete(i,j)+1:thor_beh.stimPatternComplete(i,j)+p.resp.numFrames;
%     end
% end
% 
% % get mean activity in analysis windows
% act_base = nan(numRois,numPatternsPerSequence,numStimSequences);
% act_resp = nan(numRois,numPatternsPerSequence,numStimSequences);
% for i=1:numPatternsPerSequence
%     for j=1:numStimSequences
%         try
%             act_base(:,i,j) = nanmean(act(:,squeeze(flw.frames_base(i,j,:))),2);
%             act_resp(:,i,j) = nanmean(act(:,squeeze(flw.frames_resp(i,j,:))),2);
%         catch
%             act_base(:,i,j) = NaN;
%             act_resp(:,i,j) = NaN;      
%         end
%     end
% end
% 
% % separate and sort pattern responses by cluster
% flw.act_base = nan(numRois,numClusters,numTrialsPerGroup);
% flw.act_resp = nan(numRois,numClusters,numTrialsPerGroup);
% flw.act_net = nan(numRois,numClusters,numTrialsPerGroup);
% stimTrials = nan(1,2*numTrialsPerGroup);
% stimTrials(1:length(task.type(task.var==1))) = task.type(task.var==1);
% this_stim1TrialCount = 0;
% this_stim2TrialCount = 0;
% for j=1:numStimSequences
%     if ~isnan(stimTrials(j))
% 
%         this_sequence = trg.sequenceOrder(j);
%         if (stimTrials(j)==1 | stimTrials(j)==3)
%             this_stim1TrialCount = this_stim1TrialCount+1;
%             this_trialIdx = this_stim1TrialCount;
%         elseif (stimTrials(j)==2 | stimTrials(j)==4)
%             this_stim2TrialCount = this_stim2TrialCount+1;
%             this_trialIdx = this_stim2TrialCount;
%         end
% 
%         these_clusters = trg.sequenceClusters(:,this_sequence);
%         for i=1:numPatternsPerSequence
%             flw.act_base(:,these_clusters(i),this_trialIdx) = act_base(:,i,j);
%             flw.act_resp(:,these_clusters(i),this_trialIdx) = act_resp(:,i,j);
%             flw.act_net(:,these_clusters(i),this_trialIdx) = flw.act_resp(:,these_clusters(i),this_trialIdx) - flw.act_base(:,these_clusters(i),this_trialIdx);
%         end
%     end
% end
% 
% % calculate trial-averaged activities
% flw.avgAct_base = nanmean(flw.act_base,3);
% flw.avgAct_resp = nanmean(flw.act_resp,3);
% flw.avgAct_net = nanmean(flw.act_net,3);
% 
% 
% %% Statstical tests
% 
% disp('--- Identifying responders.')
% 
% % calculate p-values
% flw.stats_p = nan(numRois,numClusters);
% for n=1:numRois
%     if iscell(n)
%         for i=1:numClusters
%             this_act_base = squeeze(flw.act_base(n,i,:));
%             this_act_resp = squeeze(flw.act_resp(n,i,:));
%             if sum(~isnan(this_act_base)) && sum(~isnan(this_act_resp))
%                 flw.stats_p(n,i) = signrank(this_act_base,this_act_resp);
%             end
%         end
%     end
% end
% if ~isempty(find(flw.stats_p==0))
%     warning('Something went wrong with p-value calculation.')
% end
% 
% % calculate fdr and q-values
% flw.stats_fdr = nan(numRois,numClusters);
% flw.stats_q = nan(numRois,numClusters);
% flw.stats_priori = nan(1,numClusters);
% for i=1:numClusters
%     [flw.stats_fdr(:,i),flw.stats_q(:,i),flw.stats_priori(i)] = mafdr(flw.stats_p(:,i)); % double-checked that it works correctly with NaNs
% end
% if any(flw.stats_priori>=1)
%     warning('The priori for at least one cluster is greater than or equal to 1.')
% end
% 
% % find significant followers/responders
% flw.stats_sig = double(flw.stats_q<p.resp.qthres);
% flw.stats_sig(isnan(flw.stats_q))=NaN;
% flw.stats_sig_pos = double(double(flw.avgAct_net>0)+flw.stats_sig==2);
% flw.stats_sig_pos(isnan(flw.stats_q))=NaN;
% flw.stats_sig_neg = double(double(flw.avgAct_net<0)+flw.stats_sig==2);
% flw.stats_sig_neg(isnan(flw.stats_q))=NaN;
% 
% 
% %% Apply criteria
% 
% for i=1:length(p.resp.amplitudeCritiera)
%     flw = applyResponderCriteria(trg,flw,s2p_meta,p.resp.amplitudeCritiera(i),0);
%     flw = applyResponderCriteria(trg,flw,s2p_meta,p.resp.amplitudeCritiera(i),1);
% end
% flw.responders_main = flw.(p.resp.mainResponderCondition);
% 
% 
% %% Photostimulation stability
% 
% n = 0;
% flw.trgRespAmps_cells = nan(size(flw.responders_main.targetedCells,1)*size(flw.responders_main.targetedCells,2),size(flw.act_net,3));
% flw.trgRespAmps_clusters = nan(size(flw.responders_main.targetedCells,2),size(flw.act_net,3));
% for j=1:size(flw.responders_main.targetedCells,2)
%     for i=1:size(flw.responders_main.targetedCells,1)
%         n=n+1;
%         if ~isnan(flw.responders_main.targetedCells(i,j))
%             flw.trgRespAmps_cells(n,:) = squeeze(flw.act_net(flw.responders_main.targetedCells(i,j),j,:));
%         end
%     end
%     flw.trgRespAmps_clusters(j,:) = nanmean(flw.trgRespAmps_cells(n-size(flw.responders_main.targetedCells,1)+1:n,:),1);
% end
% temp = reshape(flw.trgRespAmps_cells,[size(flw.trgRespAmps_cells,1),((info.task.trialsPerBlock/2)*(info.task.numStimTrials/info.task.numTrials)),info.task.numBlocks]);
% flw.trgRespAmps_cells_bw = squeeze(nanmean(temp,2));
% temp = reshape(flw.trgRespAmps_clusters,[size(flw.trgRespAmps_clusters,1),((info.task.trialsPerBlock/2)*(info.task.numStimTrials/info.task.numTrials)),info.task.numBlocks]);
% flw.trgRespAmps_clusters_bw = squeeze(nanmean(temp,2));
% 
% 
% %% Save
% 
% flw = orderfields(flw);
% save([path.filepart_out,'flw.mat'],'flw','-v7.3');
% disp(['--- Saved flw file as ',[path.filepart_out,'flw.mat'],'.'])
% 
% 
% %% Make figures
% 
% % cluster responses maps
% if ~ops.resp.skip_clusterResponseMaps
%     if ~exist([path.filepart_outX,'plots\ClusterResponseMaps'],'dir')
%         mkdir([path.filepart_outX,'plots\ClusterResponseMaps']);
%     end
%     for i=1:size(flw.responders_main.targetedCells,2)
%         F = clusterResponseMap(s2p_meta,iscell,trg,flw,i,0); 
%         savefig(F,[path.filepart_outX,'plots\ClusterResponseMaps\',info.animal,'_',info.date,'_','cluster_',num2str(i),'_flw.fig']);
%         saveas(F,[path.filepart_outX,'plots\ClusterResponseMaps\',info.animal,'_',info.date,'_','cluster_',num2str(i),'_flw.png']);
%     end
% end
% 
% % response analysis figure
% F = responseAnalysisFigure(s2p_meta,iscell,trg,flw,p,info);
% % F = responseAnalysisFigure_reduced(s2p_meta,iscell,trg,resp,p,info);
% savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','responseAnalysis_flw.fig']);
% saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','responseAnalysis_flw.png']);
% disp(['--- Saved response analysis (flw) figures to ',path.filepart_outX,'plots.'])


%% Return

if ops.close_figures
    close all;
end
end


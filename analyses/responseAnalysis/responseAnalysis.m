function [resp] = responseAnalysis(info,iscell,ops,p,path,trg,task,s2p_meta,thor_beh,act)
% act = dFFn_beh; p = get_p;

%% Preparations

resp.stimType = info.stimType;

% pre-process activity measure
if p.resp.zscore
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
try
    resp.pos_closestLaser_x = trg.pos_closestLaser_x;
    resp.pos_closestLaser_y = trg.pos_closestLaser_y;
catch
end


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

% calculate trial-averaged activities
resp.avgAct_base = nanmean(resp.act_base,3);
resp.avgAct_resp = nanmean(resp.act_resp,3);
resp.avgAct_net = nanmean(resp.act_net,3);


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
if p.resp.excludeOutsideExclusionZone
    for i=1:numClusters
        resp.stats_q(trg.dist_closestLaser(:,i) > p.resp.exclusionRadius,i) = NaN; % exclude neurons above 50 um radius
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


%% Apply criteria

for i=1:length(p.resp.amplitudeCritiera)
    resp = applyResponderCriteria(trg,resp,s2p_meta,p.resp.amplitudeCritiera(i),0);
    resp = applyResponderCriteria(trg,resp,s2p_meta,p.resp.amplitudeCritiera(i),1);
end
resp.responders_main = resp.(p.resp.mainResponderCondition);


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

n = 0;
resp.respAmps_cells = nan(size(resp.responders_main.respAll,1),size(resp.act_net,3));
for j=1:size(resp.responders_main.respAll,2)
    for i=1:size(resp.responders_main.respAll,1)
        if resp.responders_main.respAll(i,j)==1
            n=n+1;
            resp.respAmps_cells(n,:) = squeeze(resp.act_net(i,j,:));
        end
    end
end
resp.respAmps_cells = rmmissing(resp.respAmps_cells,1,'MinNumMissing',size(resp.respAmps_cells,2));
temp = reshape(resp.respAmps_cells,[size(resp.respAmps_cells,1),((info.task.trialsPerBlock/2)*(info.task.numStimTrials/info.task.numTrials)),info.task.numBlocks]);
resp.respAmps_cells_bw = squeeze(nanmean(temp,2));


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
        % fishing for good cluster response maps (has to be a first seq stim cluster, ideally of a runner)
        % - Jobs_20220721: 36 nope, 31 nope
        % - Austin_20220729: 12 maybe but bad resolution and too much inhibition, 16 apart from too much inhibition its nice
        % - Philip_20211004: 8 is nice!, 34 is nice as well!
        % - Celo_20210706: 10 okay I guess, 20 nah
        % - Turing_20220916: 38, 7
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


%% Follower analysis

if ops.do_followerAnalysis

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
                elseif false % subtracting average baseline from baseline and response
                    flw.act_base(n,i,:) = resp.act_base(n,i,:) - resp.avgAct_base(n,i);
                    flw.act_resp(n,i,:) = resp.act_resp(n,i,:) - resp.avgAct_base(n,i);
                    flw.act_net(n,i,:) = resp.act_net(n,i,:) - resp.avgAct_base(n,i);
                else
                    flw.act_base(n,i,:) = resp.act_base(n,i,:);
                    flw.act_resp(n,i,:) = resp.act_resp(n,i,:);
                    flw.act_net(n,i,:) = resp.act_net(n,i,:);
                end

                % trial exclusions
                if true
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
    
    % find significant followers/responders - negative followers
    flw.stats_q = flw.stats_q_raw;
    flw.stats_sig = double(flw.stats_q<p.resp.qthres);
    flw.stats_sig(isnan(flw.stats_q))=NaN;
    flw.stats_sig_neg = double(double(flw.avgAct_net<0)+flw.stats_sig==2);
    flw.stats_sig_neg(isnan(flw.stats_q))=NaN;
    
    % filter out 
    for i=1:numClusters
        flw.stats_q(trg.dist_closestLaser(:,i)< p.resp.exclusionRadius,i) = NaN; % exclude neurons within 50 um radius
%         %flw.stats_q(resp.responders_main.idcs_resp{i},i) = NaN; % exclude responders
%         flw.stats_q(resp.responders_main.idcs_targeted{i},i) = NaN; % exclude targeted neurons
    end

    % find significant followers/responders - positive followers
    flw.stats_sig = double(flw.stats_q<p.resp.qthres);
    flw.stats_sig(isnan(flw.stats_q))=NaN;
    flw.stats_sig_pos = double(double(flw.avgAct_net>0)+flw.stats_sig==2);
    flw.stats_sig_pos(isnan(flw.stats_q))=NaN;
    
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


    %% Photostimulation amplitude and stability

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
    
    
    %% Spatial plot


    
    
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

end


%% Make follower STAs

if ops.do_followerSTAs
    
    close all;

    act_woArtefact = act;
    %act_woArtefact(:,thor_beh.artefact) = NaN;

    frames_pre = length(p.general.frames_pre); %60
    frames_post = -length(p.general.frames_pre)+length(p.general.t_unbinned); %60
    flwsta.frames = -frames_pre+1:frames_post;
    flwsta.frame_stim = 0;
    flwsta.frame_preOffset = -1;
    flwsta.frame_preOnset = -4;
    flwsta.frame_postOnset = 4;
    flwsta.frame_postOffset = 7;

    % negative followers
    flwsta.neg_roi = [];
    flwsta.neg_clust = [];
    flwsta.neg_sta = [];
    for n=1:numRois
        for i=1:numClusters
            try
                if flw.stats_sig_neg(n,i)


                    % identify trials
                    temp = trg.sequenceClusters(:,trg.sequenceOrder)==i;
                    these_stimTrials = find(sum(temp));
                    these_includedIpsiStimTrials = find(~squeeze(flw.inExclusionZone(n,i,:)));
                    these_trials = these_stimTrials(these_includedIpsiStimTrials);

                    % identify cluster positions
                    this_seqPosition = nan(length(these_trials),1);
                    for k=1:length(these_trials)
                        this_seqPosition(k) = find(temp(:,these_trials(k)));
                    end

                    % make sta
                    sta = nan(length(these_trials),length(flwsta.frames));
                    for k=1:length(these_trials)
                        this_anchor = thor_beh.stimPatternOnset(this_seqPosition(k),these_trials(k));
                        sta(k,:) = act_woArtefact(n,this_anchor-frames_pre+1:this_anchor+frames_post);
                    end

    %                 figure(); hold on;
    %                 for k=1:length(these_trials)
    %                     plot(flwsta.frames,sta(k,:),'Color',p.col.gray)
    %                 end
    %                 shadedErrorBar(flwsta.frames,nanmean(sta,1),nansem(sta,1),'lineProps',p.col.photostim)
    %                 xline(flwsta.frame_stim,'r-');
    %                 xline(flwsta.frame_preOffset,'b-');
    %                 xline(flwsta.frame_preOnset,'b-');
    %                 xline(flwsta.frame_postOnset,'b-');
    %                 xline(flwsta.frame_postOffset,'b-');
    %                 title('negative follower')
    %                 drawnow;

                    flwsta.neg_roi = [flwsta.neg_roi;n];
                    flwsta.neg_clust = [flwsta.neg_clust;i];
                    flwsta.neg_sta = [flwsta.neg_sta;nanmean(sta,1)];
                end
            catch
            end
        end
    end

    % positive followers
    flwsta.pos_roi = [];
    flwsta.pos_clust = [];
    flwsta.pos_sta = [];
    for n=1:numRois
        for i=1:numClusters
            try
                if flw.stats_sig_pos(n,i)


                    % identify trials
                    temp = trg.sequenceClusters(:,trg.sequenceOrder)==i;
                    these_stimTrials = find(sum(temp));
                    these_includedIpsiStimTrials = find(~squeeze(flw.inExclusionZone(n,i,:)));
                    these_trials = these_stimTrials(these_includedIpsiStimTrials);

                    % identify cluster positions
                    this_seqPosition = nan(length(these_trials),1);
                    for k=1:length(these_trials)
                        this_seqPosition(k) = find(temp(:,these_trials(k)));
                    end

                    % make sta
                    sta = nan(length(these_trials),length(flwsta.frames));
                    for k=1:length(these_trials)
                        this_anchor = thor_beh.stimPatternOnset(this_seqPosition(k),these_trials(k));
                        sta(k,:) = act_woArtefact(n,this_anchor-frames_pre+1:this_anchor+frames_post);
                    end

    %                 figure(); hold on;
    %                 for k=1:length(these_trials)
    %                     plot(flwsta.frames,sta(k,:),'Color',p.col.gray)
    %                 end
    %                 shadedErrorBar(flwsta.frames,nanmean(sta,1),nansem(sta,1),'lineProps',p.col.photostim)
    %                 xline(flwsta.frame_stim,'r-');
    %                 xline(flwsta.frame_preOffset,'b-');
    %                 xline(flwsta.frame_preOnset,'b-');
    %                 xline(flwsta.frame_postOnset,'b-');
    %                 xline(flwsta.frame_postOffset,'b-');
    %                 title('positive follower')
    %                 drawnow;

                    flwsta.pos_roi = [flwsta.pos_roi;n];
                    flwsta.pos_clust = [flwsta.pos_clust;i];
                    flwsta.pos_sta = [flwsta.pos_sta;nanmean(sta,1)];
                end
            catch
            end
        end
    end

    % responders
    flwsta.resp_roi = [];
    flwsta.resp_clust = [];
    flwsta.resp_sta = [];
    for n=1:numRois
        for i=1:numClusters
            try
                if resp.responders_main.respAll  (n,i)


                    % identify trials
                    temp = trg.sequenceClusters(:,trg.sequenceOrder)==i;
                    these_stimTrials = find(sum(temp));
                    these_includedIpsiStimTrials = find(~squeeze(flw.inExclusionZone(n,i,:)));
                    these_trials = these_stimTrials(these_includedIpsiStimTrials);

                    % identify cluster positions
                    this_seqPosition = nan(length(these_trials),1);
                    for k=1:length(these_trials)
                        this_seqPosition(k) = find(temp(:,these_trials(k)));
                    end

                    % make sta
                    sta = nan(length(these_trials),length(flwsta.frames));
                    for k=1:length(these_trials)
                        this_anchor = thor_beh.stimPatternOnset(this_seqPosition(k),these_trials(k));
                        sta(k,:) = act_woArtefact(n,this_anchor-frames_pre+1:this_anchor+frames_post);
                    end

    %                 figure(); hold on;
    %                 for k=1:length(these_trials)
    %                     plot(flwsta.frames,sta(k,:),'Color',p.col.gray)
    %                 end
    %                 shadedErrorBar(flwsta.frames,nanmean(sta,1),nansem(sta,1),'lineProps',p.col.photostim)
    %                 xline(flwsta.frame_stim,'r-');
    %                 xline(flwsta.frame_preOffset,'b-');
    %                 xline(flwsta.frame_preOnset,'b-');
    %                 xline(flwsta.frame_postOnset,'b-');
    %                 xline(flwsta.frame_postOffset,'b-');
    %                 title('responder')
    %                 drawnow;

                    flwsta.resp_roi = [flwsta.resp_roi;n];
                    flwsta.resp_clust = [flwsta.resp_clust;i];
                    flwsta.resp_sta = [flwsta.resp_sta;nanmean(sta,1)];
                end
            catch
            end
        end
    end

    % % neutral followers (non-follower)
    % flwsta.neu_roi = [];
    % flwsta.neu_clust = [];
    % flwsta.neu_sta = [];
    % for n=1:numRois
    %     for i=1:numClusters
    %         try
    %             if ~flw.stats_sig_neg(n,i)
    % 
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
    %                 sta = nan(length(these_trials),length(flwsta.frames));
    %                 for k=1:length(these_trials)
    %                     this_anchor = thor_beh.stimPatternOnset(this_seqPosition(k),these_trials(k));
    %                     sta(k,:) = act_woArtefact(n,this_anchor-frames_pre+1:this_anchor+frames_post);
    %                 end
    % 
    % %                 figure(); hold on;
    % %                 for k=1:length(these_trials)
    % %                     plot(flwsta.frames,sta(k,:),'Color',p.col.gray)
    % %                 end
    % %                 shadedErrorBar(flwsta.frames,nanmean(sta,1),nansem(sta,1),'lineProps',p.col.photostim)
    % %                 xline(flwsta.frame_stim,'r-');
    % %                 xline(flwsta.frame_preOffset,'b-');
    % %                 xline(flwsta.frame_preOnset,'b-');
    % %                 xline(flwsta.frame_postOnset,'b-');
    % %                 xline(flwsta.frame_postOffset,'b-');
    % %                 drawnow;
    %                 
    %                 flwsta.neu_roi = [flwsta.neu_roi;n];
    %                 flwsta.neu_clust = [flwsta.neu_clust;i];
    %                 flwsta.neu_sta = [flwsta.neu_sta;nanmean(sta,1)];
    %             end
    %         catch
    %         end
    %     end
    % end


    %%% Save

    flwsta = orderfields(flwsta);
    save([path.filepart_out,'flwsta.mat'],'flwsta','-v7.3');
    disp(['--- Saved flwsta file as ',[path.filepart_out,'flwsta.mat'],'.'])
end


%% Return

if ops.close_figures
    close all;
end
end


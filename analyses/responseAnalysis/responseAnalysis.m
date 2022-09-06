function [resp] = responseAnalysis(info,iscell,ops,p,path,trg,task,s2p_meta,thor_beh,act)
% act = dFFn_beh; p = get_p;

%% Preparations

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
            act_base(:,i,j) = nanmean(act(:,squeeze(resp.frames_base(i,j,:))),2);
            act_resp(:,i,j) = nanmean(act(:,squeeze(resp.frames_resp(i,j,:))),2);
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
        elseif (stimTrials(j)==2 | stimTrials(j)==4)
            this_stim2TrialCount = this_stim2TrialCount+1;
            this_trialIdx = this_stim2TrialCount;
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

% find significant followers/responders
resp.stats_sig = double(resp.stats_q<p.resp.qthres);
resp.stats_sig(isnan(resp.stats_q))=NaN;
resp.stats_sig_pos = double(double(resp.avgAct_net>0)+resp.stats_sig==2);
resp.stats_sig_pos(isnan(resp.stats_q))=NaN;
resp.stats_sig_neg = double(double(resp.avgAct_net<0)+resp.stats_sig==2);
resp.stats_sig_neg(isnan(resp.stats_q))=NaN;


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

if ops.close_figures
    close all;
end
end


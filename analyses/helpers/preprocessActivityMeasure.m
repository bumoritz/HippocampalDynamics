function [prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,this_p,p,paq_beh,iscell)
% act = spks_beh; this_p = this_p;

%% Preparations

% get properties of raw data
prop.iscell = iscell;
prop.numCells = sum(prop.iscell);
prop.numRois = size(act,1);
prop.numFramesTotal = size(act,2);

% get properties of binned data
prop.meanFrames2bins = floor(p.general.binSize/2)+1:p.general.binSize:prop.numFramesTotal;
prop.numFramesTotal_binned = length(prop.meanFrames2bins) - 1;

% get properties of trial-wise data
if isfield(paq_beh,'sync_all')
    prop.sync = paq_beh.sync;
    prop.sync_all = paq_beh.sync_all;
else
    prop.sync = paq_beh.sync;
    prop.sync_all = paq_beh.sync;
end
prop.numTrials = length(prop.sync_all);

% bin sync data
temp = zeros(1,prop.numFramesTotal);
temp(prop.sync_all) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
prop.sync_all_binned = temp(:,prop.meanFrames2bins);
prop.sync_all_binned = find(prop.sync_all_binned(:,1:end-1));
if length(prop.sync_all_binned)~=prop.numTrials
    warning('Number of SniffinSyncs not equal to number of trials')
end

% identify trial bounds
these_trial_bounds(1,:) = prop.sync_all_binned-length(p.general.bins_pre);
these_trial_bounds(2,:) = these_trial_bounds(1,:)+p.general.numBins-1;
for i=1:length(prop.sync_all_binned)
    prop.trial_frames_binned(:,i) = these_trial_bounds(1,i) : these_trial_bounds(2,i);
end

prop = orderfields(prop);


%% Pre-process activity measure

% deal with problematic sessions (NaN trials, too short imaging)
% if isfield(paq_beh,'sync_all')
%     act(:,prop.sync_all(nanmin(find(isnan(prop.sync))))-2*this_p.frames_pre:end) = NaN;
% end
% if nanmax(prop.trial_frames(:))>size(act,2)
%     temp = act;
%     act = NaN(size(act,1),nanmax(prop.trial_frames(:)));
%     act(:,1:size(temp,2)) = temp;
% end

% set non-iscells to NaN
temp = act;
act = nan(size(act));
act(find(prop.iscell),:) = temp(find(prop.iscell),:);

% 1) smoothing (pre-binning)
if this_p.smoothingSd_preBinning~=0
	act = smoothdata(act,2,'gaussian',this_p.smoothingSd_preBinning*5);
end

% 2) binning
temp = movmean(act,p.general.binSize,2,'omitnan');
act = temp(:,floor(p.general.binSize/2)+1:p.general.binSize:end);
act = act(:,1:end-1);

% 3) smoothing (post-binning)
if this_p.smoothingSd_postBinning~=0
	act = smoothdata(act,2,'gaussian',this_p.smoothingSd_postBinning*5);
end

% 4) z-scoring
if this_p.zscore
     act = nanzscore(act,[],2);
end

% 5) split into trials
try
    nf_binned = act;
    nft_binned = reshape(nf_binned(:,prop.trial_frames_binned),[],size(prop.trial_frames_binned,1),prop.numTrials);
catch
    nft_binned = nan(size(nf_binned,1),size(prop.trial_frames_binned,1),prop.numTrials);
    try
        nf_binned = act;
        temp = reshape(nf_binned(:,prop.trial_frames_binned),[],size(prop.trial_frames_binned,1),size(prop.trial_frames_binned,2));
    catch
        prop.sync = prop.sync(1:end-1);
        prop.sync_all = prop.sync_all(1:end-1);
        prop.sync_all_binned = prop.sync_all_binned(1:end-1);
        prop.trial_frames_binned = prop.trial_frames_binned(:,1:end-1);
        nf_binned = act;
        temp = reshape(nf_binned(:,prop.trial_frames_binned),[],size(prop.trial_frames_binned,1),size(prop.trial_frames_binned,2));
    end
    nft_binned(:,:,1:size(temp,3)) = temp;
    warning('Dealt with the fact that number of SniffinSyncs not equal to number of trials')
end
prop.numTrials_incl = size(prop.trial_frames_binned,2);


%% Get mean and sd of traces

prop.nf_mean = nanmean(nf_binned,2);
prop.nf_median = nanmedian(nf_binned,2);
prop.nf_sd = nanstd(nf_binned,[],2);


end

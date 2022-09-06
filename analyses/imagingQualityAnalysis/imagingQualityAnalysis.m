%function [iqa] = imagingQualityAnalysis(info,iscell,ops,p,path,s2p_meta,thor_beh,act)
p=get_p; sync_beh = paq_beh;
act = dFF_beh; % spks_beh


%% Calculate dFF-based (unbinned, unsmoothed, but z-scored) noise measure

iqa.mean = nanmean(act,2);
iqa.std = nanstd(act,[],2);
iqa.snr = iqa.mean ./ iqa.std;
iqa.skew = skewness(act,[],2);

temp = nanzscore(act,[],2);
iqa.noise = nanmedian(abs(diff(temp,[],2)),2)/sqrt(info.scope.frameRate);


%%












%% Preparations


%%

absoff = hypot(s2p_meta.ops.xoff,s2p_meta.ops.yoff);
absdisp = abs(diff(absoff));  % use this to identify motion

absdisp = smoothdata(absdisp,2,'gaussian',p.tng.smoothingSd_preBinning*5);


%%

temp = movmean([NaN,absdisp]',p.general.binSize,2,'omitnan');
events_binned.absdisp = temp(floor(p.general.binSize/2)+1:p.general.binSize:end);
events_binned.absdisp = events_binned.absdisp(1:end-1);

absdisp_binned = nan(size(prop.trial_frames_binned));
for i=1:size(prop.trial_frames_binned,2)
    absdisp_binned(:,i) = events_binned.absdisp(prop.trial_frames_binned(:,i));
end

%%

figure;
plot(nanmean(absdisp_binned,2))

%%











%% Make figure

% imaging quality analysis figure
F = imagingQualityAnalysisFigure(info,p,iscell,s2p_meta,thor_beh);
savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','imagingQualityAnalysis.fig']);
saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','imagingQualityAnalysis.png']);
disp(['--- Saved imaging quality analysis figures to ',path.filepart_outX,'plots.'])

if ops.close_figures
    close all;
end
%end
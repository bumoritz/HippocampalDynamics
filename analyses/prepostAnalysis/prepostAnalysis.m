function [ppa] = prepostAnalysis(info,iscell,ops,p,path,act_pre,act_post)
% p = get_p;

%% Preparations

disp('--- Preparations')

% pre-process activity measure pre  
act = act_pre;
temp = act;
act = nan(size(act));
act(find(iscell),:) = temp(find(iscell),:);
act = smoothdata(act,2,'gaussian',p.ppa.smoothingSd_preBinning*5);
temp = movmean(act,p.general.binSize,2,'omitnan');
act = temp(:,floor(p.general.binSize/2)+1:p.general.binSize:end);
act = act(:,1:end-1);
nf_binned_pre_noZ = act;
nf_binned_pre_withinZ = nanzscore(act,[],2);

% pre-process activity measure post  
act = act_post;
temp = act;
act = nan(size(act));
act(find(iscell),:) = temp(find(iscell),:);
act = smoothdata(act,2,'gaussian',p.ppa.smoothingSd_preBinning*5);
temp = movmean(act,p.general.binSize,2,'omitnan');
act = temp(:,floor(p.general.binSize/2)+1:p.general.binSize:end);
act = act(:,1:end-1);
nf_binned_post_noZ = act;  
nf_binned_post_withinZ = nanzscore(act,[],2);

% z-score across pre and post
this_mean = nanmean([nf_binned_pre_noZ,nf_binned_post_noZ],2);
this_std = nanstd([nf_binned_pre_noZ,nf_binned_post_noZ],[],2);
nf_binned_pre_acrossZ = (nf_binned_pre_noZ - this_mean) ./ this_std;
nf_binned_post_acrossZ = (nf_binned_post_noZ - this_mean) ./ this_std;


%% Calculate correlation matrix

disp('--- Calculating correlations')

this_selection = find(~(sum(isnan(nf_binned_pre_withinZ),2)>100));
ppa.withinZ_pre = nan(length(iscell));
ppa.withinZ_post = nan(length(iscell));
ppa.acrossZ_pre = nan(length(iscell));
ppa.acrossZ_post = nan(length(iscell));
temp = corr(nf_binned_pre_withinZ(this_selection,:)',nf_binned_pre_withinZ(this_selection,:)','Type','Pearson','Rows','Complete');
ppa.withinZ_pre(this_selection,this_selection) = temp;
temp = corr(nf_binned_post_withinZ(this_selection,:)',nf_binned_post_withinZ(this_selection,:)','Type','Pearson','Rows','Complete');
ppa.withinZ_post(this_selection,this_selection) = temp;
temp = corr(nf_binned_pre_acrossZ(this_selection,:)',nf_binned_pre_acrossZ(this_selection,:)','Type','Pearson','Rows','Complete');
ppa.acrossZ_pre(this_selection,this_selection) = temp;
temp = corr(nf_binned_post_acrossZ(this_selection,:)',nf_binned_post_acrossZ(this_selection,:)','Type','Pearson','Rows','Complete');
ppa.acrossZ_post(this_selection,this_selection) = temp;


%% Save

ppa = orderfields(ppa);
save([path.filepart_out,'ppa.mat'],'ppa','-v7.3');
disp(['--- Saved ppa file as ',[path.filepart_out,'ppa.mat'],'.'])


%% Return

if ops.close_figures
    close all;
end
end







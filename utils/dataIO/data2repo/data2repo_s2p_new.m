function [path] = data2repo_s2p_new(info,p,path,ops,s2p_meta,paq_beh,thor_beh,raw)

%% Assign imaging frames to pre, beh and post epochs and sanity checks

try 
    s2pFrames_pre = s2p_meta.raw.pre.frames;
    s2pFrames_post = s2p_meta.raw.post.frames;
    s2pFrames_beh = s2p_meta.raw.beh.frames;
catch
    
    % pre
    s2pFrames_pre = raw.pre.frames;
    if raw.pre.frames > info.scope.numFrames_pre
        warning(['Mismatch in number of pre frames between imaging raw file and user input',newline...
            'imaging raw file: ', num2str(raw.pre.frames),', user input: ', num2str(info.scope.numFrames_pre)])
    elseif raw.pre.frames < info.scope.numFrames_pre
        warning(['Mismatch in number of pre frames between imaging raw file and user input',newline...
            'imaging raw file: ', num2str(raw.pre.frames),', user input: ', num2str(info.scope.numFrames_pre)])
    end

    % post
    s2pFrames_post = raw.post.frames;
    if raw.post.frames > info.scope.numFrames_post
        warning(['Mismatch in number of post frames between imaging raw file and user input',newline...
            'imaging raw file: ', num2str(raw.post.frames),', user input: ', num2str(info.scope.numFrames_post)])
    elseif raw.post.frames < info.scope.numFrames_post
        warning(['Mismatch in number of post frames between imaging raw file and user input',newline...
            'imaging raw file: ', num2str(raw.post.frames),', user input: ', num2str(info.scope.numFrames_post)])
    end

    % beh
    if info.data.numFragments~=1
        s2pFrames_beh = 0;
        for i=1:info.data.numFragments
            s2pFrames_beh = s2pFrames_beh + raw.(['beh_',num2str(i)]).frames;
        end
    else
        s2pFrames_beh = raw.beh.frames;
    end
end


%% Preparations

iscell = s2p_meta.iscell(:,1);
if isfield(paq_beh,'sync_all')
    sync_beh.sync = paq_beh.sync + s2pFrames_pre;
    sync_beh.sync_all = paq_beh.sync_all + s2pFrames_pre;
else
    sync_beh.sync = paq_beh.sync + s2pFrames_pre;
    sync_beh.sync_all = paq_beh.sync + s2pFrames_pre;
end


%% Load F (and save F_pre, F_beh, F_post)

disp('--- Loading fluorescence data...')

% load F_all data
path.file_in_s2p_F = [path.folder_data,'Imaging\suite2p\plane0\F.npy'];
F_all = readNPY(path.file_in_s2p_F);
[F_pre,F_beh,F_post] = splitPreBehPost(F_all,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);
save([path.filepart_outX,'F_pre.mat'],'F_pre','-v7.3');
disp(['--- Added F_pre file to repoX as ',[path.filepart_outX,'F_pre.mat'],'.'])
save([path.filepart_outX,'F_beh.mat'],'F_beh','-v7.3');
disp(['--- Added F_beh file to repoX as ',[path.filepart_outX,'F_beh.mat'],'.'])
save([path.filepart_outX,'F_post.mat'],'F_post','-v7.3');
disp(['--- Added F_post file to repoX as ',[path.filepart_outX,'F_post.mat'],'.'])


%% Do SVD and remove artefact components

disp('--- Doing SVD and removing artefact components...')

% do SVD
if p.s2p_new.includeOnlyIscells
    data_in = single(F_all(find(iscell==1),:));
else
    data_in = single(F_all);
end
M = data_in';
[U,S,V] = svdecon(M);
M_hat = U*S*V';
data_out = M_hat';

% calculate zero crossing rates of temporal singular vectors
if info.stimSession
    if p.s2p_new.holdOutStimData
        sync_beh.stimSequence = thor_beh.stimSequence - 9 + s2pFrames_pre;
        these_frames = 1:size(U,1);
        for i=1:length(sync_beh.stimSequence)
            these_frames = setdiff(these_frames,sync_beh.stimSequence(i):sync_beh.stimSequence(i)+8*round(info.scope.frameRate));
        end
        this_U = U(these_frames,:);
    else
        this_U = U;
    end
else
    this_U = U;
end
zcr = nansum(abs(diff(this_U>=0,[],1)),1);
zcrd = abs(diff(zcr));
zcrd_threshold = p.s2p_new.zcrdSigmaThreshold*nanstd(zcrd);
singularValuesToKeep = ones(1,length(zcr));
singularValuesToKeep([true,zcrd>zcrd_threshold]) = 0;

% reconstruct data without artefact components
S_zcrFiltered = S;
S_zcrFiltered(:,find(singularValuesToKeep==0)) = 0;
M_hat_zcrFiltered = U*S_zcrFiltered*V';
data_out_zcrFiltered = M_hat_zcrFiltered';

% reinsert non-iscells
if p.s2p_new.includeOnlyIscells
    Fn = nan(size(F_all));
    Fn(find(iscell==1),:) = data_out_zcrFiltered;
else
    Fn = data_out_zcrFiltered;
end


%% Prepare SVD visualisation

disp('--- Visualising SVD...')

% create trial-averages of SVD matrices
[prop,~,data_in_avg] = preprocessActivityMeasure(data_in,p.inh,p,sync_beh,ones(size(data_in,1),1));
data_in_avg = nanmean(data_in_avg,3);
[~,~,data_out_avg] = preprocessActivityMeasure(data_out,p.inh,p,sync_beh,ones(size(data_out,1),1));
data_out_avg = nanmean(data_out_avg,3);
[~,~,data_out_zcrFiltered_avg] = preprocessActivityMeasure(data_out_zcrFiltered,p.inh,p,sync_beh,ones(size(data_out_zcrFiltered,1),1));
data_out_zcrFiltered_avg = nanmean(data_out_zcrFiltered_avg,3);
[~,~,U_avg] = preprocessActivityMeasure(U',p.inh,p,sync_beh,ones(size(data_out_zcrFiltered_avg,1),1));
U_avg = nanmean(U_avg,3);

% create trial-averages of removed components
remU_avg = cell(1,sum(1-singularValuesToKeep));
temp = find(singularValuesToKeep==0);
for i=1:sum(1-singularValuesToKeep)
    temp0 = U(:,temp(i))';
    temp1 = movmean(temp0,p.general.binSize,2,'omitnan');
    temp2 = temp1(:,floor(p.general.binSize/2)+1:p.general.binSize:end);
    temp2 = temp2(:,1:end-1);
    temp3 = nan(1,size(prop.trial_frames_binned,1)*prop.numTrials);
    temp3(1:length(prop.trial_frames_binned(:))) = temp2(:,prop.trial_frames_binned);
    remU_avg{1,i} = nanmean(reshape(temp3,[],size(prop.trial_frames_binned,1),prop.numTrials),3);
end


%% SVD summary figure

nrows = 2; ncols = 5;
F = default_figure();

r=1; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_in_avg',[],1))
taskLines(p,info,'full','heatmap',false,0,1);
xlabel('ROI (n)')
ylabel('Bin (~m)')
title('M (z-scored)')

r=1; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_out_avg',[],1))
taskLines(p,info,'full','heatmap',false,0,1);
xlabel('ROI (n)')
ylabel('Bin (~m)')
title('M_{hat} original (z-scored)')

r=1; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_out_zcrFiltered_avg',[],1))
taskLines(p,info,'full','heatmap',false,0,1);
xlabel('ROI (n)')
ylabel('Bin (~m)')
title('M_{hat} decontaminated (z-scored)')

r=1; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_in_avg(find(iscell==1),:)',[],1))
taskLines(p,info,'full','heatmap',false,0,1);
xlabel('Cell')
ylabel('Bin (~m)')
title('M (z-scored)')

r=1; c=5; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_out_zcrFiltered_avg(find(iscell==1),:)',[],1))
taskLines(p,info,'full','heatmap',false,0,1);
xlabel('Cell')
ylabel('Bin (~m)')
title('M_{hat} decontaminated (z-scored)')

r=2; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U_avg(1:10,:)')
taskLines(p,info,'full','heatmap',false,0,1);
title('Left (temporal) singular vectors U')
xlabel('Component (m)')
ylabel('Bin (~m)')

r=2; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(S(1:10,1:10)')
title('Singular values S')
xlabel('Component (n)')
ylabel('Component (m)')

r=2; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(V',[0,0.05])
title('Right ("spatial") singular vectors V^T')
xlabel('ROI or Component (n)')
ylabel('ROI or Component (n)')

r=2; c=4:5; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
plot(zcrd)
yline(zcrd_threshold,':');
xlabel('Temporal singular vectors')
ylabel('Absolute derivative of zero crossing rate')
title([num2str(sum(1-singularValuesToKeep)),' components removed'])

suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',num2str(info.stimSession)])
savefig(F,[path.filepart_out,'plots\',info.animal,'_',info.date,'_','svdSummaryFigure.fig']);
saveas(F,[path.filepart_out,'plots\',info.animal,'_',info.date,'_','svdSummaryFigure.png']);
disp(['--- Saved SVD summary figure to ',path.filepart_out,'plots.'])
drawnow;


%% SVD removed components figure

nrows = 2; ncols = 10;
F = default_figure();

temp = find(singularValuesToKeep==0);
temp1 = sum(1-singularValuesToKeep);
if temp1<nrows*ncols
    temp2 = temp1;
else
    temp2 = nrows*ncols;
end
for i=1:temp2
    subplot(nrows,ncols,i);
    imagesc(remU_avg{i}'); hold on;
    %taskLines(p,info,'full','heatmap',false,0,1);
    xticks([])
    ylabel('Bin (~m)')
    title(['S(',num2str(temp(i)),')=',num2str(S(temp(i),temp(i)),2)])
end

suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',num2str(info.stimSession),', ',num2str(sum(1-singularValuesToKeep)),' components removed'])
savefig(F,[path.filepart_out,'plots\',info.animal,'_',info.date,'_','svdRemovedComponentsFigure.fig']);
saveas(F,[path.filepart_out,'plots\',info.animal,'_',info.date,'_','svdRemovedComponentsFigure.png']);
disp(['--- Saved SVD removed components figure to ',path.filepart_out,'plots.'])
drawnow;


%% Calculate dFF

disp('--- Calculating dFF...')

% reestablish fluorescence baseline
Fn = Fn - nanmedian(Fn,2) + nanmedian(F_all,2);

% calculate F0
Fn_trend_all = movmedian(Fn,1+2*floor(info.scope.frameRate*p.scope.detrendingWindow/2),2);
    
% calculate dFF
dFFn_all = (Fn - Fn_trend_all) ./ Fn_trend_all;

% split into epochs
writeNPY(dFFn_all,[path.filepart_outX,'dFFn_all.npy']);
disp(['--- Added dFFn_all file to repoX as ',[path.filepart_outX,'dFFn_all.npy'],'.'])
[dFFn_pre,dFFn_beh,dFFn_post] = splitPreBehPost(dFFn_all,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);
save([path.filepart_out,'dFFn_pre.mat'],'dFFn_pre','-v7.3');
disp(['--- Added dFFn_pre file to repo as ',[path.filepart_out,'dFFn_pre.mat'],'.'])
save([path.filepart_out,'dFFn_beh.mat'],'dFFn_beh','-v7.3');
disp(['--- Added dFFn_beh file to repo as ',[path.filepart_out,'dFFn_beh.mat'],'.'])
save([path.filepart_out,'dFFn_post.mat'],'dFFn_post','-v7.3');
disp(['--- Added dFFn_post file to repo as ',[path.filepart_out,'dFFn_post.mat'],'.'])


%% Deconvolution

% disp('--- Running suite2p deconvolution...')


%% Save and return

% save new s2p_meta
s2p_meta.singularValuesToKeep = singularValuesToKeep;
s2p_meta = orderfields(s2p_meta);
save([path.filepart_out,'s2p_meta.mat'],'s2p_meta','-v7.3');
disp(['--- Added s2p_meta file to repo as ',[path.filepart_out,'s2p_meta.mat'],'.'])

end

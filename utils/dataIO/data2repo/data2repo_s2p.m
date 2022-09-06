function [path,s2p_meta] = data2repo_s2p(info,p,path,raw,ops)

%% Assign imaging frames to pre, beh and post epochs and sanity checks

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


%% Load s2p_meta data

if ops.s2p.use_mat_data
    disp('--- Loading Fall.mat...')
    path.file_in_s2p = [path.folder_data,'Imaging\suite2p\plane0\Fall.mat'];
    s2p = load(path.file_in_s2p);
else
    
    if ~ops.s2p.skip_convert_ops_and_stat
        system(['python ',fileparts(mfilename('fullpath')),'\convert_ops_and_stat.py ',info.animal,' ',info.date,' ',path.root_data(1:end-1)]);
    end
    
    disp('--- Loading s2p_meta data...')

    % iscell
    path.file_in_s2p_iscell = [path.folder_data,'Imaging\suite2p\plane0\iscell.npy'];
    s2p.iscell = readNPY(path.file_in_s2p_iscell);
    s2p_meta.iscell = s2p.iscell;
    
    % ops
    path.file_in_s2p_ops = [path.folder_data,'Imaging\suite2p\plane0\_ops.mat'];
    s2p.ops = load(path.file_in_s2p_ops);
    s2p.ops = s2p.ops.ops;
    s2p_meta.ops = orderfields(s2p.ops);
    
    % stat
    path.file_in_s2p_stat = [path.folder_data,'Imaging\suite2p\plane0\_stat.mat'];
    s2p.stat = load(path.file_in_s2p_stat);
    s2p.stat = s2p.stat.stat;
    s2p_meta.stat = s2p.stat;
end

% Clear memory
clear s2p;


%% Load fluorescence data and do neuropil subtraction

% Loading Fns (or calculate it from F and Fneu)
disp('--- Loading fluorescence data...')
if ~ops.s2p.import_Fns
    if ~ops.s2p.use_mat_data
        path.file_in_s2p_F = [path.folder_data,'Imaging\suite2p\plane0\F.npy'];
        s2p.F = readNPY(path.file_in_s2p_F);
        path.file_in_s2p_Fneu = [path.folder_data,'Imaging\suite2p\plane0\Fneu.npy'];
        s2p.Fneu = readNPY(path.file_in_s2p_Fneu);
    end
    
    % save F and Fneu
    if ops.s2p.saveFandFneu
        [F_pre,F_beh,F_post] = splitPreBehPost(s2p.F,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);
        save([path.filepart_outX,'F_pre.mat'],'F_pre','-v7.3');
        disp(['--- Added F_pre file to repoX as ',[path.filepart_outX,'F_pre.mat'],'.'])
        save([path.filepart_outX,'F_beh.mat'],'F_beh','-v7.3');
        disp(['--- Added F_beh file to repoX as ',[path.filepart_outX,'F_beh.mat'],'.'])
        save([path.filepart_outX,'F_post.mat'],'F_post','-v7.3');
        disp(['--- Added F_post file to repoX as ',[path.filepart_outX,'F_post.mat'],'.'])
        clear F_pre; clear F_beh; clear F_post; 
        [Fneu_pre,Fneu_beh,Fneu_post] = splitPreBehPost(s2p.Fneu,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);
        save([path.filepart_outX,'Fneu_pre.mat'],'Fneu_pre','-v7.3');
        disp(['--- Added Fneu_pre file to repoX as ',[path.filepart_outX,'Fneu_pre.mat'],'.'])
        save([path.filepart_outX,'Fneu_beh.mat'],'Fneu_beh','-v7.3');
        disp(['--- Added Fneu_beh file to repoX as ',[path.filepart_outX,'Fneu_beh.mat'],'.'])
        save([path.filepart_outX,'Fneu_post.mat'],'Fneu_post','-v7.3');
        disp(['--- Added Fneu_post file to repoX as ',[path.filepart_outX,'Fneu_post.mat'],'.'])
        clear Fneu_pre; clear Fneu_beh; clear Fneu_post; 
    end

    if ops.s2p.skip_neuropil_subtraction
        s2p.Fns = s2p.F;
    else
        s2p.Fns = s2p.F - p.scope.neuropilSubtraction*s2p.Fneu;
    end
    s2p = rmfield(s2p,'F');
    s2p = rmfield(s2p,'Fneu');
else
    path.file_in_s2p_Fns = [path.folder_data,'Imaging\suite2p\plane0\_Fns.npy'];
    s2p.Fns = readNPY(path.file_in_s2p_Fns);
end

% Calculating baseline fluorescence
if ~ops.s2p.skip_baseline_fluorescence
    disp('--- Calculating baseline fluorescence...')
    [Fns_pre,Fns_beh,Fns_post] = splitPreBehPost(s2p.Fns,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);
    for i=1:length(p.scope.baseF_prctiles)
        s2p_meta.base.(['prctile_',num2str(p.scope.baseF_prctiles(i))]).all = prctile(rmmissing(s2p.Fns,2),p.scope.baseF_prctiles(i),2);
        s2p_meta.base.(['prctile_',num2str(p.scope.baseF_prctiles(i))]).pre = prctile(rmmissing(Fns_pre,2),p.scope.baseF_prctiles(i),2);
        s2p_meta.base.(['prctile_',num2str(p.scope.baseF_prctiles(i))]).beh = prctile(rmmissing(Fns_beh,2),p.scope.baseF_prctiles(i),2);
        s2p_meta.base.(['prctile_',num2str(p.scope.baseF_prctiles(i))]).post = prctile(rmmissing(Fns_post,2),p.scope.baseF_prctiles(i),2);
    end
    clear Fns_pre; clear Fns_beh; clear Fns_post; 
end


%% Calculate dFF (with global median as F0)

if ops.s2p.do_dFF_gm
    
    % Calculating dFF_gm (global median as F0) and saving results
    disp('--- Calculating dFF_gm (global median as F0)...')
    temp = nanmedian(s2p.Fns,2);
    s2p.dFF_gm = (s2p.Fns-temp)./temp;
    [dFF_gm_pre,dFF_gm_beh,dFF_gm_post] = splitPreBehPost(s2p.dFF_gm,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);
    save([path.filepart_out,'dFF_gm_pre.mat'],'dFF_gm_pre','-v7.3');
    disp(['--- Added dFF_gm_pre file to repo as ',[path.filepart_out,'dFF_gm_pre.mat'],'.'])
    save([path.filepart_out,'dFF_gm_beh.mat'],'dFF_gm_beh','-v7.3');
    disp(['--- Added dFF_gm_beh file to repo as ',[path.filepart_out,'dFF_gm_beh.mat'],'.'])
    save([path.filepart_out,'dFF_gm_post.mat'],'dFF_gm_post','-v7.3');
    disp(['--- Added dFF_gm_post file to repo as ',[path.filepart_out,'dFF_gm_post.mat'],'.'])

    % Calculating noise level of dFF_gm
    s2p_meta.noise.dFF_gm.all = 100*nanmedian(abs(diff(s2p.dFF_gm,[],2)),2)/sqrt(info.scope.frameRate);
    s2p_meta.noise.dFF_gm.pre = 100*nanmedian(abs(diff(dFF_gm_pre,[],2)),2)/sqrt(info.scope.frameRate);
    s2p_meta.noise.dFF_gm.beh = 100*nanmedian(abs(diff(dFF_gm_beh,[],2)),2)/sqrt(info.scope.frameRate);
    s2p_meta.noise.dFF_gm.post = 100*nanmedian(abs(diff(dFF_gm_post,[],2)),2)/sqrt(info.scope.frameRate);

    % % Z-scoring epoch-wise dFF_gm and saving results
    % disp('--- Z-scoring epoch-wise dFF_gm...')
    % zdFF_gm_pre = nanzscore(dFF_gm_pre,[],2);
    % save([path.filepart_out,'zdFF_gm_pre.mat'],'zdFF_gm_pre','-v7.3');
    % disp(['--- Added zdFF_gm_pre file to repo as ',[path.filepart_out,'zdFF_gm_pre.mat'],'.'])
    % clear zdFF_gm_pre;
    % zdFF_gm_beh = nanzscore(dFF_gm_beh,[],2);
    % save([path.filepart_out,'zdFF_gm_beh.mat'],'zdFF_gm_beh','-v7.3');
    % disp(['--- Added zdFF_gm_beh file to repo as ',[path.filepart_out,'zdFF_gm_beh.mat'],'.'])
    % clear zdFF_gm_beh;
    % zdFF_gm_post = nanzscore(dFF_gm_post,[],2);
    % save([path.filepart_out,'zdFF_gm_post.mat'],'zdFF_gm_post','-v7.3');
    % disp(['--- Added zdFF_gm_post file to repo as ',[path.filepart_out,'zdFF_gm_post.mat'],'.'])
    % clear zdFF_gm_post;

    % Clear memory
    clear dFF_gm_pre; clear dFF_gm_beh; clear dFF_gm_post; 
    s2p = rmfield(s2p,'dFF_gm');
end


%% Calculate dFF (with moving median as F0)

if ~ops.s2p.skip_dFF
    
    % Calculating F0
    disp('--- Calculating F0 (moving median)...')
    s2p.Fns_trend = movmedian(s2p.Fns,1+2*floor(info.scope.frameRate*p.scope.detrendingWindow/2),2);
    
    % Calculating dFF
    disp('--- Calculating dFF...')
    s2p.dFF = (s2p.Fns - s2p.Fns_trend) ./ s2p.Fns_trend;
    
%     % Loading dFF_mm
%     disp('--- Loading dFF_mm data...')
%     path.file_in_s2p_dFF_mm = [path.folder_data,'Imaging\suite2p\plane0\_dFF_mp_50_60.npy'];
%     s2p.dFF = readNPY(path.file_in_s2p_dFF_mm);
%     
    % Splitting dFF
    [dFF_pre,dFF_beh,dFF_post] = splitPreBehPost(s2p.dFF,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);

    % Calculating noise level of dFF
    s2p_meta.noise.dFF.all = 100*nanmedian(abs(diff(s2p.dFF,[],2)),2)/sqrt(info.scope.frameRate);
    s2p_meta.noise.dFF.pre = 100*nanmedian(abs(diff(dFF_pre,[],2)),2)/sqrt(info.scope.frameRate);
    s2p_meta.noise.dFF.beh = 100*nanmedian(abs(diff(dFF_beh,[],2)),2)/sqrt(info.scope.frameRate);
    s2p_meta.noise.dFF.post = 100*nanmedian(abs(diff(dFF_post,[],2)),2)/sqrt(info.scope.frameRate);

    % Saving dFF to repo and clear memory
    save([path.filepart_out,'dFF_pre.mat'],'dFF_pre','-v7.3');
    disp(['--- Added dFF_pre file to repo as ',[path.filepart_out,'dFF_pre.mat'],'.'])
    save([path.filepart_out,'dFF_beh.mat'],'dFF_beh','-v7.3');
    disp(['--- Added dFF_beh file to repo as ',[path.filepart_out,'dFF_beh.mat'],'.'])
    save([path.filepart_out,'dFF_post.mat'],'dFF_post','-v7.3');
    disp(['--- Added dFF_post file to repo as ',[path.filepart_out,'dFF_post.mat'],'.'])
    clear dFF_pre; clear dFF_beh; clear dFF_post; 
    dFF = s2p.dFF;
    clear s2p;
    
    % Saving dFF.npy and clear memory
%     save([path.folder_data,'Imaging\suite2p\plane0\_dFF.mat'],'dFF','-v7.3');
%     disp(['--- Added dFF file to data path as ',[path.folder_data,'Imaging\suite2p\plane0\_dFF.mat'],'.'])
    %writeNPY(dFF,[path.folder_data,'Imaging\suite2p\plane0\_dFF_',info.animal,'_',info.date,'.npy']);
    %disp(['--- Added dFF file to data path as ',[path.folder_data,'Imaging\suite2p\plane0\_dFF_',info.animal,'_',info.date,'.npy'],'.'])
    temp = regexp(path.root_data,'\','split');
    if ~exist([temp{1},'\Cascade\',info.animal],'dir')
        mkdir([temp{1},'\Cascade\',info.animal]);
    end
    writeNPY(dFF,[temp{1},'\Cascade\',info.animal,'\_dFF_',info.animal,'_',info.date,'.npy']);
    disp(['--- Added dFF file to data path as ',[temp{1},'\Cascade\',info.animal,'\_dFF_',info.animal,'_',info.date,'.npy'],'.'])
    clear dFF;
end


%% Load, process and save s2p spks data

if ~ops.s2p.skip_spks
    
    % Loading spks
    if ~ops.s2p.use_mat_data
        disp('--- Loading suite2p spks data...')
        path.file_in_s2p_spks = [path.folder_data,'Imaging\suite2p\plane0\spks.npy'];
        s2p.spks = readNPY(path.file_in_s2p_spks);
    end
    
    % Splitting spks
    [spks_pre,spks_beh,spks_post] = splitPreBehPost(s2p.spks,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);

    % Calculating signal noise of spks
    disp('--- Calculating signal noise of suite2p spks...')
    try
        s2p_meta.noise.spks.all = GetSn(rmmissing(s2p.spks,2));
        s2p_meta.noise.spks.pre = GetSn(rmmissing(spks_pre,2));
        s2p_meta.noise.spks.beh = GetSn(rmmissing(spks_beh,2));
        s2p_meta.noise.spks.post = GetSn(rmmissing(spks_post,2));
    catch
        warning('Didnt calculate signal noise of s2p spks due to out of memory error')
    end

    % Saving spks
    save([path.filepart_outX,'spks_pre.mat'],'spks_pre','-v7.3');
    disp(['--- Added spks_pre file to repoX as ',[path.filepart_outX,'spks_pre.mat'],'.'])
    save([path.filepart_outX,'spks_beh.mat'],'spks_beh','-v7.3');
    disp(['--- Added spks_beh file to repoX as ',[path.filepart_outX,'spks_beh.mat'],'.'])
    save([path.filepart_outX,'spks_post.mat'],'spks_post','-v7.3');
    disp(['--- Added spks_post file to repoX as ',[path.filepart_outX,'spks_post.mat'],'.'])

    % Clear memory
    clear spks_pre; clear spks_beh; clear spks_post; 
    clear s2p;
end


%% Save s2p_meta

disp('--- Saving s2p_meta data...')

% Saving s2p_meta
s2p_meta.raw = raw;
s2p_meta.noise = orderfields(s2p_meta.noise);
s2p_meta = orderfields(s2p_meta);
save([path.filepart_out,'s2p_meta.mat'],'s2p_meta','-v7.3');
disp(['--- Added s2p_meta file to repo as ',[path.filepart_out,'s2p_meta.mat'],'.'])

% Clear memory
clear s2p;


end

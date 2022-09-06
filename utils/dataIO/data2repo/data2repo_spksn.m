function [path] = data2repo_spksn(info,p,path,ops,s2p_meta,raw)

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


%% Load and process data

disp('--- Loading and processing spks data...')

% load spks_all data
spksn_all = readNPY([path.filepart_outX,'spksn_all.npy']);

% split into pre, beh, post
[spksn_pre,spksn_beh,spksn_post] = splitPreBehPost(spksn_all,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);

% save
save([path.filepart_out,'spksn_pre.mat'],'spksn_pre','-v7.3');
disp(['--- Added spksn_pre file to repo as ',[path.filepart_out,'spksn_pre.mat'],'.'])
save([path.filepart_out,'spksn_beh.mat'],'spksn_beh','-v7.3');
disp(['--- Added spksn_beh file to repo as ',[path.filepart_out,'spksn_beh.mat'],'.'])
save([path.filepart_out,'spksn_post.mat'],'spksn_post','-v7.3');
disp(['--- Added spksn_post file to repo as ',[path.filepart_out,'spksn_post.mat'],'.'])


%% Calculate noise level

% calculate signal noise of spks
disp('--- Calculating signal noise of spksn...')
try
    s2p_meta.noise.spksn.all = GetSn(rmmissing(spksn_all,2));
    s2p_meta.noise.spksn.pre = GetSn(rmmissing(spksn_pre,2));
    s2p_meta.noise.spksn.beh = GetSn(rmmissing(spksn_beh,2));
    s2p_meta.noise.spksn.post = GetSn(rmmissing(spksn_post,2));
catch
    warning('Didnt calculate signal noise of spksn due to out of memory error')
end

% save new s2p_meta
s2p_meta = orderfields(s2p_meta);
save([path.filepart_out,'s2p_meta.mat'],'s2p_meta','-v7.3');
disp(['--- Added s2p_meta file to repo as ',[path.filepart_out,'s2p_meta.mat'],'.'])


end





function [path] = data2repo_casc(info,path,raw,ops,s2p_meta)

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


%% Load, process and save CASCADE 50g data

if ops.casc.do_50g
    
    % Loading casc
    disp('--- Loading CASCADE 50g data...')
    temp = regexp(path.root_data,'\','split');
    path.file_in_casc_50g = [temp{1},'\Cascade\',info.animal,'_out\_casc_',info.animal,'_',info.date,'_50g.npy'];
    s2p.casc_50g = readNPY(path.file_in_casc_50g);
    
    % Splitting casc
    [casc_50g_pre,casc_50g_beh,casc_50g_post] = splitPreBehPost(s2p.casc_50g,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);

    % Calculating firing properties
    disp('--- Calculating firing properties in 50g data...')
    s2p_meta.firing.casc_50g.numSpikes.all = nansum(s2p.casc_50g,2);
    s2p_meta.firing.casc_50g.numSpikes.pre = nansum(casc_50g_pre,2);
    s2p_meta.firing.casc_50g.numSpikes.beh = nansum(casc_50g_beh,2);
    s2p_meta.firing.casc_50g.numSpikes.post = nansum(casc_50g_post,2);
    s2p_meta.firing.casc_50g.avgRate.all = s2p_meta.firing.casc_50g.numSpikes.all*info.scope.frameRate/(s2pFrames_pre+s2pFrames_beh+s2pFrames_post);
    s2p_meta.firing.casc_50g.avgRate.pre = s2p_meta.firing.casc_50g.numSpikes.pre*info.scope.frameRate/s2pFrames_pre;
    s2p_meta.firing.casc_50g.avgRate.beh = s2p_meta.firing.casc_50g.numSpikes.beh*info.scope.frameRate/s2pFrames_beh;
    s2p_meta.firing.casc_50g.avgRate.post = s2p_meta.firing.casc_50g.numSpikes.post*info.scope.frameRate/s2pFrames_post;

    % Saving casc
    save([path.filepart_outX,'casc_50g_pre.mat'],'casc_50g_pre','-v7.3');
    disp(['--- Added casc_50g_pre file to repoX as ',[path.filepart_outX,'casc_50g_pre.mat'],'.'])
    save([path.filepart_outX,'casc_50g_beh.mat'],'casc_50g_beh','-v7.3');
    disp(['--- Added casc_50g_beh file to repoX as ',[path.filepart_outX,'casc_50g_beh.mat'],'.'])
    save([path.filepart_outX,'casc_50g_post.mat'],'casc_50g_post','-v7.3');
    disp(['--- Added casc_50g_post file to repoX as ',[path.filepart_outX,'casc_50g_post.mat'],'.'])

    % Clear memory
    clear casc_50g_pre; clear casc_50g_beh; clear casc_50g_post; 
    clear s2p;
end


%% Load, process and save CASCADE 50c data

if ops.casc.do_50c
    
    % Loading casc
    disp('--- Loading CASCADE 50c data...')
    temp = regexp(path.root_data,'\','split');
    path.file_in_casc_50c = [temp{1},'\Cascade\',info.animal,'_out\_casc_',info.animal,'_',info.date,'_50c.npy'];
    s2p.casc_50c = readNPY(path.file_in_casc_50c);
    
    % Splitting casc
    [casc_50c_pre,casc_50c_beh,casc_50c_post] = splitPreBehPost(s2p.casc_50c,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);

    % Calculating firing properties
    disp('--- Calculating firing properties in 50c data...')
    s2p_meta.firing.casc_50c.numSpikes.all = nansum(s2p.casc_50c,2);
    s2p_meta.firing.casc_50c.numSpikes.pre = nansum(casc_50c_pre,2);
    s2p_meta.firing.casc_50c.numSpikes.beh = nansum(casc_50c_beh,2);
    s2p_meta.firing.casc_50c.numSpikes.post = nansum(casc_50c_post,2);
    s2p_meta.firing.casc_50c.avgRate.all = s2p_meta.firing.casc_50c.numSpikes.all*info.scope.frameRate/(s2pFrames_pre+s2pFrames_beh+s2pFrames_post);
    s2p_meta.firing.casc_50c.avgRate.pre = s2p_meta.firing.casc_50c.numSpikes.pre*info.scope.frameRate/s2pFrames_pre;
    s2p_meta.firing.casc_50c.avgRate.beh = s2p_meta.firing.casc_50c.numSpikes.beh*info.scope.frameRate/s2pFrames_beh;
    s2p_meta.firing.casc_50c.avgRate.post = s2p_meta.firing.casc_50c.numSpikes.post*info.scope.frameRate/s2pFrames_post;

    % Saving casc
    save([path.filepart_outX,'casc_50c_pre.mat'],'casc_50c_pre','-v7.3');
    disp(['--- Added casc_50c_pre file to repoX as ',[path.filepart_outX,'casc_50c_pre.mat'],'.'])
    save([path.filepart_outX,'casc_50c_beh.mat'],'casc_50c_beh','-v7.3');
    disp(['--- Added casc_50c_beh file to repoX as ',[path.filepart_outX,'casc_50c_beh.mat'],'.'])
    save([path.filepart_outX,'casc_50c_post.mat'],'casc_50c_post','-v7.3');
    disp(['--- Added casc_50c_post file to repoX as ',[path.filepart_outX,'casc_50c_post.mat'],'.'])

    % Clear memory
    clear casc_50c_pre; clear casc_50c_beh; clear casc_50c_post; 
    clear s2p;
end


%% Load, process and save CASCADE 100g data

if ops.casc.do_100g
    
    % Loading casc
    disp('--- Loading CASCADE 100g data...')
    temp = regexp(path.root_data,'\','split');
    path.file_in_casc_100g = [temp{1},'\Cascade\',info.animal,'_out\_casc_',info.animal,'_',info.date,'_100g.npy'];
    s2p.casc_100g = readNPY(path.file_in_casc_100g);
    
    % Splitting casc
    [casc_100g_pre,casc_100g_beh,casc_100g_post] = splitPreBehPost(s2p.casc_100g,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);

    % Calculating firing properties
    disp('--- Calculating firing properties in 100g data...')
    s2p_meta.firing.casc_100g.numSpikes.all = nansum(s2p.casc_100g,2);
    s2p_meta.firing.casc_100g.numSpikes.pre = nansum(casc_100g_pre,2);
    s2p_meta.firing.casc_100g.numSpikes.beh = nansum(casc_100g_beh,2);
    s2p_meta.firing.casc_100g.numSpikes.post = nansum(casc_100g_post,2);
    s2p_meta.firing.casc_100g.avgRate.all = s2p_meta.firing.casc_100g.numSpikes.all*info.scope.frameRate/(s2pFrames_pre+s2pFrames_beh+s2pFrames_post);
    s2p_meta.firing.casc_100g.avgRate.pre = s2p_meta.firing.casc_100g.numSpikes.pre*info.scope.frameRate/s2pFrames_pre;
    s2p_meta.firing.casc_100g.avgRate.beh = s2p_meta.firing.casc_100g.numSpikes.beh*info.scope.frameRate/s2pFrames_beh;
    s2p_meta.firing.casc_100g.avgRate.post = s2p_meta.firing.casc_100g.numSpikes.post*info.scope.frameRate/s2pFrames_post;

    % Saving casc
    save([path.filepart_outX,'casc_100g_pre.mat'],'casc_100g_pre','-v7.3');
    disp(['--- Added casc_100g_pre file to repoX as ',[path.filepart_outX,'casc_100g_pre.mat'],'.'])
    save([path.filepart_outX,'casc_100g_beh.mat'],'casc_100g_beh','-v7.3');
    disp(['--- Added casc_100g_beh file to repoX as ',[path.filepart_outX,'casc_100g_beh.mat'],'.'])
    save([path.filepart_outX,'casc_100g_post.mat'],'casc_100g_post','-v7.3');
    disp(['--- Added casc_100g_post file to repoX as ',[path.filepart_outX,'casc_100g_post.mat'],'.'])

    % Clear memory
    clear casc_100g_pre; clear casc_100g_beh; clear casc_100g_post; 
    clear s2p;
end


%% Load, process and save CASCADE 100c data

if ops.casc.do_100c
    
    % Loading casc
    disp('--- Loading CASCADE 100c data...')
    temp = regexp(path.root_data,'\','split');
    path.file_in_casc_100c = [temp{1},'\Cascade\',info.animal,'_out\_casc_',info.animal,'_',info.date,'_100c.npy'];
    s2p.casc_100c = readNPY(path.file_in_casc_100c);
    
    % Splitting casc
    [casc_100c_pre,casc_100c_beh,casc_100c_post] = splitPreBehPost(s2p.casc_100c,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);

    % Calculating firing properties
    disp('--- Calculating firing properties in 100c data...')
    s2p_meta.firing.casc_100c.numSpikes.all = nansum(s2p.casc_100c,2);
    s2p_meta.firing.casc_100c.numSpikes.pre = nansum(casc_100c_pre,2);
    s2p_meta.firing.casc_100c.numSpikes.beh = nansum(casc_100c_beh,2);
    s2p_meta.firing.casc_100c.numSpikes.post = nansum(casc_100c_post,2);
    s2p_meta.firing.casc_100c.avgRate.all = s2p_meta.firing.casc_100c.numSpikes.all*info.scope.frameRate/(s2pFrames_pre+s2pFrames_beh+s2pFrames_post);
    s2p_meta.firing.casc_100c.avgRate.pre = s2p_meta.firing.casc_100c.numSpikes.pre*info.scope.frameRate/s2pFrames_pre;
    s2p_meta.firing.casc_100c.avgRate.beh = s2p_meta.firing.casc_100c.numSpikes.beh*info.scope.frameRate/s2pFrames_beh;
    s2p_meta.firing.casc_100c.avgRate.post = s2p_meta.firing.casc_100c.numSpikes.post*info.scope.frameRate/s2pFrames_post;

    % Saving casc
    save([path.filepart_outX,'casc_100c_pre.mat'],'casc_100c_pre','-v7.3');
    disp(['--- Added casc_100c_pre file to repoX as ',[path.filepart_outX,'casc_100c_pre.mat'],'.'])
    save([path.filepart_outX,'casc_100c_beh.mat'],'casc_100c_beh','-v7.3');
    disp(['--- Added casc_100c_beh file to repoX as ',[path.filepart_outX,'casc_100c_beh.mat'],'.'])
    save([path.filepart_outX,'casc_100c_post.mat'],'casc_100c_post','-v7.3');
    disp(['--- Added casc_100c_post file to repoX as ',[path.filepart_outX,'casc_100c_post.mat'],'.'])

    % Clear memory
    clear casc_100c_pre; clear casc_100c_beh; clear casc_100c_post; 
    clear s2p;
end


%% Load, process and save CASCADE 200g data

if ops.casc.do_200g
    
    % Loading casc
    disp('--- Loading CASCADE 200g data...')
    temp = regexp(path.root_data,'\','split');
    path.file_in_casc_200g = [temp{1},'\Cascade\',info.animal,'_out\_casc_',info.animal,'_',info.date,'_200g.npy'];
    s2p.casc_200g = readNPY(path.file_in_casc_200g);
    
    % Splitting casc
    [casc_200g_pre,casc_200g_beh,casc_200g_post] = splitPreBehPost(s2p.casc_200g,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post);

    % Calculating firing properties
    disp('--- Calculating firing properties in 200g data...')
    s2p_meta.firing.casc_200g.numSpikes.all = nansum(s2p.casc_200g,2);
    s2p_meta.firing.casc_200g.numSpikes.pre = nansum(casc_200g_pre,2);
    s2p_meta.firing.casc_200g.numSpikes.beh = nansum(casc_200g_beh,2);
    s2p_meta.firing.casc_200g.numSpikes.post = nansum(casc_200g_post,2);
    s2p_meta.firing.casc_200g.avgRate.all = s2p_meta.firing.casc_200g.numSpikes.all*info.scope.frameRate/(s2pFrames_pre+s2pFrames_beh+s2pFrames_post);
    s2p_meta.firing.casc_200g.avgRate.pre = s2p_meta.firing.casc_200g.numSpikes.pre*info.scope.frameRate/s2pFrames_pre;
    s2p_meta.firing.casc_200g.avgRate.beh = s2p_meta.firing.casc_200g.numSpikes.beh*info.scope.frameRate/s2pFrames_beh;
    s2p_meta.firing.casc_200g.avgRate.post = s2p_meta.firing.casc_200g.numSpikes.post*info.scope.frameRate/s2pFrames_post;

    % Saving casc
    save([path.filepart_outX,'casc_200g_pre.mat'],'casc_200g_pre','-v7.3');
    disp(['--- Added casc_200g_pre file to repoX as ',[path.filepart_outX,'casc_200g_pre.mat'],'.'])
    save([path.filepart_outX,'casc_200g_beh.mat'],'casc_200g_beh','-v7.3');
    disp(['--- Added casc_200g_beh file to repoX as ',[path.filepart_outX,'casc_200g_beh.mat'],'.'])
    save([path.filepart_outX,'casc_200g_post.mat'],'casc_200g_post','-v7.3');
    disp(['--- Added casc_200g_post file to repoX as ',[path.filepart_outX,'casc_200g_post.mat'],'.'])

    % Clear memory
    clear casc_200g_pre; clear casc_200g_beh; clear casc_200g_post; 
    clear s2p;
end


%% Save s2p_meta

disp('--- Overwriting s2p_meta data...')

% Saving s2p_meta
s2p_meta = orderfields(s2p_meta);
save([path.filepart_out,'s2p_meta.mat'],'s2p_meta','-v7.3');
disp(['--- Added s2p_meta file to repo as ',[path.filepart_out,'s2p_meta.mat'],'.'])

% Clear memory
clear s2p;

end

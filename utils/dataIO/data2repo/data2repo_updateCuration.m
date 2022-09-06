function [path,s2p_meta] = data2repo_updateCuration(info,path)

%% Load scope data

disp('--- Loading s2p_meta and new curation...')
  
% s2p_meta
load([path.filepart_out,'s2p_meta.mat']);

% load new curation
path.file_in_s2p_iscell = [path.folder_data,'Imaging\suite2p\plane0\iscell.npy'];
s2p.iscell = readNPY(path.file_in_s2p_iscell);


%% Process suite2p data

s2p_meta.iscell = s2p.iscell;


%% Save suite2p file

s2p_meta = orderfields(s2p_meta);
save([path.filepart_out,'s2p_meta.mat'],'s2p_meta','-v7.3');
disp(['--- Updated and saved s2p_meta file to repo as ',[path.filepart_out,'s2p_meta.mat'],'.'])

end
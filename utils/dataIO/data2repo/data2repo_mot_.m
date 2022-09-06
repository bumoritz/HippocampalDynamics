function [path] = data2repo_mot_(info,p,path)

path.root_data                      = 'G:\Data\'; %  %'F:\Data\'; % 'Z:\WIBR_Hippos\SniffinHippo\Data\'
path.root_repo                      = 'D:\SniffinHippo\Repo\'; % 'D:\temp\';
path.root_repoX                     = 'E:\SniffinHippo\RepoX\'; % 'D:\temp\';

% gather info and p
p = get_p;
info = get_info;
info.animal = animal;
info.date = date;
info.expDay = expDay;
info.stimSession = stimSession;

% identify session folders
path.folder_data = [path.root_data,info.date(1:4),'\',info.date(1:4),'-',info.date(5:6),'\',info.date(1:4),'-',info.date(5:6),'-',info.date(7:8),'\',info.animal,'\'];
path.folder_repo = [path.root_repo,info.animal,'\',info.animal,'_',info.date,'\'];
path.filepart_out = [path.folder_repo,info.animal,'_',info.date,'_'];
if ~exist(path.folder_repo,'dir')
   mkdir(path.folder_repo);
end
path.folder_repoX = [path.root_repoX,info.animal,'\',info.animal,'_',info.date,'\'];
path.filepart_outX = [path.folder_repoX,info.animal,'_',info.date,'_'];
if ~exist(path.folder_repoX,'dir')
   mkdir(path.folder_repoX);
end

if ~exist([path.filepart_out,'plots'],'dir')
    mkdir([path.filepart_out,'plots']);
end
if ~exist([path.filepart_out,'log'],'dir')
    mkdir([path.filepart_out,'log']);
end


%% Load and process camera data

disp('--- Loading and processing camera data...')

% find videos
if ~isempty(dir([path.folder_data,'Camera\',info.animal,'_',info.date,'_*.avi']))
    input_files = dir([path.folder_data,'Camera\',info.animal,'_',info.date,'_*.avi']);
elseif ~isempty(dir([path.folder_data,'Camera\',info.animal,'-',info.date,'-*.avi']))
    input_files = dir([path.folder_data,'Camera\',info.animal,'-',info.date,'-*.avi']);
else
    disp(['Did not find any videos in ',[path.folder_data,'Camera\']])
end

% process videos one-by-one
numFiles = length(input_files);
mot_all = cell(3,numFiles);
for n=1:numFiles
    
    % load video
    mot_all{1,n} = input_files(n).name;
    videoReader = VideoReader([path.folder_data,'Camera\',input_files(n).name]);

    % calculate motion energy
    this_numFrames = videoReader.NumFrames;
    mot_all{2,n} = this_numFrames;
    this_mot = nan(1,this_numFrames);
    for k = 1:this_numFrames-1
        this_frame = rgb2gray(read(videoReader,k));
        next_frame = rgb2gray(read(videoReader,k+1));
        this_mot(k+1) = nansum(abs(next_frame(:)-this_frame(:)));
    end 
    mot_all{3,n} = this_mot;
    
    disp(['--- Processed file ',num2str(n),' / ',num2str(numFiles)])
end


%% Save results

save([path.filepart_outX,'mot_all.mat'],'mot_all','-v7.3');
disp(['--- Added mot_all file to repo as ',[path.filepart_outX,'mot_all.mat'],'.'])


end





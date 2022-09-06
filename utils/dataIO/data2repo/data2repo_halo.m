function [path] = data2repo_halo(info,p,path,raw,ops,s2p_meta,paq_beh)

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


%% Define background

this_meanImage = s2p_meta.ops.meanImg;
this_meanImage_inset = this_meanImage(p.halo.bkgd_inset+1:size(this_meanImage,1)-p.halo.bkgd_inset,p.halo.bkgd_inset+1:size(this_meanImage,1)-p.halo.bkgd_inset);
this_meanImage_inset_lin = this_meanImage_inset(:);
these_background_idcs = find(this_meanImage_inset_lin<prctile(this_meanImage_inset_lin,p.halo.bkgd_prctile));
this_backgroundImage_inset_lin = zeros(size(this_meanImage_inset_lin));
this_backgroundImage_inset_lin(these_background_idcs)=1;
this_backgroundImage_inset = reshape(this_backgroundImage_inset_lin,size(this_meanImage_inset));
this_backgroundImage = zeros(size(this_meanImage));
this_backgroundImage(p.halo.bkgd_inset+1:size(this_meanImage,1)-p.halo.bkgd_inset,p.halo.bkgd_inset+1:size(this_meanImage,1)-p.halo.bkgd_inset) = this_backgroundImage_inset;
temp = this_backgroundImage';
this_backgroundImage_lin = temp(:);

% figure;
% subplot(1,2,1)
% imshow(this_meanImage,[nanmin(this_meanImage(:)),nanmax(this_meanImage(:))/4]);
% title('average image')
% subplot(1,2,2)
% imshow(this_backgroundImage,[0,1]);
% title('background mask')


%% Define halos

iscell = s2p_meta.iscell(:,1);

% get centroids from suite2p output
centroids = {};
centroids.x = nan(length(iscell),1);
centroids.y = nan(length(iscell),1);
for i=1:length(iscell)
    centroids.x(i) = round(s2p_meta.stat{i}.med(2));
    centroids.y(i) = round(s2p_meta.stat{i}.med(1));
end

% generate halos
halos = dilateCentroids_mod(centroids,[],[],p.halo.sigma,[],false, p.halo.halo_inset);

% subtract cell masks from halos
halos_final = halos;
for j=1:length(halos)
    [~,these_pixels]=setdiff(halos{j}.coords,sub2ind([s2p_meta.ops.Lx,s2p_meta.ops.Ly],s2p_meta.stat{j}.xpix,s2p_meta.stat{j}.ypix));
    halos_final{j}.coords = halos{j}.coords(these_pixels);
    halos_final{j}.weights = halos{j}.weights(these_pixels);  
    for i=1:length(s2p_meta.stat{j}.xpix)
        halos_final{j}.image(s2p_meta.stat{j}.ypix(i),s2p_meta.stat{j}.xpix(i)) = 0;
    end
end

% create figure
% idx = 1;
% this_img = zeros(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
% for i=1:length(s2p_meta.stat{idx}.xpix)
%     this_img(s2p_meta.stat{idx}.ypix(i),s2p_meta.stat{idx}.xpix(i)) = 1;
% end
% figure;
% imagesc(halos_final{idx}.image);
% title('example halo')


%% Extracting background and halo traces from registered_movie.raw

disp('--- Extracting background and halo traces...')

path.imagingFile = [path.folder_data,'Imaging\registered_movie.raw'];

temp = dir(path.imagingFile);
numFrames = temp.bytes/(512*512*2);
Fbkgd_all = nan(1,numFrames);
Fhalo_all= nan(numel(halos_final),numFrames);
for i=1:numFrames
    
    % load images one-by-one
    this_fid = fopen(path.imagingFile,'r');
    fseek(this_fid,(i-1)*512*512*2,'bof');
    this_frame = uint16(fread(this_fid,512*512,'uint16',0,'l'));
    frewind(this_fid);
    fclose(this_fid);

    % calculate background
    Fbkgd_all(i) = nanmean(this_frame(find(this_backgroundImage_lin)));
    
    % extract halo traces
    for j = 1:numel(halos_final)
        halo_px = this_frame(halos_final{j}.coords);
        halo_w = halos_final{j}.weights;
        Fhalo_all(j,i) = sum(halo_w .* double(halo_px), 1) ./ sum(halo_w, 1);
    end

    if mod(i,50000)==0
        disp(['--- ',num2str(i),' / ',num2str(numFrames),'...'])
    end
end

% save background and halo traces
save([path.root_repoX,info.animal,'\',info.animal,'_',info.date,'\',info.animal,'_',info.date,'_Fbkgd_all.mat'],'Fbkgd_all','-v7.3');
disp(['--- Added Fbkgd_all file to repoX as ',[path.root_repoX,info.animal,'\',info.animal,'_',info.date,'\',info.animal,'_',info.date,'_Fbkgd_all.mat'],'.'])
save([path.root_repoX,info.animal,'\',info.animal,'_',info.date,'\',info.animal,'_',info.date,'_Fhalo_all.mat'],'Fhalo_all','-v7.3');
disp(['--- Added Fhalo_all file to repoX as ',[path.root_repoX,info.animal,'\',info.animal,'_',info.date,'\',info.animal,'_',info.date,'_Fhalo_all.mat'],'.'])


%% Load and process data

disp('--- Loading and processing fluorescence data...')

% load F_all data
path.file_in_s2p_F = [path.folder_data,'Imaging\suite2p\plane0\F.npy'];
F_all = readNPY(path.file_in_s2p_F);
[F_pre,F_beh,F_post] = splitPreBehPost(F_all,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_pre);
save([path.filepart_outX,'F_pre.mat'],'F_pre','-v7.3');
disp(['--- Added F_pre file to repoX as ',[path.filepart_outX,'F_pre.mat'],'.'])
save([path.filepart_outX,'F_beh.mat'],'F_beh','-v7.3');
disp(['--- Added F_beh file to repoX as ',[path.filepart_outX,'F_beh.mat'],'.'])
save([path.filepart_outX,'F_post.mat'],'F_post','-v7.3');
disp(['--- Added F_post file to repoX as ',[path.filepart_outX,'F_post.mat'],'.'])

% subtract background from data
F_bs_all = F_all - Fbkgd_all;
Fhalo_bs_all = Fhalo_all - Fbkgd_all;

% split into epochs for estimating subtraction factor
[~,F_bs_beh,~] = splitPreBehPost(F_bs_all,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_pre);
[~,Fhalo_bs_beh,~] = splitPreBehPost(Fhalo_bs_all,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_pre);

% average traces
[~,~,F_bs_beh_avg] = preprocessActivityMeasure(F_bs_beh,p.halo,p,paq_beh,iscell);
F_bs_beh_avg = nanmean(F_bs_beh_avg,3);
[~,~,Fhalo_bs_beh_avg] = preprocessActivityMeasure(Fhalo_bs_beh,p.halo,p,paq_beh,iscell);
Fhalo_bs_beh_avg = nanmean(Fhalo_bs_beh_avg,3);

% estimate optimal subtraction factor
warning('off','all'); warning;
s2p_meta.optimalSubtractionFactor = nan(size(F_bs_beh_avg,1),1);
for i=1:size(F_bs_beh_avg,1)
    if iscell(i)==1
        temp = robustfit(Fhalo_bs_beh_avg(i,:),F_bs_beh_avg(i,:));
        s2p_meta.optimalSubtractionFactor(i) = temp(2);
    end
end
warning('on','all'); warning('query','all');

% save new s2p_meta
s2p_meta = orderfields(s2p_meta);
save([path.filepart_out,'s2p_meta.mat'],'s2p_meta','-v7.3');
disp(['--- Added s2p_meta file to repo as ',[path.filepart_out,'s2p_meta.mat'],'.'])

% subtract scaled halo signals
F_hs_bs_all = F_bs_all - s2p_meta.optimalSubtractionFactor .* Fhalo_bs_all;

% reestablish fluorescence baseline
F_hs_bs_all = F_hs_bs_all - nanmedian(F_hs_bs_all,2) + nanmedian(F_all,2);

% calculate F0
F_hs_bs_trend_all = movmedian(F_hs_bs_all,1+2*floor(info.scope.frameRate*p.scope.detrendingWindow/2),2);
    
% calculate dFF
dFF_hs_bs_all = (F_hs_bs_all - F_hs_bs_trend_all) ./ F_hs_bs_trend_all;

% split into epochs
writeNPY(dFF_hs_bs_all,[path.filepart_outX,'dFF_hs_bs_all.npy']);
disp(['--- Added dFF_hs_bs_all file to data path as ',[path.filepart_outX,'dFF_hs_bs_all.npy'],'.'])
[dFF_hs_bs_pre,dFF_hs_bs_beh,dFF_hs_bs_post] = splitPreBehPost(dFF_hs_bs_all,info.scope.numFrames_pre,info.scope.numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_pre);
save([path.filepart_outX,'dFF_hs_bs_pre.mat'],'dFF_hs_bs_pre','-v7.3');
disp(['--- Added dFF_hs_bs_pre file to repoX as ',[path.filepart_outX,'dFF_hs_bs_pre.mat'],'.'])
save([path.filepart_outX,'dFF_hs_bs_beh.mat'],'dFF_hs_bs_beh','-v7.3');
disp(['--- Added dFF_hs_bs_beh file to repoX as ',[path.filepart_outX,'dFF_hs_bs_beh.mat'],'.'])
save([path.filepart_outX,'dFF_hs_bs_post.mat'],'dFF_hs_bs_post','-v7.3');
disp(['--- Added dFF_hs_bs_post file to repoX as ',[path.filepart_outX,'dFF_hs_bs_post.mat'],'.'])


%% Run deconvolution

% disp('--- Running deconvolution...')
% run deconvolution
% system(['conda activate suite2p & python C:\Users\Moritz\Documents\MATLAB\SniffinHippo\utils\dataIO\data2repo\run_suite2p_deconvolution.py']);
% system(['conda activate suite2p & python C:\Users\Moritz\Documents\MATLAB\SniffinHippo\utils\dataIO\data2repo\run_suite2p_deconvolution_TS.py']);
% system(['conda activate suite2p & runipy C:\Users\Moritz\Documents\MATLAB\SniffinHippo\utils\dataIO\data2repo\run_suite2p_deconvolution_TS.ipynb']);
% system(['conda activate suite2p & python C:\Users\Moritz\Documents\MATLAB\SniffinHippo\utils\dataIO\data2repo\run_suite2p_deconvolution_TS.html']);
% system(['conda activate suite2p & jupyter nbconvert --execute C:\Users\Moritz\Documents\MATLAB\SniffinHippo\utils\dataIO\data2repo\run_suite2p_deconvolution_TS.ipynb']);
% C:\Users\Moritz\Documents\MATLAB\SniffinHippo\utils\dataIO\data2repo\run_suite2p_deconvolution_TS.html
% system(['conda activate suite2p & ipython -c C:\Users\Moritz\Documents\MATLAB\SniffinHippo\utils\dataIO\data2repo\run_suite2p_deconvolution_TS.ipynb'])
% spks_all = readNPY([path.filepart_outX,'dFF_hs_bs_post.mat']) %[path.folder_data,'Imaging\suite2p\plane0\Fneu.npy']);
%system(['python ',fileparts(mfilename('fullpath')),'\convert_ops_and_stat.py ',info.animal,' ',info.date,' ',path.root_data(1:end-1)]);


end





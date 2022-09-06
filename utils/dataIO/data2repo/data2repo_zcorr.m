% function [path] = data2repo_zcorr(info,p,path,raw,ops,s2p_meta,paq_beh)

%% Get paths

path.stackPath      = [path.folder_data,'Expression\',info.animal,'-',info.date,'-930A-stack1\';...
    path.folder_data,'Expression\',info.animal,'-',info.date,'-930A-stack2\';...
    path.folder_data,'Expression\',info.animal,'-',info.date,'-930A-stack3\'];
path.imagingFile    = [path.folder_data,'Imaging\registered_movie.raw'];


%% Handle frame numbers for merged videos

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

% s2pFrameIdcs
s2pFrameIdcs_beh = s2pFrames_pre+1:s2pFrames_pre+s2pFrames_beh;


%% Preparations for stack correlation

p.zcorr.stackTop          = 25; % [um]  
p.zcorr.stackBottom       = -25; % [um]
p.zcorr.stackStepSize     = 1; % [um]
p.zcorr.frameRate         = 30.0; % [Hz]
p.zcorr.frameRange        = -30:300;

% Load reference frame and set up registration process
refFrame = uint16(s2p_meta.ops.meanImg);
ops = setup_registration_phasecorr(refFrame);
p.zcorr.bytesPerFrame = size(refFrame,1)*size(refFrame,2)*2; % 2 bytes per value (precision) for uint16

% Load stack
p.zcorr.numSlices = length(p.zcorr.stackBottom:p.zcorr.stackStepSize:p.zcorr.stackTop);
stack = zeros(size(refFrame,1),size(refFrame,2),p.zcorr.numSlices,'uint16');
for i=1:p.zcorr.numSlices
    temp = zeros(size(refFrame,1),size(refFrame,2),size(path.stackPath,1),'uint16');
    for j=1:size(path.stackPath,1)        
        temp2 = dir([path.stackPath(j,:),'Chan*']);
        this_slice_path = [path.stackPath(j,:),temp2(i).name];
        temp(:,:,j) = imread(this_slice_path);
    end
    stack(:,:,i) = mean(temp,3);
end

% Register reference against stack
out0.shift = zeros(p.zcorr.numSlices,2);
out0.corr = zeros(p.zcorr.numSlices,1);
for i=1:p.zcorr.numSlices
    [this_regFrame, out0.shift(i,:), ~] = return_offsets_phasecorr(single(stack(:,:,i)), ops);
    out0.corr(i) = corr2(this_regFrame,refFrame);
end
out0.corr = smoothdata(out0.corr,'gaussian',5);
[~,out0.best_slice] = max(out0.corr);
out0.best_um = p.zcorr.stackTop - (out0.best_slice-1)*p.zcorr.stackStepSize;


%% Frame-wise z-movement (averaged over 500 trials)

avgImages_timeInTrial = zeros(size(refFrame,1),size(refFrame,2),length(p.zcorr.frameRange),'uint16');
avgImages_bestSlice = zeros(length(p.zcorr.frameRange),1);
for k=1:length(p.zcorr.frameRange)
    
    % Load images for 1 specific frame in all trials
    this_fid = fopen(path.imagingFile, 'r');
    fseek(this_fid,0,'bof');
    these_frame_numbers = paq_beh.sync + s2pFrames_pre + p.zcorr.frameRange(k);
    these_frames = zeros(size(refFrame,1),size(refFrame,2),length(paq_beh.sync),'uint16');
    for i=1:length(paq_beh.sync)
        fseek(this_fid,(these_frame_numbers(i)-1)*p.zcorr.bytesPerFrame,'bof');
        temp = uint16(fread(this_fid,size(refFrame,1)*size(refFrame,2),'uint16',0,'l'));
        these_frames(:,:,i) = reshape(temp,size(refFrame,1),size(refFrame,2));
        frewind(this_fid);
    end
    fclose(this_fid);
    
    % Calculate average image and set up registration process
    this_avgImg = nanmean(these_frames,3);
    this_avgImg = this_avgImg';
    this_ops = setup_registration_phasecorr(this_avgImg);
    
    % Register average image against stack
    for i=1:p.zcorr.numSlices
        [this_regFrame, out.shift(i,:), ~] = return_offsets_phasecorr(single(stack(:,:,i)), this_ops);
        out.corr(i) = corr2(this_regFrame,this_avgImg);
    end
    out.corr = smoothdata(out.corr,'gaussian',5);
    [~,out.best_slice_max] = max(out.corr);
    out.best_um_max = p.zcorr.stackTop - (out.best_slice_max-1)*p.zcorr.stackStepSize;

    % Return
    avgImages_timeInTrial(:,:,k) = this_avgImg;
    avgImages_bestSlice(k) = out.best_slice_max;
    k
end

% make movie
for k=1:length(p.zcorr.frameRange)
    if k==1
        imwrite(avgImages_timeInTrial(:,:,k),['C:\SniffinHippo\Snaps\RunningStory\z-problem\','framewiseZmovement.tif'],'Compression','none');
    else
        imwrite(avgImages_timeInTrial(:,:,k),['C:\SniffinHippo\Snaps\RunningStory\z-problem\','framewiseZmovement.tif'],'WriteMode','append','Compression','none');
    end
end

% make z-scored movie
this_baseline_mean = nanmean(avgImages_timeInTrial(:,:,1:30),3);
this_baseline_std = nanstd(double(avgImages_timeInTrial(:,:,1:30)),[],3);
avgImages_timeInTrial_zscored = (double(avgImages_timeInTrial) - this_baseline_mean) ./ this_baseline_std;
avgImages_timeInTrial_zscored_uint16 = im2uint16(rescale(avgImages_timeInTrial_zscored));
for k=1:length(p.zcorr.frameRange)
    if k==1
        imwrite(avgImages_timeInTrial_zscored_uint16(:,:,k),['C:\SniffinHippo\Snaps\RunningStory\z-problem\','framewiseZmovement_zscore.tif'],'Compression','none');
    else
        imwrite(avgImages_timeInTrial_zscored_uint16(:,:,k),['C:\SniffinHippo\Snaps\RunningStory\z-problem\','framewiseZmovement_zscore.tif'],'WriteMode','append','Compression','none');
    end
end

%%

figure;
imshow(avgImages_timeInTrial_zscored_uint16(:,:,200),[0,500000])












%% Frame-wise z-movement (averaged over 500 trials) - plot

%p.zcorr.frameRange
figure;
plot(avgImages_bestSlice - out0.best_slice)
hold on
xline(31)
xline(41)
xline(191)
xline(201)
% ylim([-5,5])
xlabel('frame')
ylabel('z distance from reference frame (1 step = 1 um)')
title('Frame-wise z-movement (averaged over 500 trials)')


%% trial-wise z-movement (averaged over many frames)

avgImages_trialInSession = zeros(size(refFrame,1),size(refFrame,2),length(paq_beh.sync),'uint16');
bestSlice_trialInSession = zeros(length(paq_beh.sync),1);
for k=1:length(paq_beh.sync)
    
    % Load images for 1 specific frame in all trials
    this_fid = fopen(path.imagingFile, 'r');
    fseek(this_fid,0,'bof');
    these_frame_numbers = paq_beh.sync(k) + s2pFrames_pre + p.zcorr.frameRange;
    these_frames = zeros(size(refFrame,1),size(refFrame,2),length(p.zcorr.frameRange),'uint16');
    for i=1:length(p.zcorr.frameRange)
        fseek(this_fid,(these_frame_numbers(i)-1)*p.zcorr.bytesPerFrame,'bof');
        temp = uint16(fread(this_fid,size(refFrame,1)*size(refFrame,2),'uint16',0,'l'));
        these_frames(:,:,i) = reshape(temp,size(refFrame,1),size(refFrame,2));
        frewind(this_fid);
    end
    fclose(this_fid);
    
    % Calculate average image and set up registration process
    this_avgImg = nanmean(these_frames,3);
    this_avgImg = this_avgImg';
    this_ops = setup_registration_phasecorr(this_avgImg);
    
    % Register average image against stack
    for i=1:p.zcorr.numSlices
        [this_regFrame, out.shift(i,:), ~] = return_offsets_phasecorr(single(stack(:,:,i)), this_ops);
        out.corr(i) = corr2(this_regFrame,this_avgImg);
    end
    out.corr = smoothdata(out.corr,'gaussian',5);
    [~,out.best_slice_max] = max(out.corr);
    out.best_um_max = p.zcorr.stackTop - (out.best_slice_max-1)*p.zcorr.stackStepSize;

    % Return
    avgImages_trialInSession(:,:,k) = this_avgImg;
    bestSlice_trialInSession(k) = out.best_slice_max;
    k
end


%% Frame-wise z-movement (averaged over 500 trials) - plot

%p.zcorr.frameRange
figure;
plot(bestSlice_trialInSession - out0.best_slice)
hold on
%ylim([-10,10])
xlabel('frame')
ylabel('z distance from reference frame (1 step = 1 um)')
title('Trial-wise z-movement (averaged over many frames from entire trial)')







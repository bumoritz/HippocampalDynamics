function [path] = data2repo_cam(info,p,path,ops,s2p_meta)

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


%% Load and process mot data

if ops.cam.do_mot
    disp('--- Loading and processing mot data...')
    
    load([path.root_data,'MOTout\',info.animal,'\',info.animal,'_',info.date,'\',info.animal,'_',info.date,'_mot_all.mat']);
end


%% Load and process dlc data

if ops.cam.do_dlc
    disp('--- Loading and processing dlc data...')
    
    % find csv files from DLC 
    if ~isempty(dir([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',info.animal,'_',info.date,'_*.csv']))
        input_files = dir([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',info.animal,'_',info.date,'_*.csv']);
    elseif ~isempty(dir([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',info.animal,'-',info.date,'-*.csv']))
        input_files = dir([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',info.animal,'-',info.date,'-*.csv']);
    else
        disp(['Did not find any DLC csv files in ',[path.root_data,'DLCout\']])
    end
    
    % load csv files one-by-one
    numFiles = length(input_files);
    dlc_all = cell(3,numFiles);
    dlcHeader = {'fileName',...
        'nose','nose_p',...
        'pupil','pupil_1_p','pupil_2_p','pupil_3_p','pupil_4_p','pupil_5_p','pupil_6_p','pupil_7_p','pupil_8_p'};
    for n=1:numFiles
        dlc_all{1,n} = input_files(n).name;
        if n==1
            csvHeader = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','text',[1,3],[]); csvHeader = csvHeader.text;
        end
        
        demo = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','mixed',[1,15],[1,10]); demo = demo.mixed;
        
        test = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],1); test = test.num;
        
        % load nose tip data
        this_nose_tip_likelihood = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'nose_tip')+strcmp(csvHeader,'likelihood'),1)==2)); this_nose_tip_likelihood = this_nose_tip_likelihood.num;
        this_nose_tip_x = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'nose_tip')+strcmp(csvHeader,'x'),1)==2)); this_nose_tip_x = this_nose_tip_x.num;
        this_nose_tip_y = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'nose_tip')+strcmp(csvHeader,'y'),1)==2)); this_nose_tip_y = this_nose_tip_y.num;
        
        % storing nose metrics
        dlc_all{2,n} = sqrt((this_nose_tip_x.^2)+(this_nose_tip_y.^2));
        dlc_all{3,n} = this_nose_tip_likelihood;
        
        % load pupil data
        this_pupil_top_likelihood = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_top')+strcmp(csvHeader,'likelihood'),1)==2)); this_pupil_top_likelihood = this_pupil_top_likelihood.num;
        this_pupil_topright_likelihood = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_top_right')+strcmp(csvHeader,'likelihood'),1)==2)); this_pupil_topright_likelihood = this_pupil_topright_likelihood.num;
        this_pupil_right_likelihood = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_right')+strcmp(csvHeader,'likelihood'),1)==2)); this_pupil_right_likelihood = this_pupil_right_likelihood.num;
        this_pupil_bottomright_likelihood = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_bottom_right')+strcmp(csvHeader,'likelihood'),1)==2)); this_pupil_bottomright_likelihood = this_pupil_bottomright_likelihood.num;
        this_pupil_bottom_likelihood = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_bottom')+strcmp(csvHeader,'likelihood'),1)==2)); this_pupil_bottom_likelihood = this_pupil_bottom_likelihood.num;
        this_pupil_bottomleft_likelihood = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_bottom_left')+strcmp(csvHeader,'likelihood'),1)==2)); this_pupil_bottomleft_likelihood = this_pupil_bottomleft_likelihood.num;
        this_pupil_left_likelihood = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_left')+strcmp(csvHeader,'likelihood'),1)==2)); this_pupil_left_likelihood = this_pupil_left_likelihood.num;
        this_pupil_topleft_likelihood = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_top_left')+strcmp(csvHeader,'likelihood'),1)==2)); this_pupil_topleft_likelihood = this_pupil_topleft_likelihood.num;
        this_pupil_top_x = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_top')+strcmp(csvHeader,'x'),1)==2)); this_pupil_top_x = this_pupil_top_x.num;
        this_pupil_topright_x = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_top_right')+strcmp(csvHeader,'x'),1)==2)); this_pupil_topright_x = this_pupil_topright_x.num;
        this_pupil_right_x = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_right')+strcmp(csvHeader,'x'),1)==2)); this_pupil_right_x = this_pupil_right_x.num;
        this_pupil_bottomright_x = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_bottom_right')+strcmp(csvHeader,'x'),1)==2)); this_pupil_bottomright_x = this_pupil_bottomright_x.num;
        this_pupil_bottom_x = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_bottom')+strcmp(csvHeader,'x'),1)==2)); this_pupil_bottom_x = this_pupil_bottom_x.num;
        this_pupil_bottomleft_x = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_bottom_left')+strcmp(csvHeader,'x'),1)==2)); this_pupil_bottomleft_x = this_pupil_bottomleft_x.num;
        this_pupil_left_x = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_left')+strcmp(csvHeader,'x'),1)==2)); this_pupil_left_x = this_pupil_left_x.num;
        this_pupil_topleft_x = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_top_left')+strcmp(csvHeader,'x'),1)==2)); this_pupil_topleft_x = this_pupil_topleft_x.num;                        
        this_pupil_top_y = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_top')+strcmp(csvHeader,'y'),1)==2)); this_pupil_top_y = this_pupil_top_y.num;
        this_pupil_topright_y = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_top_right')+strcmp(csvHeader,'y'),1)==2)); this_pupil_topright_y = this_pupil_topright_y.num;
        this_pupil_right_y = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_right')+strcmp(csvHeader,'y'),1)==2)); this_pupil_right_y = this_pupil_right_y.num;
        this_pupil_bottomright_y = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_bottom_right')+strcmp(csvHeader,'y'),1)==2)); this_pupil_bottomright_y = this_pupil_bottomright_y.num;
        this_pupil_bottom_y = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_bottom')+strcmp(csvHeader,'y'),1)==2)); this_pupil_bottom_y = this_pupil_bottom_y.num;
        this_pupil_bottomleft_y = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_bottom_left')+strcmp(csvHeader,'y'),1)==2)); this_pupil_bottomleft_y = this_pupil_bottomleft_y.num;
        this_pupil_left_y = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_left')+strcmp(csvHeader,'y'),1)==2)); this_pupil_left_y = this_pupil_left_y.num;
        this_pupil_topleft_y = delimread([path.root_data,'DLCout\',info.animal,'\',info.animal,'_',info.date,'\',input_files(n).name],...
            ',','num',[],find(sum(strcmp(csvHeader,'pupil_top_left')+strcmp(csvHeader,'y'),1)==2)); this_pupil_topleft_y = this_pupil_topleft_y.num;                
        
        % calculate pupil area by fitting ellipse
        this_pupil_area = nan(size(this_pupil_top_likelihood));
        for i=1:length(this_pupil_top_likelihood)
            this_x = [this_pupil_top_x(i);this_pupil_topright_x(i);this_pupil_right_x(i);this_pupil_bottomright_x(i);this_pupil_bottom_x(i);this_pupil_bottomleft_x(i);this_pupil_left_x(i);this_pupil_topleft_x(i)];
            this_y = [this_pupil_top_y(i);this_pupil_topright_y(i);this_pupil_right_y(i);this_pupil_bottomright_y(i);this_pupil_bottom_y(i);this_pupil_bottomleft_y(i);this_pupil_left_y(i);this_pupil_topleft_y(i)];
            this_pupil_struct = fit_ellipse(this_x,this_y);
            try
                this_pupil_area(i) = pi*this_pupil_struct.long_axis*this_pupil_struct.short_axis;
            catch
            end
        end
        
        % storing nose metrics
        dlc_all{4,n} = this_pupil_area;
        dlc_all{5,n} = this_pupil_top_likelihood;
        dlc_all{6,n} = this_pupil_topright_likelihood;
        dlc_all{7,n} = this_pupil_right_likelihood;
        dlc_all{8,n} = this_pupil_bottomright_likelihood;
        dlc_all{9,n} = this_pupil_bottom_likelihood;
        dlc_all{10,n} = this_pupil_bottomleft_likelihood;
        dlc_all{11,n} = this_pupil_left_likelihood;
        dlc_all{12,n} = this_pupil_topleft_likelihood;

        disp(['--- Loaded file ',num2str(n),' / ',num2str(numFiles)])
    end
    
end


%% Prepare to bring camera data into pre, beh, post structure

% check if mot and dlc data match up
if ops.cam.do_mot && ops.cam.do_dlc
    if size(mot_all,2)~=size(dlc_all,2)
        warning('Different number of video files in mot_all and dlc_all structs.')
    else
        for n=1:size(mot_all,2)
            if length(mot_all{2,n})~=length(dlc_all{2,n})
                warning(['Video ',num2str(n),' has different number of frames in mot_all and dlc_all structs.'])
            end
        end
    end
end

% store file headers
if ops.cam.do_dlc
    cam_pre.csvHeader = csvHeader;
    cam_pre.dlcHeader = dlcHeader;
    cam_post.csvHeader = csvHeader;
    cam_post.dlcHeader = dlcHeader;
    cam_beh.csvHeader = csvHeader;
    cam_beh.dlcHeader = dlcHeader;
end

% compare frame numbers between s2p and cam files
if ops.cam.do_mot
    ref_all = dlc_all;
elseif ops.cam.do_dlc
    ref_all = mot_all;
end
cam_pre.videoIndices = find(contains({ref_all{1,:}},'_pre'));
cam_post.videoIndices = find(contains({ref_all{1,:}},'_post'));
cam_beh.videoIndices = find(contains({ref_all{1,:}},'_beh'));
temp = cellfun(@length,{ref_all{2,:}});
cam_pre.camFrames = nansum(temp(cam_pre.videoIndices));
cam_post.camFrames = nansum(temp(cam_post.videoIndices));
cam_beh.camFrames = nansum(temp(cam_beh.videoIndices));
cam_pre.s2pFrames = s2pFrames_pre;
cam_post.s2pFrames = s2pFrames_post;
cam_beh.s2pFrames = s2pFrames_beh;
if cam_pre.s2pFrames~=cam_pre.camFrames
    warning('Mismatch in number of pre frames between s2p and cam files.')
    if abs(cam_pre.s2pFrames-cam_pre.camFrames)<=5
        warning('They differ by only max 5 frames.')
    end
end
if cam_post.s2pFrames~=cam_post.camFrames
    warning('Mismatch in number of post frames between s2p and cam files.')
    if abs(cam_post.s2pFrames-cam_post.camFrames)<=5
        warning('They differ by only max 5 frames.')
    end
end
if cam_beh.s2pFrames~=cam_beh.camFrames
    warning('Mismatch in number of beh frames between s2p and cam files.')
    if abs(cam_beh.s2pFrames-cam_beh.camFrames)<=5
        warning('They differ by only max 5 frames.')
    end
    if abs(cam_beh.s2pFrames-nansum(temp(find(contains({ref_all{1,:}},'_beh_')))))<=5
        cam_beh.videoIndices = find(contains({ref_all{1,:}},'_beh_'));
        cam_beh.camFrames = nansum(temp(cam_beh.videoIndices));
        warning('Identified the cause for this mismatch and dealt with it.')
    elseif abs(cam_beh.s2pFrames-nansum(temp(find(contains({ref_all{1,:}},'_beh1_')))))<=5
        cam_beh.videoIndices = find(contains({ref_all{1,:}},'_beh1_'));
        cam_beh.camFrames = nansum(temp(cam_beh.videoIndices));
        warning('Identified the cause for this mismatch and dealt with it.')
    elseif abs(cam_beh.s2pFrames-nansum(temp(find(contains({ref_all{1,:}},'_beh2_')))))<=5
        cam_beh.videoIndices = find(contains({ref_all{1,:}},'_beh2_'));
        cam_beh.camFrames = nansum(temp(cam_beh.videoIndices));
        warning('Identified the cause for this mismatch and dealt with it.')
    elseif abs(cam_beh.s2pFrames-nansum(temp(find(contains({ref_all{1,:}},'_beh3_')))))<=5
        cam_beh.videoIndices = find(contains({ref_all{1,:}},'_beh3_'));
        cam_beh.camFrames = nansum(temp(cam_beh.videoIndices));
        warning('Identified the cause for this mismatch and dealt with it.')
    else
        warning('Could not automatically identify what the problem is.')
    end
end


%% Bring mot data into pre, beh, post structure

if ops.cam.do_mot
    
    % pre
    this_cumsum = 0;
    cam_pre.mot = nan(1,s2pFrames_pre);
    for n=1:length(cam_pre.videoIndices)
        this_videoIndex = cam_pre.videoIndices(n);
        this_content = mot_all{2,this_videoIndex};
        cam_pre.mot(this_cumsum+1:this_cumsum+length(this_content)) = this_content;
        this_cumsum = this_cumsum + length(this_content);
    end
    cam_pre.mot = nanzscore(cam_pre.mot);
    
    % post
    this_cumsum = 0;
    cam_post.mot = nan(1,s2pFrames_post);
    for n=1:length(cam_post.videoIndices)
        this_videoIndex = cam_post.videoIndices(n);
        this_content = mot_all{2,this_videoIndex};
        cam_post.mot(this_cumsum+1:this_cumsum+length(this_content)) = this_content;
        this_cumsum = this_cumsum + length(this_content);
    end
    cam_post.mot = nanzscore(cam_post.mot);
    
    % beh
    this_cumsum = 0;
    cam_beh.mot = nan(1,s2pFrames_beh);
    for n=1:length(cam_beh.videoIndices)
        this_videoIndex = cam_beh.videoIndices(n);
        this_content = mot_all{2,this_videoIndex};
        cam_beh.mot(this_cumsum+1:this_cumsum+length(this_content)) = this_content;
        this_cumsum = this_cumsum + length(this_content);
    end
    cam_beh.mot = nanzscore(cam_beh.mot);
end


%% Bring dlc data into pre, beh, post structure

if ops.cam.do_dlc
    
    % pre
    this_cumsum = 0;
    cam_pre.nose = nan(1,s2pFrames_pre);
    cam_pre.nose_P = nan(1,s2pFrames_pre);
    cam_pre.pupil = nan(1,s2pFrames_pre);
    cam_pre.pupil_p = nan(8,s2pFrames_pre);
    for n=1:length(cam_pre.videoIndices)
        this_videoIndex = cam_pre.videoIndices(n);
        
        this_content = dlc_all{find(strcmp(dlcHeader,'nose')),this_videoIndex};
        this_contentLength = length(this_content);
        cam_pre.nose(this_cumsum+1:this_cumsum+this_contentLength) = this_content;
        
        this_content = dlc_all{find(strcmp(dlcHeader,'nose_p')),this_videoIndex};
        cam_pre.nose_p(this_cumsum+1:this_cumsum+this_contentLength) = this_content;
        
        this_content = dlc_all{find(strcmp(dlcHeader,'pupil')),this_videoIndex};
        cam_pre.pupil(this_cumsum+1:this_cumsum+this_contentLength) = this_content;
        
        this_content = [];
        for k=1:8
            this_content = [this_content, dlc_all{find(strcmp(dlcHeader,['pupil_',num2str(k),'_p'])),this_videoIndex}];
        end
        cam_pre.pupil_p(:,this_cumsum+1:this_cumsum+this_contentLength) = this_content';
        
        this_cumsum = this_cumsum + this_contentLength;
    end
    cam_pre.sniffing = bandpass(cam_pre.nose,p.cam.sniffing_bp,info.scope.frameRate);
    cam_pre.sniffing(1:200)=NaN;
    cam_pre.sniffing(end-200+1:end)=NaN;
    
    % post
    this_cumsum = 0;
    cam_post.nose = nan(1,s2pFrames_post);
    cam_post.nose_P = nan(1,s2pFrames_post);
    cam_post.pupil = nan(1,s2pFrames_post);
    cam_post.pupil_p = nan(8,s2pFrames_post);
    for n=1:length(cam_post.videoIndices)
        this_videoIndex = cam_post.videoIndices(n);
        
        this_content = dlc_all{find(strcmp(dlcHeader,'nose')),this_videoIndex};
        this_contentLength = length(this_content);
        cam_post.nose(this_cumsum+1:this_cumsum+this_contentLength) = this_content;
        
        this_content = dlc_all{find(strcmp(dlcHeader,'nose_p')),this_videoIndex};
        cam_post.nose_p(this_cumsum+1:this_cumsum+this_contentLength) = this_content;
        
        this_content = dlc_all{find(strcmp(dlcHeader,'pupil')),this_videoIndex};
        cam_post.pupil(this_cumsum+1:this_cumsum+this_contentLength) = this_content;
        
        this_content = [];
        for k=1:8
            this_content = [this_content, dlc_all{find(strcmp(dlcHeader,['pupil_',num2str(k),'_p'])),this_videoIndex}];
        end
        cam_post.pupil_p(:,this_cumsum+1:this_cumsum+this_contentLength) = this_content';
        
        this_cumsum = this_cumsum + this_contentLength;
    end
    cam_post.sniffing = bandpass(cam_post.nose,p.cam.sniffing_bp,info.scope.frameRate);
    cam_post.sniffing(1:200)=NaN;
    cam_post.sniffing(end-200+1:end)=NaN;
    
    % beh
    this_cumsum = 0;
    cam_beh.nose = nan(1,s2pFrames_beh);
    cam_beh.nose_P = nan(1,s2pFrames_beh);
    cam_beh.pupil = nan(1,s2pFrames_beh);
    cam_beh.pupil_p = nan(8,s2pFrames_beh);
    for n=1:length(cam_beh.videoIndices)
        this_videoIndex = cam_beh.videoIndices(n);
        
        this_content = dlc_all{find(strcmp(dlcHeader,'nose')),this_videoIndex};
        this_contentLength = length(this_content);
        cam_beh.nose(this_cumsum+1:this_cumsum+this_contentLength) = this_content;
        
        this_content = dlc_all{find(strcmp(dlcHeader,'nose_p')),this_videoIndex};
        cam_beh.nose_p(this_cumsum+1:this_cumsum+this_contentLength) = this_content;
        
        this_content = dlc_all{find(strcmp(dlcHeader,'pupil')),this_videoIndex};
        cam_beh.pupil(this_cumsum+1:this_cumsum+this_contentLength) = this_content;
        
        this_content = [];
        for k=1:8
            this_content = [this_content, dlc_all{find(strcmp(dlcHeader,['pupil_',num2str(k),'_p'])),this_videoIndex}];
        end
        cam_beh.pupil_p(:,this_cumsum+1:this_cumsum+this_contentLength) = this_content';
        
        this_cumsum = this_cumsum + this_contentLength;
    end
    cam_beh.sniffing = nan(size(cam_beh.nose));
    temp = rmmissing(cam_beh.nose);
    cam_beh.sniffing(1:length(temp)) = bandpass(temp,p.cam.sniffing_bp,info.scope.frameRate);
    cam_beh.sniffing(1:200)=NaN;
    cam_beh.sniffing(end-200+1:end)=NaN;
end


%% Save and return

save([path.filepart_outX,'cam_pre.mat'],'cam_pre','-v7.3');
disp(['--- Added cam_pre file to repoX as ',[path.filepart_outX,'cam_pre.mat'],'.'])
save([path.filepart_outX,'cam_beh.mat'],'cam_beh','-v7.3');
disp(['--- Added cam_beh file to repoX as ',[path.filepart_outX,'cam_beh.mat'],'.'])
save([path.filepart_outX,'cam_post.mat'],'cam_post','-v7.3');
disp(['--- Added cam_post file to repoX as ',[path.filepart_outX,'cam_post.mat'],'.'])


end





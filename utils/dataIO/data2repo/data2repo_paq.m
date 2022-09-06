function [path,paq_beh] = data2repo_paq(info,p,path,raw,task,thor_beh)

if nargin < 6
    skipThorChecks = true;
else
    skipThorChecks = false;
end
if nargin < 5
    skipTaskChecks = true;
else
    skipTaskChecks = false;
end
if any(info.data.numSubFragments_paq>1)
    error('Multiple sub-fragments in paq files not yet implemented.')
end

%% Load and process paq data

% Loading ThorSync data and detecting events
disp('--- Loading paq data and detecting events.')
for e=1:length(info.epochs)
    if e~=2 || info.data.numFragments==1
        
        if (e==1 && info.data.missingPaqPreFile) || (e==3 && info.data.missingPaqPostFile)
        else

            % paq data -> paq
            try 
                path.(['file_in_paq_',info.epochs{e}]) = [path.folder_data,'Behaviour\',info.animal,'_',info.date,'_',info.epochs{e},'.paq'];
                paq.(info.epochs{e}) = paq2lab(path.(['file_in_paq_',info.epochs{e}])); 
            catch
                path.(['file_in_paq_',info.epochs{e}]) = [path.folder_data,'Behaviour\',info.animal,'-',info.date,'-',info.epochs{e},'.paq'];
                paq.(info.epochs{e}) = paq2lab(path.(['file_in_paq_',info.epochs{e}])); 
            end
            disp(['--- paq ',info.epochs{e},' file loaded.'])

            % extract paq data
            ts.(info.epochs{e}).cue = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.cue),'above',p.paq.voltageThreshold);
            ts.(info.epochs{e}).finalvalve_start = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.finalvalve),'above',p.paq.voltageThreshold);
            ts.(info.epochs{e}).finalvalve_stop = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.finalvalve),'below',p.paq.voltageThreshold)-1;
            ts.(info.epochs{e}).imaging = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.imaging),'above',p.paq.voltageThreshold);
            ts.(info.epochs{e}).lick_start = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.lick),'above',p.paq.voltageThreshold);
            ts.(info.epochs{e}).lick_stop = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.lick),'below',p.paq.voltageThreshold)-1;
            ts.(info.epochs{e}).responsewindow_start = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.responsewindow),'above',p.paq.voltageThreshold);
            ts.(info.epochs{e}).responsewindow_stop = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.responsewindow),'below',p.paq.voltageThreshold)-1;
            ts.(info.epochs{e}).reward_start = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.reward),'above',p.paq.voltageThreshold);
            ts.(info.epochs{e}).reward_stop = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.reward),'below',p.paq.voltageThreshold)-1;

            % handle Grid data without imaging
            if isempty(raw)
                if length(ts.(info.epochs{e}).imaging)>10000
                    warning('Thought its a Grid session but somehow there seem to be imaging frame triggers?');
                end
                ts.(info.epochs{e}).imaging = (1:round(info.paq.samplingRate/info.scope.frameRate):size(paq.(info.epochs{e}),1))';
            end

            % rotary encoder     
            RotA = paq.(info.epochs{e})(:,info.paq.config.rotA) > p.paq.voltageThreshold;
            RotB = paq.(info.epochs{e})(:,info.paq.config.rotB) > p.paq.voltageThreshold;
            temp = find([0;diff(RotA)]); % detect change in A (i.e. pulses)
            temp2 = double(RotA(temp)==RotB(temp)); % check if RotA and RotB the same or not at pulses
            temp2(temp2==0) = -1;
            position_pulses = zeros(size(RotA));
            position_pulses(temp) = temp2;
            position_pulses = cumsum(position_pulses);
            position = info.paq.rot.wheelCircumference * position_pulses/info.paq.rot.pulsesPerRevolution; % [cm]
            ts.(info.epochs{e}).position = position(ts.(info.epochs{e}).imaging)'; % [cm]
            ts.(info.epochs{e}).speed = [0,diff(ts.(info.epochs{e}).position)] * info.scope.frameRate; % [cm/s]
            ts.(info.epochs{e}).acceleration = [0,diff(ts.(info.epochs{e}).speed)] * info.scope.frameRate; % [cm/s^2]
        end
       
    else
        for i=1:info.data.numFragments
            
            % paq data -> paq
            try
                if info.data.numSubFragments_paq(i)>1
                    path.(['file_in_paq_',info.epochs{e},'_',num2str(i)]) = [path.folder_data,'Behaviour\',info.animal,'_',info.date,'_',info.epochs{e},'-',num2str(i),'-1.paq'];
                    warning('Multiple sub-fragments in paq files not yet implemented. Only using the first sub-fragment for now.')
                else
                    path.(['file_in_paq_',info.epochs{e},'_',num2str(i)]) = [path.folder_data,'Behaviour\',info.animal,'_',info.date,'_',info.epochs{e},'-',num2str(i),'.paq'];
                end
                paq.(info.epochs{e}) = paq2lab(path.(['file_in_paq_',info.epochs{e},'_',num2str(i)]));
            catch
                if info.data.numSubFragments_paq(i)>1
                    path.(['file_in_paq_',info.epochs{e},'_',num2str(i)]) = [path.folder_data,'Behaviour\',info.animal,'-',info.date,'-',info.epochs{e},'-',num2str(i),'-1.paq'];
                    warning('Multiple sub-fragments in paq files not yet implemented. Only using the first sub-fragment for now.')
                else
                    path.(['file_in_paq_',info.epochs{e},'_',num2str(i)]) = [path.folder_data,'Behaviour\',info.animal,'-',info.date,'-',info.epochs{e},'-',num2str(i),'.paq'];
                end
                paq.(info.epochs{e}) = paq2lab(path.(['file_in_paq_',info.epochs{e},'_',num2str(i)]));
            end
            disp(['--- paq ',info.epochs{e},'-',num2str(i),' file loaded.'])

            % extract paq data
            ts.([info.epochs{e},'_',num2str(i)]).cue = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.cue),'above',p.paq.voltageThreshold);
            ts.([info.epochs{e},'_',num2str(i)]).finalvalve_start = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.finalvalve),'above',p.paq.voltageThreshold);
            ts.([info.epochs{e},'_',num2str(i)]).finalvalve_stop = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.finalvalve),'below',p.paq.voltageThreshold)-1;
            ts.([info.epochs{e},'_',num2str(i)]).imaging = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.imaging),'above',p.paq.voltageThreshold);
            ts.([info.epochs{e},'_',num2str(i)]).lick_start = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.lick),'above',p.paq.voltageThreshold);
            ts.([info.epochs{e},'_',num2str(i)]).lick_stop = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.lick),'below',p.paq.voltageThreshold)-1;
            ts.([info.epochs{e},'_',num2str(i)]).responsewindow_start = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.responsewindow),'above',p.paq.voltageThreshold);
            ts.([info.epochs{e},'_',num2str(i)]).responsewindow_stop = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.responsewindow),'below',p.paq.voltageThreshold)-1;
            ts.([info.epochs{e},'_',num2str(i)]).reward_start = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.reward),'above',p.paq.voltageThreshold);
            ts.([info.epochs{e},'_',num2str(i)]).reward_stop = detectThresholdCrossing(paq.(info.epochs{e})(:,info.paq.config.reward),'below',p.paq.voltageThreshold)-1;

            % handle Grid data without imaging
            if isempty(raw)
                if length(ts.([info.epochs{e},'_',num2str(i)]).imaging)>10000
                    warning('Thought its a Grid session but somehow there seem to be imaging frame triggers?');
                end
                ts.([info.epochs{e},'_',num2str(i)]).imaging = (1:round(info.paq.samplingRate/info.scope.frameRate):size(paq.(info.epochs{e}),1))';
            end
        
            % rotary encoder     
            RotA = paq.(info.epochs{e})(:,info.paq.config.rotA) > p.paq.voltageThreshold;
            RotB = paq.(info.epochs{e})(:,info.paq.config.rotB) > p.paq.voltageThreshold;
            temp = find([0;diff(RotA)]); % detect change in A (i.e. pulses)
            temp2 = double(RotA(temp)==RotB(temp)); % check if RotA and RotB the same or not at pulses
            temp2(temp2==0) = -1;
            position_pulses = zeros(size(RotA));
            position_pulses(temp) = temp2;
            position_pulses = cumsum(position_pulses);
            position = info.paq.rot.wheelCircumference * position_pulses/info.paq.rot.pulsesPerRevolution; % [cm]
            ts.([info.epochs{e},'_',num2str(i)]).position = position(ts.([info.epochs{e},'_',num2str(i)]).imaging)'; % [cm]
            ts.([info.epochs{e},'_',num2str(i)]).speed = [0,diff(ts.([info.epochs{e},'_',num2str(i)]).position)] * info.scope.frameRate; % [cm/s]
            ts.([info.epochs{e},'_',num2str(i)]).acceleration = [0,diff(ts.([info.epochs{e},'_',num2str(i)]).speed)] * info.scope.frameRate; % [cm/s^2]
        end
    end
end
clear paq; 


%% Deal with fragmentations and additional frame triggers

% Dealing with fragmentations and additional frame triggers
disp('--- Dealing with fragmentations and additional frame triggers.')

e=2;
paq_beh.numFrameTriggers_original = [];
paq_beh.numFrameTriggers_postDropping = [];
if info.data.numFragments~=1
    ts.beh.cue = [];
    ts.beh.finalvalve_start = [];
    ts.beh.finalvalve_stop = [];
    ts.beh.imaging = [];
    ts.beh.lick_start = [];
    ts.beh.lick_stop = [];
    ts.beh.responsewindow_start = [];
    ts.beh.responsewindow_stop = [];
    ts.beh.reward_start = [];
    ts.beh.reward_stop = [];
    ts.beh.position = [];
    ts.beh.speed = [];
    ts.beh.acceleration = [];
    temp2 = 0;
    temp3 = 0;
    for i=1:info.data.numFragments
        paq_beh.numFrameTriggers_original(i) = length(ts.([info.epochs{e},'_',num2str(i)]).imaging);
        if ~isempty(raw)
            paq_beh.numFrameTriggers_postDropping(i) = raw.([info.epochs{e},'_',num2str(i)]).frames;
        else
            paq_beh.numFrameTriggers_postDropping(i) = paq_beh.numFrameTriggers_original(i);
        end
        if info.data.numSubFragments_paq(i)>1
            paq_beh.numFrameTriggers_postDropping(i) = paq_beh.numFrameTriggers_original(i);
            warning('Multiple sub-fragments in paq files not yet implemented. This line has to be deleted after implementation.')
        end
        
        try
            ts.beh.imaging = [ts.beh.imaging; temp2+ ts.([info.epochs{e},'_',num2str(i)]).imaging(1:paq_beh.numFrameTriggers_postDropping(i))];
            ts.beh.position = [ts.beh.position; temp3+ ts.([info.epochs{e},'_',num2str(i)]).position(1:paq_beh.numFrameTriggers_postDropping(i))'];
            ts.beh.speed = [ts.beh.speed; ts.([info.epochs{e},'_',num2str(i)]).speed(1:paq_beh.numFrameTriggers_postDropping(i))'];
            ts.beh.acceleration = [ts.beh.acceleration; ts.([info.epochs{e},'_',num2str(i)]).acceleration(1:paq_beh.numFrameTriggers_postDropping(i))'];
        catch % deals with the issue that the paq recording somehow crashed before the imaging / ThorSync could be stopped
            temp4 = nan(paq_beh.numFrameTriggers_postDropping(i),1);
            temp4(1:length(ts.([info.epochs{e},'_',num2str(i)]).imaging)) = ts.([info.epochs{e},'_',num2str(i)]).imaging;
            ts.beh.imaging = [ts.beh.imaging; temp2+ temp4];            
            temp4 = nan(paq_beh.numFrameTriggers_postDropping(i),1);
            temp4(1:length(ts.([info.epochs{e},'_',num2str(i)]).position)) = ts.([info.epochs{e},'_',num2str(i)]).position';
            ts.beh.position = [ts.beh.position; temp2+ temp4];           
            temp4 = nan(paq_beh.numFrameTriggers_postDropping(i),1);
            temp4(1:length(ts.([info.epochs{e},'_',num2str(i)]).speed)) = ts.([info.epochs{e},'_',num2str(i)]).speed';
            ts.beh.speed = [ts.beh.speed; temp2+ temp4];
            temp4 = nan(paq_beh.numFrameTriggers_postDropping(i),1);
            temp4(1:length(ts.([info.epochs{e},'_',num2str(i)]).acceleration)) = ts.([info.epochs{e},'_',num2str(i)]).acceleration';
            ts.beh.acceleration = [ts.beh.acceleration; temp2+ temp4];
            warning('A paq file ends before ThorSync and imaging end. Dealt with it.')
        end
        
        if info.data.numSubFragments_paq(i)>1
            ts.beh.cue = [ts.beh.cue; temp2+ ts.([info.epochs{e},'_',num2str(i)]).cue(1:end)];
            warning('Multiple sub-fragments in paq files not yet implemented. This line has to be deleted after implementation.')
        else
            ts.beh.cue = [ts.beh.cue; temp2+ ts.([info.epochs{e},'_',num2str(i)]).cue(1:info.data.fragments_numTrials(i))];
        end
        if length(ts.([info.epochs{e},'_',num2str(i)]).cue) > info.data.fragments_numTrials(i) 
            temp = ts.([info.epochs{e},'_',num2str(i)]).cue(info.data.fragments_numTrials(i)+1) - info.scope.samplingRate; % last frame for valid events
            ts.beh.finalvalve_start = [ts.beh.finalvalve_start; temp2+ ts.([info.epochs{e},'_',num2str(i)]).finalvalve_start(ts.([info.epochs{e},'_',num2str(i)]).finalvalve_start<=temp)];
            ts.beh.finalvalve_stop = [ts.beh.finalvalve_stop; temp2+ ts.([info.epochs{e},'_',num2str(i)]).finalvalve_stop(ts.([info.epochs{e},'_',num2str(i)]).finalvalve_stop<=temp)];
            ts.beh.lick_start = [ts.beh.lick_start; temp2+ ts.([info.epochs{e},'_',num2str(i)]).lick_start(ts.([info.epochs{e},'_',num2str(i)]).lick_start<=temp)];
            ts.beh.lick_stop = [ts.beh.lick_stop; temp2+ ts.([info.epochs{e},'_',num2str(i)]).lick_stop(ts.([info.epochs{e},'_',num2str(i)]).lick_stop<=temp)];
            ts.beh.responsewindow_start = [ts.beh.responsewindow_start; temp2+ ts.([info.epochs{e},'_',num2str(i)]).responsewindow_start(ts.([info.epochs{e},'_',num2str(i)]).responsewindow_start<=temp)];
            ts.beh.responsewindow_stop = [ts.beh.responsewindow_stop; temp2+ ts.([info.epochs{e},'_',num2str(i)]).responsewindow_stop(ts.([info.epochs{e},'_',num2str(i)]).responsewindow_stop<=temp)];
            ts.beh.reward_start = [ts.beh.reward_start; temp2+ ts.([info.epochs{e},'_',num2str(i)]).reward_start(ts.([info.epochs{e},'_',num2str(i)]).reward_start<=temp)];
            ts.beh.reward_stop = [ts.beh.reward_stop; temp2+ ts.([info.epochs{e},'_',num2str(i)]).reward_stop(ts.([info.epochs{e},'_',num2str(i)]).reward_stop<=temp)];
        else
            ts.beh.finalvalve_start = [ts.beh.finalvalve_start; temp2+ ts.([info.epochs{e},'_',num2str(i)]).finalvalve_start];
            ts.beh.finalvalve_stop = [ts.beh.finalvalve_stop; temp2+ ts.([info.epochs{e},'_',num2str(i)]).finalvalve_stop];
            ts.beh.lick_start = [ts.beh.lick_start; temp2+ ts.([info.epochs{e},'_',num2str(i)]).lick_start];
            ts.beh.lick_stop = [ts.beh.lick_stop; temp2+ ts.([info.epochs{e},'_',num2str(i)]).lick_stop];
            ts.beh.responsewindow_start = [ts.beh.responsewindow_start; temp2+ ts.([info.epochs{e},'_',num2str(i)]).responsewindow_start];
            ts.beh.responsewindow_stop = [ts.beh.responsewindow_stop; temp2+ ts.([info.epochs{e},'_',num2str(i)]).responsewindow_stop];
            ts.beh.reward_start = [ts.beh.reward_start; temp2+ ts.([info.epochs{e},'_',num2str(i)]).reward_start];
            ts.beh.reward_stop = [ts.beh.reward_stop; temp2+ ts.([info.epochs{e},'_',num2str(i)]).reward_stop];
        end
        temp2 = temp2 + ts.([info.epochs{e},'_',num2str(i)]).imaging(end);
        temp3 = temp3 + ts.([info.epochs{e},'_',num2str(i)]).position(end);
        
        if info.data.numSubFragments_paq(i)>1
            if i>1
                error('Multiple sub-fragments in non-first paq file are not yet implemented.')
            else
                warning('Multiple sub-fragments in paq files not yet implemented. This fraction has to be adjusted after implementation.')
                temp = ts.beh.cue; ts.beh.cue = nan(info.data.fragments_numTrials(i),1); ts.beh.cue(1:length(temp)) = temp;
                temp = ts.beh.finalvalve_start; ts.beh.finalvalve_start = nan(info.data.fragments_numTrials(i)*2,1); ts.beh.finalvalve_start(1:length(temp)) = temp;
                temp = ts.beh.finalvalve_stop; ts.beh.finalvalve_stop = nan(info.data.fragments_numTrials(i)*2,1); ts.beh.finalvalve_stop(1:length(temp)) = temp;
                temp = ts.beh.responsewindow_start; ts.beh.responsewindow_start = nan(info.data.fragments_numTrials(i),1); ts.beh.responsewindow_start(1:length(temp)) = temp;
                temp = ts.beh.responsewindow_stop; ts.beh.responsewindow_stop = nan(info.data.fragments_numTrials(i),1); ts.beh.responsewindow_stop(1:length(temp)) = temp;            
            end
        end
        
        % deals with the issue that the paq recording somehow crashed before the imaging / ThorSync could be stopped
        if length(ts.([info.epochs{e},'_',num2str(i)]).cue) > info.data.fragments_numTrials(i)
            upper = ts.([info.epochs{e},'_',num2str(i)]).cue(info.data.fragments_numTrials(i)+1);
            ts.([info.epochs{e},'_',num2str(i)]).cue = ts.([info.epochs{e},'_',num2str(i)]).cue(1:info.data.fragments_numTrials(i));
            ts.([info.epochs{e},'_',num2str(i)]).responsewindow_start = ts.([info.epochs{e},'_',num2str(i)]).responsewindow_start(1:info.data.fragments_numTrials(i));
            ts.([info.epochs{e},'_',num2str(i)]).responsewindow_stop = ts.([info.epochs{e},'_',num2str(i)]).responsewindow_stop(1:info.data.fragments_numTrials(i));
            ts.([info.epochs{e},'_',num2str(i)]).finalvalve_start = ts.([info.epochs{e},'_',num2str(i)]).finalvalve_start(1:2*info.data.fragments_numTrials(i));
            ts.([info.epochs{e},'_',num2str(i)]).finalvalve_stop = ts.([info.epochs{e},'_',num2str(i)]).finalvalve_stop(1:2*info.data.fragments_numTrials(i));
            ts.([info.epochs{e},'_',num2str(i)]).lick_start = ts.([info.epochs{e},'_',num2str(i)]).lick_start(ts.([info.epochs{e},'_',num2str(i)]).lick_start<upper);
            ts.([info.epochs{e},'_',num2str(i)]).lick_stop = ts.([info.epochs{e},'_',num2str(i)]).lick_stop(ts.([info.epochs{e},'_',num2str(i)]).lick_stop<upper);
            ts.([info.epochs{e},'_',num2str(i)]).reward_start = ts.([info.epochs{e},'_',num2str(i)]).reward_start(ts.([info.epochs{e},'_',num2str(i)]).reward_start<upper);
            ts.([info.epochs{e},'_',num2str(i)]).reward_stop = ts.([info.epochs{e},'_',num2str(i)]).reward_stop(ts.([info.epochs{e},'_',num2str(i)]).reward_stop<upper);
            warning('A paq file ends before task ends. Dealt with it.')
        end
    end
else
    paq_beh.numFrameTriggers_original = length(ts.beh.imaging);
    if ~isempty(raw)
        paq_beh.numFrameTriggers_postDropping = raw.beh.frames;
    else
        paq_beh.numFrameTriggers_postDropping = paq_beh.numFrameTriggers_original;
    end
    % MB20210630: added the catch to deal with this issue
    try
        ts.beh.imaging = ts.beh.imaging(1:paq_beh.numFrameTriggers_postDropping);
        ts.beh.position = ts.beh.position(1:paq_beh.numFrameTriggers_postDropping);
        ts.beh.speed = ts.beh.speed(1:paq_beh.numFrameTriggers_postDropping);
        ts.beh.acceleration = ts.beh.acceleration(1:paq_beh.numFrameTriggers_postDropping);
    catch
        warning('There are somehow less imaging frames in the paq file than in the scope files. Check manually.')
    end
end
if any(paq_beh.numFrameTriggers_original-paq_beh.numFrameTriggers_postDropping>2)
    warning('For at least one paq beh fragment more than 2 frame triggers were dropped.')
end


%% Identify frames with paq events in pre and post

disp('--- Aligning paq events to imaging frames.')

% identify frames with paq events in pre
if ~info.data.missingPaqPreFile
    e=1;
    temp = [-Inf, ts.(info.epochs{e}).imaging(2:end)', Inf]; % even if an event starts at the very end of a frame, it's still counted to the same frame
    paq_pre.licks = nan(2,max([length(discretize(ts.(info.epochs{e}).lick_start', temp)),length(discretize(ts.(info.epochs{e}).lick_stop', temp))]));
    temp2 = discretize(ts.(info.epochs{e}).lick_start', temp); % take care of signals that may be active before imaging starts
    temp3 = temp2((temp2~=1)&(temp2~=length(temp)-1));
    paq_pre.licks(1,1:length(temp3)) = temp3;
    temp2 = discretize(ts.(info.epochs{e}).lick_stop', temp); % take care of signals that may be active before imaging starts
    temp3 = temp2((temp2~=1)&(temp2~=length(temp)-1));
    paq_pre.licks(2,1:length(temp3)) = temp3;
    paq_pre.licks = paq_pre.licks(:,~all(isnan(paq_pre.licks)));
    paq_pre.position = ts.(info.epochs{e}).position;
    paq_pre.speed = ts.(info.epochs{e}).speed;
    paq_pre.acceleration = ts.(info.epochs{e}).acceleration;
else
    paq_pre.position = nan(1,info.scope.numFrames_pre);
    paq_pre.speed = nan(1,info.scope.numFrames_pre);
    paq_pre.acceleration = nan(1,info.scope.numFrames_pre);
    paq_pre.licks = nan(2,1);
    warning('Missing or corrupted paq_pre file. Assigned NaNs.')
end

% identify frames with paq events in post
if ~info.data.missingPaqPostFile
    e=3;
    temp = [-Inf, ts.(info.epochs{e}).imaging(2:end)', Inf]; % even if an event starts at the very end of a frame, it's still counted to the same frame
    paq_post.licks = nan(2,max([length(discretize(ts.(info.epochs{e}).lick_start', temp)),length(discretize(ts.(info.epochs{e}).lick_stop', temp))]));
    temp2 = discretize(ts.(info.epochs{e}).lick_start', temp); % take care of signals that may be active before imaging starts
    temp3 = temp2((temp2~=1)&(temp2~=length(temp)-1));
    paq_post.licks(1,1:length(temp3)) = temp3;
    temp2 = discretize(ts.(info.epochs{e}).lick_stop', temp); % take care of signals that may be active before imaging starts
    temp3 = temp2((temp2~=1)&(temp2~=length(temp)-1));
    paq_post.licks(2,1:length(temp3)) = temp3;
    paq_post.licks = paq_post.licks(:,~all(isnan(paq_post.licks)));
    paq_post.position = ts.(info.epochs{e}).position;
    paq_post.speed = ts.(info.epochs{e}).speed;
    paq_post.acceleration = ts.(info.epochs{e}).acceleration;
else
    paq_post.position = nan(1,info.scope.numFrames_post);
    paq_post.speed = nan(1,info.scope.numFrames_post);
    paq_post.acceleration = nan(1,info.scope.numFrames_post);
    paq_post.licks = nan(2,1);
    warning('Missing or corrupted paq_post file. Assigned NaNs.')
end

% SANITY CHECK: number of frames

if ~info.data.missingPaqPreFile
    if length(paq_pre.position)~=info.scope.numFrames_pre
        paq_pre.position = standardisePrePostSize(paq_pre.position,info.scope.numFrames_pre);
        paq_pre.speed = standardisePrePostSize(paq_pre.speed,info.scope.numFrames_pre);
        paq_pre.acceleration = standardisePrePostSize(paq_pre.acceleration,info.scope.numFrames_pre);
        warning(['Mismatch in number of pre frames between paq file and user input. Dealt with it.'])   
    end   
end
if ~info.data.missingPaqPostFile
    if length(paq_post.position)~=info.scope.numFrames_post
        paq_post.position = standardisePrePostSize(paq_post.position,info.scope.numFrames_post);
        paq_post.speed = standardisePrePostSize(paq_post.speed,info.scope.numFrames_post);
        paq_post.acceleration = standardisePrePostSize(paq_post.acceleration,info.scope.numFrames_post);
        warning(['Mismatch in number of post frames between paq file and user input. Dealt with it.'])   
    end
end


%% Identify frames with paq events in beh

e=2;
temp = [-Inf, ts.(info.epochs{e}).imaging(2:end)', Inf]; % even if an event starts at the very end of a frame, it's still counted to the same frame
if any(isnan(temp))
    temp2 = find(isnan(temp));
    temp(temp2) = linspace(temp(temp2(1)-1),temp(temp2(end)+1),length(temp2));
end

% Special case fixer
if isfield(info.data,'specialCaseFixer')
    if info.data.specialCaseFixer == 1
        ts.(info.epochs{e}).cue = ts.(info.epochs{e}).cue(9:end);
        ts.(info.epochs{e}).finalvalve_start = ts.(info.epochs{e}).finalvalve_start(17:end);
        ts.(info.epochs{e}).finalvalve_stop = ts.(info.epochs{e}).finalvalve_stop(17:end);
        ts.(info.epochs{e}).finalvalve_start = [ts.(info.epochs{e}).finalvalve_start;NaN;NaN];
        ts.(info.epochs{e}).finalvalve_stop = [ts.(info.epochs{e}).finalvalve_stop;NaN;NaN];
        ts.(info.epochs{e}).responsewindow_start = ts.(info.epochs{e}).responsewindow_start(9:end);
        ts.(info.epochs{e}).responsewindow_stop = ts.(info.epochs{e}).responsewindow_stop(9:end);
    end
end

paq_beh.finalvalve = nan(2,info.task.numTrials*2);
paq_beh.responsewindow = nan(2,max([length(discretize(ts.(info.epochs{e}).responsewindow_start', temp)),length(discretize(ts.(info.epochs{e}).responsewindow_stop', temp))]));
paq_beh.reward = nan(2,max([length(discretize(ts.(info.epochs{e}).reward_start', temp)),length(discretize(ts.(info.epochs{e}).reward_stop', temp))]));
paq_beh.licks = nan(2,max([length(discretize(ts.(info.epochs{e}).lick_start', temp)),length(discretize(ts.(info.epochs{e}).lick_stop', temp))]));

paq_beh.sync = discretize(ts.(info.epochs{e}).cue', temp);
temp3 = discretize(ts.(info.epochs{e}).finalvalve_start', temp);
paq_beh.finalvalve(1,1:length(temp3)) = temp3;
temp3 = discretize(ts.(info.epochs{e}).finalvalve_stop', temp);
paq_beh.finalvalve(2,1:length(temp3)) = temp3;
paq_beh.finalvalve = reshape(paq_beh.finalvalve,4,info.task.numTrials);
temp3 = discretize(ts.(info.epochs{e}).responsewindow_start', temp);
paq_beh.responsewindow(1,1:length(temp3)) = temp3;
temp3 = discretize(ts.(info.epochs{e}).responsewindow_stop', temp);
paq_beh.responsewindow(2,1:length(temp3)) = temp3;
temp3 = discretize(ts.(info.epochs{e}).reward_start', temp);
paq_beh.reward(1,1:length(temp3)) = temp3;
temp3 = discretize(ts.(info.epochs{e}).reward_stop', temp);
paq_beh.reward(2,1:length(temp3)) = temp3;
temp2 = discretize(ts.(info.epochs{e}).lick_start', temp); % take care of signals that may be active before imaging starts
temp3 = temp2((temp2~=1)&(temp2~=length(temp)-1));
paq_beh.licks(1,1:length(temp3)) = temp3;
temp2 = discretize(ts.(info.epochs{e}).lick_stop', temp); % take care of signals that may be active before imaging starts
temp3 = temp2((temp2~=1)&(temp2~=length(temp)-1));
paq_beh.licks(2,1:length(temp3)) = temp3;
paq_beh.licks = paq_beh.licks(:,~all(isnan(paq_beh.licks)));
paq_beh.position = ts.(info.epochs{e}).position;
paq_beh.speed = ts.(info.epochs{e}).speed;
paq_beh.acceleration = ts.(info.epochs{e}).acceleration;

if isfield(info.data,'prematureImagingStop_nanAfter')
    paq_beh.sync_all = paq_beh.sync;
    paq_beh.responsewindow_all = paq_beh.responsewindow;
    paq_beh.licks_all = paq_beh.licks;
    paq_beh.reward_all = paq_beh.reward;
    paq_beh.finalvalve_all = paq_beh.finalvalve;
    warning('Dealing with premature imaging stop for paq_beh file.')
end
    
% remove events that happen after imaging has stopped
if length(paq_beh.sync)~=length(unique(paq_beh.sync))
    warning('Removed events that happened after imaging has stopped from paq_beh file')
    temp2 = unique(paq_beh.sync);
    temp2 = temp2(1:end-1);
    paq_beh.sync = nan(1,info.task.numTrials);
    paq_beh.sync(1:length(temp2)) = temp2;
    temp2 = [unique(paq_beh.responsewindow(1,:));unique(paq_beh.responsewindow(2,:))];
    temp2 = temp2(:,1:end-1);
    paq_beh.responsewindow = nan(2,info.task.numTrials);
    paq_beh.responsewindow(:,1:length(temp2)) = temp2;
    try
        paq_beh.licks = [unique(paq_beh.licks(1,:));unique(paq_beh.licks(2,:))];
        paq_beh.licks = paq_beh.licks(:,1:end-1);
    catch
        warning('...but skipped licking data as that caused trouble.')
    end  
    paq_beh.reward = [unique(paq_beh.reward(1,:));unique(paq_beh.reward(2,:))];
    paq_beh.reward = paq_beh.reward(:,1:end-1);
    temp2 = [unique(paq_beh.finalvalve(1,:));unique(paq_beh.finalvalve(2,:));unique(paq_beh.finalvalve(3,:));unique(paq_beh.finalvalve(4,:))];
    temp2 = temp2(:,1:end-1);
    paq_beh.finalvalve = nan(4,info.task.numTrials);
    paq_beh.finalvalve(:,1:length(temp2)) = temp2;
end

% SANITY CHECKS number of frames
if ~skipThorChecks
    
    % number of frames
    if ~all(thor_beh.numFrameTriggers_postDropping == paq_beh.numFrameTriggers_postDropping)
        warning(['Mismatch in number of beh frames between paq file(s) and ThorSync file(s) - postDropping'])   
    end
    
    % sync happens at identical imaging frames (thor vs paq)
    if length(paq_beh.sync) ~= length(thor_beh.sync) % MB20210619: added extra SniffinSync handling
        if length(thor_beh.sync) > length(paq_beh.sync) 
            thor_beh.sync_raw = thor_beh.sync;
            thor_beh.sync = [];
            for i=1:length(thor_beh.sync_raw)
                if any(abs(thor_beh.sync_raw(i)-paq_beh.sync)<=1)
                    thor_beh.sync = [thor_beh.sync, thor_beh.sync_raw(i)];
                end
            end
            disp(['--- Deleted ',num2str(length(thor_beh.sync_raw)-length(paq_beh.sync)),' extra SniffinSyncs from thor_beh file.'])
            thor_beh = orderfields(thor_beh);
            save([path.filepart_out,'thor_beh.mat'],'thor_beh','-v7.3');
            disp(['--- Overwrote thor_beh file in repo at ',[path.filepart_out,'thor_beh.mat'],'.'])
        else
            warning(['Unequal number of SniffinSyncs between thor and paq data. But not extra SniffinSyncs.'])
            if length(thor_beh.sync) < length(paq_beh.sync)  % MB20211102: to deal with sessions like Stanage_20210925, where the sync of one thor_beh fragment is useless
                thor_beh.sync_raw = thor_beh.sync;
                thor_beh.sync = paq_beh.sync;
                thor_beh = orderfields(thor_beh);
                save([path.filepart_out,'thor_beh.mat'],'thor_beh','-v7.3');
                disp(['--- Overwrote thor_beh file in repo at ',[path.filepart_out,'thor_beh.mat'],'.'])
                warning(['Cloned thor_beh.sync from paq_beh.sync due to thor_beh.sync issues in a fragmented file.'])                
            end
        end
    end
    if length(paq_beh.sync) ~= length(thor_beh.sync)
        warning(['Somehow still unequal number of SniffinSyncs between thor and paq data.'])
    end
    if any(abs(paq_beh.sync-thor_beh.sync)>1)
        temp = paq_beh.sync-thor_beh.sync;
        if temp(1)==0 && temp(end)==0
            thor_beh.sync_raw = thor_beh.sync;
            thor_beh.sync = paq_beh.sync;
            thor_beh = orderfields(thor_beh);
            save([path.filepart_out,'thor_beh.mat'],'thor_beh','-v7.3');
            disp(['--- Overwrote thor_beh file in repo at ',[path.filepart_out,'thor_beh.mat'],'.'])
            disp(['--- Cloned thor_beh.sync from paq_beh.sync due to extra SniffinSyncs in a fragmented file.'])
        else
            warning(['SniffinSync do not perfectly correspond to the same imaging frames in thor vs paq data (differ by more than 1 frame). Possibly extra SniffinSyncs. No correction implemented yet.'])
        end  
    end
    
    % SANITY CHECK: stim trials agree between triggers and task struct
    if (~skipTaskChecks) && info.stimSession
        if length(find(rmmissing(task.var)))==length(thor_beh.stimSequence)
            temp = [];
            for i=1:length(thor_beh.sync)
                if i~=length(thor_beh.sync)
                    if sum((thor_beh.stimSequence > thor_beh.sync(i)) & (thor_beh.stimSequence < thor_beh.sync(i+1)))==1
                        temp = [temp,i];
                    elseif sum((thor_beh.stimSequence > thor_beh.sync(i)) & (thor_beh.stimSequence < thor_beh.sync(i+1)))>1
                        warning('Something weird is wrong.')
                    end
                else
                    if sum(thor_beh.stimSequence > thor_beh.sync(i))==1
                        temp = [temp,i];
                    elseif sum(thor_beh.stimSequence > thor_beh.sync(i))>1
                        warning('Something weird is wrong.')
                    end
                end
            end
            if ~all(find(task.var)==temp)
                warning('Mismatch in assignment of stim trials between task file and thor_beh stim triggers. But same number.')
            end
            if ~isempty(find((thor_beh.stimSequence-thor_beh.sync(find(task.var))<9) | (thor_beh.stimSequence-thor_beh.sync(find(task.var))>11)))
                warning('There is some wrong timing between SniffinSync and StimTrigger. Check manually.')
            end
        else
            if ~info.data.stimTriggerIssue_nanAfter
                warning('Mismatch in assignment of stim trials between task file and thor_beh stim triggers. Unequal number.')
            end
        end
    end
end

% number of trials and task events
if length(ts.beh.cue)~=info.task.numTrials
    warning(['Mismatch between number of trials: paq.cue=',num2str(length(ts.beh.cue)),', user input=',num2str(info.task.numTrials),'.'])
elseif length(ts.beh.cue)*2~=length(ts.beh.finalvalve_start)
    warning(['Mismatch between number of trials and valve openings: cue=',num2str(length(ts.beh.cue)),', finalvalve_start=',num2str(length(ts.beh.finalvalve_start)),'.'])
elseif length(ts.beh.cue)*2~=length(ts.beh.finalvalve_stop)
    warning(['Mismatch between number of trials and valve closures: cue=',num2str(length(ts.beh.cue)),', finalvalve_stop=',num2str(length(ts.beh.finalvalve_stop)),'.'])
elseif length(ts.beh.cue)~=length(ts.beh.responsewindow_start)
    warning(['Mismatch between number of trials and response window onsets: cue=',num2str(length(ts.beh.cue)),', responsewindow_start=',num2str(length(ts.beh.responsewindow_start)),'.'])
elseif length(ts.beh.cue)~=length(ts.beh.responsewindow_stop)
    warning(['Mismatch between number of trials and response window offsets: cue=',num2str(length(ts.beh.cue)),', responsewindow_stop=',num2str(length(ts.beh.responsewindow_stop)),'.'])
end

% rotary encoder data
if paq_beh.position(end)<100
    paq_beh.position(:) = NaN;
    paq_beh.speed(:) = NaN;
    paq_beh.acceleration(:) = NaN;
    warning(['Rotary encoder data corrupted as less than 1 m position change. Assigned NaN.'])
end


%% Lick patterns

ts.(info.epochs{e}).finalvalve_start1 = ts.(info.epochs{e}).finalvalve_start(1:2:end);
ts.(info.epochs{e}).finalvalve_start2 = ts.(info.epochs{e}).finalvalve_start(2:2:end);

% licks in s relative to first odour onset (od1)
paq_beh.licks_od1 = cell(1,info.task.numTrials);
for i=1:info.task.numTrials
    try
        temp = find(ts.beh.lick_start >= ts.beh.finalvalve_start1(i) & ts.beh.lick_start <= ts.beh.responsewindow_stop(i));
        if ~isempty(temp)
            paq_beh.licks_od1{i} = (ts.beh.lick_start(temp) - ts.beh.finalvalve_start1(i)) / info.paq.samplingRate;
        end
    catch
    end
end

% licks in s relative to second odour onset (od2)
paq_beh.licks_od2 = cell(1,info.task.numTrials);
for i=1:info.task.numTrials
    try
        temp = find(ts.beh.lick_start >= ts.beh.finalvalve_start2(i) & ts.beh.lick_start <= ts.beh.responsewindow_stop(i));
        if ~isempty(temp)
            paq_beh.licks_od2{i} = (ts.beh.lick_start(temp) - ts.beh.finalvalve_start2(i)) / info.paq.samplingRate;
        end
    catch
    end
end

% licks in s relative to response window onset (rw)
paq_beh.licks_rw = cell(1,info.task.numTrials);
for i=1:info.task.numTrials
    try
        temp = find(ts.beh.lick_start >= ts.beh.responsewindow_start(i) & ts.beh.lick_start <= ts.beh.responsewindow_stop(i));
        if ~isempty(temp)
            paq_beh.licks_rw{i} = (ts.beh.lick_start(temp) - ts.beh.responsewindow_start(i)) / info.paq.samplingRate;
        end
    catch
    end
end


%% Dealing with stim trigger issues or premature imaging stop

if isfield(info.data,'stimTriggerIssue_nanAfter')
    temp = fields(paq_beh);
    for i=1:length(temp)
        if size(paq_beh.(temp{i}),2)==info.task.numTrials
            paq_beh.([temp{i},'_raw']) = paq_beh.(temp{i});
            try
                paq_beh.(temp{i}) = [{paq_beh.([temp{i},'_raw']){:,1:info.data.stimTriggerIssue_nanAfter}},cell(1,info.task.numTrials-info.data.stimTriggerIssue_nanAfter)];
            catch
                paq_beh.(temp{i}) = nan(size(paq_beh.([temp{i},'_raw'])));
                paq_beh.(temp{i})(:,1:info.data.stimTriggerIssue_nanAfter) = paq_beh.([temp{i},'_raw'])(:,1:info.data.stimTriggerIssue_nanAfter);
            end
        end
    end
    disp('--- Dealt with stim trigger issue. Assigned NaN to paq_beh events after issue.')
end

if isfield(info.data,'prematureImagingStop_nanAfter')  %MB20211105
    temp = paq_beh.sync;
    paq_beh.sync = nan(1,info.task.numTrials);
    paq_beh.sync(1:info.data.prematureImagingStop_nanAfter) = temp(1:info.data.prematureImagingStop_nanAfter);
    temp = paq_beh.responsewindow;
    paq_beh.responsewindow = nan(1,info.task.numTrials);
    paq_beh.responsewindow(1:info.data.prematureImagingStop_nanAfter) = temp(1:info.data.prematureImagingStop_nanAfter);
    warning('Dealt with premature imaging stop. Assigned NaN to paq_beh events after issue.')
end


%% Save paq files

paq_pre = orderfields(paq_pre);
save([path.filepart_out,'paq_pre.mat'],'paq_pre','-v7.3');
disp(['--- Added paq_pre file to repo as ',[path.filepart_out,'paq_pre.mat'],'.'])
paq_beh = orderfields(paq_beh);
save([path.filepart_out,'paq_beh.mat'],'paq_beh','-v7.3');
disp(['--- Added paq_beh file to repo as ',[path.filepart_out,'paq_beh.mat'],'.'])
paq_post = orderfields(paq_post);
save([path.filepart_out,'paq_post.mat'],'paq_post','-v7.3');
disp(['--- Added paq_post file to repo as ',[path.filepart_out,'paq_post.mat'],'.'])


end


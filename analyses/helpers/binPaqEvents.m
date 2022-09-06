function [events_binned] = binPaqEvents(paq_beh,task,p,prop,cam_beh)

%% Bin continuous paq data

% position
temp = movmean(paq_beh.position',p.general.binSize,2,'omitnan');
events_binned.position = temp(floor(p.general.binSize/2)+1:p.general.binSize:end);
events_binned.position = events_binned.position(1:end-1);

% speed
velocity = paq_beh.speed;
if p.general.smoothingSd_preBinning_velocity~=0
	velocity = smoothdata(velocity,2,'gaussian',p.general.smoothingSd_preBinning_velocity*5);
end
temp = movmean(velocity',p.general.binSize,2,'omitnan');
events_binned.velocity = temp(floor(p.general.binSize/2)+1:p.general.binSize:end);
events_binned.velocity = events_binned.velocity(1:end-1);

% acceleration
acceleration = paq_beh.acceleration;
if p.general.smoothingSd_preBinning_acceleration~=0
	acceleration = smoothdata(acceleration,2,'gaussian',p.general.smoothingSd_preBinning_acceleration*5);
end
temp = movmean(acceleration',p.general.binSize,2,'omitnan');
events_binned.acceleration = temp(floor(p.general.binSize/2)+1:p.general.binSize:end);
events_binned.acceleration = events_binned.acceleration(1:end-1);


%% Bin discrete paq data

% sync
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.sync) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.sync = temp(:,prop.meanFrames2bins);
events_binned.sync = ceil(events_binned.sync(:,1:end-1)');
if length(find(events_binned.sync))~=prop.numTrials
    error('Something went wrong when binning the syncs in binPaqEvents.m')
end

% odourA
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.finalvalve(1,find(task.odour1=="A"))) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.odourA = temp(:,prop.meanFrames2bins);
events_binned.odourA = ceil(events_binned.odourA(:,1:end-1)');
if length(find(events_binned.odourA))~=prop.numTrials/2
    error('Something went wrong when binning the odourA in binPaqEvents.m')
end
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.finalvalve(2,find(task.odour1=="A"))) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.odourA_off = temp(:,prop.meanFrames2bins);
events_binned.odourA_off = ceil(events_binned.odourA_off(:,1:end-1)');

% odourX
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.finalvalve(1,find(task.odour1=="X"))) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.odourX = temp(:,prop.meanFrames2bins);
events_binned.odourX = ceil(events_binned.odourX(:,1:end-1)');
if length(find(events_binned.odourX))~=prop.numTrials/2
    error('Something went wrong when binning the odourX in binPaqEvents.m')
end
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.finalvalve(2,find(task.odour1=="X"))) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.odourX_off = temp(:,prop.meanFrames2bins);
events_binned.odourX_off = ceil(events_binned.odourX_off(:,1:end-1)');

% odourB
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.finalvalve(3,find(task.odour2=="B"))) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.odourB = temp(:,prop.meanFrames2bins);
events_binned.odourB = ceil(events_binned.odourB(:,1:end-1)');
if length(find(events_binned.odourB))~=prop.numTrials/2
    error('Something went wrong when binning the odourB in binPaqEvents.m')
end
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.finalvalve(4,find(task.odour2=="B"))) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.odourB_off = temp(:,prop.meanFrames2bins);
events_binned.odourB_off = ceil(events_binned.odourB_off(:,1:end-1)');

% odourY
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.finalvalve(3,find(task.odour2=="Y"))) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.odourY = temp(:,prop.meanFrames2bins);
events_binned.odourY = ceil(events_binned.odourY(:,1:end-1)');
if length(find(events_binned.odourY))~=prop.numTrials/2
    error('Something went wrong when binning the odourY in binPaqEvents.m')
end
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.finalvalve(4,find(task.odour2=="Y"))) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.odourY_off = temp(:,prop.meanFrames2bins);
events_binned.odourY_off = ceil(events_binned.odourY_off(:,1:end-1)');

% lick
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.licks(1,:)) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.lick = temp(:,prop.meanFrames2bins);
events_binned.lick = ceil(events_binned.lick(:,1:end-1)');
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.licks(2,:)) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.lick_off = temp(:,prop.meanFrames2bins);
events_binned.lick_off = ceil(events_binned.lick_off(:,1:end-1)');

% reward
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.reward(1,:)) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.reward = temp(:,prop.meanFrames2bins);
events_binned.reward = ceil(events_binned.reward(:,1:end-1)');
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.reward(2,:)) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.reward_off = temp(:,prop.meanFrames2bins);
events_binned.reward_off = ceil(events_binned.reward_off(:,1:end-1)');

% responsewindow
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.responsewindow(1,:)) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.responsewindow = temp(:,prop.meanFrames2bins);
events_binned.responsewindow = ceil(events_binned.responsewindow(:,1:end-1)');
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.responsewindow(2,:)) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.responsewindow_off = temp(:,prop.meanFrames2bins);
events_binned.responsewindow_off = ceil(events_binned.responsewindow_off(:,1:end-1)');


%% Get events for interaction terms

% odourPostA
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.finalvalve(3,find(task.odour1=="A"))) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.odourPostA = temp(:,prop.meanFrames2bins);
events_binned.odourPostA = ceil(events_binned.odourPostA(:,1:end-1)');

% odourPostX
temp = zeros(1,prop.numFramesTotal);
temp(paq_beh.finalvalve(3,find(task.odour1=="X"))) = 1;
temp = movmean(temp,p.general.binSize,2,'omitnan');
events_binned.odourPostX = temp(:,prop.meanFrames2bins);
events_binned.odourPostX = ceil(events_binned.odourPostX(:,1:end-1)');

% % odourBafterA
% temp = zeros(1,prop.numFramesTotal);
% temp(paq_beh.finalvalve(3,intersect(find(task.odour2=="B"),find(task.odour1=="A")))) = 1;
% temp = movmean(temp,p.general.binSize,2,'omitnan');
% events_binned.odourBafterA = temp(:,prop.meanFrames2bins);
% events_binned.odourBafterA = ceil(events_binned.odourBafterA(:,1:end-1)');
% 
% % odourBafterX
% temp = zeros(1,prop.numFramesTotal);
% temp(paq_beh.finalvalve(3,intersect(find(task.odour2=="B"),find(task.odour1=="X")))) = 1;
% temp = movmean(temp,p.general.binSize,2,'omitnan');
% events_binned.odourBafterX = temp(:,prop.meanFrames2bins);
% events_binned.odourBafterX = ceil(events_binned.odourBafterX(:,1:end-1)');
% 
% % odourYafterA
% temp = zeros(1,prop.numFramesTotal);
% temp(paq_beh.finalvalve(3,intersect(find(task.odour2=="Y"),find(task.odour1=="A")))) = 1;
% temp = movmean(temp,p.general.binSize,2,'omitnan');
% events_binned.odourYafterA = temp(:,prop.meanFrames2bins);
% events_binned.odourYafterA = ceil(events_binned.odourYafterA(:,1:end-1)');
% 
% % odourYafterX
% temp = zeros(1,prop.numFramesTotal);
% temp(paq_beh.finalvalve(3,intersect(find(task.odour2=="Y"),find(task.odour1=="X")))) = 1;
% temp = movmean(temp,p.general.binSize,2,'omitnan');
% events_binned.odourYafterX = temp(:,prop.meanFrames2bins);
% events_binned.odourYafterX = ceil(events_binned.odourYafterX(:,1:end-1)');


%% Get distance and duration from trial onset

% distance
events_binned.distance = nan(size(events_binned.position));
latest_sync = NaN;
position_at_trial_start = NaN; % NaN crashes GLM, so need to change to 0 on the way
for i=1:length(events_binned.distance)
    if events_binned.sync(i)==1
        latest_sync = i;
        position_at_trial_start = events_binned.position(i);
    end
    events_binned.distance(i) = events_binned.position(i)-position_at_trial_start;
end

% duration in bins
events_binned.duration = (1:length(events_binned.distance))';
latest_sync = NaN;
position_at_trial_start = NaN; % NaN crashes GLM, so need to change to 0 on the way
for i=1:length(events_binned.duration)
    if events_binned.sync(i)==1
        latest_sync = i;
        position_at_trial_start = i;
    end
    events_binned.duration(i) = events_binned.duration(i)-position_at_trial_start;
end


%% Bin cam data

if nargin==5
    
    % pupil
    temp = movmean(cam_beh.pupil',p.general.binSize,2,'omitnan');
    events_binned.pupil = temp(floor(p.general.binSize/2)+1:p.general.binSize:end);
    events_binned.pupil = events_binned.pupil(1:end-1);

    % nose
    temp = movmean(cam_beh.nose',p.general.binSize,2,'omitnan');
    events_binned.nose = temp(floor(p.general.binSize/2)+1:p.general.binSize:end);
    events_binned.nose = events_binned.nose(1:end-1);
    
    % motion
    temp = movmean(cam_beh.mot',p.general.binSize,2,'omitnan');
    events_binned.motion = temp(floor(p.general.binSize/2)+1:p.general.binSize:end);
    events_binned.motion = events_binned.motion(1:end-1);
end






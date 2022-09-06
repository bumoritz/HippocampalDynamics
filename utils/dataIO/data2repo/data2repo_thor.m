function [path,thor_beh] = data2repo_thor(info,p,path,raw)

%% Load and process Thorsync data

% Loading ThorSync data and detecting events
disp('--- Loading ThorSync data and detecting events.')
for e = 2 %1:length(info.epochs)
    if e~=2 | info.data.numFragments==1
        
        % ThorSync data -> ts
        temp = [path.folder_data,'Imaging\',info.animal,'-',info.date,'-',info.epochs{e},'\syncData\'];
        temp2 = dir(temp);
        path.(['file_in_thor_',info.epochs{e}]) = [temp,temp2(end-1).name];
        ts.(info.epochs{e}) = ThorLink_ReadThorSync(path.(['file_in_thor_',info.epochs{e}]));
        disp(['--- ThorSync ',info.epochs{e},' file loaded.'])
        
        % extract ThorSync data
        try
            ts2.(info.epochs{e}).frames = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.frames),'above',p.scope.voltageThreshold);
        catch
            ts2.(info.epochs{e}).frames = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.frames_alt),'above',p.scope.voltageThreshold);
        end
        try
            try
                ts2.(info.epochs{e}).stimSequence = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.stimSequence),'above',p.scope.voltageThreshold);
            catch
                ts2.(info.epochs{e}).stimSequence = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.stimSequence_alt),'above',p.scope.voltageThreshold);
            end
            ts2.(info.epochs{e}).stimPatternComplete = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.stimPatternComplete),'above',p.scope.voltageThreshold);
            ts2.(info.epochs{e}).stimIteration = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.stimIteration),'above',p.scope.voltageThreshold);
        catch
            if ~info.stimSession
                ts2.(info.epochs{e}).stimSequence = [];
                ts2.(info.epochs{e}).stimPatternComplete = [];
                ts2.(info.epochs{e}).stimIteration = [];
            else
                error('Some stim-related fields are missing in ThorSync files.')
            end
        end
        ts2.(info.epochs{e}).sync = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.sync),'above',p.scope.voltageThreshold);   

    else
        % ThorSync data -> ts
        for i=1:info.data.numFragments
            % ThorSync data -> ts
            temp = [path.folder_data,'Imaging\',info.animal,'-',info.date,'-',info.epochs{e},'-',num2str(i),'\syncData\'];
            temp2 = dir(temp);
            path.(['file_in_thor_',info.epochs{e},'_',num2str(i)]) = [temp,temp2(end-1).name];
            ts.(info.epochs{e}) = ThorLink_ReadThorSync(path.(['file_in_thor_',info.epochs{e},'_',num2str(i)]));
            disp(['--- ThorSync ',info.epochs{e},'-',num2str(i),' file loaded.'])
            
            % extract ThorSync data
            try
                ts2.([info.epochs{e},'_',num2str(i)]).frames = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.frames),'above',p.scope.voltageThreshold);
            catch
                ts2.([info.epochs{e},'_',num2str(i)]).frames = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.frames_alt),'above',p.scope.voltageThreshold);
            end           
            try
                try
                    ts2.([info.epochs{e},'_',num2str(i)]).stimSequence = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.stimSequence),'above',p.scope.voltageThreshold);
                catch
                    ts2.([info.epochs{e},'_',num2str(i)]).stimSequence = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.stimSequence_alt),'above',p.scope.voltageThreshold);
                end
                ts2.([info.epochs{e},'_',num2str(i)]).stimPatternComplete = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.stimPatternComplete),'above',p.scope.voltageThreshold);
                ts2.([info.epochs{e},'_',num2str(i)]).stimIteration = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.stimIteration),'above',p.scope.voltageThreshold);
            catch
                if ~info.stimSession
                    ts2.([info.epochs{e},'_',num2str(i)]).stimSequence = [];
                    ts2.([info.epochs{e},'_',num2str(i)]).stimPatternComplete = [];
                    ts2.([info.epochs{e},'_',num2str(i)]).stimIteration = [];
                else
                    error('Some stim-related fields are missing in ThorSync files.')
                end
            end
            ts2.([info.epochs{e},'_',num2str(i)]).sync = detectThresholdCrossing(ts.(info.epochs{e}).(info.scope.chan.sync),'above',p.scope.voltageThreshold);
        end
    end
end
%clear ts; 


%% Deal with fragmentations and additional frame triggers

% Dealing with fragmentations and additional frame triggers
disp('--- Dealing with fragmentations and additional frame triggers.')

thor_beh.numFrameTriggers_original = [];
thor_beh.numFrameTriggers_postDropping = [];
if info.data.numFragments~=1
    ts2.beh.frames = [];
    ts2.beh.stimSequence = [];
    ts2.beh.stimPatternComplete = [];
    ts2.beh.stimIteration = [];
    ts2.beh.sync = [];
    temp2 = 0;
    for i=1:info.data.numFragments
        thor_beh.numFrameTriggers_original(i) = length(ts2.([info.epochs{e},'_',num2str(i)]).frames);
        thor_beh.numFrameTriggers_postDropping(i) = raw.([info.epochs{e},'_',num2str(i)]).frames;
        ts2.beh.frames = [ts2.beh.frames, temp2+ ts2.([info.epochs{e},'_',num2str(i)]).frames(1:thor_beh.numFrameTriggers_postDropping(i))];
        
        % MB20210719: rectified this to make stim trigger assingment completely independent from error-prone thor_beh.sync
        try
            ts2.beh.sync = [ts2.beh.sync, temp2+ ts2.([info.epochs{e},'_',num2str(i)]).sync(1:info.data.fragments_numTrials(i))];
        catch
            ts2.beh.sync = [ts2.beh.sync, []];
            warning('Skipped thor_beh.sync for at least one fragment as it made trouble. Needs to be replaced using paq_beh.sync.')
        end
        if ~info.stimSession
            try
                ts2.beh.stimSequence = [ts2.beh.stimSequence, temp2+ ts2.([info.epochs{e},'_',num2str(i)]).stimSequence(1:info.data.fragments_numTrials(i)*(info.task.numStimTrials/info.task.numTrials))];
                ts2.beh.stimPatternComplete = [ts2.beh.stimPatternComplete, temp2+ ts2.([info.epochs{e},'_',num2str(i)]).stimPatternComplete(1:info.data.fragments_numTrials(i)*(info.task.numStimTrials/info.task.numTrials)*info.scope.patternsPerSequence)];
                ts2.beh.stimIteration = [ts2.beh.stimIteration, temp2+ ts2.([info.epochs{e},'_',num2str(i)]).stimIteration];
            catch
                ts2.beh.stimSequence = [];
                ts2.beh.stimPatternComplete = [];
                ts2.beh.stimIteration = [];
            end
        else
            ts2.beh.stimSequence = [ts2.beh.stimSequence, temp2+ ts2.([info.epochs{e},'_',num2str(i)]).stimSequence(1:info.data.fragments_numTrials(i)*(info.task.numStimTrials/info.task.numTrials))];
            ts2.beh.stimPatternComplete = [ts2.beh.stimPatternComplete, temp2+ ts2.([info.epochs{e},'_',num2str(i)]).stimPatternComplete(1:info.data.fragments_numTrials(i)*(info.task.numStimTrials/info.task.numTrials)*info.scope.patternsPerSequence)];
            ts2.beh.stimIteration = [ts2.beh.stimIteration, temp2+ ts2.([info.epochs{e},'_',num2str(i)]).stimIteration];
        end

        temp2 = temp2 + ts2.([info.epochs{e},'_',num2str(i)]).frames(end);
    end
else
    thor_beh.numFrameTriggers_original = length(ts2.beh.frames);
    thor_beh.numFrameTriggers_postDropping = raw.beh.frames;
    if thor_beh.numFrameTriggers_postDropping > thor_beh.numFrameTriggers_original
        warning('The thor_beh file has less frames than the movie. Assigned NaN to later frame times.')
        temp = ts2.beh.frames;
        ts2.beh.frames = nan(1,thor_beh.numFrameTriggers_postDropping);
        ts2.beh.frames(1:length(temp)) = temp;
    else
        ts2.beh.frames = ts2.beh.frames(1:thor_beh.numFrameTriggers_postDropping);
    end
end
if any(thor_beh.numFrameTriggers_original-thor_beh.numFrameTriggers_postDropping>2)
    if thor_beh.numFrameTriggers_postDropping==0
        warning('There are somehow no frame triggers in ThorSync beh file.')
    else
        warning('For at least one ThorSync beh fragment more than 2 frame triggers were dropped.')
    end
end


%% Identify imaging frames with scope events

disp('--- Identifying imaging frames with scope events.')

thor_beh.numFrames = sum(thor_beh.numFrameTriggers_postDropping);
temp = [-Inf, rmmissing(ts2.beh.frames(2:end)), +Inf]; % even if an event starts at the very end of a frame, it's still counted to the same frame
thor_beh.sync = discretize(ts2.beh.sync, temp);
thor_beh.stimSequence = discretize(ts2.beh.stimSequence, temp);

% SANITY CHECK: number of different stim events
if info.stimSession
    if info.task.numStimTrials ~= length(thor_beh.stimSequence);
        warning('Mismatch in number of stim trials between ThorSync file and user input')
    end
    if ~isfield(info.data,'stimTriggerIssue_nanAfter')
        if length(ts2.beh.stimPatternComplete)/info.scope.patternsPerSequence ~= length(thor_beh.stimSequence);
            warning('Mismatch in number of pattern completes and user input')
        end
    end
end

% identify frames with photoartefact
if info.stimSession
    thor_beh.stimPatternComplete = discretize(ts2.beh.stimPatternComplete, temp);
    temp2 = discretize(ts2.beh.stimIteration, temp);
    thor_beh.stimPatternOnset = zeros(size(thor_beh.stimPatternComplete));
    thor_beh.artefact = [];
    for i=1:length(thor_beh.stimPatternComplete)
        if i==1
            thor_beh.stimPatternOnset(i) = temp2(min( intersect(find(temp2<thor_beh.stimPatternComplete(1)),find(temp2>=thor_beh.stimSequence(1))) ));
        else
            thor_beh.stimPatternOnset(i) = temp2(min( intersect(find((temp2<thor_beh.stimPatternComplete(i))&(temp2>thor_beh.stimPatternComplete(i-1))),find(temp2>=thor_beh.stimSequence(fix((i-1)/info.scope.patternsPerSequence)+1))) ));
        end
        thor_beh.artefact = [thor_beh.artefact, thor_beh.stimPatternOnset(i):thor_beh.stimPatternComplete(i)];
    end
    
    % fixer
    if strcmp(info.animal,'Jobs') && strcmp(info.date,'20220721')
        thor_beh.stimPatternOnset = horzcat(thor_beh.stimPatternOnset(1:10200),379812,thor_beh.stimPatternOnset(10201:12799));
        thor_beh.stimPatternComplete = horzcat(thor_beh.stimPatternComplete(1:10200),379815,thor_beh.stimPatternComplete(10201:12799));
    end
    
    % dealing with stim trigger issues
    if ~isfield(info.data,'stimTriggerIssue_nanAfter')
        thor_beh.stimPatternOnset = reshape(thor_beh.stimPatternOnset,length(thor_beh.stimPatternOnset)/length(thor_beh.stimSequence),length(thor_beh.stimSequence));
        thor_beh.stimPatternComplete = reshape(thor_beh.stimPatternComplete,length(thor_beh.stimPatternComplete)/length(thor_beh.stimSequence),length(thor_beh.stimSequence));
    else
        temp2 = thor_beh.sync(1:info.data.stimTriggerIssue_nanAfter);
        thor_beh.sync = nan(1,length(thor_beh.sync));
        thor_beh.sync(1:length(temp2)) = temp2;
        temp2 = thor_beh.stimSequence;
        thor_beh.stimSequence = nan(1,length(thor_beh.stimSequence));
        temp3 = max(find(temp2<nanmax(thor_beh.sync)))+1;
        thor_beh.stimSequence(1:temp3) = temp2(1:temp3); 
        temp = reshape(thor_beh.stimPatternOnset(1:20*temp3),info.scope.patternsPerSequence,temp3);
        thor_beh.stimPatternOnset = nan(info.scope.patternsPerSequence,length(thor_beh.stimSequence));
        thor_beh.stimPatternOnset(:,1:temp3) = temp;
        temp = reshape(thor_beh.stimPatternComplete(1:20*temp3),info.scope.patternsPerSequence,temp3);
        thor_beh.stimPatternComplete = nan(info.scope.patternsPerSequence,length(thor_beh.stimSequence));
        thor_beh.stimPatternComplete(:,1:temp3) = temp;
        disp('--- Dealt with stim trigger issue. Assigned NaN to thor_beh events after issue.')
    end

    % SANITY CHECK: frames with stim onset
    if find(abs(thor_beh.stimPatternOnset(1,1:end)-thor_beh.stimSequence)>1)
        thor_beh.problematicStims = find(abs(thor_beh.stimPatternOnset(1,1:end)-thor_beh.stimSequence)>1);
        warning(['Mismatch in onset frame between stim trigger and first stim pattern above one frame in ',num2str(length(thor_beh.problematicStims)),' trials'])
    else
        thor_beh.problematicStims = [];
    end
end
%clear ts2;

% Special case fixer
if isfield(info.data,'specialCaseFixer')
    if info.data.specialCaseFixer == 1
        thor_beh.sync = thor_beh.sync(9:end);
    end
end

%% Save thor_beh file

% sanity check again (MB20211104)
if length(thor_beh.sync)~=info.task.numTrials
    warning('Incorrect number of SniffinSyncs in thor_beh file.')
end
if (info.stimSession && length(thor_beh.stimSequence)~=info.task.numTrials*0.8)
    warning('Incorrect number of stim triggers in thor_beh file.')
end

% remove events that happen after imaging has stopped
if length(thor_beh.sync)~=length(unique(thor_beh.sync))
    if isfield(info.data,'prematureImagingStop_nanAfter')
        thor_beh.sync_all = thor_beh.sync;
        thor_beh.stimSequence_all = thor_beh.stimSequence;
        warning('Dealt with premature imaging stop for thor_beh file.')
        if info.stimSession
            warning('Truncated imaging not yet implemented for stim sessions');
        end
    end
    temp2 = unique(thor_beh.sync);
    temp2 = temp2(1:end-1);
    thor_beh.sync = nan(1,info.task.numTrials);
    thor_beh.sync(1:length(temp2)) = temp2;
    temp2 = unique(thor_beh.stimSequence);
    temp2 = temp2(1:end-1);
    thor_beh.stimSequence = nan(1,info.task.numTrials*0.8);
    thor_beh.stimSequence(1:length(temp2)) = temp2;
    if isfield(thor_beh,'stimPatternOnset')
        temp = [];
        for i=1:size(thor_beh.stimPatternOnset,1)
            temp = [temp;unique(thor_beh.stimPatternOnset(i,:))];
        end
        thor_beh.stimPatternOnset = temp;
        temp = [];
        for i=1:size(thor_beh.stimPatternComplete,1)
            temp = [temp;unique(thor_beh.stimPatternComplete(i,:))];
        end
        thor_beh.stimPatternComplete = temp;
    end
    warning('Removed events that happened after imaging has stopped from thor_beh file.')
end


thor_beh = orderfields(thor_beh);
save([path.filepart_out,'thor_beh.mat'],'thor_beh','-v7.3');
disp(['--- Added thor_beh file to repo as ',[path.filepart_out,'thor_beh.mat'],'.'])


end

%% Online sequence analysis

clc; clear;

info.animal             = 'Python'; %'Mufasa';
info.date               = '20211210'; %'20200520';
info.recording          = 'Python-20211210-base1'; %'Mufasa-20200520-5min';
info.rootPath           = 'N:\Data\SniffinHippo\'; % 'N:\Data\SniffinHippo\';


%%% --------------------------- %%%
%%% ------ DON'T CHANGE ------- %%%
%%% --------------------------- %%%


info.numTrials                          = 60;

% sync data
p.syncData                              = 'paq_raw'; % 'paq_raw' or 'sniffinSyncs' or 'thor'
p.sync.voltageThreshold                 = 0.5; % [V], threshold for digitising
info.scope.frameRate                    = 30.04; % [Hz]
info.scope.chan.frames                  = 'FrameTrigger';
info.scope.chan.frames_alt              = 'FrameOut';
info.scope.chan.sync                    = 'SniffinSync';
info.paq.config.imaging                 = 7;
info.paq.config.cue                     = 8;

% sequence cell analysis
p.sca.zscore                            = true;
p.sca.smoothingSd_preBinning            = 6; % in frames (not in bins)
p.sca.binSize                           = 6; % number of frames per bin (should really be odd number)
p.sca.smoothingSd_postBinning           = 0; % in bins (not in frames) [was 0, 2 with no pre looks nice as well]   % in Sami's, no pre-smoothing, but post-binning smoothing of 10 bins (not *5)
p.sca.nannoniscell                      = true;
p.sca.frames_pre                        = 3*30; % (ideally dividable by binSize)
p.sca.frames_core                       = 5.3*30; % (ideally dividable by binSize)
p.sca.frames_post                       = 9+15+30+2*30; % (ideally dividable by binSize)
p.sca.frames_baselineWindow             = 3*30; % (ideally dividable by binSize)
p.sca.frames_odour2                     = 0.8*30; % (ideally dividable by binSize)   % seems to short; also consider e.g. in Changa (session with shorter delay between odour2 and response window)
p.sca.rgnSeed                           = 1234;
p.sca.numShuffles                       = 1000;
p.sca.shufflingThreshold                = 95;
p.sca.firingFieldBoundaries             = 0.5; % percent over baseline. 0.5 will lead to PF width being the FWHM
p.sca.activeInFiringFieldSd             = 1; % number of sd above baseline to be called active in firing field [default: 1] 3 getted rid von den prae Aktivitaeten
p.sca.reliabilityCriterion              = 0.2; %0.1; % fraction of trials with activity in firing field [default: 0.2] 0.3 sieht besser aus, gerade bei activeInFiringFieldSd=3
p.sca.earlyVsLateThreshold              = 1; % s after first odour onset [default: 1]
p.sca.qThreshold                        = 0.05;

% experiment structure
info.task.trialStructure.tOdour1        = 0.3; % [s]
info.task.trialStructure.tGap           = 5; % [s]
info.task.trialStructure.tOdour2        = 0.3; % [s]
info.task.trialStructure.tRespDelay     = 0.5; % [s]
info.task.trialStructure.tRespWindow    = 1; % [s]

% colours
p.col.odour                             = [252,193,2]/255; % yellow
p.col.reward                            = [128,204,223]/255; % light blue


%% Load data

disp(['- Loading data.'])

info.year = info.date(1:4);
info.month = info.date(5:6);
info.day = info.date(7:8);
path.homeFolder = [info.rootPath,info.year,'\',info.year,'-',info.month,'\',info.year,'-',info.month,'-',info.day,'\',info.animal,'\','Targeting','\'];
path.s2pPath = [path.homeFolder,info.recording,'\suite2p\plane0\Fall.mat'];
path.iscellPath = [path.homeFolder,info.recording,'\suite2p\plane0\iscell.npy']; % just as a control, is not used otherwise
temp = dir([path.homeFolder,info.animal,'_',info.date,'_*_SEQ.txt']);
path.seqPath = [path.homeFolder,temp.name];

% import seq file
seq = readmatrix(path.seqPath);
disp(['--- seq file loaded.'])
if size(seq,2)~=info.numTrials
    warning('seq file has incorrect number of trials.')
end

% import sync file
if strcmp(p.syncData,'sniffinSyncs')
    path.syncPath = [path.homeFolder,info.animal,'_',info.date,'_base1_sniffinSyncs.mat'];
    load(path.syncPath);
    disp(['--- paq_extr file loaded.'])
elseif strcmp(p.syncData,'paq_raw')
    path.syncPath = [path.homeFolder,info.animal,'_',info.date,'_base1.paq'];
    paq = paq2lab_ov20220105(path.syncPath); 
    ts.cue = detectThresholdCrossing_ov20220105(paq(:,info.paq.config.cue),'above',p.sync.voltageThreshold);
    ts.imaging = detectThresholdCrossing_ov20220105(paq(:,info.paq.config.imaging),'above',p.sync.voltageThreshold);
    temp = [-Inf, ts.imaging(2:end)', Inf]; % even if an event starts at the very end of a frame, it's still counted to the same frame
    if any(isnan(temp))
        temp2 = find(isnan(temp));
        temp(temp2) = linspace(temp(temp2(1)-1),temp(temp2(end)+1),length(temp2));
    end
    sniffinSyncs = discretize(ts.cue', temp);
    disp(['--- paq_extr file loaded.'])
elseif strcmp(p.syncData,'thor')
    path.syncPath = [path.homeFolder,info.recording,'\syncData\Episode001.h5'];
    ts = ThorLink_ReadThorSync_ov20220105(path.syncPath);
    try
        ts2.frames = detectThresholdCrossing_ov20220105(ts.(info.scope.chan.frames),'above',p.sync.voltageThreshold);
    catch
        ts2.frames = detectThresholdCrossing_ov20220105(ts.(info.scope.chan.frames_alt),'above',p.sync.voltageThreshold);
    end
    ts2.sync = detectThresholdCrossing_ov20220105(ts.(info.scope.chan.sync),'above',p.sync.voltageThreshold);
    ts2.numFrames = length(ts2.frames);
    temp = [-Inf, rmmissing(ts2.frames(2:end)), +Inf]; % even if an event starts at the very end of a frame, it's still counted to the same frame
    sniffinSyncs = discretize(ts2.sync, temp);
    disp(['--- thor file loaded.'])
end
if length(sniffinSyncs)~=info.numTrials
    warning('sync file has incorrect number of trials.')
end

% import suite2p and iscell files
s2p = load(path.s2pPath);
temp = readNPY(path.iscellPath);
iscell = find(temp(:,1));
disp(['--- suite2p files loaded.'])


%% Preparations

disp(['- Processing data.'])

% pre-process activity measure
act = s2p.spks;
if p.sca.zscore
     act = nanzscore(act,[],2);
end
if p.sca.smoothingSd_preBinning~=0
	act = smoothdata(act,2,'gaussian',p.sca.smoothingSd_preBinning*5);
end

% get basic properties
prop.numRois = size(act,1);
prop.numFramesTotal = size(act,2);
prop.numTrials = length(sniffinSyncs);
prop.numFrames_raw = p.sca.frames_pre + p.sca.frames_core + p.sca.frames_post;

% identify trial bounds
trial_bounds(1,:) = sniffinSyncs-p.sca.frames_pre;
trial_bounds(2,:) = sniffinSyncs+p.sca.frames_core+p.sca.frames_post-1;
for i=1:prop.numTrials
    prop.trial_frames(:,i) = trial_bounds(1,i) : trial_bounds(2,i);
end

% split into trials
nft = reshape(act(:,prop.trial_frames),[],prop.numFrames_raw,prop.numTrials);

% binning
temp = movmean(nft,p.sca.binSize,2,'omitnan');
nft_binned = temp(:,floor(p.sca.binSize/2)+1:p.sca.binSize:end,:);
nft_binned = nft_binned(:,1:end-1,:);
if p.sca.smoothingSd_postBinning~=0
	nft_binned = smoothdata(nft_binned,2,'gaussian',p.sca.smoothingSd_postBinning*5);
end

% get basic properties of binned data
prop.numFrames_binned = size(nft_binned,2);
prop.frames_pre_binned = p.sca.frames_pre/p.sca.binSize;
prop.frames_core_binned = p.sca.frames_core/p.sca.binSize;
prop.frames_post_binned = p.sca.frames_post/p.sca.binSize;
prop.frames_odour2_binned = p.sca.frames_odour2/p.sca.binSize;
prop.t_binned = ([1:prop.numFrames_binned]-prop.frames_pre_binned-1)*p.sca.binSize/info.scope.frameRate;
prop.frames_baselineWindow_binned = p.sca.frames_baselineWindow/p.sca.binSize;
prop.baselineWindow = prop.frames_pre_binned-prop.frames_baselineWindow_binned+1:prop.frames_pre_binned;
if floor(prop.frames_pre_binned)==prop.frames_pre_binned
    prop.analysisWindow = prop.frames_pre_binned+1:floor(prop.frames_pre_binned+prop.frames_core_binned);
    if floor(prop.frames_core_binned)==prop.frames_core_binned
        prop.odour2Window = prop.frames_pre_binned+prop.frames_core_binned+1:prop.frames_pre_binned+prop.frames_core_binned+prop.frames_odour2_binned;
    else
        prop.odour2Window = prop.frames_pre_binned+ceil(prop.frames_core_binned)+1:prop.frames_pre_binned+ceil(prop.frames_core_binned)+prop.frames_odour2_binned;
    end
else
    error('Check binning for analysis window.')
end


%% Sequence cell identifiation - all trials

disp(['- Sequence cell identification.'])

% % identify trial types
trials.trials_A = sort([find(seq(1,:)==1),find(seq(1,:)==3)]);
trials.trials_X = sort([find(seq(1,:)==2),find(seq(1,:)==4)]);
trials.trials_B = sort([find(seq(1,:)==1),find(seq(1,:)==4)]);
trials.trials_Y = sort([find(seq(1,:)==2),find(seq(1,:)==3)]);

sca_general = sequenceCellIdentification_ov20220105(nft_binned,trials,prop,iscell,p);
sca = sca_general;


%% Make snake plot

F = snakePlots_ov20220105(sca,p,info);


%% Save data

disp(['- Saving sca file.'])

savefig(F,[path.homeFolder,info.animal,'_',info.date,'_sca_8tile.fig']);
saveas(F,[path.homeFolder,info.animal,'_',info.date,'_sca_8tile.png']);
save([path.homeFolder,info.animal,'_',info.date,'_sca.mat'],'sca','-v7.3');

disp('- Done.')







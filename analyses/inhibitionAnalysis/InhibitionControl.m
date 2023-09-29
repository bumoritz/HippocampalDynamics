%% Inhibition control


%% Load data

% spks
s2p = load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\suite2p\plane0\Fall.mat');
load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\background.mat');
numFrames = size(s2p.F,2);
spks.unblocked_765 = s2p.spks(:,1:25333); % 1, 25333
spks.blocked_765 = s2p.spks(:,25333+1:25333+23676); % 2, 23676
spks.unblocked_930 = s2p.spks(:,25333+23676+1:25333+23676+22739); % 3, 22739
spks.blocked_930 = s2p.spks(:,25333+23676+22739+1:25333+23676+22739+22724); % 4, 22724
F.unblocked_765 = s2p.F(:,1:25333); % 1, 25333
F.blocked_765 = s2p.F(:,25333+1:25333+23676); % 2, 23676
F.unblocked_930 = s2p.F(:,25333+23676+1:25333+23676+22739); % 3, 22739
F.blocked_930 = s2p.F(:,25333+23676+22739+1:25333+23676+22739+22724); % 4, 22724
Fneu.unblocked_765 = s2p.Fneu(:,1:25333); % 1, 25333
Fneu.blocked_765 = s2p.Fneu(:,25333+1:25333+23676); % 2, 23676
Fneu.unblocked_930 = s2p.Fneu(:,25333+23676+1:25333+23676+22739); % 3, 22739
Fneu.blocked_930 = s2p.Fneu(:,25333+23676+22739+1:25333+23676+22739+22724); % 4, 22724
bkg.unblocked_765 = background(:,1:25333); % 1, 25333
bkg.blocked_765 = background(:,25333+1:25333+23676); % 2, 23676
bkg.unblocked_930 = background(:,25333+23676+1:25333+23676+22739); % 3, 22739
bkg.blocked_930 = background(:,25333+23676+22739+1:25333+23676+22739+22724); % 4, 22724
iscell = s2p.iscell(:,1);

% sniffinSyncs
temp = load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_765_blocked_sniffinSyncs.mat');
sniffinSyncs.blocked_765 = temp.sniffinSyncs;
temp = load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_765_unblocked_sniffinSyncs.mat');
sniffinSyncs.unblocked_765 = temp.sniffinSyncs;
temp = load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_930_blocked_sniffinSyncs.mat');
sniffinSyncs.blocked_930 = temp.sniffinSyncs;
temp = load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_930_unblocked_sniffinSyncs.mat');
sniffinSyncs.unblocked_930 = temp.sniffinSyncs;

% task
p.trialsPerBlock = 20; p.localWindow = 20;
trl = extractTrl('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_123211_Exp_NoSwitch_SEQ.txt',...
    'C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_214911_SNH.txt',...
    'C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_214917_PYB.mat',...
    'C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\out\',p);
task.odour1 = trl.odour1;
task.odour2 = trl.odour2;
task.response = trl.response;
task.autoreward = zeros(size(task.response));

% other preparations
p = get_p; 
info = get_info;
ops.tng.do_allTrials = true;
ops.tng.do_60t = false;
ops.tng.do_60t_onlyFirst = false;
ops.tng.do_100t = false;
ops.tng.do_100t_onlyFirst = false;
ops.tng.do_eventWiseAnalysis = false;
ops.tng.skip_boringSnakePlots = true;


%% blocked_765

act = spks.blocked_765; sync_beh.sync = sniffinSyncs.blocked_765; 

[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
[traces_all,avgTraces_all,temp] = createTracesStructs(nft_binned,trials_all);

these_avgTraces = avgTraces_all.A;

this_title = 'blocked, 765 nm, A';


%% unblocked_765

act = spks.unblocked_765; sync_beh.sync = sniffinSyncs.unblocked_765; 

[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
[traces_all,avgTraces_all,temp] = createTracesStructs(nft_binned,trials_all);

these_avgTraces = avgTraces_all.A;

this_title = 'unblocked, 765 nm, A';


%% blocked_930

act = spks.blocked_930; sync_beh.sync = sniffinSyncs.blocked_930; 

[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
[traces_all,avgTraces_all,temp] = createTracesStructs(nft_binned,trials_all);

these_avgTraces = avgTraces_all.A;

this_title = 'blocked, 930 nm, A';


%% unblocked_930

act = spks.unblocked_930; sync_beh.sync = sniffinSyncs.unblocked_930; 

[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
[traces_all,avgTraces_all,temp] = createTracesStructs(nft_binned,trials_all);

these_avgTraces = avgTraces_all.A;

this_title = 'unblocked, 930 nm, A';


%% popActSuppPanel

% assign windows
this_baselineWindow = 1:10;
this_sortingWindow = 16:20;
this_plottingWindow = 1:25;

% normalise traces
% first take:
% this_baseline_mean = nanmean(these_avgTraces(:,this_baselineWindow),2);
% this_baseline_std = nanstd(these_avgTraces(:,this_baselineWindow),[],2);
% these_normAvgTraces = (these_avgTraces - this_baseline_mean) ./ this_baseline_std;
% new take:
these_normAvgTraces = these_avgTraces;

% sort traces
these_idcs_unsorted = find(iscell==1);
[~,temp] = sort(nanmean(these_normAvgTraces(these_idcs_unsorted,this_sortingWindow),2));
these_idcs_sorted = these_idcs_unsorted(flip(temp));

figure;
plt.clim = [-1,1]; plt.colormap = redblue;
heatMap_task(these_normAvgTraces(these_idcs_unsorted,this_plottingWindow),NaN,[],p,info,plt);
title(this_title)



%% --- Testing preprocessing ---


%% Figure

act = F.unblocked_765; sync_beh.sync = sniffinSyncs.unblocked_765; 

p.inh.smoothingSd_preBinning = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('765 nm, unblocked, F, 40 trial avg')

%

act = F.blocked_765; sync_beh.sync = sniffinSyncs.blocked_765; 

p.inh.smoothingSd_preBinning = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('765 nm, blocked, F, 40 trial avg')

%

act = F.unblocked_930; sync_beh.sync = sniffinSyncs.unblocked_930; 

p.inh.smoothingSd_preBinning = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('930 nm, unblocked, F, 40 trial avg')

%

act = F.blocked_930; sync_beh.sync = sniffinSyncs.blocked_930; 

p.inh.smoothingSd_preBinning = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('930 nm, blocked, F, 40 trial avg')


%% Identifying background area

p.bkg.prctile = 8;
p.bkg.inset = 20;

this_meanImage = s2p.ops.meanImg;

this_meanImage_inset = this_meanImage(p.bkg.inset+1:size(this_meanImage,1)-p.bkg.inset,p.bkg.inset+1:size(this_meanImage,1)-p.bkg.inset);
this_meanImage_inset_lin = this_meanImage_inset(:);
these_background_idcs = find(this_meanImage_inset_lin<prctile(this_meanImage_inset_lin,p.bkg.prctile));
this_backgroundImage_inset_lin = zeros(size(this_meanImage_inset_lin));
this_backgroundImage_inset_lin(these_background_idcs)=1;
this_backgroundImage_inset = reshape(this_backgroundImage_inset_lin,size(this_meanImage_inset));
this_backgroundImage = zeros(size(this_meanImage));
this_backgroundImage(p.bkg.inset+1:size(this_meanImage,1)-p.bkg.inset,p.bkg.inset+1:size(this_meanImage,1)-p.bkg.inset) = this_backgroundImage_inset;
this_backgroundImage_lin = this_backgroundImage(:);

figure;
subplot(1,2,1)
imshow(this_meanImage,[nanmin(this_meanImage(:)),nanmax(this_meanImage(:))/4]);
title('average image')
subplot(1,2,2)
imshow(this_backgroundImage,[0,1]);
title('background mask')


%% Applying background mask on registered video

% path.imagingFile    = 'C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\registered_movie.raw';
% 
% background = nan(1,numFrames);
% for i=1:numFrames
%     
%     % load images one-by-one
%     this_fid = fopen(path.imagingFile,'r');
%     fseek(this_fid,(i-1)*512*512*2,'bof');
%     this_frame = uint16(fread(this_fid,512*512,'uint16',0,'l'));
%     %temp = reshape(temp,512,512);
%     frewind(this_fid);
%     fclose(this_fid);
%     
%     % calculate background
%     background(i) = nanmean(this_frame(find(this_backgroundImage_lin)));
%     
%     if mod(i,10000)==0
%         disp(num2str(i))
%     end
% end


%% Trying out background subtraction

act = F.unblocked_930; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.unblocked_930; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('930 nm unblocked, F')

%

act = F.blocked_930; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.blocked_930; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('930 nm blocked, F')

%

act = F.unblocked_765; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.unblocked_765; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('765 nm unblocked, F')

%

act = F.blocked_765; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.blocked_765; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('765 nm blocked, F')


%% Trying out background subtraction

act = Fneu.unblocked_930; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.unblocked_930; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('930 nm unblocked, Fneu')

%

act = Fneu.blocked_930; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.blocked_930; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('930 nm blocked, Fneu')

%

act = Fneu.unblocked_765; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.unblocked_765; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('765 nm unblocked, Fneu')

%

act = Fneu.blocked_765; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.blocked_765; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('765 nm blocked, Fneu')



%% Trying out background subtraction

act = F.unblocked_930 - Fneu.unblocked_930; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.unblocked_930; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('930 nm unblocked, F-Fneu')

%

act = F.blocked_930 - Fneu.blocked_930; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.blocked_930; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('930 nm blocked, F-Fneu')

%

act = F.unblocked_765 - Fneu.unblocked_765; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.unblocked_765; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('765 nm unblocked, F-Fneu')

%

act = F.blocked_765 - Fneu.blocked_765; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.blocked_765; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('765 nm blocked, F-Fneu')


%% Trying out background subtraction

act = F.unblocked_930 - 0.7*Fneu.unblocked_930; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.unblocked_930; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('930 nm unblocked, F-0.7*Fneu')

%

act = F.blocked_930 - 0.7*Fneu.blocked_930; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.blocked_930; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('930 nm blocked, F-0.7*Fneu')

%

act = F.unblocked_765 - 0.7*Fneu.unblocked_765; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.unblocked_765; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('765 nm unblocked, F-0.7*Fneu')

%

act = F.blocked_765 - 0.7*Fneu.blocked_765; %- bkg.blocked_930;  % repmat(bkg.blocked_930,2000,1);
sync_beh.sync = sniffinSyncs.blocked_765; 

p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('765 nm blocked, F-0.7*Fneu')


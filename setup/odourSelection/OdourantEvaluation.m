%% Source data and parameters

clear

basepath = 'D:\SniffinHippo\2019\2019-09\2019-09-29\';
basename = 'PID_Troubleshoot_pear_every10min_thenBaselineCorrect (hGain).paq';

chan_FinalValve = 9;
chan_PID = 6;

t_pre = 1; % s
t_stim = 1; % s
t_post = 5; % s
sd = 0.02; % s

frameRate = 10000; % Hz
downsampleFactor = 10;

ID = (1:4)';
Name = {'IsobutylPropionate_L','IsobutylPropionate_S','EthylValerate_L','EthylValerate_S'}';
Set = [1,1,1,1]';
File = {'','','',''}';
Pulses = {1:5,6:10,11:15,16:20}';
% ID = (1:32)';
% Name = {'Eucalyptol','IsoamylAcetate','Pinene','EthylButyrate',...
%     'Heptanone','Octanal','Terpinene','IsobutylPropionate',...
%     'BenzylAcetate','Benzaldehyde','Nonanone','Dimethoxybenzene',...
%     'MethylBenzoate','Allylanisole','Acetophenone','EthylValerate',...
%     'Limonene',...
%     'EthylButyrate_05_20_H_10','EthylButyrate_05_20_L_5','EthylButyrate_05_10_H_5','EthylButyrate_05_20_H_5',...
%     'EthylButyrate_1_20_H_10','EthylButyrate_1_20_L_5','EthylButyrate_1_10_H_5','EthylButyrate_1_20_H_5',...
%     'EthylValerate_10_H','EthylValerate_10_L','IsobutylPropionate_10_L',...
%     'Acetophenone_10_H','Acetophenone_10_L','BenzylAcetate_10_H','BenzylAcetate_10_L'}';
% Set = [1,1,1,1,...
%     2,2,2,2,...
%     3,3,3,3,...
%     4,4,4,4,...
%     5,...
%     6,6,6,6,...
%     7,7,7,7,...
%     8,8,8,...
%     9,9,9,9]';
% File = {'set1','set1','set1','set1',...
%     'set2','set2','set2_v34','set2_v4',...
%     'set3','set3','set3','set3_v4',...
%     'set4','set4','set4','set4',...
%     'set5',...
%     'set6_v1','set6_v2','set6_v3','set6_v4',...
%     'set7_v1','set7_v2_diedAfter4th','set7_v3_diedAfter4th','set7_v4_diedAfter4th',...
%     'set8_v12','set8_v12','set8_v4',...
%     'set9','set9','set9','set9'}';
% Pulses = {1:5,6:10,11:15,16:20,...
%     1:5,6:10,1:5,1:5,...
%     1:5,6:10,11:15,1:5,...
%     1:5,6:10,11:15,16:20,...
%     1:5,...
%     1:5,1:5,1:5,1:5,...
%     1:5,1:4,1:4,1:4,...
%     1:4,5:8,1:4,...
%     1:4,5:8,9:12,13:16}';
dat = table(ID,Name,Set,File,Pulses);
clear ID
clear Name
clear Set
clear File
clear Pulses

disp([num2str(frameRate/downsampleFactor),' frames correspond to 1 s.'])

%% Load data and save traces

numOdourants = size(dat,1);
numPulses = 4;%size(dat.Pulses{1},2);

pid = [];
PID = {};
for f=1:numOdourants
    paq = paq2lab([basepath,dat.File{f},basename]);
    paqd = paq(1:downsampleFactor:end,:);
    fr_onset = thresholdDetect(paqd(:,chan_FinalValve),'above',1);
    fr_onset = fr_onset(dat.Pulses{f});
    fr_start = fr_onset-t_pre*frameRate/downsampleFactor;
    fr_end = fr_onset+(t_stim+t_post)*frameRate/downsampleFactor;
    
    for p=1:numPulses
        pid(:,p) = smoothdata(paqd(fr_start(p):fr_end(p),chan_PID)','Gaussian',5*sd*frameRate/downsampleFactor);
    end
    PID{f,1} = pid;
end
dat = addvars(dat,PID);
numFrames = size(dat.PID{1},1);
clear PID   


%% Calculate metrics

% mean and sd traces
PID_mean = {};
PID_sd = {};
for f=1:numOdourants
    PID_mean{f,1} = nanmean(dat.PID{f},2);
    PID_sd{f,1} = nanstd(dat.PID{f},[],2);
end
dat = addvars(dat,PID_mean);
dat = addvars(dat,PID_sd);
clear PID_mean
clear PID_sd

% norm. traces
PID_norm = {};
PID_mean_norm = {};
PID_sd_norm = {};
for f=1:numOdourants
    PID_norm{f,1} = dat.PID{f};
    PID_mean_norm{f,1} = nanmean(dat.PID{f},2);
    PID_sd_norm{f,1} = nanstd(dat.PID{f},[],2);
end
dat = addvars(dat,PID_norm);
dat = addvars(dat,PID_mean_norm);
dat = addvars(dat,PID_sd_norm);
clear PID_norm
clear PID_mean_norm
clear PID_sd_norm


%% Baseline-subtracted metrics

% peak amplitude (V)
PeakAmplitude = {};
for f=1:numOdourants
    PeakAmplitude{f,1} = nanmax(dat.PID_mean{f}) - nanmean(dat.PID_mean{f}(1:t_pre*frameRate/downsampleFactor));
end
dat = addvars(dat,PeakAmplitude);
clear PeakAmplitude

% rise time (95% of peak; ms)
RiseTime = {};
for f=1:numOdourants
    RiseTime{f,1} = nanmin(find(dat.PID_mean{f}-nanmean(dat.PID_mean{f}(1:t_pre*frameRate/downsampleFactor))>0.95*dat.PeakAmplitude{f})) - t_pre*frameRate/downsampleFactor;
end
dat = addvars(dat,RiseTime);
clear RiseTime

% decay time (1/e of peak; ms)
criterion = 1/exp(1);
DecayTimeE = {};
for f=1:numOdourants
    temp = find(dat.PID_mean{f}-nanmean(dat.PID_mean{f}(1:t_pre*frameRate/downsampleFactor))<criterion*dat.PeakAmplitude{f});
    DecayTimeE{f,1} = nanmin(temp(find(temp>(t_pre+0.5*t_stim)*frameRate/downsampleFactor))) - (t_pre+t_stim)*frameRate/downsampleFactor;
end
dat = addvars(dat,DecayTimeE);
clear DecayTimeE
 
% decay time (5% of peak; ms)
criterion = 0.05;
DecayTime5 = {};
for f=1:numOdourants
    temp = find(dat.PID_mean{f}-nanmean(dat.PID_mean{f}(1:t_pre*frameRate/downsampleFactor))<criterion*dat.PeakAmplitude{f});
    DecayTime5{f,1} = nanmin(temp(find(temp>(t_pre+0.5*t_stim)*frameRate/downsampleFactor))) - (t_pre+t_stim)*frameRate/downsampleFactor;
end
dat = addvars(dat,DecayTime5);
clear DecayTime5
 
% decay time (1% of peak; ms)
criterion = 0.01;
DecayTime1 = {};
for f=1:numOdourants
    temp = find(dat.PID_mean{f}-nanmean(dat.PID_mean{f}(1:t_pre*frameRate/downsampleFactor))<criterion*dat.PeakAmplitude{f});
    DecayTime1{f,1} = nanmin(temp(find(temp>(t_pre+0.5*t_stim)*frameRate/downsampleFactor))) - (t_pre+t_stim)*frameRate/downsampleFactor;
end
dat = addvars(dat,DecayTime1);
clear DecayTime1


%% Plots

t = -t_pre:downsampleFactor/frameRate:(t_stim+t_post);

for f=1:numOdourants
    F = figure;
    for p=1:numPulses
        plot(t,dat.PID{f}(:,p),'color',[1,0,p/numPulses])
        hold on
    end
    plot(t,dat.PID_mean{f},'k-','LineWidth',3) %shadedErrorBar(t,dat.PID_mean{f},dat.PID_sd{f}) 
    hold off     
    lower = nanmin(nanmin(dat.PID{f})) - 0.05*nanmax(nanmax(dat.PID{f})); % 0;
    upper = nanmax(nanmax(dat.PID{f})) + 0.05*nanmax(nanmax(dat.PID{f}));
    ylim([lower,upper])
    xlabel('time (s)')
    ylabel('amplitude (V)')
    title(dat.Name{f})
    savefig(F,[basepath,'\',dat.Name{f},'.fig']);
end



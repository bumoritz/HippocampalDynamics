%% Select parameters

clc
clear

%title_string = 'x5_45cmOld';
paqPath     = 'C:\Data\SniffinHippo\2020\2020-06\2020-06-26\PID\PID_20200625_PA_o_2_4_7.paq';%strcat('C:\Data\PID_20200122\PID_',title_string,'.paq');
seqPath     = 'C:\Data\SniffinHippo\2020\2020-06\2020-06-26\PID\SEQ_121314_21it.txt';
savePath    = 'C:\Data\SniffinHippo\2020\2020-06\2020-06-26\PID\PA';
mkdir(savePath)

paqChannels = {'PID','FinalValve'};%{'Trigger','TrialOn','ResponseWindow','Odour','MFC1','MFC2','MFC3','PID'};
frameRate = 1000;

downsampleFactor = 1;


%% Load data

seq = importdata(seqPath)';
numTrials = size(seq,1);
stim = seq(:,1);
var = seq(:,2);
clear seq;

paq = paq2lab(paqPath);
paqd = paq(1:downsampleFactor:end,:);
numFrames = size(paqd,1);

% let's do it the stupid way
% Trigger         = paqd(:,find(contains(paqChannels,'Trigger')));
% TrialOn         = paqd(:,find(contains(paqChannels,'TrialOn')));
% ResponseWindow  = paqd(:,find(contains(paqChannels,'ResponseWindow')));
% Odour           = paqd(:,find(contains(paqChannels,'Odour')));
% MFC1            = paqd(:,find(contains(paqChannels,'MFC1')));
% MFC2            = paqd(:,find(contains(paqChannels,'MFC2')));
% MFC3            = paqd(:,find(contains(paqChannels,'MFC3')));
% PID             = paqd(:,find(contains(paqChannels,'PID')));
Odour           = paqd(:,9);
PID             = paqd(:,4);%paqd(:,4);%paqd(:,6);


%% if only one odour is presented each

pulses = thresholdDetect_MB(Odour,'above',0.5);
stim = stim(1:length(pulses));
vial1 = pulses(find(stim(1:length(pulses))==1));
vial2 = pulses(find(stim(1:length(pulses))==2));
vial3 = pulses(find(stim(1:length(pulses))==3));
vial4 = pulses(find(stim(1:length(pulses))==4));

%% Get frames with odour presentation onset
% 
% pulses = thresholdDetect(Odour,'above',0.5);
% 
% % fixer
% %stim = stim(1:(197/2)-1);
% 
% odourA_type1 = pulses(find(stim==1)*2-1);
% odourA_type3 = pulses(find(stim==3)*2-1);
% odourA_all = sort([odourA_type1;odourA_type3]);
% odourX_type2 = pulses(find(stim==2)*2-1);
% odourX_type4 = pulses(find(stim==4)*2-1);
% odourX_all = sort([odourX_type2;odourX_type4]);
% odourB_type1 = pulses(find(stim==1)*2);
% odourB_type4 = pulses(find(stim==4)*2);
% odourB_all = sort([odourB_type1;odourB_type4]);
% odourY_type2 = pulses(find(stim==2)*2);
% odourY_type3 = pulses(find(stim==3)*2);
% odourY_all = sort([odourY_type2;odourY_type3]);


%% Select data

plotTitle = 'R3-7-0d1%' % set1: 'A-IsobutylPropionate', 'B-Limonene', 'X-Prenol', 'Y-Benzaldehyde'
onsets = vial4(2:end);
onset_subset2 = []; % leave empty otherwise -> green
onset_subset1 = []; % leave empty otherwise -> red

stimSeconds = 0.3;

% %% raw and bls plots

%lower = -1;
%upper = 11;
preSeconds = 1;
%stimSeconds = 1;
postSeconds = 5;
filterOrder = 10;

F = figure;
for i=1:length(onsets)
    data = PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)+postSeconds*frameRate/downsampleFactor);%- median(PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)-1));
    data = medfilt1(data,filterOrder);
    
    if (~isempty(onset_subset1) && ~isempty(onset_subset2))
        if ismember(onsets(i),onset_subset1)
            plot((1:length(data))*downsampleFactor/frameRate,data,'color',[1,0,i/length(onsets)])
        elseif ismember(onsets(i),onset_subset2)
            plot((1:length(data))*downsampleFactor/frameRate,data,'color',[0,1,i/length(onsets)])
        end
    else
        plot((1:length(data))*downsampleFactor/frameRate,data,'color',[1,0,i/length(onsets)])
    end
    hold on
end
x = (1:length(data))*downsampleFactor/frameRate;
hold off
xlim([x(filterOrder),x(length(x)-filterOrder)])
%ylim([lower,upper])
xlabel('time [s]')
ylabel('PID voltage')
title([plotTitle])
savefig(F,[savePath,'\',plotTitle,'_raw.fig']);
close

% %% Generate full plot

preSeconds = 1;

postSeconds = stimSeconds+5;
postSeconds_plot = stimSeconds+5;
insetSeconds = 0.1;
filterOrder = 100;

%%% --- Make traces ---

data_raw = zeros(length(onsets(1)-preSeconds*frameRate/downsampleFactor:onsets(1)+postSeconds*frameRate/downsampleFactor),length(onsets));
for i=1:length(onsets)
    data_raw(:,i) = PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)+postSeconds*frameRate/downsampleFactor);
    %data_raw(:,i) = smoothdata(data_raw(:,i),'gaussian',10);
    %data_raw(:,i) = medfilt1(data_raw(:,i),filterOrder);
end

data_bls = zeros(length(onsets(1)-preSeconds*frameRate/downsampleFactor:onsets(1)+postSeconds*frameRate/downsampleFactor),length(onsets));
for i=1:length(onsets)
    
    GAUSS = 50;
    
    data_bls(:,i) = PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)+postSeconds*frameRate/downsampleFactor) - median( smoothdata(PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)-1),'gaussian',GAUSS));
    %data_bls(:,i) = smoothdata(data_bls(:,i),'gaussian',GAUSS);
    %data_bls(:,i) = medfilt1(data_bls(:,i),filterOrder);
end

data_raw_avg = nanmean(data_raw,2);
data_bls_avg = nanmean(data_bls,2);

% DEBUG - for odourA
% [~,type1_idcs]=intersect(sort([find(stim==1);find(stim==3)]),find(stim==1));
% [~,type3_idcs]=intersect(sort([find(stim==1);find(stim==3)]),find(stim==3));
% data_raw_avg_type1 = nanmean(data_raw(:,type1_idcs),2);
% data_raw_avg_type3 = nanmean(data_raw(:,type3_idcs),2);



%%% --- Extract properties and generate full plot ---

F = default_figure([25,3.5,6,4]);
x = (1:length(onsets(1)-preSeconds*frameRate/downsampleFactor:onsets(1)+postSeconds_plot*frameRate/downsampleFactor)-1)*downsampleFactor/frameRate;
yline(0,'k:');
hold on


% %% DEBUG - for odourA
% figure;
% plot(data_raw_avg_type1);
% hold on
% plot(data_raw_avg_type3);
% %shadedErrorBar(x,nanmean(data_raw(:,type1_idcs),2),sem(data_raw(:,type1_idcs),2),'r')

% wdw = 900:1040;
% figure;
% for i=1:length(type1_idcs);  
%     hold on
%     plot(data_raw(wdw,type1_idcs(i)),'g')
%     plot(data_raw(wdw,type3_idcs(i)),'r')
% end


% peak amplitude (V)
% out.PeakAmplitude = nanmean(data_bls_avg); %nanmax(data_bls_avg);
% yline(out.PeakAmplitude,'k:');

% core plot
for i=1:length(onsets)
    if (~isempty(onset_subset1) && ~isempty(onset_subset2))
        if ismember(onsets(i),onset_subset1)
            plot(x,data_bls(1:length(x),i),'color',[1,0,i/length(onsets)])
        elseif ismember(onsets(i),onset_subset2)
            plot(x,data_bls(1:length(x),i),'color',[0,1,i/length(onsets)])
        end
    else
        plot(x,data_bls(1:length(x),i),'color',[1,0,i/length(onsets)])
    end
end
plot(x,data_bls_avg(1:length(x)),'k','LineWidth',2)

% rise time (90% of peak; ms)
% try
%     out.RiseTime = nanmin(find(data_bls_avg>0.90*nanmax(data_bls_avg))) - preSeconds*frameRate/downsampleFactor;
%     plot(x(nanmin(find(data_bls_avg>0.90*nanmax(data_bls_avg)))),0.90*nanmax(data_bls_avg),'bx','LineWidth',3,'MarkerSize',10)
% catch
% end

% decay time (1/e of peak; ms)
try
    criterion = 1/exp(1);
    temp = find(data_bls_avg<criterion*nanmax(data_bls_avg));
    out.DecayTime_1e = nanmin(temp(find(temp>(preSeconds+0.5*stimSeconds)*frameRate/downsampleFactor))) - (preSeconds+stimSeconds)*frameRate/downsampleFactor;
    plot(x(nanmin(temp(find(temp>(preSeconds+0.5*stimSeconds)*frameRate/downsampleFactor)))),nanmax(data_bls_avg)*criterion,'bx','LineWidth',3,'MarkerSize',10)
catch
end

% decay time (5% of peak; ms)
try
    criterion = 0.05;
    temp = find(data_bls_avg<criterion*nanmax(data_bls_avg));
    out.DecayTime_5 = nanmin(temp(find(temp>(preSeconds+0.5*stimSeconds)*frameRate/downsampleFactor))) - (preSeconds+stimSeconds)*frameRate/downsampleFactor;
    plot(x(nanmin(temp(find(temp>(preSeconds+0.5*stimSeconds)*frameRate/downsampleFactor)))),nanmax(data_bls_avg)*criterion,'bx','LineWidth',3,'MarkerSize',10)
catch
end

% decay time (1% of peak; ms)
try
    criterion = 0.01;
    temp = find(data_bls_avg<criterion*nanmax(data_bls_avg));
    out.DecayTime_1 = nanmin(temp(find(temp>(preSeconds+0.5*stimSeconds)*frameRate/downsampleFactor))) - (preSeconds+stimSeconds)*frameRate/downsampleFactor;
    plot(x(nanmin(temp(find(temp>(preSeconds+0.5*stimSeconds)*frameRate/downsampleFactor)))),nanmax(data_bls_avg)*criterion,'bx','LineWidth',3,'MarkerSize',10)
catch
end

% depletion (minimum/maximum)
trialResponses = zeros(length(onsets),1);
for i=1:length(onsets)
    trialResponses(i) = nanmean(data_bls((preSeconds+insetSeconds)*frameRate/downsampleFactor:(preSeconds+stimSeconds-insetSeconds)*frameRate/downsampleFactor,i));
end
out.Depletion_fraction = min(trialResponses)/max(trialResponses);

% cosmetics
hold off
xlim([x(filterOrder),x(length(x)-filterOrder)])
xlabel('time (s)')
ylabel('PID voltage')
title([plotTitle])
savefig(F,[savePath,'\',plotTitle,'_full.fig']);

out


%% 

F = figure;
hold on
yline(0,'k:');
yline(0.01*max(nanmean(data_bls(1:length(x),:),2)),'r:');
shadedErrorBar(x,nanmean(data_bls(1:length(x),:),2),nanstd(data_bls(1:length(x),:),[],2)/sqrt(50))

xlim([0,7])
ylim([-0.2,2])
xlabel('time (s)')
ylabel('PID voltage')
xticks([1,2,3,4,5,6])
yticks([0,0.5,1,1.5,2])

savefig(F,[savePath,'\',plotTitle,'_full_sem.fig']);

%%

F = figure;
hold on
yline(0,'k:');
yline(0.01*max(nanmean(data_bls(1:length(x),:),2)),'r:');
shadedErrorBar(x,nanmean(data_bls(1:length(x),:),2),nanstd(data_bls(1:length(x),:),[],2)/sqrt(50))

xlim([0,7])
ylim([-0.01*max(nanmean(data_bls(1:length(x),:),2)),0.06*max(nanmean(data_bls(1:length(x),:),2))])
xlabel('time (s)')
ylabel('PID voltage')
xticks([1,2,3,4,5,6])
yticks([0,0.01*max(nanmean(data_bls(1:length(x),:),2)),0.05*max(nanmean(data_bls(1:length(x),:),2))])

savefig(F,[savePath,'\',plotTitle,'_full_sem_1perc.fig']);


%% Trial-type plot

plotTitle = 'all trials' % set1: 'A-IsobutylPropionate', 'B-Limonene', 'X-Prenol', 'Y-Benzaldehyde'
onsets = [odourA_all;odourX_all];
onset_subset1 = odourA_type1; % leave empty otherwise -> green
onset_subset2 = odourA_type3; % leave empty otherwise -> cyan
onset_subset3 = odourX_type2; % leave empty otherwise -> red
onset_subset4 = odourX_type4; % leave empty otherwise -> pink


% %% Generate full plot

preSeconds = 1;
stimSeconds = 1;
postSeconds = stimSeconds+17;
postSeconds_plot = stimSeconds+17;
insetSeconds = 0.1;
filterOrder = 100;

%%% --- Make traces ---

data_raw = zeros(length(onsets(1)-preSeconds*frameRate/downsampleFactor:onsets(1)+postSeconds*frameRate/downsampleFactor),length(onsets));
for i=1:length(onsets)
    data_raw(:,i) = PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)+postSeconds*frameRate/downsampleFactor);
    data_raw(:,i) = medfilt1(data_raw(:,i),filterOrder);
end

data_bls = zeros(length(onsets(1)-preSeconds*frameRate/downsampleFactor:onsets(1)+postSeconds*frameRate/downsampleFactor),length(onsets));
for i=1:length(onsets)
    data_bls(:,i) = PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)+postSeconds*frameRate/downsampleFactor);% - median(PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)-1));
    data_bls(:,i) = medfilt1(data_bls(:,i),filterOrder);
end

data_raw_avg = nanmean(data_raw,2);
data_bls_avg = nanmean(data_bls,2);

%%% --- Extract properties and generate full plot ---

F = default_figure([20,1.5,19,9]);
x = (1:length(onsets(1)-preSeconds*frameRate/downsampleFactor:onsets(1)+postSeconds_plot*frameRate/downsampleFactor)-1)*downsampleFactor/frameRate;
yline(0,'k:');
hold on

% core plot
for i=1:length(onsets)
    if ismember(onsets(i),onset_subset1)
        plot(x,data_bls(1:length(x),i),'color',[0,1,0])
    elseif ismember(onsets(i),onset_subset2)
        plot(x,data_bls(1:length(x),i),'color',[0,1,1])
    elseif ismember(onsets(i),onset_subset3)
        plot(x,data_bls(1:length(x),i),'color',[1,0,0])
    elseif ismember(onsets(i),onset_subset4)
        plot(x,data_bls(1:length(x),i),'color',[1,0,1])
    end
end

% cosmetics
hold off
xlim([x(filterOrder),x(length(x)-filterOrder)])
xlabel('time (s)')
ylabel('PID voltage')
title([plotTitle])
% savefig(F,[savePath,'\',plotTitle,'_full.fig']);

%% Trial-type plot

plotTitle = 'all trials' % set1: 'A-IsobutylPropionate', 'B-Limonene', 'X-Prenol', 'Y-Benzaldehyde'
onsets = [odourA_all;odourX_all];
onset_subset1 = odourA_type1; % leave empty otherwise -> green
onset_subset2 = odourA_type3; % leave empty otherwise -> cyan
onset_subset3 = odourX_type2; % leave empty otherwise -> red
onset_subset4 = odourX_type4; % leave empty otherwise -> pink


% %% Generate full plot

preSeconds = 1;
stimSeconds = 1;
postSeconds = stimSeconds+17;
postSeconds_plot = stimSeconds+17;
insetSeconds = 0.1;
filterOrder = 100;

%%% --- Make traces ---

data_raw = zeros(length(onsets(1)-preSeconds*frameRate/downsampleFactor:onsets(1)+postSeconds*frameRate/downsampleFactor),length(onsets));
for i=1:length(onsets)
    data_raw(:,i) = PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)+postSeconds*frameRate/downsampleFactor);
    data_raw(:,i) = medfilt1(data_raw(:,i),filterOrder);
end

data_bls = zeros(length(onsets(1)-preSeconds*frameRate/downsampleFactor:onsets(1)+postSeconds*frameRate/downsampleFactor),length(onsets));
for i=1:length(onsets)
    data_bls(:,i) = PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)+postSeconds*frameRate/downsampleFactor);% - median(PID(onsets(i)-preSeconds*frameRate/downsampleFactor:onsets(i)-1));
    data_bls(:,i) = medfilt1(data_bls(:,i),filterOrder);
end

data_raw_avg = nanmean(data_raw,2);
data_bls_avg = nanmean(data_bls,2);

%%% --- Extract properties and generate full plot ---

F = default_figure([20,1.5,19,9]);
x = (1:length(onsets(1)-preSeconds*frameRate/downsampleFactor:onsets(1)+postSeconds_plot*frameRate/downsampleFactor)-1)*downsampleFactor/frameRate;
yline(0,'k:');
hold on

% core plot
for i=1:length(onsets)
    if ismember(onsets(i),onset_subset1)
        plot(x,data_bls(1:length(x),i),'color',[0,1,0])
    elseif ismember(onsets(i),onset_subset2)
        plot(x,data_bls(1:length(x),i),'color',[0,1,1])
    elseif ismember(onsets(i),onset_subset3)
        plot(x,data_bls(1:length(x),i),'color',[1,0,0])
    elseif ismember(onsets(i),onset_subset4)
        plot(x,data_bls(1:length(x),i),'color',[1,0,1])
    end
end

% cosmetics
hold off
xlim([x(filterOrder),x(length(x)-filterOrder)])
xlabel('time (s)')
ylabel('PID voltage')
title([plotTitle])
% savefig(F,[savePath,'\',plotTitle,'_full.fig']);




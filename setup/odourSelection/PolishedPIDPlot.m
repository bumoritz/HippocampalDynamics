function [norm_H,odour_diff_norm_mean_H,odour_diff_norm_sem_H] = PolishedPIDPlot(name,norm_H,odour_diff_norm_mean_H,odour_diff_norm_sem_H)

%% Polished PID Plot

% final for
% set iii: BF_low, PF_low

%clear;
%clc;

p.name              = name; %'MethylButyrate';

p.numRecordings     = 5;
p.numTrials         = 20;
p.samplingRate      = 10000;

p.downsampleFactor  = 1;
p.smoothingType     = 'Gaussian'; % 'none', 'median', 'Gaussian'
p.smoothingAmount   = 100; %300; % [frames], sd for Gaussian
p.doBaselineSubtraction = true;
p.baselineSubtractionWindow = 1; % [s]

p.preSeconds        = 1;
p.stimSeconds       = 0.3;
p.postSeconds       = 5;

pid = {};
onsets_odour = [];
onsets_oil = [];
for i=1:p.numRecordings
    %temp = load(['Z:\WIBR_Hippos\SniffinHippo\Analysis\Odourants\EthylPropionate_low\EP_low_',num2str(i),'.mat']);

    % final recordings
    temp = load(['Z:\WIBR_Hippos\SniffinHippo\Analysis\Odourants\',p.name,'\',p.name,'_',num2str(i),'.mat']);
    
    pid{1,i} = temp.pid.data;
    onsets_odour(:,i) = temp.pid.odour;
    onsets_oil(:,i) = temp.pid.oil;
end

numSamples = length(onsets_odour(1,1)-p.preSeconds*p.samplingRate/p.downsampleFactor:onsets_odour(1,1)+p.postSeconds*p.samplingRate/p.downsampleFactor);


%% Filtering

data = {};
if strcmp(p.smoothingType,'none')
    data = pid;
elseif strcmp(p.smoothingType,'median')
    for i=1:p.numRecordings
        data{1,i} = medfilt1(pid{1,i},p.smoothingAmount);
    end
elseif strcmp(p.smoothingType,'Gaussian')
    for i=1:p.numRecordings
        data{1,i} = smoothdata(pid{1,i},'Gaussian',p.smoothingAmount*5);
    end
end


%% Extracting traces

t = (1:numSamples)*p.downsampleFactor/p.samplingRate - p.preSeconds-p.downsampleFactor/p.samplingRate;

odour = zeros(p.numTrials,numSamples,p.numRecordings);
oil = zeros(p.numTrials,numSamples,p.numRecordings);
for i=1:p.numRecordings
    for j=1:p.numTrials
        odour(j,:,i) = data{1,i}(onsets_odour(j,i)-p.preSeconds*p.samplingRate/p.downsampleFactor:onsets_odour(j,i)+p.postSeconds*p.samplingRate/p.downsampleFactor);
        oil(j,:,i) = data{1,i}(onsets_oil(j,i)-p.preSeconds*p.samplingRate/p.downsampleFactor:onsets_oil(j,i)+p.postSeconds*p.samplingRate/p.downsampleFactor);
    end
end

if p.doBaselineSubtraction
    for i=1:p.numRecordings
        for j=1:p.numTrials
            odour(j,:,i) = odour(j,:,i) - median(odour(j,1:p.baselineSubtractionWindow*p.samplingRate/p.downsampleFactor,i));
            oil(j,:,i) = oil(j,:,i) - median(oil(j,1:p.baselineSubtractionWindow*p.samplingRate/p.downsampleFactor,i));
        end
    end
end

odour_diff = odour-oil;
odour_diff_all = [];
for i=1:p.numRecordings
    odour_diff_all = [odour_diff_all; odour_diff(:,:,i)];
end

% %%% MB20201118
% odour_diff = zeros(p.numTrials,numSamples,p.numRecordings);
% for i=1:p.numRecordings
%     for j=1:p.numTrials
%         if (j==1)
%             odour_diff(j,:,i) = odour(j,:,i)-oil(j,:,i);
%         else
%             odour_diff(j,:,i) = odour(j,:,i)-nanmean([oil(j-1,:,i);oil(j,:,i)],1);
%         end
%     end
% end

odour_avg = squeeze(nanmean(odour,1))';
oil_avg = squeeze(nanmean(oil,1))';
odour_diff_avg = squeeze(nanmean(odour_diff,1))';

odour_avg_avg = nanmean(odour_avg,1);
oil_avg_avg = nanmean(oil_avg,1);
odour_diff_avg_avg = nanmean(odour_diff_avg,1);

odour_diff_all_mean = nanmean(odour_diff_all,1);
odour_diff_all_std = nanstd(odour_diff_all,1);
odour_diff_all_sem = nanstd(odour_diff_all,1)/sqrt(p.numRecordings);

%%

% F = figure;
% hold on
% yline(0,'k:');
% for i=1:p.numRecordings
%     plot(t,odour_diff_avg(i,:),'color',[1,0,i/p.numRecordings])
% end
% plot(t,odour_diff_avg_avg,'k','LineWidth',2)
% title(p.name)
% savefig(F,['Z:\WIBR_Hippos\SniffinHippo\Analysis\Odourants\',p.name,'\',p.name,'.fig']);
% close

%%

% figure;
% hold on
% yline(0,'k:');
% shadedErrorBar(t,odour_diff_avg_avg,nanstd(odour_diff_avg,1)/sqrt(p.numRecordings));
% close

%%

% figure;
% hold on
% yline(0,'k:');
% for i=[1,2,3,5]%1:p.numRecordings
%     for j=1:p.numTrials
%         plot(t,odour_diff(j,:,i),'color',[1,0,i/p.numRecordings])
%     end
% end
% plot(t,odour_diff_avg_avg,'k','LineWidth',2)
% close


%% For upgrade - standard figure

odour_diff_norm_mean = odour_diff_all_mean / nanmax(odour_diff_all_mean);
odour_diff_norm_sem = odour_diff_all_sem / nanmax(odour_diff_all_mean);

F = figure;
hold on
yline(0,'k:');
h=shadedErrorBar(t,odour_diff_norm_mean,odour_diff_norm_sem);
ylim([-0.2,1.4])

% temp = ['Z:\WIBR_Hippos\SniffinHippo\Analysis\Odourants\',p.name,'\',p.name,'_upgrade.png'];
% saveas(F,temp);
% 
% 
% G = figure;
% hold on
% h=shadedErrorBar(t,odour_diff_norm_mean,odour_diff_norm_sem);
% ylim([-0.2,1.4])
% 
% temp = ['Z:\WIBR_Hippos\SniffinHippo\Analysis\Odourants\',p.name,'\',p.name,'_upgrade2.png'];
% saveas(F,temp);


%% For upgrade - ctrl figure

if ~strcmp(p.name(end-3:end),'_low')
    norm_H = nanmax(odour_diff_all_mean);
    odour_diff_norm_mean_H = odour_diff_norm_mean;
    odour_diff_norm_sem_H = odour_diff_norm_sem;
else
    odour_diff_norm_mean = odour_diff_all_mean / norm_H;
    odour_diff_norm_sem = odour_diff_all_sem / norm_H;

    F = figure;
    hold on
    yline(0,'k:');
    h=shadedErrorBar(t,odour_diff_norm_mean_H*100,odour_diff_norm_sem_H*100,'patchSaturation',0.05);
    hold on
    h=shadedErrorBar(t,odour_diff_norm_mean*100,odour_diff_norm_sem*100,'lineProps','-b','patchSaturation',0.05);
    yline(max(odour_diff_norm_mean)*100,'b:');
    hold off
    ylim([-0.005*100,0.02*100])
    ytickformat(gca, 'percentage');
    temp = ['Z:\WIBR_Hippos\SniffinHippo\Analysis\Odourants\',p.name,'\',p.name,'_upgrade_conc.png'];
    saveas(F,temp);
    
    norm_H = [];
    odour_diff_norm_mean_H = [];
    odour_diff_norm_sem_H = [];
end

end


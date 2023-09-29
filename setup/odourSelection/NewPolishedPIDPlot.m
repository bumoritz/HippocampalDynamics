function [norm_H,odour_diff_norm_mean_H,odour_diff_norm_sem_H] = PolishedPIDPlot(name,norm_H,odour_diff_norm_mean_H,odour_diff_norm_sem_H)

%% Polished PID Plot

% final for
% set iii: BF_low, PF_low

% clear;
% clc;

p.name              = 'BF_low' %name; %'MethylButyrate';

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
    temp = load(['E:\SniffinHippo\RepoO\PID_in\',p.name,'_',num2str(i),'.mat']);
    
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


%% For upgrade - standard figure

odour_diff_norm_mean = odour_diff_all_mean / nanmax(odour_diff_all_mean);
odour_diff_norm_sem = odour_diff_all_sem / nanmax(odour_diff_all_mean);

% F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
% yline(0,'k:');
% h=shadedErrorBar(t(1:10:end),odour_diff_norm_mean(1:10:end),odour_diff_norm_sem(1:10:end),'lineProps',[0,0,0]);
% xlim([min(t),max(t)])
% ylim([-0.2,1.4])
% 
% savefig(F,['E:\SniffinHippo\RepoO\PID_out\',p.name,'.fig']);
% saveas(F,['E:\SniffinHippo\RepoO\PID_out\',p.name,'.png']);
% set(gcf,'Color','none'); saveas(F,['E:\SniffinHippo\RepoO\PID_out\',p.name,'.pdf']); set(gcf,'Color',[1,1,1])


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
    %h=shadedErrorBar(t,odour_diff_norm_mean_H*100,odour_diff_norm_sem_H*100,'patchSaturation',0.05);
    h=shadedErrorBar(t,odour_diff_norm_mean_H*100,odour_diff_norm_sem_H*100,'lineProps',[0.5,0.5,0.5]);
    hold on
    %h=shadedErrorBar(t,odour_diff_norm_mean*100,odour_diff_norm_sem*100,'lineProps','-b','patchSaturation',0.05);
    h=shadedErrorBar(t,odour_diff_norm_mean*100,odour_diff_norm_sem*100,'lineProps',[0,0,1]);
    yline(max(odour_diff_norm_mean)*100,'b:');
    hold off
    ylim([-0.005*100,0.02*100])
    ytickformat(gca, 'percentage');
%     temp = ['Z:\WIBR_Hippos\SniffinHippo\Analysis\Odourants\',p.name,'\',p.name,'_upgrade_conc.png'];
%     saveas(F,temp);
    
    norm_H = [];
    odour_diff_norm_mean_H = [];
    odour_diff_norm_sem_H = [];
end


%% Paper figures

pp = get_p;
path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig1_fig\'];
save_root_png = [path.root_summary,'figures\Fig1_png\'];
save_root_pdf = [path.root_summary,'figures\Fig1_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig1_txt\'];

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;
yline(0,'k:');
h1=shadedErrorBar(t(1:10:end),odour_diff_norm_mean_H(1:10:end)*100,odour_diff_norm_sem_H(1:10:end)*100,'lineProps',pp.col.odour);
h2=shadedErrorBar(t(1:10:end),odour_diff_norm_mean(1:10:end)*100,odour_diff_norm_sem(1:10:end)*100,'lineProps',pp.col.gray);
ylim([-0.5,1.1]) % iA: ylim([-0.5,2.6]); iX: ylim([-0.25,1.1]); iiA: ylim([-0.15,0.7]); iiX: ylim([-0.2,0.85]); iiiA: ylim([-0.5,1.1]); iiiX: ylim([-0.3,0.6])
ytickformat('percentage')
%yticks([0,0.5]); yticklabels({'0%','5%'});
yticks([0,1]) %yticks([0,1,2])
ylabel({'PID signal','(rel. to peak at normal conc.)'})
%xticks([0,5])
xlabel({'Time (s)'})

this_title = 'Fig1_LocConcPID_iiiA';
savefig(F,[save_root_fig,'\',this_title,'.fig']);
saveas(F,[save_root_png,'\',this_title,'.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\',this_title,'.pdf']); set(gcf,'Color',[1,1,1])




%%


% this_data_1 = all(info.set==3&info.lowConc==0.001,:)*100;
% this_data_2 = all(info.set==3&info.lowConc==0,:)*100;
% v = bar(these_labels,[nanmean(this_data_1,1),nanmean(this_data_2,1)]);
% v.FaceColor = 'flat'; v.CData(1,:) = p.col.odour; v.CData(2,:) = mean([p.col.odour;p.col.white]); v.CData(3,:) = p.col.odour; v.CData(4,:) = p.col.gray; v.EdgeColor = 'none'; 
% yline(50,':');
% for i=1:size(this_data_1,1)
%     plot(these_labels(1:2),this_data_1(i,:),'-k','LineWidth',1)
% end
% for i=1:size(this_data_2,1)
%     plot(these_labels(3:4),this_data_2(i,:),'-k','LineWidth',1)
% end
% xlabel({'Concentration of 1^{st} odor '})
% ylim([35,100])
% yticks([50,100])
% ylabel({'Performance'})
% 
% savefig(F,[save_root_fig,'\Fig1_ResidualOdourControl_iii.fig']);
% saveas(F,[save_root_png,'\Fig1_ResidualOdourControl_iii.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_ResidualOdourControl_iii.pdf']); set(gcf,'Color',[1,1,1])
% 
% [temp1,~,temp2] = signrank(this_data_1(:,2)-50)
% [temp1,~,temp2] = signrank(this_data_2(:,2)-50)
% [temp1,~,temp2] = signrank(this_data_1(:,2)-50,this_data_2(:,2)-50)
% 
% 


end


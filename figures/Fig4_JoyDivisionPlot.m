%% Fig4_JoyDivisionPlot

% import data using Analyses_Master with ops.do_responseAnalysis = true;
% use Jobs_20220719 (ctrl) and Jobs_20220721 (seq) as examples (new PMT, good photostim responses, physiological amplitude)

% lines that are commented out need to be toggled to switch between seq and ctrl stim

save_root_fig = [path.root_summary,'figures\Fig4_fig\'];
save_root_png = [path.root_summary,'figures\Fig4_png\'];
save_root_pdf = [path.root_summary,'figures\Fig4_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig4_txt\'];


%% Preparations

% overkill but okay
[prop,~,~] = preprocessActivityMeasure(act,p.genPrep,p,thor_beh,iscell);
trials = createTrialsStruct(task,1:prop.numTrials_incl);

% smooth act
act_smoothed = smoothdata(act,2,'gaussian',3*5);

% split unbinned data into trials: act -> actt
these_trial_bounds(1,:) = prop.sync_all-length(p.general.frames_pre);
these_trial_bounds(2,:) = these_trial_bounds(1,:)+p.general.numBins*p.general.binSize-1;
for i=1:length(prop.sync_all_binned)
    prop.trial_frames_unbinned(:,i) = these_trial_bounds(1,i) : these_trial_bounds(2,i);
end
actt = reshape(act_smoothed(:,prop.trial_frames_unbinned),[],size(prop.trial_frames_unbinned,1),prop.numTrials);


%% Select data and make figure

% seq: Jobs_20220721
% - selected: A 26,17,55,(72),(25)
% ctrl: Jobs_20220719
% - selected: X 15,23,37 (26,6,42,36)

stimType = 'seq';
trialType = 'A';
this_trial = 1;

for this_trial=[26,17,55] %51:100 %[13,73, 70, 46,42,38,35,29,93,62]
try

% get data
% this_actt = actt(:,:,trials.stimuli.(trialType));
this_actt = actt(:,:,intersect(find(trg.seq(2,:)),trials.stimuli.(trialType)));
sta = this_actt(:,:,this_trial);
if strcmp(trialType,'A')
    % clusters = trg.grouping(:,1);
    clusters = trg.sequenceClusters(:,1);
elseif strcmp(trialType,'X')
    % clusters = trg.grouping(:,2);
    clusters = trg.sequenceClusters(:,2);
end

% get colours
if strcmp(stimType,'seq')
    cols = p.col.seq;
elseif strcmp(stimType,'ctrl')
    cols = p.col.ctrl;
end

F = paper_figure([0,0.5,mm2inch(1.9*34),mm2inch(1.3*34)]); hold on;
%F = paper_figure([0,0.5,mm2inch(1.9*34*3),mm2inch(1.3*34*3)]); hold on;

distanceScaling = 2.5; ymax = 52.5;
this_t = p.general.t_unbinned;
patch([interp1(this_t,1:length(this_t),0),interp1(this_t,1:length(this_t),0),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1)],...
    [0,ymax,ymax,0],mean([p.col.odour;p.col.white]),'EdgeColor','none');
patch([interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2)],...
    [0,ymax,ymax,0],mean([p.col.odour;p.col.white]),'EdgeColor','none');
patch([interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow)],...
    [0,ymax,ymax,0],mean([p.col.reward;p.col.white]),'EdgeColor','none');
for i=1:length(clusters)
    temp = sta(rmmissing(trg.idcs_targetedCells(:,clusters(i))),:);
    for j=1:size(temp,1)
        plot(1:length(p.general.t_unbinned),temp(j,:)+(length(clusters)+1-i)*distanceScaling,'Color',cols,'LineWidth',0.5)
    end
end
temp = 100:(250-100)/20:250;
for j=1:length(temp)-1
    if strcmp(stimType,'seq')
        this_j = j;
    elseif strcmp(stimType,'ctrl')
        this_actualStimTrial = find(find(trg.seq(2,:)) == trials.stimuli.(trialType)(this_trial));
        this_clusterOrder = trg.sequenceClusters(:,trg.sequenceOrder(this_actualStimTrial));
        this_j = find(clusters==this_clusterOrder(j));
    end
    line([temp(j),temp(j)],[(length(temp)-this_j-0.5)*distanceScaling,(length(temp)-this_j+0.5)*distanceScaling],'Color',p.col.photostim,'LineWidth',1,'LineStyle','-')
end
xlim([0,length(p.general.t_unbinned)])
ylim([0,ymax+1])
box off; axis off;

% save plot
savefig(F,[save_root_fig,'\Fig4_JoyDivisionPlot_',stimType,trialType,num2str(this_trial),'.fig']);
saveas(F,[save_root_png,'\Fig4_JoyDivisionPlot_',stimType,trialType,num2str(this_trial),'.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_JoyDivisionPlot_',stimType,trialType,num2str(this_trial),'.pdf']); set(gcf,'Color',[1,1,1])

%title(num2str(this_trial))
drawnow;
catch
end
end

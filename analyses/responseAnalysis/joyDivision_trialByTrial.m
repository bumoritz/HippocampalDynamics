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


%% Select data

stimType = 'seq';
trialType = 'A';
this_trial = 2;

% get data
this_actt = actt(:,:,trials.stimuli.(trialType));
sta = this_actt(:,:,this_trial);

% get colours
if strcmp(stimType,'seq')
    cols = p.col.seq_rainbow;
elseif strcmp(stimType,'ctrl')
    cols = p.col.ctrl_rainbow;
end

% get stim clusters
if strcmp(trialType,'A')
    clusters = trg.grouping(:,1);
elseif strcmp(trialType,'X')
    clusters = trg.grouping(:,2);
end


%% Figure

distanceScaling = 0.5;

F = default_figure(); hold on;

for i=1:length(clusters)
    temp = sta(rmmissing(trg.idcs_targetedCells(:,clusters(i))),:);
    for j=1:size(temp,1)
        plot(1:length(p.general.t_unbinned),temp(j,:)+(length(clusters)+1-i)*distanceScaling,'Color',cols(i,:))
    end
end
taskLines(p,info,'full','traces',true,1,0,p.general.t_unbinned,1);
temp = 100:(250-100)/20:250;
for j=1:length(temp)-1
    if strcmp(stimType,'seq')
        this_j = j;
    elseif strcmp(stimType,'ctrl')
        this_actualStimTrial = find(find(trg.seq(2,:)) == trials.stimuli.(trialType)(this_trial));
        this_clusterOrder = trg.sequenceClusters(:,trg.sequenceOrder(this_actualStimTrial));
        this_j = find(clusters==this_clusterOrder(j));
    end
    line([temp(j),temp(j)],[(length(temp)-this_j-0.5)*distanceScaling,(length(temp)-this_j+0.5)*distanceScaling],'Color',p.col.photostim,'LineWidth',3,'LineStyle','-')
end
%xlim([60,332])
%ylim([0,length(clusters)+3])
% title('STA sequence')
% ylabel('z-scored \DeltaF/F (stacked by stim cluster)')

savefig(F,'C:\SniffinHippo\Summary\plots\responseSummary\joyDivision.fig');
saveas(F,'C:\SniffinHippo\Summary\plots\responseSummary\joyDivision.png');






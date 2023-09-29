%% Preparations
% need resp, flw, task

% overkill but okay
[prop,~,~] = preprocessActivityMeasure(act,p.genPrep,p,thor_beh,resp.iscell_used);
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


%% Get data

this_actt = actt(:,:,trials.stimuli.A);
stimType = 'seq';

% select trial type
if strcmp(stimType,'seq')
    sta = nanmean(this_actt,3); %resp.avgSta_seq;
    cols = p.col.seq_rainbow;
    clusters = trg.grouping(:,1);
elseif strcmp(stimType,'ctrl')
    sta = nanmean(this_actt,3); %resp.avgSta_ctrl;
    cols = p.col.ctrl_rainbow;
    clusters = trg.grouping(:,2);
end


%% Figure

F = default_figure(); hold on;

for i=1:length(clusters)
    temp = sta(rmmissing(trg.idcs_targetedCells(:,clusters(i))),:);
    for j=1:size(temp,1)
        plot(1:length(p.general.t_unbinned),temp(j,:)+length(clusters)+1-i,'Color',cols(i,:))
    end
end
taskLines(p,info,'full','traces',true,1,0,p.general.t_unbinned)
temp = 100:(250-100)/20:250;
for j=1:length(temp)-1
    line([temp(j),temp(j)],[length(temp)-j-0.5,length(temp)-j+1],'Color',p.col.photostim,'LineWidth',3,'LineStyle','-')
end
%xlim([60,332])
ylim([0,length(clusters)+3])
% title('STA sequence')
% ylabel('z-scored \DeltaF/F (stacked by stim cluster)')

savefig(F,'C:\SniffinHippo\Summary\plots\responseSummary\joyDivision.fig');
saveas(F,'C:\SniffinHippo\Summary\plots\responseSummary\joyDivision.png');






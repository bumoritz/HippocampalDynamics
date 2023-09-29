%% Preparations

blocks = [5,8,5,8,5];
numBlocks = sum(blocks);
numBlocks_switches = 13;


%% Get average correlation within stim cluster by block and animal (averaged over clusters)

conditions = fields(d{15,4}.impr_100t{1}.avgCorr_pw);
for c=1:length(conditions)
    pwcorr_withinCluster.(conditions{c}) = nan(d_info.numAnimals,numBlocks); % [animal, block]
    
    avgAct.(conditions{c}) = nan(d_info.numAnimals,numBlocks); % [animal, block]
    
end
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if j==1
            temp0=0;
        else
            temp0 = cumsum(blocks);
            temp0 = temp0(j-1);
        end
        if ~isempty(d{i,j})
            for k=1:length(d{i,j}.impr_100t)
                this_block = temp0+k;
                
                % get performance metrics
                for c=1:length(conditions)
                    pwcorr_withinCluster.(conditions{c})(i,this_block) = nanmean(d{i,j}.impr_100t{k}.avgCorr_pw.(conditions{c}));
                    
                    avgAct.(conditions{c})(i,this_block) = nanmean(d{i,j}.impr_100t{k}.avgAct.(conditions{c}));
                    
                end
            end
        end
    end
end
% for c=1:length(conditions)
%     pwcorr_withinCluster_avg.(conditions{c}) = nanmean(pwcorr_withinCluster.(conditions{c}),1);
% end


%% Combine conditions

pwcorr_withinCluster_final = {};
pwcorr_withinCluster_final.ipsi_stim = (pwcorr_withinCluster.Aseq_Atrials_stim + pwcorr_withinCluster.Xseq_Xtrials_stim)/2;
pwcorr_withinCluster_final.contra_stim = (pwcorr_withinCluster.Aseq_Xtrials_stim + pwcorr_withinCluster.Xseq_Atrials_stim)/2;
pwcorr_withinCluster_final.none_stim = (pwcorr_withinCluster.none_Atrials_stim + pwcorr_withinCluster.none_Xtrials_stim)/2;
pwcorr_withinCluster_final.ipsi_catch = (pwcorr_withinCluster.Aseq_Atrials_catch + pwcorr_withinCluster.Xseq_Xtrials_catch)/2;
pwcorr_withinCluster_final.contra_catch = (pwcorr_withinCluster.Aseq_Xtrials_catch + pwcorr_withinCluster.Xseq_Atrials_catch)/2;
pwcorr_withinCluster_final.none_catch = (pwcorr_withinCluster.none_Atrials_catch + pwcorr_withinCluster.none_Xtrials_catch)/2;

pwcorr_withinCluster_final_byStim = {};
these_conditions = fields(pwcorr_withinCluster_final);
for c=1:length(these_conditions)
    pwcorr_withinCluster_final_byStim.([these_conditions{c},'_seq']) = nan(size(pwcorr_withinCluster_final.(these_conditions{c})));
    pwcorr_withinCluster_final_byStim.([these_conditions{c},'_ctrl']) = nan(size(pwcorr_withinCluster_final.(these_conditions{c})));
    pwcorr_withinCluster_final_byStim_avg.([these_conditions{c},'_seq']) = nan(size(pwcorr_withinCluster_final.(these_conditions{c})(:,6:18)));
    pwcorr_withinCluster_final_byStim_avg.([these_conditions{c},'_ctrl']) = nan(size(pwcorr_withinCluster_final.(these_conditions{c})(:,19:31)));
    for i=1:d_info.numAnimals
        if d_info.group(i)==7
            pwcorr_withinCluster_final_byStim.([these_conditions{c},'_seq'])(i,6:18) = pwcorr_withinCluster_final.(these_conditions{c})(i,6:18);
            pwcorr_withinCluster_final_byStim.([these_conditions{c},'_ctrl'])(i,19:31) = pwcorr_withinCluster_final.(these_conditions{c})(i,19:31);
            pwcorr_withinCluster_final_byStim_avg.([these_conditions{c},'_seq'])(i,:) = pwcorr_withinCluster_final.(these_conditions{c})(i,6:18);
            pwcorr_withinCluster_final_byStim_avg.([these_conditions{c},'_ctrl'])(i,:) = pwcorr_withinCluster_final.(these_conditions{c})(i,19:31);
        elseif d_info.group(i)==8
            pwcorr_withinCluster_final_byStim.([these_conditions{c},'_ctrl'])(i,6:18) = pwcorr_withinCluster_final.(these_conditions{c})(i,6:18);
            pwcorr_withinCluster_final_byStim.([these_conditions{c},'_seq'])(i,19:31) = pwcorr_withinCluster_final.(these_conditions{c})(i,19:31);
            pwcorr_withinCluster_final_byStim_avg.([these_conditions{c},'_ctrl'])(i,:) = pwcorr_withinCluster_final.(these_conditions{c})(i,6:18);
            pwcorr_withinCluster_final_byStim_avg.([these_conditions{c},'_seq'])(i,:) = pwcorr_withinCluster_final.(these_conditions{c})(i,19:31);
        end
    end
end
pwcorr_withinCluster_final_byStim_avg = orderfields(pwcorr_withinCluster_final_byStim_avg);



avgAct_final = {};
avgAct_final.ipsi_stim = (avgAct.Aseq_Atrials_stim + avgAct.Xseq_Xtrials_stim)/2;
avgAct_final.contra_stim = (avgAct.Aseq_Xtrials_stim + avgAct.Xseq_Atrials_stim)/2;
avgAct_final.none_stim = (avgAct.none_Atrials_stim + avgAct.none_Xtrials_stim)/2;
avgAct_final.ipsi_catch = (avgAct.Aseq_Atrials_catch + avgAct.Xseq_Xtrials_catch)/2;
avgAct_final.contra_catch = (avgAct.Aseq_Xtrials_catch + avgAct.Xseq_Atrials_catch)/2;
avgAct_final.none_catch = (avgAct.none_Atrials_catch + avgAct.none_Xtrials_catch)/2;

pwcorr_withinCluster_final_byStim = {};
these_conditions = fields(avgAct_final);
for c=1:length(these_conditions)
    avgAct_final_byStim.([these_conditions{c},'_seq']) = nan(size(avgAct_final.(these_conditions{c})));
    avgAct_final_byStim.([these_conditions{c},'_ctrl']) = nan(size(avgAct_final.(these_conditions{c})));
    avgAct_final_byStim_avg.([these_conditions{c},'_seq']) = nan(size(avgAct_final.(these_conditions{c})(:,6:18)));
    avgAct_final_byStim_avg.([these_conditions{c},'_ctrl']) = nan(size(avgAct_final.(these_conditions{c})(:,19:31)));
    for i=1:d_info.numAnimals
        if d_info.group(i)==7
            avgAct_final_byStim.([these_conditions{c},'_seq'])(i,6:18) = avgAct_final.(these_conditions{c})(i,6:18);
            avgAct_final_byStim.([these_conditions{c},'_ctrl'])(i,19:31) = avgAct_final.(these_conditions{c})(i,19:31);
            avgAct_final_byStim_avg.([these_conditions{c},'_seq'])(i,:) = avgAct_final.(these_conditions{c})(i,6:18);
            avgAct_final_byStim_avg.([these_conditions{c},'_ctrl'])(i,:) = avgAct_final.(these_conditions{c})(i,19:31);
        elseif d_info.group(i)==8
            avgAct_final_byStim.([these_conditions{c},'_ctrl'])(i,6:18) = avgAct_final.(these_conditions{c})(i,6:18);
            avgAct_final_byStim.([these_conditions{c},'_seq'])(i,19:31) = avgAct_final.(these_conditions{c})(i,19:31);
            avgAct_final_byStim_avg.([these_conditions{c},'_ctrl'])(i,:) = avgAct_final.(these_conditions{c})(i,6:18);
            avgAct_final_byStim_avg.([these_conditions{c},'_seq'])(i,:) = avgAct_final.(these_conditions{c})(i,19:31);
        end
    end
end
avgAct_final_byStim_avg = orderfields(avgAct_final_byStim_avg);


%% --- SUMMARY FIGURE --- 

% select sessions
these_animals = 'all animals'; these_sessions = ~d_info.excl & d_info.group==7; % 7 or 8 for different stim groups

% select data
this_data = pwcorr_withinCluster;

% preparations
these_rgbs = discretisedColourMap('jet',true,length(conditions));
this_n_lower = min(sum(these_sessions,1));
this_n_upper = max(sum(these_sessions,1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 1; ncols = 1; m=0;
F = default_figure();

m = m+1; subplot(nrows,ncols,m); hold on;
for c=1:length(conditions)
    temp = nan(size(these_blocks));
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0
            temp(idx,find(these_blocks(idx,:))) = this_data.(conditions{c})(idx,find(these_blocks(idx,:)));
        end
    end
    %shadedErrorBar(1:numBlocks,nanmean(temp,1),nansem(temp,1),'lineProps',these_rgbs(c,:));
    plot(1:numBlocks,nanmean(temp,1),'Color',these_rgbs(c,:))
end
temp = cumsum(blocks(1:end-1))+0.5;
for i=1:length(temp)
    if i==1 || i==3
        xline(temp(i),'k-');
    else
        xline(temp(i),'k:');
    end
end
xlim([0,numBlocks])
xlabel('Trial block')
% temp = nanmean(this_data(:,:,:),3);
% ylim([0,ceil(nanmax(temp(:))*100)])
ylim([0,ceil(nanmax(this_data.(conditions{c})(:))*100)])
ylabel('Average pairwise correlation between cells within the same stim cluster')
legend('Aseq_Atrials_catch','Aseq_Atrials_stim','Aseq_Xtrials_catch','Aseq_Xtrials_stim',...
    'Xseq_Atrials_catch','Xseq_Atrials_stim','Xseq_Xtrials_catch','Xseq_Xtrials_stim',...
    'none_Atrials_catch','none_Atrials_stim','none_Xtrials_catch','none_Xtrials_stim')
title('Seq stim first, ctrl stim second')


%% --- SUMMARY FIGURE --- bad

% select sessions
these_animals = 'all animals'; these_sessions = ~d_info.excl & d_info.group==7; % 7 or 8 for different stim groups

% select data
this_data = {};
this_data.Aseq_Atrials_catch = pwcorr_withinCluster.Aseq_Atrials_catch;
this_data.Aseq_Xtrials_catch = pwcorr_withinCluster.Aseq_Xtrials_catch;
this_data.Xseq_Xtrials_catch = pwcorr_withinCluster.Xseq_Xtrials_catch;
this_data.Xseq_Atrials_catch = pwcorr_withinCluster.Xseq_Atrials_catch;
this_data.none_Atrials_catch = pwcorr_withinCluster.none_Atrials_catch;
this_data.none_Xtrials_catch = pwcorr_withinCluster.none_Xtrials_catch;

% preparations
these_conditions = fields(this_data);
these_rgbs = discretisedColourMap('jet',true,length(these_conditions));
this_n_lower = min(sum(these_sessions,1));
this_n_upper = max(sum(these_sessions,1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 1; ncols = 1; m=0;
F = default_figure();

m = m+1; subplot(nrows,ncols,m); hold on;
for c=1:length(these_conditions)
    temp = nan(size(these_blocks));
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0
            temp(idx,find(these_blocks(idx,:))) = this_data.(these_conditions{c})(idx,find(these_blocks(idx,:)));
        end
    end
    %shadedErrorBar(1:numBlocks,nanmean(temp,1),nansem(temp,1),'lineProps',these_rgbs(c,:));
    plot(1:numBlocks,nanmean(temp,1),'Color',these_rgbs(c,:))
end
temp = cumsum(blocks(1:end-1))+0.5;
for i=1:length(temp)
    if i==1 || i==3
        xline(temp(i),'k-');
    else
        xline(temp(i),'k:');
    end
end
xlim([0,numBlocks])
xlabel('Trial block')
ylabel('Average pairwise correlation between cells within the same stim cluster')
legend('Aseq_Atrials_catch','Aseq_Xtrials_catch',...
    'Xseq_Xtrials_catch','Xseq_Atrials_catch',...
    'none_Atrials_catch','none_Xtrials_catch')
title('Seq stim first, ctrl stim second')
title('Ctrl stim first, seq stim second')



%% --- SUMMARY FIGURE --- good

% select sessions
these_animals = 'all animals'; these_sessions = ~d_info.excl & d_info.group==7; % 7 or 8 for different stim groups

% select data
this_data = pwcorr_withinCluster_final;

% preparations
these_conditions = fields(this_data);
these_rgbs = discretisedColourMap('jet',true,length(these_conditions));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];

nrows = 1; ncols = 1; m=0;
F = default_figure();

m = m+1; subplot(nrows,ncols,m); hold on;
for c=1:length(these_conditions)
    temp = nan(size(these_blocks));
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0
            temp(idx,find(these_blocks(idx,:))) = this_data.(these_conditions{c})(idx,find(these_blocks(idx,:)));
        end
    end
    %shadedErrorBar(1:numBlocks,nanmean(temp,1),nansem(temp,1),'lineProps',these_rgbs(c,:));
    plot(1:numBlocks,nanmean(temp,1),'Color',these_rgbs(c,:))
end
temp = cumsum(blocks(1:end-1))+0.5;
for i=1:length(temp)
    if i==1 || i==3
        xline(temp(i),'k-');
    else
        xline(temp(i),'k:');
    end
end
xlim([0,numBlocks])
xlabel('Trial block')
ylabel('Average pairwise correlation between cells within the same stim cluster')
legend('ipsi_stim','contra_stim','none_stim',...
    'ipsi_catch','contra_catch','none_catch')
title('Seq stim first, ctrl stim second')
%title('Ctrl stim first, seq stim second')



%% --- SUMMARY FIGURE --- best

% select sessions
these_animals = 'all animals'; these_sessions = ~d_info.excl;

% select data
this_data = avgAct_final_byStim_avg; %pwcorr_withinCluster_final_byStim_avg;
this_metric = 'Average activity'; % 'Average pairwise correlation'

% preparations
these_conditions = fields(this_data);
these_rgbs = discretisedColourMap('jet',true,length(these_conditions));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];

nrows = 3; ncols = 4; m=0;
F = default_figure();

for c=1:length(these_conditions)

    m = m+1; subplot(nrows,ncols,m); hold on;
    temp = this_data.(these_conditions{c});
    plot(repmat(1:numBlocks_switches,d_info.numAnimals,1)',temp','Color','y') %these_rgbs(c,:))
    shadedErrorBar(1:numBlocks_switches,nanmean(temp,1),nansem(temp,1),'lineProps','k') %'lineProps',these_rgbs(c,:));

    xlim([0,numBlocks_switches])
    xlabel('Trial block')
    ylabel(this_metric)
    title(strrep(these_conditions{c},'_','-'))
end










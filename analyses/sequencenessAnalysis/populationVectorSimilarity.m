function povSim = populationVectorSimilarity(prop,traces,avgTraces,tng,nem_cmpr,p)
%traces = traces_all; avgTraces = avgTraces_all; tng = tng_all;

%% Create templates

% full
% these_trialTypes = {'AB','AY','XY','XB'};
these_trialTypes = {'A_H_O','A_CR_O','X_H_O','X_CR_O'};
for i=1:length(these_trialTypes)
    
    % iscells
    povSim.dynTemplates.(['iscells_full_',these_trialTypes{i},'trials']) = nanmean(traces.(these_trialTypes{i})(find(prop.iscell==1),:,:),3);
    
    % sequence cells
    povSim.dynTemplates.(['A_full_',these_trialTypes{i},'trials']) = nanmean(traces.(these_trialTypes{i})(find(tng.passed.AW.A==1),:,:),3);
    povSim.dynTemplates.(['X_full_',these_trialTypes{i},'trials']) = nanmean(traces.(these_trialTypes{i})(find(tng.passed.AW.X==1),:,:),3);
    
    % LR and nLR cells
    povSim.dynTemplates.(['LR_full_',these_trialTypes{i},'trials']) = nanmean(traces.(these_trialTypes{i})(find(nem_cmpr.sigM_sigTG{16}==1),:,:),3);
    povSim.dynTemplates.(['nLR_full_',these_trialTypes{i},'trials']) = nanmean(traces.(these_trialTypes{i})(find(nem_cmpr.sigM_sigTG{16}==0),:,:),3);
end

% AW
% these_trialTypes = {'A','X'};
these_trialTypes = {'A_corr_O','X_corr_O'};
for i=1:length(these_trialTypes)
    
    % iscells
    povSim.dynTemplates.(['iscells_AW_',these_trialTypes{i},'trials']) = nanmean(traces.(these_trialTypes{i})(find(prop.iscell==1),p.general.bins_analysisWindow,:),3);

    % sequence cells
    povSim.dynTemplates.(['A_AW_',these_trialTypes{i},'trials']) = nanmean(traces.(these_trialTypes{i})(find(tng.passed.AW.A==1),p.general.bins_analysisWindow,:),3);
    povSim.dynTemplates.(['X_AW_',these_trialTypes{i},'trials']) = nanmean(traces.(these_trialTypes{i})(find(tng.passed.AW.X==1),p.general.bins_analysisWindow,:),3);
end

% iscells, base
% these_trialTypes = {'A','X','B','Y','AB','AY','XY','XB'};
these_trialTypes = {'A_corr_O','X_corr_O','B_corr_O','Y_corr_O','A_H_O','A_CR_O','X_H_O','X_CR_O'};
for i=1:length(these_trialTypes)
    povSim.statTemplates.(['iscells_base1s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(prop.iscell==1),p.general.bins_base1s,:),3),2);
    povSim.statTemplates.(['iscells_base2s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(prop.iscell==1),p.general.bins_base2s,:),3),2);
end

% iscells, 1st
% these_trialTypes = {'A','X'};
these_trialTypes = {'A_corr_O','X_corr_O'};
for i=1:length(these_trialTypes)
    
    % iscells
    povSim.statTemplates.(['iscells_1st1s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(prop.iscell==1),p.general.bins_1st1s,:),3),2);
    povSim.statTemplates.(['iscells_1st2s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(prop.iscell==1),p.general.bins_1st2s,:),3),2);
    povSim.statTemplates.(['iscells_1st3s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(prop.iscell==1),p.general.bins_1st3s,:),3),2);


    % sequence cells
    povSim.statTemplates.(['A_1st1s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(tng.passed.AW.A==1),p.general.bins_1st1s,:),3),2);
    povSim.statTemplates.(['A_1st2s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(tng.passed.AW.A==1),p.general.bins_1st2s,:),3),2);
    povSim.statTemplates.(['A_1st3s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(tng.passed.AW.A==1),p.general.bins_1st3s,:),3),2);
    povSim.statTemplates.(['X_1st1s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(tng.passed.AW.X==1),p.general.bins_1st1s,:),3),2);
    povSim.statTemplates.(['X_1st2s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(tng.passed.AW.X==1),p.general.bins_1st2s,:),3),2);
    povSim.statTemplates.(['X_1st3s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(tng.passed.AW.X==1),p.general.bins_1st3s,:),3),2);
end

% iscells, 2nd
% these_trialTypes = {'B','Y','AB','AY','XY','XB'};
these_trialTypes = {'B_corr_O','Y_corr_O','A_H_O','A_CR_O','X_H_O','X_CR_O'};
for i=1:length(these_trialTypes)
    povSim.statTemplates.(['iscells_2nd1s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(prop.iscell==1),p.general.bins_2nd1s,:),3),2);
    povSim.statTemplates.(['iscells_2nd2s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(prop.iscell==1),p.general.bins_2nd2s,:),3),2);
    povSim.statTemplates.(['iscells_2nd3s_',these_trialTypes{i},'trials']) = nanmean(nanmean(traces.(these_trialTypes{i})(find(prop.iscell==1),p.general.bins_2nd3s,:),3),2);
end


%% Calculate population vector similarities - dynTemplatesAuto

% iscells, full
% these_trialTypes = {'AB','AY','XY','XB'};
these_trialTypes = {'A_H_O','A_CR_O','X_H_O','X_CR_O'};
for i=1:length(these_trialTypes)
    for j=1:length(these_trialTypes)
        povSim.dynTemplatesAuto.Pearson.iscells_full.([these_trialTypes{i},'_',these_trialTypes{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_full_',these_trialTypes{i},'trials']),povSim.dynTemplates.(['iscells_full_',these_trialTypes{j},'trials']),'Pearson');
        povSim.dynTemplatesAuto.cosTheta.iscells_full.([these_trialTypes{i},'_',these_trialTypes{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_full_',these_trialTypes{i},'trials']),povSim.dynTemplates.(['iscells_full_',these_trialTypes{j},'trials']),'cosTheta');
    end
end

% iscells, AW
% these_trialTypes = {'A','X'};
these_trialTypes = {'A_corr_O','X_corr_O'};
for i=1:length(these_trialTypes)
    for j=1:length(these_trialTypes)
        povSim.dynTemplatesAuto.Pearson.iscells_AW.([these_trialTypes{i},'_',these_trialTypes{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_AW_',these_trialTypes{i},'trials']),povSim.dynTemplates.(['iscells_AW_',these_trialTypes{j},'trials']),'Pearson');
        povSim.dynTemplatesAuto.cosTheta.iscells_AW.([these_trialTypes{i},'_',these_trialTypes{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_AW_',these_trialTypes{i},'trials']),povSim.dynTemplates.(['iscells_AW_',these_trialTypes{j},'trials']),'cosTheta');
    end
end


%% Calculate population vector similarities - dynTrialwise

% % iscells, full
% these_trialTypes_1 = {'AB','AY','XY','XB'};
% these_trialTypes_2 = {'A','A_H','A_M','A_CR','A_FA','X','X_H','X_M','X_CR','X_FA'};
% for i=1:length(these_trialTypes_1)
%     for j=1:length(these_trialTypes_2)
%         povSim.dynTrialwise.Pearson.iscells_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(prop.iscell==1),:,:),'Pearson');
%         povSim.dynTrialwise.cosTheta.iscells_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(prop.iscell==1),:,:),'cosTheta');
%     end
% end
% 
% % iscells, AW
% these_trialTypes_1 = {'A','X'};
% these_trialTypes_2 = {'A','A_H','A_M','A_CR','A_FA','X','X_H','X_M','X_CR','X_FA'};
% for i=1:length(these_trialTypes_1)
%     for j=1:length(these_trialTypes_2)
%         povSim.dynTrialwise.Pearson.iscells_AW.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_AW_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_analysisWindow,:),'Pearson');
%         povSim.dynTrialwise.cosTheta.iscells_AW.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_AW_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_analysisWindow,:),'cosTheta');
%     end
% end


%% Calculate population vector similarities - statTrialwise

% % iscells, 1st
% these_trialTypes_1 = {'A','X'};
% these_trialTypes_2 = {'A','A_H','A_M','A_CR','A_FA','X','X_H','X_M','X_CR','X_FA'};
% for i=1:length(these_trialTypes_1)
%     for j=1:length(these_trialTypes_2)
%         povSim.statTrialwise.Pearson.iscells_1st1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st1s,:),2),'Pearson');
%         povSim.statTrialwise.cosTheta.iscells_1st1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st1s,:),2),'cosTheta');
%         povSim.statTrialwise.Pearson.iscells_1st2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st2s,:),2),'Pearson');
%         povSim.statTrialwise.cosTheta.iscells_1st2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st2s,:),2),'cosTheta');
%         povSim.statTrialwise.Pearson.iscells_1st3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st3s,:),2),'Pearson');
%         povSim.statTrialwise.cosTheta.iscells_1st3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st3s,:),2),'cosTheta');
%     end
% end
% 
% % iscells, 2nd
% these_trialTypes_1 = {'B','Y','AB','AY','XY','XB'};
% these_trialTypes_2 = {'A','A_H','A_M','A_CR','A_FA','X','X_H','X_M','X_CR','X_FA'};
% for i=1:length(these_trialTypes_1)
%     for j=1:length(these_trialTypes_2)
%         povSim.statTrialwise.Pearson.iscells_2nd1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd1s,:),2),'Pearson');
%         povSim.statTrialwise.cosTheta.iscells_2nd1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd1s,:),2),'cosTheta');
%         povSim.statTrialwise.Pearson.iscells_2nd2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd2s,:),2),'Pearson');
%         povSim.statTrialwise.cosTheta.iscells_2nd2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd2s,:),2),'cosTheta');
%         povSim.statTrialwise.Pearson.iscells_2nd3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd3s,:),2),'Pearson');
%         povSim.statTrialwise.cosTheta.iscells_2nd3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd3s,:),2),'cosTheta');
%     end
% end


%% Calculate population vector similarities - dynTrialwiseBalanced

% full
% these_trialTypes_1 = {'AB','AY','XY','XB'};
these_trialTypes_1 = {'A_H_O','A_CR_O','X_H_O','X_CR_O'};
these_trialTypes_2 = {'A_H_M','A_M_M','A_CR_M','A_FA_M','X_H_M','X_M_M','X_CR_M','X_FA_M'};
for i=1:length(these_trialTypes_1)
    for j=1:length(these_trialTypes_2)
        
        % iscells
        povSim.dynTrialwiseBalanced.Pearson.iscells_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(prop.iscell==1),:,:),'Pearson');
        povSim.dynTrialwiseBalanced.cosTheta.iscells_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(prop.iscell==1),:,:),'cosTheta');
    
        % sequence cells
        povSim.dynTrialwiseBalanced.Pearson.A_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['A_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(tng.passed.AW.A==1),:,:),'Pearson');
        povSim.dynTrialwiseBalanced.cosTheta.A_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['A_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(tng.passed.AW.A==1),:,:),'cosTheta');
        povSim.dynTrialwiseBalanced.Pearson.X_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['X_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(tng.passed.AW.X==1),:,:),'Pearson');
        povSim.dynTrialwiseBalanced.cosTheta.X_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['X_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(tng.passed.AW.X==1),:,:),'cosTheta');    

        % sequence cells
        povSim.dynTrialwiseBalanced.Pearson.LR_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['LR_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(nem_cmpr.sigM_sigTG{16}==1),:,:),'Pearson');
        povSim.dynTrialwiseBalanced.cosTheta.LR_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['LR_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(nem_cmpr.sigM_sigTG{16}==1),:,:),'cosTheta');
        povSim.dynTrialwiseBalanced.Pearson.nLR_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['nLR_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(nem_cmpr.sigM_sigTG{16}==0),:,:),'Pearson');
        povSim.dynTrialwiseBalanced.cosTheta.nLR_full.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['nLR_full_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(nem_cmpr.sigM_sigTG{16}==0),:,:),'cosTheta');
    end
end

% AW
% these_trialTypes_1 = {'A','X'};
these_trialTypes_1 = {'A_corr_O','X_corr_O'};
these_trialTypes_2 = {'A_H_M','A_M_M','A_CR_M','A_FA_M','X_H_M','X_M_M','X_CR_M','X_FA_M'};
for i=1:length(these_trialTypes_1)
    for j=1:length(these_trialTypes_2)
        
        % iscells
        povSim.dynTrialwiseBalanced.Pearson.iscells_AW.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_AW_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_analysisWindow,:),'Pearson');
        povSim.dynTrialwiseBalanced.cosTheta.iscells_AW.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['iscells_AW_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_analysisWindow,:),'cosTheta');
        
        % sequence cells
        povSim.dynTrialwiseBalanced.Pearson.A_AW.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['A_AW_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(tng.passed.AW.A==1),p.general.bins_analysisWindow,:),'Pearson');
        povSim.dynTrialwiseBalanced.cosTheta.A_AW.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['A_AW_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(tng.passed.AW.A==1),p.general.bins_analysisWindow,:),'cosTheta');
        povSim.dynTrialwiseBalanced.Pearson.X_AW.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['X_AW_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(tng.passed.AW.X==1),p.general.bins_analysisWindow,:),'Pearson');
        povSim.dynTrialwiseBalanced.cosTheta.X_AW.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.dynTemplates.(['X_AW_',these_trialTypes_1{i},'trials']),traces.(these_trialTypes_2{j})(find(tng.passed.AW.X==1),p.general.bins_analysisWindow,:),'cosTheta');
    end
end


%% Calculate population vector similarities - statTrialwiseBalanced

% iscells, 1st
% these_trialTypes_1 = {'A','X'};
these_trialTypes_1 = {'A_corr_O','X_corr_O'};
these_trialTypes_2 = {'A_H_M','A_M_M','A_CR_M','A_FA_M','X_H_M','X_M_M','X_CR_M','X_FA_M'};
for i=1:length(these_trialTypes_1)
    for j=1:length(these_trialTypes_2)
        
        % iscells
        povSim.statTrialwiseBalanced.Pearson.iscells_1st1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st1s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.iscells_1st1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st1s,:),2),'cosTheta');
        povSim.statTrialwiseBalanced.Pearson.iscells_1st2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st2s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.iscells_1st2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st2s,:),2),'cosTheta');
        povSim.statTrialwiseBalanced.Pearson.iscells_1st3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st3s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.iscells_1st3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_1st3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_1st3s,:),2),'cosTheta');
    
        % sequence cells
        povSim.statTrialwiseBalanced.Pearson.A_1st1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['A_1st1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.A==1),p.general.bins_1st1s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.A_1st1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['A_1st1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.A==1),p.general.bins_1st1s,:),2),'cosTheta');
        povSim.statTrialwiseBalanced.Pearson.A_1st2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['A_1st2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.A==1),p.general.bins_1st2s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.A_1st2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['A_1st2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.A==1),p.general.bins_1st2s,:),2),'cosTheta');
        povSim.statTrialwiseBalanced.Pearson.A_1st3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['A_1st3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.A==1),p.general.bins_1st3s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.A_1st3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['A_1st3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.A==1),p.general.bins_1st3s,:),2),'cosTheta');
        povSim.statTrialwiseBalanced.Pearson.X_1st1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['X_1st1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.X==1),p.general.bins_1st1s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.X_1st1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['X_1st1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.X==1),p.general.bins_1st1s,:),2),'cosTheta');
        povSim.statTrialwiseBalanced.Pearson.X_1st2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['X_1st2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.X==1),p.general.bins_1st2s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.X_1st2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['X_1st2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.X==1),p.general.bins_1st2s,:),2),'cosTheta');
        povSim.statTrialwiseBalanced.Pearson.X_1st3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['X_1st3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.X==1),p.general.bins_1st3s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.X_1st3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['X_1st3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(tng.passed.AW.X==1),p.general.bins_1st3s,:),2),'cosTheta');
    end
end

% iscells, 2nd
% these_trialTypes_1 = {'B','Y','AB','AY','XY','XB'};
these_trialTypes_1 = {'B_corr_O','Y_corr_O','A_H_O','A_CR_O','X_H_O','X_CR_O'};
these_trialTypes_2 = {'A_H_M','A_M_M','A_CR_M','A_FA_M','X_H_M','X_M_M','X_CR_M','X_FA_M'};
for i=1:length(these_trialTypes_1)
    for j=1:length(these_trialTypes_2)
        povSim.statTrialwiseBalanced.Pearson.iscells_2nd1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd1s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.iscells_2nd1s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd1s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd1s,:),2),'cosTheta');
        povSim.statTrialwiseBalanced.Pearson.iscells_2nd2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd2s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.iscells_2nd2s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd2s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd2s,:),2),'cosTheta');
        povSim.statTrialwiseBalanced.Pearson.iscells_2nd3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd3s,:),2),'Pearson');
        povSim.statTrialwiseBalanced.cosTheta.iscells_2nd3s.([these_trialTypes_1{i},'temp_',these_trialTypes_2{j}]) = calculatePopulationVectorSimilarities(povSim.statTemplates.(['iscells_2nd3s_',these_trialTypes_1{i},'trials']),nanmean(traces.(these_trialTypes_2{j})(find(prop.iscell==1),p.general.bins_2nd3s,:),2),'cosTheta');
    end
end
















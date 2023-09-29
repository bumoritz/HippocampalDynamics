function dyn = getDynamicECA(prop,traces,avgTraces,nem_cmpr,p)
%traces = traces_all; avgTraces = avgTraces_all; nem_cmpr = nem_all_cmpr;

numTestGroups = length(nem_cmpr.testGroup.label);


%% All individual test groups (+/-, +, -)

these_trialTypes = {'AB','AY','XY','XB',...
    'A_H','A_CR','X_H','X_CR',...
    'A_H_M','A_CR_M','X_H_M','X_CR_M','A_M_M','A_FA_M','X_M_M','X_FA_M'};
for k=1:length(these_trialTypes)
    dyn.(these_trialTypes{k}).sigM_sigTG = {};
    dyn.(these_trialTypes{k}).sigM_sigTG_pos = {};
    dyn.(these_trialTypes{k}).sigM_sigTG_neg = {};
    for i=1:numTestGroups
        dyn.(these_trialTypes{k}).sigM_sigTG{i} = nanmean(avgTraces.(these_trialTypes{k})(find(nem_cmpr.sigM_sigTG{i}==1),:),1);
        dyn.(these_trialTypes{k}).sigM_sigTG_pos{i} = nanmean(avgTraces.(these_trialTypes{k})(find(nem_cmpr.sigM_sigTG_pos{i}==1),:),1);
        dyn.(these_trialTypes{k}).sigM_sigTG_neg{i} = nanmean(avgTraces.(these_trialTypes{k})(find(nem_cmpr.sigM_sigTG_neg{i}==1),:),1);
    end
end


%% Specific combinations of test groups

these_trialTypes = {'AB','AY','XY','XB',...
    'A_H','A_CR','X_H','X_CR',...
    'A_H_M','A_CR_M','X_H_M','X_CR_M','A_M_M','A_FA_M','X_M_M','X_FA_M'};

dyn.(these_trialTypes{k}).specificCombs_labels = {'A+_nLR','X+_nLR','B+_nLR','Y+_nLR','A_{sens}+_nLR','X_{sens}+_nLR',...
    'int+_nLR','int-_nLR','A_{int}+_nLR','A_{int}-_nLR','X_{int}+_nLR','X_{int}-_nLR',...
    'R+_nL','R-_nL',...
    'A+_nX+','X+_nA+','B+_nY+','Y+_B+','A_{sens}+_nX_{sens}+','X_{sens}+_nA_{sens}+'};

for k=1:length(these_trialTypes)
    dyn.(these_trialTypes{k}).specificCombs = {};
    dyn.(these_trialTypes{k}).specificCombs{1} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{2}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{2} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{3}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{3} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{11}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{4} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{12}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{5} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{6}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{6} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{8}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);

    dyn.(these_trialTypes{k}).specificCombs{7} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{13}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{8} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_neg{13}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{9} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{14}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{10} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_neg{14}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{11} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{15}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{12} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_neg{15}==1),find(nem_cmpr.sigM_sigTG{16}==1)),:),1);

    dyn.(these_trialTypes{k}).specificCombs{13} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{17}==1),find(nem_cmpr.sigM_sigTG{18}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{14} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_neg{17}==1),find(nem_cmpr.sigM_sigTG{18}==1)),:),1);

    dyn.(these_trialTypes{k}).specificCombs{15} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{2}==1),find(nem_cmpr.sigM_sigTG_pos{3}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{16} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{3}==1),find(nem_cmpr.sigM_sigTG_pos{2}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{17} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{11}==1),find(nem_cmpr.sigM_sigTG_pos{12}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{18} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{12}==1),find(nem_cmpr.sigM_sigTG_pos{11}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{19} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{6}==1),find(nem_cmpr.sigM_sigTG_pos{8}==1)),:),1);
    dyn.(these_trialTypes{k}).specificCombs{20} = nanmean(avgTraces.(these_trialTypes{k})(setdiff(find(nem_cmpr.sigM_sigTG_pos{8}==1),find(nem_cmpr.sigM_sigTG_pos{6}==1)),:),1);
end

%% Return

end
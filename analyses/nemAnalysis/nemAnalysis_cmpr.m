function nemAnalysis_cmpr(path,nem_all,nem_100t)

%% nem_all

nem_all_cmpr.testGroup = nem_all.testGroup;
nem_all_cmpr.sigM = nem_all.shuffled.significant;
for i=1:nem_all.numTestGroups
    nem_all_cmpr.sigTG{i} = nem_all.testGroupShuffled{i}.significant;
    nem_all_cmpr.sigM_sigTG{i} = floor((nem_all_cmpr.sigM + nem_all.testGroupShuffled{i}.significant)/2);
    
    these_predictors = cell2mat(nem_all.predictorGroups.idcs(nem_all.testGroup.group{i}));
    this_influence = nanmean(nem_all.testGroupResidual{i}.coefs(:,these_predictors+1),2);
    nem_all_cmpr.sigTG_pos{i} = floor((nem_all_cmpr.sigTG{i} + (this_influence>0))/2);
    nem_all_cmpr.sigTG_neg{i} = floor((nem_all_cmpr.sigTG{i} + (this_influence<0))/2);
    nem_all_cmpr.sigM_sigTG_pos{i} = floor((nem_all_cmpr.sigM_sigTG{i} + (this_influence>0))/2);
    nem_all_cmpr.sigM_sigTG_neg{i} = floor((nem_all_cmpr.sigM_sigTG{i} + (this_influence<0))/2);
end
nem_all_cmpr = orderfields(nem_all_cmpr);


%% nem_100t

nem_100t_cmpr{1}.testGroup = nem_100t{1}.testGroup;
for j=1:length(nem_100t)
    nem_100t_cmpr{j}.sigM = nem_100t{j}.shuffled.significant;
    for i=1:nem_all.numTestGroups
        nem_100t_cmpr{j}.sigTG{i} = nem_100t{j}.testGroupShuffled{i}.significant;
        nem_100t_cmpr{j}.sigM_sigTG{i} = floor((nem_100t_cmpr{j}.sigM + nem_100t{j}.testGroupShuffled{i}.significant)/2);
        
        these_predictors = cell2mat(nem_100t{j}.predictorGroups.idcs(nem_100t{j}.testGroup.group{i}));
        this_influence = nanmean(nem_100t{j}.testGroupResidual{i}.coefs(:,these_predictors+1),2);
        nem_100t_cmpr{j}.sigTG_pos{i} = floor((nem_100t_cmpr{j}.sigTG{i} + (this_influence>0))/2);
        nem_100t_cmpr{j}.sigTG_neg{i} = floor((nem_100t_cmpr{j}.sigTG{i} + (this_influence<0))/2);
        nem_100t_cmpr{j}.sigM_sigTG_pos{i} = floor((nem_100t_cmpr{j}.sigM_sigTG{i} + (this_influence>0))/2);
        nem_100t_cmpr{j}.sigM_sigTG_neg{i} = floor((nem_100t_cmpr{j}.sigM_sigTG{i} + (this_influence<0))/2);
    end
    nem_100t_cmpr{j} = orderfields(nem_100t_cmpr{j});
end


%% Save and return

save([path.filepart_out,'nem_all_cmpr.mat'],'nem_all_cmpr','-v7.3');
disp(['--- Saved nem_all_cmpr file as ',[path.filepart_out,'nem_all_cmpr.mat'],'.'])

save([path.filepart_out,'nem_100t_cmpr.mat'],'nem_100t_cmpr','-v7.3');
disp(['--- Saved nem_100t_cmpr file as ',[path.filepart_out,'nem_100t_cmpr.mat'],'.'])

end
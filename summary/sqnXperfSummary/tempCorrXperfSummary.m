%function respXperfSummary(d_info,d,ops,p,path)

this_corrType = 'maxBinCorr'; % 'rankCorr', 'maxBinCorr'
this_corrMetric = 'Spearman'; % 'Pearson', 'Spearman', 'Kendall'


%% Extract data

% tempCorr
tempCorr.A_Atemp.AB = {}; tempCorr.A_Atemp.AY = {}; tempCorr.A_Atemp.XY = {}; tempCorr.A_Atemp.XB = {};
tempCorr.A_Xtemp.AB = {}; tempCorr.A_Xtemp.AY = {}; tempCorr.A_Xtemp.XY = {}; tempCorr.A_Xtemp.XB = {}; 
tempCorr.X_Atemp.AB = {}; tempCorr.X_Atemp.AY = {}; tempCorr.X_Atemp.XY = {}; tempCorr.X_Atemp.XB = {};
tempCorr.X_Xtemp.AB = {}; tempCorr.X_Xtemp.AY = {}; tempCorr.X_Xtemp.XY = {}; tempCorr.X_Xtemp.XB = {}; 
tempCorr.Aonly_Atemp.AB = {}; tempCorr.Aonly_Atemp.AY = {}; tempCorr.Aonly_Atemp.XY = {}; tempCorr.Aonly_Atemp.XB = {};
tempCorr.Aonly_Xtemp.AB = {}; tempCorr.Aonly_Xtemp.AY = {}; tempCorr.Aonly_Xtemp.XY = {}; tempCorr.Aonly_Xtemp.XB = {}; 
tempCorr.Xonly_Atemp.AB = {}; tempCorr.Xonly_Atemp.AY = {}; tempCorr.Xonly_Atemp.XY = {}; tempCorr.Xonly_Atemp.XB = {};
tempCorr.Xonly_Xtemp.AB = {}; tempCorr.Xonly_Xtemp.AY = {}; tempCorr.Xonly_Xtemp.XY = {}; tempCorr.Xonly_Xtemp.XB = {}; 
tempCorr.iscells_Atemp.AB = {}; tempCorr.iscells_Atemp.AY = {}; tempCorr.iscells_Atemp.XY = {}; tempCorr.iscells_Atemp.XB = {};
tempCorr.iscells_Xtemp.AB = {}; tempCorr.iscells_Xtemp.AY = {}; tempCorr.iscells_Xtemp.XY = {}; tempCorr.iscells_Xtemp.XB = {}; 
for i=1:d_info.numAnimals
    try    
        % A_Atemp
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="B"));
        tempCorr.A_Atemp.AB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).A_Atemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.A_Atemp.AY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).A_Atemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.A_Atemp.XY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).A_Atemp_Xtrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="B"));
        tempCorr.A_Atemp.XB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).A_Atemp_Xtrials(temp);
        
        % A_Xtemp
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="B"));
        tempCorr.A_Xtemp.AB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).A_Xtemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.A_Xtemp.AY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).A_Xtemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.A_Xtemp.XY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).A_Xtemp_Xtrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="B"));
        tempCorr.A_Xtemp.XB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).A_Xtemp_Xtrials(temp);
        
        % X_Atemp
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="B"));
        tempCorr.X_Atemp.AB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).X_Atemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.X_Atemp.AY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).X_Atemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.X_Atemp.XY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).X_Atemp_Xtrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="B"));
        tempCorr.X_Atemp.XB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).X_Atemp_Xtrials(temp);
        
        % X_Xtemp
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="B"));
        tempCorr.X_Xtemp.AB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).X_Xtemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.X_Xtemp.AY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).X_Xtemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.X_Xtemp.XY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).X_Xtemp_Xtrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="B"));
        tempCorr.X_Xtemp.XB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).X_Xtemp_Xtrials(temp);
        
        % Aonly_Atemp
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="B"));
        tempCorr.Aonly_Atemp.AB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Aonly_Atemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.Aonly_Atemp.AY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Aonly_Atemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.Aonly_Atemp.XY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Aonly_Atemp_Xtrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="B"));
        tempCorr.Aonly_Atemp.XB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Aonly_Atemp_Xtrials(temp);
        
        % Aonly_Xtemp
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="B"));
        tempCorr.Aonly_Xtemp.AB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Aonly_Xtemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.Aonly_Xtemp.AY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Aonly_Xtemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.Aonly_Xtemp.XY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Aonly_Xtemp_Xtrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="B"));
        tempCorr.Aonly_Xtemp.XB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Aonly_Xtemp_Xtrials(temp);
        
        % Xonly_Atemp
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="B"));
        tempCorr.Xonly_Atemp.AB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Xonly_Atemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.Xonly_Atemp.AY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Xonly_Atemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.Xonly_Atemp.XY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Xonly_Atemp_Xtrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="B"));
        tempCorr.Xonly_Atemp.XB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Xonly_Atemp_Xtrials(temp);
        
        % Xonly_Xtemp
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="B"));
        tempCorr.Xonly_Xtemp.AB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Xonly_Xtemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.Xonly_Xtemp.AY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Xonly_Xtemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.Xonly_Xtemp.XY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Xonly_Xtemp_Xtrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="B"));
        tempCorr.Xonly_Xtemp.XB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).Xonly_Xtemp_Xtrials(temp);
        
        % iscells_Atemp
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="B"));
        tempCorr.iscells_Atemp.AB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).iscells_Atemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.iscells_Atemp.AY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).iscells_Atemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.iscells_Atemp.XY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).iscells_Atemp_Xtrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="B"));
        tempCorr.iscells_Atemp.XB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).iscells_Atemp_Xtrials(temp);
        
        % iscells_Xtemp
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="B"));
        tempCorr.iscells_Xtemp.AB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).iscells_Xtemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="A"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.iscells_Xtemp.AY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).iscells_Xtemp_Atrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="Y"));
        tempCorr.iscells_Xtemp.XY{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).iscells_Xtemp_Xtrials(temp);
        [~,temp,~] = intersect(find(d{i,1}.task.odour1=="X"),find(d{i,1}.task.odour2=="B"));
        tempCorr.iscells_Xtemp.XB{i,1} = d{i,1}.sqn_all.tempCorr.(this_corrType).(this_corrMetric).iscells_Xtemp_Xtrials(temp);
    catch
    end
end

% correct
correct.AB = {}; correct.AY = {}; correct.XY = {}; correct.XB = {}; 
for i=1:d_info.numAnimals
    try
        temp =  d{i,1}.task.response';
        correct.AB{i,1} = temp(find(d{i,1}.task.odour1=="A"&d{i,1}.task.odour2=="B"));
        correct.AB{i,1} = correct.AB{i,1} == "H";
        correct.AY{i,1} = temp(find(d{i,1}.task.odour1=="A"&d{i,1}.task.odour2=="Y"));
        correct.AY{i,1} = correct.AY{i,1} == "CR";
        correct.XY{i,1} = temp(find(d{i,1}.task.odour1=="X"&d{i,1}.task.odour2=="Y"));
        correct.XY{i,1} = correct.XY{i,1} == "H";
        correct.XB{i,1} = temp(find(d{i,1}.task.odour1=="X"&d{i,1}.task.odour2=="B"));
        correct.XB{i,1} = correct.XB{i,1} == "CR";
    catch
    end
end


%% Calculate correlations (tempCorr x correct)

% tempCorrXcorrect
tempCorrXcorrect.A_Atemp.AB = nan(d_info.numAnimals,1); tempCorrXcorrect.A_Atemp.AY = nan(d_info.numAnimals,1); tempCorrXcorrect.A_Atemp.XY = nan(d_info.numAnimals,1); tempCorrXcorrect.A_Atemp.XB = nan(d_info.numAnimals,1);
tempCorrXcorrect.A_Xtemp.AB = nan(d_info.numAnimals,1); tempCorrXcorrect.A_Xtemp.AY = nan(d_info.numAnimals,1); tempCorrXcorrect.A_Xtemp.XY = nan(d_info.numAnimals,1); tempCorrXcorrect.A_Xtemp.XB = nan(d_info.numAnimals,1); 
tempCorrXcorrect.X_Atemp.AB = nan(d_info.numAnimals,1); tempCorrXcorrect.X_Atemp.AY = nan(d_info.numAnimals,1); tempCorrXcorrect.X_Atemp.XY = nan(d_info.numAnimals,1); tempCorrXcorrect.X_Atemp.XB = nan(d_info.numAnimals,1);
tempCorrXcorrect.X_Xtemp.AB = nan(d_info.numAnimals,1); tempCorrXcorrect.X_Xtemp.AY = nan(d_info.numAnimals,1); tempCorrXcorrect.X_Xtemp.XY = nan(d_info.numAnimals,1); tempCorrXcorrect.X_Xtemp.XB = nan(d_info.numAnimals,1); 
tempCorrXcorrect.Aonly_Atemp.AB = nan(d_info.numAnimals,1); tempCorrXcorrect.Aonly_Atemp.AY = nan(d_info.numAnimals,1); tempCorrXcorrect.Aonly_Atemp.XY = nan(d_info.numAnimals,1); tempCorrXcorrect.Aonly_Atemp.XB = nan(d_info.numAnimals,1);
tempCorrXcorrect.Aonly_Xtemp.AB = nan(d_info.numAnimals,1); tempCorrXcorrect.Aonly_Xtemp.AY = nan(d_info.numAnimals,1); tempCorrXcorrect.Aonly_Xtemp.XY = nan(d_info.numAnimals,1); tempCorrXcorrect.Aonly_Xtemp.XB = nan(d_info.numAnimals,1); 
tempCorrXcorrect.Xonly_Atemp.AB = nan(d_info.numAnimals,1); tempCorrXcorrect.Xonly_Atemp.AY = nan(d_info.numAnimals,1); tempCorrXcorrect.Xonly_Atemp.XY = nan(d_info.numAnimals,1); tempCorrXcorrect.Xonly_Atemp.XB = nan(d_info.numAnimals,1);
tempCorrXcorrect.Xonly_Xtemp.AB = nan(d_info.numAnimals,1); tempCorrXcorrect.Xonly_Xtemp.AY = nan(d_info.numAnimals,1); tempCorrXcorrect.Xonly_Xtemp.XY = nan(d_info.numAnimals,1); tempCorrXcorrect.Xonly_Xtemp.XB = nan(d_info.numAnimals,1); 
tempCorrXcorrect.iscells_Atemp.AB = nan(d_info.numAnimals,1); tempCorrXcorrect.iscells_Atemp.AY = nan(d_info.numAnimals,1); tempCorrXcorrect.iscells_Atemp.XY = nan(d_info.numAnimals,1); tempCorrXcorrect.iscells_Atemp.XB = nan(d_info.numAnimals,1);
tempCorrXcorrect.iscells_Xtemp.AB = nan(d_info.numAnimals,1); tempCorrXcorrect.iscells_Xtemp.AY = nan(d_info.numAnimals,1); tempCorrXcorrect.iscells_Xtemp.XY = nan(d_info.numAnimals,1); tempCorrXcorrect.iscells_Xtemp.XB = nan(d_info.numAnimals,1); 
for i=1:d_info.numAnimals
    try
        % A_Atemp
        tempCorrXcorrect.A_Atemp.AB(i) = corr(tempCorr.A_Atemp.AB{i,1},correct.AB{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.A_Atemp.AY(i) = corr(tempCorr.A_Atemp.AY{i,1},correct.AY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.A_Atemp.XY(i) = corr(tempCorr.A_Atemp.XY{i,1},correct.XY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.A_Atemp.XB(i) = corr(tempCorr.A_Atemp.XB{i,1},correct.XB{i,1},'Type','Pearson','Rows','Complete');
        
        % A_Xtemp
        tempCorrXcorrect.A_Xtemp.AB(i) = corr(tempCorr.A_Xtemp.AB{i,1},correct.AB{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.A_Xtemp.AY(i) = corr(tempCorr.A_Xtemp.AY{i,1},correct.AY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.A_Xtemp.XY(i) = corr(tempCorr.A_Xtemp.XY{i,1},correct.XY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.A_Xtemp.XB(i) = corr(tempCorr.A_Xtemp.XB{i,1},correct.XB{i,1},'Type','Pearson','Rows','Complete');       
        
        % X_Atemp
        tempCorrXcorrect.X_Atemp.AB(i) = corr(tempCorr.X_Atemp.AB{i,1},correct.AB{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.X_Atemp.AY(i) = corr(tempCorr.X_Atemp.AY{i,1},correct.AY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.X_Atemp.XY(i) = corr(tempCorr.X_Atemp.XY{i,1},correct.XY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.X_Atemp.XB(i) = corr(tempCorr.X_Atemp.XB{i,1},correct.XB{i,1},'Type','Pearson','Rows','Complete');
        
        % X_Xtemp
        tempCorrXcorrect.X_Xtemp.AB(i) = corr(tempCorr.X_Xtemp.AB{i,1},correct.AB{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.X_Xtemp.AY(i) = corr(tempCorr.X_Xtemp.AY{i,1},correct.AY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.X_Xtemp.XY(i) = corr(tempCorr.X_Xtemp.XY{i,1},correct.XY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.X_Xtemp.XB(i) = corr(tempCorr.X_Xtemp.XB{i,1},correct.XB{i,1},'Type','Pearson','Rows','Complete');  
        
        % Aonly_Atemp
        tempCorrXcorrect.Aonly_Atemp.AB(i) = corr(tempCorr.Aonly_Atemp.AB{i,1},correct.AB{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Aonly_Atemp.AY(i) = corr(tempCorr.Aonly_Atemp.AY{i,1},correct.AY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Aonly_Atemp.XY(i) = corr(tempCorr.Aonly_Atemp.XY{i,1},correct.XY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Aonly_Atemp.XB(i) = corr(tempCorr.Aonly_Atemp.XB{i,1},correct.XB{i,1},'Type','Pearson','Rows','Complete');
        
        % Aonly_Xtemp
        tempCorrXcorrect.Aonly_Xtemp.AB(i) = corr(tempCorr.Aonly_Xtemp.AB{i,1},correct.AB{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Aonly_Xtemp.AY(i) = corr(tempCorr.Aonly_Xtemp.AY{i,1},correct.AY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Aonly_Xtemp.XY(i) = corr(tempCorr.Aonly_Xtemp.XY{i,1},correct.XY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Aonly_Xtemp.XB(i) = corr(tempCorr.Aonly_Xtemp.XB{i,1},correct.XB{i,1},'Type','Pearson','Rows','Complete');       
        
        % Xonly_Atemp
        tempCorrXcorrect.Xonly_Atemp.AB(i) = corr(tempCorr.Xonly_Atemp.AB{i,1},correct.AB{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Xonly_Atemp.AY(i) = corr(tempCorr.Xonly_Atemp.AY{i,1},correct.AY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Xonly_Atemp.XY(i) = corr(tempCorr.Xonly_Atemp.XY{i,1},correct.XY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Xonly_Atemp.XB(i) = corr(tempCorr.Xonly_Atemp.XB{i,1},correct.XB{i,1},'Type','Pearson','Rows','Complete');
        
        % Xonly_Xtemp
        tempCorrXcorrect.Xonly_Xtemp.AB(i) = corr(tempCorr.Xonly_Xtemp.AB{i,1},correct.AB{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Xonly_Xtemp.AY(i) = corr(tempCorr.Xonly_Xtemp.AY{i,1},correct.AY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Xonly_Xtemp.XY(i) = corr(tempCorr.Xonly_Xtemp.XY{i,1},correct.XY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.Xonly_Xtemp.XB(i) = corr(tempCorr.Xonly_Xtemp.XB{i,1},correct.XB{i,1},'Type','Pearson','Rows','Complete');  
        
        % iscells_Atemp
        tempCorrXcorrect.iscells_Atemp.AB(i) = corr(tempCorr.iscells_Atemp.AB{i,1},correct.AB{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.iscells_Atemp.AY(i) = corr(tempCorr.iscells_Atemp.AY{i,1},correct.AY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.iscells_Atemp.XY(i) = corr(tempCorr.iscells_Atemp.XY{i,1},correct.XY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.iscells_Atemp.XB(i) = corr(tempCorr.iscells_Atemp.XB{i,1},correct.XB{i,1},'Type','Pearson','Rows','Complete');
        
        % iscells_Xtemp
        tempCorrXcorrect.iscells_Xtemp.AB(i) = corr(tempCorr.iscells_Xtemp.AB{i,1},correct.AB{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.iscells_Xtemp.AY(i) = corr(tempCorr.iscells_Xtemp.AY{i,1},correct.AY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.iscells_Xtemp.XY(i) = corr(tempCorr.iscells_Xtemp.XY{i,1},correct.XY{i,1},'Type','Pearson','Rows','Complete');
        tempCorrXcorrect.iscells_Xtemp.XB(i) = corr(tempCorr.iscells_Xtemp.XB{i,1},correct.XB{i,1},'Type','Pearson','Rows','Complete');
    catch
    end
end



%% Figure

nrows = 2;
ncols = 5;
F = default_figure([20,0.5,20,9.9]);


% a) Template correlation of A cells with A template

these_labels = {'AB','XY','AY','XB'};
these_data = [tempCorrXcorrect.A_Atemp.AB,tempCorrXcorrect.A_Atemp.XY,tempCorrXcorrect.A_Atemp.AY,tempCorrXcorrect.A_Atemp.XB];
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

subplot(nrows,ncols,1)
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of A cells with A template')


% b) Template correlation of X cells with A template

these_labels = {'AB','XY','AY','XB'};
these_data = [tempCorrXcorrect.X_Atemp.AB,tempCorrXcorrect.X_Atemp.XY,tempCorrXcorrect.X_Atemp.AY,tempCorrXcorrect.X_Atemp.XB];
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

subplot(nrows,ncols,2)
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of X cells with A template')


% c) Template correlation of Aonly cells with A template

these_labels = {'AB','XY','AY','XB'};
these_data = [tempCorrXcorrect.Aonly_Atemp.AB,tempCorrXcorrect.Aonly_Atemp.XY,tempCorrXcorrect.Aonly_Atemp.AY,tempCorrXcorrect.Aonly_Atemp.XB];
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

subplot(nrows,ncols,3)
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of Aonly cells with A template')


% d) Template correlation of Xonly cells with A template

these_labels = {'AB','XY','AY','XB'};
these_data = [tempCorrXcorrect.Xonly_Atemp.AB,tempCorrXcorrect.Xonly_Atemp.XY,tempCorrXcorrect.Xonly_Atemp.AY,tempCorrXcorrect.Xonly_Atemp.XB];
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

subplot(nrows,ncols,4)
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of Xonly cells with A template')


% e) Template correlation of all cells with A template

these_labels = {'AB','XY','AY','XB'};
these_data = [tempCorrXcorrect.iscells_Atemp.AB,tempCorrXcorrect.iscells_Atemp.XY,tempCorrXcorrect.iscells_Atemp.AY,tempCorrXcorrect.iscells_Atemp.XB];
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

subplot(nrows,ncols,5)
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of all cells with A template')


% f) Template correlation of A cells with X template

these_labels = {'AB','XY','AY','XB'};
these_data = [tempCorrXcorrect.A_Xtemp.AB,tempCorrXcorrect.A_Xtemp.XY,tempCorrXcorrect.A_Xtemp.AY,tempCorrXcorrect.A_Xtemp.XB];
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

subplot(nrows,ncols,6)
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of A cells with X template')


% g) Template correlation of X cells with X template

these_labels = {'AB','XY','AY','XB'};
these_data = [tempCorrXcorrect.X_Xtemp.AB,tempCorrXcorrect.X_Xtemp.XY,tempCorrXcorrect.X_Xtemp.AY,tempCorrXcorrect.X_Xtemp.XB];
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

subplot(nrows,ncols,7)
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of X cells with X template')


% h) Template correlation of Aonly cells with X template

these_labels = {'AB','XY','AY','XB'};
these_data = [tempCorrXcorrect.Aonly_Xtemp.AB,tempCorrXcorrect.Aonly_Xtemp.XY,tempCorrXcorrect.Aonly_Xtemp.AY,tempCorrXcorrect.Aonly_Xtemp.XB];
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

subplot(nrows,ncols,8)
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of Aonly cells with X template')


% i) Template correlation of Xonly cells with X template

these_labels = {'AB','XY','AY','XB'};
these_data = [tempCorrXcorrect.Xonly_Xtemp.AB,tempCorrXcorrect.Xonly_Xtemp.XY,tempCorrXcorrect.Xonly_Xtemp.AY,tempCorrXcorrect.Xonly_Xtemp.XB];
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

subplot(nrows,ncols,9)
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of Xonly cells with X template')


% j) Template correlation of all cells with X template

these_labels = {'AB','XY','AY','XB'};
these_data = [tempCorrXcorrect.iscells_Xtemp.AB,tempCorrXcorrect.iscells_Xtemp.XY,tempCorrXcorrect.iscells_Xtemp.AY,tempCorrXcorrect.iscells_Xtemp.XB];
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

subplot(nrows,ncols,10)
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of all cells with X template')


%% Figure - combined into ipsi and contra

nrows = 1;
ncols = 3;
F = default_figure([20,0.5,20,9.9]);


% a) Template correlation (A or X cells)

these_labels = {'ipsi-go','contra-go','ipsi-nogo','contra-nogo'};
these_data = [nanmean([tempCorrXcorrect.A_Atemp.AB,tempCorrXcorrect.X_Xtemp.XY],2),...
    nanmean([tempCorrXcorrect.A_Atemp.XY,tempCorrXcorrect.X_Xtemp.AB],2),...
    nanmean([tempCorrXcorrect.A_Atemp.AY,tempCorrXcorrect.X_Xtemp.XB],2),...
    nanmean([tempCorrXcorrect.A_Atemp.XB,tempCorrXcorrect.X_Xtemp.AY],2)];
these_cols = {p.col.gray,p.col.gray,p.col.gray,p.col.gray};

subplot(nrows,ncols,1)
barsWithLines(these_labels,these_data,these_cols);
sigstar({[1,2],[3,4]},[signrank(these_data(:,1),these_data(:,2)),signrank(these_data(:,3),these_data(:,4))]);
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['1-2: n = ',num2str(length(rmmissing(these_data(:,1:2)))),', p = ',num2str(signrank(these_data(:,1),these_data(:,2)),2),newline,...
    '3-4: n = ',num2str(length(rmmissing(these_data(:,3:4)))),', p = ',num2str(signrank(these_data(:,3),these_data(:,4)),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation (A or X cells)')


% b) Template correlation (Aonly or Xonly cells)

these_labels = {'ipsi-go','contra-go','ipsi-nogo','contra-nogo'};
these_data = [nanmean([tempCorrXcorrect.Aonly_Atemp.AB,tempCorrXcorrect.Xonly_Xtemp.XY],2),...
    nanmean([tempCorrXcorrect.Aonly_Atemp.XY,tempCorrXcorrect.Xonly_Xtemp.AB],2),...
    nanmean([tempCorrXcorrect.Aonly_Atemp.AY,tempCorrXcorrect.Xonly_Xtemp.XB],2),...
    nanmean([tempCorrXcorrect.Aonly_Atemp.XB,tempCorrXcorrect.Xonly_Xtemp.AY],2)];
these_cols = {p.col.gray,p.col.gray,p.col.gray,p.col.gray};

subplot(nrows,ncols,2)
barsWithLines(these_labels,these_data,these_cols);
sigstar({[1,2],[3,4]},[signrank(these_data(:,1),these_data(:,2)),signrank(these_data(:,3),these_data(:,4))]);
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['1-2: n = ',num2str(length(rmmissing(these_data(:,1:2)))),', p = ',num2str(signrank(these_data(:,1),these_data(:,2)),2),newline,...
    '3-4: n = ',num2str(length(rmmissing(these_data(:,3:4)))),', p = ',num2str(signrank(these_data(:,3),these_data(:,4)),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation (Aonly or Xonly cells)')


% c) Template correlation (all cells)

these_labels = {'ipsi-go','contra-go','ipsi-nogo','contra-nogo'};
these_data = [nanmean([tempCorrXcorrect.iscells_Atemp.AB,tempCorrXcorrect.iscells_Xtemp.XY],2),...
    nanmean([tempCorrXcorrect.iscells_Atemp.XY,tempCorrXcorrect.iscells_Xtemp.AB],2),...
    nanmean([tempCorrXcorrect.iscells_Atemp.AY,tempCorrXcorrect.iscells_Xtemp.XB],2),...
    nanmean([tempCorrXcorrect.iscells_Atemp.XB,tempCorrXcorrect.iscells_Xtemp.AY],2)];
these_cols = {p.col.gray,p.col.gray,p.col.gray,p.col.gray};

subplot(nrows,ncols,3)
barsWithLines(these_labels,these_data,these_cols);
sigstar({[1,2],[3,4]},[signrank(these_data(:,1),these_data(:,2)),signrank(these_data(:,3),these_data(:,4))]);
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['1-2: n = ',num2str(length(rmmissing(these_data(:,1:2)))),', p = ',num2str(signrank(these_data(:,1),these_data(:,2)),2),newline,...
    '3-4: n = ',num2str(length(rmmissing(these_data(:,3:4)))),', p = ',num2str(signrank(these_data(:,3),these_data(:,4)),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
ylabel('Correlation (tempCorr x %correct)')
title('Template correlation of all cells')


%% Figure - iscells, by percentile

this_quantile = 5;

nrows = 1;
ncols = 3;
F = default_figure([20,0.5,20,9.9]);


% a) Template correlation (all cells with A template)

these_labels = {'AB','XY','AY','XB'};
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

these_data_x_raw = [tempCorr.iscells_Atemp.AB,tempCorr.iscells_Atemp.XY,tempCorr.iscells_Atemp.AY,tempCorr.iscells_Atemp.XB];
these_data_y_raw = [correct.AB,correct.XY,correct.AY,correct.XB];
these_data_x = {nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile)};
these_data_y = {nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile)};
for i=1:d_info.numAnimals
    for j=1:length(these_labels)
        if ~isempty(these_data_x_raw{i,j})
            these_ranks = quantileranks(these_data_x_raw{i,j},10);
            for k=1:this_quantile
                these_data_x{j}(i,k) = nanmean(these_data_x_raw{i,j}(find(these_ranks==k)));
                these_data_y{j}(i,k) = nanmean(these_data_y_raw{i,j}(find(these_ranks==k)));
            end
        end
    end
end

subplot(nrows,ncols,1)
yline(50,':');
hold on
for j=1:length(these_labels)
    temp=shadedErrorBar(1:this_quantile,nanmean(these_data_y{j},1)*100,nansem(these_data_y{j},1)*100,'lineProps',these_cols{j}); temp.mainLine.LineWidth = 2;  
end
xlim([0,this_quantile+1])
ylim([0,100])
ytickformat('percentage')
xlabel(['Template correlation (',num2str(this_quantile),'-quantiles)'])
ylabel('Performance')
title('Template correlation (all cells with A template)')


% b) Template correlation (all cells with X template)

these_labels = {'AB','XY','AY','XB'};
these_cols = {p.col.AB,p.col.XY,p.col.AY,p.col.XB};

these_data_x_raw = [tempCorr.iscells_Xtemp.AB,tempCorr.iscells_Xtemp.XY,tempCorr.iscells_Xtemp.AY,tempCorr.iscells_Xtemp.XB];
these_data_y_raw = [correct.AB,correct.XY,correct.AY,correct.XB];
these_data_x = {nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile)};
these_data_y = {nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile)};
for i=1:d_info.numAnimals
    for j=1:length(these_labels)
        if ~isempty(these_data_x_raw{i,j})
            these_ranks = quantileranks(these_data_x_raw{i,j},10);
            for k=1:this_quantile
                these_data_x{j}(i,k) = nanmean(these_data_x_raw{i,j}(find(these_ranks==k)));
                these_data_y{j}(i,k) = nanmean(these_data_y_raw{i,j}(find(these_ranks==k)));
            end
        end
    end
end

subplot(nrows,ncols,2)
yline(50,':');
hold on
for j=1:length(these_labels)
    temp=shadedErrorBar(1:this_quantile,nanmean(these_data_y{j},1)*100,nansem(these_data_y{j},1)*100,'lineProps',these_cols{j}); temp.mainLine.LineWidth = 2;  
end
xlim([0,this_quantile+1])
ylim([0,100])
ytickformat('percentage')
xlabel(['Template correlation (',num2str(this_quantile),'-quantiles)'])
ylabel('Performance')
title('Template correlation (all cells with X template)')


% c) Template correlation (all cells)

these_labels = {'AB','XY','AY','XB'};
these_cols = {p.col.darkGray,p.col.gray,p.col.darkGray,p.col.gray};

these_data_x_raw = {};
these_data_y_raw = {};
for i=1:d_info.numAnimals
    try
        these_data_x_raw{i,1} = [tempCorr.iscells_Atemp.AB{i};tempCorr.iscells_Xtemp.XY{i}];
        these_data_x_raw{i,2} = [tempCorr.iscells_Atemp.XY{i};tempCorr.iscells_Xtemp.AB{i}];
        these_data_x_raw{i,3} = [tempCorr.iscells_Atemp.AY{i};tempCorr.iscells_Xtemp.XB{i}];
        these_data_x_raw{i,4} = [tempCorr.iscells_Atemp.XB{i};tempCorr.iscells_Xtemp.AY{i}];
        these_data_y_raw{i,1} = [correct.AB{i};correct.XY{i}];
        these_data_y_raw{i,2} = [correct.XY{i};correct.AB{i}];
        these_data_y_raw{i,3} = [correct.AY{i};correct.XB{i}];
        these_data_y_raw{i,4} = [correct.XB{i};correct.AY{i}];        
    catch
    end
end
these_data_x = {nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile)};
these_data_y = {nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile),nan(d_info.numAnimals,this_quantile)};
for i=1:d_info.numAnimals
    for j=1:length(these_labels)
        if ~isempty(these_data_x_raw{i,j})
            these_ranks = quantileranks(these_data_x_raw{i,j},10);
            for k=1:this_quantile
                these_data_x{j}(i,k) = nanmean(these_data_x_raw{i,j}(find(these_ranks==k)));
                these_data_y{j}(i,k) = nanmean(these_data_y_raw{i,j}(find(these_ranks==k)));
            end
        end
    end
end

subplot(nrows,ncols,3)
yline(50,':');
hold on
for j=1:length(these_labels)
    temp=shadedErrorBar(1:this_quantile,nanmean(these_data_y{j},1)*100,nansem(these_data_y{j},1)*100,'lineProps',these_cols{j}); temp.mainLine.LineWidth = 2;  
end
xlim([0,this_quantile+1])
ylim([0,100])
ytickformat('percentage')
xlabel(['Template correlation (',num2str(this_quantile),'-quantiles)'])
ylabel('Performance')
title('Template correlation (all cells)')





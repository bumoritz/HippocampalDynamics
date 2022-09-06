%% Select data
clear; clc;

%% Running options

ops.do_learningCurveSummary         = false;
ops.do_lickingSummary               = false;
ops.do_respXperfSummary             = false;

ops.do_sequenceCellSummary          = false;
ops.do_sqnXperfSummary              = false;
ops.do_tngEncXPerfSummary           = false;
ops.do_ecaXperfSummary              = false;
ops.do_nemSummary                   = false;
ops.do_flexibleSummary              = true;

ops.skipIncompletelyProcessed       = true;

ops.addBackInTroubleshooter         = false;
ops.runOnSubsetOfSessions           = 1; % 0 for all sessions, 1 for the 11 expert imaging pilots

ops.lcs.alternativeLickMetrics      = true; % delete
ops.lcs.runningMetrics              = true; % re-organise

ops.close_figures                   = false;

if ispc
    path.root_repo                  = 'D:\SniffinHippo\Repo\';
    path.root_repoX                 = 'E:\SniffinHippo\RepoX\'; 
    path.root_analysis              = 'C:\SniffinHippo\Analysis\';
    path.root_analysisX             = 'F:\SniffinHippo\AnalysisX\';
    path.root_summary               = 'C:\SniffinHippo\Summary\';
else
    path.root_repo                  = '/Users/Moritz/Documents/MATLAB/PhD/Repo/';
    path.root_repoX                 = '/Users/Moritz/Documents/MATLAB/PhD/RepoX/';
    path.root_analysis              = '/Users/Moritz/Documents/MATLAB/PhD/Analysis/';
    path.root_summary               = '/Users/Moritz/Documents/MATLAB/PhD/Summary/';
end


%% Load dataset Google sheet

out = summary(path,ops);


%% Extract

temp = fieldnames(out);
for i=1:length(temp)
    eval([temp{i},'=out.',temp{i},';']);
end
clear('out');





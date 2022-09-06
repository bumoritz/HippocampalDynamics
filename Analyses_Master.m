%% Select data
clear; clc;

animal = 'Arasaka'; date = '20210323';


%% Running options

% behaviour data
ops.do_lickingAnalysis              = false;

% imaging data
ops.do_imagingQualityAnalysis       = false;
ops.do_tuningAnalysis               = false;
ops.do_inhibitionAnalysis           = false;
ops.do_errorTrialTuningAnalysis     = false;
ops.do_sequencenessAnalysis         = false;
ops.do_decodingAnalysis             = true;
ops.do_nemAnalysis                  = false;
ops.do_nemAnalysis_cmpr             = false;
ops.do_ecaAnalysis                  = false;
ops.do_osipAnalysis                 = false;

% stim data
ops.do_responseAnalysis             = false;
ops.do_imprintingAnalysis           = false;
ops.do_reactivationAnalysis         = false;

% other options
ops.close_figures                   = false;
ops.resp.skip_clusterResponseMaps   = false;
ops.dec.do_allTrials                = true;
ops.eca.do_allTrials                = true;
ops.nem.do_allTrials                = true;
ops.nem.do_100t_onlyFirst           = false;
ops.nem.do_100t                     = false;
ops.tng.do_allTrials                = true;
ops.tng.do_60t_onlyFirst            = false;
ops.tng.do_60t                      = false;
ops.tng.do_100t_onlyFirst           = false;
ops.tng.do_100t                     = false;
ops.tng.do_eventWiseAnalysis        = false;
ops.tng.skip_boringSnakePlots       = true;
ops.sqn.do_allTrials                = true;
ops.sqn.do_100t                     = false;
ops.ett.do_eventWiseAnalysis        = false; % not yet implemented
ops.impr.do_allTrials               = true;
ops.impr.do_60t_onlyFirst           = false;
ops.impr.do_60t                     = false;
ops.impr.do_100t_onlyFirst          = false;
ops.impr.do_100t                    = false;
ops.impr.skip_boringSnakePlots      = true;

if ispc
    path.root_repo                  = 'D:\SniffinHippo\Repo\';
    path.root_repoX                 = 'E:\SniffinHippo\RepoX\';
    path.root_repoXX                = 'E:\SniffinHippo\RepoXX\';
    path.root_analysis              = 'C:\SniffinHippo\Analysis\';
    path.root_analysisX             = 'F:\SniffinHippo\AnalysisX\';
else
    path.root_repo                  = '/Users/Moritz/Documents/MATLAB/PhD/Repo/';
    path.root_repoX                 = '/Users/Moritz/Documents/MATLAB/PhD/RepoX/';
    path.root_analysis              = '/Users/Moritz/Documents/MATLAB/PhD/Analysis/';
end


%% Run analyses

out = analyses(animal,date,path,ops);


%% Extract

temp = fieldnames(out);
for i=1:length(temp)
    eval([temp{i},'=out.',temp{i},';']);
end
clear('out');


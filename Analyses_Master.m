%% Select data
clear; clc;

%animal = 'Ao'; date = '20220402';
animal = 'Faramir'; date = '20210608';


%% Running options

% behaviour data
ops.do_lickingAnalysis              = false;
ops.do_behConfoundAnalysis          = false;

% imaging data
ops.do_imagingQualityAnalysis       = false;
ops.do_tuningAnalysis               = true;
ops.do_saveTraces                   = true;
ops.do_placeCellAnalysis            = false;
ops.do_tuningAnalysisDff            = false;
ops.do_tuningAnalysisDffS           = false;
ops.do_inhibitionAnalysis           = false;
ops.do_errorTrialTuningAnalysis     = false;
ops.do_sequenceWarpAnalysis         = false;
ops.do_sequencenessAnalysis         = false;
ops.do_decodingAnalysis             = false;
ops.do_nemAnalysis                  = false;
ops.do_nemAnalysis_cmpr             = false;
ops.do_ecaAnalysis                  = false;
ops.do_osipAnalysis                 = false;
ops.do_trnsAnalysis                 = false;
ops.do_prepostAnalysis              = false;

% stim data
ops.do_responseAnalysis             = false;
ops.do_followerAnalysis             = false;
ops.do_followerSTAs                 = false;
ops.do_imprintingAnalysis           = false;
ops.do_reactivationAnalysis         = false;

% other options
ops.close_figures                   = false;
ops.resp.skip_clusterResponseMaps   = true;
ops.dec.do_allTrials                = true;
ops.dec.do_100t                     = false;
ops.eca.do_allTrials                = true;
ops.nem.do_allTrials                = true;
ops.nem.do_100t_onlyFirst           = false;
ops.nem.do_100t                     = true;
ops.tng.do_allTrials                = true;
ops.tng.do_allTrials_stimVersion    = true;
ops.tng.do_60t_onlyFirst            = false;
ops.tng.do_60t                      = false;
ops.tng.do_100t_onlyFirst           = false;
ops.tng.do_100t                     = false;
ops.tng.do_100t_stimVersion         = false;
ops.tng.do_eventWiseAnalysis        = false;
ops.tng.skip_boringSnakePlots       = true;
ops.pca.do_allTrials                = true;
ops.tngn.do_allTrials                = true;
ops.tngn.do_allTrials_stimVersion    = false;
ops.tngn.do_60t_onlyFirst            = false;
ops.tngn.do_60t                      = false;
ops.tngn.do_100t_onlyFirst           = false;
ops.tngn.do_100t                     = true;
ops.tngn.do_100t_stimVersion         = false;
ops.tngn.do_eventWiseAnalysis        = false;
ops.tngn.skip_boringSnakePlots       = true;
ops.tngnn.do_allTrials                = false;
ops.tngnn.do_allTrials_stimVersion    = false;
ops.tngnn.do_60t_onlyFirst            = false;
ops.tngnn.do_60t                      = false;
ops.tngnn.do_100t_onlyFirst           = false;
ops.tngnn.do_100t                     = false;
ops.tngnn.do_100t_stimVersion         = false;
ops.tngnn.do_eventWiseAnalysis        = false;
ops.tngnn.skip_boringSnakePlots       = true;
ops.sqn.do_allTrials                = false;
ops.sqn.do_100t                     = false;
ops.warp.do_allTrials               = false;
ops.warp.do_allTrials_stimVersion   = false;
ops.warp.do_100t                    = true;
ops.warp.do_100t_stimVersion        = false;
ops.warp.input                      = "tng";
ops.ett.do_eventWiseAnalysis        = false; % not yet implemented
ops.impr.do_allTrials               = true;
ops.impr.do_60t_onlyFirst           = false;
ops.impr.do_60t                     = false;
ops.impr.do_100t_onlyFirst          = false;
ops.impr.do_100t                    = false;
ops.impr.skip_boringSnakePlots      = false;

if ispc
    path.root_repo                  = 'D:\SniffinHippo\Repo\';
    path.root_repoX                 = 'E:\SniffinHippo\RepoX\';
    path.root_repoXX                = 'E:\SniffinHippo\RepoXX\';
    path.root_analysis              = 'C:\SniffinHippo\Analysis\';
    path.root_analysisX             = 'F:\SniffinHippo\AnalysisX\';
    path.root_summary               = 'C:\SniffinHippo\Summary\';
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


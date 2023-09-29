function p = get_p

%% data2repo

p.scope.neuropilSubtraction         = 0.7;
p.scope.detrendingWindow            = 60; % [s]
p.scope.baseF_prctiles              = [8,10,20,30,40,50];
p.scope.voltageThreshold            = 0.5; % [V], threshold for digitising

p.paq.voltageThreshold              = 0.5; % [V], threshold for digitising

p.cam.sniffing_bp                   = [2,12]; % [Hz] refs for 2-12: Uchida & Mainen, review by Grimaud and Murthy https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6139454/

p.s2p_new.zscoreBeforeSVD           = true;
p.s2p_new.includeOnlyIscells        = true;
p.s2p_new.holdOutStimData           = false;
p.s2p_new.zcrdSigmaThreshold        = 3;
p.s2p_new.fovCorrThreshold          = 5;
p.s2p_new.smoothingSd_preBinning   	= 0; % in frames (not in bins)
p.s2p_new.binSize                  	= 6; % number of frames per bin
p.s2p_new.smoothingSd_postBinning  	= 0; % in bins (not in frames)
p.s2p_new.zscore                 	= false;

p.halo.sigma                        = 10;
p.halo.bkgd_prctile                 = 8;
p.halo.bkgd_inset                   = 15;
p.halo.halo_inset                   = 15;
p.halo.smoothingSd_preBinning       = 0; % in frames (not in bins)
p.halo.binSize                      = 6; % number of frames per bin
p.halo.smoothingSd_postBinning      = 0; % in bins (not in frames)
p.halo.zscore                       = false;


%% repo2repo

p.trg.minCellProb                  = -1; % suite2p classifier iscell probability to include in registration process (-1 is all Rois), DEFAULT: -1
p.trg.use_parallel_processing      = true; % either true or false
p.trg.maximal_rotation             = 30; % in degrees - only relevant if 'Translations and Rotations' is used
p.trg.transformation_smoothness    = 1; % levels of non-rigid FOV transformation smoothness (range 0.5-3) %MB20210716: tested for Faramir and Carlo, non-rigid wins, best is 1, 2 is okay, 0.5 and 3 are worse
p.trg.reference_session_index      = 2; % numeric or 'middle'
p.trg.maximal_distance             = 10; % cell-pairs that are more than 12 micrometers apart are assumed to be different cells
p.trg.p_same_certainty_threshold   = 0.95; % certain cells are those with p_same>threshld or <1-threshold
p.trg.initial_registration_type    = 'Spatial correlation'; % either 'Spatial correlation', 'Centroid distance', or 'best_model_string';
p.trg.registration_approach        = 'Probabilistic'; % either 'Probabilistic' or 'Simple threshold'
p.trg.model_type                   = 'Spatial correlation'; % either 'Spatial correlation', 'Centroid distance', or 'best_model_string';
p.trg.p_same_threshold             = 0.5; % only relevant if probabilistic approach is used

% p.trck.minCellProb                  = 1000; % suite2p classifier iscell probability to include in registration process, if 1000 then only curated, DEFAULT: 1000
p.trck.use_parallel_processing      = true; % either true or false
p.trck.maximal_rotation             = 30; % in degrees - only relevant if 'Translations and Rotations' is used
p.trck.transformation_smoothness    = 1; % levels of non-rigid FOV transformation smoothness (range 0.5-3)   %MB20210716: tested for Faramir and Carlo, non-rigid wins, best is 1, 2 is okay, 0.5 and 3 are worse
p.trck.reference_session_index      = 'middle'; % numeric or 'middle'
p.trck.maximal_distance             = 10; % cell-pairs that are more than 12 micrometers apart are assumed to be different cells
p.trck.p_same_certainty_threshold   = 0.95; % certain cells are those with p_same>threshld or <1-threshold
p.trck.initial_registration_type    = 'Spatial correlation'; % either 'Spatial correlation', 'Centroid distance', or 'best_model_string';
p.trck.registration_approach        = 'Probabilistic'; % either 'Probabilistic' or 'Simple threshold'
p.trck.model_type                   = 'Spatial correlation'; % either 'Spatial correlation', 'Centroid distance', or 'best_model_string';
p.trck.p_same_threshold             = 0.5; % only relevant if probabilistic approach is used


%% analyses - pre-processing

% general pre-processing
p.genPrep.activityMeasure               = "spksn_beh";
p.genPrep.smoothingSd_preBinning        = 12; % in frames (not in bins)
p.genPrep.binSize                       = 6; % number of frames per bin
p.genPrep.smoothingSd_postBinning       = 0; % in bins (not in frames)
p.genPrep.zscore                        = true;

% specific analyses - inherits
p.tng = p.genPrep;
p.pca = p.genPrep;
p.inh = p.genPrep;
p.ett = p.genPrep;
p.iqa = p.genPrep;
p.eca = p.genPrep;
p.impr = p.genPrep;
p.sqn = p.genPrep;
p.warp = p.genPrep;
p.dec = p.genPrep;
p.nem = p.genPrep;
p.osip = p.genPrep;
p.bcon = p.genPrep;
p.ppa = p.genPrep;

% other
p.resp.activityMeasure                  = "dFFn_beh";
p.resp.zscore                           = true;
p.resp.previousClustersForExclusionZone = 12; %8;
p.resp.excludeOutsideExclusionZone      = true;
p.resp.exclusionRadius                  = 50; % in um
p.resp.excludeInSelfCluster             = false;
p.resp.minNumTrials                     = 0;
p.resp.subtractBlockwiseTuning          = false;
p.resp.subtractOverallTuning            = false;
p.resp.ampBins_edges                    = [-5.5:1:15.5];
p.resp.ampBins_x                        = [-5:1:15]';
p.resp.ampBins_minNumDataPoints         = 20;
p.impr.activityMeasure                  = "dFFn_beh";
p.trns.activityMeasure                  = "dFFn_beh";
p.ppa.activityMeasure_pre               = "spksn_pre";
p.ppa.activityMeasure_post              = "spksn_post";

% temp
if false
    p.resp.activityMeasure = "dFF_beh";
    disp('Watch out. Running with non-default parameters')
end

%% analyses - general

p.general.rgnSeed                       = 1234;
p.general.numBins                       = 67;
p.general.bins_pre                      = 1:15; % (15 bins*6/30 = 3 s)
p.general.bins_baselineWindow           = 1:10;
p.general.bins_analysisWindow           = 16:41; % (26 bins*6/30 = 5.2 s ? 5.3s)
p.general.bins_odour1Window_bl          = 1:10;
p.general.bins_odour1Window             = 16:25;
p.general.bins_odour2Window_bl          = 27:36;
p.general.bins_odour2Window             = 42:51;
% p.general.bins_odour2Window_woRW        = 42:45; % also consider e.g. in Changa (session with shorter delay between odour2 and response window)
p.general.bins_rewardWindow_bl          = 31:40;
p.general.bins_rewardWindow             = 46:55;
p.general.bins_base1s                   = 6:10;
p.general.bins_base2s                   = 1:10;
p.general.bins_1st1s                    = 16:20;
p.general.bins_1st2s                    = 16:25;
p.general.bins_1st3s                    = 16:30;
p.general.bins_2nd1s                    = 42:46;
p.general.bins_2nd2s                    = 42:51;
p.general.bins_2nd3s                    = 42:56;
p.general.binSize                       = 6; % number of frames per bin
p.general.t_binned = ([1:p.general.numBins]-length(p.general.bins_pre)-1)*p.general.binSize/30.04; % info.scope.frameRate == 30.04
p.general.numFrames                     = p.general.numBins*p.general.binSize;
p.general.frames_pre                      = 1:15*p.general.binSize;
p.general.t_unbinned = ([1:p.general.numFrames]-length(p.general.frames_pre)-1)/30.04; % info.scope.frameRate == 30.04
% p.general.smoothingSd_preBinning_velocity = 30; % in frames (not in bins)
% p.general.smoothingSd_preBinning_acceleration = 30; % in frames (not in bins)
p.general.smoothingSd_preBinning_velocity = 12;
p.general.smoothingSd_preBinning_acceleration = 12;


%% analyses - tng

% shuffling
p.tng.analysisBlockRestrictedNf     = true;
p.tng.numShuffles                   = 1000;
p.tng.shufflingThreshold            = 95;
p.tng.pfdrCorrectionForShuffle      = false;

% firing field properties
p.tng.sdAboveBaselineType           = 'sd_of_trialwise_baseline'; % 'sd_of_entire_recording' or 'sd_of_trialwise_baseline'
p.tng.activeInFiringFieldSd         = 1; % number of sd above baseline to be called active in firing field [default: 1] 3 getted rid von den prae Aktivitaeten
p.tng.firingFieldBoundaries         = 0.5; % percent over baseline. 0.5 will lead to PF width being the FWHM
p.tng.earlyVsLateThreshold          = 1.3; % s after first odour onset [default: 1]
p.tng.reliabilityThreshold          = 0.1; %0.2; % fraction of trials with activity in firing field [default: 0.2] 0.3 sieht besser aus, gerade bei activeInFiringFieldSd=3
p.tng.selectivityThreshold          = 0;
p.tng.qThreshold                    = 0.05;
p.tng.amplitudeThreshold            = 0.2;

% tuning criteria
p.tng.criteria_shuffling            = true;
p.tng.criteria_peakInWindow         = true;
p.tng.criteria_aboveBaseline        = false;
p.tng.criteria_reliability          = true; % true;
p.tng.criteria_amplitude            = false;

% tng dff
p.tngn = p.tng;
p.tngn.activityMeasure = "dFFn_beh";
p.tngn.zscore = true;
p.tngn.smoothingSd_preBinning = 9;

% tng dff sig
p.tngnn = p.tng;
p.tngnn.activityMeasure = "dFFns_beh";
p.tngnn.zscore = true;
p.tngnn.smoothingSd_preBinning = 9;


%% analyses - nem

% model predictors
p.nem.predictorStd                      = 3; % in bins
p.nem.predictorEveryXBins               = 5; % in bins (has to be uneven, predictor will be centered)

% model fitting
p.nem.family                            = 'gaussian';
p.nem.alpha                             = 0.5;
p.nem.testSet                           = 0.2; % training set = 1-p.nem.testSet (assigned trial block-wise)
p.nem.cvFolds                           = 5;
p.nem.cvLoss                            = 'deviance';
p.nem.lambda                            = 'lambda_min';

% model versions
p.nem.do_shuffled                       = true;
p.nem.do_testGroupShuffled              = true;
p.nem.do_testGroupResidual              = true;
p.nem.do_shuffle1group                  = false;
p.nem.do_leave1out                  	= false;
p.nem.do_leave1outr                  	= false;

% significance
p.nem.numShuffles                       = 20;
p.nem.modelSignificanceStd              = 2;
p.nem.predictorGroupSignificanceStd     = 2;


%% analyses - other

p.resp.numFrames                    = 3;
p.resp.qthres                       = 0.05;
p.resp.sta_framesPre                = 30;
p.resp.sta_framesPost               = 30;
p.resp.maxTargetCentroidDistance_um = 25;
p.resp.amplitudeCritiera            = [NaN,0.5,1]; % in z-scores
p.resp.mainAmplitudeCritierion      = 0.5;
p.resp.mainResponderCondition       = 'responders_0d5z_ext';

p.sqn.shuffleOrderForSameBins       = true;
p.dec.cvFold                        = 5; %10;
p.dec.trainOnCompleteSet            = true; % used to be: false;
p.ett.cvFolds                       = 10;
p.ett.numIncorrectTrialsAX          = 'same'; %'same','max'

p.warp.histBins                     = [0:0.2:5.2];

%% summary

p.excl.engagement                   = 0.2; % fraction of trials with licks to be called engaged
p.excl.running                      = 10.01; % [cm/s], min. avg. running speed to be called a runner
p.excl.presponsive                  = 0.5; % min. responders_main.proportion_RespTargeted_AllTargeted
p.lcs.smoothing_sd                  = 1;
p.lick.smoothing_sd                 = 1;
p.lick.crit_engagement              = 0.25;
p.lick.crit_bias                    = 0.6;


%% colours

% standard colours
p.col.black                         = [0,0,0]/255; % black
p.col.gray                          = [170,170,170]/255; % gray
p.col.darkGray                      = [117,114,115]/255; % darkGray
p.col.white                         = [255,255,255]/255; % white

% optics
p.col.imaging                       = [90,188,173]/255; % turquoise
p.col.photostim                     = [229,49,18]/255; % red

% task events
p.col.odour                         = [252,193,2]/255; % yellow
p.col.reward                        = [128,204,223]/255; % light blue

% stimulation types
p.col.seq                           = [136,19,80]/255; % wine red
p.col.seq_rainbow                   = [linspace(136,241,20);linspace(19,143,20);linspace(80,177,20)]'/255; % wine red (gradient to magenta)
p.col.ctrl                          = [4,88,156]/255; % blue
p.col.ctrl_rainbow                  = [linspace(4,133,20);linspace(88,209,20);linspace(156,245,20)]'/255; % blue (gradient to cyan)
p.col.noStim                        = p.col.darkGray; % darkGray    %%%[28,95,43]/255; %green
p.col.pos                           = [1,0.3,0.3];
p.col.neg                           = [0.3,0.3,1];

% % odours
% p.col.A                             = [231,126,33]/255; % orange
% p.col.AB                            = [231,126,33]/255; % orange
% p.col.AB_rainbow                    = [linspace(231,147,20);linspace(126,82,20);linspace(33,23,20)]'/255; % orange (gradient to brown)
% p.col.B                             = [147,82,23]/255; % brown
% p.col.AY                            = [147,82,23]/255; % brown
% p.col.X                             = [192,202,50]/255; % light green
% p.col.XY                            = [192,202,50]/255; % light green
% p.col.XY_rainbow                    = [linspace(192,130,20);linspace(202,120,20);linspace(50,29,20)]'/255; % light green (gradient to olive green)
% p.col.XB                            = [130,120,29]/255; % olive green
% p.col.Y                             = [130,120,29]/255; % olive green

% odours
p.col.AB                            = [194,79,30]/255; % ocre
p.col.AY                            = [106,43,17]/255; % brown
p.col.XY                            = [184,170,40]/255; % light green
p.col.XB                            = [130,120,29]/255; % dark green
p.col.A                             = p.col.AB;
p.col.X                             = p.col.XY;
p.col.B                             = p.col.AY;
p.col.Y                             = p.col.XB;

% splits
p.col.runner                        = [250,131,52]/255; % orange
p.col.nonrunner                     = [28,95,43]/255; % green
% temp=copper;
% p.col.runner                        = temp((size(temp,1)-1)/3,:); % dark brown
% p.col.nonrunner                     = temp(2*(size(temp,1)-1)/3,:); % light brown
p.col.learning                      = 'g';
p.col.nonlearning                   = 'r';

% GLM
p.col.lick                          = p.col.darkGray;
p.col.velocity                      = p.col.runner;
p.col.acceleration                  = p.col.nonrunner;
p.col.distance                      = p.col.black;
p.col.duration                      = p.col.black;
p.col.motion                        = p.col.darkGray;

end


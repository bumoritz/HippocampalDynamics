function info = get_info

info.epochs                         = {'pre','beh','post'};

info.data.numFragments              = 1;

info.task.trialsPerBlock            = 20;
info.task.trialStructure.tOdour1    = 0.3; % [s]
info.task.trialStructure.tGap       = 5; % [s]
info.task.trialStructure.tOdour2    = 0.3; % [s]
info.task.trialStructure.tRespDelay = 0.5; % [s]
info.task.trialStructure.tRespWindow= 1; % [s]
info.task.vials.vialNumber          = [1,2,3,4];
info.task.vials.role                = ["A","X","B","Y"];
info.task.contingencies.type        = [1,2,3,4];
info.task.contingencies.odour1      = ["A","X","A","X"];
info.task.contingencies.odour2      = ["B","Y","Y","B"];
info.task.contingencies.requirement = ["GO","GO","NOGO","NOGO"];
info.task.roles.odour1              = ["A","X"];
info.task.roles.odour2              = ["B","Y"];
info.task.var.catch                 = 0;
info.task.var.stim                  = 1;

info.scope.numFrames_pre            = 108000;
info.scope.numFrames_post           = 108000;
info.scope.samplingRate             = 10000; % [Hz]
info.scope.frameRate                = 30.04; % [Hz]
info.scope.fovSize_um               = 1000.78; % [um]
info.scope.fovSize_pix              = 512; % pixels
info.scope.bytesPerNumber           = 2; % pixels
info.scope.chan.frames              = 'FrameTrigger';
info.scope.chan.frames_alt          = 'FrameOut';
info.scope.chan.stimSequence        = 'StimTrigger';
info.scope.chan.stimSequence_alt    = 'BleachOut';
info.scope.chan.stimPatternComplete = 'PatternComplete';
info.scope.chan.stimIteration       = 'IterationOutput';
info.scope.chan.sync                = 'SniffinSync';
info.scope.patternsPerSequence      = 20;

info.paq.samplingRate               = 10000; % [Hz]
info.paq.rot.pulsesPerRevolution    = 1024;
info.paq.rot.wheelDiameter_Grid     = 15.8; % [cm]
info.paq.rot.wheelDiameter_Valhalla = 19.8; % [cm]
info.paq.config.lick                = 2;
info.paq.config.imaging             = 7;
info.paq.config.cue                 = 8;
info.paq.config.finalvalve          = 9;
info.paq.config.responsewindow      = 10;
info.paq.config.reward              = 11;
info.paq.config.rotA                = 15;
info.paq.config.rotB                = 16;

%


end
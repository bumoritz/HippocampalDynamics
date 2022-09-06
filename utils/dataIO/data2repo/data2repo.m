function [out] = data2repo(animal,date,expDay,stimSession,path,ops)

%% Preparations

% gather info and p
p = get_p;
info = get_info;
info.animal = animal;
info.date = date;
info.expDay = expDay;
info.stimSession = stimSession;

% identify session folders
path.folder_data = [path.root_data,info.date(1:4),'\',info.date(1:4),'-',info.date(5:6),'\',info.date(1:4),'-',info.date(5:6),'-',info.date(7:8),'\',info.animal,'\'];
path.folder_repo = [path.root_repo,info.animal,'\',info.animal,'_',info.date,'\'];
path.filepart_out = [path.folder_repo,info.animal,'_',info.date,'_'];
if ~exist(path.folder_repo,'dir')
   mkdir(path.folder_repo);
end
path.folder_repoX = [path.root_repoX,info.animal,'\',info.animal,'_',info.date,'\'];
path.filepart_outX = [path.folder_repoX,info.animal,'_',info.date,'_'];
if ~exist(path.folder_repoX,'dir')
   mkdir(path.folder_repoX);
end

if ~exist([path.filepart_out,'plots'],'dir')
    mkdir([path.filepart_out,'plots']);
end
if ~exist([path.filepart_out,'log'],'dir')
    mkdir([path.filepart_out,'log']);
end

% start log
t_start = tic;
log.executionDate = datetime(now,'ConvertFrom','datenum');
log.done = false;
temp = [path.filepart_out,'log\',info.animal,'_',info.date,'_',datestr(log.executionDate,'yyyymmdd_HHMMSS'),'.log'];
diary(temp);
diary on;

disp([log.executionDate])
disp(['Running data2repo.m for ',info.animal,'_',info.date,'.'])

% get task information
if info.expDay==1
    info.task.numTrials = 500;
    info.task.odourSet = 'i';
    info.task.sessionType = 'expert';
elseif info.expDay==2
    info.task.numTrials = 800;
    info.task.numStimTrials = 640;
    info.task.odourSet = 'ii';
    info.task.sessionType = 'switch';
elseif info.expDay==3
    info.task.numTrials = 500;
    info.task.numStimTrials = 400;
    info.task.odourSet = 'ii';
    info.task.sessionType = 'expert';
elseif info.expDay==4
    info.task.numTrials = 800;
    info.task.numStimTrials = 640;
    info.task.odourSet = 'iii';
    info.task.sessionType = 'switch';
elseif info.expDay==5
    info.task.numTrials = 500;
    info.task.numStimTrials = 400;
    info.task.odourSet = 'iii';
    info.task.sessionType = 'expert';
elseif info.expDay==101
    info.task.numTrials_stim = 300;
    info.task.trialsPerBlock_stim = 20;
    info.task.odourSet = 'i';
    info.task.sessionType = 'ExpertStim1';
elseif info.expDay==102
    info.task.numTrials_stim = 300;
    info.task.trialsPerBlock_stim = 20;
    info.task.odourSet = 'i';
    info.task.sessionType = 'ExpertStim2';
elseif info.expDay==103
    info.task.numTrials_stim = 448;
    info.task.trialsPerBlock_stim = 56;
    info.task.odourSet = 'i';
    info.task.sessionType = 'ExpertStim3';
elseif info.expDay==104
    info.task.numTrials_stim = 360;
    info.task.trialsPerBlock_stim = 24;
    info.task.odourSet = 'i';
    info.task.sessionType = 'ExpertStim4';
elseif info.expDay==105
    info.task.numTrials_stim = 360;
    info.task.trialsPerBlock_stim = 24;
    info.task.odourSet = 'i';
    info.task.sessionType = 'ExpertStim5';
elseif info.expDay==106
    info.task.numTrials_stim = 300;
    info.task.trialsPerBlock_stim = 20;
    info.task.odourSet = 'i';
    info.task.sessionType = 'ExpertStim6';
elseif info.expDay==107
    info.task.numTrials_stim = 300;
    info.task.trialsPerBlock_stim = 20;
    info.task.odourSet = 'i';
    info.task.sessionType = 'ExpertStim7';
end
if info.expDay>100
    info.epochs = {'base1','base2','stim'};
    info.task.numTrials_base1 = 60;
    info.task.numTrials_base2 = 60;
    info.task.numTrials = info.task.numTrials_stim; 
    info.task.trialsPerBlock_base1 = 20;
    info.task.trialsPerBlock_base2 = 20;
    info.task.trialsPerBlock = info.task.trialsPerBlock_stim;
    info.task.numBlocks_base1 = info.task.numTrials_base1/info.task.trialsPerBlock_base1;
    info.task.numBlocks_base2 = info.task.numTrials_base2/info.task.trialsPerBlock_base2;
    info.task.numBlocks_stim = info.task.numTrials_stim/info.task.trialsPerBlock_stim;
end
info.task.numBlocks = info.task.numTrials/info.task.trialsPerBlock;
if info.expDay>=102
    info.task.vials.vialNumber          = [1,2,3,4,5,6];
    info.task.vials.role                = ["A","X","B","Y","O","E"];
    info.task.contingencies.type        = [1,2,3,4,5,6,7,8];
    info.task.contingencies.odour1      = ["A","X","A","X","O","O","O","O"];
    info.task.contingencies.odour2      = ["B","Y","Y","B","B","Y","Y","B"];
    info.task.contingencies.requirement = ["GO","GO","NOGO","NOGO","GO","GO","NOGO","NOGO"];
    info.task.roles.odour1              = ["A","X","O"];
end
if strcmp(info.task.odourSet,'i')
    info.task.vials.odourant = ["MethylButyrate","EthylPropionate","Pinene","Benzaldehyde"];
elseif strcmp(info.task.odourSet,'ii')
    info.task.vials.odourant = ["Pentanone","Hexanone","Eucalyptol","Limonene"];
elseif strcmp(info.task.odourSet,'iii')
    info.task.vials.odourant = ["ButylFormate","PropylFormate","Guaiacol","IsobutylPropionate"];
end


%% Troubleshooting

% Expert stim experiment specific things
if str2num(info.date)==20220116
    info.task.numTrials_stim = 400;
    info.task.numTrials = info.task.numTrials_stim; 
    info.task.numBlocks_stim = info.task.numTrials_stim/info.task.trialsPerBlock_stim;
    info.task.numBlocks = info.task.numTrials/info.task.trialsPerBlock;
end

% First two animals
if strcmp(info.animal,'Lancelot') || ...
        strcmp(info.animal,'Changa')
    info.task.numTrials = 500;
    info.task.numStimTrials = 400;
    info.task.numBlocks = info.task.numTrials/info.task.trialsPerBlock;
    info.task.trialStructure.tRespDelay = 0;
end

% Thorsync frame rate
if str2num(info.date)>20210501
    info.scope.samplingRate = 20000;
end
if strcmp(info.animal,'Lancelot') || ...
        (strcmp(info.animal,'Changa') && strcmp(info.date,'20201118'))
    info.scope.samplingRate = 30000;
end
if (strcmp(info.animal,'Changa') && strcmp(info.date,'20201116')) || ...
        (strcmp(info.animal,'Changa') && strcmp(info.date,'20201117'))
    info.scope.samplingRate = 100000;
end

% Valhalla vs Grid
if strcmp(info.animal,'Pfizer') || ...
        strcmp(info.animal,'Oxford') || ...
        strcmp(info.animal,'Asakusa') || ...
        strcmp(info.animal,'Arwen') || ...
        strcmp(info.animal,'MamTor') || ...
        strcmp(info.animal,'Diana') || ...
        strcmp(info.animal,'Java')
    info.task.room = 'Grid';
    info.paq.rot.wheelCircumference = pi*info.paq.rot.wheelDiameter_Grid; % [cm]
else
    info.task.room = 'Valhalla';
    info.paq.rot.wheelCircumference = pi*info.paq.rot.wheelDiameter_Valhalla; % [cm]
end

% File fragmentations
if (strcmp(info.animal,'Carlo') && strcmp(info.date,'20210314'))
    info.data.numFragments = 3;
    info.data.fragments_numTrials = [180,100,520];
    info.data.fragments_seq = [1:180,1:100,1:520];
end
if (strcmp(info.animal,'Carlo') && strcmp(info.date,'20210316'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [460,340];
    info.data.fragments_seq = [1:460,1:340];
end
if (strcmp(info.animal,'Cardano') && strcmp(info.date,'20210401'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [220,580];
    info.data.fragments_seq = [1:220,1:580];
end
if (strcmp(info.animal,'Frida') && strcmp(info.date,'20210410'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [400,400];
    info.data.fragments_seq = [1:400,1:400];
end
if (strcmp(info.animal,'Frida') && strcmp(info.date,'20210411'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [440,60];
    info.data.fragments_seq = [1:440,1:60];
end
if (strcmp(info.animal,'Frida') && strcmp(info.date,'20210412'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [60,740];
    info.data.fragments_seq = [1:60,1:740];
end
if (strcmp(info.animal,'Legolas') && strcmp(info.date,'20210423'))
    info.data.numFragments = 2;
    info.data.fragments_taskContinuous = true;
    info.data.fragments_numTrials = [250,250];
    info.data.fragments_seq = [1:250,1:250];
end
if (strcmp(info.animal,'Ripple') && strcmp(info.date,'20210525'))
    info.data.numFragments = 3;
    info.data.fragments_numTrials = [340,160,300];
    info.data.fragments_seq = [1:340,1:160,1:300];
    info.data.numSubFragments_paq = [2,1,1];
end
if (strcmp(info.animal,'Ripple') && strcmp(info.date,'20210526'))
    info.data.numFragments = 4;
    info.data.fragments_numTrials = [340,60,60,40];
    info.data.fragments_seq = [1:340,1:60,1:60,1:40];
end
if (strcmp(info.animal,'Sauron') && strcmp(info.date,'20210621'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [80,420];
    info.data.fragments_seq = [1:80,1:420];
end
if (strcmp(info.animal,'Merry') && strcmp(info.date,'20210702'))
    info.data.numFragments = 3;
    info.data.fragments_numTrials = [60,640,100];
    info.data.fragments_seq = [1:60,1:640,1:100];
end
if (strcmp(info.animal,'Celo') && strcmp(info.date,'20210705'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [600,200];
    info.data.fragments_seq = [1:600,1:200];
end
if (strcmp(info.animal,'Celo') && strcmp(info.date,'20210706'))
    info.data.numFragments = 3;
    info.data.fragments_numTrials = [80,380,340];
    info.data.fragments_seq = [1:80,1:380,1:340]; % last block has 20 trials more
end
if (strcmp(info.animal,'Shaw') && strcmp(info.date,'20210813'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [480,20];
    info.data.fragments_seq = [1:480,1:20];
end
if (strcmp(info.animal,'Sterling') && strcmp(info.date,'20210812'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [440,360];
    info.data.fragments_seq = [1:440,1:360];
end
if (strcmp(info.animal,'Maguire') && strcmp(info.date,'20210821'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [380,420];
    info.data.fragments_seq = [1:380,1:420];
end
if (strcmp(info.animal,'Pickford') && strcmp(info.date,'20210827'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [680,120];
    info.data.fragments_seq = [1:680,1:120];
end
if (strcmp(info.animal,'Pickford') && strcmp(info.date,'20210828'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [60,440];
    info.data.fragments_seq = [1:60,1:440];
end
if (strcmp(info.animal,'Stanage') && strcmp(info.date,'20210925'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [120,680];
    info.data.fragments_seq = [1:120,1:680];
end
if (strcmp(info.animal,'Philip') && strcmp(info.date,'20211006'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [760,40];
    info.data.fragments_seq = [1:760,1:40];
end
if (strcmp(info.animal,'Philip') && strcmp(info.date,'20211007'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [180,320];
    info.data.fragments_seq = [1:180,1:320];
end
if (strcmp(info.animal,'Diana') && strcmp(info.date,'20211024'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [600,200];
    info.data.fragments_seq = [1:600,1:200];
end
if (strcmp(info.animal,'Meghan') && strcmp(info.date,'20211124'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [120,380];
    info.data.fragments_seq = [1:120,1:380];
end
if (strcmp(info.animal,'Meghan') && strcmp(info.date,'20211125'))
    info.data.numFragments = 3;
    info.data.fragments_numTrials = [100,520,180];
    info.data.fragments_seq = [1:100,1:520,1:180];
end
if (strcmp(info.animal,'Meghan') && strcmp(info.date,'20211124'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [120,380];
    info.data.fragments_seq = [1:120,1:380];
end
if (strcmp(info.animal,'Java') && strcmp(info.date,'20211201'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [460,40];
    info.data.fragments_seq = [1:460,1:40];
end
if (strcmp(info.animal,'Python') && strcmp(info.date,'20211212'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [440,60];
    info.data.fragments_seq = [1:440,1:60];
end
if (strcmp(info.animal,'Python') && strcmp(info.date,'20211213'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [460,340];
    info.data.fragments_seq = [1:460,1:340];
end
if (strcmp(info.animal,'Hope') && strcmp(info.date,'20220130'))
    info.data.numFragments = 3;
    info.data.fragments_numTrials = [200,200,400];
    info.data.fragments_seq = [1:200,1:200,1:400];
end
if (strcmp(info.animal,'BullyBoy') && strcmp(info.date,'20220217'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [680,120];
    info.data.fragments_seq = [1:680,1:120];
end
if (strcmp(info.animal,'Ao') && strcmp(info.date,'20220403'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [440,360];
    info.data.fragments_seq = [1:440,1:360];
end
if (strcmp(info.animal,'Ao') && strcmp(info.date,'20220405'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [520,280];
    info.data.fragments_seq = [1:520,1:280];
end
if (strcmp(info.animal,'Jobs') && strcmp(info.date,'20220719'))
    info.data.numFragments = 2;
    info.data.fragments_numTrials = [740,60];
    info.data.fragments_seq = [1:740,1:60];
end

% Stim trigger issues and premature imaging stops
if (strcmp(info.animal,'Cardano') && strcmp(info.date,'20210329'))
    info.data.prematureImagingStop_nanAfter = 480;
end
if (strcmp(info.animal,'Legolas') && strcmp(info.date,'20210423'))
    info.data.prematureImagingStop_nanAfter = 497;
end
if (strcmp(info.animal,'Arwen') && strcmp(info.date,'20210610'))
    info.data.prematureImagingStop_nanAfter = 218; % note there is a second paq fragment which is 240 trials
end
if (strcmp(info.animal,'Shaw') && strcmp(info.date,'20210810'))
    info.data.stimTriggerIssue_nanAfter = 300; % last correct stim is stim 250 (trial 311), stim 312 (trial 769) is missed, all following stims are wrong
end
if (strcmp(info.animal,'Philip') && strcmp(info.date,'20211004'))
    info.data.stimTriggerIssue_nanAfter = 760; % last correct stim is stim 614 (trial 768), stim 615 (trial 769) is missed, all following stims are wrong
end
if (strcmp(info.animal,'Lizzy') && strcmp(info.date,'20211102'))
    info.data.prematureImagingStop_nanAfter = 418;
end

% Missing files
if (strcmp(info.animal,'Pfizer') && strcmp(info.date,'20210208'))
    info.data.missingPaqPreFile = true;
end
if (strcmp(info.animal,'Pfizer') && strcmp(info.date,'20210212'))
    info.data.missingPaqPostFile = true;
end
if (strcmp(info.animal,'Oxford') && strcmp(info.date,'20210208'))
    info.data.missingPaqPreFile = true;
end
if (strcmp(info.animal,'Oxford') && strcmp(info.date,'20210212'))
    info.data.missingPaqPostFile = true;
end
if (strcmp(info.animal,'Cardano') && strcmp(info.date,'20210401'))
    info.data.missingThorBeh = true; % not strictly missing, but somehow corrupted and can't load it into Matlab
end
if (strcmp(info.animal,'Faramir') && strcmp(info.date,'20210610'))
    info.data.missingThorBeh = true;
end
if (strcmp(info.animal,'Merry') && strcmp(info.date,'20210702'))
    info.data.missingMatFile = [1,0,0];
end
if (strcmp(info.animal,'Celo') && strcmp(info.date,'20210706'))
    info.data.missingMatFile = [1,0,0];
    info.data.missingPaqPreFile = true;
    info.data.missingThorBeh = true; % second ThorSync file is corrupted
end
if (strcmp(info.animal,'Sterling') && strcmp(info.date,'20210812'))
    info.data.missingMatFile = [1,0];
end
if (strcmp(info.animal,'Shaw') && strcmp(info.date,'20210813'))
    info.data.missingMatFile = [1,0];
end
if (strcmp(info.animal,'Kane') && strcmp(info.date,'20210911'))
    info.data.missingPaqPostFile = true;
end
if (strcmp(info.animal,'Philip') && strcmp(info.date,'20211006'))
    info.data.missingMatFile = [1,0];
end
if (strcmp(info.animal,'Hope') && strcmp(info.date,'20220130'))
    info.data.missingPaqPostFile = true;
end
if (strcmp(info.animal,'Jobs') && strcmp(info.date,'20220721'))
    info.data.missingPaqPreFile = true;
end

% Special case fixers
if (strcmp(info.animal,'Lizzy') && strcmp(info.date,'20211030'))
    info.data.specialCaseFixer = 1; % behaviour control was restarted 8 trials in while everything else was kept running
end

% Append fields for uncorrupted data
if ~isfield(info.data,'numSubFragments_paq')
    info.data.numSubFragments_paq = ones(1,info.data.numFragments);
end
if ~isfield(info.data,'missingThorBeh')
    info.data.missingThorBeh = false;
end
if ~isfield(info.data,'missingMatFile')
    info.data.missingMatFile = zeros(1,info.data.numFragments);
end
if ~isfield(info.data,'missingPaqPreFile')
    info.data.missingPaqPreFile = false;
end
if ~isfield(info.data,'missingPaqPostFile')
    info.data.missingPaqPostFile = false;
end
if ~isfield(info.data,'fragments_taskContinuous')
    info.data.fragments_taskContinuous = false;
end


%% Add data to repo

out = struct();
out.info = info;

if ops.do_taskData
    disp('- [task] module')
    if info.expDay<100
        [path,task,perf] = data2repo_task(info,path);
        out.path = path; out.task = task; out.perf = perf;
    else
        data2repo_task(info,path,'base1');
        data2repo_task(info,path,'base2');
        [path,task,perf] = data2repo_task(info,path,'stim');
        out.path = path; out.task = task; out.perf = perf;
    end
end

if ops.do_thorData || ops.do_s2pData || ops.do_paqData || ops.do_cascData || ops.do_zcorrData %|| ops.do_haloData || ops.do_spksnData || ops.do_s2pData_new
    if strcmp(info.task.room,'Valhalla')
        [path,raw] = data2repo_gatherRawDataInfo(info,path);
        if ops.do_thorData && (~info.data.missingThorBeh)
            disp('- [thor] module')
            [path,thor_beh] = data2repo_thor(info,p,path,raw);
            out.thor_beh = thor_beh;
        end
        if ops.do_s2pData
            disp('- [s2p] module')
            [path,s2p_meta] = data2repo_s2p(info,p,path,raw,ops);
            out.s2p_meta = s2p_meta;
        end
        if ops.do_cascData
            disp('- [casc] module')
            if ~exist('s2p_meta','var')
                disp('--- Loading s2p_meta...')
                load([path.filepart_out,'s2p_meta.mat']);
            end
            [path] = data2repo_casc(info,path,raw,ops,s2p_meta);
        end
    end
    if ops.do_paqData
        disp('- [paq] module')
        if ~exist('task','var')
            disp('--- Loading task...')
            load([path.filepart_out,'task.mat']);
        end
        if exist('thor_beh','var')
            [path,paq_beh] = data2repo_paq(info,p,path,raw,task,thor_beh);
        elseif exist([path.filepart_out,'thor_beh.mat'],'file')
            disp('--- Loading thor_beh...')
            load([path.filepart_out,'thor_beh.mat']);
            [path,paq_beh] = data2repo_paq(info,p,path,raw,task,thor_beh);
        elseif strcmp(info.task.room,'Grid')
            [path,paq_beh] = data2repo_paq(info,p,path,[],task);
        else
            warning('Running [paq] module without thor_beh as thor_beh data was not available')
            [path,paq_beh] = data2repo_paq(info,p,path,raw,task);
        end
        if info.data.missingThorBeh
            disp('--- Loading paq_beh...')
            load([path.filepart_out,'paq_beh.mat']);
            if ~exist('task','var')
                disp('--- Loading task...')
                load([path.filepart_out,'task.mat']);
            end
            [path] = data2repo_thor_mock(info,path,paq_beh,task);
        end
        out.paq_beh = paq_beh;
    end
    if strcmp(info.task.room,'Valhalla')
        if ops.do_zcorrData
            disp('- [zcorr] module')
            if ~exist('s2p_meta','var')
                disp('--- Loading s2p_meta...')
                load([path.filepart_out,'s2p_meta.mat']);
            end
            if ~exist('paq_beh','var')
                disp('--- Loading paq_beh...')
                load([path.filepart_out,'paq_beh.mat']);
            end
            [path] = data2repo_zcorr(info,path,raw,ops,s2p_meta,paq_beh);
        end
        if ops.do_haloData
            disp('- [halo] module')
            if ~exist('s2p_meta','var')
                disp('--- Loading s2p_meta...')
                load([path.filepart_out,'s2p_meta.mat']);
            end
            if ~exist('paq_beh','var')
                disp('--- Loading paq_beh...')
                load([path.filepart_out,'paq_beh.mat']);
            end
            [path] = data2repo_halo(info,p,path,raw,ops,s2p_meta,paq_beh);
        end
    end
    out.path = path;
end

if strcmp(info.task.room,'Valhalla')
    if ops.do_s2pData_new
        disp('- [s2p_new] module')
        if ~exist('s2p_meta','var')
            disp('--- Loading s2p_meta...')
            load([path.filepart_out,'s2p_meta.mat']);
        end
        if ~exist('paq_beh','var')
            disp('--- Loading paq_beh...')
            load([path.filepart_out,'paq_beh.mat']);
        end
        if ~exist('thor_beh','var') && info.stimSession
            disp('--- Loading thor_beh...')
            load([path.filepart_out,'thor_beh.mat']);
            try
                [path,raw] = data2repo_gatherRawDataInfo(info,path);
                [path] = data2repo_s2p_new(info,p,path,ops,s2p_meta,paq_beh,thor_beh,raw);
            catch
                [path] = data2repo_s2p_new(info,p,path,ops,s2p_meta,paq_beh,thor_beh);
            end
        else
            try
                [path,raw] = data2repo_gatherRawDataInfo(info,path);
                [path] = data2repo_s2p_new(info,p,path,ops,s2p_meta,paq_beh,[],raw);
            catch
                [path] = data2repo_s2p_new(info,p,path,ops,s2p_meta,paq_beh);
            end
        end
    end
    if ops.do_spksnData
        disp('- [spksn] module')
        if ~exist('s2p_meta','var')
            disp('--- Loading s2p_meta...')
            load([path.filepart_out,'s2p_meta.mat']);
        end
        try
            [path,raw] = data2repo_gatherRawDataInfo(info,path);
            [path] = data2repo_spksn(info,p,path,ops,s2p_meta,raw);
        catch
            [path] = data2repo_spksn(info,p,path,ops,s2p_meta);
        end
    end
    if ops.do_camData
        disp('- [cam] module')
        if ~exist('s2p_meta','var')
            disp('--- Loading s2p_meta...')
            load([path.filepart_out,'s2p_meta.mat']);
        end
        try
            [path,raw] = data2repo_gatherRawDataInfo(info,path);
            [path] = data2repo_cam(info,p,path,ops,s2p_meta,raw);
        catch
            [path] = data2repo_cam(info,p,path,ops,s2p_meta);
        end
    end
    out.path = path;
end

if ops.s2p.update_curation
    disp('- [s2p] curation update')
    [path,s2p_meta] = data2repo_updateCuration(info,path);
    out.path = path; out.s2p_meta = s2p_meta;
end

if ops.do_trgData && info.stimSession
    disp('- [trg] module')
    if ~exist('task','var')
        disp('--- Loading task...')
        load([path.filepart_out,'task.mat']);
    end
    [path] = data2repo_trg(info,path,task);
    out.path = path;
end


%% Complete execution

log.done = true;
log.runTime = toc(t_start);

if ~exist([path.filepart_out,'meta.mat'],'file')
    disp('- [meta] module')
    [path] = data2repo_meta(info,log,ops,p,path);
    out.path = path;
elseif ops.s2p.update_meta
    disp('- [meta] module')
    [path] = data2repo_meta(info,log,ops,p,path);
    out.path = path;
end

disp(['- Done in ',num2str(log.runTime/60,3),' min.'])
diary off;

end
function [out] = analyses(animal,date,path,ops)
% info.animal = animal; info.date = date;

%% Preparations

% gather p
p = get_p;

% identify session folders
path.folder_repo = [path.root_repo,animal,'/',animal,'_',date,'/'];
path.filepart_in = [path.folder_repo,animal,'_',date,'_'];
path.folder_repoX = [path.root_repoX,animal,'/',animal,'_',date,'/'];
path.folder_repoXX = [path.root_repoXX,'dFF/'];
path.filepart_in_repoX = [path.folder_repoX,animal,'_',date,'_'];
path.filepart_in_repoXX = [path.folder_repoXX,animal,'_',date,'_'];
path.folder_analysis = [path.root_analysis,animal,'/',animal,'_',date,'/'];
path.filepart_in_analysis = [path.folder_analysis,animal,'_',date,'_'];
path.folder_analysisX = [path.root_analysisX,animal,'/',animal,'_',date,'/'];
path.filepart_in_analysisX = [path.folder_analysisX,animal,'_',date,'_'];
path.filepart_out = [path.folder_analysis,animal,'_',date,'_'];
path.filepart_outX = [path.folder_analysisX,animal,'_',date,'_'];
if ~exist(path.folder_analysis,'dir')
   mkdir(path.folder_analysis);
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
temp = [path.filepart_out,'log/',animal,'_',date,'_',datestr(log.executionDate,'yyyymmdd_HHMMSS'),'.log'];
diary(temp);
diary on;

disp([log.executionDate])
disp(['Running analyses.m for ',animal,'_',date,'.'])


%% Load data

disp('- Loading data.')

% load meta
load([path.filepart_in,'meta.mat']);
info = meta.info;
if info.stimSession
    try
        info.stimType = meta.stimType;
    catch
        info.stimType = 'missing';
        warning('Info about stim type not in meta file')
    end
else
    info.stimType = 'nostim';
end

% identify relevant data
cmpr_list = [];
% if ops.do_learningCurveAnalysis
%     cmpr_list = [cmpr_list,"perf"];
% end
if ops.do_responseAnalysis
	cmpr_list = [cmpr_list,meta.reg.trg,"task","s2p_meta","thor_beh",p.resp.activityMeasure];
end
if ops.do_imprintingAnalysis
	cmpr_list = [cmpr_list,meta.reg.trg,"task","s2p_meta","thor_beh",p.impr.activityMeasure];
end
if ops.do_imagingQualityAnalysis
    cmpr_list = [cmpr_list,"s2p_meta","paq_beh",p.iqa.activityMeasure];
end
if ops.do_tuningAnalysis
    cmpr_list = [cmpr_list,"task","s2p_meta","paq_beh",p.tng.activityMeasure];
end
if ops.do_inhibitionAnalysis
    cmpr_list = [cmpr_list,"task","s2p_meta","paq_beh",p.tng.activityMeasure,"tng_all_cmpr","nem_all_cmpr"];
end
if ops.do_errorTrialTuningAnalysis
    cmpr_list = [cmpr_list,"task","s2p_meta","paq_beh",p.ett.activityMeasure];
end
if ops.do_sequencenessAnalysis
    cmpr_list = [cmpr_list,"task","s2p_meta","paq_beh",p.sqn.activityMeasure];
    if ops.sqn.do_allTrials
        cmpr_list = [cmpr_list,"tng_all","nem_all_cmpr"];
    end
    if ops.sqn.do_100t
        cmpr_list = [cmpr_list,"tng_100t"];
    end
end
if ops.do_decodingAnalysis
    cmpr_list = [cmpr_list,"task","paq_beh",p.sqn.activityMeasure];
end
if ops.do_nemAnalysis
    cmpr_list = [cmpr_list,"task","s2p_meta","paq_beh",p.nem.activityMeasure];
end
if ops.do_nemAnalysis_cmpr
    cmpr_list = [cmpr_list,"nem_all","nem_100t"];
    %cmpr_list = [cmpr_list,"nem_all"];
end
if ops.do_osipAnalysis
    cmpr_list = [cmpr_list,"task","perf","s2p_meta","paq_beh",p.nem.activityMeasure,"tng_all","nem_all"];
end
if ops.do_ecaAnalysis
    cmpr_list = [cmpr_list,"task","s2p_meta","paq_beh",p.eca.activityMeasure];
    if ops.eca.do_allTrials
        cmpr_list = [cmpr_list,"nem_all_cmpr"];
    end
end
if ops.do_lickingAnalysis
    cmpr_list = [cmpr_list,"task","perf","paq_beh"];
end
cmpr_list = unique(cmpr_list);

% load data
cmpr = repo2ws(path,cmpr_list);
temp = fieldnames(cmpr);
for i=1:length(temp)
    if strcmp(temp{i},'F_beh') || strcmp(temp{i},'Fneu_beh') || strcmp(temp{i},'dFF_beh') || strcmp(temp{i},'spks_beh') || ...
            strcmp(temp{i},'dFFn_beh') || strcmp(temp{i},'spksn_beh')
        act_struct.(temp{i}) = cmpr.(temp{i});
    else
        eval([temp{i},'=cmpr.',temp{i},';']);
    end
end
if exist('trg_rigid','var')
    trg = trg_rigid;
    clear('trg_rigid');
elseif exist('trg_nonrigid','var')
    trg = trg_nonrigid;
    clear('trg_nonrigid');
end
clear('cmpr');

try
    try 
        iscell = meta.iscell;
    catch
        iscell = s2p_meta.iscell(:,1);
        warning('Info about iscell not in meta file')
    end
catch
    warning('Did not find iscell information')
end

% quick sanity checks
if exist('paq_beh','var') && exist('thor_beh','var')
    if ~( (length(paq_beh.sync)==info.task.numTrials) && (length(thor_beh.sync)==info.task.numTrials) && (sum(paq_beh.sync-thor_beh.sync)==0) )
        warning('There is something wrong with the SniffinSyncs of this session.')
    end
end
if exist('thor_beh','var')
    if info.stimSession
        if ~(length(thor_beh.stimSequence)==info.task.numStimTrials)
            warning('There is something wrong with the StimTriggers of this session.')
        end
    end
end
if exist('task','var') && exist('thor_beh','var')
    if any(diff(thor_beh.sync(task.var==1)-thor_beh.stimSequence)>=2)
        warning('There is something wrong with the timing of the StimTriggers of this session.')
    end    
end

% NaN frames with massive lateral offsets (presumably incorrect registration)
if exist('s2p_meta','var')
    absoff = hypot(s2p_meta.ops.xoff,s2p_meta.ops.yoff); % use this to identify frames with registration faults
    absdisp = diff(absoff);  % use this to identify motion
end


%% Do analyses

out = struct();
out.info = info;
out.meta = meta;
out.path = path;
out.p = p;
if exist('paq_beh','var')
    out.paq_beh = paq_beh;
end

if ops.do_responseAnalysis
    disp('- [resp] module')
        act = getActivityMeasure(p.resp,act_struct);
        [resp] = responseAnalysis(info,iscell,ops,p,path,trg,task,s2p_meta,thor_beh,act);
        out.resp = resp;
end

if ops.do_imagingQualityAnalysis
    disp('- [iqa] module')
    act = getActivityMeasure(p.iqa,act_struct);
    [iqa] = imagingQualityAnalysis(info,iscell,ops,p,path,s2p_meta,thor_beh,act);
    out.iqa = iqa;
end

if ops.do_tuningAnalysis
    disp('- [tng] module')
    act = getActivityMeasure(p.tng,act_struct);
    [tng,nft_binned] = tuningAnalysis(info,iscell,ops,p,path,task,s2p_meta,paq_beh,act);
    out.tng = tng; out.nft_binned = nft_binned;
end

if ops.do_inhibitionAnalysis
    disp('- [inh] module')
    act = getActivityMeasure(p.inh,act_struct);
    inh = inhibitionAnalysis(info,iscell,ops,p,path,task,s2p_meta,paq_beh,act);
    out.inh = inh;
end

if ops.do_errorTrialTuningAnalysis
    disp('- [ett] module')
    act = getActivityMeasure(p.ett,act_struct);
    [ett] = errorTrialTuningAnalysis(info,iscell,ops,p,path,task,s2p_meta,paq_beh,act);
    out.ett = ett;
end

if ops.do_sequencenessAnalysis
    disp('- [sqn] module')
    act = getActivityMeasure(p.sqn,act_struct);
    if ops.sqn.do_allTrials && ops.sqn.do_100t
        [sqn] = sequencenessAnalysis(info,iscell,ops,p,path,task,s2p_meta,paq_beh,act,tng_all,tng_100t,nem_all_cmpr);
    elseif ops.sqn.do_allTrials
        [sqn] = sequencenessAnalysis(info,iscell,ops,p,path,task,s2p_meta,paq_beh,act,tng_all,[],nem_all_cmpr);
    elseif ops.sqn.do_100t
        [sqn] = sequencenessAnalysis(info,iscell,ops,p,path,task,s2p_meta,paq_beh,act,[],tng_100t);
    end
    out.sqn = sqn;
end

if ops.do_ecaAnalysis
    disp('- [eca] module')
    act = getActivityMeasure(p.eca,act_struct);
    if ops.eca.do_allTrials
        encodingCellActivityAnalysis(info,iscell,ops,p,path,task,paq_beh,act,nem_all_cmpr,[]);
    end
end

if ops.do_decodingAnalysis
    disp('- [dec] module')
    act = getActivityMeasure(p.dec,act_struct);
    decodingAnalysis(info,iscell,ops,p,path,task,paq_beh,act);
end

if ops.do_nemAnalysis
    disp('- [nem] module')
    act = getActivityMeasure(p.nem,act_struct);
    nemAnalysis(info,iscell,ops,p,path,task,paq_beh,act);
end
if ops.do_nemAnalysis_cmpr
    disp('- [nem] module (cmpr)')
    nemAnalysis_cmpr(path,nem_all,nem_100t);
end

if ops.do_osipAnalysis
    disp('- [osip] module')
    act = getActivityMeasure(p.osip,act_struct);
    osipAnalysis(info,iscell,ops,p,path,task,paq_beh,act,tng_nem,nem_all);
end

if ops.do_imprintingAnalysis
    disp('- [impr] module')
    act = getActivityMeasure(p.impr,act_struct);
    [impr] = imprintingAnalysis(info,iscell,ops,p,path,task,s2p_meta,paq_beh,act);
    out.impr = impr;
end

if ops.do_lickingAnalysis
    disp('- [lick] module')
    [lick] = lickingAnalysis(info,ops,p,path,task,perf,paq_beh);
    out.lick = lick;
end


%% Complete execution

log.done = true;
log.runTime = toc(t_start);

disp(['- Done in ',num2str(log.runTime/60,3),' min.'])
diary off;

end


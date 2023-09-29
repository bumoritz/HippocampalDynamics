%% Select data
clear; clc;

animal = 'Guzman'; date = '20200805'; expDay = 1; stimSession = 0;


%% Running options

% core
ops.do_taskData                     = true;
ops.do_thorData                     = false;
ops.do_s2pData                      = false;
ops.do_paqData                      = false;
ops.do_trgData                      = false;

% new additional core 
ops.do_s2pData_new                  = false;
ops.do_spksnData                    = false;

%abandoned
ops.do_s2pData_newn                 = false;
ops.do_spksnnData                   = false;

% add-ons
ops.do_cascData                     = false;
ops.do_zcorrData                    = false;
ops.do_haloData                     = false;
ops.do_camData                      = false;

ops.s2p.update_curation             = false;
ops.s2p.update_meta                 = false; % keep at false!

% options for base (should be all false)
ops.s2p.saveFandFneu                = false;
ops.s2p.do_dFF_gm                   = false;
ops.s2p.use_mat_data                = false;
ops.s2p.import_Fns                  = false;
ops.s2p.skip_convert_ops_and_stat   = false;
ops.s2p.skip_neuropil_subtraction   = false;
ops.s2p.skip_baseline_fluorescence  = false;
ops.s2p.skip_dFF                    = false;
ops.s2p.skip_spks                   = false;

% options for add-ons
ops.casc.do_50g                     = false;
ops.casc.do_50c                     = false;
ops.casc.do_100g                    = false;
ops.casc.do_100c                    = false;
ops.casc.do_200g                    = false;
ops.cam.do_mot                      = false;
ops.cam.do_dlc                      = false;

path.root_data                      = 'Z:\WIBR_Hippos\SniffinHippo\Data\'; %  %'F:\Data\'; % 'E:\SniffinHippo\CamData\'
path.root_repo                      = 'D:\SniffinHippo\RepoO\'; %'D:\SniffinHippo\Repo\'; % 'D:\temp\';
path.root_repoX                     = 'E:\SniffinHippo\RepoX\'; % 'D:\temp\';


%% Run data2repo

out = data2repo(animal,date,expDay,stimSession,path,ops);


%% Extract

temp = fieldnames(out);
for i=1:length(temp)
    eval([temp{i},'=out.',temp{i},';']);
end
clear('out');


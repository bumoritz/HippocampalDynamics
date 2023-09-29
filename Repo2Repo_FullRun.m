%% Select data
clear; clc;

% animal = 'Turing'; date = ['20220913';'20220914';'20220915';'20220916';'20220917'];
% animal = 'Zuse'; date = ['20220922';'20220923';'20220924';'20220926';'20220927'];
% animal = 'Austin'; date = ['20220726';'20220727';'20220728';'20220729';'20220730'];
% animal = 'Berners'; date = ['20221010';'20221011';'20221012'];
% animal = 'Dickinson'; date = ['20220824';'20220825';'20220826';'20220827';'20220828'];
% animal = 'Musk'; date = ['20220617';'20220618';'20220619';'20220620';'20220622'];
% animal = 'Biontech'; date = ['20210208';'20210209';'20210210';'20210211';'20210212'];
% animal = 'Arasaka'; date = ['20210323';'20210324';'20210325';'20210326';'20210327'];
% animal = 'Elrond'; date = ['20210531';'20210601';'20210602'];
% animal = 'Merry'; date = ['20210629';'20210630';'20210701';'20210702';'20210703'];
% animal = 'Pippin'; date = ['20210523';'20210524';'20210525';'20210526';'20210527'];
% animal = 'Legolas'; date = ['20210419';'20210420';'20210421';'20210422';'20210423'];
% animal = 'Sauron'; date = ['20210617';'20210618';'20210619';'20210620';'20210621'];
% animal = 'Abdi'; date = ['20210705';'20210706';'20210707';'20210708';'20210709'];
% animal = 'Maguire'; date = ['20210820';'20210821';'20210822';'20210823';'20210824'];
% animal = 'Maria'; date = ['20211021';'20211022';'20211023';'20211024';'20211025'];
% animal = 'Lizzy'; date = ['20211029';'20211030';'20211031';'20211101';'20211102'];
% animal = 'William'; date = ['20211108';'20211109';'20211110';'20211111';'20211112'];
% animal = 'Ao'; date = ['20220402';'20220403';'20220404';'20220405';'20220406'];
% animal = 'Stanage'; date = ['20210924';'20210925';'20210926';'20210927';'20210928'];
% animal = 'Shaw'; date = ['20210809';'20210810';'20210811';'20210812';'20210813'];
% animal = 'Pickford'; date = ['20210825';'20210827';'20210828';'20210829';'20210830'];
% animal = 'Celo'; date = ['20210704';'20210706';'20210707';'20210708';'20210709'];
% animal = 'Margaret'; date = ['20211011';'20211012';'20211013';'20211014'];
% animal = 'Philip'; date = ['20211003';'20211004';'20211005';'20211006';'20211007'];
% animal = 'Faramir'; date = ['20210607';'20210608';'20210609';'20210610';'20210611'];
% animal = 'BullyBoy'; date = ['20220214';'20220215';'20220216';'20220217';'20220218'];
% animal = 'Hope'; date = ['20220129';'20220130';'20220131';'20220201';'20220202'];
% animal = 'Kura'; date = ['20220524';'20220525';'20220526';'20220527';'20220528'];
% animal = 'Jobs'; date = ['20220718';'20220719';'20220720';'20220721';'20220722'];
%animal = 'Kane'; date = ['20210911';'20210912'];


do_cbalData     = true;
do_trgData      = true;
do_iscellData   = true;
do_trckData     = true;


%%% --- %%% --- %%% --- %%% --- %%% --- %%%

%% Running options

ops.winnerReg2meta                  = true;

ops.trck.do_all                     = true;
ops.trck.do_123                     = true;
ops.trck.do_345                     = true;
ops.trck.do_rigid                   = true;
ops.trck.do_nonrigid                = true;
ops.trck.showFigures                = false;

ops.trg.do_rigid                    = true;
ops.trg.do_nonrigid                 = true;
ops.trg.showFigures                 = false;

ops.cbal.showFigures                = true;

path.root_repo                      = 'D:\SniffinHippo\Repo\';
path.root_repoX                     = 'E:\SniffinHippo\RepoX\';


%% Get animal info

disp(['Repo2Repo Full Run for ',animal,'.',newline])

% get animal info from first potential stim session
numSessions = size(date,1);
if numSessions>=2
    temp_date = date(2,:);
    path.folder_repo = [path.root_repo,animal,'\',animal,'_',temp_date,'\'];
    path.filepart_out = [path.folder_repo,animal,'_',temp_date,'_'];
    if exist([path.filepart_out,'meta.mat'])==2
        load([path.filepart_out,'meta.mat']);
    end
    stimSession = meta.info.stimSession;
else
    stimSession = false;
end

if stimSession && numSessions<5
    error('Will not run properly because less than 5 sessions.')
end


%% [cbal] block

if stimSession
    ops.do_cbalData                     = do_cbalData;
    ops.do_trgData                      = false;
    ops.do_iscellData                   = false;
    ops.do_trckData                     = false;
    if ops.do_cbalData
        repo2repo(animal,[date(2,:);date(4,:)],path,ops);
    end
end


%% [trg] block

if stimSession
    ops.do_cbalData                     = false;
    ops.do_trgData                      = do_trgData;
    ops.do_iscellData                   = false;
    ops.do_trckData                     = false;
    if ops.do_trgData
        for i=2:5
            repo2repo(animal,date(i,:),path,ops);
        end
    end
end


%% NEW INFO block

ops.do_cbalData                         = false;
ops.do_trgData                          = false;
ops.do_iscellData                       = do_iscellData;
ops.do_trckData                         = false;
if ops.do_iscellData
    disp('- [NEW INFO] block')
    for i=1:numSessions
        
        % load meta
        temp_date = date(i,:);
        path.folder_repo = [path.root_repo,animal,'\',animal,'_',temp_date,'\'];
        path.filepart_out = [path.folder_repo,animal,'_',temp_date,'_'];
        if exist([path.filepart_out,'meta.mat'])==2
            load([path.filepart_out,'meta.mat']);
        end
        
        % load and extract s2p_meta to iscell
        if exist([path.filepart_out,'s2p_meta.mat'])==2
            load([path.filepart_out,'s2p_meta.mat']);
        end        
        meta.iscell_raw = s2p_meta.iscell(:,1);
        
        % imaging only
        if ~meta.info.stimSession
            meta.iscell = meta.iscell_raw;
        
        % stim
        else
            % stimType
            if exist([path.filepart_out,char(meta.reg.trg),'.mat'])==2
                load([path.filepart_out,char(meta.reg.trg),'.mat']);
            else
                error('A meta file did not have a reg field.')
            end 
            if exist('trg_rigid','var')
                trg = trg_rigid;
                clear('trg_rigid');
            elseif exist('trg_nonrigid','var')
                trg = trg_nonrigid;
                clear('trg_nonrigid');
            end
            meta.stimType = trg.p.stimType;
            
            % iscell
            meta.iscell = meta.iscell_raw;
            meta.iscell(rmmissing(trg.idcs_targetedCells(:))) = 1;
        end
        
        % save meta
        save([path.filepart_out,'meta.mat'],'meta','-v7.3');
        disp(['--- Overwrote meta file in repo at ',[path.filepart_out,'meta.mat'],'.'])    
    end
end


%% [trck] block

if numSessions>=3
    ops.do_cbalData                     = false;
    ops.do_trgData                      = false;
    ops.do_iscellData                   = false;
    ops.do_trckData                     = do_trckData;
    if ops.do_trckData
        repo2repo(animal,date,path,ops);
    end
end







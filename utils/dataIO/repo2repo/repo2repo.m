function [out] = repo2repo(animal,date,path,ops)

%% Preparations

% gather info and p
p = get_p;
info.animal = animal;

if size(date,1)>1
    
    % mutli-session processing
    info.date = date;
    path.folder_repo = [path.root_repo,info.animal,'\',info.animal,'_combined\'];
    path.folder_repoX = [path.root_repoX,info.animal,'\',info.animal,'_combined\'];
    path.filepart_out = [path.folder_repo,info.animal,'_combined_'];
    path.filepart_outX = [path.folder_repoX,info.animal,'_combined_'];
    if ~exist(path.folder_repo,'dir')
       mkdir(path.folder_repo);
    end
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
    temp = [path.filepart_out,'log\',info.animal,'_combined_',datestr(log.executionDate,'yyyymmdd_HHMMSS'),'.log'];
    diary(temp);
    diary on;

    disp([log.executionDate])
    disp(['Running repo2repo.m for ',info.animal,'_combined.'])
else
    
    % single session processing
    info.date = date;
    path.folder_repo = [path.root_repo,info.animal,'\',info.animal,'_',info.date,'\'];
    path.folder_repoX = [path.root_repoX,info.animal,'\',info.animal,'_',info.date,'\'];
    path.filepart_out = [path.folder_repo,info.animal,'_',info.date,'_'];
    path.filepart_outX = [path.folder_repoX,info.animal,'_',info.date,'_'];
    if ~exist(path.folder_repo,'dir')
       mkdir(path.folder_repo);
    end
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
    disp(['Running repo2repo.m for ',info.animal,'_',info.date,'.'])
end


%% Repo to repo

out = struct();
out.info = info;

if ops.do_trckData
    disp('- [trck] module')
    if ops.trck.do_all
        if ops.trck.do_rigid
            [path,trck_15_rigid] = repo2repo_trck(info,ops,p,path,[1:size(date,1)],'rigid');
        end
        if ops.trck.do_nonrigid
            [path,trck_15_nonrigid] = repo2repo_trck(info,ops,p,path,[1:size(date,1)],'nonrigid');
        end
    end
    if ops.trck.do_123 && size(date,1)==5
       if ops.trck.do_rigid
           [path,trck_13_rigid] = repo2repo_trck(info,ops,p,path,[1:3],'rigid');
       end
       if ops.trck.do_nonrigid
           [path,trck_13_nonrigid] = repo2repo_trck(info,ops,p,path,[1:3],'nonrigid');
       end
    end
    if ops.trck.do_345 && size(date,1)==5
        if ops.trck.do_rigid
            [path,trck_35_rigid] = repo2repo_trck(info,ops,p,path,[3:5],'rigid');
        end
        if ops.trck.do_nonrigid
            [path,trck_35_nonrigid] = repo2repo_trck(info,ops,p,path,[3:5],'nonrigid');
        end
    end
    out.path = path;
end

if ops.do_trgData
    disp('- [trg] module')
    if ops.trg.do_rigid
        [path,trg_rigid] = repo2repo_trg(info,ops,p,path,'rigid');
    end
    if ops.trg.do_nonrigid
        [path,trg_nonrigid] = repo2repo_trg(info,ops,p,path,'nonrigid');
    end
    out.path = path;
end

if ops.winnerReg2meta
    disp('- Updating meta with registration info.')
    
    % load required data
    if size(date,1)>1
        if exist([path.filepart_out,'meta.mat'])==2
            load([path.filepart_out,'meta.mat']);
        else
            meta = {};
        end
        try
            if (~exist('trck_15_rigid','var')) | (~exist('trck_15_nonrigid','var'))
                trck_15_rigid = load([path.filepart_out,'trck_15_rigid.mat']);
                trck_15_rigid = trck_15_rigid.trck;
                trck_15_nonrigid = load([path.filepart_out,'trck_15_nonrigid.mat']);
                trck_15_nonrigid = trck_15_nonrigid.trck;
            end
            if sum(sum(trck_15_rigid.cell_to_index_map>0,2)==size(trck_15_rigid.cell_to_index_map,2)) > sum(sum(trck_15_nonrigid.cell_to_index_map>0,2)==size(trck_15_nonrigid.cell_to_index_map,2))
                meta.reg.trck_15 = "trck_15_rigid";
            else
                meta.reg.trck_15 = "trck_15_nonrigid";
            end
        catch
            disp('--- Skipped trck_15.')
        end
        try
            if (~exist('trck_13_rigid','var')) | (~exist('trck_13_nonrigid','var'))
                trck_13_rigid = load([path.filepart_out,'trck_13_rigid.mat']);
                trck_13_rigid = trck_13_rigid.trck;
                trck_13_nonrigid = load([path.filepart_out,'trck_13_nonrigid.mat']);
                trck_13_nonrigid = trck_13_nonrigid.trck;
            end
            if sum(sum(trck_13_rigid.cell_to_index_map(:,2:3)>0,2)==2) > sum(sum(trck_13_nonrigid.cell_to_index_map(:,2:3)>0,2)==2)
                meta.reg.trck_13 = "trck_13_rigid";
            else
                meta.reg.trck_13 = "trck_13_nonrigid";
            end
        catch
            disp('--- Skipped trck_13.')
        end
        try
            if (~exist('trck_35_rigid','var')) | (~exist('trck_35_nonrigid','var'))
                trck_35_rigid = load([path.filepart_out,'trck_35_rigid.mat']);
                trck_35_rigid = trck_35_rigid.trck;
                trck_35_nonrigid = load([path.filepart_out,'trck_35_nonrigid.mat']);
                trck_35_nonrigid = trck_35_nonrigid.trck;
            end
            if sum(sum(trck_35_rigid.cell_to_index_map(:,2:3)>0,2)==2) > sum(sum(trck_35_nonrigid.cell_to_index_map(:,2:3)>0,2)==2)
                meta.reg.trck_35 = "trck_35_rigid";
            else
                meta.reg.trck_35 = "trck_35_nonrigid";
            end
        catch
            disp('--- Skipped trck_35.')
        end
    else
        load([path.filepart_out,'meta.mat']);
        if (~exist('trg_rigid','var')) | (~exist('trg_nonrigid','var'))
            trg_rigid = load([path.filepart_out,'trg_rigid.mat']);
            trg_rigid = trg_rigid.trg;
            trg_nonrigid = load([path.filepart_out,'trg_nonrigid.mat']);
            trg_nonrigid = trg_nonrigid.trg;
        end
        if trg_rigid.numIdentifiedTargeted > trg_nonrigid.numIdentifiedTargeted
            meta.reg.trg = "trg_rigid";
        else
            meta.reg.trg = "trg_nonrigid";
        end
    end
    
    % update meta
    save([path.filepart_out,'meta.mat'],'meta','-v7.3');
    disp(['--- Overwrote meta file in repo at ',[path.filepart_out,'meta.mat'],'.'])
    out.meta = meta;
end

if ops.do_cbalData
    disp('- [cbal] module')
    path = repo2repo_cbal(info,ops,p,path);
    out.path = path;
end


%% Complete execution

log.done = true;
log.runTime = toc(t_start);
disp(['- Done in ',num2str(log.runTime/60,3),' min.'])
diary off;

end
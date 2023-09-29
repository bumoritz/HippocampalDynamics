%% Select data
clear; clc;

%animal = 'Biontech'; date = ['20210208';'20210209';'20210210';'20210211';'20210212'];
animal = 'Philip'; date = '20211006';


%% Running options

ops.do_cbalData                     = false;
ops.do_trgData                      = true;
% ops.do_iscellData                   = false;
ops.do_trckData                     = false;

ops.winnerReg2meta                  = false;

ops.trck.do_all                     = false;
ops.trck.do_123                     = false;
ops.trck.do_345                     = false;
ops.trck.do_rigid                   = false;
ops.trck.do_nonrigid                = false;
ops.trck.showFigures                = false;

ops.trg.do_rigid                    = true;
ops.trg.do_nonrigid                 = true;
ops.trg.showFigures                 = false;

ops.cbal.showFigures                = true;

path.root_repo                      = 'D:\SniffinHippo\Repo\';
path.root_repoX                     = 'E:\SniffinHippo\RepoX\';


%% Run data2repo

out = repo2repo(animal,date,path,ops);


%% Extract

temp = fieldnames(out);
for i=1:length(temp)
    eval([temp{i},'=out.',temp{i},';']);
end
clear('out');


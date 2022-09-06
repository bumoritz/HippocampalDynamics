function [path] = data2repo_trg(info,path,task)

%% Load trg data

disp('--- Loading trg data...')

temp = dir([path.folder_data,'Targeting\TargetSelection\Targeting_',info.animal,'_',info.date,'_*.mat']);
path.file_in_trg = [path.folder_data,'Targeting\TargetSelection\',temp.name];
trg = load(path.file_in_trg);
trg = trg.trg;

% troubleshoot session fragmentations
if info.data.numFragments>1
    trg.sequenceOrder_raw = trg.sequenceOrder;
    temp = [];
    for i=1:info.data.numFragments
        temp2 = sum(task.var(1:info.data.fragments_numTrials(i)));
        temp = [temp,trg.sequenceOrder_raw(1:temp2)];
    end
    trg.sequenceOrder = temp;
end

% save reg_spont file
trg_online = orderfields(trg);
save([path.filepart_outX,'trg_online.mat'],'trg_online','-v7.3');
disp(['--- Saved trg_online file to repoX as ',[path.filepart_outX,'trg_online.mat'],'.'])

end




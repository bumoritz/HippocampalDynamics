function [path] = data2repo_meta(info,log,ops,p,path)

disp('--- Saving meta data...')

path.file_out_meta = [path.folder_repo,info.animal,'_',info.date,'_meta.mat'];

meta.info = orderfields(info);
meta.log = orderfields(log);
meta.ops = orderfields(ops);
meta.p = orderfields(p);
meta.path = orderfields(path);
meta = orderfields(meta);

save(path.file_out_meta,'meta','-v7.3');
disp(['--- Added meta file to repo as ',path.file_out_meta,'.'])

end
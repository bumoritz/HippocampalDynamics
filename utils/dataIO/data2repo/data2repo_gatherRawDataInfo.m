function [path,raw] = data2repo_gatherRawDataInfo(info,path)

%% Gather raw data info

% imaging raw data info -> movie size in bytes
for e = 1:length(info.epochs)
    if e~=2 | info.data.numFragments==1
        temp = [path.folder_data,'Imaging\',info.animal,'-',info.date,'-',info.epochs{e},'\'];
        temp2 = dir([temp,'Image_00*001.raw']);
        path.(['file_in_raw_',info.epochs{e}]) = [temp,temp2.name];
        raw.(info.epochs{e}).bytes = temp2.bytes;
    else
        for i=1:info.data.numFragments
            temp = [path.folder_data,'Imaging\',info.animal,'-',info.date,'-',info.epochs{e},'-',num2str(i),'\'];
            temp2 = dir([temp,'Image_00*001.raw']);
            path.(['file_in_raw_',info.epochs{e},'_',num2str(i)]) = [temp,temp2.name];
            raw.([info.epochs{e},'_',num2str(i)]).bytes = temp2.bytes;
        end
    end
end

% movie size in bytes -> movie size in frames
temp = fields(raw);
for j=1:length(temp)
    raw.(temp{j}).frames = raw.(temp{j}).bytes /(info.scope.bytesPerNumber*info.scope.fovSize_pix^2);
    if raw.(temp{j}).frames==0
        error(['There is a raw video file without any frames']);
    end
end

% MB20210720: troubleshoot for example Cardano_20210401
temp = fields(raw);
i=0;
for j=1:length(temp)
    if raw.(temp{j}).frames==600000 | raw.(temp{j}).frames==350000
        i=i+1;
        if i==1
            raw.(temp{j}).frames_original = raw.(temp{j}).frames;
            temp1 = [path.folder_data,'Imaging\'];
            temp2 = dir([temp1,'*-merged.raw']);
            total_frames = temp2.bytes /(info.scope.bytesPerNumber*info.scope.fovSize_pix^2);
            temp3 = total_frames;
            for k=1:length(temp)
                if k~=j
                    temp3 = temp3 - raw.(temp{k}).frames;
                end
            end
            raw.(temp{j}).frames = temp3;
            i=0;
        end
    end
end

end
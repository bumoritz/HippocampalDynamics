function data = ThorLink_ReadThorSync_ov20220105(filePath, varargin)
% read thorlabs thorsync h5 file into a matlab struct - fieldnames are dataset names
% varargin = channels to load
% Lloyd Russell 2017-04-12


if nargin > 1
    channels = varargin{1};
else
    channels = [];
end

fileInfo = h5info(filePath);
groupInfo = fileInfo.Groups;
numGroups = numel(groupInfo);

data = [];
if numGroups
    for i = 1:numGroups
        groupName = groupInfo(i).Name;
        datasetInfo = groupInfo(i).Datasets;
        numDatasets = numel(datasetInfo);
        if numDatasets
            for j = 1: numDatasets
                datasetName = groupInfo(i).Datasets(j).Name;
                validDatasetName = genvarname(datasetName);
                
                if ~isempty(channels)
                    if any(contains(channels, validDatasetName))
                        completeDatasetName = [groupName '/' datasetName];
                        data.(validDatasetName) = h5read(filePath, completeDatasetName);
                    end
                else
                    completeDatasetName = [groupName '/' datasetName];
                    data.(validDatasetName) = h5read(filePath, completeDatasetName);
                end
            end
        end
    end
end

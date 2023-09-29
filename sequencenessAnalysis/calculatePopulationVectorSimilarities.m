function data_out = calculatePopulationVectorSimilarities(data_in_1,data_in_2,metric) % data_in_1 is template, data_in_2 is test
%data_in_1 = povSim.dynTemplates.iscells_full_ABtrials; data_in_2 = povSim.dynTemplates.iscells_full_AYtrials; metric = 'cosTheta'; % (or 'Pearson')
%data_in_1 = povSim.dynTemplates.iscells_full_ABtrials; data_in_2 = povSim.dynTemplates.iscells_full_ABtrials; metric = 'Pearson';

numCells = size(data_in_2,1);
numBins = size(data_in_2,2);
numTrials = size(data_in_2,3);

if size(data_in_1,1)~=size(data_in_2,1)
    error('numCells of input data not the same.');
end
if size(data_in_1,2)~=size(data_in_2,2)
    error('numBins of input data not the same.');
end

% core
if strcmp(metric,'cosTheta')
    if numTrials==1
        data_out = nan(numBins,numBins);
        for i=1:numBins
            for j = 1:numBins
                u = data_in_1(:,i);
                v = data_in_2(:,j);
                data_out(i,j) = dot(u,v)/(norm(u)*norm(v)); % = cos(theta)
            end
        end
    else
%         data_out = nan(numBins,numBins,numTrials);
%         for k=1:numTrials
%             for i=1:numBins
%                 for j = 1:numBins
%                     u = data_in_1(:,i);
%                     v = data_in_2(:,j,k);
%                     data_out(i,j,k) = dot(u,v)/(norm(u)*norm(v)); % = cos(theta)
%                 end
%             end
%         end
        data_out = nan(numTrials,numBins);
        for k=1:numTrials
            for i=1:numBins
                u = data_in_1(:,i);
                v = data_in_2(:,i,k);
                data_out(k,i) = dot(u,v)/(norm(u)*norm(v)); % = cos(theta)
            end
        end
    end
elseif strcmp(metric,'Pearson')
    if numTrials==1
        data_out = nan(numBins,numBins);
        for i=1:numBins
            for j = 1:numBins
                u = data_in_1(:,i);
                v = data_in_2(:,j);
                data_out(i,j) = corr(u,v,'Type','Pearson','Rows','Complete');
            end
        end
    else
%         data_out = nan(numBins,numBins,numTrials);
%         for k=1:numTrials
%             for i=1:numBins
%                 for j = 1:numBins
%                     u = data_in_1(:,i);
%                     v = data_in_2(:,j,k);
%                     data_out(i,j,k) = corr(u,v,'Type','Pearson','Rows','Complete');
%                 end
%             end
%         end
        data_out = nan(numTrials,numBins);
        for k=1:numTrials
            for i=1:numBins
                u = data_in_1(:,i);
                v = data_in_2(:,i,k);
                data_out(k,i) = corr(u,v,'Type','Pearson','Rows','Complete');
            end
        end
    end
end




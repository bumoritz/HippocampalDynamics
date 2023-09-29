function data_out = getPairwiseCorrelations(data_in)

data_out = nan(size(data_in,1),size(data_in,1),size(data_in,3));
for k=1:size(data_in,3)
    data_out(:,:,k) = corr(data_in(:,:,k)',data_in(:,:,k)','Type','Pearson','Rows','Complete');
end
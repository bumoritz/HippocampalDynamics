function data_out = standardisePrePostSize(data_in,target_size)
% data_in = paq_post.position;

if size(data_in,2)<target_size
    data_out = nan(1,target_size);
    data_out(1,1:size(data_in,2)) = data_in;
elseif size(data_in,2)>target_size
    data_out = data_in(1,1:target_size);
else
    warning('Something went wrong')
end

end
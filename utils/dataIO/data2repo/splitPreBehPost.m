function [data_out_pre,data_out_beh,data_out_post] = splitPreBehPost(data_in,numFrames_pre,numFrames_post,s2pFrames_pre,s2pFrames_beh,s2pFrames_post)
% data_in = s2p.dFF_gm; numFrames_pre = info.scope.numFrames_pre; numFrames_post = info.scope.numFrames_post;

data_out_pre = nan(size(data_in,1),numFrames_pre,'single');
data_out_beh = nan(size(data_in,1),s2pFrames_beh,'single');
if s2pFrames_pre == numFrames_pre
    data_out_pre(:,1:s2pFrames_pre) = data_in(:,1:s2pFrames_pre);
elseif s2pFrames_pre > numFrames_pre
    data_out_pre(:,1:numFrames_pre) = data_in(:,1:numFrames_pre);
elseif s2pFrames_pre < numFrames_pre
    data_out_pre(:,1:s2pFrames_pre) = data_in(:,1:s2pFrames_pre);
end
data_out_beh(:,1:s2pFrames_beh) = data_in(:,s2pFrames_pre+1:s2pFrames_pre+s2pFrames_beh);

data_out_post = nan(size(data_in,1),numFrames_post,'single');
if s2pFrames_post == numFrames_post
    data_out_post(:,1:s2pFrames_post) = data_in(:,end-(s2pFrames_post-1):end);
elseif s2pFrames_post > numFrames_post
    temp = s2pFrames_post - numFrames_post;
    data_out_post(:,1:numFrames_post) = data_in(:,end-(numFrames_post-1)-temp:end-temp);
elseif s2pFrames_post < numFrames_post
    data_out_post(:,1:s2pFrames_post) = data_in(:,end-(s2pFrames_post-1):end);
end

% for troubleshooting:
% [data_out_post(1,1:3);data_out_post(1,end-2:end)]
% [data_in(1,s2pFrames_pre+s2pFrames_beh+1:s2pFrames_pre+s2pFrames_beh+3);data_in(1,end-2:end)]

end




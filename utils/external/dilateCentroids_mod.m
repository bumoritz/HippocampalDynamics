function halo = dilateCentroids_mod(centroids, roi_radius, watershed_width, halo_multiplier, remove_overlap, show_plot, halo_inset)
% LR 2017
% centroids: is a struct, with fields .x and .y
% halo_multiplier: relative to roi_radius, the alpha of a guassian kernal to approxiamte neuropil contamination.

% prepare arrays
img_size = [512,512];
num_points = numel(centroids.x);
all_rois = zeros(img_size(1), img_size(2), num_points);
all_halos = all_rois;
% all_watershed = all_rois;
all_centroids = all_rois;

% % make ROI mask
% [x,y] = meshgrid(1:(roi_radius*2)+1, 1:(roi_radius*2)+1);
% roi_mask = double((x - (roi_radius+1)).^2 + (y - (roi_radius+1)).^2 < roi_radius^2);

% make halo and watershed mask
halo_mask = fspecial('gaussian', halo_multiplier*6, halo_multiplier);
% [x2,y2] = meshgrid(1:size(halo_mask,1), 1:size(halo_mask,2));
% watershed_mask = (x2 - (round(size(halo_mask,1)/2))-1).^2 + (y2 - (round(size(halo_mask,2)/2))-1).^2 < (roi_radius+watershed_width)^2;
% halo_mask(watershed_mask) = 0;
halo_mask = halo_mask./max(max(halo_mask));

% convolve all the centroids with the ROI/halo/watershed kernels
for i = 1:num_points
    % make a centroid point image
    this_img = zeros(img_size(1), img_size(2));
    this_img(centroids.y(i), centroids.x(i)) = 1;
    all_centroids(:,:,i) = this_img;
    
%     % now convolve the point with the kernels
%     all_rois(:,:,i)      = conv2(this_img, roi_mask, 'same'); 
%     all_watershed(:,:,i) = conv2(this_img, double(watershed_mask), 'same'); 
    all_halos(:,:,i)     = conv2(this_img, halo_mask, 'same');
end

% % remove all ROI pixels (plus watershed) from halos
% max_watershed = sum(all_watershed, 3) > 0;
% for i = 1:num_points
%     temp = all_halos(:,:,i);
%     temp(max_watershed) = 0;
%     all_halos(:,:,i) = temp;
% end
% 
% % remove overlapping pixels from all rois
% if remove_overlap
%     overlap = sum(all_rois, 3) > 1;
%     for i = 1:num_points
%         temp = all_rois(:,:,i);
%         temp(overlap) = 0;
%         all_rois(:,:,i) = temp;
%     end
% end

% plot the results
if show_plot
%     max_roi = max(all_rois,[],3);
    max_halo = max(all_halos,[],3);
    max_centroid = max(all_centroids,[],3);
    figure; imagesc(max_roi + (max_halo/2)); axis square; axis off
end

% save the results
% roi = cell(num_points,1);
halo = cell(num_points,1);
for i = 1:num_points
    % roi
%     [coords_y, coords_x, weights] = find(all_rois(:,:,i) > 0);
%     roi{i}.coords = sub2ind(img_size, coords_y, coords_x);
%     roi{i}.weights = weights;
%     roi{i}.image = all_rois(:,:,i);
    
    % halo
    [coords_x, coords_y, weights] = find(all_halos(:,:,i) > 0.05);
    temp =  find((coords_x > halo_inset) & (coords_x < 512-halo_inset) & (coords_y > halo_inset) & (coords_y < 512-halo_inset));
    coords_x_filt = coords_x(temp);
    coords_y_filt = coords_y(temp);
    weights_filt = weights(temp);
    
    halo{i}.coords = sub2ind(img_size, coords_y_filt, coords_x_filt);
    halo{i}.image = all_halos(:,:,i);
    halo{i}.image(1:halo_inset,:) = 0;
    halo{i}.image(512-halo_inset+1:512,:) = 0;
    halo{i}.image(:,1:halo_inset) = 0;
    halo{i}.image(:,512-halo_inset+1:512) = 0;
    temp2 = halo{i}.image';
    temp = temp2(:);
    halo{i}.weights = temp(halo{i}.coords);
end











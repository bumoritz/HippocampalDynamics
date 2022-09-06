function [path,trg] = repo2repo_trg(info,ops,p,path,regmode)

%% Load s2p_meta files

disp('--- Loading s2p_meta and trg_online files...')

number_of_sessions = 2;
file_names = {};
spatial_footprints = {};
includedRois = {};
for i=1:2
    if i==1
        load([path.filepart_out,'meta.mat']);
    end
    if i==1
        file_names{i} = [path.filepart_out,'s2p_meta.mat'];
    elseif i==2
        file_names{i} = [path.filepart_outX,'trg_online.mat'];
    end
    load(file_names{i});
    if i==1
        includedRois{i} = s2p_meta.iscell(:,2)>p.trg.minCellProb;
    elseif i==2
        includedRois{i} = trg_online.s2p.iscell(:,1);
    end
    spatial_footprints{i} = zeros(sum(includedRois{i}),meta.info.scope.fovSize_pix,meta.info.scope.fovSize_pix);
    n=0;
    if i==1
        for j=1:length(s2p_meta.stat)
            if includedRois{i}(j)
                n=n+1;
                for k=1:length(s2p_meta.stat{j}.lam)
                    spatial_footprints{i}(n,s2p_meta.stat{j}.ypix(k),s2p_meta.stat{j}.xpix(k)) = s2p_meta.stat{j}.lam(k);
                end
            end
        end
    elseif i==2
        for j=1:length(trg_online.s2p.stat)
            if includedRois{i}(j)
                n=n+1;
                for k=1:length(trg_online.s2p.stat{j}.lam)
                    spatial_footprints{i}(n,trg_online.s2p.stat{j}.ypix(k),trg_online.s2p.stat{j}.xpix(k)) = trg_online.s2p.stat{j}.lam(k);
                end
            end
        end
    end
end


%% --- CellReg ---

disp('--- Running CellReg...')


%% Setting paths for the cell registration procedure:

% Defining the results_directory and creating the figures_directory:
figures_directory           = [path.filepart_out,'plots/trck_',regmode];
if ops.trg.showFigures
    figures_visibility      = 'on'; % either 'on' or 'off' (in any case figures are saved)
else
    figures_visibility      = 'off'; % either 'on' or 'off' (in any case figures are saved)
end
if ~exist(figures_directory,'dir')
    mkdir(figures_directory);
end


%% Stage 1 - Loading the spatial footprints of cellular activity:
% This stage loads a new data set which includes several sessions with the
% identified spatial footprints.

% Defining the parameters:
microns_per_pixel           = meta.info.scope.fovSize_um/meta.info.scope.fovSize_pix;

% Loading the data:
disp('Stage 1 - Loading sessions')
%[spatial_footprints,number_of_sessions]=load_multiple_sessions(file_names);
[footprints_projections]=compute_footprints_projections(spatial_footprints);
plot_all_sessions_projections(footprints_projections,figures_directory,figures_visibility)
disp('Done')


%% Stage 2 - Aligning all the sessions to a reference coordinate system:
% A rigid-body transfomration is applied to all the sessions
% according to a chosen reference ssseion. The alignment includes:
% 1. Preparing the data for alignment
% 2. Aligning all the sessions according to a reference coordinate system
% 3. Evaluating how suitable the data is for longitudinal analysis

% Defining the parameters for image alignment:
if strcmp(regmode,'rigid')
    alignment_type          = 'Translations and Rotations';
elseif strcmp(regmode,'nonrigid')
    alignment_type          = 'Non-rigid';
end
use_parallel_processing     = p.trg.use_parallel_processing; % either true or false
maximal_rotation            = p.trg.maximal_rotation; % in degrees - only relevant if 'Translations and Rotations' is used
transformation_smoothness   = p.trg.transformation_smoothness; % levels of non-rigid FOV transformation smoothness (range 0.5-3)
if isnumeric(p.trg.reference_session_index)
    reference_session_index = p.trg.reference_session_index; 
elseif strcmp(p.trg.reference_session_index,'middle')
    reference_session_index = ceil(number_of_sessions/2); 
end

% Preparing the data for alignment:
disp('Stage 2 - Aligning sessions')
[normalized_spatial_footprints]=normalize_spatial_footprints(spatial_footprints);
[adjusted_spatial_footprints,adjusted_FOV,adjusted_x_size,adjusted_y_size,adjustment_zero_padding]=...
    adjust_FOV_size(normalized_spatial_footprints);
[adjusted_footprints_projections]=compute_footprints_projections(adjusted_spatial_footprints);
[centroid_locations]=compute_centroid_locations(adjusted_spatial_footprints,microns_per_pixel); 
[centroid_projections]=compute_centroids_projections(centroid_locations,adjusted_spatial_footprints);

% Aligning the cells according to the tranlations/rotations that maximize their similarity:
sufficient_correlation_centroids=0.2; % smaller correlation imply no similarity between sessions
sufficient_correlation_footprints=0.3; % smaller correlation imply no similarity between sessions
if strcmp(alignment_type,'Translations and Rotations')
    [spatial_footprints_corrected,centroid_locations_corrected,footprints_projections_corrected,centroid_projections_corrected,maximal_cross_correlation,alignment_translations,overlapping_FOV]=...
        align_images(adjusted_spatial_footprints,centroid_locations,adjusted_footprints_projections,centroid_projections,adjusted_FOV,microns_per_pixel,reference_session_index,alignment_type,sufficient_correlation_centroids,sufficient_correlation_footprints,use_parallel_processing,maximal_rotation);
elseif strcmp(alignment_type,'Non-rigid')
    [spatial_footprints_corrected,centroid_locations_corrected,footprints_projections_corrected,centroid_projections_corrected,maximal_cross_correlation,alignment_translations,overlapping_FOV,displacement_fields]=...
        align_images(adjusted_spatial_footprints,centroid_locations,adjusted_footprints_projections,centroid_projections,adjusted_FOV,microns_per_pixel,reference_session_index,alignment_type,sufficient_correlation_centroids,sufficient_correlation_footprints,use_parallel_processing,transformation_smoothness);
else
    [spatial_footprints_corrected,centroid_locations_corrected,footprints_projections_corrected,centroid_projections_corrected,maximal_cross_correlation,alignment_translations,overlapping_FOV]=...
        align_images(adjusted_spatial_footprints,centroid_locations,adjusted_footprints_projections,centroid_projections,adjusted_FOV,microns_per_pixel,reference_session_index,alignment_type,sufficient_correlation_centroids,sufficient_correlation_footprints,use_parallel_processing);
end

% Evaluating data quality:
[all_projections_correlations,number_of_cells_per_session]=...
    evaluate_data_quality(spatial_footprints_corrected,centroid_projections_corrected,footprints_projections_corrected,maximal_cross_correlation,alignment_translations,reference_session_index,sufficient_correlation_footprints,alignment_type);

% plotting alignment results:
if strcmp(alignment_type,'Non-rigid')
    plot_alignment_results(adjusted_spatial_footprints,centroid_locations,spatial_footprints_corrected,centroid_locations_corrected,adjusted_footprints_projections,footprints_projections_corrected,reference_session_index,all_projections_correlations,maximal_cross_correlation,alignment_translations,overlapping_FOV,alignment_type,number_of_cells_per_session,figures_directory,figures_visibility,displacement_fields)
else
    plot_alignment_results(adjusted_spatial_footprints,centroid_locations,spatial_footprints_corrected,centroid_locations_corrected,adjusted_footprints_projections,footprints_projections_corrected,reference_session_index,all_projections_correlations,maximal_cross_correlation,alignment_translations,overlapping_FOV,alignment_type,number_of_cells_per_session,figures_directory,figures_visibility)
end

if use_parallel_processing
    delete(gcp);
end
disp('Done')


%% Stage 3 (part a) - Calculating the similarities distributions from the data:
% This stage uses the ditribtuions of centroid distance and spatial correlations
% to compute the probabilities of neighboring cell-pairs to be the same cell (P_same).

% part a includes the calculation of the distributions of centroid distances and spatial
% correlations from the data.

% Defining the parameters for the probabilstic modeling:
maximal_distance            = p.trg.maximal_distance; % cell-pairs that are more than 12 micrometers apart are assumed to be different cells
p_same_certainty_threshold  = p.trg.p_same_certainty_threshold; % certain cells are those with p_same>threshld or <1-threshold
normalized_maximal_distance=maximal_distance/microns_per_pixel;

% Computing correlations and distances across days:
disp('Stage 3 - Calculating a probabilistic model of the data')
[number_of_bins,centers_of_bins]=estimate_number_of_bins(spatial_footprints,normalized_maximal_distance);
[all_to_all_indexes,all_to_all_spatial_correlations,all_to_all_centroid_distances,neighbors_spatial_correlations,neighbors_centroid_distances,neighbors_x_displacements,neighbors_y_displacements,NN_spatial_correlations,NNN_spatial_correlations,NN_centroid_distances,NNN_centroid_distances]=...
    compute_data_distribution(spatial_footprints_corrected,centroid_locations_corrected,normalized_maximal_distance);

% Plotting the (x,y) displacements:
plot_x_y_displacements(neighbors_x_displacements,neighbors_y_displacements,microns_per_pixel,normalized_maximal_distance,number_of_bins,centers_of_bins,figures_directory,figures_visibility);
disp('Part a done')


%% Stage 3 (part b) - Compute a probabilistic model:
% Modeling the data as a weighted sum of same cells and different cells,
% and estimating the attainable registration accuracy:

disp('Calculating a probabilistic model of the data')
% Modeling the distribution of centroid distances:
[centroid_distances_model_parameters,p_same_given_centroid_distance,centroid_distances_distribution,centroid_distances_model_same_cells,centroid_distances_model_different_cells,centroid_distances_model_weighted_sum,MSE_centroid_distances_model,centroid_distance_intersection]=...
    compute_centroid_distances_model(neighbors_centroid_distances,microns_per_pixel,centers_of_bins);

% Modeling the distribution of spatial correlations:
[spatial_correlations_model_parameters,p_same_given_spatial_correlation,spatial_correlations_distribution,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,spatial_correlations_model_weighted_sum,MSE_spatial_correlations_model,spatial_correlation_intersection]=...
    compute_spatial_correlations_model(neighbors_spatial_correlations,centers_of_bins);

% estimating registration accuracy:
[p_same_centers_of_bins,uncertain_fraction_centroid_distances,cdf_p_same_centroid_distances,false_positive_per_distance_threshold,true_positive_per_distance_threshold,uncertain_fraction_spatial_correlations,cdf_p_same_spatial_correlations,false_positive_per_correlation_threshold,true_positive_per_correlation_threshold]=...
    estimate_registration_accuracy(p_same_certainty_threshold,neighbors_centroid_distances,centroid_distances_model_same_cells,centroid_distances_model_different_cells,p_same_given_centroid_distance,centers_of_bins,neighbors_spatial_correlations,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,p_same_given_spatial_correlation);
% Checking which model is better according to a defined cost function:
[best_model_string]=choose_best_model(MSE_centroid_distances_model,centroid_distances_model_same_cells,centroid_distances_model_different_cells,p_same_given_centroid_distance,MSE_spatial_correlations_model,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,p_same_given_spatial_correlation);

% Plotting the probabilistic models and estimated registration accuracy:
plot_models(centroid_distances_model_parameters,NN_centroid_distances,NNN_centroid_distances,centroid_distances_distribution,centroid_distances_model_same_cells,centroid_distances_model_different_cells,centroid_distances_model_weighted_sum,centroid_distance_intersection,centers_of_bins,microns_per_pixel,normalized_maximal_distance,figures_directory,figures_visibility,spatial_correlations_model_parameters,NN_spatial_correlations,NNN_spatial_correlations,spatial_correlations_distribution,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,spatial_correlations_model_weighted_sum,spatial_correlation_intersection)
plot_estimated_registration_accuracy(p_same_centers_of_bins,p_same_certainty_threshold,p_same_given_centroid_distance,centroid_distances_distribution,cdf_p_same_centroid_distances,uncertain_fraction_centroid_distances,true_positive_per_distance_threshold,false_positive_per_distance_threshold,centers_of_bins,normalized_maximal_distance,microns_per_pixel,figures_directory,figures_visibility,p_same_given_spatial_correlation,spatial_correlations_distribution,cdf_p_same_spatial_correlations,uncertain_fraction_spatial_correlations,true_positive_per_correlation_threshold,false_positive_per_correlation_threshold)

% Computing the P_same for each neighboring cell-pair according to the different models:
[all_to_all_p_same_centroid_distance_model,all_to_all_p_same_spatial_correlation_model]=...
    compute_p_same(all_to_all_centroid_distances,p_same_given_centroid_distance,centers_of_bins,all_to_all_spatial_correlations,p_same_given_spatial_correlation);

disp('Done')


%% Stage 4 - Initial cell registration
% This stage performs an initial cell registration according to an
% optimized threshold of either spatial correlations or centroid distances.

% Defining the parameters for initial registration:
initial_registration_type   = p.trg.initial_registration_type; % either 'Spatial correlation', 'Centroid distance', or 'best_model_string';
% The threshold that corresponds to p_same=0.5 is automatically chosen.
% if a specific distance/correlation threshold is to be used - change the
% initial threshold manually in the next few lines.

% Computing the initial registration according to a simple threshold:
disp('Stage 4 - Performing initial registration')
if strcmp(initial_registration_type,'Spatial correlation') % if spatial correlations are used
    if exist('spatial_correlation_intersection','var')
        initial_threshold=spatial_correlation_intersection; % the threshold for p_same=0.5;
    else
        initial_threshold=0.65; % a fixed correlation threshold not based on the model
    end   
        [cell_to_index_map,registered_cells_spatial_correlations,non_registered_cells_spatial_correlations]=...
            initial_registration_spatial_correlations(normalized_maximal_distance,initial_threshold,spatial_footprints_corrected,centroid_locations_corrected);
        plot_initial_registration(cell_to_index_map,number_of_bins,spatial_footprints_corrected,initial_registration_type,figures_directory,figures_visibility,registered_cells_spatial_correlations,non_registered_cells_spatial_correlations)
else % if centroid distances are used
    if exist('centroid_distance_intersection','var')
        initial_threshold=centroid_distance_intersection; % the threshold for p_same=0.5;
    else
        initial_threshold=5; % a fixed distance threshold not based on the model
    end
    normalized_distance_threshold=initial_threshold/microns_per_pixel;
    [cell_to_index_map,registered_cells_centroid_distances,non_registered_cells_centroid_distances]=...
        initial_registration_centroid_distances(normalized_maximal_distance,normalized_distance_threshold,centroid_locations_corrected);
    plot_initial_registration(cell_to_index_map,number_of_bins,spatial_footprints_corrected,initial_registration_type,figures_directory,figures_visibility,registered_cells_centroid_distances,non_registered_cells_centroid_distances,microns_per_pixel,normalized_maximal_distance)
end

disp([num2str(size(cell_to_index_map,1)) ' cells were found'])
disp('Done')


%% Stage 5 - Final cell registration:
% This stage performs the final cell registration with a clustering algorithm 
% that is based on the probability model for same cells and different cells. 
% P_same can be either according to centroid distances or spatial
% correlations.

% Defining the parameters for final registration:
registration_approach       = p.trg.registration_approach; % either 'Probabilistic' or 'Simple threshold'
model_type                  = p.trg.model_type; % either 'Spatial correlation' or 'Centroid distance'
p_same_threshold            = p.trg.p_same_threshold; % only relevant if probabilistic approach is used

% Deciding on the registration threshold:
transform_data=false;
if strcmp(registration_approach,'Simple threshold') % only relevant if a simple threshold is used
    if strcmp(model_type,'Spatial correlation')
        if exist('spatial_correlation_intersection','var')
            final_threshold=spatial_correlation_intersection; % the threshold for p_same=0.5;
        else
            final_threshold=0.65; % a fixed correlation threshold not based on the model
        end
    elseif strcmp(model_type,'Centroid distance')
        if exist('centroid_distance_intersection','var')
            final_threshold=centroid_distance_intersection; % the threshold for p_same=0.5;
        else
            final_threshold=5; % a fixed distance threshold not based on the model
        end
        normalized_distance_threshold=(maximal_distance-final_threshold)/maximal_distance;
        transform_data=true;
    end
else
    final_threshold=p_same_threshold;
end

% Registering the cells with the clustering algorithm:
disp('Stage 5 - Performing final registration')
if strcmp(registration_approach,'Probabilistic')    
    if strcmp(model_type,'Spatial correlation')
        [optimal_cell_to_index_map,registered_cells_centroids,cell_scores,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,p_same_registered_pairs]=...
            cluster_cells(cell_to_index_map,all_to_all_p_same_spatial_correlation_model,all_to_all_indexes,normalized_maximal_distance,p_same_threshold,centroid_locations_corrected,registration_approach,transform_data);
    elseif strcmp(model_type,'Centroid distance')
        [optimal_cell_to_index_map,registered_cells_centroids,cell_scores,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,p_same_registered_pairs]=...
            cluster_cells(cell_to_index_map,all_to_all_p_same_centroid_distance_model,all_to_all_indexes,normalized_maximal_distance,p_same_threshold,centroid_locations_corrected,registration_approach,transform_data);
    end
    plot_cell_scores(cell_scores_positive,cell_scores_negative,cell_scores_exclusive,cell_scores,p_same_registered_pairs,figures_directory,figures_visibility)
elseif strcmp(registration_approach,'Simple threshold')
    if strcmp(model_type,'Spatial correlation')
        [optimal_cell_to_index_map,registered_cells_centroids]=...
            cluster_cells(cell_to_index_map,all_to_all_spatial_correlations,all_to_all_indexes,normalized_maximal_distance,final_threshold,centroid_locations_corrected,registration_approach,transform_data);
    elseif strcmp(model_type,'Centroid distance')
        [optimal_cell_to_index_map,registered_cells_centroids]=...
            cluster_cells(cell_to_index_map,all_to_all_centroid_distances,all_to_all_indexes,normalized_maximal_distance,normalized_distance_threshold,centroid_locations_corrected,registration_approach,transform_data);
    end
end
[is_in_overlapping_FOV]=check_if_in_overlapping_FOV(registered_cells_centroids,overlapping_FOV);

% Plotting the registration results with the cell maps from all sessions:
plot_all_registered_projections(spatial_footprints_corrected,optimal_cell_to_index_map,figures_directory,figures_visibility)


%% Saving

% saving the final registration results:
disp('Saving the results')
cell_registered_struct=struct;
cell_registered_struct.cell_to_index_map=optimal_cell_to_index_map;
if strcmp(registration_approach,'Probabilistic')
    cell_registered_struct.cell_scores=cell_scores';
    cell_registered_struct.true_positive_scores=cell_scores_positive';
    cell_registered_struct.true_negative_scores=cell_scores_negative';
    cell_registered_struct.exclusivity_scores=cell_scores_exclusive';
    cell_registered_struct.p_same_registered_pairs=p_same_registered_pairs';
end
cell_registered_struct.is_cell_in_overlapping_FOV=is_in_overlapping_FOV';
cell_registered_struct.registered_cells_centroids=registered_cells_centroids';
cell_registered_struct.centroid_locations_corrected=centroid_locations_corrected';
cell_registered_struct.spatial_footprints_corrected=spatial_footprints_corrected';
cell_registered_struct.spatial_footprints_corrected=spatial_footprints_corrected';
cell_registered_struct.alignment_x_translations=alignment_translations(1,:);
cell_registered_struct.alignment_y_translations=alignment_translations(2,:);
if strcmp(alignment_type,'Translations and Rotations')
    cell_registered_struct.alignment_rotations=alignment_translations(3,:);
end
cell_registered_struct.adjustment_x_zero_padding=adjustment_zero_padding(1,:);
cell_registered_struct.adjustment_y_zero_padding=adjustment_zero_padding(2,:);
%save(fullfile(results_directory,['cellRegistered_' datestr(clock,'yyyymmdd_HHMMss') '.mat']),'cell_registered_struct','-v7.3')

% Saving a log file with all the chosen parameters:
comments=''; % anything written here will be added to the log file
if strcmp(registration_approach,'Probabilistic')
    if strcmp(model_type,'Spatial correlation')
        save_log_file([path.filepart_out,'log'],file_names,microns_per_pixel,adjusted_x_size,adjusted_y_size,alignment_type,reference_session_index,maximal_distance,number_of_bins,initial_registration_type,initial_threshold,registration_approach,model_type,final_threshold,optimal_cell_to_index_map,cell_registered_struct,comments,uncertain_fraction_spatial_correlations,false_positive_per_correlation_threshold,true_positive_per_correlation_threshold,MSE_spatial_correlations_model)
    elseif strcmp(model_type,'Centroid distance')
        save_log_file([path.filepart_out,'log'],file_names,microns_per_pixel,adjusted_x_size,adjusted_y_size,alignment_type,reference_session_index,maximal_distance,number_of_bins,initial_registration_type,initial_threshold,registration_approach,model_type,final_threshold,optimal_cell_to_index_map,cell_registered_struct,comments,uncertain_fraction_centroid_distances,false_positive_per_distance_threshold,true_positive_per_distance_threshold,MSE_centroid_distances_model)
    end
elseif strcmp(registration_approach,'Simple threshold')
    if strcmp(model_type,'Spatial correlation')
        save_log_file([path.filepart_out,'log'],file_names,microns_per_pixel,adjusted_x_size,adjusted_y_size,alignment_type,reference_session_index,maximal_distance,number_of_bins,initial_registration_type,initial_threshold,registration_approach,model_type,final_threshold,optimal_cell_to_index_map,cell_registered_struct,comments)
    elseif strcmp(model_type,'Centroid distance')
        save_log_file([path.filepart_out,'log'],file_names,microns_per_pixel,adjusted_x_size,adjusted_y_size,alignment_type,reference_session_index,maximal_distance,number_of_bins,initial_registration_type,initial_threshold,registration_approach,model_type,final_threshold,optimal_cell_to_index_map,cell_registered_struct,comments)
    end
end
disp([num2str(size(optimal_cell_to_index_map,1)) ' cells were found'])
disp('End of cell registration procedure')


%% ---------------

%% Start creating trg file from CellReg outputs and save detailed CellReg output

cell_registered_struct.cell_to_index_map_original = cell_registered_struct.cell_to_index_map;
cell_registered_struct.includedRois = includedRois;
for i=1:size(cell_registered_struct.cell_to_index_map_original,2)
    temp = find(includedRois{i});
    for j=1:size(cell_registered_struct.cell_to_index_map_original,1)
        if cell_registered_struct.cell_to_index_map_original(j,i)~=0
            cell_registered_struct.cell_to_index_map(j,i) = temp(cell_registered_struct.cell_to_index_map_original(j,i));
        end
    end
end

% save trg_trck_details file
trg_trck_details = orderfields(cell_registered_struct);
save([path.filepart_outX,'trg_trck_details_',regmode,'.mat'],'trg_trck_details','-v7.3');
disp(['--- Saved trg_trck_details_',regmode,' file to repoX as ',[path.filepart_outX,'trg_trck_details_',regmode,'.mat'],'.'])

% start creating trg file
trg = trg_online;
trg = rmfield(trg,'s2p'); trg = rmfield(trg,'xml'); trg = rmfield(trg,'cmp_across'); trg = rmfield(trg,'cmp_within');
trg.trck.cell_to_index_map = cell_registered_struct.cell_to_index_map;
trg.trck.cell_to_index_map_original = cell_registered_struct.cell_to_index_map_original;
trg.trck.cell_scores = cell_registered_struct.cell_scores;
trg.trck.registered_cells_centroids = cell_registered_struct.registered_cells_centroids;
trg.trck = orderfields(trg.trck);


%% Identify targets

% identify laser spots
trg.laser_x = zeros(size(trg.clustering));
trg.laser_y = zeros(size(trg.clustering));
for j=1:size(trg.clustering,2)
    for i=1:size(trg.clustering,1)
        trg.laser_y(i,j) = trg.f.med2(trg.clustering(i,j));
        trg.laser_x(i,j) = trg.f.med1(trg.clustering(i,j));
    end
end

% identify targets
trg.idcs_targetedCells_spontIndexing = zeros(size(trg.clustering));
trg.idcs_targetedCells = zeros(size(trg.clustering));
trg.targeted = false(length(s2p_meta.stat),size(trg.clustering,2));
trg.trck.targeted = false(size(trg.trck.cell_to_index_map,1),size(trg.clustering,2));
for j=1:size(trg.clustering,2)
    for i=1:size(trg.clustering,1)
        trg.idcs_targetedCells_spontIndexing(i,j) = trg.iscells(trg.clustering(i,j));
        
        if ismember(trg.idcs_targetedCells_spontIndexing(i,j),trg.trck.cell_to_index_map(:,2))
            temp = find(trg.trck.cell_to_index_map(:,2)==trg.idcs_targetedCells_spontIndexing(i,j));
            trg.trck.targeted(temp,j) = true;
            if trg.trck.cell_to_index_map(temp,1)~=0
                trg.idcs_targetedCells(i,j) = trg.trck.cell_to_index_map(temp,1);
                trg.targeted(trg.idcs_targetedCells(i,j),j) = true;
            else
                trg.idcs_targetedCells(i,j) = NaN;
            end
        else
            trg.idcs_targetedCells(i,j) = NaN;
        end
    end
end
trg.numIdentifiedTargeted = sum(~isnan(trg.idcs_targetedCells(:)));
disp(['--- Re-identified targets: ',num2str(trg.numIdentifiedTargeted),' / ',num2str(length(trg.idcs_targetedCells(:))),'.'])

plot_all_registered_projections_MB(spatial_footprints_corrected,optimal_cell_to_index_map,figures_directory,figures_visibility,trg);


%% Calculate distance to closest laser and closest target

% EXAMPLE HOW TO INDEX
% idx = trg.idcs_targetedCells(1,1);
% temp = find(trg.trck.cell_to_index_map(:,1)==idx);
% % WRONG: [trg.f.med2(trg.iscells(trg.trck.cell_to_index_map(temp,2))),trg.f.med1(trg.iscells(trg.trck.cell_to_index_map(temp,2)))]
% [trg.f.med2(trg.trck.cell_to_index_map_original(temp,2)),trg.f.med1(trg.trck.cell_to_index_map_original(temp,2))]
% trg.trck.registered_cells_centroids(temp,:)

% calculate distance to closest laser
trg.dist_closestLaser = zeros(length(s2p_meta.iscell(:,1)),size(trg.idcs_targetedCells,2));
for j=1:size(trg.idcs_targetedCells,2)
    temp3 = find(trg.trck.cell_to_index_map(:,1)~=0);
    if length(temp3)~=length(s2p_meta.iscell(:,1))
        error('Indexing issue')
    end
    temp = pdist2(trg.trck.registered_cells_centroids(temp3,:),[trg.laser_y(:,j),trg.laser_x(:,j)],'euclidean');
    temp4 = sortrows([trg.trck.cell_to_index_map(temp3,1),nanmin(temp,[],2)],1);
    trg.dist_closestLaser(:,j) = temp4(:,2) * (meta.info.scope.fovSize_um/meta.info.scope.fovSize_pix); 
end

% calculate distance to closest target
trg.dist_closestTargeted = nan(length(s2p_meta.iscell(:,1)),size(trg.idcs_targetedCells,2));
for j=1:size(trg.idcs_targetedCells,2)
    temp1 = trg.idcs_targetedCells(~isnan(trg.idcs_targetedCells(:,j)),j);
    temp2 = [];
    for i=1:length(temp1)
        temp2 = [temp2; find(trg.trck.cell_to_index_map(:,1)==temp1(i))];
    end
    temp3 = find(trg.trck.cell_to_index_map(:,1)~=0);
    if length(temp3)~=length(s2p_meta.iscell(:,1))
        error('Indexing issue')
    end
    temp = pdist2(trg.trck.registered_cells_centroids(temp3,:),trg.trck.registered_cells_centroids(temp2,:),'euclidean');
    temp4 = sortrows([trg.trck.cell_to_index_map(temp3,1),nanmin(temp,[],2)],1);
    try
        trg.dist_closestTargeted(:,j) = temp4(:,2) * (meta.info.scope.fovSize_um/meta.info.scope.fovSize_pix);
    catch
    end
end


%% Save

trg = orderfields(trg);
save([path.filepart_out,'trg_',regmode,'.mat'],'trg','-v7.3');
disp(['--- Saved trg_',regmode,' file to repo as ',[path.filepart_out,'trg_',regmode,'.mat'],'.'])

end

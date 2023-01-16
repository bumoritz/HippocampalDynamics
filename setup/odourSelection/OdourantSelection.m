%% Load and descriptors

path_descriptors = 'C:\Users\Moritz\Desktop\Odourant Similarity Analysis\FinalOdourSpace\OdourSpace_DESCRIPTORS.xls';

temp = sheetnames(path_descriptors);
temp2 = [];
for i=1:length(temp)
    temp2 = [temp2,readmatrix(path_descriptors,'Sheet',temp(i))];
end
data_raw = temp2(2:end,:);


%% Load and extract local odour info

path_local = 'C:\Users\Moritz\Desktop\Odourant Similarity Analysis\FinalOdourSpace\OdourSpace_local.xls';

local = readtable(path_local);

idcs_local = local.Global;
idcs_candidates = local.Global(~isnan(local.Candidates));
labels_candidates = local.ShortName(~isnan(local.Candidates));

[~,sorter] = sort(local.SampleOdour(~isnan(local.SampleOdour)));
idcs_sampleOdours = local.Global(~isnan(local.SampleOdour));
idcs_sampleOdours = idcs_sampleOdours(sorter);
[~,~,idcs_sampleOdours_rel_canditates] = intersect(idcs_sampleOdours,idcs_candidates);
idcs_sampleOdours_rel_canditates = sort(idcs_sampleOdours_rel_canditates);
idcs_sampleOdours_rel_canditates = idcs_sampleOdours_rel_canditates(sorter);

labels_candidates = local.ShortName(~isnan(local.Candidates));
labels_sampleOdours = local.ShortName(~isnan(local.SampleOdour));
labels_sampleOdours = labels_sampleOdours(sorter);


%% Global odour space

[numGlobalOdours, numDescriptors] = size(data_raw);
data_zscored = zscore(data_raw,0,1);
[pca_coeff,pca_score,pca_latent,pca_tsquared,pca_explained,pca_mu] = pca(data_zscored);

corr_global = zeros(numGlobalOdours);
for i=1:numGlobalOdours
    for j=1:numGlobalOdours
        corr_global(i,j) = corr(data_zscored(i,:)',data_zscored(j,:)','Type','Pearson','Rows','Complete');
    end
end
dist_global = 1-corr_global;

numCandidates = length(idcs_candidates);
corr_candidates = zeros(numCandidates);
for i=1:numCandidates
    for j=1:numCandidates
        corr_candidates(i,j) = corr(data_zscored(idcs_candidates(i),:)',data_zscored(idcs_candidates(j),:)','Type','Pearson','Rows','Complete');
    end
end
dist_candidates = 1-corr_candidates;

numSampleOdours = length(idcs_sampleOdours);
corr_sampleOdours = zeros(numSampleOdours);
for i=1:numSampleOdours
    for j=1:numSampleOdours
        corr_sampleOdours(i,j) = corr(data_zscored(idcs_sampleOdours(i),:)',data_zscored(idcs_sampleOdours(j),:)','Type','Pearson','Rows','Complete');
    end
end
dist_sampleOdours = 1-corr_sampleOdours;


%% Plot - Global odour space - PCA

figure;
scatter(pca_score(:,1),-pca_score(:,2),'.','MarkerEdgeColor',[0.8,0.8,0.8])
hold on
scatter(pca_score(idcs_local,1),-pca_score(idcs_local,2),'ko','LineWidth',2)
scatter(pca_score(idcs_candidates,1),-pca_score(idcs_candidates,2),'ro','LineWidth',2)
scatter(pca_score(idcs_sampleOdours,1),-pca_score(idcs_sampleOdours,2),'gx','LineWidth',2)

xlabel('1st PC (33% of variance)')
ylabel('2nd PC (10% of variance)')
annotation('textbox',[0.25 0 0.3 0.3],'String',[num2str(numGlobalOdours),' odourants',newline,num2str(numDescriptors),' descriptors'],'FitBoxToText','on');
title('Global odour space')


%% Plot - Global odour space - UMAP embedding 

[reduction,umap,clusterIdentifiers,extras]=run_umap(data_zscored,'n_components',2,'metric','correlation');


%% Plot - Global odour space - Clustered heatmap

cgo = clustergram(dist_global,'Colormap','pink','Symmetric',false);


%% Plot - Candidate odours - Clustered heatmap

cgo = clustergram(dist_candidates,'Colormap','pink','Symmetric',false,...
    'ColumnLabels',labels_candidates,'RowLabels',labels_candidates);


%% Plot - Sample odours - Heatmap

figure;
heatmap(dist_sampleOdours);
colormap('pink')
caxis([0,0.2])
set(gca,'XData',{'Ai','Xi','Aii','Xii','Aiii','Xiii'})
set(gca,'YData',labels_sampleOdours)
title('Chemical distances between sample odourants')



%% k-means clustering

% k = 9;
% [kmeans_idx,kmeans_C,kmeans_sumd,kmeans_D] = kmeans(data,k,'Replicates',1000,'Options',statset('Display','final'));


%% Maximise differences within and across set

numSelectedOdours = 12;

data = data_zscored(idcs_candidates,:);

% list all possible combinations
numAdditionalOdours = numSelectedOdours - length(idcs_sampleOdours_rel_canditates);
possibleAdditionalOdours = setdiff(1:numCandidates,idcs_sampleOdours_rel_canditates);
combs.combinationList = nchoosek(possibleAdditionalOdours,numAdditionalOdours);
combs.numCombinations = size(combs.combinationList,1);

% go through every possible combination
combs.combinationMeanDistance_all = zeros(combs.numCombinations,1);
combs.combinationMeanDistance_woDiagonal = zeros(combs.numCombinations,1);
combs.combinationMeanDistance_woFirstOdours = zeros(combs.numCombinations,1);
combs.combinationMeanDistance_woFirstOdoursAcross = zeros(combs.numCombinations,1);
combs.combinationMinDistance_all = zeros(combs.numCombinations,1);
combs.combinationMinDistance_woDiagonal = zeros(combs.numCombinations,1);
combs.combinationMinDistance_woFirstOdours = zeros(combs.numCombinations,1);
combs.combinationMinDistance_woFirstOdoursAcross = zeros(combs.numCombinations,1);

for c=1:combs.numCombinations

    this_selection = combs.combinationList(c,:);
    these_odours = [idcs_sampleOdours_rel_canditates',this_selection];
    this_data = data(these_odours,:);
    
    % calculate correlation matrix
    these_correlations = zeros(numSelectedOdours);
    for i=1:numSelectedOdours
        for j=1:numSelectedOdours
            these_correlations(i,j) = corr(this_data(i,:)',this_data(j,:)','Type','Pearson','Rows','Complete');
        end
    end
    these_distances = 1-these_correlations;
   
%     figure;
%     heatmap(these_distances);
%     colormap('pink');
%     title('chemical distance');

    % remove irrelevant entries
    these_distances_woDiagonal = these_distances;
    these_distances_woDiagonal(1:numSelectedOdours+1:end) = diag(NaN(numSelectedOdours));
    these_distances_woFirstOdours = these_distances_woDiagonal;
    for k=0:length(idcs_sampleOdours_rel_canditates)/2-1
        these_distances_woFirstOdours(2*k+1,2*k+2) = NaN;
        these_distances_woFirstOdours(2*k+2,2*k+1) = NaN;
    end
    these_distances_woFirstOdoursAcross = these_distances_woDiagonal;
    these_distances_woFirstOdoursAcross(1:length(idcs_sampleOdours_rel_canditates),1:length(idcs_sampleOdours_rel_canditates)) = NaN;

    combs.combinationMeanDistance_all(c) = nanmean(these_distances(:));
    combs.combinationMeanDistance_woDiagonal(c) = nanmean(these_distances_woDiagonal(:));
    combs.combinationMeanDistance_woFirstOdours(c) = nanmean(these_distances_woFirstOdours(:));
    combs.combinationMeanDistance_woFirstOdoursAcross(c) = nanmean(these_distances_woFirstOdoursAcross(:));
    
    combs.combinationMinDistance_all(c) = nanmin(these_distances(:));
    combs.combinationMinDistance_woDiagonal(c) = nanmin(these_distances_woDiagonal(:));
    combs.combinationMinDistance_woFirstOdours(c) = nanmin(these_distances_woFirstOdours(:));
    combs.combinationMinDistance_woFirstOdoursAcross(c) = nanmin(these_distances_woFirstOdoursAcross(:));
end


%% Maximise differences for MAGNUM OPUS

numSelectedOdours = 6;

data = data_zscored(idcs_candidates,:);

% list all possible combinations
numAdditionalOdours = numSelectedOdours;
possibleAdditionalOdours = 1:numCandidates;
combs.combinationList = nchoosek(possibleAdditionalOdours,numAdditionalOdours);
combs.numCombinations = size(combs.combinationList,1);

% go through every possible combination
combs.combinationMeanDistance_all = zeros(combs.numCombinations,1);
combs.combinationMeanDistance_woDiagonal = zeros(combs.numCombinations,1);
combs.combinationMeanDistance_woFirstOdours = zeros(combs.numCombinations,1);
combs.combinationMeanDistance_woFirstOdoursAcross = zeros(combs.numCombinations,1);
combs.combinationMinDistance_all = zeros(combs.numCombinations,1);
combs.combinationMinDistance_woDiagonal = zeros(combs.numCombinations,1);
combs.combinationMinDistance_woFirstOdours = zeros(combs.numCombinations,1);
combs.combinationMinDistance_woFirstOdoursAcross = zeros(combs.numCombinations,1);

for c=1:combs.numCombinations

    this_selection = combs.combinationList(c,:);
    these_odours = this_selection;
    this_data = data(these_odours,:);
    
    % calculate correlation matrix
    these_correlations = zeros(numSelectedOdours);
    for i=1:numSelectedOdours
        for j=1:numSelectedOdours
            these_correlations(i,j) = corr(this_data(i,:)',this_data(j,:)','Type','Pearson','Rows','Complete');
        end
    end
    these_distances = 1-these_correlations;
   
%     figure;
%     heatmap(these_distances);
%     colormap('pink');
%     title('chemical distance');

    % remove irrelevant entries
    these_distances_woDiagonal = these_distances;
    these_distances_woDiagonal(1:numSelectedOdours+1:end) = diag(NaN(numSelectedOdours));
%     these_distances_woFirstOdours = these_distances_woDiagonal;
%     for k=0:length(idcs_sampleOdours_rel_canditates)/2-1
%         these_distances_woFirstOdours(2*k+1,2*k+2) = NaN;
%         these_distances_woFirstOdours(2*k+2,2*k+1) = NaN;
%     end
%     these_distances_woFirstOdoursAcross = these_distances_woDiagonal;
%     these_distances_woFirstOdoursAcross(1:length(idcs_sampleOdours_rel_canditates),1:length(idcs_sampleOdours_rel_canditates)) = NaN;

    combs.combinationMeanDistance_all(c) = nanmean(these_distances(:));
    combs.combinationMeanDistance_woDiagonal(c) = nanmean(these_distances_woDiagonal(:));
%     combs.combinationMeanDistance_woFirstOdours(c) = nanmean(these_distances_woFirstOdours(:));
%     combs.combinationMeanDistance_woFirstOdoursAcross(c) = nanmean(these_distances_woFirstOdoursAcross(:));
    
    combs.combinationMinDistance_all(c) = nanmin(these_distances(:));
    combs.combinationMinDistance_woDiagonal(c) = nanmin(these_distances_woDiagonal(:));
%     combs.combinationMinDistance_woFirstOdours(c) = nanmin(these_distances_woFirstOdours(:));
%     combs.combinationMinDistance_woFirstOdoursAcross(c) = nanmin(these_distances_woFirstOdoursAcross(:));
end

%%%
% find(combs.combinationMinDistance_woDiagonal>0.4045) -> 8 X 0.4045, 
% 5109, 6078, 6894, 8278, 27264, 28233, 29049, 30433
% of those, 6078 maximised the mean distance
% combs.combinationList(6078,:): 1     2     7    10    11    19
% Eucalyptol, IsoamylAcetate, MethylBenzoate, Limonene, Pentanol, Pentanedione


%% Test winning combinations

% [a,b]=max(combs.combinationMeanDistance_woFirstOdoursAcross) 
% -> 12 selected, 26 candidates, global zscore -> 0.6664 at 3738, same for across

% [a,b]=max(combs.combinationMinDistance_woFirstOdoursAcross)
% -> 12 selected, 26 candidates, global zscore -> e.g. 0.0932 at 1, 0.2747 at 228
% 4 equally good solutions: 228, 265, 3623, 3660

% THE FINAL WINNER IS 3660

data = data_zscored(idcs_candidates,:);

selectedCombination = 3738 %3660;

firstOdours = idcs_sampleOdours_rel_canditates;

this_selection = combs.combinationList(selectedCombination,:);
these_odours = [idcs_sampleOdours_rel_canditates(1:2)',this_selection(1:2),...
    idcs_sampleOdours_rel_canditates(3:4)',this_selection(3:4),...
    idcs_sampleOdours_rel_canditates(5:6)',this_selection(5:6)]; %temp
this_data = data(these_odours,:);

% calculate correlation matrix
these_correlations = zeros(numSelectedOdours);
for i=1:numSelectedOdours
    for j=1:numSelectedOdours
        these_correlations(i,j) = corr(this_data(i,:)',this_data(j,:)','Type','Pearson','Rows','Complete');
    end
end
these_distances = 1-these_correlations;
    
% plot heatmap
figure;
heatmap(these_distances);
colormap('pink');
%caxis([0,0.8])
set(gca,'XData',{'Ai','Xi','Bi','Yi','Aii','Xii','Bii','Yii','Aiii','Xiii','Biii','Yiii'})
set(gca,'YData',labels_candidates(these_odours))
title(['chemical distance, combination ',num2str(selectedCombination)]);

% plot heatmap on log scale
figure;
heatmap(these_distances);
colormap('pink');
set(gca,'ColorScaling','log')
temp = these_distances(:);
caxis([log(min(temp(temp>0.001))),log(max(these_distances(:)))])
set(gca,'XData',{'Ai','Xi','Bi','Yi','Aii','Xii','Bii','Yii','Aiii','Xiii','Biii','Yiii'})
set(gca,'YData',labels_candidates(these_odours))
title(['chemical distance, combination ',num2str(selectedCombination)]);


% plot heatmap with cheeky colormap
figure;
heatmap(these_distances);
set(gca,'XData',{'Ai','Xi','Bi','Yi','Aii','Xii','Bii','Yii','Aiii','Xiii','Biii','Yiii'})
set(gca,'YData',labels_candidates(these_odours))
title(['chemical distance, combination ',num2str(selectedCombination)]);

cMap_fundament = colormap('pink');
cMap_max = max(these_distances(:))+max(these_distances(:))/100;
cMap_min = min(these_distances(:));
cMap_center = 0.3;
cMap_scalingIntensity = 4;
temp = 1:length(cMap_fundament); 
temp = temp - (cMap_center-cMap_min)*length(temp)/(cMap_max-cMap_min);
temp = cMap_scalingIntensity * temp/max(abs(temp));
temp = sign(temp).* exp(abs(temp));
temp = temp - min(temp); 
temp = temp*511/max(temp)+1; 
cMap = interp1(temp,cMap_fundament,1:512);
colormap(cMap);


%% Plot - Selection in PC space

figure;
hold on
scatter(pca_score(:,1),-pca_score(:,2),'.','MarkerEdgeColor',[0.8,0.8,0.8])
scatter(pca_score(idcs_local,1),-pca_score(idcs_local,2),'o','MarkerEdgeColor',[0.5,0.5,0.5],'LineWidth',2)
scatter(pca_score(idcs_candidates,1),-pca_score(idcs_candidates,2),'o','MarkerEdgeColor',[0.2,0.2,0.2],'LineWidth',2)

scatter(pca_score(idcs_candidates(these_odours(1)),1),-pca_score(idcs_candidates(these_odours(1)),2),100,'x','MarkerEdgeColor',[1,0,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(2)),1),-pca_score(idcs_candidates(these_odours(2)),2),100,'x','MarkerEdgeColor',[1,0,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(3)),1),-pca_score(idcs_candidates(these_odours(3)),2),100,'d','MarkerEdgeColor',[1,0,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(4)),1),-pca_score(idcs_candidates(these_odours(4)),2),100,'d','MarkerEdgeColor',[1,0,0],'LineWidth',2)

scatter(pca_score(idcs_candidates(these_odours(5)),1),-pca_score(idcs_candidates(these_odours(5)),2),100,'x','MarkerEdgeColor',[0,1,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(6)),1),-pca_score(idcs_candidates(these_odours(6)),2),100,'x','MarkerEdgeColor',[0,1,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(7)),1),-pca_score(idcs_candidates(these_odours(7)),2),100,'d','MarkerEdgeColor',[0,1,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(8)),1),-pca_score(idcs_candidates(these_odours(8)),2),100,'d','MarkerEdgeColor',[0,1,0],'LineWidth',2)

scatter(pca_score(idcs_candidates(these_odours(9)),1),-pca_score(idcs_candidates(these_odours(9)),2),100,'x','MarkerEdgeColor',[0,0,1],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(10)),1),-pca_score(idcs_candidates(these_odours(10)),2),100,'x','MarkerEdgeColor',[0,0,1],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(11)),1),-pca_score(idcs_candidates(these_odours(11)),2),100,'d','MarkerEdgeColor',[0,0,1],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(12)),1),-pca_score(idcs_candidates(these_odours(12)),2),100,'d','MarkerEdgeColor',[0,0,1],'LineWidth',2)

xlabel('1st PC (33% of variance)')
ylabel('2nd PC (10% of variance)')
annotation('textbox',[0.25 0 0.3 0.3],'String',[num2str(numGlobalOdours),' odourants',newline,num2str(numDescriptors),' descriptors'],'FitBoxToText','on');
title('Global odour space')


%% Maximise differences within set

data = data_zscored(idcs_candidates(these_odours),:);

permus.permutationList = perms([1:6]);
permus.numPermutations = size(permus.permutationList,1);
permus.permutationMinDistance = zeros(permus.numPermutations,1);
temp = [3,4,7,8,11,12];

for p=1:permus.numPermutations

    this_data = data([1,2,temp(permus.permutationList(p,1)),temp(permus.permutationList(p,2)),...
        5,6,temp(permus.permutationList(p,3)),temp(permus.permutationList(p,4)),...
        9,10,temp(permus.permutationList(p,5)),temp(permus.permutationList(p,6))],:);
    
    % calculate correlation matrix
    these_correlations = zeros(numSelectedOdours);
    for i=1:numSelectedOdours
        for j=1:numSelectedOdours
            these_correlations(i,j) = corr(this_data(i,:)',this_data(j,:)','Type','Pearson','Rows','Complete');
        end
    end
    these_distances = 1-these_correlations;
   
    these_distances_woDiagonal = these_distances;
    these_distances_woDiagonal(1:numSelectedOdours+1:end) = diag(NaN(numSelectedOdours));

    these_dist_i = these_distances_woDiagonal(1:4,1:4);
    these_dist_ii = these_distances_woDiagonal(5:8,5:8);
    these_dist_iii = these_distances_woDiagonal(9:12,9:12);
    these_dist_i(1:2,1:2) = NaN;
    these_dist_ii(1:2,1:2) = NaN;
    these_dist_iii(1:2,1:2) = NaN;

    this_minCorr_i = nanmin(these_dist_i(:));
    this_minCorr_ii = nanmin(these_dist_ii(:));
    this_minCorr_iii = nanmin(these_dist_iii(:));

    permus.permutationMinDistance(p) = min([this_minCorr_i,this_minCorr_ii,this_minCorr_iii]);
end


%% Test winning permutations

% [a,b]=max(permus.permutationMinDistance)
% find(permus.permutationMinDistance>0.3972)
% permus.permutationList(find(permus.permutationMinDistance>0.3972),:)
% 

% THE FINAL WINNER IS 539

data = data_zscored(idcs_candidates(these_odours),:);

selectedPermutation = 549;

%old_these_odours = these_odours;
temp = [3,4,7,8,11,12];
these_odours = old_these_odours([1,2,temp(permus.permutationList(selectedPermutation,1)),temp(permus.permutationList(selectedPermutation,2)),...
        5,6,temp(permus.permutationList(selectedPermutation,3)),temp(permus.permutationList(selectedPermutation,4)),...
        9,10,temp(permus.permutationList(selectedPermutation,5)),temp(permus.permutationList(selectedPermutation,6))]);
this_data = data([1,2,temp(permus.permutationList(selectedPermutation,1)),temp(permus.permutationList(selectedPermutation,2)),...
        5,6,temp(permus.permutationList(selectedPermutation,3)),temp(permus.permutationList(selectedPermutation,4)),...
        9,10,temp(permus.permutationList(selectedPermutation,5)),temp(permus.permutationList(selectedPermutation,6))],:);

% calculate correlation matrix
these_correlations = zeros(numSelectedOdours);
for i=1:numSelectedOdours
    for j=1:numSelectedOdours
        these_correlations(i,j) = corr(this_data(i,:)',this_data(j,:)','Type','Pearson','Rows','Complete');
    end
end
these_distances = 1-these_correlations;
    
% plot heatmap
figure;
heatmap(these_distances);
colormap('pink');
%caxis([0,0.8])
set(gca,'XData',{'Ai','Xi','Bi','Yi','Aii','Xii','Bii','Yii','Aiii','Xiii','Biii','Yiii'})
set(gca,'YData',labels_candidates(these_odours))
title(['chemical distance, permutation ',num2str(selectedPermutation)]);

% plot heatmap on log scale
figure;
heatmap(these_distances);
colormap('pink');
set(gca,'ColorScaling','log')
temp = these_distances(:);
caxis([log(min(temp(temp>0.001))),log(max(these_distances(:)))])
set(gca,'XData',{'Ai','Xi','Bi','Yi','Aii','Xii','Bii','Yii','Aiii','Xiii','Biii','Yiii'})
set(gca,'YData',labels_candidates(these_odours))
title(['chemical distance, permutation ',num2str(selectedPermutation)]);


% plot heatmap with cheeky colormap
figure;
heatmap(these_distances);
set(gca,'XData',{'Ai','Xi','Bi','Yi','Aii','Xii','Bii','Yii','Aiii','Xiii','Biii','Yiii'})
set(gca,'YData',labels_candidates(these_odours))
title(['chemical distance, permutation ',num2str(selectedPermutation)]);
cMap_fundament = colormap('pink');
cMap_max = max(these_distances(:))+max(these_distances(:))/100;
cMap_min = min(these_distances(:));
cMap_center = 0.3;
cMap_scalingIntensity = 4;
temp = 1:length(cMap_fundament); 
temp = temp - (cMap_center-cMap_min)*length(temp)/(cMap_max-cMap_min);
temp = cMap_scalingIntensity * temp/max(abs(temp));
temp = sign(temp).* exp(abs(temp));
temp = temp - min(temp); 
temp = temp*511/max(temp)+1; 
cMap = interp1(temp,cMap_fundament,1:512);
colormap(cMap);


%% Test all winning permutations

% [a,b]=max(permus.permutationMinDistance)
% find(permus.permutationMinDistance>0.3972)
% permus.permutationList(find(permus.permutationMinDistance>0.3972),:)
% 

data = data_zscored(idcs_candidates(these_odours),:);

temp2 = find(permus.permutationMinDistance>0.3972);
for p=1:size(temp2,1)

    selectedPermutation = temp2(p);

    %old_these_odours = these_odours;
    temp = [3,4,7,8,11,12];
    these_odours = old_these_odours([1,2,temp(permus.permutationList(selectedPermutation,1)),temp(permus.permutationList(selectedPermutation,2)),...
            5,6,temp(permus.permutationList(selectedPermutation,3)),temp(permus.permutationList(selectedPermutation,4)),...
            9,10,temp(permus.permutationList(selectedPermutation,5)),temp(permus.permutationList(selectedPermutation,6))]);
    this_data = data([1,2,temp(permus.permutationList(selectedPermutation,1)),temp(permus.permutationList(selectedPermutation,2)),...
            5,6,temp(permus.permutationList(selectedPermutation,3)),temp(permus.permutationList(selectedPermutation,4)),...
            9,10,temp(permus.permutationList(selectedPermutation,5)),temp(permus.permutationList(selectedPermutation,6))],:);

    % calculate correlation matrix
    these_correlations = zeros(numSelectedOdours);
    for i=1:numSelectedOdours
        for j=1:numSelectedOdours
            these_correlations(i,j) = corr(this_data(i,:)',this_data(j,:)','Type','Pearson','Rows','Complete');
        end
    end
    these_distances = 1-these_correlations;

    % % plot heatmap
    % figure;
    % heatmap(these_distances);
    % colormap('pink');
    % %caxis([0,0.8])
    % set(gca,'XData',{'Ai','Xi','Bi','Yi','Aii','Xii','Bii','Yii','Aiii','Xiii','Biii','Yiii'})
    % set(gca,'YData',labels_candidates(these_odours))
    % title(['chemical distance, permutation ',num2str(selectedPermutation)]);

    % plot heatmap on log scale
    figure;
    heatmap(these_distances);
    colormap('pink');
    set(gca,'ColorScaling','log')
    temp = these_distances(:);
    caxis([log(min(temp(temp>0.001))),log(max(these_distances(:)))])
    set(gca,'XData',{'Ai','Xi','Bi','Yi','Aii','Xii','Bii','Yii','Aiii','Xiii','Biii','Yiii'})
    set(gca,'YData',labels_candidates(these_odours))
    title(['chemical distance, permutation ',num2str(selectedPermutation)]);


    % % plot heatmap with cheeky colormap
    % figure;
    % heatmap(these_distances);
    % set(gca,'XData',{'Ai','Xi','Bi','Yi','Aii','Xii','Bii','Yii','Aiii','Xiii','Biii','Yiii'})
    % set(gca,'YData',labels_candidates(these_odours))
    % title(['chemical distance, permutation ',num2str(selectedPermutation)]);
    % cMap_fundament = colormap('pink');
    % cMap_max = max(these_distances(:))+max(these_distances(:))/100;
    % cMap_min = min(these_distances(:));
    % cMap_center = 0.3;
    % cMap_scalingIntensity = 4;
    % temp = 1:length(cMap_fundament); 
    % temp = temp - (cMap_center-cMap_min)*length(temp)/(cMap_max-cMap_min);
    % temp = cMap_scalingIntensity * temp/max(abs(temp));
    % temp = sign(temp).* exp(abs(temp));
    % temp = temp - min(temp); 
    % temp = temp*511/max(temp)+1; 
    % cMap = interp1(temp,cMap_fundament,1:512);
    % colormap(cMap);
end


%% Checking results

% 1) checking for which second odourants come together in same set
% --- combinations to avoid: any two of Pinene-Eucalyptol-Guaiacol, or Guaiacol-Limonene
% --- only remaining option: LE,BP,GI

% 2) checking which first odourants to combine the second odourants with
% --- 2 options: ME-LE,PH-BP,(BF-GI) (e.g. 645) or ME-BP,PH-LE,(BF-GI) (e.g.323)
% --- comparing in PC space: ME-BP,PH-LE,(BF-GI) (e.g.323) wins

% 3) checking which second odourant is B or Y
% --- check in PCA plot that all 4 odours of a set form NOT always the same a cross -> NOT 539
% => 549

%% Plot - Selection in PC space

selectedPermutation = 549;

temp = [3,4,7,8,11,12];
these_odours = old_these_odours([1,2,temp(permus.permutationList(selectedPermutation,1)),temp(permus.permutationList(selectedPermutation,2)),...
        5,6,temp(permus.permutationList(selectedPermutation,3)),temp(permus.permutationList(selectedPermutation,4)),...
        9,10,temp(permus.permutationList(selectedPermutation,5)),temp(permus.permutationList(selectedPermutation,6))]);
this_data = data([1,2,temp(permus.permutationList(selectedPermutation,1)),temp(permus.permutationList(selectedPermutation,2)),...
        5,6,temp(permus.permutationList(selectedPermutation,3)),temp(permus.permutationList(selectedPermutation,4)),...
        9,10,temp(permus.permutationList(selectedPermutation,5)),temp(permus.permutationList(selectedPermutation,6))],:);

figure;
hold on
scatter(pca_score(:,1),-pca_score(:,2),'.','MarkerEdgeColor',[0.8,0.8,0.8])
scatter(pca_score(idcs_local,1),-pca_score(idcs_local,2),'o','MarkerEdgeColor',[0.5,0.5,0.5],'LineWidth',2)
scatter(pca_score(idcs_candidates,1),-pca_score(idcs_candidates,2),'o','MarkerEdgeColor',[0.2,0.2,0.2],'LineWidth',2)

scatter(pca_score(idcs_candidates(these_odours(1)),1),-pca_score(idcs_candidates(these_odours(1)),2),100,'x','MarkerEdgeColor',[1,0,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(2)),1),-pca_score(idcs_candidates(these_odours(2)),2),100,'x','MarkerEdgeColor',[0.5,0,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(3)),1),-pca_score(idcs_candidates(these_odours(3)),2),100,'d','MarkerEdgeColor',[1,0,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(4)),1),-pca_score(idcs_candidates(these_odours(4)),2),100,'d','MarkerEdgeColor',[0.5,0,0],'LineWidth',2)

scatter(pca_score(idcs_candidates(these_odours(5)),1),-pca_score(idcs_candidates(these_odours(5)),2),100,'x','MarkerEdgeColor',[0,1,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(6)),1),-pca_score(idcs_candidates(these_odours(6)),2),100,'x','MarkerEdgeColor',[0,0.5,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(7)),1),-pca_score(idcs_candidates(these_odours(7)),2),100,'d','MarkerEdgeColor',[0,1,0],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(8)),1),-pca_score(idcs_candidates(these_odours(8)),2),100,'d','MarkerEdgeColor',[0,0.5,0],'LineWidth',2)

scatter(pca_score(idcs_candidates(these_odours(9)),1),-pca_score(idcs_candidates(these_odours(9)),2),100,'x','MarkerEdgeColor',[0,0,1],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(10)),1),-pca_score(idcs_candidates(these_odours(10)),2),100,'x','MarkerEdgeColor',[0,0,0.5],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(11)),1),-pca_score(idcs_candidates(these_odours(11)),2),100,'d','MarkerEdgeColor',[0,0,1],'LineWidth',2)
scatter(pca_score(idcs_candidates(these_odours(12)),1),-pca_score(idcs_candidates(these_odours(12)),2),100,'d','MarkerEdgeColor',[0,0,0.5],'LineWidth',2)

xlabel('1st PC (33% of variance)')
ylabel('2nd PC (10% of variance)')
annotation('textbox',[0.25 0 0.3 0.3],'String',[num2str(numGlobalOdours),' odourants',newline,num2str(numDescriptors),' descriptors'],'FitBoxToText','on');
title(['Global odour space, selected permutation: ',num2str(selectedPermutation)])


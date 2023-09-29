function resp = applyResponderCriteria_firstCluster(trg,resp,s2p_meta,amplitudeCriterion,extendedTargets)

if extendedTargets
    if isnan(amplitudeCriterion)
        fieldname = char('responders_ext');
    else
        fieldname = char(['responders_',num2str(amplitudeCriterion),'z_ext']);
        fieldname(find(fieldname=='.'))='d';
    end
else
    if isnan(amplitudeCriterion)
        fieldname = char('responders');
    else
        fieldname = char(['responders_',num2str(amplitudeCriterion),'z']);
        fieldname(find(fieldname=='.'))='d';
    end
end
fieldname_neg = char([fieldname,'_neg']);


%% Identify responders

% find all responsive cells
resp.firstCluster.(fieldname).respAll = nan(size(resp.firstCluster.stats_sig_pos));
for i=1:size(resp.firstCluster.stats_sig_pos,2)
    for j=1:size(resp.firstCluster.stats_sig_pos,1)
        if isnan(amplitudeCriterion)
            resp.firstCluster.(fieldname).respAll(j,i) = resp.firstCluster.stats_sig_pos(j,i);
        elseif ~isnan(resp.firstCluster.stats_sig_pos(j,i))
            resp.firstCluster.(fieldname).respAll(j,i) = resp.firstCluster.stats_sig_pos(j,i) && (resp.firstCluster.avgAct_net(j,i)>amplitudeCriterion);
        end
    end
end

% find all responsive cells (for negative ones)
resp.firstCluster.(fieldname_neg).respAll = nan(size(resp.firstCluster.stats_sig_neg));
for i=1:size(resp.firstCluster.stats_sig_neg,2)
    for j=1:size(resp.firstCluster.stats_sig_neg,1)
        if isnan(amplitudeCriterion)
            resp.firstCluster.(fieldname_neg).respAll(j,i) = resp.firstCluster.stats_sig_neg(j,i);
        elseif ~isnan(resp.firstCluster.stats_sig_pos(j,i))
            resp.firstCluster.(fieldname_neg).respAll(j,i) = resp.firstCluster.stats_sig_neg(j,i) && (resp.firstCluster.avgAct_net(j,i)<amplitudeCriterion);
        end
    end
end


%% Re-identify targets that could not be matched so far

if extendedTargets
    resp.firstCluster.(fieldname).targetedCells = trg.idcs_targetedCells;
    resp.firstCluster.(fieldname).targeted = trg.targeted;
    for j=1:size(resp.firstCluster.(fieldname).targetedCells,2)
        for i=1:size(resp.firstCluster.(fieldname).targetedCells,1)
            if isnan(trg.idcs_targetedCells(i,j))
                % can happen that the same cell will be assigned as target
                % for multiple clusters (otherwise, do sum(trg.targeted,2)==0)
                
                % get responding iscells that overlap with the laser beamlet
                candidateCellList = [];
                temp = find(resp.iscell_used==1 & resp.firstCluster.(fieldname).respAll(:,j)==1);
                for k=1:length(temp)
                    if ismember(round(trg.laser_y(i,j)),s2p_meta.stat{temp(k)}.xpix) && ismember(round(trg.laser_x(i,j)),s2p_meta.stat{temp(k)}.ypix)
                        candidateCellList = [candidateCellList, temp(k)];
                    end
                end
                
                % if 1 candidate, take it
                if length(candidateCellList)==1
                    resp.firstCluster.(fieldname).targetedCells(i,j) = candidateCellList;
                    resp.firstCluster.(fieldname).targeted(resp.firstCluster.(fieldname).targetedCells(i,j),j) = true;
                    
                % if more than 1 candidate, take the one with closer median    
                elseif length(candidateCellList)>1
                    [~,temp] = nanmin(trg.dist_closestLaser(candidateCellList,j));
                    resp.firstCluster.(fieldname).targetedCells(i,j) = candidateCellList(temp);
                    resp.firstCluster.(fieldname).targeted(resp.firstCluster.(fieldname).targetedCells(i,j),j) = true;
               
                % if no candidate,then try finding one with centroid distance of max 20um from laser target  
                else
                    candidateCellList = [];
                    candidateCellDist = [];
                    temp = find(resp.iscell_used==1 & resp.firstCluster.(fieldname).respAll(:,j)==1 & trg.targeted(:,j)==0);
                    for k=1:length(temp)
                        temp0 = pdist2(s2p_meta.stat{temp(k)}.med,[trg.laser_x(i,j),trg.laser_y(i,j)],'euclidean');
                        if temp0 <= resp.p.maxTargetCentroidDistance_pix
                            candidateCellList = [candidateCellList, temp(k)];
                            candidateCellDist = [candidateCellDist, temp0];
                        end
                    end
                    if ~isempty(candidateCellList)
                        [~,temp] = nanmin(candidateCellDist);
                        resp.firstCluster.(fieldname).targetedCells(i,j) = candidateCellList(temp);
                        resp.firstCluster.(fieldname).targeted(resp.firstCluster.(fieldname).targetedCells(i,j),j) = true;
                    end
                end
            end
        end
    end
else
    resp.firstCluster.(fieldname).targetedCells = trg.idcs_targetedCells;
    resp.firstCluster.(fieldname).targeted = trg.targeted;
end


%% Identify responders (continued)

% find responsive targets
resp.firstCluster.(fieldname).respTargeted = nan(size(resp.firstCluster.(fieldname).targetedCells));
for i=1:size(resp.firstCluster.(fieldname).targetedCells,2)
    for j=1:size(resp.firstCluster.(fieldname).targetedCells,1)
        temp = resp.firstCluster.(fieldname).targetedCells(j,i);
        if ~isnan(temp)
            if isnan(amplitudeCriterion)
                resp.firstCluster.(fieldname).respTargeted(j,i) = resp.firstCluster.stats_sig_pos(temp,i);
            elseif ~isnan(resp.firstCluster.stats_sig_pos(temp,i))
                resp.firstCluster.(fieldname).respTargeted(j,i) = resp.firstCluster.stats_sig_pos(temp,i) && (resp.firstCluster.avgAct_net(temp,i)>amplitudeCriterion);
            end
        end
    end
end

% find responsive non-targets
resp.firstCluster.(fieldname).respNonTargeted = nan(size(resp.firstCluster.(fieldname).targeted));
for i=1:size(resp.firstCluster.(fieldname).targeted,2)
    for j=1:size(resp.firstCluster.(fieldname).targeted,1)
        if ~resp.firstCluster.(fieldname).targeted(j,i)
            if isnan(amplitudeCriterion)
                resp.firstCluster.(fieldname).respNonTargeted(j,i) = resp.firstCluster.stats_sig_pos(j,i);
            elseif ~isnan(resp.firstCluster.stats_sig_pos(j,i))
                resp.firstCluster.(fieldname).respNonTargeted(j,i) = resp.firstCluster.stats_sig_pos(j,i) && (resp.firstCluster.avgAct_net(j,i)>amplitudeCriterion);
            end
        end
    end
end


%% Get responder numbers

temp = resp.firstCluster.(fieldname).respAll(:); 
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numAll = length(temp);
resp.firstCluster.(fieldname).numRespAll = sum(temp);
resp.firstCluster.(fieldname).proportion_RespAll_All = resp.firstCluster.(fieldname).numRespAll / resp.firstCluster.(fieldname).numAll;

temp = resp.firstCluster.(fieldname).respTargeted(:);
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numAllTargeted = size(resp.firstCluster.(fieldname).respTargeted,1)*size(resp.firstCluster.(fieldname).respTargeted,2);
resp.firstCluster.(fieldname).numIdentifiedTargeted = length(temp);
resp.firstCluster.(fieldname).numRespTargeted = sum(temp);
resp.firstCluster.(fieldname).proportion_RespTargeted_IdentifiedTargeted = resp.firstCluster.(fieldname).numRespTargeted / resp.firstCluster.(fieldname).numIdentifiedTargeted;
resp.firstCluster.(fieldname).proportion_RespTargeted_AllTargeted = resp.firstCluster.(fieldname).numRespTargeted / resp.firstCluster.(fieldname).numAllTargeted;
resp.firstCluster.(fieldname).proportion_RespTargeted_All = resp.firstCluster.(fieldname).numRespTargeted / resp.firstCluster.(fieldname).numAll;

temp = resp.firstCluster.(fieldname).respNonTargeted(:); 
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numAllNonTargeted = length(temp);
resp.firstCluster.(fieldname).numRespNonTargeted = sum(temp);
resp.firstCluster.(fieldname).proportion_RespNonTargeted_AllNonTargeted = resp.firstCluster.(fieldname).numRespNonTargeted / resp.firstCluster.(fieldname).numAllNonTargeted;
resp.firstCluster.(fieldname).proportion_RespNonTargeted_All = resp.firstCluster.(fieldname).numRespNonTargeted / resp.firstCluster.(fieldname).numAll;

resp.firstCluster.(fieldname).proportion_RespTargeted_RespNonTargeted = resp.firstCluster.(fieldname).numRespTargeted / resp.firstCluster.(fieldname).numRespNonTargeted;
resp.firstCluster.(fieldname).proportion_RespNonTargeted_RespTargeted = resp.firstCluster.(fieldname).numRespNonTargeted / resp.firstCluster.(fieldname).numRespTargeted;
resp.firstCluster.(fieldname).proportion_RespTargeted_RespAll = resp.firstCluster.(fieldname).numRespTargeted / resp.firstCluster.(fieldname).numRespAll;
resp.firstCluster.(fieldname).proportion_RespNonTargeted_RespAll = resp.firstCluster.(fieldname).numRespNonTargeted / resp.firstCluster.(fieldname).numRespAll;


%% Photostimulation group specificity

% responsive targeted A cells in A stim
temp = resp.firstCluster.(fieldname).respTargeted(:,trg.grouping(:,1));
temp = temp(:);
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numAllTargeted_A = size(resp.firstCluster.(fieldname).respTargeted,1)*size(resp.firstCluster.(fieldname).respTargeted,2)/2;
resp.firstCluster.(fieldname).numIdentifiedTargeted_A = length(temp);
resp.firstCluster.(fieldname).numRespTargeted_A = sum(temp);
resp.firstCluster.(fieldname).proportion_RespTargeted_IdentifiedTargeted_A = resp.firstCluster.(fieldname).numRespTargeted_A / resp.firstCluster.(fieldname).numIdentifiedTargeted_A;
resp.firstCluster.(fieldname).proportion_RespTargeted_AllTargeted_A = resp.firstCluster.(fieldname).numRespTargeted_A / resp.firstCluster.(fieldname).numAllTargeted_A;

% responsive targeted X cells in X stim
temp = resp.firstCluster.(fieldname).respTargeted(:,trg.grouping(:,2));
temp = temp(:);
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numAllTargeted_X = size(resp.firstCluster.(fieldname).respTargeted,1)*size(resp.firstCluster.(fieldname).respTargeted,2)/2;
resp.firstCluster.(fieldname).numIdentifiedTargeted_X = length(temp);
resp.firstCluster.(fieldname).numRespTargeted_X = sum(temp);
resp.firstCluster.(fieldname).proportion_RespTargeted_IdentifiedTargeted_X = resp.firstCluster.(fieldname).numRespTargeted_X / resp.firstCluster.(fieldname).numIdentifiedTargeted_X;
resp.firstCluster.(fieldname).proportion_RespTargeted_AllTargeted_X = resp.firstCluster.(fieldname).numRespTargeted_X / resp.firstCluster.(fieldname).numAllTargeted_X;

% (all) responsive stim A cells in A stim
temp2 = resp.firstCluster.(fieldname).targetedCells(:,trg.grouping(:,1));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.firstCluster.(fieldname).respAll(temp2,trg.grouping(:,1)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numRespIpsi_A = sum(temp);

% (all) responsive stim X cells in X stim
temp2 = resp.firstCluster.(fieldname).targetedCells(:,trg.grouping(:,2));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.firstCluster.(fieldname).respAll(temp2,trg.grouping(:,2)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numRespIpsi_X = sum(temp);
resp.firstCluster.(fieldname).numRespIpsi = resp.firstCluster.(fieldname).numRespIpsi_A + resp.firstCluster.(fieldname).numRespIpsi_X;

% responsive stim X cells in A stim
temp2 = resp.firstCluster.(fieldname).targetedCells(:,trg.grouping(:,2));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.firstCluster.(fieldname).respAll(temp2,trg.grouping(:,1)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numRespContra_A = sum(temp);

% responsive stim A cells in X stim
temp2 = resp.firstCluster.(fieldname).targetedCells(:,trg.grouping(:,1));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.firstCluster.(fieldname).respAll(temp2,trg.grouping(:,2)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numRespContra_X = sum(temp);
resp.firstCluster.(fieldname).numRespContra = resp.firstCluster.(fieldname).numRespContra_A + resp.firstCluster.(fieldname).numRespContra_X;

% responsive no-group cells in A stim
temp2 = 1:size(resp.firstCluster.(fieldname).respAll,1);
temp2 = setdiff(temp2,sort(rmmissing(resp.firstCluster.(fieldname).targetedCells(:))));
temp = resp.firstCluster.(fieldname).respAll(temp2,trg.grouping(:,1)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numRespNoGroupCells_A = sum(temp);

% responsive no-group cells in X stim
temp2 = 1:size(resp.firstCluster.(fieldname).respAll,1);
temp2 = setdiff(temp2,sort(rmmissing(resp.firstCluster.(fieldname).targetedCells(:))));
temp = resp.firstCluster.(fieldname).respAll(temp2,trg.grouping(:,2)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.firstCluster.(fieldname).numRespNoGroupCells_X = sum(temp);
resp.firstCluster.(fieldname).numRespNoGroupCells = resp.firstCluster.(fieldname).numRespNoGroupCells_A + resp.firstCluster.(fieldname).numRespNoGroupCells_X;

% specificity index
resp.firstCluster.(fieldname).specificity_respTargeted = resp.firstCluster.(fieldname).numRespTargeted / resp.firstCluster.(fieldname).numRespAll;
resp.firstCluster.(fieldname).specificity_respIpsiGroup = resp.firstCluster.(fieldname).numRespIpsi / resp.firstCluster.(fieldname).numRespAll;
resp.firstCluster.(fieldname).specificity_respContraGroup = resp.firstCluster.(fieldname).numRespContra / resp.firstCluster.(fieldname).numRespAll;
resp.firstCluster.(fieldname).specificity_respNoGroupCells = resp.firstCluster.(fieldname).numRespNoGroupCells / resp.firstCluster.(fieldname).numRespAll;

resp.firstCluster.(fieldname) = orderfields(resp.firstCluster.(fieldname));
end
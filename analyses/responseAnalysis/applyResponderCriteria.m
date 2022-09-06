function resp = applyResponderCriteria(trg,resp,s2p_meta,amplitudeCriterion,extendedTargets)

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


%% Identify responders

% find all responsive cells
resp.(fieldname).respAll = nan(size(resp.stats_sig_pos));
for i=1:size(resp.stats_sig_pos,2)
    for j=1:size(resp.stats_sig_pos,1)
        if isnan(amplitudeCriterion)
            resp.(fieldname).respAll(j,i) = resp.stats_sig_pos(j,i);
        elseif ~isnan(resp.stats_sig_pos(j,i))
            resp.(fieldname).respAll(j,i) = resp.stats_sig_pos(j,i) && (resp.avgAct_net(j,i)>amplitudeCriterion);
        end
    end
end


%% Re-identify targets that could not be matched so far

if extendedTargets
    resp.(fieldname).targetedCells = trg.idcs_targetedCells;
    resp.(fieldname).targeted = trg.targeted;
    for j=1:size(resp.(fieldname).targetedCells,2)
        for i=1:size(resp.(fieldname).targetedCells,1)
            if isnan(trg.idcs_targetedCells(i,j))
                % can happen that the same cell will be assigned as target
                % for multiple clusters (otherwise, do sum(trg.targeted,2)==0)
                
                % get responding iscells that overlap with the laser beamlet
                candidateCellList = [];
                temp = find(resp.iscell_used==1 & resp.(fieldname).respAll(:,j)==1);
                for k=1:length(temp)
                    if ismember(round(trg.laser_y(i,j)),s2p_meta.stat{temp(k)}.xpix) && ismember(round(trg.laser_x(i,j)),s2p_meta.stat{temp(k)}.ypix)
                        candidateCellList = [candidateCellList, temp(k)];
                    end
                end
                
                % if 1 candidate, take it
                if length(candidateCellList)==1
                    resp.(fieldname).targetedCells(i,j) = candidateCellList;
                    resp.(fieldname).targeted(resp.(fieldname).targetedCells(i,j),j) = true;
                    
                % if more than 1 candidate, take the one with closer median    
                elseif length(candidateCellList)>1
                    [~,temp] = nanmin(trg.dist_closestLaser(candidateCellList,j));
                    resp.(fieldname).targetedCells(i,j) = candidateCellList(temp);
                    resp.(fieldname).targeted(resp.(fieldname).targetedCells(i,j),j) = true;
               
                % if no candidate,then try finding one with centroid distance of max 20um from laser target  
                else
                    candidateCellList = [];
                    candidateCellDist = [];
                    temp = find(resp.iscell_used==1 & resp.(fieldname).respAll(:,j)==1 & trg.targeted(:,j)==0);
                    for k=1:length(temp)
                        temp0 = pdist2(s2p_meta.stat{temp(k)}.med,[trg.laser_x(i,j),trg.laser_y(i,j)],'euclidean');
                        if temp0 <= resp.p.maxTargetCentroidDistance_pix
                            candidateCellList = [candidateCellList, temp(k)];
                            candidateCellDist = [candidateCellDist, temp0];
                        end
                    end
                    if ~isempty(candidateCellList)
                        [~,temp] = nanmin(candidateCellDist);
                        resp.(fieldname).targetedCells(i,j) = candidateCellList(temp);
                        resp.(fieldname).targeted(resp.(fieldname).targetedCells(i,j),j) = true;
                    end
                end
            end
        end
    end
else
    resp.(fieldname).targetedCells = trg.idcs_targetedCells;
    resp.(fieldname).targeted = trg.targeted;
end


%% Identify responders (continued)

% find responsive targets
resp.(fieldname).respTargeted = nan(size(resp.(fieldname).targetedCells));
for i=1:size(resp.(fieldname).targetedCells,2)
    for j=1:size(resp.(fieldname).targetedCells,1)
        temp = resp.(fieldname).targetedCells(j,i);
        if ~isnan(temp)
            if isnan(amplitudeCriterion)
                resp.(fieldname).respTargeted(j,i) = resp.stats_sig_pos(temp,i);
            elseif ~isnan(resp.stats_sig_pos(temp,i))
                resp.(fieldname).respTargeted(j,i) = resp.stats_sig_pos(temp,i) && (resp.avgAct_net(temp,i)>amplitudeCriterion);
            end
        end
    end
end

% find responsive non-targets
resp.(fieldname).respNonTargeted = nan(size(resp.(fieldname).targeted));
for i=1:size(resp.(fieldname).targeted,2)
    for j=1:size(resp.(fieldname).targeted,1)
        if ~resp.(fieldname).targeted(j,i)
            if isnan(amplitudeCriterion)
                resp.(fieldname).respNonTargeted(j,i) = resp.stats_sig_pos(j,i);
            elseif ~isnan(resp.stats_sig_pos(j,i))
                resp.(fieldname).respNonTargeted(j,i) = resp.stats_sig_pos(j,i) && (resp.avgAct_net(j,i)>amplitudeCriterion);
            end
        end
    end
end


%% Get responder numbers

temp = resp.(fieldname).respAll(:); 
temp = temp(~isnan(temp));
resp.(fieldname).numAll = length(temp);
resp.(fieldname).numRespAll = sum(temp);
resp.(fieldname).proportion_RespAll_All = resp.(fieldname).numRespAll / resp.(fieldname).numAll;

temp = resp.(fieldname).respTargeted(:);
temp = temp(~isnan(temp));
resp.(fieldname).numAllTargeted = size(resp.(fieldname).respTargeted,1)*size(resp.(fieldname).respTargeted,2);
resp.(fieldname).numIdentifiedTargeted = length(temp);
resp.(fieldname).numRespTargeted = sum(temp);
resp.(fieldname).proportion_RespTargeted_IdentifiedTargeted = resp.(fieldname).numRespTargeted / resp.(fieldname).numIdentifiedTargeted;
resp.(fieldname).proportion_RespTargeted_AllTargeted = resp.(fieldname).numRespTargeted / resp.(fieldname).numAllTargeted;
resp.(fieldname).proportion_RespTargeted_All = resp.(fieldname).numRespTargeted / resp.(fieldname).numAll;

temp = resp.(fieldname).respNonTargeted(:); 
temp = temp(~isnan(temp));
resp.(fieldname).numAllNonTargeted = length(temp);
resp.(fieldname).numRespNonTargeted = sum(temp);
resp.(fieldname).proportion_RespNonTargeted_AllNonTargeted = resp.(fieldname).numRespNonTargeted / resp.(fieldname).numAllNonTargeted;
resp.(fieldname).proportion_RespNonTargeted_All = resp.(fieldname).numRespNonTargeted / resp.(fieldname).numAll;

resp.(fieldname).proportion_RespTargeted_RespNonTargeted = resp.(fieldname).numRespTargeted / resp.(fieldname).numRespNonTargeted;
resp.(fieldname).proportion_RespNonTargeted_RespTargeted = resp.(fieldname).numRespNonTargeted / resp.(fieldname).numRespTargeted;
resp.(fieldname).proportion_RespTargeted_RespAll = resp.(fieldname).numRespTargeted / resp.(fieldname).numRespAll;
resp.(fieldname).proportion_RespNonTargeted_RespAll = resp.(fieldname).numRespNonTargeted / resp.(fieldname).numRespAll;


%% Photostimulation group specificity

% responsive targeted A cells in A stim
temp = resp.(fieldname).respTargeted(:,trg.grouping(:,1));
temp = temp(:);
temp = temp(~isnan(temp));
resp.(fieldname).numAllTargeted_A = size(resp.(fieldname).respTargeted,1)*size(resp.(fieldname).respTargeted,2)/2;
resp.(fieldname).numIdentifiedTargeted_A = length(temp);
resp.(fieldname).numRespTargeted_A = sum(temp);
resp.(fieldname).proportion_RespTargeted_IdentifiedTargeted_A = resp.(fieldname).numRespTargeted_A / resp.(fieldname).numIdentifiedTargeted_A;
resp.(fieldname).proportion_RespTargeted_AllTargeted_A = resp.(fieldname).numRespTargeted_A / resp.(fieldname).numAllTargeted_A;

% responsive targeted X cells in X stim
temp = resp.(fieldname).respTargeted(:,trg.grouping(:,2));
temp = temp(:);
temp = temp(~isnan(temp));
resp.(fieldname).numAllTargeted_X = size(resp.(fieldname).respTargeted,1)*size(resp.(fieldname).respTargeted,2)/2;
resp.(fieldname).numIdentifiedTargeted_X = length(temp);
resp.(fieldname).numRespTargeted_X = sum(temp);
resp.(fieldname).proportion_RespTargeted_IdentifiedTargeted_X = resp.(fieldname).numRespTargeted_X / resp.(fieldname).numIdentifiedTargeted_X;
resp.(fieldname).proportion_RespTargeted_AllTargeted_X = resp.(fieldname).numRespTargeted_X / resp.(fieldname).numAllTargeted_X;

% (all) responsive stim A cells in A stim
temp2 = resp.(fieldname).targetedCells(:,trg.grouping(:,1));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.(fieldname).respAll(temp2,trg.grouping(:,1)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.(fieldname).numRespIpsi_A = sum(temp);

% (all) responsive stim X cells in X stim
temp2 = resp.(fieldname).targetedCells(:,trg.grouping(:,2));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.(fieldname).respAll(temp2,trg.grouping(:,2)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.(fieldname).numRespIpsi_X = sum(temp);
resp.(fieldname).numRespIpsi = resp.(fieldname).numRespIpsi_A + resp.(fieldname).numRespIpsi_X;

% responsive stim X cells in A stim
temp2 = resp.(fieldname).targetedCells(:,trg.grouping(:,2));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.(fieldname).respAll(temp2,trg.grouping(:,1)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.(fieldname).numRespContra_A = sum(temp);

% responsive stim A cells in X stim
temp2 = resp.(fieldname).targetedCells(:,trg.grouping(:,1));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.(fieldname).respAll(temp2,trg.grouping(:,2)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.(fieldname).numRespContra_X = sum(temp);
resp.(fieldname).numRespContra = resp.(fieldname).numRespContra_A + resp.(fieldname).numRespContra_X;

% responsive no-group cells in A stim
temp2 = 1:size(resp.(fieldname).respAll,1);
temp2 = setdiff(temp2,sort(rmmissing(resp.(fieldname).targetedCells(:))));
temp = resp.(fieldname).respAll(temp2,trg.grouping(:,1)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.(fieldname).numRespNoGroupCells_A = sum(temp);

% responsive no-group cells in X stim
temp2 = 1:size(resp.(fieldname).respAll,1);
temp2 = setdiff(temp2,sort(rmmissing(resp.(fieldname).targetedCells(:))));
temp = resp.(fieldname).respAll(temp2,trg.grouping(:,2)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.(fieldname).numRespNoGroupCells_X = sum(temp);
resp.(fieldname).numRespNoGroupCells = resp.(fieldname).numRespNoGroupCells_A + resp.(fieldname).numRespNoGroupCells_X;

% specificity index
resp.(fieldname).specificity_respTargeted = resp.(fieldname).numRespTargeted / resp.(fieldname).numRespAll;
resp.(fieldname).specificity_respIpsiGroup = resp.(fieldname).numRespIpsi / resp.(fieldname).numRespAll;
resp.(fieldname).specificity_respContraGroup = resp.(fieldname).numRespContra / resp.(fieldname).numRespAll;
resp.(fieldname).specificity_respNoGroupCells = resp.(fieldname).numRespNoGroupCells / resp.(fieldname).numRespAll;

resp.(fieldname) = orderfields(resp.(fieldname));
end
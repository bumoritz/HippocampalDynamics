function resp = applyResponderCriteria_noStim(trg,resp,s2p_meta,amplitudeCriterion,extendedTargets)

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
resp.noStim.(fieldname).respAll = nan(size(resp.noStim.stats_sig_pos));
for i=1:size(resp.noStim.stats_sig_pos,2)
    for j=1:size(resp.noStim.stats_sig_pos,1)
        if isnan(amplitudeCriterion)
            resp.noStim.(fieldname).respAll(j,i) = resp.noStim.stats_sig_pos(j,i);
        elseif ~isnan(resp.noStim.stats_sig_pos(j,i))
            resp.noStim.(fieldname).respAll(j,i) = resp.noStim.stats_sig_pos(j,i) && (resp.noStim.avgAct_net(j,i)>amplitudeCriterion);
        end
    end
end

% find all responsive cells (for negative ones)
resp.noStim.(fieldname_neg).respAll = nan(size(resp.noStim.stats_sig_neg));
for i=1:size(resp.noStim.stats_sig_neg,2)
    for j=1:size(resp.noStim.stats_sig_neg,1)
        if isnan(amplitudeCriterion)
            resp.noStim.(fieldname_neg).respAll(j,i) = resp.noStim.stats_sig_neg(j,i);
        elseif ~isnan(resp.noStim.stats_sig_pos(j,i))
            resp.noStim.(fieldname_neg).respAll(j,i) = resp.noStim.stats_sig_neg(j,i) && (resp.noStim.avgAct_net(j,i)<amplitudeCriterion);
        end
    end
end


%% Re-identify targets that could not be matched so far

if extendedTargets
    resp.noStim.(fieldname).targetedCells = trg.idcs_targetedCells;
    resp.noStim.(fieldname).targeted = trg.targeted;
    for j=1:size(resp.noStim.(fieldname).targetedCells,2)
        for i=1:size(resp.noStim.(fieldname).targetedCells,1)
            if isnan(trg.idcs_targetedCells(i,j))
                % can happen that the same cell will be assigned as target
                % for multiple clusters (otherwise, do sum(trg.targeted,2)==0)
                
                % get responding iscells that overlap with the laser beamlet
                candidateCellList = [];
                temp = find(resp.iscell_used==1 & resp.noStim.(fieldname).respAll(:,j)==1);
                for k=1:length(temp)
                    if ismember(round(trg.laser_y(i,j)),s2p_meta.stat{temp(k)}.xpix) && ismember(round(trg.laser_x(i,j)),s2p_meta.stat{temp(k)}.ypix)
                        candidateCellList = [candidateCellList, temp(k)];
                    end
                end
                
                % if 1 candidate, take it
                if length(candidateCellList)==1
                    resp.noStim.(fieldname).targetedCells(i,j) = candidateCellList;
                    resp.noStim.(fieldname).targeted(resp.noStim.(fieldname).targetedCells(i,j),j) = true;
                    
                % if more than 1 candidate, take the one with closer median    
                elseif length(candidateCellList)>1
                    [~,temp] = nanmin(trg.dist_closestLaser(candidateCellList,j));
                    resp.noStim.(fieldname).targetedCells(i,j) = candidateCellList(temp);
                    resp.noStim.(fieldname).targeted(resp.noStim.(fieldname).targetedCells(i,j),j) = true;
               
                % if no candidate,then try finding one with centroid distance of max 20um from laser target  
                else
                    candidateCellList = [];
                    candidateCellDist = [];
                    temp = find(resp.iscell_used==1 & resp.noStim.(fieldname).respAll(:,j)==1 & trg.targeted(:,j)==0);
                    for k=1:length(temp)
                        temp0 = pdist2(s2p_meta.stat{temp(k)}.med,[trg.laser_x(i,j),trg.laser_y(i,j)],'euclidean');
                        if temp0 <= resp.p.maxTargetCentroidDistance_pix
                            candidateCellList = [candidateCellList, temp(k)];
                            candidateCellDist = [candidateCellDist, temp0];
                        end
                    end
                    if ~isempty(candidateCellList)
                        [~,temp] = nanmin(candidateCellDist);
                        resp.noStim.(fieldname).targetedCells(i,j) = candidateCellList(temp);
                        resp.noStim.(fieldname).targeted(resp.noStim.(fieldname).targetedCells(i,j),j) = true;
                    end
                end
            end
        end
    end
else
    resp.noStim.(fieldname).targetedCells = trg.idcs_targetedCells;
    resp.noStim.(fieldname).targeted = trg.targeted;
end


%% Identify responders (continued)

% find responsive targets
resp.noStim.(fieldname).respTargeted = nan(size(resp.noStim.(fieldname).targetedCells));
for i=1:size(resp.noStim.(fieldname).targetedCells,2)
    for j=1:size(resp.noStim.(fieldname).targetedCells,1)
        temp = resp.noStim.(fieldname).targetedCells(j,i);
        if ~isnan(temp)
            if isnan(amplitudeCriterion)
                resp.noStim.(fieldname).respTargeted(j,i) = resp.noStim.stats_sig_pos(temp,i);
            elseif ~isnan(resp.noStim.stats_sig_pos(temp,i))
                resp.noStim.(fieldname).respTargeted(j,i) = resp.noStim.stats_sig_pos(temp,i) && (resp.noStim.avgAct_net(temp,i)>amplitudeCriterion);
            end
        end
    end
end

% find responsive non-targets
resp.noStim.(fieldname).respNonTargeted = nan(size(resp.noStim.(fieldname).targeted));
for i=1:size(resp.noStim.(fieldname).targeted,2)
    for j=1:size(resp.noStim.(fieldname).targeted,1)
        if ~resp.noStim.(fieldname).targeted(j,i)
            if isnan(amplitudeCriterion)
                resp.noStim.(fieldname).respNonTargeted(j,i) = resp.noStim.stats_sig_pos(j,i);
            elseif ~isnan(resp.noStim.stats_sig_pos(j,i))
                resp.noStim.(fieldname).respNonTargeted(j,i) = resp.noStim.stats_sig_pos(j,i) && (resp.noStim.avgAct_net(j,i)>amplitudeCriterion);
            end
        end
    end
end


%% Get responder numbers

temp = resp.noStim.(fieldname).respAll(:); 
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numAll = length(temp);
resp.noStim.(fieldname).numRespAll = sum(temp);
resp.noStim.(fieldname).proportion_RespAll_All = resp.noStim.(fieldname).numRespAll / resp.noStim.(fieldname).numAll;

temp = resp.noStim.(fieldname).respTargeted(:);
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numAllTargeted = size(resp.noStim.(fieldname).respTargeted,1)*size(resp.noStim.(fieldname).respTargeted,2);
resp.noStim.(fieldname).numIdentifiedTargeted = length(temp);
resp.noStim.(fieldname).numRespTargeted = sum(temp);
resp.noStim.(fieldname).proportion_RespTargeted_IdentifiedTargeted = resp.noStim.(fieldname).numRespTargeted / resp.noStim.(fieldname).numIdentifiedTargeted;
resp.noStim.(fieldname).proportion_RespTargeted_AllTargeted = resp.noStim.(fieldname).numRespTargeted / resp.noStim.(fieldname).numAllTargeted;
resp.noStim.(fieldname).proportion_RespTargeted_All = resp.noStim.(fieldname).numRespTargeted / resp.noStim.(fieldname).numAll;

temp = resp.noStim.(fieldname).respNonTargeted(:); 
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numAllNonTargeted = length(temp);
resp.noStim.(fieldname).numRespNonTargeted = sum(temp);
resp.noStim.(fieldname).proportion_RespNonTargeted_AllNonTargeted = resp.noStim.(fieldname).numRespNonTargeted / resp.noStim.(fieldname).numAllNonTargeted;
resp.noStim.(fieldname).proportion_RespNonTargeted_All = resp.noStim.(fieldname).numRespNonTargeted / resp.noStim.(fieldname).numAll;

resp.noStim.(fieldname).proportion_RespTargeted_RespNonTargeted = resp.noStim.(fieldname).numRespTargeted / resp.noStim.(fieldname).numRespNonTargeted;
resp.noStim.(fieldname).proportion_RespNonTargeted_RespTargeted = resp.noStim.(fieldname).numRespNonTargeted / resp.noStim.(fieldname).numRespTargeted;
resp.noStim.(fieldname).proportion_RespTargeted_RespAll = resp.noStim.(fieldname).numRespTargeted / resp.noStim.(fieldname).numRespAll;
resp.noStim.(fieldname).proportion_RespNonTargeted_RespAll = resp.noStim.(fieldname).numRespNonTargeted / resp.noStim.(fieldname).numRespAll;


%% Photostimulation group specificity

% responsive targeted A cells in A stim
temp = resp.noStim.(fieldname).respTargeted(:,trg.grouping(:,1));
temp = temp(:);
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numAllTargeted_A = size(resp.noStim.(fieldname).respTargeted,1)*size(resp.noStim.(fieldname).respTargeted,2)/2;
resp.noStim.(fieldname).numIdentifiedTargeted_A = length(temp);
resp.noStim.(fieldname).numRespTargeted_A = sum(temp);
resp.noStim.(fieldname).proportion_RespTargeted_IdentifiedTargeted_A = resp.noStim.(fieldname).numRespTargeted_A / resp.noStim.(fieldname).numIdentifiedTargeted_A;
resp.noStim.(fieldname).proportion_RespTargeted_AllTargeted_A = resp.noStim.(fieldname).numRespTargeted_A / resp.noStim.(fieldname).numAllTargeted_A;

% responsive targeted X cells in X stim
temp = resp.noStim.(fieldname).respTargeted(:,trg.grouping(:,2));
temp = temp(:);
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numAllTargeted_X = size(resp.noStim.(fieldname).respTargeted,1)*size(resp.noStim.(fieldname).respTargeted,2)/2;
resp.noStim.(fieldname).numIdentifiedTargeted_X = length(temp);
resp.noStim.(fieldname).numRespTargeted_X = sum(temp);
resp.noStim.(fieldname).proportion_RespTargeted_IdentifiedTargeted_X = resp.noStim.(fieldname).numRespTargeted_X / resp.noStim.(fieldname).numIdentifiedTargeted_X;
resp.noStim.(fieldname).proportion_RespTargeted_AllTargeted_X = resp.noStim.(fieldname).numRespTargeted_X / resp.noStim.(fieldname).numAllTargeted_X;

% (all) responsive stim A cells in A stim
temp2 = resp.noStim.(fieldname).targetedCells(:,trg.grouping(:,1));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.noStim.(fieldname).respAll(temp2,trg.grouping(:,1)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numRespIpsi_A = sum(temp);

% (all) responsive stim X cells in X stim
temp2 = resp.noStim.(fieldname).targetedCells(:,trg.grouping(:,2));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.noStim.(fieldname).respAll(temp2,trg.grouping(:,2)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numRespIpsi_X = sum(temp);
resp.noStim.(fieldname).numRespIpsi = resp.noStim.(fieldname).numRespIpsi_A + resp.noStim.(fieldname).numRespIpsi_X;

% responsive stim X cells in A stim
temp2 = resp.noStim.(fieldname).targetedCells(:,trg.grouping(:,2));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.noStim.(fieldname).respAll(temp2,trg.grouping(:,1)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numRespContra_A = sum(temp);

% responsive stim A cells in X stim
temp2 = resp.noStim.(fieldname).targetedCells(:,trg.grouping(:,1));
temp2 = sort(rmmissing(temp2(:)));
temp = resp.noStim.(fieldname).respAll(temp2,trg.grouping(:,2)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numRespContra_X = sum(temp);
resp.noStim.(fieldname).numRespContra = resp.noStim.(fieldname).numRespContra_A + resp.noStim.(fieldname).numRespContra_X;

% responsive no-group cells in A stim
temp2 = 1:size(resp.noStim.(fieldname).respAll,1);
temp2 = setdiff(temp2,sort(rmmissing(resp.noStim.(fieldname).targetedCells(:))));
temp = resp.noStim.(fieldname).respAll(temp2,trg.grouping(:,1)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numRespNoGroupCells_A = sum(temp);

% responsive no-group cells in X stim
temp2 = 1:size(resp.noStim.(fieldname).respAll,1);
temp2 = setdiff(temp2,sort(rmmissing(resp.noStim.(fieldname).targetedCells(:))));
temp = resp.noStim.(fieldname).respAll(temp2,trg.grouping(:,2)); 
temp = temp(:);
temp = temp(~isnan(temp));
resp.noStim.(fieldname).numRespNoGroupCells_X = sum(temp);
resp.noStim.(fieldname).numRespNoGroupCells = resp.noStim.(fieldname).numRespNoGroupCells_A + resp.noStim.(fieldname).numRespNoGroupCells_X;

% specificity index
resp.noStim.(fieldname).specificity_respTargeted = resp.noStim.(fieldname).numRespTargeted / resp.noStim.(fieldname).numRespAll;
resp.noStim.(fieldname).specificity_respIpsiGroup = resp.noStim.(fieldname).numRespIpsi / resp.noStim.(fieldname).numRespAll;
resp.noStim.(fieldname).specificity_respContraGroup = resp.noStim.(fieldname).numRespContra / resp.noStim.(fieldname).numRespAll;
resp.noStim.(fieldname).specificity_respNoGroupCells = resp.noStim.(fieldname).numRespNoGroupCells / resp.noStim.(fieldname).numRespAll;

resp.noStim.(fieldname) = orderfields(resp.noStim.(fieldname));
end
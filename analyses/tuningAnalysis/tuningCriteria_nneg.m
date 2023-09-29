function [passed_nneg,passed_stats_nneg] = tuningCriteria_nneg(shuffling,firingField,timewindow,ipsi,contra,p,prop)
%shuffling = shuffling_all; firingField = firingField_all; timewindow = 'AW'; ipsi = 'A'; contra = 'X';

%% Preparations

shuffling_ipsi = shuffling.([ipsi,'_',timewindow]);
shuffling_contra = shuffling.([contra,'_',timewindow]);
firingField_ipsi = firingField.([ipsi,'_',timewindow]);
firingField_contra = firingField.([contra,'_',timewindow]);


%% Apply basic tuning critera

temp1 = ones(size(prop.iscell,1),1);
temp2 = ones(size(prop.iscell,1),1);

if p.tng.criteria_shuffling
    temp1 = floor((shuffling_ipsi.significant + temp1)/2);
    temp2 = floor((shuffling_contra.significant + temp2)/2);
end
if p.tng.criteria_peakInWindow
    temp1 = floor((firingField_ipsi.binary.peakInWindow + temp1)/2);
    temp2 = floor((firingField_contra.binary.peakInWindow + temp2)/2);
end
if p.tng.criteria_aboveBaseline
    temp1 = floor((firingField_ipsi.diffBl_pos + temp1)/2);
    temp2 = floor((firingField_contra.diffBl_pos + temp2)/2);
end
if p.tng.criteria_reliability
    temp1 = floor((firingField_ipsi.binary.reliability_ipsi + temp1)/2);
    temp2 = floor((firingField_contra.binary.reliability_ipsi + temp2)/2); % firingField_contra.binary.reliability_ipsi is indeed correct
end
if p.tng.criteria_amplitude
    temp1 = floor((firingField_ipsi.binary.meanAmplitude_blSub_ipsi + temp1)/2);
    temp2 = floor((firingField_contra.binary.meanAmplitude_blSub_ipsi + temp2)/2); % firingField_contra.binary.reliability_ipsi is indeed correct
end

% --- nneg section ---
nneg_ipsi = firingField.([ipsi,'_1s']);
nneg_contra = firingField.([contra,'_1s']);
temp1 = floor(((1-nneg_ipsi.negativeGoing) + temp1)/2);
temp2 = floor(((1-nneg_contra.negativeGoing) + temp2)/2);
% --------------------

passed_nneg.(ipsi) = temp1;
passed_nneg.(contra) = temp2;


%% Apply all tuning critera

% exclusive vs non-exclusive
passed_nneg.([ipsi,'only']) = floor((passed_nneg.(ipsi) + 1-passed_nneg.(contra))/2);
passed_nneg.([contra,'only']) = floor((1-passed_nneg.(ipsi) + passed_nneg.(contra))/2);
passed_nneg.([ipsi,'and',contra]) = floor((passed_nneg.(ipsi) + passed_nneg.(contra))/2);
passed_nneg.([ipsi,'or',contra]) = round((passed_nneg.(ipsi) + passed_nneg.(contra))/2);

% selectivity threshold
passed_nneg.([ipsi,'_sel']) = floor((passed_nneg.(ipsi) + double(firingField_ipsi.selectivity_ipsi>p.tng.selectivityThreshold))/2); 
passed_nneg.([contra,'_sel']) = floor((passed_nneg.(contra) + double(firingField_contra.selectivity_ipsi>p.tng.selectivityThreshold))/2); % firingField_contra.selectivity_ipsi is indeed correct
passed_nneg.([ipsi,'only_sel']) = floor((passed_nneg.([ipsi,'only']) + double(firingField_ipsi.selectivity_ipsi>p.tng.selectivityThreshold))/2); 
passed_nneg.([contra,'only_sel']) = floor((passed_nneg.([contra,'only']) + double(firingField_contra.selectivity_ipsi>p.tng.selectivityThreshold))/2); % firingField_contra.selectivity_ipsi is indeed correct

% early vs late
if strcmp(timewindow,'AW')
    passed_nneg.([ipsi,'_early']) = floor((passed_nneg.(ipsi) + firingField_ipsi.binary.early)/2);
    passed_nneg.([ipsi,'_late']) = floor((passed_nneg.(ipsi) + firingField_ipsi.binary.late)/2);
    passed_nneg.([contra,'_early']) = floor((passed_nneg.(contra) + firingField_contra.binary.early)/2);
    passed_nneg.([contra,'_late']) = floor((passed_nneg.(contra) + firingField_contra.binary.late)/2);
    passed_nneg.early = ceil((passed_nneg.([ipsi,'_early']) + passed_nneg.([contra,'_early']))/2);
    passed_nneg.late = ceil((passed_nneg.([ipsi,'_late']) + passed_nneg.([contra,'_late']))/2);
end


%% Summarise results

% basic
passed_stats_nneg.([ipsi,'_num']) = nansum(passed_nneg.(ipsi));
passed_stats_nneg.([contra,'_num']) = nansum(passed_nneg.(contra));
passed_stats_nneg.([ipsi,'_fractionOfCells']) = passed_stats_nneg.([ipsi,'_num'])/nansum(prop.iscell);
passed_stats_nneg.([contra,'_fractionOfCells']) = passed_stats_nneg.([contra,'_num'])/nansum(prop.iscell);

% exclusive vs non-exclusive
passed_stats_nneg.([ipsi,'or',contra,'_num']) = nansum(passed_nneg.([ipsi,'or',contra]));
passed_stats_nneg.([ipsi,'and',contra,'_num']) = nansum(passed_nneg.([ipsi,'and',contra]));
passed_stats_nneg.([ipsi,'or',contra,'_fractionOfCells']) = passed_stats_nneg.([ipsi,'or',contra,'_num'])/nansum(prop.iscell);
passed_stats_nneg.([ipsi,'and',contra,'_fractionOfCells']) = passed_stats_nneg.([ipsi,'and',contra,'_num'])/nansum(prop.iscell);
passed_stats_nneg.([ipsi,'and',contra,'_fractionOf',ipsi,'or',contra,'Cells']) = passed_stats_nneg.([ipsi,'and',contra,'_num'])/passed_stats_nneg.([ipsi,'or',contra,'_num']);

% only
passed_stats_nneg.([ipsi,'only_num']) = nansum(passed_nneg.([ipsi,'only']));
passed_stats_nneg.([contra,'only_num']) = nansum(passed_nneg.([contra,'only']));
passed_stats_nneg.([ipsi,'only_fractionOfCells']) = passed_stats_nneg.([ipsi,'only_num'])/nansum(prop.iscell);
passed_stats_nneg.([contra,'only_fractionOfCells']) = passed_stats_nneg.([contra,'only_num'])/nansum(prop.iscell);

% sel
passed_stats_nneg.([ipsi,'_sel_num']) = nansum(passed_nneg.([ipsi,'_sel']));
passed_stats_nneg.([contra,'_sel_num']) = nansum(passed_nneg.([contra,'_sel']));
passed_stats_nneg.([ipsi,'only_sel_num']) = nansum(passed_nneg.([ipsi,'only_sel']));
passed_stats_nneg.([contra,'only_sel_num']) = nansum(passed_nneg.([contra,'only_sel']));
passed_stats_nneg.([ipsi,'_sel_fractionOfCells']) = passed_stats_nneg.([ipsi,'_sel_num'])/nansum(prop.iscell);
passed_stats_nneg.([contra,'_sel_fractionOfCells']) = passed_stats_nneg.([contra,'_sel_num'])/nansum(prop.iscell);
passed_stats_nneg.([ipsi,'only_sel_fractionOfCells']) = passed_stats_nneg.([ipsi,'only_sel_num'])/nansum(prop.iscell);
passed_stats_nneg.([contra,'only_sel_fractionOfCells']) = passed_stats_nneg.([contra,'only_sel_num'])/nansum(prop.iscell);

% early vs late
if strcmp(timewindow,'AW')
    passed_stats_nneg.([ipsi,'_early_num']) = nansum(passed_nneg.([ipsi,'_early']));
    passed_stats_nneg.([ipsi,'_late_num']) = nansum(passed_nneg.([ipsi,'_late']));
    passed_stats_nneg.([contra,'_early_num']) = nansum(passed_nneg.([contra,'_early']));
    passed_stats_nneg.([contra,'_late_num']) = nansum(passed_nneg.([contra,'_late']));
    passed_stats_nneg.early_num = nansum(passed_nneg.early);
    passed_stats_nneg.late_num = nansum(passed_nneg.late);
end


%% Return

passed_nneg = orderfields(passed_nneg);
passed_stats_nneg = orderfields(passed_stats_nneg);
end

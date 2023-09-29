function [passed,passed_stats] = tuningCriteria(shuffling,firingField,timewindow,ipsi,contra,p,prop)
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

passed.(ipsi) = temp1;
passed.(contra) = temp2;


%% Apply all tuning critera

% exclusive vs non-exclusive
passed.([ipsi,'only']) = floor((passed.(ipsi) + 1-passed.(contra))/2);
passed.([contra,'only']) = floor((1-passed.(ipsi) + passed.(contra))/2);
passed.([ipsi,'and',contra]) = floor((passed.(ipsi) + passed.(contra))/2);
passed.([ipsi,'or',contra]) = round((passed.(ipsi) + passed.(contra))/2);

% selectivity threshold
passed.([ipsi,'_sel']) = floor((passed.(ipsi) + double(firingField_ipsi.selectivity_ipsi>p.tng.selectivityThreshold))/2); 
passed.([contra,'_sel']) = floor((passed.(contra) + double(firingField_contra.selectivity_ipsi>p.tng.selectivityThreshold))/2); % firingField_contra.selectivity_ipsi is indeed correct
passed.([ipsi,'only_sel']) = floor((passed.([ipsi,'only']) + double(firingField_ipsi.selectivity_ipsi>p.tng.selectivityThreshold))/2); 
passed.([contra,'only_sel']) = floor((passed.([contra,'only']) + double(firingField_contra.selectivity_ipsi>p.tng.selectivityThreshold))/2); % firingField_contra.selectivity_ipsi is indeed correct

% early vs late
if strcmp(timewindow,'AW')
    passed.([ipsi,'_early']) = floor((passed.(ipsi) + firingField_ipsi.binary.early)/2);
    passed.([ipsi,'_late']) = floor((passed.(ipsi) + firingField_ipsi.binary.late)/2);
    passed.([contra,'_early']) = floor((passed.(contra) + firingField_contra.binary.early)/2);
    passed.([contra,'_late']) = floor((passed.(contra) + firingField_contra.binary.late)/2);
    passed.early = ceil((passed.([ipsi,'_early']) + passed.([contra,'_early']))/2);
    passed.late = ceil((passed.([ipsi,'_late']) + passed.([contra,'_late']))/2);
end


%% Summarise results

% basic
passed_stats.([ipsi,'_num']) = nansum(passed.(ipsi));
passed_stats.([contra,'_num']) = nansum(passed.(contra));
passed_stats.([ipsi,'_fractionOfCells']) = passed_stats.([ipsi,'_num'])/nansum(prop.iscell);
passed_stats.([contra,'_fractionOfCells']) = passed_stats.([contra,'_num'])/nansum(prop.iscell);

% exclusive vs non-exclusive
passed_stats.([ipsi,'or',contra,'_num']) = nansum(passed.([ipsi,'or',contra]));
passed_stats.([ipsi,'and',contra,'_num']) = nansum(passed.([ipsi,'and',contra]));
passed_stats.([ipsi,'or',contra,'_fractionOfCells']) = passed_stats.([ipsi,'or',contra,'_num'])/nansum(prop.iscell);
passed_stats.([ipsi,'and',contra,'_fractionOfCells']) = passed_stats.([ipsi,'and',contra,'_num'])/nansum(prop.iscell);
passed_stats.([ipsi,'and',contra,'_fractionOf',ipsi,'or',contra,'Cells']) = passed_stats.([ipsi,'and',contra,'_num'])/passed_stats.([ipsi,'or',contra,'_num']);

% only
passed_stats.([ipsi,'only_num']) = nansum(passed.([ipsi,'only']));
passed_stats.([contra,'only_num']) = nansum(passed.([contra,'only']));
passed_stats.([ipsi,'only_fractionOfCells']) = passed_stats.([ipsi,'only_num'])/nansum(prop.iscell);
passed_stats.([contra,'only_fractionOfCells']) = passed_stats.([contra,'only_num'])/nansum(prop.iscell);

% sel
passed_stats.([ipsi,'_sel_num']) = nansum(passed.([ipsi,'_sel']));
passed_stats.([contra,'_sel_num']) = nansum(passed.([contra,'_sel']));
passed_stats.([ipsi,'only_sel_num']) = nansum(passed.([ipsi,'only_sel']));
passed_stats.([contra,'only_sel_num']) = nansum(passed.([contra,'only_sel']));
passed_stats.([ipsi,'_sel_fractionOfCells']) = passed_stats.([ipsi,'_sel_num'])/nansum(prop.iscell);
passed_stats.([contra,'_sel_fractionOfCells']) = passed_stats.([contra,'_sel_num'])/nansum(prop.iscell);
passed_stats.([ipsi,'only_sel_fractionOfCells']) = passed_stats.([ipsi,'only_sel_num'])/nansum(prop.iscell);
passed_stats.([contra,'only_sel_fractionOfCells']) = passed_stats.([contra,'only_sel_num'])/nansum(prop.iscell);

% early vs late
if strcmp(timewindow,'AW')
    passed_stats.([ipsi,'_early_num']) = nansum(passed.([ipsi,'_early']));
    passed_stats.([ipsi,'_late_num']) = nansum(passed.([ipsi,'_late']));
    passed_stats.([contra,'_early_num']) = nansum(passed.([contra,'_early']));
    passed_stats.([contra,'_late_num']) = nansum(passed.([contra,'_late']));
    passed_stats.early_num = nansum(passed.early);
    passed_stats.late_num = nansum(passed.late);
end


%% Return

passed = orderfields(passed);
passed_stats = orderfields(passed_stats);
end

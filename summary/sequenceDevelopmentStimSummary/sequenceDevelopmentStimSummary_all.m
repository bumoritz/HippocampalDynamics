%% Get data

% get performance metrics
correct_catch = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if ~isempty(d{i,j})
            correct_catch(i,j) = nanmean(d{i,j}.perf.general.correct_catch);
        end
    end
end

% get warping
warp_exp1_a = nan(d_info.numAnimals,d_info.numDays);
warp_exp1_b = nan(d_info.numAnimals,d_info.numDays); 
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if isfield(d{i,j},'warp_all_stimVersion')
            warp_exp1_a(i,j) = d{i,j}.warp_all_stimVersion.exp1.a;
            warp_exp1_b(i,j) = d{i,j}.warp_all_stimVersion.exp1.b;
        end
    end
end


%% Simplify correct_catch and warp_exp1 data

correct_catch_switch1 = correct_catch(:,2:3);
correct_catch_switch2 = correct_catch(:,4:5);
correct_catch_all = [correct_catch_switch1;correct_catch_switch1];
correct_catch_switch1_seq = correct_catch_switch1(find(d_info.group==7),:);
correct_catch_switch1_ctrl = correct_catch_switch1(find(d_info.group==8),:);
correct_catch_switch2_seq = correct_catch_switch2(find(d_info.group==8),:);
correct_catch_switch2_ctrl = correct_catch_switch2(find(d_info.group==7),:);
correct_catch_seq = [correct_catch_switch1_seq;correct_catch_switch2_seq];
correct_catch_ctrl = [correct_catch_switch1_ctrl;correct_catch_switch2_ctrl];

warp_exp1_a_switch1 = warp_exp1_a(:,2:3);
warp_exp1_a_switch2 = warp_exp1_a(:,4:5);
warp_exp1_a_all = [warp_exp1_a_switch1;warp_exp1_a_switch1];
warp_exp1_a_switch1_seq = warp_exp1_a_switch1(find(d_info.group==7),:);
warp_exp1_a_switch1_ctrl = warp_exp1_a_switch1(find(d_info.group==8),:);
warp_exp1_a_switch2_seq = warp_exp1_a_switch2(find(d_info.group==8),:);
warp_exp1_a_switch2_ctrl = warp_exp1_a_switch2(find(d_info.group==7),:);
warp_exp1_a_seq = [warp_exp1_a_switch1_seq;warp_exp1_a_switch2_seq];
warp_exp1_a_ctrl = [warp_exp1_a_switch1_ctrl;warp_exp1_a_switch2_ctrl];

warp_exp1_b_switch1 = warp_exp1_b(:,2:3);
warp_exp1_b_switch2 = warp_exp1_b(:,4:5);
warp_exp1_b_all = [warp_exp1_b_switch1;warp_exp1_b_switch1];
warp_exp1_b_switch1_seq = warp_exp1_b_switch1(find(d_info.group==7),:);
warp_exp1_b_switch1_ctrl = warp_exp1_b_switch1(find(d_info.group==8),:);
warp_exp1_b_switch2_seq = warp_exp1_b_switch2(find(d_info.group==8),:);
warp_exp1_b_switch2_ctrl = warp_exp1_b_switch2(find(d_info.group==7),:);
warp_exp1_b_seq = [warp_exp1_b_switch1_seq;warp_exp1_b_switch2_seq];
warp_exp1_b_ctrl = [warp_exp1_b_switch1_ctrl;warp_exp1_b_switch2_ctrl];


%%

F = default_figure();

this_data = nan(d_info.numAnimals,4);
this_data(1:length(warp_exp1_b_seq(:,1)),1) = warp_exp1_b_seq(:,1);
this_data(1:length(warp_exp1_b_seq(:,2)),2) = warp_exp1_b_seq(:,2);
this_data(1:length(warp_exp1_b_ctrl(:,1)),3) = warp_exp1_b_ctrl(:,1);
this_data(1:length(warp_exp1_b_ctrl(:,2)),4) = warp_exp1_b_ctrl(:,2);
v = violinplot(this_data,{'seq stim (d1)','seq stim (d2)','ctrl stim (d1)','ctrl stim (d2)'});


%%

ranksum(warp_exp1_b_seq(:,1),warp_exp1_b_seq(:,2))
ranksum(warp_exp1_b_ctrl(:,1),warp_exp1_b_ctrl(:,2))










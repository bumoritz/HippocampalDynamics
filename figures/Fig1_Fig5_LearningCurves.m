%% Fig5_GrandFinale

% import data using Summary_Master with ops.do_learningCurveSummary = true;
% d_info = selectDataset(d_info,'-g0123456789101112-d12345-e01-r01-p1-l01-i00',sheet,path,ops);

save_root_fig = [path.root_summary,'figures\Fig5_fig\'];
save_root_png = [path.root_summary,'figures\Fig5_png\'];
save_root_pdf = [path.root_summary,'figures\Fig5_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig5_txt\'];

blocks_20t = [25,40,25,40,25];
blocks_100t = [5,8,5,8,5];


%% Get data

correct_20t = nan(d_info.numAnimals,sum(blocks_20t)); % [animal, block]
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if j==1
            temp0=0;
        else
            temp0 = cumsum(blocks_20t);
            temp0 = temp0(j-1);
        end
        if ~isempty(d{i,j})
            for k=1:blocks_20t(j)
                this_block = temp0+k;
                correct_20t(i,this_block) = d{i,j}.perf.blocks_general.correct(k);
            end
        end
    end
end

go_20t = nan(d_info.numAnimals,sum(blocks_20t)); % [animal, block]
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if j==1
            temp0=0;
        else
            temp0 = cumsum(blocks_20t);
            temp0 = temp0(j-1);
        end
        if ~isempty(d{i,j})
            for k=1:blocks_20t(j)
                this_block = temp0+k;
                go_20t(i,this_block) = d{i,j}.perf.blocks_general.H(k);
            end
        end
    end
end

nogo_20t = nan(d_info.numAnimals,sum(blocks_20t)); % [animal, block]
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if j==1
            temp0=0;
        else
            temp0 = cumsum(blocks_20t);
            temp0 = temp0(j-1);
        end
        if ~isempty(d{i,j})
            for k=1:blocks_20t(j)
                this_block = temp0+k;
                nogo_20t(i,this_block) = d{i,j}.perf.blocks_general.CR(k);
            end
        end
    end
end

correct_100t = nan(d_info.numAnimals,sum(blocks_100t)); % [animal, block]
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if j==1
            temp0=0;
        else
            temp0 = cumsum(blocks_100t);
            temp0 = temp0(j-1);
        end
        if ~isempty(d{i,j})
            for k=1:blocks_100t(j)
                this_block = temp0+k;
                these_blocks_20t = (k-1)*5+1:k*5;
                correct_100t(i,this_block) = nanmean(d{i,j}.perf.blocks_general.correct(these_blocks_20t));
            end
        end
    end
end

correct_general = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if ~isempty(d{i,j})
            correct_general(i,j) = d{i,j}.perf.general.correct;
        end
    end
end


%% Fit mixed-effect model - 100t

% this_data = correct_100t(:,blocks_100t(1)+1:end);
% this_data(~(d_info.group==7|d_info.group==8),:) = NaN;
% 
% tbl_in_100t.y = this_data(:);
% temp = repmat([1:13,1:13],d_info.numAnimals,1);
% tbl_in_100t.block = temp(:);
% temp = repmat([1:8,1:5,1:8,1:5],d_info.numAnimals,1);
% tbl_in_100t.block_sess = temp(:);
% tbl_in_100t.animal = nominal(repmat((1:d_info.numAnimals)',size(this_data,2),1));
% temp = repmat([1*ones(1,13),2*ones(1,13)],d_info.numAnimals,1);
% tbl_in_100t.switch = nominal(temp(:));
% temp = repmat([1*ones(1,8),2*ones(1,5),1*ones(1,8),2*ones(1,5)],d_info.numAnimals,1);
% tbl_in_100t.switchday = nominal(temp(:));
% temp = repmat([1*ones(1,13),2*ones(1,13)],d_info.numAnimals,1);
% tbl_in_100t.switch = nominal(temp(:));
% temp = nan(size(this_data));
% for j=2:5
%     if j==2
%         temp(:,1:8) = repmat(d_info.running(:,j),1,8);
%     elseif j==3
%         temp(:,9:13) = repmat(d_info.running(:,j),1,5);
%     elseif j==4
%         temp(:,14:21) = repmat(d_info.running(:,j),1,8);
%     elseif j==5
%         temp(:,22:26) = repmat(d_info.running(:,j),1,5);
%     end
% end
% tbl_in_100t.running = nominal(temp(:));
% % temp = zeros(size(this_data));
% % for i=1:d_info.numAnimals
% %     if d_info.group(i)==7
% %         temp(i,:) = [1*ones(1,13),2*ones(1,13)]; % 1=seq
% %     elseif d_info.group(i)==8
% %         temp(i,:) = [2*ones(1,13),1*ones(1,13)]; % 2=ctrl
% %     end
% % end
% % tbl_in.stim = nominal(temp(:));
% temp = nan(size(this_data));
% for i=1:d_info.numAnimals
%     if d_info.group(i)==7
%         temp(i,:) = [1*ones(1,13),2*ones(1,13)]; % 1=seq
%     elseif d_info.group(i)==8
%         temp(i,:) = [2*ones(1,13),1*ones(1,13)]; % 2=ctrl
%     end
% end
% tbl_in_100t.stim = nominal(temp(:));
% 
% tbl_100t = table(tbl_in_100t.animal,tbl_in_100t.switch,tbl_in_100t.switchday,tbl_in_100t.block,tbl_in_100t.block_sess,tbl_in_100t.running,tbl_in_100t.stim,tbl_in_100t.y,...
%     'VariableNames',{'animal','switch','switchday','block','block_sess','running','stim','y'});
% 
% lme_100t = fitlme(tbl_100t,'y ~ (running+stim+switchday)*block_sess + (1|animal) + (1|switch)')
% 
% %lme = fitlme(tbl,'y ~ stim*running*block_sess + stim*block_sess + running*block_sess + switchday*block_sess + (1|animal) + (1|switch)')
% 
% %lme = fitlme(tbl,'y ~ (running*stim*switchday*block) + (1|animal) + (1|switch)') % faulty
% 
% 
% 
% %lme = fitlme(tbl,'y ~ (running*stim*block) + (1|animal) + (1|switch)')
% 
% %lme = fitlme(tbl,'y ~ block + switchday + running:block + stim:block + running:stim:block + (1|animal) + (1|switch)')


%% Fit mixed-effect model - 20t

this_data = correct_20t(:,blocks_20t(1)+1:end);
this_data(~(d_info.group==7|d_info.group==8),:) = NaN;

tbl_in_20t.y = this_data(:) - 0.5;
temp = repmat([1:65,1:65],d_info.numAnimals,1);
tbl_in_20t.block = temp(:);
temp = repmat([1:40,1:25,1:40,1:25],d_info.numAnimals,1);
tbl_in_20t.block_sess = temp(:);
tbl_in_20t.animal = nominal(repmat((1:d_info.numAnimals)',size(this_data,2),1));
temp = repmat([1*ones(1,65),2*ones(1,65)],d_info.numAnimals,1);
tbl_in_20t.switch = nominal(temp(:));
temp = repmat([1*ones(1,40),2*ones(1,25),1*ones(1,40),2*ones(1,25)],d_info.numAnimals,1);
tbl_in_20t.switchday = nominal(temp(:));
temp = repmat([1*ones(1,65),2*ones(1,65)],d_info.numAnimals,1);
tbl_in_20t.switch = nominal(temp(:));
temp = nan(size(this_data));
for j=2:5
    if j==2
        temp(:,1:40) = repmat(d_info.running(:,j),1,40);
    elseif j==3
        temp(:,41:65) = repmat(d_info.running(:,j),1,25);
    elseif j==4
        temp(:,66:105) = repmat(d_info.running(:,j),1,40);
    elseif j==5
        temp(:,106:130) = repmat(d_info.running(:,j),1,25);
    end
end
tbl_in_20t.running = nominal(temp(:));
temp = nan(size(this_data));
for i=1:d_info.numAnimals
    if d_info.group(i)==7
        temp(i,:) = [2*ones(1,65),1*ones(1,65)]; % 2=seq
    elseif d_info.group(i)==8
        temp(i,:) = [1*ones(1,65),2*ones(1,65)]; % 1=ctrl
    end
end
tbl_in_20t.stim = nominal(temp(:));

tbl_20t = table(tbl_in_20t.animal,tbl_in_20t.switch,tbl_in_20t.switchday,tbl_in_20t.block,tbl_in_20t.block_sess,tbl_in_20t.running,tbl_in_20t.stim,tbl_in_20t.y,...
    'VariableNames',{'animal','switch','switchday','block','block_sess','running','stim','y'});

tbl_20t_nostim = tbl_20t;
tbl_20t_nostim(find(~isnan(double(tbl_20t_nostim.stim))),:) = [];

% NOT USED
% lme_20t = fitlme(tbl_20t,'y ~ (running+stim+switchday)*block_sess + (1|animal) + (1|switch)') % stim not significant
% lme_20t = fitlme(tbl_20t,'y ~ (running+stim+switchday)*block + (1|animal) + (1|switch)') % not ideal
% lme_20t = fitlme(tbl_20t,'y ~ (running*stim*switchday*block_sess) + (1|animal) + (1|switch)') % not ideal (not ideal, overfitting triple interaction)
% lme_20t = fitlme(tbl_20t,'y ~ (running*stim*switchday*block) + (1|animal) + (1|switch)') % not ideal (almost good but gives huge weight to switch day intercept)
% lme_20t = fitlme(tbl_20t,'y ~ (running*stim*block) + (1|animal) + (1|switch)') % best
% lme_20t = fitlme(tbl_20t,'y ~ (running*stim*block) + (1|animal) + (1|switch)')
% lme_20t = fitlme(tbl,'y ~ (running*stim*switchday*block) + (1|animal) + (1|switch)') % faulty
% lme_20t = fitlme(tbl,'y ~ (running*stim*block) + (1|animal) + (1|switch)')

% USED - stim data
%lme_20t = fitlme(tbl_20t,'y ~ (running+stim+switchday)*block_sess + (1|animal) + (1|switch)') % not ideal
lme_20t = fitlme(tbl_20t,'y ~ (running+stim)*block_sess + (1|animal) + (1|switch) + (1|switchday)')

% USED - no-stim data
%lme_20t = fitlme(tbl_20t,'y ~ (running+switchday)*block_sess + (1|animal) + (1|switch)')
%lme_20t = fitlme(tbl_20t,'y ~ running*switchday*block_sess + (1|animal) + (1|switch)')
lme_20t = fitlme(tbl_20t,'y ~ running*block_sess + (1|animal) + (1|switch) + (1|switchday)')


%%

lme_20t = fitlme(tbl_20t,'y ~ (running*stim*block) + (1|animal) + (1|switch)') % best
 
these_labels = categorical({'(offset)','(slope)','Running (offset)','Running (slope)','Seq stim (initial)','Seq stim (slope)','Running * seq stim (initial)','Running * seq stim (slope)'});
these_labels = reordercats(these_labels,{'(offset)','(slope)','Running (offset)','Running (slope)','Seq stim (initial)','Seq stim (slope)','Running * seq stim (initial)','Running * seq stim (slope)'});

slope_scaling = 5; % 100t
this_order = [1,2,3,5,4,6,7,8];
this_data = double(lme_20t.Coefficients(:,2));
this_data = this_data(this_order);
this_ci_lower = double(lme_20t.Coefficients(:,7));
this_ci_lower = this_data - this_ci_lower(this_order);
this_ci_lower([2,4,6,8]) = this_ci_lower([2,4,6,8])*slope_scaling;
this_ci_upper = double(lme_20t.Coefficients(:,8));
this_ci_upper = this_data - this_ci_upper(this_order);
this_ci_upper([2,4,6,8]) = this_ci_upper([2,4,6,8])*slope_scaling;
this_data([2,4,6,8]) = this_data([2,4,6,8])*slope_scaling;

figure; hold on;
%F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
v=bar(these_labels,this_data*100);
h=errorbar(these_labels,this_data*100,this_ci_upper*100,this_ci_lower*100);
v.FaceColor = 'flat'; v.EdgeColor = 'none'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.darkGray; v.CData(3,:) = p.col.runner; v.CData(4,:) = p.col.runner; v.CData(5,:) = p.col.seq; v.CData(6,:) = p.col.seq; v.CData(7,:) = mean([p.col.runner;p.col.seq]); v.CData(8,:) = mean([p.col.runner;p.col.seq]);
h.Color = p.col.black; h.LineStyle = 'none'; h.LineWidth = 1;  
                          

ytickformat('percentage')
ylabel('Performance regression coefficient')



%% ---- %% ---- %% ---- %%

%% Fig5_PSF

save_root_fig = [path.root_summary,'figures\Fig5_fig\'];
save_root_png = [path.root_summary,'figures\Fig5_png\'];
save_root_pdf = [path.root_summary,'figures\Fig5_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig5_txt\'];


%% Extract data

this_data = correct_20t(:,blocks_20t(1)+1:end);
this_data_go = go_20t(:,blocks_20t(1)+1:end);
this_data_nogo = nogo_20t(:,blocks_20t(1)+1:end);

this_data_2 = correct_20t(:,1:blocks_20t(1));
this_data_2_go = go_20t(:,1:blocks_20t(1));
this_data_2_nogo = nogo_20t(:,1:blocks_20t(1));

this_running = nan(size(this_data));
for j=2:5
    if j==2
        this_running(:,1:40) = repmat(d_info.running(:,j),1,40);
    elseif j==3
        this_running(:,41:65) = repmat(d_info.running(:,j),1,25);
    elseif j==4
        this_running(:,66:105) = repmat(d_info.running(:,j),1,40);
    elseif j==5
        this_running(:,106:130) = repmat(d_info.running(:,j),1,25);
    end
end

this_stim = nan(size(this_data));
for i=1:d_info.numAnimals
    if d_info.group(i)==7
        this_stim(i,:) = [2*ones(1,65),1*ones(1,65)]; % 2=seq
    elseif d_info.group(i)==8
        this_stim(i,:) = [1*ones(1,65),2*ones(1,65)]; % 1=ctrl
    else
        this_stim(i,:) = zeros(1,130);
    end
end

temp = nan(size(this_data));
for i=1:d_info.numAnimals
    for j=1:size(this_data,2)
        this_condition = this_stim==2;
        if this_condition(i,j)
            temp(i,j) = this_data(i,j);
        end
    end
end
this_data_seq = nanmean(cat(3,temp(:,1:65),temp(:,66:end)),3);
this_data_seq = [smoothdata(this_data_seq(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_seq(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(this_data));
for i=1:d_info.numAnimals
    for j=1:size(this_data,2)
        this_condition = this_running==1 & this_stim==2;
        if this_condition(i,j)
            temp(i,j) = this_data(i,j);
        end
    end
end
this_data_runner_seq = nanmean(cat(3,temp(:,1:65),temp(:,66:end)),3);
this_data_runner_seq = [smoothdata(this_data_runner_seq(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_runner_seq(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(this_data));
for i=1:d_info.numAnimals
    for j=1:size(this_data,2)
        this_condition = this_running==0 & this_stim==2;
        if this_condition(i,j)
            temp(i,j) = this_data(i,j);
        end
    end
end
this_data_nonrunner_seq = nanmean(cat(3,temp(:,1:65),temp(:,66:end)),3);
this_data_nonrunner_seq = [smoothdata(this_data_nonrunner_seq(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nonrunner_seq(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(this_data));
for i=1:d_info.numAnimals
    for j=1:size(this_data,2)
        this_condition = this_stim==1;
        if this_condition(i,j)
            temp(i,j) = this_data(i,j);
        end
    end
end
this_data_ctrl = nanmean(cat(3,temp(:,1:65),temp(:,66:end)),3);
this_data_ctrl = [smoothdata(this_data_ctrl(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_ctrl(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(this_data));
for i=1:d_info.numAnimals
    for j=1:size(this_data,2)
        this_condition = this_running==1 & this_stim==1;
        if this_condition(i,j)
            temp(i,j) = this_data(i,j);
        end
    end
end
this_data_runner_ctrl = nanmean(cat(3,temp(:,1:65),temp(:,66:end)),3);
this_data_runner_ctrl = [smoothdata(this_data_runner_ctrl(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_runner_ctrl(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(this_data));
for i=1:d_info.numAnimals
    for j=1:size(this_data,2)
        this_condition = this_running==0 & this_stim==1;
        if this_condition(i,j)
            temp(i,j) = this_data(i,j);
        end
    end
end
this_data_nonrunner_ctrl = nanmean(cat(3,temp(:,1:65),temp(:,66:end)),3);
this_data_nonrunner_ctrl = [smoothdata(this_data_nonrunner_ctrl(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nonrunner_ctrl(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(this_data));
for i=1:d_info.numAnimals
    for j=1:size(this_data,2)
        this_condition = this_stim==0;
        if this_condition(i,j)
            temp(i,j) = this_data(i,j);
        end
    end
end
this_data_nostim_sep = temp;
this_data_nostim_sep = [smoothdata(this_data_nostim_sep(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nostim_sep(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nostim_sep(:,66:105),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nostim_sep(:,106:130),2,'gaussian',p.lcs.smoothing_sd*5)];
this_data_nostim = nanmean(cat(3,temp(:,1:65),temp(:,66:end)),3);
this_data_nostim = [smoothdata(this_data_nostim(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nostim(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(this_data_go));
for i=1:d_info.numAnimals
    for j=1:size(this_data_go,2)
        this_condition = this_stim==0;
        if this_condition(i,j)
            temp(i,j) = this_data_go(i,j);
        end
    end
end
this_data_nostim_sep_go = temp;
this_data_nostim_sep_go = [smoothdata(this_data_nostim_sep_go(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nostim_sep_go(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nostim_sep_go(:,66:105),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nostim_sep_go(:,106:130),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(this_data_nogo));
for i=1:d_info.numAnimals
    for j=1:size(this_data_nogo,2)
        this_condition = this_stim==0;
        if this_condition(i,j)
            temp(i,j) = this_data_nogo(i,j);
        end
    end
end
this_data_nostim_sep_nogo = temp;
this_data_nostim_sep_nogo = [smoothdata(this_data_nostim_sep_nogo(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nostim_sep_nogo(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nostim_sep_nogo(:,66:105),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nostim_sep_nogo(:,106:130),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(this_data_2));
for i=1:d_info.numAnimals
    for j=1:size(this_data_2,2)
        temp(i,j) = this_data_2(i,j);
    end
end
this_data_nostim_exp = temp;
this_data_nostim_exp = smoothdata(this_data_nostim_exp,2,'gaussian',p.lcs.smoothing_sd*5);

temp = nan(size(this_data_2));
for i=1:d_info.numAnimals
    if d_info.running(i,1)==1
        for j=1:size(this_data_2,2)
            temp(i,j) = this_data_2(i,j);
        end
    end
end
this_data_runner_nostim_exp = temp;
this_data_runner_nostim_exp = smoothdata(this_data_runner_nostim_exp,2,'gaussian',p.lcs.smoothing_sd*5);

temp = nan(size(this_data_2));
for i=1:d_info.numAnimals
    if d_info.running(i,1)==0
        for j=1:size(this_data_2,2)
            temp(i,j) = this_data_2(i,j);
        end
    end
end
this_data_nonrunner_nostim_exp = temp;
this_data_nonrunner_nostim_exp = smoothdata(this_data_nonrunner_nostim_exp,2,'gaussian',p.lcs.smoothing_sd*5);

temp = nan(size(this_data_2_go));
for i=1:d_info.numAnimals
    for j=1:size(this_data_2_go,2)
        temp(i,j) = this_data_2_go(i,j);
    end
end
this_data_nostim_exp_go = temp;
this_data_nostim_exp_go = smoothdata(this_data_nostim_exp_go,2,'gaussian',p.lcs.smoothing_sd*5);

temp = nan(size(this_data_2_nogo));
for i=1:d_info.numAnimals
    for j=1:size(this_data_2_nogo,2)
        temp(i,j) = this_data_2_nogo(i,j);
    end
end
this_data_nostim_exp_nogo = temp;
this_data_nostim_exp_nogo = smoothdata(this_data_nostim_exp_nogo,2,'gaussian',p.lcs.smoothing_sd*5);

temp = nan(size(this_data));
for i=1:d_info.numAnimals
    for j=1:size(this_data,2)
        this_condition = this_running==1 & this_stim==0;
        if this_condition(i,j)
            temp(i,j) = this_data(i,j);
        end
    end
end
this_data_runner_nostim_sep = temp;
this_data_runner_nostim_sep = [smoothdata(this_data_runner_nostim_sep(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_runner_nostim_sep(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_runner_nostim_sep(:,66:105),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_runner_nostim_sep(:,106:130),2,'gaussian',p.lcs.smoothing_sd*5)];
this_data_runner_nostim = nanmean(cat(3,temp(:,1:65),temp(:,66:end)),3);
this_data_runner_nostim = [smoothdata(this_data_runner_nostim(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_runner_nostim(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(this_data));
for i=1:d_info.numAnimals
    for j=1:size(this_data,2)
        this_condition = this_running==0 & this_stim==0;
        if this_condition(i,j)
            temp(i,j) = this_data(i,j);
        end
    end
end
this_data_nonrunner_nostim_sep = temp;
this_data_nonrunner_nostim_sep = [smoothdata(this_data_nonrunner_nostim_sep(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nonrunner_nostim_sep(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nonrunner_nostim_sep(:,66:105),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nonrunner_nostim_sep(:,106:130),2,'gaussian',p.lcs.smoothing_sd*5)];
this_data_nonrunner_nostim = nanmean(cat(3,temp(:,1:65),temp(:,66:end)),3);
this_data_nonrunner_nostim = [smoothdata(this_data_nonrunner_nostim(:,1:40),2,'gaussian',p.lcs.smoothing_sd*5),smoothdata(this_data_nonrunner_nostim(:,41:65),2,'gaussian',p.lcs.smoothing_sd*5)];

temp = nan(size(correct_general));
for i=1:d_info.numAnimals
    for j=1:size(correct_general,2)
        temp(i,j) = correct_general(i,j);
    end
end
this_data_general = temp;

temp = nan(size(correct_general));
for i=1:d_info.numAnimals
    for j=1:size(correct_general,2)
        if d_info.running(i,j)==1
            temp(i,j) = correct_general(i,j);
        end
    end
end
this_data_general_runner = temp;

temp = nan(size(correct_general));
for i=1:d_info.numAnimals
    for j=1:size(correct_general,2)
        if d_info.running(i,j)==0
            temp(i,j) = correct_general(i,j);
        end
    end
end
this_data_general_nonrunner = temp;


%% Fig5_LearningCurves_all_ind

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

xline(40.5,'k-');
yline(50,'k:');

plot(1:size(this_data_ctrl,2),this_data_ctrl*100,'Color',p.col.ctrl,'LineWidth',0.5);
plot(1:size(this_data_seq,2),this_data_seq*100,'Color',p.col.seq,'LineWidth',0.5);

xlim([0,65])
xticks([0:5:65])
xticklabels({'0','','200','','400','','600','','800','','1000','','1200',''})
xlabel('Trial')
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig5_LearningCurves_all_ind.fig']);
saveas(F,[save_root_png,'\Fig5_LearningCurves_all_ind.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_LearningCurves_all_ind.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_LearningCurves_all_avg

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

xline(40.5,'k-');
yline(50,'k:');
shadedErrorBar(1:size(this_data_nostim,2),nanmean(this_data_nostim,1)*100,nansem(this_data_nostim,1)*100,'lineProps',p.col.darkGray);
shadedErrorBar(1:size(this_data_ctrl,2),nanmean(this_data_ctrl,1)*100,nansem(this_data_ctrl,1)*100,'lineProps',p.col.ctrl);
shadedErrorBar(1:size(this_data_seq,2),nanmean(this_data_seq,1)*100,nansem(this_data_seq,1)*100,'lineProps',p.col.seq);
xlim([0,65])
xticks([0:5:65])
xticklabels({'0','','200','','400','','600','','800','','1000','','1200',''})
xlabel('Trial')
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig5_LearningCurves_all_avg.fig']);
saveas(F,[save_root_png,'\Fig5_LearningCurves_all_avg.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_LearningCurves_all_avg.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_LearningCurves_R_ind

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

xline(40.5,'k-');
yline(50,'k:');

plot(1:size(this_data_runner_ctrl,2),this_data_runner_ctrl*100,'Color',p.col.ctrl,'LineWidth',0.5);
plot(1:size(this_data_runner_seq,2),this_data_runner_seq*100,'Color',p.col.seq,'LineWidth',0.5);

xlim([0,65])
xticks([0:5:65])
xticklabels({'0','','200','','400','','600','','800','','1000','','1200',''})
xlabel('Trial')
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig5_LearningCurves_R_ind.fig']);
saveas(F,[save_root_png,'\Fig5_LearningCurves_R_ind.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_LearningCurves_R_ind.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_LearningCurves_R_avg

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

xline(40.5,'k-');
yline(50,'k:');
shadedErrorBar(1:size(this_data_runner_nostim,2),nanmean(this_data_runner_nostim,1)*100,nansem(this_data_runner_nostim,1)*100,'lineProps',p.col.darkGray);
shadedErrorBar(1:size(this_data_runner_ctrl,2),nanmean(this_data_runner_ctrl,1)*100,nansem(this_data_runner_ctrl,1)*100,'lineProps',p.col.ctrl);
shadedErrorBar(1:size(this_data_runner_seq,2),nanmean(this_data_runner_seq,1)*100,nansem(this_data_runner_seq,1)*100,'lineProps',p.col.seq);
xlim([0,65])
xticks([0:5:65])
xticklabels({'0','','200','','400','','600','','800','','1000','','1200',''})
xlabel('Trial')
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig5_LearningCurves_R_avg.fig']);
saveas(F,[save_root_png,'\Fig5_LearningCurves_R_avg.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_LearningCurves_R_avg.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_LearningCurves_NR_ind

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

xline(40.5,'k-');
yline(50,'k:');

plot(1:size(this_data_nonrunner_ctrl,2),this_data_nonrunner_ctrl*100,'Color',p.col.ctrl,'LineWidth',0.5);
plot(1:size(this_data_nonrunner_seq,2),this_data_nonrunner_seq*100,'Color',p.col.seq,'LineWidth',0.5);

xlim([0,65])
xticks([0:5:65])
xticklabels({'0','','200','','400','','600','','800','','1000','','1200',''})
xlabel('Trial')
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig5_LearningCurves_NR_ind.fig']);
saveas(F,[save_root_png,'\Fig5_LearningCurves_NR_ind.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_LearningCurves_NR_ind.pdf']); set(gcf,'Color',[1,1,1])


%% Fig5_LearningCurves_NR_avg

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

xline(40.5,'k-');
yline(50,'k:');
shadedErrorBar(1:size(this_data_nonrunner_nostim,2),nanmean(this_data_nonrunner_nostim,1)*100,nansem(this_data_nonrunner_nostim,1)*100,'lineProps',p.col.darkGray);
shadedErrorBar(1:size(this_data_nonrunner_ctrl,2),nanmean(this_data_nonrunner_ctrl,1)*100,nansem(this_data_nonrunner_ctrl,1)*100,'lineProps',p.col.ctrl);
shadedErrorBar(1:size(this_data_nonrunner_seq,2),nanmean(this_data_nonrunner_seq,1)*100,nansem(this_data_nonrunner_seq,1)*100,'lineProps',p.col.seq);
xlim([0,65])
xticks([0:5:65])
xticklabels({'0','','200','','400','','600','','800','','1000','','1200',''})
xlabel('Trial')
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig5_LearningCurves_NR_avg.fig']);
saveas(F,[save_root_png,'\Fig5_LearningCurves_NR_avg.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_LearningCurves_NR_avg.pdf']); set(gcf,'Color',[1,1,1])





%% Fig1_LearningCurves

save_root_fig = [path.root_summary,'figures\Fig1_fig\'];
save_root_png = [path.root_summary,'figures\Fig1_png\'];
save_root_pdf = [path.root_summary,'figures\Fig1_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig1_txt\'];


%% Fig1_LearningCurves_all_d15_avg

spacer = 1;
F = paper_figure([0,0.5,mm2inch(3.25*34),mm2inch(0.85*34)]); hold on;

for i=1:5
    if i==1
        subplot(1,5+8+5+8+5+4*spacer,1:5)
    elseif i==2
        subplot(1,5+8+5+8+5+4*spacer,(6:13)+spacer)
    elseif i==3
        subplot(1,5+8+5+8+5+4*spacer,(14:18)+2*spacer)
    elseif i==4
        subplot(1,5+8+5+8+5+4*spacer,(19:26)+3*spacer)
    elseif i==5
        subplot(1,5+8+5+8+5+4*spacer,(27:31)+4*spacer)
    end
    yline(50,'k:');
    if i==1
        shadedErrorBar(1:size(this_data_nostim_exp,2),nanmean(this_data_nostim_exp,1)*100,nansem(this_data_nostim_exp,1)*100,'lineProps',p.col.black);
        xlabel('')
        ylabel('Performance')
    elseif i==2
        shadedErrorBar(1:size(this_data_nostim_sep(:,1:40),2),nanmean(this_data_nostim_sep(:,1:40),1)*100,nansem(this_data_nostim_sep(:,1:40),1)*100,'lineProps',p.col.black);
        xlabel('')
    elseif i==3
        shadedErrorBar(1:size(this_data_nostim_sep(:,41:65),2),nanmean(this_data_nostim_sep(:,41:65),1)*100,nansem(this_data_nostim_sep(:,41:65),1)*100,'lineProps',p.col.black);
        xlabel('Trial')
    elseif i==4
        shadedErrorBar(1:size(this_data_nostim_sep(:,66:105),2),nanmean(this_data_nostim_sep(:,66:105),1)*100,nansem(this_data_nostim_sep(:,66:105),1)*100,'lineProps',p.col.black);
        xlabel('')
    elseif i==5
        shadedErrorBar(1:size(this_data_nostim_sep(:,106:130),2),nanmean(this_data_nostim_sep(:,106:130),1)*100,nansem(this_data_nostim_sep(:,106:130),1)*100,'lineProps',p.col.black);
        xlabel('')
    end
    ylim([35,100])
    yticks([50,100])
    ytickformat('percentage')
    if i~=1
        yticklabels({'',''})
    end
    if i==1 || i==3 || i==5
        xlim([0,25])
        xticks([0:5:25])
        xticklabels({'0','','','','','500'})
    elseif i==2 || i==4
        xlim([0,40])
        xticks([0:5:40])
        xticklabels({'0','','','','','','','','800'})
    end
end

savefig(F,[save_root_fig,'\Fig1_LearningCurves_all_d15_avg.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_all_d15_avg.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_all_d15_avg.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_R_NR_d1_avg

labels = categorical({'Non-runners','Runners'});
labels = reordercats(labels,{'Non-runners','Runners'});

this_data = nan(d_info.numAnimals,2);
this_data_nonrunner = this_data_general_nonrunner(:,1);
this_data_runner = this_data_general_runner(:,1);
this_data(1:length(this_data_nonrunner),1) = this_data_nonrunner;
this_data(1:length(this_data_runner),2) = this_data_runner;

%[pval,~,stats] = ranksum(this_data_nonrunner,this_data_runner);

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
yline(50,'k:');
   
v = violinplot(this_data*100,labels);%,'bandwidth',4); % to find out what the default bandwidth for each violin is [density, value, bandwidth] = ksdensity(this_data(:,2)*100);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;

ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig1_LearningCurves_R_NR_d1_avg.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_R_NR_d1_avg.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_R_NR_d1_avg.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_R_NR_d3_avg

labels = categorical({'Non-runners','Runners'});
labels = reordercats(labels,{'Non-runners','Runners'});

this_data = nan(d_info.numAnimals,2);
this_data_nonrunner = this_data_general_nonrunner(:,3);
this_data_runner = this_data_general_runner(:,3);
this_data(1:length(this_data_nonrunner),1) = this_data_nonrunner;
this_data(1:length(this_data_runner),2) = this_data_runner;

[pval,~,stats] = ranksum(this_data_nonrunner,this_data_runner)

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
yline(50,'k:');
   
v = violinplot(this_data*100,labels);%,'bandwidth',4); % to find out what the default bandwidth for each violin is [density, value, bandwidth] = ksdensity(this_data(:,2)*100);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;

ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig1_LearningCurves_R_NR_d3_avg.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_R_NR_d3_avg.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_R_NR_d3_avg.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_R_NR_d5_avg

labels = categorical({'Non-runners','Runners'});
labels = reordercats(labels,{'Non-runners','Runners'});

this_data = nan(d_info.numAnimals,2);
this_data_nonrunner = this_data_general_nonrunner(:,5);
this_data_runner = this_data_general_runner(:,5);
this_data(1:length(this_data_nonrunner),1) = this_data_nonrunner;
this_data(1:length(this_data_runner),2) = this_data_runner;

[pval,~,stats] = ranksum(this_data_nonrunner,this_data_runner)

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
yline(50,'k:');
   
v = violinplot(this_data*100,labels);%,'bandwidth',4); % to find out what the default bandwidth for each violin is [density, value, bandwidth] = ksdensity(this_data(:,2)*100);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;

ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig1_LearningCurves_R_NR_d5_avg.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_R_NR_d5_avg.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_R_NR_d5_avg.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_R_NR_d35_avg

labels = categorical({'Non-runners','Runners'});
labels = reordercats(labels,{'Non-runners','Runners'});

this_data = nan(d_info.numAnimals,2);
this_data_nonrunner = nanmean([this_data_general_nonrunner(:,3),this_data_general_nonrunner(:,5)],2);
this_data_runner = nanmean([this_data_general_runner(:,3),this_data_general_runner(:,5)],2);
this_data(1:length(this_data_nonrunner),1) = this_data_nonrunner;
this_data(1:length(this_data_runner),2) = this_data_runner;

[pval,~,stats] = ranksum(this_data_nonrunner,this_data_runner)

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
yline(50,'k:');
   
v = violinplot(this_data*100,labels);%,'bandwidth',4); % to find out what the default bandwidth for each violin is [density, value, bandwidth] = ksdensity(this_data(:,2)*100);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;

ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig1_LearningCurves_R_NR_d35_avg.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_R_NR_d35_avg.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_R_NR_d35_avg.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_R_NR_d25_avg

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

xline(40.5,'k-');
yline(50,'k:');
shadedErrorBar(1:size(this_data_nonrunner_nostim,2),nanmean(this_data_nonrunner_nostim,1)*100,nansem(this_data_nonrunner_nostim,1)*100,'lineProps',p.col.nonrunner);
shadedErrorBar(1:size(this_data_runner_nostim,2),nanmean(this_data_runner_nostim,1)*100,nansem(this_data_runner_nostim,1)*100,'lineProps',p.col.runner);
xlim([0,65])
xticks([0:5:65])
xticklabels({'0','','200','','400','','600','','800','','1000','','1200',''})
xlabel('Trial')
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig1_LearningCurves_R_NR_d25_avg.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_R_NR_d25_avg.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_R_NR_d25_avg.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_all_d15_avg_gonogo

spacer = 1;
F = paper_figure([0,0.5,mm2inch(3.25*34),mm2inch(0.85*34)]); hold on;

for i=1:5
    if i==1
        subplot(1,5+8+5+8+5+4*spacer,1:5)
    elseif i==2
        subplot(1,5+8+5+8+5+4*spacer,(6:13)+spacer)
    elseif i==3
        subplot(1,5+8+5+8+5+4*spacer,(14:18)+2*spacer)
    elseif i==4
        subplot(1,5+8+5+8+5+4*spacer,(19:26)+3*spacer)
    elseif i==5
        subplot(1,5+8+5+8+5+4*spacer,(27:31)+4*spacer)
    end
    yline(50,'k:');
    if i==1
        shadedErrorBar(1:size(this_data_nostim_exp_nogo,2),nanmean(this_data_nostim_exp_nogo,1)*100,nansem(this_data_nostim_exp_nogo,1)*100,'lineProps',p.col.darkGray);
        shadedErrorBar(1:size(this_data_nostim_exp_go,2),nanmean(this_data_nostim_exp_go,1)*100,nansem(this_data_nostim_exp_go,1)*100,'lineProps',p.col.reward);
        shadedErrorBar(1:size(this_data_nostim_exp,2),nanmean(this_data_nostim_exp,1)*100,nansem(this_data_nostim_exp,1)*100,'lineProps',p.col.black);
        xlabel('')
        ylabel('Performance')
    elseif i==2
        shadedErrorBar(1:size(this_data_nostim_sep_nogo(:,1:40),2),nanmean(this_data_nostim_sep_nogo(:,1:40),1)*100,nansem(this_data_nostim_sep_nogo(:,1:40),1)*100,'lineProps',p.col.darkGray);
        shadedErrorBar(1:size(this_data_nostim_sep_go(:,1:40),2),nanmean(this_data_nostim_sep_go(:,1:40),1)*100,nansem(this_data_nostim_sep_go(:,1:40),1)*100,'lineProps',p.col.reward);
        shadedErrorBar(1:size(this_data_nostim_sep(:,1:40),2),nanmean(this_data_nostim_sep(:,1:40),1)*100,nansem(this_data_nostim_sep(:,1:40),1)*100,'lineProps',p.col.black);
        xlabel('')
    elseif i==3
        shadedErrorBar(1:size(this_data_nostim_sep_nogo(:,41:65),2),nanmean(this_data_nostim_sep_nogo(:,41:65),1)*100,nansem(this_data_nostim_sep_nogo(:,41:65),1)*100,'lineProps',p.col.darkGray);
        shadedErrorBar(1:size(this_data_nostim_sep_go(:,41:65),2),nanmean(this_data_nostim_sep_go(:,41:65),1)*100,nansem(this_data_nostim_sep_go(:,41:65),1)*100,'lineProps',p.col.reward);
        shadedErrorBar(1:size(this_data_nostim_sep(:,41:65),2),nanmean(this_data_nostim_sep(:,41:65),1)*100,nansem(this_data_nostim_sep(:,41:65),1)*100,'lineProps',p.col.black);
        xlabel('Trial')
    elseif i==4
        shadedErrorBar(1:size(this_data_nostim_sep_nogo(:,66:105),2),nanmean(this_data_nostim_sep_nogo(:,66:105),1)*100,nansem(this_data_nostim_sep_nogo(:,66:105),1)*100,'lineProps',p.col.darkGray);
        shadedErrorBar(1:size(this_data_nostim_sep_go(:,66:105),2),nanmean(this_data_nostim_sep_go(:,66:105),1)*100,nansem(this_data_nostim_sep_go(:,66:105),1)*100,'lineProps',p.col.reward);
        shadedErrorBar(1:size(this_data_nostim_sep(:,66:105),2),nanmean(this_data_nostim_sep(:,66:105),1)*100,nansem(this_data_nostim_sep(:,66:105),1)*100,'lineProps',p.col.black);
        xlabel('')
    elseif i==5
        shadedErrorBar(1:size(this_data_nostim_sep_nogo(:,106:130),2),nanmean(this_data_nostim_sep_nogo(:,106:130),1)*100,nansem(this_data_nostim_sep_nogo(:,106:130),1)*100,'lineProps',p.col.darkGray);
        shadedErrorBar(1:size(this_data_nostim_sep_go(:,106:130),2),nanmean(this_data_nostim_sep_go(:,106:130),1)*100,nansem(this_data_nostim_sep_go(:,106:130),1)*100,'lineProps',p.col.reward);
        shadedErrorBar(1:size(this_data_nostim_sep(:,106:130),2),nanmean(this_data_nostim_sep(:,106:130),1)*100,nansem(this_data_nostim_sep(:,106:130),1)*100,'lineProps',p.col.black);
        xlabel('')
    end
    ylim([0,100])
    yticks([0,50,100])
    ytickformat('percentage')
    if i~=1
        yticklabels({'',''})
    end
    if i==1 || i==3 || i==5
        xlim([0,25])
        xticks([0:5:25])
        xticklabels({'0','','','','','500'})
    elseif i==2 || i==4
        xlim([0,40])
        xticks([0:5:40])
        xticklabels({'0','','','','','','','','800'})
    end
end

savefig(F,[save_root_fig,'\Fig1_LearningCurves_all_d15_avg_gonogo.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_all_d15_avg_gonogo.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_all_d15_avg_gonogo.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_R_NR_d15_ind

spacer = 1;
F = paper_figure([0,0.5,mm2inch(3.25*34),mm2inch(0.85*34)]); hold on;

for i=1:5
    if i==1
        subplot(1,5+8+5+8+5+4*spacer,1:5); hold on;
    elseif i==2
        subplot(1,5+8+5+8+5+4*spacer,(6:13)+spacer); hold on;
    elseif i==3
        subplot(1,5+8+5+8+5+4*spacer,(14:18)+2*spacer); hold on;
    elseif i==4
        subplot(1,5+8+5+8+5+4*spacer,(19:26)+3*spacer); hold on;
    elseif i==5
        subplot(1,5+8+5+8+5+4*spacer,(27:31)+4*spacer); hold on;
    end
    yline(50,'k:');
    if i==1
        plot(1:size(this_data_nonrunner_nostim_exp,2),this_data_nonrunner_nostim_exp*100,'Color',p.col.nonrunner,'LineWidth',0.5);
        plot(1:size(this_data_runner_nostim_exp,2),this_data_runner_nostim_exp*100,'Color',p.col.runner,'LineWidth',0.5);
        xlabel('')
        ylabel('Performance')
    elseif i==2
        plot(1:size(this_data_nonrunner_nostim_sep(:,1:40),2),this_data_nonrunner_nostim_sep(:,1:40)*100,'Color',p.col.nonrunner,'LineWidth',0.5);
        plot(1:size(this_data_runner_nostim_sep(:,1:40),2),this_data_runner_nostim_sep(:,1:40)*100,'Color',p.col.runner,'LineWidth',0.5);
        xlabel('')
    elseif i==3
        plot(1:size(this_data_nonrunner_nostim_sep(:,41:65),2),this_data_nonrunner_nostim_sep(:,41:65)*100,'Color',p.col.nonrunner,'LineWidth',0.5);
        plot(1:size(this_data_runner_nostim_sep(:,41:65),2),this_data_runner_nostim_sep(:,41:65)*100,'Color',p.col.runner,'LineWidth',0.5);
        xlabel('Trial')
    elseif i==4
        plot(1:size(this_data_nonrunner_nostim_sep(:,66:105),2),this_data_nonrunner_nostim_sep(:,66:105)*100,'Color',p.col.nonrunner,'LineWidth',0.5);
        plot(1:size(this_data_runner_nostim_sep(:,66:105),2),this_data_runner_nostim_sep(:,66:105)*100,'Color',p.col.runner,'LineWidth',0.5);
        xlabel('')
    elseif i==5
        plot(1:size(this_data_nonrunner_nostim_sep(:,106:130),2),this_data_nonrunner_nostim_sep(:,106:130)*100,'Color',p.col.nonrunner,'LineWidth',0.5);
        plot(1:size(this_data_runner_nostim_sep(:,106:130),2),this_data_runner_nostim_sep(:,106:130)*100,'Color',p.col.runner,'LineWidth',0.5);
        xlabel('')
    end
    ylim([35,100])
    yticks([50,100])
    ytickformat('percentage')
    if i~=1
        yticklabels({'',''})
    end
    if i==1 || i==3 || i==5
        xlim([0,25])
        xticks([0:5:25])
        xticklabels({'0','','','','','500'})
    elseif i==2 || i==4
        xlim([0,40])
        xticks([0:5:40])
        xticklabels({'0','','','','','','','','800'})
    end
end

savefig(F,[save_root_fig,'\Fig1_LearningCurves_R_NR_d15_ind.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_R_NR_d15_ind.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_R_NR_d15_ind.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_R_NR_d15_speed

these_speeds = d_info.speed_sess; % d_info.speed_sess or d_info.speed_AW
this_col = 'copper';
this_col_nan = [1,1,1]*2/3;
this_clim = [0,100];

spacer = 1;
F = paper_figure([0,0.5,mm2inch(3.25*34),mm2inch(0.85*34)]); hold on;
%     plt.cbar_numTicks = 7;
%     plt.clabel = 'Average running speed (cm/s)';
    
for i=1:5
    if i==1
        subplot(1,5+8+5+8+5+4*spacer,1:5); hold on;
    elseif i==2
        subplot(1,5+8+5+8+5+4*spacer,(6:13)+spacer); hold on;
    elseif i==3
        subplot(1,5+8+5+8+5+4*spacer,(14:18)+2*spacer); hold on;
    elseif i==4
        subplot(1,5+8+5+8+5+4*spacer,(19:26)+3*spacer); hold on;
    elseif i==5
        subplot(1,5+8+5+8+5+4*spacer,(27:31)+4*spacer); hold on;
    end
    yline(50,'k:');
    if i==1
        for k=1:d_info.numAnimals
            if ~isnan(these_speeds(k,i))
                temp = flipud(eval('copper'));
                [~,temp2] = discretize(this_clim,size(temp,1));
                plot(1:size(this_data_nostim_exp(k,:),2),this_data_nostim_exp(k,:)*100,'Color',temp(discretize(these_speeds(k,i),temp2),:),'LineWidth',0.5)
            else
                plot(1:size(this_data_nostim_exp(k,:),2),this_data_nostim_exp(k,:)*100,'Color',this_col_nan,'LineWidth',0.5);
            end
        end
        xlabel('')
        ylabel('Performance')
    elseif i==2
        for k=1:d_info.numAnimals
            if ~isnan(these_speeds(k,i))
                temp = flipud(eval('copper'));
                [~,temp2] = discretize(this_clim,size(temp,1));
                plot(1:size(this_data_nostim_sep(k,1:40),2),this_data_nostim_sep(k,1:40)*100,'Color',temp(discretize(these_speeds(k,i),temp2),:),'LineWidth',0.5)
            else
                plot(1:size(this_data_nostim_sep(k,1:40),2),this_data_nostim_sep(k,1:40)*100,'Color',this_col_nan,'LineWidth',0.5);
            end
        end
        xlabel('')
    elseif i==3
        for k=1:d_info.numAnimals
            if ~isnan(these_speeds(k,i))
                temp = flipud(eval('copper'));
                [~,temp2] = discretize(this_clim,size(temp,1));
                plot(1:size(this_data_nostim_sep(k,41:65),2),this_data_nostim_sep(k,41:65)*100,'Color',temp(discretize(these_speeds(k,i),temp2),:),'LineWidth',0.5)
            else
                plot(1:size(this_data_nostim_sep(k,41:65),2),this_data_nostim_sep(k,41:65)*100,'Color',this_col_nan,'LineWidth',0.5);
            end
        end
        xlabel('Trial')
    elseif i==4
        for k=1:d_info.numAnimals
            if ~isnan(these_speeds(k,i))
                temp = flipud(eval('copper'));
                [~,temp2] = discretize(this_clim,size(temp,1));
                plot(1:size(this_data_nostim_sep(k,66:105),2),this_data_nostim_sep(k,66:105)*100,'Color',temp(discretize(these_speeds(k,i),temp2),:),'LineWidth',0.5)
            else
                plot(1:size(this_data_nostim_sep(k,66:105),2),this_data_nostim_sep(k,66:105)*100,'Color',this_col_nan,'LineWidth',0.5);
            end
        end
        xlabel('')
    elseif i==5
        for k=1:d_info.numAnimals
            if ~isnan(these_speeds(k,i))
                temp = flipud(eval('copper'));
                [~,temp2] = discretize(this_clim,size(temp,1));
                plot(1:size(this_data_nostim_sep(k,106:130),2),this_data_nostim_sep(k,106:130)*100,'Color',temp(discretize(these_speeds(k,i),temp2),:),'LineWidth',0.5)
            else
                plot(1:size(this_data_nostim_sep(k,106:130),2),this_data_nostim_sep(k,106:130)*100,'Color',this_col_nan,'LineWidth',0.5);
            end
        end
        xlabel('')
    end
    ylim([35,100])
    yticks([50,100])
    ytickformat('percentage')
    if i~=1
        yticklabels({'',''})
    end
    if i==1 || i==3 || i==5
        xlim([0,25])
        xticks([0:5:25])
        xticklabels({'0','','','','','500'})
    elseif i==2 || i==4
        xlim([0,40])
        xticks([0:5:40])
        xticklabels({'0','','','','','','','','800'})
    end
end

savefig(F,[save_root_fig,'\Fig1_LearningCurves_R_NR_d15_speed.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_R_NR_d15_speed.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_R_NR_d15_speed.pdf']); set(gcf,'Color',[1,1,1])

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
imagesc([0,1,100,NaN],this_clim)
colormap(flipud(eval('copper')))
h=colorbar;
savefig(F,[save_root_fig,'\Fig1_LearningCurves_R_NR_d15_speed_colorbar.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_R_NR_d15_speed_colorbar.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_R_NR_d15_speed_colorbar.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_l100t % mean+-std numbers

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

% all no-stim mice
[nanmean(correct_100t(:,5)),nanstd(correct_100t(:,5))]
temp = correct_100t(~(d_info.group==7 | d_info.group==8),:);
[nanmean(temp(:,5+8)),nanstd(temp(:,5+8))]
[nanmean(temp(:,5+8+5)),nanstd(temp(:,5+8+5))]
[nanmean(temp(:,5+8+5+8)),nanstd(temp(:,5+8+5+8))]
[nanmean(temp(:,5+8+5+8+5)),nanstd(temp(:,5+8+5+8+5))]

% only runners
temp_R = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,2)==1,:);
[nanmean(temp_R(:,5+8)),nanstd(temp_R(:,5+8))]
temp_R = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,3)==1,:);
[nanmean(temp_R(:,5+8+5)),nanstd(temp_R(:,5+8+5))]
temp_R = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,4)==1,:);
[nanmean(temp_R(:,5+8+5+8)),nanstd(temp_R(:,5+8+5+8))]
temp_R = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,5)==1,:);
[nanmean(temp_R(:,5+8+5+8+5)),nanstd(temp_R(:,5+8+5+8+5))]

% only non-runners
temp_NR = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,2)==0,:);
[nanmean(temp_NR(:,5+8)),nanstd(temp_NR(:,5+8))]
temp_NR = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,3)==0,:);
[nanmean(temp_NR(:,5+8+5)),nanstd(temp_NR(:,5+8+5))]
temp_NR = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,4)==0,:);
[nanmean(temp_NR(:,5+8+5+8)),nanstd(temp_NR(:,5+8+5+8))]
temp_NR = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,5)==0,:);
[nanmean(temp_NR(:,5+8+5+8+5)),nanstd(temp_NR(:,5+8+5+8+5))]

% p-values: d2=0.79, d2=0.13, ...
ranksum(temp_R(:,5+8),temp_NR(:,5+8))
ranksum(temp_R(:,5+8+5),temp_NR(:,5+8+5))

temp_R = nan(d_info.numAnimals,2);
temp_R((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,2)==1,1) = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,2)==1,5+8);
temp_R((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,4)==1,2) = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,4)==1,5+8+5+8);
temp_R = nanmean(temp_R,2);
[nanmean(temp_R),nanstd(temp_R)]
temp_NR = nan(d_info.numAnimals,2);
temp_NR((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,2)==0,1) = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,2)==0,5+8);
temp_NR((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,4)==0,2) = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,4)==0,5+8+5+8);
temp_NR = nanmean(temp_NR,2);
[nanmean(temp_NR),nanstd(temp_NR)]
ranksum(temp_NR,temp_R)

temp_R = nan(d_info.numAnimals,2);
temp_R((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,3)==1,1) = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,3)==1,5+8+5);
temp_R((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,5)==1,2) = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,5)==1,5+8+5+8+5);
temp_R = nanmean(temp_R,2);
[nanmean(temp_R),nanstd(temp_R)]
temp_NR = nan(d_info.numAnimals,2);
temp_NR((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,3)==0,1) = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,3)==0,5+8+5);
temp_NR((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,5)==0,2) = correct_100t((~(d_info.group==7 | d_info.group==8)) & d_info.running(:,5)==0,5+8+5+8+5);
temp_NR = nanmean(temp_NR,2);
[nanmean(temp_NR),nanstd(temp_NR)]
ranksum(temp_NR,temp_R)


%% Fig1_LearningCurves_speedHistogram

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

h=histogram(d_info.speed_sess(:,1),0:10:100,'Normalization','count','FaceColor',p.col.darkGray,'FaceAlpha',1,'EdgeColor','none');
xline(10,':','Color',p.col.black,'LineWidth',1,'Alpha',1);

xlim([0,100])
xticks([0,50,100])
xlabel('Running speed (cm/s)')
ylim([0,25])
yticks([0,25])
ylabel('Number of mice')

savefig(F,[save_root_fig,'\Fig1_LearningCurves_speedHistogram.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_speedHistogram.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_speedHistogram.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_switchCorrelation

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

temp = correct_general(~(d_info.group==7 | d_info.group==8),:); %correct_100t(~(d_info.group==7 | d_info.group==8),:);
this_data_1 = temp(:,3);
this_data_2 = temp(:,5);

xline(50,'k:');
yline(50,'k:');
scatter(this_data_1*100,this_data_2*100,'.','SizeData',100,'MarkerEdgeColor',p.col.black);
[this_corr_r,this_corr_p] = fitLine(this_data_1*100,this_data_2*100,p.col.black);

xlim([35,100])
xticks([50,100])
xtickformat('percentage')
xlabel({'Performance on day 3'})
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel({'Performance on day 5'})


savefig(F,[save_root_fig,'\Fig1_LearningCurves_switchCorrelation.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_switchCorrelation.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_switchCorrelation.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_runningCorrelation

% these_labels = categorical({'NR -> NR','R -> R','NR -> R','R -> NR'});
% these_labels = reordercats(these_labels,{'NR -> NR','R -> R','NR -> R','R -> NR'});
these_labels = categorical({'Same','Different'});
these_labels = reordercats(these_labels,{'Same','Different'});
this_data = d_info.running(~(d_info.group==7 | d_info.group==8),:);

runningStatusChanges = nan(4,4);
for k=1:4
    temp = rmmissing(this_data(:,k:k+1));
    runningStatusChanges(1,k) = sum(temp(:,1)==0 & temp(:,2)==0) / length(temp);
    runningStatusChanges(2,k) = sum(temp(:,1)==1 & temp(:,2)==1) / length(temp);
    runningStatusChanges(3,k) = sum(temp(:,1)==0 & temp(:,2)==1) / length(temp);
    runningStatusChanges(4,k) = sum(temp(:,1)==1 & temp(:,2)==0) / length(temp);
end

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_plotting_data = [nanmean(runningStatusChanges(1:2,:),2)';nanmean([runningStatusChanges(4,:);runningStatusChanges(3,:)],2)']*100;
v = bar(these_labels,this_plotting_data,'stacked','DisplayName','runningStatusChanges');
% for i=1:4
%     plot(these_labels,[runningStatusChanges(i),runningStatusChanges(i),runningStatusChanges(i),runningStatusChanges(i)],'-k','LineWidth',1)
% end
v(1).FaceColor = 'flat'; v(2).FaceColor = 'flat'; v(1).EdgeColor = 'none'; v(2).EdgeColor = 'none'; 
v(1).CData(1,:) = p.col.nonrunner; v(1).CData(2,:) = p.col.nonrunner; v(2).CData(1,:) = p.col.runner; v(2).CData(2,:) = p.col.runner; 
ylim([0,100])
yticks([0,50,100])
ytickformat('percentage')
ylabel('Proportion of sessions')

savefig(F,[save_root_fig,'\Fig1_LearningCurves_runningCorrelation.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_runningCorrelation.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_runningCorrelation.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_speedCorrelation_d1

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_1 = d_info.speed_sess(:,1);
this_data_2 = correct_general(:,1);

xline(10,'k:');
yline(50,'k:');
scatter(this_data_1(this_data_1<10),this_data_2(this_data_1<10)*100,'.','SizeData',100,'MarkerEdgeColor',p.col.nonrunner);
scatter(this_data_1(this_data_1>10),this_data_2(this_data_1>10)*100,'.','SizeData',100,'MarkerEdgeColor',p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_1(this_data_1>10),this_data_2(this_data_1>10)*100,p.col.black);
%[this_corr_r,this_corr_p] = fitLine(this_data_1,this_data_2*100,p.col.black);

[a,b,c] = ranksum(this_data_2(this_data_1>10),this_data_2(this_data_1<10))

xlim([-10,100])
xticks([0,50,100])
xlabel({'Running speed (cm/s)'})
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel({'Performance on day 1'})

savefig(F,[save_root_fig,'\Fig1_LearningCurves_speedCorrelation_d1.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_speedCorrelation_d1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_speedCorrelation_d1.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_speedCorrelation_d35

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

temp = d_info.speed_sess(~(d_info.group==7 | d_info.group==8),[3,5]);
this_data_1 = temp(:);
temp = correct_general(~(d_info.group==7 | d_info.group==8),[3,5]);
this_data_2 = temp(:);

xline(10,'k:');
yline(50,'k:');
scatter(this_data_1(this_data_1<10),this_data_2(this_data_1<10)*100,'.','SizeData',100,'MarkerEdgeColor',p.col.nonrunner);
scatter(this_data_1(this_data_1>10),this_data_2(this_data_1>10)*100,'.','SizeData',100,'MarkerEdgeColor',p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_1(this_data_1>10),this_data_2(this_data_1>10)*100,p.col.black);
%[this_corr_r,this_corr_p] = fitLine(this_data_1,this_data_2*100,p.col.black);
[a,b,c]=ranksum(this_data_2(this_data_1<10),this_data_2(this_data_1>10))

xlim([-10,100])
xticks([0,50,100])
xlabel({'Running speed (cm/s)'})
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel({'Performance on 2nd switch days'})

savefig(F,[save_root_fig,'\Fig1_LearningCurves_speedCorrelation_d35.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_speedCorrelation_d35.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_speedCorrelation_d35.pdf']); set(gcf,'Color',[1,1,1])


%% Fig1_LearningCurves_speedCorrelation_d24

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

temp = d_info.speed_sess(~(d_info.group==7 | d_info.group==8),[2,4]);
this_data_1 = temp(:);
temp = correct_general(~(d_info.group==7 | d_info.group==8),[2,4]);
this_data_2 = temp(:);

xline(10,'k:');
yline(50,'k:');
scatter(this_data_1(this_data_1<10),this_data_2(this_data_1<10)*100,'.','SizeData',100,'MarkerEdgeColor',p.col.nonrunner);
scatter(this_data_1(this_data_1>10),this_data_2(this_data_1>10)*100,'.','SizeData',100,'MarkerEdgeColor',p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_1(this_data_1>10),this_data_2(this_data_1>10)*100,p.col.black);
%[this_corr_r,this_corr_p] = fitLine(this_data_1,this_data_2*100,p.col.black);

[a,b,c] = ranksum(this_data_2(this_data_1<10),this_data_2(this_data_1>10))

xlim([-10,100])
xticks([0,50,100])
xlabel({'Running speed (cm/s)'})
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel({'Performance on 1st switch days'})

savefig(F,[save_root_fig,'\Fig1_LearningCurves_speedCorrelation_d24.fig']);
saveas(F,[save_root_png,'\Fig1_LearningCurves_speedCorrelation_d24.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_LearningCurves_speedCorrelation_d24.pdf']); set(gcf,'Color',[1,1,1])




%% --- ExpertStim Quick Summary ---

clear; clc; close all;
disp(['Running ExpertStim_QuickSummary.m.'])


%% Load data

numTypes = 8;
numVars = 7;

% define root folders
if ispc
    path.root_repo                  = 'D:\SniffinHippo\Repo\';
    path.root_repoX                 = 'E:\SniffinHippo\RepoX\';
    path.root_analysis              = 'C:\SniffinHippo\Analysis\';
    path.root_summary               = 'C:\SniffinHippo\Summary\';
else
    path.root_repo                  = '/Users/Moritz/Documents/MATLAB/PhD/Repo/';
    path.root_repoX                 = '/Users/Moritz/Documents/MATLAB/PhD/RepoX/';
    path.root_analysis              = '/Users/Moritz/Documents/MATLAB/PhD/Analysis/';
    path.root_summary               = '/Users/Moritz/Documents/MATLAB/PhD/Summary/';
end

% load dataset sheet
path.file_in_Gsheet = [path.root_summary,'ExpertStim Pilot DATASET.xlsx'];
sheet = readtable(path.file_in_Gsheet,detectImportOptions(path.file_in_Gsheet,...
    'Sheet','Sheet1',...
    'NumHeaderLines',0,'ReadVariableNames',true,'DatetimeType','datetime'));
numSessions = nanmax(sheet.ID);
sheet = sheet(1:numSessions,:);

% load data
cmpr_list_repo = [];
cmpr_list_repo = [cmpr_list_repo,"perf_stim"];
d = {};
for i=1:numSessions
    for j=1:1
        this_session = sheet.ID(i);
        this_animal = sheet.Animal{i};
        this_date = [sheet.Date{i}(7:end),sheet.Date{i}(4:5),sheet.Date{i}(1:2)];
        path.folder_repo = [path.root_repo,this_animal,'/',this_animal,'_',this_date,'/'];
        path.filepart_in = [path.folder_repo,this_animal,'_',this_date,'_'];
        path.folder_analysis = [path.root_analysis,this_animal,'/',this_animal,'_',this_date,'/'];
        path.filepart_in_analysis = [path.folder_analysis,this_animal,'_',this_date,'_'];
        d{i,j} = repo2ws(path,cmpr_list_repo,0);
    end
end


%% Organise data in matrices - performance

performance_raw = nan(numTypes,numVars,numSessions);
for t=1:numTypes
    for v=1:numVars
        temp = extractVariable(d,['perf_stim.type',num2str(t),'.correct_var',num2str(v-1)],'cell','all');
        temp(cellfun('isempty',temp)) = {NaN};
        performance_raw(t,v,1:length(cell2mat(temp))) = cell2mat(temp);
    end
end

performance_tvs = nan(4,numVars,numSessions);
for t=1:4
    for v=1:numVars
        for s=1:numSessions
            
            if sheet.FirstOdours(s)
                t_source = t;
            else
                if t==1
                    t_source = 5;
                elseif t==2
                    t_source = 6;
                elseif t==3
                    t_source = 7;
                elseif t==4
                    t_source = 8;
                else
                    t_source = t;
                end
            end
            
            if sheet.Code(s)==101 || sheet.Code(s)==102
                if v==3
                    v_source = NaN;
                elseif v==4
                    v_source = 3;
                elseif v==5
                    v_source = 4;
                elseif v==6
                    v_source = NaN;
                elseif v==7
                    v_source = 5;
                else
                    v_source = v;
                end
            elseif sheet.Code(s)==103
                v_source = v;
            elseif sheet.Code(s)==104 || sheet.Code(s)==105
                if v==3
                    v_source = NaN;
                elseif v==4
                    v_source = NaN;
                elseif v==5
                    v_source = 3;
                elseif v==6
                    v_source = NaN;
                elseif v==7
                    v_source = NaN;
                else
                    v_source = v;
                end
            elseif sheet.Code(s)==106 || sheet.Code(s)==107
                if v==4
                    v_source = NaN;
                elseif v==5
                    v_source = 4;
                elseif v==6
                    v_source = 5;
                elseif v==7
                    v_source = NaN;
                else
                    v_source = v;
                end
            end
            
            if isnan(t_source) || isnan(v_source)
                performance_tvs(t,v,s) = NaN;
            else
                performance_tvs(t,v,s) = performance_raw(t_source,v_source,s);
            end
        end
    end
end

performance_vs = squeeze(nanmean(performance_tvs,1));
performance_sv = performance_vs';


%% Organise data in matrices - Hit rate

hitrate_raw = nan(numTypes,numVars,numSessions);
for t=1:numTypes
    for v=1:numVars
        temp = extractVariable(d,['perf_stim.type',num2str(t),'.H_var',num2str(v-1)],'cell','all');
        temp(cellfun('isempty',temp)) = {NaN};
        hitrate_raw(t,v,1:length(cell2mat(temp))) = cell2mat(temp);
    end
end

hitrate_tvs = nan(4,numVars,numSessions);
for t=1:4
    for v=1:numVars
        for s=1:numSessions
            
            if sheet.FirstOdours(s)
                t_source = t;
            else
                if t==1
                    t_source = 5;
                elseif t==2
                    t_source = 6;
                elseif t==3
                    t_source = 7;
                elseif t==4
                    t_source = 8;
                else
                    t_source = t;
                end
            end
            
            if sheet.Code(s)==101 || sheet.Code(s)==102
                if v==3
                    v_source = NaN;
                elseif v==4
                    v_source = 3;
                elseif v==5
                    v_source = 4;
                elseif v==6
                    v_source = NaN;
                elseif v==7
                    v_source = 5;
                else
                    v_source = v;
                end
            elseif sheet.Code(s)==103
                v_source = v;
            elseif sheet.Code(s)==104 || sheet.Code(s)==105
                if v==3
                    v_source = NaN;
                elseif v==4
                    v_source = NaN;
                elseif v==5
                    v_source = 3;
                elseif v==6
                    v_source = NaN;
                elseif v==7
                    v_source = NaN;
                else
                    v_source = v;
                end
            elseif sheet.Code(s)==106 || sheet.Code(s)==107
                if v==4
                    v_source = NaN;
                elseif v==5
                    v_source = 4;
                elseif v==6
                    v_source = 5;
                elseif v==7
                    v_source = NaN;
                else
                    v_source = v;
                end
            end
            
            if isnan(t_source) || isnan(v_source)
                hitrate_tvs(t,v,s) = NaN;
            else
                hitrate_tvs(t,v,s) = hitrate_raw(t_source,v_source,s);
            end
        end
    end
end

hitrate_vs = squeeze(nanmean(hitrate_tvs,1));
hitrate_sv = hitrate_vs';


%% Organise data in matrices - Correct rejection rate

crrate_raw = nan(numTypes,numVars,numSessions);
for t=1:numTypes
    for v=1:numVars
        temp = extractVariable(d,['perf_stim.type',num2str(t),'.CR_var',num2str(v-1)],'cell','all');
        temp(cellfun('isempty',temp)) = {NaN};
        crrate_raw(t,v,1:length(cell2mat(temp))) = cell2mat(temp);
    end
end

crrate_tvs = nan(4,numVars,numSessions);
for t=1:4
    for v=1:numVars
        for s=1:numSessions
            
            if sheet.FirstOdours(s)
                t_source = t;
            else
                if t==1
                    t_source = 5;
                elseif t==2
                    t_source = 6;
                elseif t==3
                    t_source = 7;
                elseif t==4
                    t_source = 8;
                else
                    t_source = t;
                end
            end
            
            if sheet.Code(s)==101 || sheet.Code(s)==102
                if v==3
                    v_source = NaN;
                elseif v==4
                    v_source = 3;
                elseif v==5
                    v_source = 4;
                elseif v==6
                    v_source = NaN;
                elseif v==7
                    v_source = 5;
                else
                    v_source = v;
                end
            elseif sheet.Code(s)==103
                v_source = v;
            elseif sheet.Code(s)==104 || sheet.Code(s)==105
                if v==3
                    v_source = NaN;
                elseif v==4
                    v_source = NaN;
                elseif v==5
                    v_source = 3;
                elseif v==6
                    v_source = NaN;
                elseif v==7
                    v_source = NaN;
                else
                    v_source = v;
                end
            elseif sheet.Code(s)==106 || sheet.Code(s)==107
                if v==4
                    v_source = NaN;
                elseif v==5
                    v_source = 4;
                elseif v==6
                    v_source = 5;
                elseif v==7
                    v_source = NaN;
                else
                    v_source = v;
                end
            end
            
            if isnan(t_source) || isnan(v_source)
                crrate_tvs(t,v,s) = NaN;
            else
                crrate_tvs(t,v,s) = crrate_raw(t_source,v_source,s);
            end
        end
    end
end

crrate_vs = squeeze(nanmean(crrate_tvs,1));
crrate_sv = crrate_vs';


%% --- Meeting with Michael and Arnd ---


%% 1) All expert stim experiments - Performance

default_figure; close; figure; 
subplot(1,3,1);
yline(50,':');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = performance_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        elseif sheet.NumPatterns(s)==21
            this_color = 'c';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if false % (s==1 || s==2 || s==16)
        plot(these_vars(temp),this_data(temp),[this_color,':',this_marker])
    else
        plot(these_vars(temp),this_data(temp),[this_color,'-',this_marker])
    end
end

xlim([0,8])
ylim([30,100])
ytickformat('percentage')
xticks(1:numVars)
xticklabels({'no-stim','ipsi-seq','ipsi-ctrl-rel','ipsi-ctrl-rnd','contra-seq','contra-ctrl-rel','contra-ctrl-rnd'})
xtickangle(45)
xlabel('Stim condition')
ylabel('Performance')
title('All expert stim experiments')


%%% 2) With first odour expert stim experiments - Performance

%default_figure; close; figure;
subplot(1,3,2);
yline(50,':');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = performance_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if sheet.FirstOdours(s)
        if false % (s==1 || s==2 || s==16)
            plot(these_vars(temp),this_data(temp),[this_color,':',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        else
            plot(these_vars(temp),this_data(temp),[this_color,'-',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        end
    end
end
% for v=1:numVars
%     this_data = performance_sv(setdiff(find(sheet.FirstOdours),[1,2,16]),v)*100;
%     plot(v,nanmean(this_data,1),'mo','LineWidth',2)
%     errorbar(v,nanmean(this_data,1),nansem(this_data,1),'m:','LineWidth',2)
% end
for v=1:numVars
    this_data = performance_sv(sheet.FirstOdours,v)*100;
    if v==3 || v==4 || v==6 || v==7
        plot(v,nanmean(this_data,1),'ko','LineWidth',0.5)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',0.5)
    else
        plot(v,nanmean(this_data,1),'ko','LineWidth',2)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',2)
    end
end

xlim([0,8])
ylim([30,100])
ytickformat('percentage')
xticks(1:numVars)
xticklabels({'no-stim','ipsi-seq','ipsi-ctrl-rel','ipsi-ctrl-rnd','contra-seq','contra-ctrl-rel','contra-ctrl-rnd'})
xtickangle(45)
xlabel('Stim condition')
ylabel('Performance')
title('Expert stim experiments with first odour')


%%% 3) Without first odour expert stim experiments - Performance

%default_figure; close; figure;
subplot(1,3,3);
yline(50,':');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = performance_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if ~sheet.FirstOdours(s)
        if false % (s==1 || s==2 || s==16)
            plot(these_vars(temp),this_data(temp),[this_color,':',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        else
            plot(these_vars(temp),this_data(temp),[this_color,'-',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        end
    end
end
% for v=1:numVars
%     this_data = performance_sv(setdiff(find(~sheet.FirstOdours),[1,2,16]),v)*100;
%     plot(v,nanmean(this_data,1),'mo','LineWidth',2)
%     errorbar(v,nanmean(this_data,1),nansem(this_data,1),'m:','LineWidth',2)
% end
for v=1:numVars
    this_data = performance_sv(~sheet.FirstOdours,v)*100;
    if v==3 || v==4 || v==6 || v==7
        plot(v,nanmean(this_data,1),'ko','LineWidth',0.5)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',0.5)
    else
        plot(v,nanmean(this_data,1),'ko','LineWidth',2)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',2)
    end
end

xlim([0,8])
ylim([30,100])
ytickformat('percentage')
xticks(1:numVars)
xticklabels({'no-stim','ipsi-seq','ipsi-ctrl-rel','ipsi-ctrl-rnd','contra-seq','contra-ctrl-rel','contra-ctrl-rnd'})
xtickangle(45)
xlabel('Stim condition')
ylabel('Performance')
title('Expert stim experiments without first odour')


%% 4) All expert stim experiments - Hit rate

default_figure; close; figure;
subplot(2,3,1);
yline(50,':');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = hitrate_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if false % (s==1 || s==2 || s==16)
        plot(these_vars(temp),this_data(temp),[this_color,':',this_marker])
    else
        plot(these_vars(temp),this_data(temp),[this_color,'-',this_marker])
    end
end

xlim([0,8])
ylim([0,100])
ytickformat('percentage')
xticks(1:numVars)
xticklabels({'no-stim','ipsi-seq','ipsi-ctrl-rel','ipsi-ctrl-rnd','contra-seq','contra-ctrl-rel','contra-ctrl-rnd'})
xtickangle(45)
xlabel('Stim condition')
ylabel('Hit rate')
title('All expert stim experiments')


%%% 5) With first odour expert stim experiments - Hit rate

%default_figure; close; figure;
subplot(2,3,2);
yline(50,':');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = hitrate_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if sheet.FirstOdours(s)
        if false % (s==1 || s==2 || s==16)
            plot(these_vars(temp),this_data(temp),[this_color,':',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        else
            plot(these_vars(temp),this_data(temp),[this_color,'-',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        end
    end
end
% for v=1:numVars
%     this_data = hitrate_sv(setdiff(find(sheet.FirstOdours),[1,2,16]),v)*100;
%     plot(v,nanmean(this_data,1),'mo','LineWidth',2)
%     errorbar(v,nanmean(this_data,1),nansem(this_data,1),'m:','LineWidth',2)
% end
for v=1:numVars
    this_data = hitrate_sv(sheet.FirstOdours,v)*100;
    if v==3 || v==4 || v==6 || v==7
        plot(v,nanmean(this_data,1),'ko','LineWidth',0.5)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',0.5)
    else
        plot(v,nanmean(this_data,1),'ko','LineWidth',2)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',2)
    end
end

xlim([0,8])
ylim([0,100])
ytickformat('percentage')
xticks(1:numVars)
xticklabels({'no-stim','ipsi-seq','ipsi-ctrl-rel','ipsi-ctrl-rnd','contra-seq','contra-ctrl-rel','contra-ctrl-rnd'})
xtickangle(45)
xlabel('Stim condition')
ylabel('Hit rate')
title('Expert stim experiments with first odour')


%%% 6) Without first odour expert stim experiments - Hit rate

%default_figure; close; figure;
subplot(2,3,3);
yline(50,':');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = hitrate_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if ~sheet.FirstOdours(s)
        if false % (s==1 || s==2 || s==16)
            plot(these_vars(temp),this_data(temp),[this_color,':',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        else
            plot(these_vars(temp),this_data(temp),[this_color,'-',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        end
    end
end
% for v=1:numVars
%     this_data = hitrate_sv(setdiff(find(~sheet.FirstOdours),[1,2,16]),v)*100;
%     plot(v,nanmean(this_data,1),'mo','LineWidth',2)
%     errorbar(v,nanmean(this_data,1),nansem(this_data,1),'m:','LineWidth',2)
% end
for v=1:numVars
    this_data = hitrate_sv(~sheet.FirstOdours,v)*100;
    if v==3 || v==4 || v==6 || v==7
        plot(v,nanmean(this_data,1),'ko','LineWidth',0.5)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',0.5)
    else
        plot(v,nanmean(this_data,1),'ko','LineWidth',2)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',2)
    end
end

xlim([0,8])
ylim([0,100])
ytickformat('percentage')
xticks(1:numVars)
xticklabels({'no-stim','ipsi-seq','ipsi-ctrl-rel','ipsi-ctrl-rnd','contra-seq','contra-ctrl-rel','contra-ctrl-rnd'})
xtickangle(45)
xlabel('Stim condition')
ylabel('Hit rate')
title('Expert stim experiments without first odour')


%%% 7) All expert stim experiments - Correct rejection rate

%default_figure; close; figure;
subplot(2,3,4);
yline(50,':');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = crrate_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if false % (s==1 || s==2 || s==16)
        plot(these_vars(temp),this_data(temp),[this_color,':',this_marker])
    else
        plot(these_vars(temp),this_data(temp),[this_color,'-',this_marker])
    end
end

xlim([0,8])
ylim([0,100])
ytickformat('percentage')
xticks(1:numVars)
xticklabels({'no-stim','ipsi-seq','ipsi-ctrl-rel','ipsi-ctrl-rnd','contra-seq','contra-ctrl-rel','contra-ctrl-rnd'})
xtickangle(45)
xlabel('Stim condition')
ylabel('Correct rejection rate')
title('All expert stim experiments')


%%% 8) With first odour expert stim experiments - Correct rejection rate

%default_figure; close; figure;
subplot(2,3,5);
yline(50,':');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = crrate_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if sheet.FirstOdours(s)
        if false % (s==1 || s==2 || s==16)
            plot(these_vars(temp),this_data(temp),[this_color,':',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        else
            plot(these_vars(temp),this_data(temp),[this_color,'-',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        end
    end
end
% for v=1:numVars
%     this_data = crrate_sv(setdiff(find(sheet.FirstOdours),[1,2,16]),v)*100;
%     plot(v,nanmean(this_data,1),'mo','LineWidth',2)
%     errorbar(v,nanmean(this_data,1),nansem(this_data,1),'m:','LineWidth',2)
% end
for v=1:numVars
    this_data = crrate_sv(sheet.FirstOdours,v)*100;
    if v==3 || v==4 || v==6 || v==7
        plot(v,nanmean(this_data,1),'ko','LineWidth',0.5)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',0.5)
    else
        plot(v,nanmean(this_data,1),'ko','LineWidth',2)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',2)
    end
end

xlim([0,8])
ylim([0,100])
ytickformat('percentage')
xticks(1:numVars)
xticklabels({'no-stim','ipsi-seq','ipsi-ctrl-rel','ipsi-ctrl-rnd','contra-seq','contra-ctrl-rel','contra-ctrl-rnd'})
xtickangle(45)
xlabel('Stim condition')
ylabel('Correct rejection rate')
title('Expert stim experiments with first odour')


%%% 9) Without first odour expert stim experiments - Correct rejection rate

%default_figure; close; figure;
subplot(2,3,6);
yline(50,':');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = crrate_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if ~sheet.FirstOdours(s)
        if false % (s==1 || s==2 || s==16)
            plot(these_vars(temp),this_data(temp),[this_color,':',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        else
            plot(these_vars(temp),this_data(temp),[this_color,'-',this_marker],'Color',nanmean([1,1,1;0.5,0.5,0.5],1))
        end
    end
end
% for v=1:numVars
%     this_data = crrate_sv(setdiff(find(~sheet.FirstOdours),[1,2,16]),v)*100;
%     plot(v,nanmean(this_data,1),'mo','LineWidth',2)
%     errorbar(v,nanmean(this_data,1),nansem(this_data,1),'m:','LineWidth',2)
% end
for v=1:numVars
    this_data = crrate_sv(~sheet.FirstOdours,v)*100;
    if v==3 || v==4 || v==6 || v==7
        plot(v,nanmean(this_data,1),'ko','LineWidth',0.5)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',0.5)
    else
        plot(v,nanmean(this_data,1),'ko','LineWidth',2)
        errorbar(v,nanmean(this_data,1),nansem(this_data,1),'k-','LineWidth',2)
    end
end

xlim([0,8])
ylim([0,100])
ytickformat('percentage')
xticks(1:numVars)
xticklabels({'no-stim','ipsi-seq','ipsi-ctrl-rel','ipsi-ctrl-rnd','contra-seq','contra-ctrl-rel','contra-ctrl-rnd'})
xtickangle(45)
xlabel('Stim condition')
ylabel('Correct rejection rate')
title('Expert stim experiments without first odour')


%% Baseline-relative performance change (with first odour sessions)

p = get_p;

default_figure; close; figure;
subplot(1,3,1);
plot(linspace(0,100),linspace(0,100),'k:');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = performance_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if sheet.FirstOdours(s)
        line([this_data(1),this_data(1)],[this_data(2),this_data(5)],'Color',p.col.darkGray);
        scatter(this_data(1),this_data(2),'MarkerFaceColor',p.col.seq,'MarkerEdgeColor',p.col.seq)
        scatter(this_data(1),this_data(5),'MarkerFaceColor',p.col.ctrl,'MarkerEdgeColor',p.col.ctrl)
    end
end

xlim([0,100])
ylim([0,100])
pbaspect([1,1,1])
xtickformat('percentage')
ytickformat('percentage')
xlabel('Performance in no-stim trials')
ylabel('Performance in stim trials')
title('Expert stim experiments with first odour')


subplot(1,3,2);
plot(linspace(0,100),linspace(0,100),'k:');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = hitrate_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if sheet.FirstOdours(s)
        line([this_data(1),this_data(1)],[this_data(2),this_data(5)],'Color',p.col.darkGray);
        scatter(this_data(1),this_data(2),'MarkerFaceColor',p.col.seq,'MarkerEdgeColor',p.col.seq)
        scatter(this_data(1),this_data(5),'MarkerFaceColor',p.col.ctrl,'MarkerEdgeColor',p.col.ctrl)
    end
end

xlim([0,100])
ylim([0,100])
pbaspect([1,1,1])
xtickformat('percentage')
ytickformat('percentage')
xlabel('Hit rate in no-stim trials')
ylabel('Hit rate in stim trials')
title('Expert stim experiments with first odour')


subplot(1,3,3);
plot(linspace(0,100),linspace(0,100),'k:');
hold on

for s=1:numSessions
    these_vars = 1:numVars;
    this_data = crrate_sv(s,:)*100;

    if sheet.FirstOdours(s)
        if sheet.NumPatterns(s)==20
            this_color = 'k';
        elseif sheet.NumPatterns(s)==7
            this_color = 'b';
        end
    else
        if sheet.NumPatterns(s)==20
            this_color = 'r';
        elseif sheet.NumPatterns(s)==7
            this_color = 'g';
        elseif sheet.NumPatterns(s)==11
            this_color = 'y';
        end
    end
    if sheet.NumPatterns(s)==20
        this_marker = 's';
    elseif sheet.NumPatterns(s)==7
        this_marker = 'd';
    elseif sheet.NumPatterns(s)==11
        this_marker = 'o';
    else
        this_marker = 'x';
    end

    temp = ~any(isnan(this_data),1);
    if sheet.FirstOdours(s)
        line([this_data(1),this_data(1)],[this_data(2),this_data(5)],'Color',p.col.darkGray);
        scatter(this_data(1),this_data(2),'MarkerFaceColor',p.col.seq,'MarkerEdgeColor',p.col.seq)
        scatter(this_data(1),this_data(5),'MarkerFaceColor',p.col.ctrl,'MarkerEdgeColor',p.col.ctrl)
    end
end

xlim([0,100])
ylim([0,100])
pbaspect([1,1,1])
xtickformat('percentage')
ytickformat('percentage')
xlabel('Correct rejection rate in no-stim trials')
ylabel('Correct rejection rate in stim trials')
title('Expert stim experiments with first odour')



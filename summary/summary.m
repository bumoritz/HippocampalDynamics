function [out] = summary(path,ops)

%% Preparations

% gather p
p = get_p;

% identify session folders
path.file_in_Gsheet = [path.root_summary,'defaults\SniffinHippo SeqStim DATASET.xlsx'];
if ~exist([path.root_summary,'log'],'dir')
    mkdir([path.root_summary,'log']);
end
if ~exist([path.root_summary,'plots'],'dir')
    mkdir([path.root_summary,'plots']);
end

% start log
t_start = tic;
log.executionDate = datetime(now,'ConvertFrom','datenum');
log.done = false;
temp = [path.root_summary,'log/',datestr(log.executionDate,'yyyymmdd_HHMMSS'),'.log'];
diary(temp);
diary on;

disp([log.executionDate])
disp(['Running summary.m.'])


%% Load dataset Google sheet

if ops.loadDefault_sheet
    load([path.root_summary,'defaults\sheet.mat']);
else
    sheet = readtable(path.file_in_Gsheet,detectImportOptions(path.file_in_Gsheet,...
        'Sheet','ANIMALS',...
        'NumHeaderLines',1,'ReadVariableNames',true,'DatetimeType','datetime'));
end
d_info.numAnimals = nanmax(sheet.ID);
d_info.numDays = 5;
sheet = sheet(1:d_info.numAnimals,:);
d_info.group = sheet.ID_1;
d_info.animals = sheet.Name;
d_info.dates_day1 = char();
for i=1:d_info.numAnimals
    try
        try
            d_info.dates_day1(i,1:8) = datestr(sheet.Day1{i},'yyyymmdd');
        catch
            d_info.dates_day1(i,1:8) = datestr(sheet.Day1(i),'yyyymmdd');
        end
    catch 
        d_info.dates_day1(i,1:8) = '00000000';
    end
end


%% Select dataset

%d_info = selectDataset(d_info,'-g0123456789101112-d1-e1-r01-p01-l01-i00',sheet,path,ops);
%d_info = selectDataset(d_info,'-g2-d12345-e01-r01-p01-l01-i00',sheet,path,ops);

%d_info = selectDataset(d_info,'-g78-d12345-e1-r01-p1-l01-i00',sheet,path,ops);

% --

%d_info = selectDataset(d_info,'-g12-d12345-e1-r1-p01-l01-i00',sheet,path,ops);

% all data
%d_info = selectDataset(d_info,'-g0123456789101112-d12345-e01-r01-p01-l01-i00',sheet,path,ops);

% all data - but only engaged and responding
%d_info = selectDataset(d_info,'-g0123456789101112-d12345-e1-r01-p01-l01-i00',sheet,path,ops);

% all data - but only responding for stim sessions
d_info = selectDataset(d_info,'-g012345678-d12345-e01-r01-p1-l01-i00',sheet,path,ops);

% all stim data
%d_info = selectDataset(d_info,'-g78-d2345-e01-r01-p01-l01-i00',sheet,path,ops);

% responsive stim data
d_info = selectDataset(d_info,'-g78-d2345-e01-r01-p1-l01-i00',sheet,path,ops);

% all switch data with imaging
%d_info = selectDataset(d_info,'-g278-d2345-e01-r01-p01-l01-i00',sheet,path,ops);

% all switch data
%d_info = selectDataset(d_info,'-g1278-d2345-e01-r01-p01-l01-i00',sheet,path,ops);

% all imaging during switches
d_info = selectDataset(d_info,'-g2-d12345-e01-r01-p01-l01-i00',sheet,path,ops);

% all switch data without stim
%d_info = selectDataset(d_info,'-g12-d12345-e01-r01-p01-l01-i00',sheet,path,ops);

% all day 1 imaging data
d_info = selectDataset(d_info,'-g02345678-d1-e01-r01-p01-l01-i00',sheet,path,ops);


%% Load data

disp('- Loading data.')

% identify relevant data
cmpr_list_repo = [];
if ops.do_learningCurveSummary
    cmpr_list_repo = [cmpr_list_repo,"perf"];
%     if ops.lcs.alternativeLickMetrics || ops.lcs.runningMetrics
%         cmpr_list_repo = [cmpr_list_repo,"paq_beh"];
%     end
end
if ops.do_lickingSummary
    cmpr_list_repo = [cmpr_list_repo,"task","paq_beh"];
end
if ops.do_respXperfSummary
	cmpr_list_repo = [cmpr_list_repo,"resp","perf"];
end
if ops.do_sequenceCellSummary
    cmpr_list_repo = [cmpr_list_repo,"tng_all","warp_all","pca_all","paq_beh","perf","bcon","task","dec_all"];
    %cmpr_list_repo = [cmpr_list_repo,"tngn_all","warpn_all","paq_beh","perf"];
end
if ops.do_sequenceDevelopmentSummary
    cmpr_list_repo = [cmpr_list_repo,"tng_100t","warp_100t","bcon","perf","dec_100T"];
end
if ops.do_sequenceDevelopmentStimSummary
    %cmpr_list_repo = [cmpr_list_repo,"warp_all_stimVersion","tng_all_stimVersion","resp","perf"];
    cmpr_list_repo = [cmpr_list_repo,"warp_all_stimVersion","tng_all_stimVersion","warp_100t_stimVersion","tng_100t_stimVersion","ppa","trg","resp","perf"];
    %cmpr_list_repo = [cmpr_list_repo,"warp_100t_stimVersion","tng_100t_stimVersion","perf"];
    %cmpr_list_repo = [cmpr_list_repo,"warpn_100t_stimVersion","tngn_100t_stimVersion","perf"];
end
if ops.do_tngEncXperfSummary
	cmpr_list_repo = [cmpr_list_repo,"tng_all","perf","paq_beh"];
end
if ops.do_sqnXperfSummary
	cmpr_list_repo = [cmpr_list_repo,"sqn_all","task"];
end
if ops.do_ecaXperfSummary
	cmpr_list_repo = [cmpr_list_repo,"eca_all"];
end
if ops.do_nemSummary
	cmpr_list_repo = [cmpr_list_repo,"nem_all_cmpr","nem_100t_cmpr","paq_beh"]; %[cmpr_list_repo,"nem_all","perf"]; %[cmpr_list_repo,"nem_all"];
end
if ops.do_imprSummary
	cmpr_list_repo = [cmpr_list_repo,"tng_all_stimVersion","resp","flw","trg","perf","str"];
    %cmpr_list_repo = [cmpr_list_repo,"impr_all","impr_100t"]; %[cmpr_list_repo,"nem_all","perf"]; %[cmpr_list_repo,"nem_all"];
end
if ops.do_responseSummary
    cmpr_list_repo = [cmpr_list_repo,"resp","flw","flwsta"];
end
if ops.do_flexibleSummary
    cmpr_list_repo = [cmpr_list_repo,"s2p_meta","paq_beh"];
end
cmpr_list_repo = unique(cmpr_list_repo);

% load data
d = {};
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        this_animal = d_info.animals{i};
        if (~isempty(sheet.(['Day',num2str(j)])(i))) && (~d_info.excl(i,j))  %%%%%%%%%%%%%%%%%%%(~isempty(sheet.(['Day',num2str(j)]){i})) && (~d_info.excl(i,j))
            try
                this_date = datestr(sheet.(['Day',num2str(j)])(i),'yyyymmdd'); %%%%%%%%%%%%%%datestr(sheet.(['Day',num2str(j)]){i},'yyyymmdd');

                % load repo data
                path.folder_repo = [path.root_repo,this_animal,'/',this_animal,'_',this_date,'/'];
                path.filepart_in = [path.folder_repo,this_animal,'_',this_date,'_'];
                path.folder_analysis = [path.root_analysis,this_animal,'/',this_animal,'_',this_date,'/'];
                path.filepart_in_analysis = [path.folder_analysis,this_animal,'_',this_date,'_'];
                path.folder_analysisX = [path.root_analysisX,this_animal,'/',this_animal,'_',this_date,'/'];
                path.filepart_in_analysisX = [path.folder_analysisX,this_animal,'_',this_date,'_'];
                
                this_cmpr_list_repo = cmpr_list_repo;
                if any(cmpr_list_repo=="trg")
                    try
                        load([path.filepart_in,'meta.mat']);
                        if strcmp(meta.reg.trg,"trg_rigid")
                            this_cmpr_list_repo(find(cmpr_list_repo=="trg")) = "trg_rigid";
                        elseif strcmp(meta.reg.trg,"trg_nonrigid")
                            this_cmpr_list_repo(find(cmpr_list_repo=="trg")) = "trg_nonrigid";
                        else
                            warning('Something went wrong')
                        end
                    catch
                    end
                end
                d{i,j} = repo2ws(path,this_cmpr_list_repo,ops.skipIncompletelyProcessed);

                %analyses(this_animal,this_date,path,ops);
            catch
            end
        else
            d{i,j} = {};
        end
    end
end


%% Do summary analyses

out = struct();
out.ops = ops;
out.p = p;
out.d_info = d_info;

if ops.do_learningCurveSummary
    disp('- [lcs] module')
    [lcs] = learningCurveSummary(d_info,d,ops,p,path);
    out.lcs = lcs;
end

if ops.do_lickingSummary
    disp('- [licks] module')
    [licks] = lickingSummary(d_info,d,ops,p,path);
    out.licks = licks;
end

if ops.do_respXperfSummary
    disp('- [respXperf] module')
    respXperfSummary(d_info,d,ops,p,path);
end

if ops.do_sequenceCellSummary
    disp('- [scs] module')
    sequenceCellSummary(d_info,d,ops,p,path);
end

if ops.do_tngEncXperfSummary
    disp('- [tngEncXperf] module')
    tngEncXperfSummary(d_info,d,ops,p,path);
end

if ops.do_sqnXperfSummary
    disp('- [sqnXperf] module')
    sqnXperfSummary(d_info,d,ops,p,path);
end

if ops.do_ecaXperfSummary
    disp('- [ecaXperf] module')
    ecaXperfSummary(d_info,d,ops,p,path);
end

if ops.do_nemSummary
    disp('- [nem] module')
    nemSummary(d_info,d,ops,p,path);
end

if ops.do_imprSummary
    disp('- [impr] module')
    imprSummary(d_info,d,ops,p,path);
end

out.d = d;


%% Complete execution

log.done = true;
log.runTime = toc(t_start);

disp(['- Done in ',num2str(log.runTime/60,3),' min.'])
diary off;

end


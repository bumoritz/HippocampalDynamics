function [path,task,perf] = data2repo_task(info,path,epoch)
% WHENEVER USING PERFORMANCE IN EXPERT STIM SESSION: ONLY USE SPECIFIC TYPE+VAR COMBINATIONS (not odourO or _stim fields)

if nargin < 3
    epoch = 'beh';
end

if strcmp(epoch,'beh')
    these_numTrials = info.task.numTrials;
    these_numBlocks = info.task.numBlocks;
    these_trialsPerBlock = info.task.trialsPerBlock;
else
    these_numTrials = info.task.(['numTrials_',epoch]);
    these_numBlocks = info.task.(['numBlocks_',epoch]);
    these_trialsPerBlock = info.task.(['trialsPerBlock_',epoch]);
end

disp('--- Loading and processing task data...')

%% Load data

% seq file
if strcmp(epoch,'beh')
    temp = dir([path.folder_data,'Behaviour\',info.animal,'_',info.date,'_*_SEQ.txt']);
    if size(temp,1)>1
        path.file_in_seq = [path.folder_data,'Behaviour\',temp(1).name];
    else
        path.file_in_seq = [path.folder_data,'Behaviour\',temp.name];
    end
else
    temp = dir([path.folder_data,'Behaviour\',epoch,'\',info.animal,'_',info.date,'_*_SEQ.txt']);
    if size(temp,1)>1
        path.file_in_seq = [path.folder_data,'Behaviour\',epoch,'\',temp(1).name];
    else
        path.file_in_seq = [path.folder_data,'Behaviour\',epoch,'\',temp.name];
    end
end
if isfile(path.file_in_seq)
    temp = readmatrix(path.file_in_seq);
    tmp.raw.seq.stim = temp(1,:);
    tmp.raw.seq.var = temp(2,:);
else
    error(['TrialSeq file at ',path.file_in_seq,' does not exist.'])
end
if (info.data.numFragments~=1) && (~info.data.fragments_taskContinuous)
    tmp.raw.seq.stim = tmp.raw.seq.stim(info.data.fragments_seq);
    tmp.raw.seq.var = tmp.raw.seq.var(info.data.fragments_seq);
end

% snh file
if strcmp(epoch,'beh')
    temp = dir([path.folder_data,'Behaviour\',info.animal,'_',info.date,'_*_SNH.txt']);
    if info.data.numFragments==1 || info.data.fragments_taskContinuous
        path.file_in_snh = [path.folder_data,'Behaviour\',temp.name];
    else
        path.file_in_snh = [path.folder_data,'Behaviour\',temp(1).name];
    end
else
    temp = dir([path.folder_data,'Behaviour\',epoch,'\',info.animal,'_',info.date,'_*_SNH.txt']);
    if info.data.numFragments==1 || info.data.fragments_taskContinuous
        path.file_in_snh = [path.folder_data,'Behaviour\',epoch,'\',temp.name];
    else
        path.file_in_snh = [path.folder_data,'Behaviour\',epoch,'\',temp(1).name];
    end
end
if isfile(path.file_in_snh)
    temp = fopen(path.file_in_snh);
    temp2 = textscan(temp,'%c');
    tmp.raw.snh.text = convertCharsToStrings(temp2{1,1});
    fclose(temp);
else
    error(['SniffinHippo file at ',path.file_in_snh,' does not exist.'])
end

% pyb file
if info.data.numFragments==1 || info.data.fragments_taskContinuous
    if strcmp(epoch,'beh')
        temp = dir([path.folder_data,'Behaviour\',info.animal,'_',info.date,'_*_PYB.mat']);
        path.file_in_pyb = [path.folder_data,'Behaviour\',temp.name];
    else
        temp = dir([path.folder_data,'Behaviour\',epoch,'\',info.animal,'_',info.date,'_*_PYB.mat']);
        path.file_in_pyb = [path.folder_data,'Behaviour\',epoch,'\',temp.name];
    end
    if isfile(path.file_in_pyb)
        tmp.raw.pyb = load(path.file_in_pyb);
    else
        error(['PyBehaviour mat file at ',path.file_in_pyb,' does not exist.'])
    end
else
    if ~strcmp(epoch,'beh')
        disp('Fragmented files not yet implemented for expert stim experiment.')
    end
    tmp.raw.pyb.results = {};
    for i=1:info.data.numFragments
        temp = dir([path.folder_data,'Behaviour\',info.animal,'_',info.date,'_*_PYB.mat']);
        if ~info.data.missingMatFile(i)
            path.(['file_in_pyb_',num2str(i)]) = [path.folder_data,'Behaviour\',temp(i).name];
            if isfile(path.(['file_in_pyb_',num2str(i)]))
                temp2 = load(path.(['file_in_pyb_',num2str(i)]));
            else
                error(['PyBehaviour mat file at ',path.(['file_in_pyb_',num2str(i)]),' does not exist.'])
            end
        else
            temp = dir([path.folder_data,'Behaviour\',info.animal,'_',info.date,'_*_PYB.pkl']);
            if length(temp)==info.data.numFragments
                warning('Mat file missing for at least one fragment. Creating it from pkl file.')
                path.(['file_in_pyb_',num2str(i)]) = [path.folder_data,'Behaviour\',temp(i).name];
                temp2 = pickle2mat(path.(['file_in_pyb_',num2str(i)]),0);
            else
                warning('Mat file missing for at least one fragment. No pkl file either. Check manually.')
            end
        end
        tmp.raw.pyb.results = [tmp.raw.pyb.results,temp2.results{1:info.data.fragments_numTrials(i)}];
    end
    tmp.raw.pyb.running_score = NaN; tmp.raw.pyb.reaction_time = NaN; tmp.raw.pyb.correct_tally = NaN;
end

if strcmp(epoch,'beh')
    path.file_out_task = [path.folder_repo,info.animal,'_',info.date,'_task.mat'];
    path.file_out_perf = [path.folder_repo,info.animal,'_',info.date,'_perf.mat'];
else
    path.file_out_task = [path.folder_repo,info.animal,'_',info.date,'_task_',epoch,'.mat'];
    path.file_out_perf = [path.folder_repo,info.animal,'_',info.date,'_perf_',epoch,'.mat'];
end


%% Sanity checks

% SANITY CHECK: trial structure
tmp.trialStructure.tOdour1 = str2num(extractBetween(tmp.raw.snh.text,strfind(tmp.raw.snh.text,'T_ODOUR1:')+9,strfind(tmp.raw.snh.text,'T_GAP:')-2))/1000; % [s]
tmp.trialStructure.tGap = str2num(extractBetween(tmp.raw.snh.text,strfind(tmp.raw.snh.text,'T_GAP:')+6,strfind(tmp.raw.snh.text,'T_ODOUR2:')-2))/1000; % [s]
tmp.trialStructure.tOdour2 = str2num(extractBetween(tmp.raw.snh.text,strfind(tmp.raw.snh.text,'T_ODOUR2:')+9,strfind(tmp.raw.snh.text,'TYPE_CONTS:')-2))/1000; % [s]
if (tmp.trialStructure.tOdour1~=info.task.trialStructure.tOdour1) | ...
        (tmp.trialStructure.tGap~=info.task.trialStructure.tGap) | ...
        (tmp.trialStructure.tOdour2~=info.task.trialStructure.tOdour2)
    warning('Mismatch in trial structure between SniffinHippo file and user input.')
end

% SANITY CHECK: vials
temp = extractBetween(tmp.raw.snh.text,strfind(tmp.raw.snh.text,'VIAL_CONTS:')+11,strfind(tmp.raw.snh.text,'NUM_TRIALS:')-2);
temp2 = strfind(temp,',');
tmp.vials.vialNumber = [];
tmp.vials.role = strings;
for j=0:length(temp2)
    if j==0
        tmp.vials.vialNumber(j+1) = str2num(extractBetween(temp,1,1));
        tmp.vials.role(j+1) = extractBetween(temp,2,2);
    else
        tmp.vials.vialNumber(j+1) = str2num(extractBetween(temp,temp2(j)+1,temp2(j)+1));
        tmp.vials.role(j+1) = extractBetween(temp,temp2(j)+2,temp2(j)+2);
    end
end
if (length(tmp.vials.vialNumber)~=length(info.task.vials.vialNumber)) | ...
        ~all(strcmp(tmp.vials.role,info.task.vials.role))
%     warning('Mismatch in vials between SniffinHippo file and user input.')
end

% SANITY CHECK: contingencies
temp = extractBetween(tmp.raw.snh.text,strfind(tmp.raw.snh.text,'TYPE_CONTS:')+11,strfind(tmp.raw.snh.text,'VIAL_CONTS:')-2);
tmp.contingencies.stim = [];
tmp.contingencies.var = [];
tmp.contingencies.odour1 = strings;
tmp.contingencies.odour2 = strings;
for j=1:ceil(strlength(temp)/2)
    if mod(j,5)==1
        tmp.contingencies.stim = [tmp.contingencies.stim; str2num(extractBetween(temp,j,j))];
    elseif mod(j,5)==2
        tmp.contingencies.var = [tmp.contingencies.var; str2num(extractBetween(temp,j,j))];
    elseif mod(j,5)==3
        tmp.contingencies.odour1 = [tmp.contingencies.odour1; extractBetween(temp,j,j)];
    elseif mod(j,5)==4
        tmp.contingencies.odour2 = [tmp.contingencies.odour2; extractBetween(temp,j,j)];
    end
end
tmp.contingencies.odour1 = tmp.contingencies.odour1(2:end);
tmp.contingencies.odour2 = tmp.contingencies.odour2(2:end);
temp = min(length(tmp.contingencies.stim),length(info.task.contingencies.type));
% if (~all(strcmp(tmp.contingencies.odour1(1:temp)',info.task.contingencies.odour1(1:temp))) | ...
%     ~all(strcmp(tmp.contingencies.odour2(1:temp)',info.task.contingencies.odour2(1:temp))))
%     warning('Mismatch in contingencies between SniffinHippo file and user input.')
% end

% SANITY CHECK: number of trials
if these_numTrials ~= size(tmp.raw.pyb.results,2)
    warning('Mismatch in number of trials between PyBehaviour file and user input.')
end


%% Create task struct

task.type = zeros(1,these_numTrials);
task.autoreward = zeros(1,these_numTrials);
task.var = zeros(1,these_numTrials);
task.requirement = strings(1,these_numTrials);
task.response = strings(1,these_numTrials);
task.odour1 = strings(1,these_numTrials);
task.odour2 = strings(1,these_numTrials);
for j=1:these_numTrials
    
    % stim
    task.type(j) = tmp.raw.seq.stim(j);
    
    % autoreward
    if tmp.raw.pyb.results{1,j}.stim_type==5 | tmp.raw.pyb.results{1,j}.stim_type==6
        task.autoreward(j) = 1;
    end
    
    % var
    task.var(j) = tmp.raw.seq.var(j); 

    % requirement
    if tmp.raw.pyb.results{1,j}.response_required
        task.requirement(j) = "GO";
    else
        task.requirement(j) = "NOGO";
    end
    
    % response
    if tmp.raw.pyb.results{1,j}.cheated
        task.response(j) = "CHEAT";
    elseif tmp.raw.pyb.results{1,j}.incorrect
        task.response(j) = "FA";
    elseif tmp.raw.pyb.results{1,j}.miss
        task.response(j) = "M";
    elseif (tmp.raw.pyb.results{1,j}.correct & tmp.raw.pyb.results{1,j}.response_required)
        task.response(j) = "H";
    elseif (tmp.raw.pyb.results{1,j}.correct & ~tmp.raw.pyb.results{1,j}.response_required)
        task.response(j) = "CR";
    else
        task.response(j) = "NaN";
    end
    
    % odours
    task.odour1(j) = info.task.contingencies.odour1(find(info.task.contingencies.type==task.type(j)));
    task.odour2(j) = info.task.contingencies.odour2(find(info.task.contingencies.type==task.type(j))); 
end


%% Sanity checks

% contingency requirements
% for j=1:length(info.task.contingencies.type)
%     if ~isempty(min(find(task.type==info.task.contingencies.type(j))))
%         tmp.contingencies.requirement(j) = task.requirement(min(find(task.type==info.task.contingencies.type(j))));
%     else
%         tmp.contingencies.requirement(j) = NaN;
%     end
% end
% % SANITY CHECK: contingency requirements
% if ~all(strcmp(tmp.contingencies.requirement,info.task.contingencies.requirement))
%     warning('Mismatch in contingency requirements between SniffinHippo/PyBehaviour files and user input.')
% end

% stim and var
if info.data.numFragments==1 || info.data.fragments_taskContinuous
    temp = extractBetween(tmp.raw.snh.text,strfind(tmp.raw.snh.text,'TRIAL_ORDER:')+13,strfind(tmp.raw.snh.text,'TRIAL_VARIATIONS:')-3);
    for j=1:strlength(temp)
        tmp.raw.snh.stim(j) = str2num([extractBetween(temp,j,j)]);
    end
    temp = extractBetween(tmp.raw.snh.text,strfind(tmp.raw.snh.text,'TRIAL_VARIATIONS:')+18,strfind(tmp.raw.snh.text,';>ARDUINO:')-2);
    for j=1:strlength(temp)
        tmp.raw.snh.var(j) = str2num([extractBetween(temp,j,j)]);
    end
    
    % SANITY CHECK: stim and var
    if ~all(task.type==tmp.raw.snh.stim(1:length(task.type))) | ...
            ~all(task.var==tmp.raw.snh.var(1:length(task.var)))
        warning('Mismatch in stim sequence between SniffinHippo file and TrialSeq file.')
    end
    if ~all(task.type==tmp.raw.snh.stim(1:length(task.type))) | ...
            ~all(task.var==tmp.raw.snh.var(1:length(task.var)))
        warning('Mismatch in stim sequence between SniffinHippo file and TrialSeq file.')
    end
end


%% Dealing with stim trigger issues

if isfield(info.data,'stimTriggerIssue_nanAfter')
    temp = fields(task);
    for i=1:length(temp)
        task.([temp{i},'_raw']) = task.(temp{i});
        task.(temp{i})(info.data.stimTriggerIssue_nanAfter+1:end) = NaN;
    end
    disp('--- Dealt with stim trigger issue. Assigned NaN to task events after issue.')
end


%% Save task file

task = orderfields(task);
save(path.file_out_task,'task','-v7.3');
disp(['--- Added task file to repo as ',path.file_out_task,'.'])


%% Calculate perf struct (general)

% catch, stim, all
temp = [fields(info.task.var);'all'];
for i=1:length(temp)
    this_condition = temp{i};
    
    if strcmp(this_condition,'catch')
        this_var = task.var==0;
        this_string = '_catch';
    elseif strcmp(this_condition,'stim')
        this_var = task.var==1;
        this_string = '_stim';
    elseif strcmp(this_condition,'all')
        this_var = true(1,these_numTrials);
        this_string = '';
    end

    % general
    perf.general.(['correct',this_string]) = (length(find(task.response=="H" & this_var)) + length(find(task.response=="CR" & this_var))) / (length(find(task.response=="H" & this_var)) + length(find(task.response=="CR" & this_var)) + length(find(task.response=="M" & this_var)) + length(find(task.response=="FA" & this_var)));
    perf.general.(['incorrect',this_string]) = 1 - perf.general.(['correct',this_string]);
    perf.general.(['H',this_string]) = length(find(task.response=="H" & this_var)) / (length(find(task.response=="H" & this_var)) + length(find(task.response=="M" & this_var)));
    perf.general.(['M',this_string]) = 1 - perf.general.(['H',this_string]);
    perf.general.(['CR',this_string]) = length(find(task.response=="CR" & this_var)) / (length(find(task.response=="CR" & this_var)) + length(find(task.response=="FA" & this_var)));
    perf.general.(['FA',this_string]) = 1 - perf.general.(['CR',this_string]);
    if ((perf.general.(['H',this_string])>0 && perf.general.(['H',this_string])<1) && (perf.general.(['FA',this_string])>0 && perf.general.(['FA',this_string])<1))
        perf.general.(['dprime',this_string]) = norminv(perf.general.(['H',this_string]))-norminv(perf.general.(['FA',this_string]));
        perf.general.(['crit',this_string]) = -0.5*(norminv(perf.general.(['H',this_string]))+ norminv(perf.general.(['FA',this_string])));
    else
        perf.general.(['dprime',this_string]) = NaN;
        perf.general.(['crit',this_string]) = NaN;
    end
    perf.general = orderfields(perf.general);
    
    % by trial type
    for j=1:length(info.task.contingencies.type)
        if strcmp(info.task.contingencies.requirement(info.task.contingencies.type(j)),"GO")
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['H',this_string]) = length(find((task.response=="H")& this_var&(task.type==info.task.contingencies.type(j)))) / (length(find((task.response=="H")& this_var&(task.type==info.task.contingencies.type(j)))) + length(find((task.response=="M")& this_var&(task.type==info.task.contingencies.type(j)))));
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['M',this_string]) = 1 - perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['H',this_string]);
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['CR',this_string]) = NaN;
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['FA',this_string]) = NaN;
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['correct',this_string]) = perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['H',this_string]);
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['incorrect',this_string]) = 1 - perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['correct',this_string]);
        elseif strcmp(info.task.contingencies.requirement(info.task.contingencies.type(j)),"NOGO")
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['H',this_string]) = NaN;
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['M',this_string]) = NaN;
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['CR',this_string]) = length(find((task.response=="CR")& this_var&(task.type==info.task.contingencies.type(j)))) / (length(find((task.response=="CR")& this_var&(task.type==info.task.contingencies.type(j)))) + length(find((task.response=="FA")& this_var&(task.type==info.task.contingencies.type(j)))));
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['FA',this_string]) = 1 - perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['CR',this_string]);
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['correct',this_string]) = perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['CR',this_string]);
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['incorrect',this_string]) = 1 - perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['correct',this_string]);        
        end
        perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['dprime',this_string]) = NaN;
        perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['crit',this_string]) = NaN;
        perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))) = orderfields(perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))));
    end

    % by first odour
    for j=1:length(info.task.roles.odour1)
        perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['correct',this_string]) = ( length(find((task.response=="H")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="CR")& this_var&(task.odour1==info.task.roles.odour1(j)))) ) / ( length(find((task.response=="H")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="CR")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="M")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="FA")& this_var&(task.odour1==info.task.roles.odour1(j)))) );
        perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['incorrect',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['correct',this_string]);
        perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string]) = length(find((task.response=="H")& this_var&(task.odour1==info.task.roles.odour1(j)))) / (length(find((task.response=="H")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="M")& this_var&(task.odour1==info.task.roles.odour1(j)))));
        perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['M',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string]);
        perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['CR',this_string]) = length(find((task.response=="CR")& this_var&(task.odour1==info.task.roles.odour1(j)))) / (length(find((task.response=="CR")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="FA")& this_var&(task.odour1==info.task.roles.odour1(j)))));
        perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['FA',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['CR',this_string]);
        if ((perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string])>0 && perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string])<1) && (perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['FA',this_string])>0 && perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['FA',this_string])<1))
            perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['dprime',this_string]) = norminv(perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string]))-norminv(perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['FA',this_string]));
            perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['crit',this_string]) = -0.5*(norminv(perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string]))+ norminv(perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['FA',this_string])));
        else
            perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['dprime',this_string]) = NaN;
            perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['crit',this_string]) = NaN;
        end
        perf.(char(strcat('odour',info.task.roles.odour1(j)))) = orderfields(perf.(char(strcat('odour',info.task.roles.odour1(j)))));
    end

    % by second odour
    for j=1:length(info.task.roles.odour2)
        perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['correct',this_string]) = ( length(find((task.response=="H")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="CR")& this_var&(task.odour2==info.task.roles.odour2(j)))) ) / ( length(find((task.response=="H")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="CR")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="M")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="FA")& this_var&(task.odour2==info.task.roles.odour2(j)))) );
        perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['incorrect',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['correct',this_string]);
        perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string]) = length(find((task.response=="H")& this_var&(task.odour2==info.task.roles.odour2(j)))) / (length(find((task.response=="H")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="M")& this_var&(task.odour2==info.task.roles.odour2(j)))));
        perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['M',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string]);
        perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['CR',this_string]) = length(find((task.response=="CR")& this_var&(task.odour2==info.task.roles.odour2(j)))) / (length(find((task.response=="CR")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="FA")& this_var&(task.odour2==info.task.roles.odour2(j)))));
        perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['FA',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['CR',this_string]);
        if ((perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string])>0 && perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string])<1) && (perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['FA',this_string])>0 && perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['FA',this_string])<1))
            perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['dprime',this_string]) = norminv(perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string]))-norminv(perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['FA',this_string]));
            perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['crit',this_string]) = -0.5*(norminv(perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string]))+ norminv(perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['FA',this_string])));
        else
            perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['dprime',this_string]) = NaN;
            perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['crit',this_string]) = NaN;
        end
        perf.(char(strcat('odour',info.task.roles.odour2(j)))) = orderfields(perf.(char(strcat('odour',info.task.roles.odour2(j)))));
    end
end


%% Calculate perf struct (by block)

% catch, stim, all
temp = [fields(info.task.var);'all'];
for i=1:length(temp)
    this_condition = temp{i};
    
    if strcmp(this_condition,'catch')
        this_var = task.var==0;
        this_string = '_catch';
    elseif strcmp(this_condition,'stim')
        this_var = task.var==1;
        this_string = '_stim';
    elseif strcmp(this_condition,'all')
        this_var = true(1,these_numTrials);
        this_string = '';
    end

    for b=1:these_numBlocks

        % general
        perf.blocks_general.(['correct',this_string])(b) = (length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR"))) / (length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")));
        perf.blocks_general.(['incorrect',this_string])(b) = 1 - perf.blocks_general.(['correct',this_string])(b);
        perf.blocks_general.(['H',this_string])(b) = length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")) / (length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")));
        perf.blocks_general.(['M',this_string])(b) = 1 - perf.blocks_general.(['H',this_string])(b);
        perf.blocks_general.(['CR',this_string])(b) = length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")) / (length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")));
        perf.blocks_general.(['FA',this_string])(b) = 1 - perf.blocks_general.(['CR',this_string])(b);
        if ((perf.blocks_general.(['H',this_string])(b)>0 && perf.blocks_general.(['H',this_string])(b)<1) && (perf.blocks_general.(['FA',this_string])(b)>0 && perf.blocks_general.(['FA',this_string])(b)<1))
            perf.blocks_general.(['dprime',this_string])(b) = norminv(perf.blocks_general.(['H',this_string])(b))-norminv(perf.blocks_general.(['FA',this_string])(b));
            perf.blocks_general.(['crit',this_string])(b) = -0.5*(norminv(perf.blocks_general.(['H',this_string])(b))+ norminv(perf.blocks_general.(['FA',this_string])(b)));
        else
            perf.blocks_general.(['dprime',this_string])(b) = NaN;
            perf.blocks_general.(['crit',this_string])(b) = NaN;
        end
        perf.blocks_general = orderfields(perf.blocks_general);

        % by trial type
        for j=1:length(info.task.contingencies.type)
            if strcmp(info.task.contingencies.requirement(info.task.contingencies.type(j)),"GO")
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['H',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))));
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['M',this_string])(b) = 1 - perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['H',this_string])(b);
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['CR',this_string])(b) = NaN;
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['FA',this_string])(b) = NaN;
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['correct',this_string])(b) = perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['H',this_string])(b);
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['incorrect',this_string])(b) = 1 - perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['correct',this_string])(b);
            elseif strcmp(info.task.contingencies.requirement(info.task.contingencies.type(j)),"NOGO")
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['H',this_string])(b) = NaN;
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['M',this_string])(b) = NaN;
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['CR',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))));
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['FA',this_string])(b) = 1 - perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['CR',this_string])(b);
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['correct',this_string])(b) = perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['CR',this_string])(b);
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['incorrect',this_string])(b) = 1 - perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['correct',this_string])(b);        
            end
            perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['dprime',this_string])(b) = NaN;
            perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['crit',this_string])(b) = NaN;
            perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))) = orderfields(perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))));
        end

        % by first odour
        for j=1:length(info.task.roles.odour1)
            perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['correct',this_string])(b) = ( length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) ) / ( length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) );
            perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['incorrect',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['correct',this_string])(b);
            perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))));
            perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['M',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b);
            perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['CR',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))));
            perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['FA',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['CR',this_string])(b);
            if ((perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b)>0 && perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b)<1) && (perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['FA',this_string])(b)>0 && perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['FA',this_string])(b)<1))
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['dprime',this_string])(b) = norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b))-norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['FA',this_string])(b));
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['crit',this_string])(b) = -0.5*(norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b))+ norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['FA',this_string])(b)));
            else
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['dprime',this_string])(b) = NaN;
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['crit',this_string])(b) = NaN;
            end
            perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))) = orderfields(perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))));
        end

        % by second odour
        for j=1:length(info.task.roles.odour2)
            perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['correct',this_string])(b) = ( length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) ) / ( length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) );
            perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['incorrect',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['correct',this_string])(b);
            perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))));
            perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['M',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b);
            perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['CR',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))));
            perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['FA',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['CR',this_string])(b);
            if ((perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b)>0 && perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b)<1) && (perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['FA',this_string])(b)>0 && perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['FA',this_string])(b)<1))
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['dprime',this_string])(b) = norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b))-norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['FA',this_string])(b));
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['crit',this_string])(b) = -0.5*(norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b))+ norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['FA',this_string])(b)));
            else
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['dprime',this_string])(b) = NaN;
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['crit',this_string])(b) = NaN;
            end
            perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))) = orderfields(perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))));
        end
    end
end


%% Calculate perf struct (general) - for different variations (e.g. in expert stim experiments)

if strcmp(epoch,'stim')
    
    temp = unique(task.var);
    for i=1:length(temp)
        this_condition = temp(i);
        this_var = task.var==temp(i);
        this_string = ['_var',num2str(temp(i))];

        % general
        perf.general.(['correct',this_string]) = (length(find(task.response=="H" & this_var)) + length(find(task.response=="CR" & this_var))) / (length(find(task.response=="H" & this_var)) + length(find(task.response=="CR" & this_var)) + length(find(task.response=="M" & this_var)) + length(find(task.response=="FA" & this_var)));
        perf.general.(['incorrect',this_string]) = 1 - perf.general.(['correct',this_string]);
        perf.general.(['H',this_string]) = length(find(task.response=="H" & this_var)) / (length(find(task.response=="H" & this_var)) + length(find(task.response=="M" & this_var)));
        perf.general.(['M',this_string]) = 1 - perf.general.(['H',this_string]);
        perf.general.(['CR',this_string]) = length(find(task.response=="CR" & this_var)) / (length(find(task.response=="CR" & this_var)) + length(find(task.response=="FA" & this_var)));
        perf.general.(['FA',this_string]) = 1 - perf.general.(['CR',this_string]);
        if ((perf.general.(['H',this_string])>0 && perf.general.(['H',this_string])<1) && (perf.general.(['FA',this_string])>0 && perf.general.(['FA',this_string])<1))
            perf.general.(['dprime',this_string]) = norminv(perf.general.(['H',this_string]))-norminv(perf.general.(['FA',this_string]));
            perf.general.(['crit',this_string]) = -0.5*(norminv(perf.general.(['H',this_string]))+ norminv(perf.general.(['FA',this_string])));
        else
            perf.general.(['dprime',this_string]) = NaN;
            perf.general.(['crit',this_string]) = NaN;
        end
        perf.general = orderfields(perf.general);

        % by trial type
        for j=1:length(info.task.contingencies.type)
            if strcmp(info.task.contingencies.requirement(info.task.contingencies.type(j)),"GO")
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['H',this_string]) = length(find((task.response=="H")& this_var&(task.type==info.task.contingencies.type(j)))) / (length(find((task.response=="H")& this_var&(task.type==info.task.contingencies.type(j)))) + length(find((task.response=="M")& this_var&(task.type==info.task.contingencies.type(j)))));
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['M',this_string]) = 1 - perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['H',this_string]);
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['CR',this_string]) = NaN;
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['FA',this_string]) = NaN;
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['correct',this_string]) = perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['H',this_string]);
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['incorrect',this_string]) = 1 - perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['correct',this_string]);
            elseif strcmp(info.task.contingencies.requirement(info.task.contingencies.type(j)),"NOGO")
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['H',this_string]) = NaN;
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['M',this_string]) = NaN;
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['CR',this_string]) = length(find((task.response=="CR")& this_var&(task.type==info.task.contingencies.type(j)))) / (length(find((task.response=="CR")& this_var&(task.type==info.task.contingencies.type(j)))) + length(find((task.response=="FA")& this_var&(task.type==info.task.contingencies.type(j)))));
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['FA',this_string]) = 1 - perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['CR',this_string]);
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['correct',this_string]) = perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['CR',this_string]);
                perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['incorrect',this_string]) = 1 - perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['correct',this_string]);        
            end
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['dprime',this_string]) = NaN;
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))).(['crit',this_string]) = NaN;
            perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))) = orderfields(perf.(char(strcat('type',num2str(info.task.contingencies.type(j))))));
        end

        % by first odour
        for j=1:length(info.task.roles.odour1)
            perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['correct',this_string]) = ( length(find((task.response=="H")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="CR")& this_var&(task.odour1==info.task.roles.odour1(j)))) ) / ( length(find((task.response=="H")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="CR")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="M")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="FA")& this_var&(task.odour1==info.task.roles.odour1(j)))) );
            perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['incorrect',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['correct',this_string]);
            perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string]) = length(find((task.response=="H")& this_var&(task.odour1==info.task.roles.odour1(j)))) / (length(find((task.response=="H")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="M")& this_var&(task.odour1==info.task.roles.odour1(j)))));
            perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['M',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string]);
            perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['CR',this_string]) = length(find((task.response=="CR")& this_var&(task.odour1==info.task.roles.odour1(j)))) / (length(find((task.response=="CR")& this_var&(task.odour1==info.task.roles.odour1(j)))) + length(find((task.response=="FA")& this_var&(task.odour1==info.task.roles.odour1(j)))));
            perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['FA',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['CR',this_string]);
            if ((perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string])>0 && perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string])<1) && (perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['FA',this_string])>0 && perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['FA',this_string])<1))
                perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['dprime',this_string]) = norminv(perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string]))-norminv(perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['FA',this_string]));
                perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['crit',this_string]) = -0.5*(norminv(perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['H',this_string]))+ norminv(perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['FA',this_string])));
            else
                perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['dprime',this_string]) = NaN;
                perf.(char(strcat('odour',info.task.roles.odour1(j)))).(['crit',this_string]) = NaN;
            end
            perf.(char(strcat('odour',info.task.roles.odour1(j)))) = orderfields(perf.(char(strcat('odour',info.task.roles.odour1(j)))));
        end

        % by second odour
        for j=1:length(info.task.roles.odour2)
            perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['correct',this_string]) = ( length(find((task.response=="H")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="CR")& this_var&(task.odour2==info.task.roles.odour2(j)))) ) / ( length(find((task.response=="H")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="CR")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="M")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="FA")& this_var&(task.odour2==info.task.roles.odour2(j)))) );
            perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['incorrect',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['correct',this_string]);
            perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string]) = length(find((task.response=="H")& this_var&(task.odour2==info.task.roles.odour2(j)))) / (length(find((task.response=="H")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="M")& this_var&(task.odour2==info.task.roles.odour2(j)))));
            perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['M',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string]);
            perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['CR',this_string]) = length(find((task.response=="CR")& this_var&(task.odour2==info.task.roles.odour2(j)))) / (length(find((task.response=="CR")& this_var&(task.odour2==info.task.roles.odour2(j)))) + length(find((task.response=="FA")& this_var&(task.odour2==info.task.roles.odour2(j)))));
            perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['FA',this_string]) = 1 - perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['CR',this_string]);
            if ((perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string])>0 && perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string])<1) && (perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['FA',this_string])>0 && perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['FA',this_string])<1))
                perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['dprime',this_string]) = norminv(perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string]))-norminv(perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['FA',this_string]));
                perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['crit',this_string]) = -0.5*(norminv(perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['H',this_string]))+ norminv(perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['FA',this_string])));
            else
                perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['dprime',this_string]) = NaN;
                perf.(char(strcat('odour',info.task.roles.odour2(j)))).(['crit',this_string]) = NaN;
            end
            perf.(char(strcat('odour',info.task.roles.odour2(j)))) = orderfields(perf.(char(strcat('odour',info.task.roles.odour2(j)))));
        end
    end
    
    
    %% Calculate perf struct (by block)

    temp = unique(task.var);
    for i=1:length(temp)
        this_condition = temp(i);
        this_var = task.var==temp(i);
        this_string = ['_var',num2str(temp(i))];

        for b=1:these_numBlocks

            % general
            perf.blocks_general.(['correct',this_string])(b) = (length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR"))) / (length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")));
            perf.blocks_general.(['incorrect',this_string])(b) = 1 - perf.blocks_general.(['correct',this_string])(b);
            perf.blocks_general.(['H',this_string])(b) = length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")) / (length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")));
            perf.blocks_general.(['M',this_string])(b) = 1 - perf.blocks_general.(['H',this_string])(b);
            perf.blocks_general.(['CR',this_string])(b) = length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")) / (length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")) + length(find(this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")));
            perf.blocks_general.(['FA',this_string])(b) = 1 - perf.blocks_general.(['CR',this_string])(b);
            if ((perf.blocks_general.(['H',this_string])(b)>0 && perf.blocks_general.(['H',this_string])(b)<1) && (perf.blocks_general.(['FA',this_string])(b)>0 && perf.blocks_general.(['FA',this_string])(b)<1))
                perf.blocks_general.(['dprime',this_string])(b) = norminv(perf.blocks_general.(['H',this_string])(b))-norminv(perf.blocks_general.(['FA',this_string])(b));
                perf.blocks_general.(['crit',this_string])(b) = -0.5*(norminv(perf.blocks_general.(['H',this_string])(b))+ norminv(perf.blocks_general.(['FA',this_string])(b)));
            else
                perf.blocks_general.(['dprime',this_string])(b) = NaN;
                perf.blocks_general.(['crit',this_string])(b) = NaN;
            end
            perf.blocks_general = orderfields(perf.blocks_general);

            % by trial type
            for j=1:length(info.task.contingencies.type)
                if strcmp(info.task.contingencies.requirement(info.task.contingencies.type(j)),"GO")
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['H',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))));
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['M',this_string])(b) = 1 - perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['H',this_string])(b);
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['CR',this_string])(b) = NaN;
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['FA',this_string])(b) = NaN;
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['correct',this_string])(b) = perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['H',this_string])(b);
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['incorrect',this_string])(b) = 1 - perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['correct',this_string])(b);
                elseif strcmp(info.task.contingencies.requirement(info.task.contingencies.type(j)),"NOGO")
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['H',this_string])(b) = NaN;
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['M',this_string])(b) = NaN;
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['CR',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")&(task.type((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.contingencies.type(j)))));
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['FA',this_string])(b) = 1 - perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['CR',this_string])(b);
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['correct',this_string])(b) = perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['CR',this_string])(b);
                    perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['incorrect',this_string])(b) = 1 - perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['correct',this_string])(b);        
                end
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['dprime',this_string])(b) = NaN;
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))).(['crit',this_string])(b) = NaN;
                perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))) = orderfields(perf.(char(strcat('blocks_type',num2str(info.task.contingencies.type(j))))));
            end

            % by first odour
            for j=1:length(info.task.roles.odour1)
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['correct',this_string])(b) = ( length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) ) / ( length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) );
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['incorrect',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['correct',this_string])(b);
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))));
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['M',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b);
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['CR',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")&(task.odour1((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour1(j)))));
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['FA',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['CR',this_string])(b);
                if ((perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b)>0 && perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b)<1) && (perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['FA',this_string])(b)>0 && perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['FA',this_string])(b)<1))
                    perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['dprime',this_string])(b) = norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b))-norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['FA',this_string])(b));
                    perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['crit',this_string])(b) = -0.5*(norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['H',this_string])(b))+ norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['FA',this_string])(b)));
                else
                    perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['dprime',this_string])(b) = NaN;
                    perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))).(['crit',this_string])(b) = NaN;
                end
                perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))) = orderfields(perf.(char(strcat('blocks_odour',info.task.roles.odour1(j)))));
            end

            % by second odour
            for j=1:length(info.task.roles.odour2)
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['correct',this_string])(b) = ( length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) ) / ( length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) );
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['incorrect',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['correct',this_string])(b);
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="H")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="M")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))));
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['M',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b);
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['CR',this_string])(b) = length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) / (length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="CR")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))) + length(find((this_var((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock) & task.response((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)=="FA")&(task.odour2((b-1)*these_trialsPerBlock+1:b*these_trialsPerBlock)==info.task.roles.odour2(j)))));
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['FA',this_string])(b) = 1 - perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['CR',this_string])(b);
                if ((perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b)>0 && perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b)<1) && (perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['FA',this_string])(b)>0 && perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['FA',this_string])(b)<1))
                    perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['dprime',this_string])(b) = norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b))-norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['FA',this_string])(b));
                    perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['crit',this_string])(b) = -0.5*(norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['H',this_string])(b))+ norminv(perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['FA',this_string])(b)));
                else
                    perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['dprime',this_string])(b) = NaN;
                    perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))).(['crit',this_string])(b) = NaN;
                end
                perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))) = orderfields(perf.(char(strcat('blocks_odour',info.task.roles.odour2(j)))));
            end
        end
    end
end

%% Save perf file

perf = orderfields(perf);
save(path.file_out_perf,'perf','-v7.3');
disp(['--- Added perf file to repo as ',path.file_out_perf,'.'])


end
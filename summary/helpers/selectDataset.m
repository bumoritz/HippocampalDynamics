function d_info = selectDataset(d_info,dataset,sheet,path,ops)
% dataset = '-g78-d12345-e1-r01-p1-l01-i00';
% dataset ='-g0123456789101112-d12345-e01-r01-p01-l01-i00';
% dataset ='-g0123456789101112-d12345-e1-r01-p01-l01-i00'

% format:
% - group (g1=beh, g2=img, g3-6=old_stim, g7=stim1, g8=stim2, g9=inh1, g10=inh2, g11=inhc1, g12=inhc2, g0123456789101112=all)
% - days (d1=day1, d12345=all_days, d2345=stim_days)
% - engagement (e0=disengaged, e1=engaged, e01=both)
% - running (r0=non_runners, r1=runners, r01=both)
% - photostim-responsive [2P stim only] (p0=bad_responses, p1=good_responses, p01=both)
% - learning (l0=non_learners, l1=learners, l01=both)
% - incomplete animal (i00=no_filter, i11=complete, i22=cmplt_within_switches, i33=cmplt_across_switches)


%% Preparations

disp('- Selecting dataset.')

d_info.excl = true(d_info.numAnimals,d_info.numDays);

% get instructions which dataset to load
instructions = split(dataset,'-');

% -g
instructions_g = instructions(find(contains(instructions,'g')));
instructions_g = instructions_g{1}(2:end);
if strcmp(instructions_g,'1')
	this_g = d_info.group==1;
end
if strcmp(instructions_g,'2')
	this_g = d_info.group==2;
end
if strcmp(instructions_g,'7')
	this_g = d_info.group==7;
end
if strcmp(instructions_g,'8')
	this_g = d_info.group==8;
end
if strcmp(instructions_g,'12')
	this_g = d_info.group==1 | d_info.group==2;
end
if strcmp(instructions_g,'78')
	this_g = d_info.group==7 | d_info.group==8;
end
if strcmp(instructions_g,'1278')
	this_g = d_info.group==1 | d_info.group==2 | d_info.group==7 | d_info.group==8;
end
if strcmp(instructions_g,'278')
	this_g = d_info.group==2 | d_info.group==7 | d_info.group==8;
end
if strcmp(instructions_g,'02345678')
    this_g = d_info.group==0 | d_info.group==2 | d_info.group==3 | d_info.group==4 | d_info.group==5 | d_info.group==6 | d_info.group==7 | d_info.group==8;
end
if strcmp(instructions_g,'0123456789101112')
	this_g = d_info.group==0 | d_info.group==1 | d_info.group==2 | d_info.group==3 | d_info.group==4 | d_info.group==5 | d_info.group==6 | d_info.group==7 | d_info.group==8 | d_info.group==9 | d_info.group==10 | d_info.group==11 | d_info.group==12;
end
if strcmp(instructions_g,'012345678')
	this_g = d_info.group==0 | d_info.group==1 | d_info.group==2 | d_info.group==3 | d_info.group==4 | d_info.group==5 | d_info.group==6 | d_info.group==7 | d_info.group==8;
end

% -d
instructions_d = instructions(find(contains(instructions,'d')));
instructions_d = instructions_d{1}(2:end);
if strcmp(instructions_d,'1')
	this_d = 1;
end
if strcmp(instructions_d,'12345')
	this_d = 1:5;
end
if strcmp(instructions_d,'2345')
	this_d = 2:5;
end

% apply -g and -d
d_info.excl(this_g,this_d) = false;


%% Exclude sessions that should never be included (because they are corrupted for whatever reason)

% exclude sessions with old stim design
d_info.excl(7,2:5) = true;      % Carlo (old 2P stim II)
d_info.excl(9,2:5) = true;      % Cardano (old 2P stim I)
d_info.excl(10,2:5) = true;     % Frida (old 2P stim III, super-biased)
d_info.excl(13,2:5) = true;     % Ripple (old 2P stim IV)

% exclude sessions with fucked-up targeting (excluding only the seq session)
d_info.excl(15,3) = true;       % Faramir (2P stim I, second switch days)
d_info.excl(18,5) = true;       % Merry (2P stim II, second switch days)
d_info.excl(19,3) = true;       % Celo (2P stim I, second switch days)
d_info.excl(22,3) = true;       % Sterling (2P stim I, second switch day fucked up, non-engaged)

% exclude sessions with fucked-up software-alignment
d_info.excl(34,2:5) = true;     % Meghan (2P stim II, fucked-up software-alignment)
d_info.excl(36,2:5) = true;     % Python (2P stim II, fucked-up software-alignment)

% exclude sessions with stim trigger issue
d_info.excl(21,2:3) = true;     % Shaw (2P stim II, missed a stim trigger issue early on during first switch session)
d_info.excl(38,5) = true;       % BullyBoy (2P stim II,, missed a stim trigger issue early on)

% exclude sessions with corrupted files
d_info.excl(22,2) = true;       % Sterling (2P stim I)

% remove sessions that have never been recorded
d_info.excl(1,2:5) = true;
d_info.excl(2,2:5) = true;
d_info.excl(14,4:5) = true;
d_info.excl(22,4:5) = true;
d_info.excl(26,3:5) = true;
d_info.excl(29,5) = true;
d_info.excl(37,3:5) = true;
d_info.excl(53,4:5) = true;


%% Get engagement

instructions_e = instructions(find(contains(instructions,'e')));
instructions_e = instructions_e{1}(2:end);
if ops.excl.loadDefault_e
    load([path.root_summary,'defaults\engagement.mat']);
    d_info.engagement = engagement;
else
    disp('--- Calculating engagement.')
    cmpr_list_repo = ["task"];
    d_info.engagement = nan(size(d_info.excl));
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            this_animal = d_info.animals{i};
            try
                try
                    temp = (~isempty(sheet.(['Day',num2str(j)]){i})) && (~d_info.excl(i,j));
                    this_date = datestr(sheet.(['Day',num2str(j)]){i},'yyyymmdd');
                catch
                    temp = (~isempty(sheet.(['Day',num2str(j)])(i))) && (~d_info.excl(i,j));
                    this_date = datestr(sheet.(['Day',num2str(j)])(i),'yyyymmdd');
                end
                if temp
                    % load repo data
                    path.folder_repo = [path.root_repo,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in = [path.folder_repo,this_animal,'_',this_date,'_'];
                    path.folder_analysis = [path.root_analysis,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in_analysis = [path.folder_analysis,this_animal,'_',this_date,'_'];
                    path.folder_analysisX = [path.root_analysisX,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in_analysisX = [path.folder_analysisX,this_animal,'_',this_date,'_'];
                    this_d = repo2ws(path,cmpr_list_repo,ops.skipIncompletelyProcessed);
                    temp = (length(find(this_d.task.response=='H'))+length(find(this_d.task.response=='FA'))) / length(this_d.task.response) >= p.excl.engagement;
                    if ~isnan(temp)
                        d_info.engagement(i,j) = temp;
                    end
                end
            catch
                warning(['Had to skip ',this_animal,' when calculating engagement.'])
            end
        end
    end
    if ops.excl.updateDefault_e
        engagement = d_info.engagement;
        save([path.root_summary,'defaults\engagement.mat'],'engagement');
    end
end


%% Get running

instructions_r = instructions(find(contains(instructions,'r')));
instructions_r = instructions_r{1}(2:end);
if ops.excl.loadDefault_r
    load([path.root_summary,'defaults\running.mat']);
    d_info.running = running;
else
    disp('--- Calculating running.')
    cmpr_list_repo = ["paq_beh"];
    d_info.running = nan(size(d_info.excl));
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            this_animal = d_info.animals{i};
            try
                try
                    temp = (~isempty(sheet.(['Day',num2str(j)]){i})) && (~d_info.excl(i,j));
                    this_date = datestr(sheet.(['Day',num2str(j)]){i},'yyyymmdd');
                catch
                    temp = (~isempty(sheet.(['Day',num2str(j)])(i))) && (~d_info.excl(i,j));
                    this_date = datestr(sheet.(['Day',num2str(j)])(i),'yyyymmdd');
                end
                if temp

                    % load repo data
                    path.folder_repo = [path.root_repo,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in = [path.folder_repo,this_animal,'_',this_date,'_'];
                    path.folder_analysis = [path.root_analysis,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in_analysis = [path.folder_analysis,this_animal,'_',this_date,'_'];
                    path.folder_analysisX = [path.root_analysisX,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in_analysisX = [path.folder_analysisX,this_animal,'_',this_date,'_'];
                    this_d = repo2ws(path,cmpr_list_repo,ops.skipIncompletelyProcessed);
                    temp = nanmean(this_d.paq_beh.speed);
                    if ~isnan(temp)
                        d_info.running(i,j) = temp >= p.excl.running;
                    end
                end
            catch
                warning(['Had to skip ',this_animal,'_',this_date,' when calculating running.'])
            end
        end
    end
    % add info for animals without rotary encoder or where it's NaN for other reasons
    d_info.running(3,4) = 0; % Biontech
    d_info.running(15,3) = 0; % Faramir
    d_info.running(16,1:5) = 1; % Arwen
    d_info.running(18,5) = 1; % Merry
    d_info.running(19,[3,5]) = 1; % Celo
	d_info.running(27,1:5) = 0; % Stanage
    d_info.running(47,1) = 1; % Austin
    d_info.running(50,4) = 0; % Turing
    d_info.running(53,3) = 1; % Berners
    if ops.excl.updateDefault_r
        running = d_info.running;
        save([path.root_summary,'defaults\running.mat'],'running');
    end
end


%% Get speed

if ops.excl.loadDefault_s
    load([path.root_summary,'defaults\speed_sess.mat']);
    load([path.root_summary,'defaults\speed_AW.mat']);
    d_info.speed_sess = speed_sess;
    d_info.speed_AW = speed_AW;
else
    disp('--- Calculating speed.')
    cmpr_list_repo = ["bcon"];
    d_info.speed_sess = nan(size(d_info.excl));
    d_info.speed_AW = nan(size(d_info.excl));
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            this_animal = d_info.animals{i};
            try
                try
                    temp = (~isempty(sheet.(['Day',num2str(j)]){i})) && (~d_info.excl(i,j));
                    this_date = datestr(sheet.(['Day',num2str(j)]){i},'yyyymmdd');
                catch
                    temp = (~isempty(sheet.(['Day',num2str(j)])(i))) && (~d_info.excl(i,j));
                    this_date = datestr(sheet.(['Day',num2str(j)])(i),'yyyymmdd');
                end
                if temp

                    % load repo data
                    path.folder_repo = [path.root_repo,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in = [path.folder_repo,this_animal,'_',this_date,'_'];
                    path.folder_analysis = [path.root_analysis,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in_analysis = [path.folder_analysis,this_animal,'_',this_date,'_'];
                    path.folder_analysisX = [path.root_analysisX,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in_analysisX = [path.folder_analysisX,this_animal,'_',this_date,'_'];
                    this_d = repo2ws(path,cmpr_list_repo,ops.skipIncompletelyProcessed);
                    d_info.speed_sess(i,j) = nanmean(this_d.bcon.trialwise_full.velocity);
                    d_info.speed_AW(i,j) = nanmean(this_d.bcon.trialwise_AW.velocity);
                end
            catch
                warning(['Had to skip ',this_animal,'_',this_date,' when calculating speed.'])
            end
        end
    end
    if ops.excl.updateDefault_s
        speed_sess = d_info.speed_sess;
        save([path.root_summary,'defaults\speed_sess.mat'],'speed_sess');
        speed_AW = d_info.speed_AW;
        save([path.root_summary,'defaults\speed_AW.mat'],'speed_AW');
    end
end


%% Get learning (above chance)

instructions_l = instructions(find(contains(instructions,'l')));
instructions_l = instructions_l{1}(2:end);
if ops.excl.loadDefault_l
    load([path.root_summary,'defaults\learning.mat']);
    d_info.learning = learning;
else
    disp('--- Calculating learning.')
    cmpr_list_repo = ["perf"];
    d_info.learning = nan(size(d_info.excl));
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            this_animal = d_info.animals{i};
            try
                try
                    temp = (~isempty(sheet.(['Day',num2str(j)]){i})) && (~d_info.excl(i,j));
                    this_date = datestr(sheet.(['Day',num2str(j)]){i},'yyyymmdd');
                catch
                    temp = (~isempty(sheet.(['Day',num2str(j)])(i))) && (~d_info.excl(i,j));
                    this_date = datestr(sheet.(['Day',num2str(j)])(i),'yyyymmdd');
                end
                if temp

                    % load repo data
                    path.folder_repo = [path.root_repo,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in = [path.folder_repo,this_animal,'_',this_date,'_'];
                    path.folder_analysis = [path.root_analysis,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in_analysis = [path.folder_analysis,this_animal,'_',this_date,'_'];
                    path.folder_analysisX = [path.root_analysisX,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in_analysisX = [path.folder_analysisX,this_animal,'_',this_date,'_'];
                    this_d = repo2ws(path,cmpr_list_repo,ops.skipIncompletelyProcessed);
                    temp = signrank(this_d.perf.blocks_general.correct-0.5);
                    if ~isnan(temp)
                        d_info.learning(i,j) = (temp < 0.05) && nanmean(this_d.perf.blocks_general.correct) > 0.5;
                    end
                end
            catch
                warning(['Had to skip ',this_animal,'_',this_date,' when calculating learning.'])
            end
        end
    end
    if ops.excl.updateDefault_l
        learning = d_info.learning;
        save([path.root_summary,'defaults\learning.mat'],'learning');
    end
end


%% Get presponsive

instructions_p = instructions(find(contains(instructions,'p')));
instructions_p = instructions_p{1}(2:end);
if ops.excl.loadDefault_p
    load([path.root_summary,'defaults\presponsive.mat']);
    d_info.presponsive = presponsive;
else
    disp('--- Calculating presponsive.')
    cmpr_list_repo = ["resp"];
    d_info.presponsive = nan(size(d_info.excl));
    for i=1:d_info.numAnimals
        for j=2:d_info.numDays
            this_animal = d_info.animals{i};
            try
                try
                    temp = (~isempty(sheet.(['Day',num2str(j)]){i})) && (~d_info.excl(i,j));
                    this_date = datestr(sheet.(['Day',num2str(j)]){i},'yyyymmdd');
                catch
                    temp = (~isempty(sheet.(['Day',num2str(j)])(i))) && (~d_info.excl(i,j));
                    this_date = datestr(sheet.(['Day',num2str(j)])(i),'yyyymmdd');
                end
                if temp

                    % load repo data
                    path.folder_repo = [path.root_repo,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in = [path.folder_repo,this_animal,'_',this_date,'_'];
                    path.folder_analysis = [path.root_analysis,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in_analysis = [path.folder_analysis,this_animal,'_',this_date,'_'];
                    path.folder_analysisX = [path.root_analysisX,this_animal,'/',this_animal,'_',this_date,'/'];
                    path.filepart_in_analysisX = [path.folder_analysisX,this_animal,'_',this_date,'_'];
                    this_d = repo2ws(path,cmpr_list_repo,ops.skipIncompletelyProcessed);
                    temp = this_d.resp.responders_main.proportion_RespTargeted_AllTargeted;
                    if ~isnan(temp)
                        d_info.presponsive(i,j) = temp > p.excl.presponsive;
                    end
                end
            catch
                warning(['Had to skip ',this_animal,'_',this_date,' when calculating presponsive.'])
            end
        end
    end
    if ops.excl.updateDefault_p
        presponsive = d_info.presponsive;
        save([path.root_summary,'defaults\presponsive.mat'],'presponsive');
    end
end


%% Apply engagement, running, presponsive, learning and incomplete filters

% apply -e
if ~strcmp(instructions_e,'01')
    if strcmp(instructions_e,'0')
        this_target = 0;
    elseif strcmp(instructions_e,'1')
        this_target = 1;
    end
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if isnan(d_info.engagement(i,j)) || d_info.excl(i,j)==true
                d_info.excl(i,j) = true;
            elseif d_info.engagement(i,j) == this_target
                d_info.excl(i,j) = false;
            elseif d_info.engagement(i,j) ~= this_target
                d_info.excl(i,j) = true;
            end
        end
    end
end

% apply -r
if ~strcmp(instructions_r,'01')
    if strcmp(instructions_r,'0')
        this_target = 0;
    elseif strcmp(instructions_r,'1')
        this_target = 1;
    end
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if isnan(d_info.running(i,j)) || d_info.excl(i,j)==true
                d_info.excl(i,j) = true;
            elseif d_info.running(i,j) == this_target
                d_info.excl(i,j) = false;
            elseif d_info.running(i,j) ~= this_target
                d_info.excl(i,j) = true;
            end
        end
    end
end

% apply -l
if ~strcmp(instructions_l,'01')
    if strcmp(instructions_l,'0')
        this_target = 0;
    elseif strcmp(instructions_l,'1')
        this_target = 1;
    end
    for i=1:d_info.numAnimals
        for j=1:d_info.numDays
            if isnan(d_info.learning(i,j)) || d_info.excl(i,j)==true
                d_info.excl(i,j) = true;
            elseif d_info.learning(i,j) == this_target
                d_info.excl(i,j) = false;
            elseif d_info.learning(i,j) ~= this_target
                d_info.excl(i,j) = true;
            end
        end
    end
end

% apply -p
if ~strcmp(instructions_p,'01')
    if strcmp(instructions_p,'0')
        this_target = 0;
    elseif strcmp(instructions_p,'1')
        this_target = 1;
    end
    for i=1:d_info.numAnimals
        if  (d_info.group(i)==7 | d_info.group(i)==8)
            for j=2:d_info.numDays
                if isnan(d_info.presponsive(i,j)) || d_info.excl(i,j)==true
                    d_info.excl(i,j) = true;
                elseif d_info.presponsive(i,j) == this_target
                    d_info.excl(i,j) = false;
                elseif d_info.presponsive(i,j) ~= this_target
                    d_info.excl(i,j) = true;
                end
            end
        end
    end
end

% apply -i
instructions_i = instructions(find(contains(instructions,'i')));
instructions_i = instructions_i{1}(2:end);
if ~strcmp(instructions_i,'00')
    if strcmp(instructions_p,'11')
        this_i = find(sum(d_info.excl(:,2:5),2)==0);
        d_info.excl(setdiff(1:d_info.numAnimals,this_i),1:5)=1;
    elseif strcmp(instructions_p,'22')
       % this_i = find(sum(d_info.excl(:,2:3),2)==0)  find(sum(d_info.excl(:,4:5),2)==0);
    elseif strcmp(instructions_p,'33')
        % this_i = find(sum(d_info.excl(:,[2,4]),2)==0);
    end
end
% - incomplete animal (i00=no_filter, i11=complete, i22=cmplt_within_switches, i33=cmplt_across_switches)

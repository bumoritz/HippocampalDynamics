function trl = extractTrl(loadPath_SEQ,loadPath_SNH,loadPath_PYB,savePath,p)

% loadPath_SEQ = path.trlDataPath_SEQ
% loadPath_SNH = path.trlDataPath_SNH
% loadPath_PYB = path.trlDataPath_PYB
% savePath = path.trlAnalysisPath
% p = p.trl

% MB20200808: modified to work with 6-port olfactometer

%% Load data

% seq file
if isfile(loadPath_SEQ)
    temp = readmatrix(loadPath_SEQ);
    trl.raw.seq.stim = temp(1,:);
    trl.raw.seq.var = temp(2,:);
else
    error(['TrialSeq file at ',loadPath_SEQ,' does not exist.'])
end

% snh file
if isfile(loadPath_SNH)
    fileID = fopen(loadPath_SNH);
    C = textscan(fileID,'%c');
    trl.raw.snh.text = convertCharsToStrings(C{1,1});
    fclose(fileID);
else
    error(['SniffinHippo file at ',loadPath_SNH,' does not exist.'])
end

% pyb file
if isfile(loadPath_PYB)
    trl.raw.pyb = load(loadPath_PYB);
else
    temp = dir([fileparts(loadPath_PYB),'\*_PYB.txt']);
    loadPath_PYB_alt = [fileparts(loadPath_PYB),'\',temp.name];
    if isfile(loadPath_PYB_alt)
        warning(['PyBehaviour mat file at ',loadPath_PYB,' does not exist. Using the txt file instead.'])
        trl.raw.pyb.results = generateMatFromTxtFile(loadPath_PYB_alt);
    else
        error(['Both PyBehaviour files at ',loadPath_PYB,' do not exist.'])
    end
end

%% General properties

trl.p = p;
trl.extractionDate = datetime(now,'ConvertFrom','datenum');

%% Extract information from SniffinHippo

% stimStructure
trl.stimStructure.tOdour1 = str2num(extractBetween(trl.raw.snh.text,strfind(trl.raw.snh.text,'T_ODOUR1:')+9,strfind(trl.raw.snh.text,'T_GAP:')-2))/1000; % [s]
trl.stimStructure.tGap = str2num(extractBetween(trl.raw.snh.text,strfind(trl.raw.snh.text,'T_GAP:')+6,strfind(trl.raw.snh.text,'T_ODOUR2:')-2))/1000; % [s]
trl.stimStructure.tOdour2 = str2num(extractBetween(trl.raw.snh.text,strfind(trl.raw.snh.text,'T_ODOUR2:')+9,strfind(trl.raw.snh.text,'TYPE_CONTS:')-2))/1000; % [s]

% vials
temp = extractBetween(trl.raw.snh.text,strfind(trl.raw.snh.text,'VIAL_CONTS:')+11,strfind(trl.raw.snh.text,'NUM_TRIALS:')-2);
temp2 = strfind(temp,',');
trl.vials.vialNumber = [];
trl.vials.role = strings;
trl.vials.odourant = strings;
for i=0:length(temp2)
    if i==0
        trl.vials.vialNumber(i+1) = str2num(extractBetween(temp,1,1));
        trl.vials.role(i+1) = extractBetween(temp,2,2);
        trl.vials.odourant(i+1) = extractBetween(temp,3,temp2(1)-1)
    else
        trl.vials.vialNumber(i+1) = str2num(extractBetween(temp,temp2(i)+1,temp2(i)+1));
        trl.vials.role(i+1) = extractBetween(temp,temp2(i)+2,temp2(i)+2);
        if i==length(temp2)
            trl.vials.odourant(i+1) = extractBetween(temp,temp2(i)+3,strlength(temp));
        else
            trl.vials.odourant(i+1) = extractBetween(temp,temp2(i)+3,temp2(i+1)-1);
        end
    end
end
    
% contingencies
temp = extractBetween(trl.raw.snh.text,strfind(trl.raw.snh.text,'TYPE_CONTS:')+11,strfind(trl.raw.snh.text,'VIAL_CONTS:')-2);
trl.contingencies.stim = [];
trl.contingencies.var = [];
trl.contingencies.odour1 = strings;
trl.contingencies.odour2 = strings;
for i=1:strlength(temp)
    if mod(i,5)==1
        trl.contingencies.stim = [trl.contingencies.stim; str2num(extractBetween(temp,i,i))];
    elseif mod(i,5)==2
        trl.contingencies.var = [trl.contingencies.var; str2num(extractBetween(temp,i,i))];
    elseif mod(i,5)==3
        trl.contingencies.odour1 = [trl.contingencies.odour1; extractBetween(temp,i,i)];
    elseif mod(i,5)==4
        trl.contingencies.odour2 = [trl.contingencies.odour2; extractBetween(temp,i,i)];
    end
end
trl.contingencies.odour1 = trl.contingencies.odour1(2:end);
trl.contingencies.odour2 = trl.contingencies.odour2(2:end);

% roles
trl.roles.odour1 = unique(trl.contingencies.odour1);
trl.roles.odour2 = unique(trl.contingencies.odour2);

% stim and var
temp = extractBetween(trl.raw.snh.text,strfind(trl.raw.snh.text,'TRIAL_ORDER:')+13,strfind(trl.raw.snh.text,'TRIAL_VARIATIONS:')-3);
for i=1:strlength(temp)
    trl.raw.snh.stim(i) = str2num([extractBetween(temp,i,i)]);
end
temp = extractBetween(trl.raw.snh.text,strfind(trl.raw.snh.text,'TRIAL_VARIATIONS:')+18,strfind(trl.raw.snh.text,';>ARDUINO:')-2);
for i=1:strlength(temp)
    trl.raw.snh.var(i) = str2num([extractBetween(temp,i,i)]);
end

%% Extract information from PyBehaviour

if isfield(trl.raw.pyb.results{1,size(trl.raw.pyb.results,2)},'stim_type')
    trl.numTrials = size(trl.raw.pyb.results,2);
else
    trl.numTrials = size(trl.raw.pyb.results,2)-1;
end

for t=1:trl.numTrials
    
    % stim type
    trl.stim(t) = trl.raw.seq.stim(t); %MB20200519, trl.raw.pyb.results{1,t}.stim_type;
    trl.var(t) = trl.raw.pyb.results{1,t}.stim_var;

    % requirement
    if trl.raw.pyb.results{1,t}.response_required % | trl.stim(t)==5 | trl.stim(t)==6)
        trl.requirement(t) = "GO";
    else
        trl.requirement(t) = "NOGO";
    end
    
    % response
    if trl.raw.pyb.results{1,t}.cheated
        trl.response(t) = "CHEAT";
    elseif trl.raw.pyb.results{1,t}.incorrect
        trl.response(t) = "FA";
    elseif trl.raw.pyb.results{1,t}.miss
        trl.response(t) = "M";
    elseif (trl.raw.pyb.results{1,t}.correct & trl.raw.pyb.results{1,t}.response_required)
        trl.response(t) = "H";
    elseif (trl.raw.pyb.results{1,t}.correct & ~trl.raw.pyb.results{1,t}.response_required)
        trl.response(t) = "CR";
    else
        trl.response(t) = "NaN";
    end
    
%     % Rusviet 20201127
%     if trl.raw.pyb.results{1,t}.cheated
%         trl.response(t) = "CHEAT";
%     elseif (trl.raw.pyb.results{1,t}.incorrect & (trl.stim(t)~=5 & trl.stim(t)~=6))
%         trl.response(t) = "FA";
%     elseif (trl.raw.pyb.results{1,t}.miss | (trl.raw.pyb.results{1,t}.correct & (trl.stim(t)==5 | trl.stim(t)==6)))
%         trl.response(t) = "M";
%     elseif (trl.raw.pyb.results{1,t}.correct & trl.raw.pyb.results{1,t}.response_required) | (trl.raw.pyb.results{1,t}.incorrect & (trl.stim(t)==5 | trl.stim(t)==6))
%         trl.response(t) = "H";
%     elseif ((trl.raw.pyb.results{1,t}.correct & ~trl.raw.pyb.results{1,t}.response_required) & (trl.stim(t)~=5 & trl.stim(t)~=6))
%         trl.response(t) = "CR";
%     else
%         trl.response(t) = "NaN";
%     end
end

%% Combine information from SniffinHippo and PyBehaviour

% odours
for t=1:trl.numTrials
    trl.odour1(t) = trl.contingencies.odour1(find(unique(trl.contingencies.stim)==trl.stim(t))); %MB20201027: unique
    trl.odour2(t) = trl.contingencies.odour2(find(unique(trl.contingencies.stim)==trl.stim(t))); %MB20201027: unique  
end

% put requirements into contingencies field
for s=1:length(unique(trl.contingencies.stim)) %MB20201027: unique
    if ~isempty(min(find(trl.stim==trl.contingencies.stim(s))))
        trl.contingencies.requirement(s) = trl.requirement(min(find(trl.stim==trl.contingencies.stim(s))));
    else
        trl.contingencies.requirement(s) = NaN;
    end
end
trl.contingencies.requirement = trl.contingencies.requirement';


%% Sanity checks

% Are all stim and var sequences matching?

%% Calculate performance metrics

%trl.contingencies.stim = [1:8]; %MB20201202

% general

% for i=1:3
%     if i==1
%         this_varLabel = '_catch';
%         this_varVect = trl.var==0;
%     elseif i==2
%         this_varLabel = '_stim';
%         this_varVect = trl.var==1;
%     else
%         this_varLabel = '';
%         this_varVect = ones(1,trl.numTrials);
%     end

    trl.performance.general.correct = (length(find(trl.response=="H")) + length(find(trl.response=="CR"))) / (length(find(trl.response=="H")) + length(find(trl.response=="CR")) + length(find(trl.response=="M")) + length(find(trl.response=="FA")));
    trl.performance.general.incorrect = 1 - trl.performance.general.correct;
    trl.performance.general.H = length(find(trl.response=="H")) / (length(find(trl.response=="H")) + length(find(trl.response=="M")));
    trl.performance.general.M = 1 - trl.performance.general.H;
    trl.performance.general.CR = length(find(trl.response=="CR")) / (length(find(trl.response=="CR")) + length(find(trl.response=="FA")));
    trl.performance.general.FA = 1 - trl.performance.general.CR;
    if ((trl.performance.general.H>0 && trl.performance.general.H<1) && (trl.performance.general.FA>0 && trl.performance.general.FA<1))
        [trl.performance.general.dprime,trl.performance.general.bias] = dprime_simple(trl.performance.general.H,trl.performance.general.FA);
    else
        trl.performance.general.dprime = NaN;
        trl.performance.general.bias = NaN;
    end
        
    % by trial type
    for s=1:length(trl.contingencies.stim)
        if strcmp(trl.contingencies.requirement(trl.contingencies.stim(s)),"GO")
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).H = length(find((trl.response=="H")&(trl.stim==trl.contingencies.stim(s)))) / (length(find((trl.response=="H")&(trl.stim==trl.contingencies.stim(s)))) + length(find((trl.response=="M")&(trl.stim==trl.contingencies.stim(s)))));
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).M = 1 - trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).H;
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).CR = NaN;
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).FA = NaN;
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).correct = trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).H;
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).incorrect = 1 - trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).correct;
        elseif strcmp(trl.contingencies.requirement(trl.contingencies.stim(s)),"NOGO")
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).H = NaN;
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).M = NaN;
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).CR = length(find((trl.response=="CR")&(trl.stim==trl.contingencies.stim(s)))) / (length(find((trl.response=="CR")&(trl.stim==trl.contingencies.stim(s)))) + length(find((trl.response=="FA")&(trl.stim==trl.contingencies.stim(s)))));
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).FA = 1 - trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).CR;
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).correct = trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).CR;
            trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).incorrect = 1 - trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).correct;        
        end
        trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).dprime = NaN;
        trl.performance.(char(strcat('type',num2str(trl.contingencies.stim(s))))).bias = NaN;
    end

    % by first odour
    for r=1:length(trl.roles.odour1)
        trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).correct = ( length(find((trl.response=="H")&(trl.odour1==trl.roles.odour1(r)))) + length(find((trl.response=="CR")&(trl.odour1==trl.roles.odour1(r)))) ) / ( length(find((trl.response=="H")&(trl.odour1==trl.roles.odour1(r)))) + length(find((trl.response=="CR")&(trl.odour1==trl.roles.odour1(r)))) + length(find((trl.response=="M")&(trl.odour1==trl.roles.odour1(r)))) + length(find((trl.response=="FA")&(trl.odour1==trl.roles.odour1(r)))) );
        trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).incorrect = 1 - trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).correct;
        trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).H = length(find((trl.response=="H")&(trl.odour1==trl.roles.odour1(r)))) / (length(find((trl.response=="H")&(trl.odour1==trl.roles.odour1(r)))) + length(find((trl.response=="M")&(trl.odour1==trl.roles.odour1(r)))));
        trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).M = 1 - trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).H;
        trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).CR = length(find((trl.response=="CR")&(trl.odour1==trl.roles.odour1(r)))) / (length(find((trl.response=="CR")&(trl.odour1==trl.roles.odour1(r)))) + length(find((trl.response=="FA")&(trl.odour1==trl.roles.odour1(r)))));
        trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).FA = 1 - trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).CR;
        if ((trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).H>0 && trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).H<1) && (trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).FA>0 && trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).FA<1))
            [trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).dprime,trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).bias] = dprime_simple(trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).H,trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).FA);
        else
            trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).dprime = NaN;
            trl.performance.(char(strcat('odour',trl.roles.odour1(r)))).bias = NaN;
        end
    end

    % by second odour
    for r=1:length(trl.roles.odour2)
        trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).correct = ( length(find((trl.response=="H")&(trl.odour2==trl.roles.odour2(r)))) + length(find((trl.response=="CR")&(trl.odour2==trl.roles.odour2(r)))) ) / ( length(find((trl.response=="H")&(trl.odour2==trl.roles.odour2(r)))) + length(find((trl.response=="CR")&(trl.odour2==trl.roles.odour2(r)))) + length(find((trl.response=="M")&(trl.odour2==trl.roles.odour2(r)))) + length(find((trl.response=="FA")&(trl.odour2==trl.roles.odour2(r)))) );
        trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).incorrect = 1 - trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).correct;
        trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).H = length(find((trl.response=="H")&(trl.odour2==trl.roles.odour2(r)))) / (length(find((trl.response=="H")&(trl.odour2==trl.roles.odour2(r)))) + length(find((trl.response=="M")&(trl.odour2==trl.roles.odour2(r)))));
        trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).M = 1 - trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).H;
        trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).CR = length(find((trl.response=="CR")&(trl.odour2==trl.roles.odour2(r)))) / (length(find((trl.response=="CR")&(trl.odour2==trl.roles.odour2(r)))) + length(find((trl.response=="FA")&(trl.odour2==trl.roles.odour2(r)))));
        trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).FA = 1 - trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).CR;
        if ((trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).H>0 && trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).H<1) && (trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).FA>0 && trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).FA<1))
            [trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).dprime,trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).bias] = dprime_simple(trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).H,trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).FA);
        else
            trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).dprime = NaN;
            trl.performance.(char(strcat('odour',trl.roles.odour2(r)))).bias = NaN;
        end
    end

% end

%% Calculate performance metrics by block

trl.blocks.trialsPerBlock = p.trialsPerBlock;
trl.blocks.numBlocks = floor(trl.numTrials / p.trialsPerBlock);

for b=1:trl.blocks.numBlocks
    
    % general
    trl.performance.blocks_general.correct(b) = (length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")) + length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR"))) / (length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")) + length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")) + length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="M")) + length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="FA")));
    trl.performance.blocks_general.incorrect(b) = 1 - trl.performance.blocks_general.correct(b);
    trl.performance.blocks_general.H(b) = length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")) / (length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")) + length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="M")));
    trl.performance.blocks_general.M(b) = 1 - trl.performance.blocks_general.H(b);
    trl.performance.blocks_general.CR(b) = length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")) / (length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")) + length(find(trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="FA")));
    trl.performance.blocks_general.FA(b) = 1 - trl.performance.blocks_general.CR(b);
    if ((trl.performance.blocks_general.H(b)>0 && trl.performance.blocks_general.H(b)<1) && (trl.performance.blocks_general.FA(b)>0 && trl.performance.blocks_general.FA(b)<1))
        [trl.performance.blocks_general.dprime(b),trl.performance.blocks_general.bias(b)] = dprime_simple(trl.performance.blocks_general.H(b),trl.performance.blocks_general.FA(b));
    else
        trl.performance.blocks_general.dprime(b) = NaN;
        trl.performance.blocks_general.bias(b) = NaN;
    end
    
    % by trial type
    for s=1:length(trl.contingencies.stim)
        if strcmp(trl.contingencies.requirement(trl.contingencies.stim(s)),"GO")
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).H(b) = length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")&(trl.stim((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.contingencies.stim(s)))) / (length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")&(trl.stim((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.contingencies.stim(s)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="M")&(trl.stim((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.contingencies.stim(s)))));
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).M(b) = 1 - trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).H(b);
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).CR(b) = NaN;
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).FA(b) = NaN;
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).correct(b) = trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).H(b);
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).incorrect(b) = 1 - trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).correct(b);
        elseif strcmp(trl.contingencies.requirement(trl.contingencies.stim(s)),"NOGO")
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).H(b) = NaN;
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).M(b) = NaN;
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).CR(b) = length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")&(trl.stim((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.contingencies.stim(s)))) / (length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")&(trl.stim((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.contingencies.stim(s)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="FA")&(trl.stim((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.contingencies.stim(s)))));
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).FA(b) = 1 - trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).CR(b);
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).correct(b) = trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).CR(b);
            trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).incorrect(b) = 1 - trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).correct(b);        
        end
        trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).dprime(b) = NaN;
        trl.performance.(char(strcat('blocks_type',num2str(trl.contingencies.stim(s))))).bias(b) = NaN;
    end
    
    % by first odour
    for r=1:length(trl.roles.odour1)
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).correct(b) = ( length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))) ) / ( length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="M")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="FA")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))) );
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).incorrect(b) = 1 - trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).correct(b);
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).H(b) = length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))) / (length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="M")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))));
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).M(b) = 1 - trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).H(b);
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).CR(b) = length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))) / (length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="FA")&(trl.odour1((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour1(r)))));
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).FA(b) = 1 - trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).CR(b);
        if ((trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).H(b)>0 && trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).H(b)<1) && (trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).FA(b)>0 && trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).FA(b)<1))
            [trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).dprime(b),trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).bias(b)] = dprime_simple(trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).H(b),trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).FA(b));
        else
            trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).dprime(b) = NaN;
            trl.performance.(char(strcat('blocks_odour',trl.roles.odour1(r)))).bias(b) = NaN;
        end
    end

    % by second odour
    for r=1:length(trl.roles.odour2)
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).correct(b) = ( length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))) ) / ( length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="M")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="FA")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))) );
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).incorrect(b) = 1 - trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).correct(b);
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).H(b) = length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))) / (length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="H")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="M")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))));
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).M(b) = 1 - trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).H(b);
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).CR(b) = length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))) / (length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="CR")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))) + length(find((trl.response((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)=="FA")&(trl.odour2((b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock)==trl.roles.odour2(r)))));
        trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).FA(b) = 1 - trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).CR(b);
        if ((trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).H(b)>0 && trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).H(b)<1) && (trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).FA(b)>0 && trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).FA(b)<1))
            [trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).dprime(b),trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).bias(b)] = dprime_simple(trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).H(b),trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).FA(b));
        else
            trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).dprime(b) = NaN;
            trl.performance.(char(strcat('blocks_odour',trl.roles.odour2(r)))).bias(b) = NaN;
        end
    end
end

%% Calculate local performance metrics

% numWindows = length(p.localWindow);
% for c=1:numWindows
%     for b=1:trl.numTrials-(p.localWindow(c) -1 )
% 
%         % general
%         trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).correct(b) = (length(find(trl.response(b:b+p.localWindow(c)-1)=="H")) + length(find(trl.response(b:b+p.localWindow(c)-1)=="CR"))) / (length(find(trl.response(b:b+p.localWindow(c)-1)=="H")) + length(find(trl.response(b:b+p.localWindow(c)-1)=="CR")) + length(find(trl.response(b:b+p.localWindow(c)-1)=="M")) + length(find(trl.response(b:b+p.localWindow(c)-1)=="FA")));
%         trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).incorrect(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).correct(b);
%         trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).H(b) = length(find(trl.response(b:b+p.localWindow(c)-1)=="H")) / (length(find(trl.response(b:b+p.localWindow(c)-1)=="H")) + length(find(trl.response(b:b+p.localWindow(c)-1)=="M")));
%         trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).M(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).H(b);
%         trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).CR(b) = length(find(trl.response(b:b+p.localWindow(c)-1)=="CR")) / (length(find(trl.response(b:b+p.localWindow(c)-1)=="CR")) + length(find(trl.response(b:b+p.localWindow(c)-1)=="FA")));
%         trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).FA(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).CR(b);
%         if ((trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).H(b)>0 && trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).H(b)<1) && (trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).FA(b)>0 && trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).FA(b)<1))
%             [trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).dprime(b),trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).bias(b)] = dprime_simple(trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).H(b),trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).FA(b));
%         else
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).dprime(b) = NaN;
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_general'))).bias(b) = NaN;
%         end
% 
%         % by trial type
%         for s=1:length(trl.contingencies.stim)
%             if strcmp(trl.contingencies.requirement(trl.contingencies.stim(s)),"GO")
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).H(b) = length(find((trl.response(b:b+p.localWindow(c)-1)=="H")&(trl.stim(b:b+p.localWindow(c)-1)==trl.contingencies.stim(s)))) / (length(find((trl.response(b:b+p.localWindow(c)-1)=="H")&(trl.stim(b:b+p.localWindow(c)-1)==trl.contingencies.stim(s)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="M")&(trl.stim(b:b+p.localWindow(c)-1)==trl.contingencies.stim(s)))));
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).M(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).H(b);
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).CR(b) = NaN;
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).FA(b) = NaN;
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).correct(b) = trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).H(b);
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).incorrect(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).correct(b);
%             elseif strcmp(trl.contingencies.requirement(trl.contingencies.stim(s)),"NOGO")
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).H(b) = NaN;
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).M(b) = NaN;
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).CR(b) = length(find((trl.response(b:b+p.localWindow(c)-1)=="CR")&(trl.stim(b:b+p.localWindow(c)-1)==trl.contingencies.stim(s)))) / (length(find((trl.response(b:b+p.localWindow(c)-1)=="CR")&(trl.stim(b:b+p.localWindow(c)-1)==trl.contingencies.stim(s)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="FA")&(trl.stim(b:b+p.localWindow(c)-1)==trl.contingencies.stim(s)))));
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).FA(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).CR(b);
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).correct(b) = trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).CR(b);
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).incorrect(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).correct(b);        
%             end
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).dprime(b) = NaN;
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_type',num2str(trl.contingencies.stim(s))))).bias(b) = NaN;
%         end
% 
%         % by first odour
%         for r=1:length(trl.roles.odour1)
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).correct(b) = ( length(find((trl.response(b:b+p.localWindow(c)-1)=="H")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="CR")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))) ) / ( length(find((trl.response(b:b+p.localWindow(c)-1)=="H")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="CR")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="M")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="FA")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))) );
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).incorrect(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).correct(b);
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).H(b) = length(find((trl.response(b:b+p.localWindow(c)-1)=="H")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))) / (length(find((trl.response(b:b+p.localWindow(c)-1)=="H")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="M")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))));
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).M(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).H(b);
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).CR(b) = length(find((trl.response(b:b+p.localWindow(c)-1)=="CR")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))) / (length(find((trl.response(b:b+p.localWindow(c)-1)=="CR")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="FA")&(trl.odour1(b:b+p.localWindow(c)-1)==trl.roles.odour1(r)))));
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).FA(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).CR(b);
%             if ((trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).H(b)>0 && trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).H(b)<1) && (trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).FA(b)>0 && trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).FA(b)<1))
%                 [trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).dprime(b),trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).bias(b)] = dprime_simple(trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).H(b),trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).FA(b));
%             else
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).dprime(b) = NaN;
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).bias(b) = NaN;
%             end
%         end
% 
%         % by second odour
%         for r=1:length(trl.roles.odour2)
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).correct(b) = ( length(find((trl.response(b:b+p.localWindow(c)-1)=="H")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="CR")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))) ) / ( length(find((trl.response(b:b+p.localWindow(c)-1)=="H")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="CR")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="M")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="FA")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))) );
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).incorrect(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).correct(b);
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).H(b) = length(find((trl.response(b:b+p.localWindow(c)-1)=="H")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))) / (length(find((trl.response(b:b+p.localWindow(c)-1)=="H")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="M")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))));
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).M(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).H(b);
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).CR(b) = length(find((trl.response(b:b+p.localWindow(c)-1)=="CR")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))) / (length(find((trl.response(b:b+p.localWindow(c)-1)=="CR")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))) + length(find((trl.response(b:b+p.localWindow(c)-1)=="FA")&(trl.odour2(b:b+p.localWindow(c)-1)==trl.roles.odour2(r)))));
%             trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).FA(b) = 1 - trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).CR(b);
%             if ((trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).H(b)>0 && trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).H(b)<1) && (trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).FA(b)>0 && trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).FA(b)<1))
%                 [trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).dprime(b),trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).bias(b)] = dprime_simple(trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).H(b),trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).FA(b));
%             else
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).dprime(b) = NaN;
%                 trl.performance.(char(strcat('local',num2str(p.localWindow(c)),'_odour',trl.roles.odour1(r)))).bias(b) = NaN;
%             end
%         end
%     end
% end

%% Global performance metrics during blocks with good performance only

% trl.engaged.criterion_FA = p.criterion_FA;
% trl.engaged.criterion_M = p.criterion_M;
% 
% trl.engaged.includedBlocks = [];
% trl.engaged.includedTrials = [];
% for b=1:trl.blocks.numBlocks
%     if ((trl.performance.blocks_general.FA(b) <= trl.engaged.criterion_FA) & ...
%             (trl.performance.blocks_general.M(b) <= trl.engaged.criterion_M))
%         trl.engaged.includedBlocks = [trl.engaged.includedBlocks,b];
%         trl.engaged.includedTrials = [trl.engaged.includedTrials,(b-1)*trl.blocks.trialsPerBlock+1:b*trl.blocks.trialsPerBlock];
%     end
% end
% 
% % general
% trl.performance.engaged_general.correct = (length(find(trl.response(trl.engaged.includedTrials)=="H")) + length(find(trl.response(trl.engaged.includedTrials)=="CR"))) / (length(find(trl.response(trl.engaged.includedTrials)=="H")) + length(find(trl.response(trl.engaged.includedTrials)=="CR")) + length(find(trl.response(trl.engaged.includedTrials)=="M")) + length(find(trl.response(trl.engaged.includedTrials)=="FA")));
% trl.performance.engaged_general.incorrect = 1 - trl.performance.engaged_general.correct;
% trl.performance.engaged_general.H = length(find(trl.response(trl.engaged.includedTrials)=="H")) / (length(find(trl.response(trl.engaged.includedTrials)=="H")) + length(find(trl.response(trl.engaged.includedTrials)=="M")));
% trl.performance.engaged_general.M = 1 - trl.performance.engaged_general.H;
% trl.performance.engaged_general.CR = length(find(trl.response(trl.engaged.includedTrials)=="CR")) / (length(find(trl.response(trl.engaged.includedTrials)=="CR")) + length(find(trl.response(trl.engaged.includedTrials)=="FA")));
% trl.performance.engaged_general.FA = 1 - trl.performance.engaged_general.CR;
% if ((trl.performance.engaged_general.H>0 && trl.performance.engaged_general.H<1) && (trl.performance.engaged_general.FA>0 && trl.performance.engaged_general.FA<1))
%     [trl.performance.engaged_general.dprime,trl.performance.engaged_general.bias] = dprime_simple(trl.performance.engaged_general.H,trl.performance.engaged_general.FA);
% else
%     trl.performance.engaged_general.dprime = NaN;
%     trl.performance.engaged_general.bias = NaN;
% end
% 
% % by trial type
% for s=1:length(trl.contingencies.stim)
%     if strcmp(trl.contingencies.requirement(trl.contingencies.stim(s)),"GO")
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).H = length(find((trl.response(trl.engaged.includedTrials)=="H")&(trl.stim(trl.engaged.includedTrials)==trl.contingencies.stim(s)))) / (length(find((trl.response(trl.engaged.includedTrials)=="H")&(trl.stim(trl.engaged.includedTrials)==trl.contingencies.stim(s)))) + length(find((trl.response(trl.engaged.includedTrials)=="M")&(trl.stim(trl.engaged.includedTrials)==trl.contingencies.stim(s)))));
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).M = 1 - trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).H;
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).CR = NaN;
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).FA = NaN;
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).correct = trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).H;
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).incorrect = 1 - trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).correct;
%     elseif strcmp(trl.contingencies.requirement(trl.contingencies.stim(s)),"NOGO")
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).H = NaN;
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).M = NaN;
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).CR = length(find((trl.response(trl.engaged.includedTrials)=="CR")&(trl.stim(trl.engaged.includedTrials)==trl.contingencies.stim(s)))) / (length(find((trl.response(trl.engaged.includedTrials)=="CR")&(trl.stim(trl.engaged.includedTrials)==trl.contingencies.stim(s)))) + length(find((trl.response(trl.engaged.includedTrials)=="FA")&(trl.stim(trl.engaged.includedTrials)==trl.contingencies.stim(s)))));
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).FA = 1 - trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).CR;
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).correct = trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).CR;
%         trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).incorrect = 1 - trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).correct;        
%     end
%     trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).dprime = NaN;
%     trl.performance.(char(strcat('engaged_type',num2str(trl.contingencies.stim(s))))).bias = NaN;
% end
% 
% % by first odour
% for r=1:length(trl.roles.odour1)
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).correct = ( length(find((trl.response(trl.engaged.includedTrials)=="H")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="CR")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))) ) / ( length(find((trl.response(trl.engaged.includedTrials)=="H")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="CR")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="M")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="FA")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))) );
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).incorrect = 1 - trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).correct;
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).H = length(find((trl.response(trl.engaged.includedTrials)=="H")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))) / (length(find((trl.response(trl.engaged.includedTrials)=="H")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="M")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))));
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).M = 1 - trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).H;
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).CR = length(find((trl.response(trl.engaged.includedTrials)=="CR")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))) / (length(find((trl.response(trl.engaged.includedTrials)=="CR")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="FA")&(trl.odour1(trl.engaged.includedTrials)==trl.roles.odour1(r)))));
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).FA = 1 - trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).CR;
%     if ((trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).H>0 && trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).H<1) && (trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).FA>0 && trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).FA<1))
%         [trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).dprime,trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).bias] = dprime_simple(trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).H,trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).FA);
%     else
%         trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).dprime = NaN;
%         trl.performance.(char(strcat('engaged_odour',trl.roles.odour1(r)))).bias = NaN;
%     end
% end
% 
% % by second odour
% for r=1:length(trl.roles.odour2)
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).correct = ( length(find((trl.response(trl.engaged.includedTrials)=="H")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="CR")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))) ) / ( length(find((trl.response(trl.engaged.includedTrials)=="H")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="CR")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="M")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="FA")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))) );
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).incorrect = 1 - trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).correct;
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).H = length(find((trl.response(trl.engaged.includedTrials)=="H")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))) / (length(find((trl.response(trl.engaged.includedTrials)=="H")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="M")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))));
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).M = 1 - trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).H;
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).CR = length(find((trl.response(trl.engaged.includedTrials)=="CR")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))) / (length(find((trl.response(trl.engaged.includedTrials)=="CR")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))) + length(find((trl.response(trl.engaged.includedTrials)=="FA")&(trl.odour2(trl.engaged.includedTrials)==trl.roles.odour2(r)))));
%     trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).FA = 1 - trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).CR;
%     if ((trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).H>0 && trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).H<1) && (trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).FA>0 && trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).FA<1))
%         [trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).dprime,trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).bias] = dprime_simple(trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).H,trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).FA);
%     else
%         trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).dprime = NaN;
%         trl.performance.(char(strcat('engaged_odour',trl.roles.odour2(r)))).bias = NaN;
%     end
% end

%% Save mat files

if ~exist(fileparts(savePath),'dir')
   mkdir(fileparts(savePath));
end
save(savePath,'trl','-v7.3');

end





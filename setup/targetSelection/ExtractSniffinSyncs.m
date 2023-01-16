%% Extract sniffinSyncs

clc; clear;

info.animal             = 'Python'; % 'Python';
info.date               = '20211210'; % '20211210';
info.rootPath           = 'N:\Data\SniffinHippo\'; % 'C:\Data\SniffinHippo\';


%%% --------------------------- %%%
%%% ------ DON'T CHANGE ------- %%%
%%% --------------------------- %%%

p.sync.voltageThreshold                 = 0.5; % [V], threshold for digitising
info.paq.config.imaging                 = 7;
info.paq.config.cue                     = 8;


%% Load data

disp(['- Extract sniffinSyncs.'])

disp(['- Loading paq data.'])

info.year = info.date(1:4);
info.month = info.date(5:6);
info.day = info.date(7:8);
path.homeFolder = [info.rootPath,info.year,'\',info.year,'-',info.month,'\',info.year,'-',info.month,'-',info.day,'\',info.animal,'\','Behaviour','\'];
path.syncPath = [path.homeFolder,info.animal,'_',info.date,'_base1.paq'];

paq = paq2lab_ov20220105(path.syncPath); 
ts.cue = detectThresholdCrossing_ov20220105(paq(:,info.paq.config.cue),'above',p.sync.voltageThreshold);
ts.imaging = detectThresholdCrossing_ov20220105(paq(:,info.paq.config.imaging),'above',p.sync.voltageThreshold);
temp = [-Inf, ts.imaging(2:end)', Inf]; % even if an event starts at the very end of a frame, it's still counted to the same frame
if any(isnan(temp))
    temp2 = find(isnan(temp));
    temp(temp2) = linspace(temp(temp2(1)-1),temp(temp2(end)+1),length(temp2));
end
sniffinSyncs = discretize(ts.cue', temp);


%% Save data

disp(['- Saving sniffinSyncs file.'])
save([path.homeFolder,info.animal,'_',info.date,'_base1_sniffinSyncs.mat'],'sniffinSyncs','-v7.3');

disp('- Done.')







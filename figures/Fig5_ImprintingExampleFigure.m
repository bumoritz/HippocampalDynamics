%% Imprinting Example Figure

path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig5_fig\'];
save_root_png = [path.root_summary,'figures\Fig5_png\'];
save_root_pdf = [path.root_summary,'figures\Fig5_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig5_txt\'];


%% Make single cell figure

load('F:\SniffinHippo\AnalysisX\Turing\Turing_20220916\Turing_20220916_tng_all_stimVersion.mat');
load('F:\SniffinHippo\AnalysisX\Kane\Kane_20210912\Kane_20210912_tng_all_stimVersion.mat');
load('F:\SniffinHippo\AnalysisX\Philip\Philip_20211004\Philip_20211004_tng_all_stimVersion.mat');
load('F:\SniffinHippo\AnalysisX\Austin\Austin_20220729\Austin_20220729_tng_all_stimVersion.mat');
load('F:\SniffinHippo\AnalysisX\Celo\Celo_20210706\Celo_20210706_tng_all_stimVersion.mat');
load('F:\SniffinHippo\AnalysisX\Faramir\Faramir_20210608\Faramir_20210608_tng_all_stimVersion.mat');
load('F:\SniffinHippo\AnalysisX\Jobs\Jobs_20220721\Jobs_20220721_tng_all_stimVersion.mat');
load('F:\SniffinHippo\AnalysisX\Stanage\Stanage_20210925\Stanage_20210925_tng_all_stimVersion.mat');


%%

idx = 1185; %67; %381; %91; % 425; 587;

tng_all.prop = prop;
tng_all.trials = trials_all_stimVersion;
tng_all.traces = traces_all_stimVersion;
tng_all.avgTraces = avgTraces_all_stimVersion;
tng_all.firingField = tng_all_stimVersion.firingField;

numBlocks = 40; %38;
tng_all.trials.stimuli_inA_var0.AB = tng_all.trials.stimuli_inA_var0.AB(1:1*numBlocks);
tng_all.trials.stimuli_inA_var1.AB = tng_all.trials.stimuli_inA_var1.AB(1:4*numBlocks);
tng_all.trials.stimuli_inA_var0.AY = tng_all.trials.stimuli_inA_var0.AY(1:1*numBlocks);
tng_all.trials.stimuli_inA_var1.AY = tng_all.trials.stimuli_inA_var1.AY(1:4*numBlocks);
tng_all.trials.stimuli_inX_var0.XY = tng_all.trials.stimuli_inX_var0.XY(1:1*numBlocks);
tng_all.trials.stimuli_inX_var1.XY = tng_all.trials.stimuli_inX_var1.XY(1:4*numBlocks);
tng_all.trials.stimuli_inX_var0.XB = tng_all.trials.stimuli_inX_var0.XB(1:1*numBlocks);
tng_all.trials.stimuli_inX_var1.XB = tng_all.trials.stimuli_inX_var1.XB(1:4*numBlocks);

limits = [];

% IN
%limits = [...]; % Stanage_20210925_idx532, A
%limits = [...]; % Stanage_20210925_idx420, A
%limits = [-2.3739,1.7940,-0.2919,3.1980]; % Philip_20211004_idx321, A, cluster_position = 14
%limits = [...]; % Kane_20210912_idx173, X
%limits = [-2.8997,2.1720,-0.5126,3.8194]; % Philip_20211004_idx134, X, cluster_position = 2
%limits = [...]; % Faramir_20210608_idx60, X

% MAYBE
%limits = [...]; % Celo_20210706_idx235, X

% PROBABLY NOT
%limits = [...]; % Stanage_20210925_idx37, A
%limits = [...]; % Jobs_20220721_idx587, X
%limits = [-2.6918,2.1667,-0.4571,5.5480]; % Kane_20210912_idx67, A
%limits = [-2.7563,2.2162,-0.4346,4.2032]; % Turing_20220916_idx91, A
%limits = [...]; % Austin_20220729_idx13, A

% good middle cells
% Philip 134
% Philip 321 (almost late)

% good late cells
% Faramir 60
% Jobs 587

%limits = [-2.6638,2.1072,-0.3080,2.7333]; % Kane_20210912_idx381, X

% F = singleCellFigure(info,p,tng_all,idx,paq_beh,'_var1',limits); % '_var0' or '_var1'
% title([info.animal,'-',info.date,'-idx',num2str(idx),'-stim'])
% 
% F = singleCellFigure(info,p,tng_all,idx,paq_beh,'_var0',limits); % '_var0' or '_var1'
% title([info.animal,'-',info.date,'-idx',num2str(idx),'-nostim'])

%title([info.animal,'-',info.date,'-idx',num2str(idx),'-nostim'])

%F = singleCellFigure(info,p,tng_all,idx,paq_beh,'_var0',limits,'A_X'); % VAR0 IS HARDCODED NOW
F = singleCellFigure(info,p,tng_all,idx,paq_beh,'_var0',limits,'X_A'); % VAR0 IS HARDCODED NOW

% cluster_position = 2;
% 
% stim_duration = 0.09;
% this_t = p.general.t_binned;
% temp = 0.3:0.25:5.3;   
% patch([interp1(this_t,1:length(this_t),temp(cluster_position)),interp1(this_t,1:length(this_t),temp(cluster_position)),interp1(this_t,1:length(this_t),temp(cluster_position)+stim_duration),interp1(this_t,1:length(this_t),temp(cluster_position)+stim_duration)],...
%         [0,1,1,0],p.col.photostim,'EdgeColor','none');

    
% %%

set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');

title(num2str(idx))

% savefig(F,[save_root_fig,'\Fig5_ImprintingAnalysis_ExampleCell_',num2str(idx),'.fig']);
% saveas(F,[save_root_png,'\Fig5_ImprintingAnalysis_ExampleCell_',num2str(idx),'.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_ImprintingAnalysis_ExampleCell_',num2str(idx),'.pdf']); set(gcf,'Color',[1,1,1])


%%


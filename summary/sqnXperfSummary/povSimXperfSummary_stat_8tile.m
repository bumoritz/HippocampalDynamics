%% Select data

this_field = 'statTrialwiseBalanced'; % 'statTrialwise','statTrialwiseBalanced'
this_metric = 'Pearson'; % 'Pearson','cosTheta'
these_idcs = 'iscells'; % 'iscells','A','X'
this_window = '1st2s'; % '1st1s','1st2s','1st3s'
this_minNumTrials = 5; % 1, 3, 5, 10

nrows = 2; ncols = 4;
F = default_figure([20,0.5,20,9.9]);

if strcmp(this_field,'statTrialwiseBalanced')
    try
        conditions = {'A_corr_Otemp_A_H_M','A_corr_Otemp_A_M_M','A_corr_Otemp_A_CR_M','A_corr_Otemp_A_FA_M','A_corr_Otemp_X_H_M','A_corr_Otemp_X_M_M','A_corr_Otemp_X_CR_M','A_corr_Otemp_X_FA_M',...
        'X_corr_Otemp_A_H_M','X_corr_Otemp_A_M_M','X_corr_Otemp_A_CR_M','X_corr_Otemp_A_FA_M','X_corr_Otemp_X_H_M','X_corr_Otemp_X_M_M','X_corr_Otemp_X_CR_M','X_corr_Otemp_X_FA_M'};
    catch
        conditions = {'Atemp_A_H_M','Atemp_A_M_M','Atemp_A_CR_M','Atemp_A_FA_M','Atemp_X_H_M','Atemp_X_M_M','Atemp_X_CR_M','Atemp_X_FA_M',...
        'Xtemp_A_H_M','Xtemp_A_M_M','Xtemp_A_CR_M','Xtemp_A_FA_M','Xtemp_X_H_M','Xtemp_X_M_M','Xtemp_X_CR_M','Xtemp_X_FA_M'};
        warning('Did not find version based on all trials (as opposed to corr_O trials')
    end
elseif strcmp(this_field,'statTrialwise')
    conditions = {'Atemp_A_H','Atemp_A_M','Atemp_A_CR','Atemp_A_FA','Atemp_X_H','Atemp_X_M','Atemp_X_CR','Atemp_X_FA',...
        'Xtemp_A_H','Xtemp_A_M','Xtemp_A_CR','Xtemp_A_FA','Xtemp_X_H','Xtemp_X_M','Xtemp_X_CR','Xtemp_X_FA'};
end
cols = {p.col.AB,mean([p.col.AB;p.col.white]),p.col.AY,mean([p.col.AY;p.col.white]),p.col.XY,mean([p.col.XY;p.col.white]),p.col.XB,mean([p.col.XB;p.col.white]),...
    p.col.AB,mean([p.col.AB;p.col.white]),p.col.AY,mean([p.col.AY;p.col.white]),p.col.XY,mean([p.col.XY;p.col.white]),p.col.XB,mean([p.col.XB;p.col.white])};


%% Extract data

povSim = {};
for k=1:length(conditions)
    povSim.(conditions{k}) = nan(d_info.numAnimals,1);
end
for i=1:d_info.numAnimals
    try    
        this_struct = d{i,1}.sqn_all.povSim.(this_field).(this_metric).([these_idcs,'_',this_window]);
        for k=1:length(conditions)
            if length(this_struct.(conditions{k}))>=this_minNumTrials
                povSim.(conditions{k})(i) = nanmean(this_struct.(conditions{k}));
            end
        end
    catch
    end
end
        

%% Make figure

for k=1:length(conditions)/2
    these_labels = {'correct','incorrect'};
    these_data = [povSim.(conditions{k*2-1}),povSim.(conditions{k*2})];
    these_cols = {cols{k*2-1},cols{k*2}};
    
    subplot(nrows,ncols,k)
    hold on
    barsWithLines(these_labels,these_data,these_cols);
    ylim([0,1])
    ylabel('Population vector correlation')
 
    this_p = signrank(these_data(:,1),these_data(:,2));
    temp = conditions{k*2}(7:end);
    if strcmp(temp,'A_H') | strcmp(temp,'A_M')
        temp = 'AB';
    elseif strcmp(temp,'A_CR') | strcmp(temp,'A_FA')
        temp = 'AY';
    elseif strcmp(temp,'X_H') | strcmp(temp,'X_M')
        temp = 'XY';
    elseif strcmp(temp,'X_CR') | strcmp(temp,'X_FA')
        temp = 'XB';
    end
    title([temp,' trials with ',conditions{k*2}(1:5),'. p=',num2str(this_p,2)])
end
suptitle([these_idcs,', ',this_window,', min. numTrials ',num2str(this_minNumTrials)])

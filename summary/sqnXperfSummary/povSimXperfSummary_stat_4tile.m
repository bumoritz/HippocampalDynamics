%% Select data

this_field = 'statTrialwiseBalanced'; % 'statTrialwise','statTrialwiseBalanced'
this_metric = 'Pearson'; % 'Pearson','cosTheta'
these_idcs = 'X'; % 'iscells','A','X'
this_window = '1st2s'; % '1st1s','1st2s','1st3s'
this_minNumTrials = 5; % 1, 3, 5, 10

nrows = 2; ncols = 2;
F = default_figure([20,0.5,20,9.9]);

conditions = {'ipsi_H','ipsi_M','ipsi_CR','ipsi_FA','contra_H','contra_M','contra_CR','contra_FA'};
cols = {p.col.darkGray,p.col.gray,p.col.darkGray,p.col.gray,p.col.darkGray,p.col.gray,p.col.darkGray,p.col.gray};


%% Extract data

povSim = {};
for k=1:length(conditions)
    povSim.(conditions{k}) = nan(d_info.numAnimals,1);
end
for i=1:d_info.numAnimals
    try    
        this_struct = d{i,1}.sqn_all.povSim.(this_field).(this_metric).([these_idcs,'_',this_window]);
        
        if strcmp(this_field,'statTrialwiseBalanced')
            temp1 = [this_struct.Atemp_A_H_M;this_struct.Xtemp_X_H_M];
            temp2 = [this_struct.Atemp_A_M_M;this_struct.Xtemp_X_M_M];
        elseif strcmp(this_field,'statTrialwise')
            temp1 = [this_struct.Atemp_A_H;this_struct.Xtemp_X_H];
            temp2 = [this_struct.Atemp_A_M;this_struct.Xtemp_X_M];
        end        
        if length(temp1)>=this_minNumTrials && length(temp2)>=this_minNumTrials
            povSim.ipsi_H(i) = nanmean(temp1);
            povSim.ipsi_M(i) = nanmean(temp2);
        end

        if strcmp(this_field,'statTrialwiseBalanced')
            temp1 = [this_struct.Atemp_A_CR_M;this_struct.Xtemp_X_CR_M];
            temp2 = [this_struct.Atemp_A_FA_M;this_struct.Xtemp_X_FA_M];
        elseif strcmp(this_field,'statTrialwise')
            temp1 = [this_struct.Atemp_A_CR;this_struct.Xtemp_X_CR];
            temp2 = [this_struct.Atemp_A_FA;this_struct.Xtemp_X_FA];
        end   
        if length(temp1)>=this_minNumTrials && length(temp2)>=this_minNumTrials
            povSim.ipsi_CR(i) = nanmean(temp1);
            povSim.ipsi_FA(i) = nanmean(temp2);
        end

        if strcmp(this_field,'statTrialwiseBalanced')
            temp1 = [this_struct.Xtemp_A_H_M;this_struct.Atemp_X_H_M];
            temp2 = [this_struct.Xtemp_A_M_M;this_struct.Atemp_X_M_M];
        elseif strcmp(this_field,'statTrialwise')
            temp1 = [this_struct.Xtemp_A_H;this_struct.Atemp_X_H];
            temp2 = [this_struct.Xtemp_A_M;this_struct.Atemp_X_M];
        end   
        if length(temp1)>=this_minNumTrials && length(temp2)>=this_minNumTrials
            povSim.contra_H(i) = nanmean(temp1);
            povSim.contra_M(i) = nanmean(temp2);
        end

        if strcmp(this_field,'statTrialwiseBalanced')
            temp1 = [this_struct.Xtemp_A_CR_M;this_struct.Atemp_X_CR_M];
            temp2 = [this_struct.Xtemp_A_FA_M;this_struct.Atemp_X_FA_M];
        elseif strcmp(this_field,'statTrialwise')
            temp1 = [this_struct.Xtemp_A_CR;this_struct.Atemp_X_CR];
            temp2 = [this_struct.Xtemp_A_FA;this_struct.Atemp_X_FA];
        end   
        if length(temp1)>=this_minNumTrials && length(temp2)>=this_minNumTrials
            povSim.contra_CR(i) = nanmean(temp1);
            povSim.contra_FA(i) = nanmean(temp2);
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
    if strcmp(conditions{k*2},'ipsi_M')
    	title(['Go trials with ipsi template. p=',num2str(this_p,2)])
    elseif strcmp(conditions{k*2},'ipsi_FA')
        title(['Nogo trials with ipsi template. p=',num2str(this_p,2)])
    elseif strcmp(conditions{k*2},'contra_M')
        title(['Go trials with contra template. p=',num2str(this_p,2)])
    elseif strcmp(conditions{k*2},'contra_FA')
        title(['Nogo trials with contra template. p=',num2str(this_p,2)])
    end
end
suptitle([these_idcs,', ',this_window,', min. numTrials ',num2str(this_minNumTrials)])

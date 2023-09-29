%% AW

this_metric = 'Pearson'; % 'Pearson','cosTheta'

nrows = 2; ncols = 4;
F = default_figure([20,0.5,20,9.9]);

these_conditions = {'Atemp_A_H','Atemp_A_M','Atemp_A_CR','Atemp_A_FA','Atemp_X_H','Atemp_X_M','Atemp_X_CR','Atemp_X_FA',...
    'Xtemp_A_H','Xtemp_A_M','Xtemp_A_CR','Xtemp_A_FA','Xtemp_X_H','Xtemp_X_M','Xtemp_X_CR','Xtemp_X_FA'};
these_cols = {p.col.AB,mean([p.col.AB;p.col.white]),p.col.AY,mean([p.col.AY;p.col.white]),p.col.XY,mean([p.col.XY;p.col.white]),p.col.XB,mean([p.col.XB;p.col.white]),...
    p.col.AB,mean([p.col.AB;p.col.white]),p.col.AY,mean([p.col.AY;p.col.white]),p.col.XY,mean([p.col.XY;p.col.white]),p.col.XB,mean([p.col.XB;p.col.white])};

this_data = {};
for k=1:length(these_conditions)
    this_data = [this_data, povSim.dynTrialwise.(this_metric).iscells_AW.(these_conditions{k})];
end

for k=1:length(these_conditions)/2
    subplot(nrows,ncols,k)
    hold on
    
    temp=shadedErrorBar(1:size(this_data{k*2-1},2),nanmean(this_data{k*2-1},1),nansem(this_data{k*2-1},1),'lineProps',these_cols{k*2-1}); temp.mainLine.LineWidth = 2;  
    temp=shadedErrorBar(1:size(this_data{k*2},2),nanmean(this_data{k*2},1),nansem(this_data{k*2},1),'lineProps',these_cols{k*2}); temp.mainLine.LineWidth = 2;  

    ylim([0,0.6])
    %taskLines(p,info,'AW','traces')
    ylabel('Population vector correlation')
    temp = these_conditions{k*2}(7:end);
    if strcmp(temp,'A_H') | strcmp(temp,'A_M')
        temp = 'AB';
    elseif strcmp(temp,'A_CR') | strcmp(temp,'A_FA')
        temp = 'AY';
    elseif strcmp(temp,'X_H') | strcmp(temp,'X_M')
        temp = 'XY';
    elseif strcmp(temp,'X_CR') | strcmp(temp,'X_FA')
        temp = 'XB';
    end
    title([temp,' trials with ',these_conditions{k*2}(1:5),'. n=',num2str(size(this_data{k*2},1)),' error trials.'])
end


%% full

this_metric = 'Pearson'; % 'Pearson','cosTheta'

nrows = 4; ncols = 4;
F = default_figure([20,0.5,20,9.9]);

these_conditions = {'ABtemp_A_H','ABtemp_A_M','ABtemp_A_CR','ABtemp_A_FA','ABtemp_X_H','ABtemp_X_M','ABtemp_X_CR','ABtemp_X_FA',...
    'AYtemp_A_H','AYtemp_A_M','AYtemp_A_CR','AYtemp_A_FA','AYtemp_X_H','AYtemp_X_M','AYtemp_X_CR','AYtemp_X_FA',...
    'XYtemp_A_H','XYtemp_A_M','XYtemp_A_CR','XYtemp_A_FA','XYtemp_X_H','XYtemp_X_M','XYtemp_X_CR','XYtemp_X_FA',...
    'XBtemp_A_H','XBtemp_A_M','XBtemp_A_CR','XBtemp_A_FA','XBtemp_X_H','XBtemp_X_M','XBtemp_X_CR','XBtemp_X_FA'};
these_cols = {p.col.AB,mean([p.col.AB;p.col.white]),p.col.AY,mean([p.col.AY;p.col.white]),p.col.XY,mean([p.col.XY;p.col.white]),p.col.XB,mean([p.col.XB;p.col.white]),...
    p.col.AB,mean([p.col.AB;p.col.white]),p.col.AY,mean([p.col.AY;p.col.white]),p.col.XY,mean([p.col.XY;p.col.white]),p.col.XB,mean([p.col.XB;p.col.white]),...
    p.col.AB,mean([p.col.AB;p.col.white]),p.col.AY,mean([p.col.AY;p.col.white]),p.col.XY,mean([p.col.XY;p.col.white]),p.col.XB,mean([p.col.XB;p.col.white]),...
    p.col.AB,mean([p.col.AB;p.col.white]),p.col.AY,mean([p.col.AY;p.col.white]),p.col.XY,mean([p.col.XY;p.col.white]),p.col.XB,mean([p.col.XB;p.col.white])};

this_data = {};
for k=1:length(these_conditions)
    this_data = [this_data, povSim.dynTrialwise.(this_metric).iscells_full.(these_conditions{k})];
end

for k=1:length(these_conditions)/2
    subplot(nrows,ncols,k)
    hold on

    temp=shadedErrorBar(1:size(this_data{k*2-1},2),nanmean(this_data{k*2-1},1),nansem(this_data{k*2-1},1),'lineProps',these_cols{k*2-1}); temp.mainLine.LineWidth = 2;  
    temp=shadedErrorBar(1:size(this_data{k*2},2),nanmean(this_data{k*2},1),nansem(this_data{k*2},1),'lineProps',these_cols{k*2}); temp.mainLine.LineWidth = 2;  
    
    ylim([0,0.6])
    taskLines(p,info,'full','traces');
    ylabel('Population vector correlation')
    temp = these_conditions{k*2}(8:end);
    if strcmp(temp,'A_H') | strcmp(temp,'A_M')
        temp = 'AB';
    elseif strcmp(temp,'A_CR') | strcmp(temp,'A_FA')
        temp = 'AY';
    elseif strcmp(temp,'X_H') | strcmp(temp,'X_M')
        temp = 'XY';
    elseif strcmp(temp,'X_CR') | strcmp(temp,'X_FA')
        temp = 'XB';
    end
    title([temp,' trials with ',these_conditions{k*2}(1:6),'. n=',num2str(size(this_data{k*2},1)),' error trials.'])
end


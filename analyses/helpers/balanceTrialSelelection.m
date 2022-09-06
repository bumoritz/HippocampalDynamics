function trials_all = balanceTrialSelelection(trials_all,info,task)
% M for matched, O for other

trials_all.balanced.A_H_M = [];
trials_all.balanced.A_H_O = [];
trials_all.balanced.A_M_M = [];
trials_all.balanced.A_M_O = [];
trials_all.balanced.A_CR_M = [];
trials_all.balanced.A_CR_O = [];
trials_all.balanced.A_FA_M = [];
trials_all.balanced.A_FA_O = [];
trials_all.balanced.X_H_M = [];
trials_all.balanced.X_H_O = [];
trials_all.balanced.X_M_M = [];
trials_all.balanced.X_M_O = [];
trials_all.balanced.X_CR_M = [];
trials_all.balanced.X_CR_O = [];
trials_all.balanced.X_FA_M = [];
trials_all.balanced.X_FA_O = [];
trials_all.balanced.A_corr_M = [];
trials_all.balanced.A_corr_O = [];
trials_all.balanced.X_corr_M = [];
trials_all.balanced.X_corr_O = [];
trials_all.balanced.B_corr_M = [];
trials_all.balanced.B_corr_O = [];
trials_all.balanced.Y_corr_M = [];
trials_all.balanced.Y_corr_O = [];
for i=1:info.task.numTrials/info.task.trialsPerBlock
    these_trials = (i-1)*info.task.trialsPerBlock+1:i*info.task.trialsPerBlock;

    % these trial types
    these_A = intersect( find(task.odour1=="A") ,these_trials);
    these_X = intersect( find(task.odour1=="X") ,these_trials);
    these_B = intersect( find(task.odour2=="B") ,these_trials);
    these_Y = intersect( find(task.odour2=="Y") ,these_trials);
    these_H = intersect( find(task.response=="H") ,these_trials);
    these_M = intersect( find(task.response=="M") ,these_trials);
    these_CR = intersect( find(task.response=="CR") ,these_trials);
    these_FA = intersect( find(task.response=="FA") ,these_trials);
    these_AB = intersect(these_A,these_B);
	these_AY = intersect(these_A,these_Y);
	these_XY = intersect(these_X,these_Y);
	these_XB = intersect(these_X,these_B);

    % these trial types with outcome
    these_A_H = intersect(these_AB,these_H);
	these_A_CR = intersect(these_AY,these_CR);
	these_X_H = intersect(these_XY,these_H);
	these_X_CR = intersect(these_XB,these_CR);
    these_A_M = intersect(these_AB,these_M);
	these_A_FA = intersect(these_AY,these_FA);
	these_X_M = intersect(these_XY,these_M);
	these_X_FA = intersect(these_XB,these_FA);    
 
    % A-go trials
    this_numMatched = nanmin([length(these_A_H),length(these_A_M)]);
    temp = datasample(these_A_H,this_numMatched,'Replace',false);
    trials_all.balanced.A_H_M = [trials_all.balanced.A_H_M, temp];
    trials_all.balanced.A_H_O = [trials_all.balanced.A_H_O, setdiff(these_A_H,temp)];
    temp = datasample(these_A_M,this_numMatched,'Replace',false);    
    trials_all.balanced.A_M_M = [trials_all.balanced.A_M_M, temp];
    trials_all.balanced.A_M_O = [trials_all.balanced.A_M_O, setdiff(these_A_M,temp)];
    
    % X-go trials
    this_numMatched = nanmin([length(these_X_H),length(these_X_M)]);
    temp = datasample(these_X_H,this_numMatched,'Replace',false);
    trials_all.balanced.X_H_M = [trials_all.balanced.X_H_M, temp];
    trials_all.balanced.X_H_O = [trials_all.balanced.X_H_O, setdiff(these_X_H,temp)];
    temp = datasample(these_X_M,this_numMatched,'Replace',false);    
    trials_all.balanced.X_M_M = [trials_all.balanced.X_M_M, temp];
    trials_all.balanced.X_M_O = [trials_all.balanced.X_M_O, setdiff(these_X_M,temp)];
    
    % A-nogo trials
    this_numMatched = nanmin([length(these_A_CR),length(these_A_FA)]);
    temp = datasample(these_A_CR,this_numMatched,'Replace',false);
    trials_all.balanced.A_CR_M = [trials_all.balanced.A_CR_M, temp];
    trials_all.balanced.A_CR_O = [trials_all.balanced.A_CR_O, setdiff(these_A_CR,temp)];
    temp = datasample(these_A_FA,this_numMatched,'Replace',false);    
    trials_all.balanced.A_FA_M = [trials_all.balanced.A_FA_M, temp];
    trials_all.balanced.A_FA_O = [trials_all.balanced.A_FA_O, setdiff(these_A_FA,temp)];
    
    % X-nogo trials
    this_numMatched = nanmin([length(these_X_CR),length(these_X_FA)]);
    temp = datasample(these_X_CR,this_numMatched,'Replace',false);
    trials_all.balanced.X_CR_M = [trials_all.balanced.X_CR_M, temp];
    trials_all.balanced.X_CR_O = [trials_all.balanced.X_CR_O, setdiff(these_X_CR,temp)];
    temp = datasample(these_X_FA,this_numMatched,'Replace',false);    
    trials_all.balanced.X_FA_M = [trials_all.balanced.X_FA_M, temp];
    trials_all.balanced.X_FA_O = [trials_all.balanced.X_FA_O, setdiff(these_X_FA,temp)];
    
    % combined
    trials_all.balanced.A_corr_M = sort([trials_all.balanced.A_H_M,trials_all.balanced.A_CR_M]);
    trials_all.balanced.A_corr_O = sort([trials_all.balanced.A_H_O,trials_all.balanced.A_CR_O]);
    trials_all.balanced.X_corr_M = sort([trials_all.balanced.X_H_M,trials_all.balanced.X_CR_M]);
    trials_all.balanced.X_corr_O = sort([trials_all.balanced.X_H_O,trials_all.balanced.X_CR_O]);
    trials_all.balanced.B_corr_M = sort([trials_all.balanced.A_H_M,trials_all.balanced.X_CR_M]);
    trials_all.balanced.B_corr_O = sort([trials_all.balanced.A_H_O,trials_all.balanced.X_CR_O]);
    trials_all.balanced.Y_corr_M = sort([trials_all.balanced.X_H_M,trials_all.balanced.A_CR_M]);
    trials_all.balanced.Y_corr_O = sort([trials_all.balanced.X_H_O,trials_all.balanced.A_CR_O]);
end





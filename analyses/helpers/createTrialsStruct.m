function [trials] = createTrialsStruct(task,these_trials,do_by_stimType)
%these_trials = 1:info.task.numTrials, do_by_stimType = true;

if nargin<3
    do_by_stimType = false;
end


%% Trial types by stimuli

% odours
trials.stimuli.A = intersect( find(task.odour1=="A") ,these_trials);
trials.stimuli.X = intersect( find(task.odour1=="X") ,these_trials);
trials.stimuli.B = intersect( find(task.odour2=="B") ,these_trials);
trials.stimuli.Y = intersect( find(task.odour2=="Y") ,these_trials);

% odour combinations
[trials.stimuli.AB,trials.stimuli_inA.AB,trials.stimuli_inB.AB] = intersect(trials.stimuli.A,trials.stimuli.B);
[trials.stimuli.AY,trials.stimuli_inA.AY,trials.stimuli_inY.AY] = intersect(trials.stimuli.A,trials.stimuli.Y);
[trials.stimuli.XY,trials.stimuli_inX.XY,trials.stimuli_inY.XY] = intersect(trials.stimuli.X,trials.stimuli.Y);
[trials.stimuli.XB,trials.stimuli_inX.XB,trials.stimuli_inB.XB] = intersect(trials.stimuli.X,trials.stimuli.B);


%% Trial types by outcome

% outcomes
trials.outcome.H = intersect( find(task.response=="H") ,these_trials);
trials.outcome.M = intersect( find(task.response=="M") ,these_trials);
trials.outcome.CR = intersect( find(task.response=="CR") ,these_trials);
trials.outcome.FA = intersect( find(task.response=="FA") ,these_trials);
trials.outcome.auto = intersect( find(task.autoreward) ,these_trials);
trials.outcome.correct = sort([trials.outcome.H,trials.outcome.CR]);
trials.outcome.incorrect = sort([trials.outcome.M,trials.outcome.FA]);
trials.outcome.rew = sort([trials.outcome.H,trials.outcome.auto]);
trials.outcome.norew = sort([trials.outcome.M,trials.outcome.CR,trials.outcome.FA]);

% outcomes - trial type combinations (A)
[trials.outcome.A_H,trials.outcome_inA.A_H,trials.outcome_inH.A_H] = intersect(trials.stimuli.A,trials.outcome.H);
[trials.outcome.A_M,trials.outcome_inA.A_M,trials.outcome_inM.A_M] = intersect(trials.stimuli.A,trials.outcome.M);
[trials.outcome.A_CR,trials.outcome_inA.A_CR,trials.outcome_inCR.A_CR] = intersect(trials.stimuli.A,trials.outcome.CR);
[trials.outcome.A_FA,trials.outcome_inA.A_FA,trials.outcome_inFA.A_FA] = intersect(trials.stimuli.A,trials.outcome.FA);
[trials.outcome.A_correct,trials.outcome_inA.A_correct,trials.outcome_incorrect.A_correct] = intersect(trials.stimuli.A,trials.outcome.correct);
[trials.outcome.A_incorrect,trials.outcome_inA.A_incorrect,trials.outcome_inincorrect.A_incorrect] = intersect(trials.stimuli.A,trials.outcome.incorrect);
[trials.outcome.A_rew,trials.outcome_inA.A_rew,trials.outcome_inrew.A_rew] = intersect(trials.stimuli.A,trials.outcome.rew);
[trials.outcome.A_norew,trials.outcome_inA.A_norew,trials.outcome_inrew.A_norew] = intersect(trials.stimuli.A,trials.outcome.norew);

% outcomes - trial type combinations (X)
[trials.outcome.X_H,trials.outcome_inX.X_H,trials.outcome_inH.X_H] = intersect(trials.stimuli.X,trials.outcome.H);
[trials.outcome.X_M,trials.outcome_inX.X_M,trials.outcome_inM.X_M] = intersect(trials.stimuli.X,trials.outcome.M);
[trials.outcome.X_CR,trials.outcome_inX.X_CR,trials.outcome_inCR.X_CR] = intersect(trials.stimuli.X,trials.outcome.CR);
[trials.outcome.X_FA,trials.outcome_inX.X_FA,trials.outcome_inFA.X_FA] = intersect(trials.stimuli.X,trials.outcome.FA);
[trials.outcome.X_correct,trials.outcome_inX.X_correct,trials.outcome_incorrect.X_correct] = intersect(trials.stimuli.X,trials.outcome.correct);
[trials.outcome.X_incorrect,trials.outcome_inX.X_incorrect,trials.outcome_inincorrect.X_incorrect] = intersect(trials.stimuli.X,trials.outcome.incorrect);
[trials.outcome.X_rew,trials.outcome_inX.X_rew,trials.outcome_inrew.X_rew] = intersect(trials.stimuli.X,trials.outcome.rew);
[trials.outcome.X_norew,trials.outcome_inX.X_norew,trials.outcome_inrew.X_norew] = intersect(trials.stimuli.X,trials.outcome.norew);

% outcomes - trial type combinations (B)
[trials.outcome.B_H,trials.outcome_inB.B_H,trials.outcome_inH.B_H] = intersect(trials.stimuli.B,trials.outcome.H);
[trials.outcome.B_M,trials.outcome_inB.B_M,trials.outcome_inM.B_M] = intersect(trials.stimuli.B,trials.outcome.M);
[trials.outcome.B_CR,trials.outcome_inB.B_CR,trials.outcome_inCR.B_CR] = intersect(trials.stimuli.B,trials.outcome.CR);
[trials.outcome.B_FA,trials.outcome_inB.B_FA,trials.outcome_inFA.B_FA] = intersect(trials.stimuli.B,trials.outcome.FA);
[trials.outcome.B_correct,trials.outcome_inB.B_correct,trials.outcome_incorrect.B_correct] = intersect(trials.stimuli.B,trials.outcome.correct);
[trials.outcome.B_incorrect,trials.outcome_inB.B_incorrect,trials.outcome_inincorrect.B_incorrect] = intersect(trials.stimuli.B,trials.outcome.incorrect);
[trials.outcome.B_rew,trials.outcome_inB.B_rew,trials.outcome_inrew.B_rew] = intersect(trials.stimuli.B,trials.outcome.rew);
[trials.outcome.B_norew,trials.outcome_inB.B_norew,trials.outcome_inrew.B_norew] = intersect(trials.stimuli.B,trials.outcome.norew);

% outcomes - trial type combinations (Y)
[trials.outcome.Y_H,trials.outcome_inY.Y_H,trials.outcome_inH.Y_H] = intersect(trials.stimuli.Y,trials.outcome.H);
[trials.outcome.Y_M,trials.outcome_inY.Y_M,trials.outcome_inM.Y_M] = intersect(trials.stimuli.Y,trials.outcome.M);
[trials.outcome.Y_CR,trials.outcome_inY.Y_CR,trials.outcome_inCR.Y_CR] = intersect(trials.stimuli.Y,trials.outcome.CR);
[trials.outcome.Y_FA,trials.outcome_inY.Y_FA,trials.outcome_inFA.Y_FA] = intersect(trials.stimuli.Y,trials.outcome.FA);
[trials.outcome.Y_correct,trials.outcome_inY.Y_correct,trials.outcome_incorrect.Y_correct] = intersect(trials.stimuli.Y,trials.outcome.correct);
[trials.outcome.Y_incorrect,trials.outcome_inY.Y_incorrect,trials.outcome_inincorrect.Y_incorrect] = intersect(trials.stimuli.Y,trials.outcome.incorrect);
[trials.outcome.Y_rew,trials.outcome_inY.Y_rew,trials.outcome_inrew.Y_rew] = intersect(trials.stimuli.Y,trials.outcome.rew);
[trials.outcome.Y_norew,trials.outcome_inY.Y_norew,trials.outcome_inrew.Y_norew] = intersect(trials.stimuli.Y,trials.outcome.norew);


%% Trial types by stimuli - by stim type

if do_by_stimType

    % odours - var0
    trials.stimuli_var0.A = intersect( intersect( find(task.odour1=="A") ,these_trials) ,find(task.var==0));
    trials.stimuli_var0.X = intersect( intersect( find(task.odour1=="X") ,these_trials) ,find(task.var==0));
    trials.stimuli_var0.B = intersect( intersect( find(task.odour2=="B") ,these_trials) ,find(task.var==0));
    trials.stimuli_var0.Y = intersect( intersect( find(task.odour2=="Y") ,these_trials) ,find(task.var==0));

    % odour combinations - var0
    [trials.stimuli_var0.AB,trials.stimuli_inA_var0.AB,trials.stimuli_inB_var0.AB] = intersect(trials.stimuli_var0.A,trials.stimuli_var0.B);
    [trials.stimuli_var0.AY,trials.stimuli_inA_var0.AY,trials.stimuli_inY_var0.AY] = intersect(trials.stimuli_var0.A,trials.stimuli_var0.Y);
    [trials.stimuli_var0.XY,trials.stimuli_inX_var0.XY,trials.stimuli_inY_var0.XY] = intersect(trials.stimuli_var0.X,trials.stimuli_var0.Y);
    [trials.stimuli_var0.XB,trials.stimuli_inX_var0.XB,trials.stimuli_inB_var0.XB] = intersect(trials.stimuli_var0.X,trials.stimuli_var0.B);

    % odours - var1
    trials.stimuli_var1.A = intersect( intersect( find(task.odour1=="A") ,these_trials) ,find(task.var==1));
    trials.stimuli_var1.X = intersect( intersect( find(task.odour1=="X") ,these_trials) ,find(task.var==1));
    trials.stimuli_var1.B = intersect( intersect( find(task.odour2=="B") ,these_trials) ,find(task.var==1));
    trials.stimuli_var1.Y = intersect( intersect( find(task.odour2=="Y") ,these_trials) ,find(task.var==1));

    % odour combinations - var1
    [trials.stimuli_var1.AB,trials.stimuli_inA_var1.AB,trials.stimuli_inB_var1.AB] = intersect(trials.stimuli_var1.A,trials.stimuli_var1.B);
    [trials.stimuli_var1.AY,trials.stimuli_inA_var1.AY,trials.stimuli_inY_var1.AY] = intersect(trials.stimuli_var1.A,trials.stimuli_var1.Y);
    [trials.stimuli_var1.XY,trials.stimuli_inX_var1.XY,trials.stimuli_inY_var1.XY] = intersect(trials.stimuli_var1.X,trials.stimuli_var1.Y);
    [trials.stimuli_var1.XB,trials.stimuli_inX_var1.XB,trials.stimuli_inB_var1.XB] = intersect(trials.stimuli_var1.X,trials.stimuli_var1.B);

end


%% Trial types by outcome - by stim type

if do_by_stimType

    % outcomes - var0
    trials.outcome_var0.H = intersect( intersect( find(task.response=="H") ,these_trials) ,find(task.var==0));
    trials.outcome_var0.M = intersect( intersect( find(task.response=="M") ,these_trials) ,find(task.var==0));
    trials.outcome_var0.CR = intersect( intersect( find(task.response=="CR") ,these_trials) ,find(task.var==0));
    trials.outcome_var0.FA = intersect( intersect( find(task.response=="FA") ,these_trials) ,find(task.var==0));
    trials.outcome_var0.auto = intersect( intersect( find(task.autoreward) ,these_trials) ,find(task.var==0));
    trials.outcome_var0.correct = sort([trials.outcome_var0.H,trials.outcome_var0.CR]);
    trials.outcome_var0.incorrect = sort([trials.outcome_var0.M,trials.outcome_var0.FA]);
    trials.outcome_var0.rew = sort([trials.outcome_var0.H,trials.outcome_var0.auto]);
    trials.outcome_var0.norew = sort([trials.outcome_var0.M,trials.outcome_var0.CR,trials.outcome_var0.FA]);

    % outcomes - trial type combinations (A) - var0
    [trials.outcome_var0.A_H,trials.outcome_inA_var0.A_H,trials.outcome_inH_var0.A_H] = intersect(trials.stimuli_var0.A,trials.outcome_var0.H);
    [trials.outcome_var0.A_M,trials.outcome_inA_var0.A_M,trials.outcome_inM_var0.A_M] = intersect(trials.stimuli_var0.A,trials.outcome_var0.M);
    [trials.outcome_var0.A_CR,trials.outcome_inA_var0.A_CR,trials.outcome_inCR_var0.A_CR] = intersect(trials.stimuli_var0.A,trials.outcome_var0.CR);
    [trials.outcome_var0.A_FA,trials.outcome_inA_var0.A_FA,trials.outcome_inFA_var0.A_FA] = intersect(trials.stimuli_var0.A,trials.outcome_var0.FA);
    [trials.outcome_var0.A_correct,trials.outcome_inA_var0.A_correct,trials.outcome_incorrect_var0.A_correct] = intersect(trials.stimuli_var0.A,trials.outcome_var0.correct);
    [trials.outcome_var0.A_incorrect,trials.outcome_inA_var0.A_incorrect,trials.outcome_inincorrect_var0.A_incorrect] = intersect(trials.stimuli_var0.A,trials.outcome_var0.incorrect);
    [trials.outcome_var0.A_rew,trials.outcome_inA_var0.A_rew,trials.outcome_inrew_var0.A_rew] = intersect(trials.stimuli_var0.A,trials.outcome_var0.rew);
    [trials.outcome_var0.A_norew,trials.outcome_inA_var0.A_norew,trials.outcome_inrew_var0.A_norew] = intersect(trials.stimuli_var0.A,trials.outcome_var0.norew);

    % outcomes - trial type combinations (X) - var0
    [trials.outcome_var0.X_H,trials.outcome_inX_var0.X_H,trials.outcome_inH_var0.X_H] = intersect(trials.stimuli_var0.X,trials.outcome_var0.H);
    [trials.outcome_var0.X_M,trials.outcome_inX_var0.X_M,trials.outcome_inM_var0.X_M] = intersect(trials.stimuli_var0.X,trials.outcome_var0.M);
    [trials.outcome_var0.X_CR,trials.outcome_inX_var0.X_CR,trials.outcome_inCR_var0.X_CR] = intersect(trials.stimuli_var0.X,trials.outcome_var0.CR);
    [trials.outcome_var0.X_FA,trials.outcome_inX_var0.X_FA,trials.outcome_inFA_var0.X_FA] = intersect(trials.stimuli_var0.X,trials.outcome_var0.FA);
    [trials.outcome_var0.X_correct,trials.outcome_inX_var0.X_correct,trials.outcome_incorrect_var0.X_correct] = intersect(trials.stimuli_var0.X,trials.outcome_var0.correct);
    [trials.outcome_var0.X_incorrect,trials.outcome_inX_var0.X_incorrect,trials.outcome_inincorrect_var0.X_incorrect] = intersect(trials.stimuli_var0.X,trials.outcome_var0.incorrect);
    [trials.outcome_var0.X_rew,trials.outcome_inX_var0.X_rew,trials.outcome_inrew_var0.X_rew] = intersect(trials.stimuli_var0.X,trials.outcome_var0.rew);
    [trials.outcome_var0.X_norew,trials.outcome_inX_var0.X_norew,trials.outcome_inrew_var0.X_norew] = intersect(trials.stimuli_var0.X,trials.outcome_var0.norew);

    % outcomes - trial type combinations (B) - var0
    [trials.outcome_var0.B_H,trials.outcome_inB_var0.B_H,trials.outcome_inH_var0.B_H] = intersect(trials.stimuli_var0.B,trials.outcome_var0.H);
    [trials.outcome_var0.B_M,trials.outcome_inB_var0.B_M,trials.outcome_inM_var0.B_M] = intersect(trials.stimuli_var0.B,trials.outcome_var0.M);
    [trials.outcome_var0.B_CR,trials.outcome_inB_var0.B_CR,trials.outcome_inCR_var0.B_CR] = intersect(trials.stimuli_var0.B,trials.outcome_var0.CR);
    [trials.outcome_var0.B_FA,trials.outcome_inB_var0.B_FA,trials.outcome_inFA_var0.B_FA] = intersect(trials.stimuli_var0.B,trials.outcome_var0.FA);
    [trials.outcome_var0.B_correct,trials.outcome_inB_var0.B_correct,trials.outcome_incorrect_var0.B_correct] = intersect(trials.stimuli_var0.B,trials.outcome_var0.correct);
    [trials.outcome_var0.B_incorrect,trials.outcome_inB_var0.B_incorrect,trials.outcome_inincorrect_var0.B_incorrect] = intersect(trials.stimuli_var0.B,trials.outcome_var0.incorrect);
    [trials.outcome_var0.B_rew,trials.outcome_inB_var0.B_rew,trials.outcome_inrew_var0.B_rew] = intersect(trials.stimuli_var0.B,trials.outcome_var0.rew);
    [trials.outcome_var0.B_norew,trials.outcome_inB_var0.B_norew,trials.outcome_inrew_var0.B_norew] = intersect(trials.stimuli_var0.B,trials.outcome_var0.norew);

    % outcomes - trial type combinations (Y) - var0
    [trials.outcome_var0.Y_H,trials.outcome_inY_var0.Y_H,trials.outcome_inH_var0.Y_H] = intersect(trials.stimuli_var0.Y,trials.outcome_var0.H);
    [trials.outcome_var0.Y_M,trials.outcome_inY_var0.Y_M,trials.outcome_inM_var0.Y_M] = intersect(trials.stimuli_var0.Y,trials.outcome_var0.M);
    [trials.outcome_var0.Y_CR,trials.outcome_inY_var0.Y_CR,trials.outcome_inCR_var0.Y_CR] = intersect(trials.stimuli_var0.Y,trials.outcome_var0.CR);
    [trials.outcome_var0.Y_FA,trials.outcome_inY_var0.Y_FA,trials.outcome_inFA_var0.Y_FA] = intersect(trials.stimuli_var0.Y,trials.outcome_var0.FA);
    [trials.outcome_var0.Y_correct,trials.outcome_inY_var0.Y_correct,trials.outcome_incorrect_var0.Y_correct] = intersect(trials.stimuli_var0.Y,trials.outcome_var0.correct);
    [trials.outcome_var0.Y_incorrect,trials.outcome_inY_var0.Y_incorrect,trials.outcome_inincorrect_var0.Y_incorrect] = intersect(trials.stimuli_var0.Y,trials.outcome_var0.incorrect);
    [trials.outcome_var0.Y_rew,trials.outcome_inY_var0.Y_rew,trials.outcome_inrew_var0.Y_rew] = intersect(trials.stimuli_var0.Y,trials.outcome_var0.rew);
    [trials.outcome_var0.Y_norew,trials.outcome_inY_var0.Y_norew,trials.outcome_inrew_var0.Y_norew] = intersect(trials.stimuli_var0.Y,trials.outcome_var0.norew);

    % outcomes - var1
    trials.outcome_var1.H = intersect( intersect( find(task.response=="H") ,these_trials) ,find(task.var==1));
    trials.outcome_var1.M = intersect( intersect( find(task.response=="M") ,these_trials) ,find(task.var==1));
    trials.outcome_var1.CR = intersect( intersect( find(task.response=="CR") ,these_trials) ,find(task.var==1));
    trials.outcome_var1.FA = intersect( intersect( find(task.response=="FA") ,these_trials) ,find(task.var==1));
    trials.outcome_var1.auto = intersect( intersect( find(task.autoreward) ,these_trials) ,find(task.var==1));
    trials.outcome_var1.correct = sort([trials.outcome_var1.H,trials.outcome_var1.CR]);
    trials.outcome_var1.incorrect = sort([trials.outcome_var1.M,trials.outcome_var1.FA]);
    trials.outcome_var1.rew = sort([trials.outcome_var1.H,trials.outcome_var1.auto]);
    trials.outcome_var1.norew = sort([trials.outcome_var1.M,trials.outcome_var1.CR,trials.outcome_var1.FA]);

    % outcomes - trial type combinations (A) - var1
    [trials.outcome_var1.A_H,trials.outcome_inA_var1.A_H,trials.outcome_inH_var1.A_H] = intersect(trials.stimuli_var1.A,trials.outcome_var1.H);
    [trials.outcome_var1.A_M,trials.outcome_inA_var1.A_M,trials.outcome_inM_var1.A_M] = intersect(trials.stimuli_var1.A,trials.outcome_var1.M);
    [trials.outcome_var1.A_CR,trials.outcome_inA_var1.A_CR,trials.outcome_inCR_var1.A_CR] = intersect(trials.stimuli_var1.A,trials.outcome_var1.CR);
    [trials.outcome_var1.A_FA,trials.outcome_inA_var1.A_FA,trials.outcome_inFA_var1.A_FA] = intersect(trials.stimuli_var1.A,trials.outcome_var1.FA);
    [trials.outcome_var1.A_correct,trials.outcome_inA_var1.A_correct,trials.outcome_incorrect_var1.A_correct] = intersect(trials.stimuli_var1.A,trials.outcome_var1.correct);
    [trials.outcome_var1.A_incorrect,trials.outcome_inA_var1.A_incorrect,trials.outcome_inincorrect_var1.A_incorrect] = intersect(trials.stimuli_var1.A,trials.outcome_var1.incorrect);
    [trials.outcome_var1.A_rew,trials.outcome_inA_var1.A_rew,trials.outcome_inrew_var1.A_rew] = intersect(trials.stimuli_var1.A,trials.outcome_var1.rew);
    [trials.outcome_var1.A_norew,trials.outcome_inA_var1.A_norew,trials.outcome_inrew_var1.A_norew] = intersect(trials.stimuli_var1.A,trials.outcome_var1.norew);

    % outcomes - trial type combinations (X) - var1
    [trials.outcome_var1.X_H,trials.outcome_inX_var1.X_H,trials.outcome_inH_var1.X_H] = intersect(trials.stimuli_var1.X,trials.outcome_var1.H);
    [trials.outcome_var1.X_M,trials.outcome_inX_var1.X_M,trials.outcome_inM_var1.X_M] = intersect(trials.stimuli_var1.X,trials.outcome_var1.M);
    [trials.outcome_var1.X_CR,trials.outcome_inX_var1.X_CR,trials.outcome_inCR_var1.X_CR] = intersect(trials.stimuli_var1.X,trials.outcome_var1.CR);
    [trials.outcome_var1.X_FA,trials.outcome_inX_var1.X_FA,trials.outcome_inFA_var1.X_FA] = intersect(trials.stimuli_var1.X,trials.outcome_var1.FA);
    [trials.outcome_var1.X_correct,trials.outcome_inX_var1.X_correct,trials.outcome_incorrect_var1.X_correct] = intersect(trials.stimuli_var1.X,trials.outcome_var1.correct);
    [trials.outcome_var1.X_incorrect,trials.outcome_inX_var1.X_incorrect,trials.outcome_inincorrect_var1.X_incorrect] = intersect(trials.stimuli_var1.X,trials.outcome_var1.incorrect);
    [trials.outcome_var1.X_rew,trials.outcome_inX_var1.X_rew,trials.outcome_inrew_var1.X_rew] = intersect(trials.stimuli_var1.X,trials.outcome_var1.rew);
    [trials.outcome_var1.X_norew,trials.outcome_inX_var1.X_norew,trials.outcome_inrew_var1.X_norew] = intersect(trials.stimuli_var1.X,trials.outcome_var1.norew);

    % outcomes - trial type combinations (B) - var1
    [trials.outcome_var1.B_H,trials.outcome_inB_var1.B_H,trials.outcome_inH_var1.B_H] = intersect(trials.stimuli_var1.B,trials.outcome_var1.H);
    [trials.outcome_var1.B_M,trials.outcome_inB_var1.B_M,trials.outcome_inM_var1.B_M] = intersect(trials.stimuli_var1.B,trials.outcome_var1.M);
    [trials.outcome_var1.B_CR,trials.outcome_inB_var1.B_CR,trials.outcome_inCR_var1.B_CR] = intersect(trials.stimuli_var1.B,trials.outcome_var1.CR);
    [trials.outcome_var1.B_FA,trials.outcome_inB_var1.B_FA,trials.outcome_inFA_var1.B_FA] = intersect(trials.stimuli_var1.B,trials.outcome_var1.FA);
    [trials.outcome_var1.B_correct,trials.outcome_inB_var1.B_correct,trials.outcome_incorrect_var1.B_correct] = intersect(trials.stimuli_var1.B,trials.outcome_var1.correct);
    [trials.outcome_var1.B_incorrect,trials.outcome_inB_var1.B_incorrect,trials.outcome_inincorrect_var1.B_incorrect] = intersect(trials.stimuli_var1.B,trials.outcome_var1.incorrect);
    [trials.outcome_var1.B_rew,trials.outcome_inB_var1.B_rew,trials.outcome_inrew_var1.B_rew] = intersect(trials.stimuli_var1.B,trials.outcome_var1.rew);
    [trials.outcome_var1.B_norew,trials.outcome_inB_var1.B_norew,trials.outcome_inrew_var1.B_norew] = intersect(trials.stimuli_var1.B,trials.outcome_var1.norew);

    % outcomes - trial type combinations (Y) - var1
    [trials.outcome_var1.Y_H,trials.outcome_inY_var1.Y_H,trials.outcome_inH_var1.Y_H] = intersect(trials.stimuli_var1.Y,trials.outcome_var1.H);
    [trials.outcome_var1.Y_M,trials.outcome_inY_var1.Y_M,trials.outcome_inM_var1.Y_M] = intersect(trials.stimuli_var1.Y,trials.outcome_var1.M);
    [trials.outcome_var1.Y_CR,trials.outcome_inY_var1.Y_CR,trials.outcome_inCR_var1.Y_CR] = intersect(trials.stimuli_var1.Y,trials.outcome_var1.CR);
    [trials.outcome_var1.Y_FA,trials.outcome_inY_var1.Y_FA,trials.outcome_inFA_var1.Y_FA] = intersect(trials.stimuli_var1.Y,trials.outcome_var1.FA);
    [trials.outcome_var1.Y_correct,trials.outcome_inY_var1.Y_correct,trials.outcome_incorrect_var1.Y_correct] = intersect(trials.stimuli_var1.Y,trials.outcome_var1.correct);
    [trials.outcome_var1.Y_incorrect,trials.outcome_inY_var1.Y_incorrect,trials.outcome_inincorrect_var1.Y_incorrect] = intersect(trials.stimuli_var1.Y,trials.outcome_var1.incorrect);
    [trials.outcome_var1.Y_rew,trials.outcome_inY_var1.Y_rew,trials.outcome_inrew_var1.Y_rew] = intersect(trials.stimuli_var1.Y,trials.outcome_var1.rew);
    [trials.outcome_var1.Y_norew,trials.outcome_inY_var1.Y_norew,trials.outcome_inrew_var1.Y_norew] = intersect(trials.stimuli_var1.Y,trials.outcome_var1.norew);

end


%% Return

trials = orderfields(trials);
end

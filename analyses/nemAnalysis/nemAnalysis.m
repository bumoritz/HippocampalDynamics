function nemAnalysis(info,iscell,ops,p,path,task,paq_beh,act)
% act = spks_beh;

%% Preparations

disp('--- Preparations.')

% start parallel pool (if not started yet)
try
    parpool; % can be closed with: delete(gcp('nocreate'));
catch
end
rng(p.general.rgnSeed);

% pre-process activity measure and paq events
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.nem,p,paq_beh,iscell);
[events_binned] = binPaqEvents(paq_beh,task,p,prop);
% [events_binned] = binPaqEvents(paq_beh,task,p,prop,cam_beh);

%% Core (incl saving)

disp('--- Core. Fitting models and predicting activity.')

these_cells = 1:length(prop.iscell);

% run nemCore function
if ops.nem.do_allTrials
    disp('--- Fitting models and predicting activity for all trials...')
    [nem_all,nemd_all] = nemCore(nf_binned,events_binned,task,info,p,prop,these_cells); 
    %nemdp_all = nemMakePartialPredictions(nem_all,nemd_all,prop);
    
    % save
    save([path.filepart_outX,'nem_all.mat'],'nem_all','-v7.3');
    disp(['--- Saved nem_all file as ',[path.filepart_outX,'nem_all.mat'],'.'])
    save([path.filepart_outX,'nemd_all.mat'],'nemd_all','-v7.3');
    disp(['--- Saved nemd_all file as ',[path.filepart_outX,'nemd_all.mat'],'.'])
        
    % create trials and traces structs
    trials_all = createTrialsStruct(task,1:prop.numTrials);
    [traces_all,~,~] = createTracesStructs(nft_binned,trials_all);

    % model setup figure
    F = nemModelSetupFigure(nem_all,events_binned,50:55,info,p);
    saveas(F,[path.filepart_outX,'plots/',info.animal,'_',info.date,'_nem_modelSetup.png']);
    
    % population-wide encoding figure
    F = nemPopulationWideEncodingFigure(nem_all,info,prop);
    saveas(F,[path.filepart_outX,'plots/',info.animal,'_',info.date,'_nem_populationWideEncoding.png']);
    
    % example single-cell figure
    F = nemSingleCellFigure(nem_all,nemd_all,[],[],[],trials_all,traces_all,nanmin(find(prop.iscell)),info,p,prop);
    saveas(F,[path.filepart_outX,'plots/',info.animal,'_',info.date,'_nem_exampleSingleCellFigure.png']);

    clear nem_all; 
    clear nemd_all;
    clear traces_all;
end

if ops.nem.do_100t || ops.nem.do_100t_onlyFirst
    nem_100t = {};
    if ~ops.nem.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        disp(['--- Fitting models and predicting activity for 100t block number ,',num2str(i),', ...'])
        these_trials = (i-1)*100+1:i*100;
        
        temp1 = find(events_binned.sync);       
        these_bins = zeros(1,prop.numFramesTotal_binned);
        if min(these_trials)==1
            these_bins(1:temp1(max(these_trials)+1)-length(p.general.bins_pre)) = 1;
        elseif max(these_trials)==prop.numTrials
            these_bins(temp1(min(these_trials))-length(p.general.bins_pre)+1:end) = 1;
        else
            these_bins(temp1(min(these_trials))-length(p.general.bins_pre)+1:temp1(max(these_trials)+1)-length(p.general.bins_pre)) = 1;
        end
        these_bins = find(these_bins);
        
        this_nf_binned = nf_binned(:,these_bins);
        
        these_events_binned = {};
        temp = fields(events_binned);
        for j=1:length(temp)
            these_events_binned.(temp{j}) = events_binned.(temp{j})(these_bins);
        end
        
        this_task = {};
        temp = fields(task);
        for j=1:length(temp)
            this_task.(temp{j}) = task.(temp{j})(these_trials);
        end
        
        [nem_100t{i},nemd_100t{i}] = nemCore(this_nf_binned,these_events_binned,this_task,info,p,prop,these_cells);
        %nemdp_100t{i} = nemMakePartialPredictions(nem_100t{i},nemd_100t{i},prop);
    end
    
    % save
    save([path.filepart_outX,'nem_100t.mat'],'nem_100t','-v7.3');
    disp(['--- Saved nem_100t file as ',[path.filepart_outX,'nem_100t.mat'],'.'])
    save([path.filepart_outX,'nemd_100t.mat'],'nemd_100t','-v7.3');
    disp(['--- Saved nemd_100t file as ',[path.filepart_outX,'nemd_100t.mat'],'.'])
    
    % block-wise encoding figure
    F = nemBlockWiseEncodingFigure(nem_100t,info,p,prop);
    saveas(F,[path.filepart_outX,'plots/',info.animal,'_',info.date,'_nem_blockWiseEncoding.png']);
    
    clear nem_100t; 
    clear nemd_100t;
end


%% Return

if ops.close_figures
    close all;
end
end

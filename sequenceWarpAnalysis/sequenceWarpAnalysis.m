function [warp] = sequenceWarpAnalysis(info,iscell,ops,p,path,tng,input_type)
%tng = tng_100t_stimVersion; input_type = '100t_stimVersion';
%tng = tng_100t; input_type = '100t';
%tng = tng_all; input_type = 'all';
%tng = tng_all_stimVersion; input_type = 'all_stimVersion';

disp(['--- Running sequence warp analysis on ',input_type])


%% Core - warp_all and warp_all_stimVersion

if (strcmp(input_type,'all') || strcmp(input_type,'all_stimVersion'))
    this_tng = tng;

    % get data
    warp.p = p;
    warp.input.numCells = length(find(tng.prop.iscell==1));
    if strcmp(input_type,'all_stimVersion')
        warp.input.idcs_Aonly = find(tng.passed_catch.AW.Acatchonly==1);
        warp.input.idcs_Xonly = find(tng.passed_catch.AW.Xcatchonly==1);
    else
        warp.input.idcs_Aonly = find(tng.passed.AW.Aonly==1);
        warp.input.idcs_Xonly = find(tng.passed.AW.Xonly==1);
    end

    % get peak locations
    if strcmp(input_type,'all_stimVersion')
        warp.input.peakTimes_Aonly = tng.firingField.Acatch_AW.peakLocation_s(warp.input.idcs_Aonly);
        warp.input.peakTimes_Xonly = tng.firingField.Xcatch_AW.peakLocation_s(warp.input.idcs_Xonly);
    else
        warp.input.peakTimes_Aonly = tng.firingField.A_AW.peakLocation_s(warp.input.idcs_Aonly);
        warp.input.peakTimes_Xonly = tng.firingField.X_AW.peakLocation_s(warp.input.idcs_Xonly);
    end
    warp.input.peakTimes = [warp.input.peakTimes_Aonly; warp.input.peakTimes_Xonly];

    [warp.peak.hist.prob_Aonly,~] = histcounts(warp.input.peakTimes_Aonly,p.warp.histBins,'Normalization','probability');
    [warp.peak.hist.prob_Xonly,~] = histcounts(warp.input.peakTimes_Xonly,p.warp.histBins,'Normalization','probability');
    [warp.peak.hist.prob,warp.peak.hist.edges] = histcounts(warp.input.peakTimes,p.warp.histBins,'Normalization','probability');

    [warp.peak.poly1_Aonly,warp.peak.gof.poly1_Aonly] = fit(warp.peak.hist.edges(1:end-1)'+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob_Aonly','poly1');
    [warp.peak.exp1_Aonly,warp.peak.gof.exp1_Aonly] = fit(warp.peak.hist.edges(1:end-1)'+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob_Aonly','exp1');
    [warp.peak.power2_Aonly,warp.peak.gof.power2_Aonly] = fit(warp.peak.hist.edges(1:end-1)'+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob_Aonly','power2');
    [warp.peak.poly1_Xonly,warp.peak.gof.poly1_Xonly] = fit(warp.peak.hist.edges(1:end-1)'+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob_Xonly','poly1');
    [warp.peak.exp1_Xonly,warp.peak.gof.exp1_Xonly] = fit(warp.peak.hist.edges(1:end-1)'+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob_Xonly','exp1');
    [warp.peak.power2_Xonly,warp.peak.gof.power2_Xonly] = fit(warp.peak.hist.edges(1:end-1)'+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob_Xonly','power2');
    [warp.peak.poly1,warp.peak.gof.poly1] = fit(warp.peak.hist.edges(1:end-1)'+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob','poly1');
    [warp.peak.exp1,warp.peak.gof.exp1] = fit(warp.peak.hist.edges(1:end-1)'+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob','exp1');
    [warp.peak.power2,warp.peak.gof.power2] = fit(warp.peak.hist.edges(1:end-1)'+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob','power2');
    
    % calculate linearity index
%     warp.peak.linIdx1_Aonly = warp.peak.gof.poly1.rsquare_Aonly / warp.peak.gof.exp1.rsquare_Aonly;
%     warp.peak.linIdx2_Aonly = (warp.peak.gof.poly1.rsquare_Aonly - warp.peak.gof.exp1.rsquare_Aonly) / (warp.peak.gof.poly1.rsquare_Aonly + warp.peak.gof.exp1.rsquare_Aonly);
%     warp.peak.linIdx1_Xonly = warp.peak.gof.poly1.rsquare_Xonly / warp.peak.gof.exp1.rsquare_Xonly;
%     warp.peak.linIdx2_Xonly = (warp.peak.gof.poly1.rsquare_Xonly - warp.peak.gof.exp1.rsquare_Xonly) / (warp.peak.gof.poly1.rsquare_Xonly + warp.peak.gof.exp1.rsquare_Xonly);
    warp.peak.linIdx1 = warp.peak.gof.poly1.rsquare / warp.peak.gof.exp1.rsquare;
    warp.peak.linIdx2 = (warp.peak.gof.poly1.rsquare - warp.peak.gof.exp1.rsquare) / warp.peak.gof.exp1.rsquare;
    warp.peak.linIdx3 = (warp.peak.gof.poly1.rsquare - warp.peak.gof.exp1.rsquare) / (warp.peak.gof.poly1.rsquare + warp.peak.gof.exp1.rsquare);
    warp.peak.linIdx4 = warp.peak.gof.poly1.rmse / warp.peak.gof.exp1.rmse;
    warp.peak.linIdx5 = (warp.peak.gof.poly1.rmse - warp.peak.gof.exp1.rmse) / warp.peak.gof.exp1.rmse;
    warp.peak.linIdx6 = (warp.peak.gof.poly1.rmse - warp.peak.gof.exp1.rmse) / (warp.peak.gof.poly1.rmse + warp.peak.gof.exp1.rmse);
    warp.peak.linIdx7 = warp.peak.gof.poly1.sse / warp.peak.gof.exp1.sse;
    warp.peak.linIdx8 = (warp.peak.gof.poly1.sse - warp.peak.gof.exp1.sse) / warp.peak.gof.exp1.sse;
    warp.peak.linIdx9 = (warp.peak.gof.poly1.sse - warp.peak.gof.exp1.sse) / (warp.peak.gof.poly1.sse + warp.peak.gof.exp1.sse);
    this_data_lin = warp.peak.gof.poly1.rmse.^2; this_data_nonlin = warp.peak.gof.exp1.rmse.^2;
    warp.peak.linIdx10 = this_data_lin / this_data_nonlin;
    warp.peak.linIdx11 = (this_data_lin - this_data_nonlin) / this_data_nonlin;
    warp.peak.linIdx12 = (this_data_lin - this_data_nonlin) / (this_data_lin + this_data_nonlin);
    warp.peak.nonLinIdx1 = warp.peak.gof.exp1.rsquare / warp.peak.gof.poly1.rsquare;
    warp.peak.nonLinIdx2 = (warp.peak.gof.exp1.rsquare - warp.peak.gof.poly1.rsquare) / warp.peak.gof.poly1.rsquare;
    warp.peak.nonLinIdx3 = (warp.peak.gof.exp1.rsquare - warp.peak.gof.poly1.rsquare) / (warp.peak.gof.poly1.rsquare + warp.peak.gof.exp1.rsquare);
    warp.peak.nonLinIdx4 = warp.peak.gof.exp1.rmse / warp.peak.gof.poly1.rmse;
    warp.peak.nonLinIdx5 = (warp.peak.gof.exp1.rmse - warp.peak.gof.poly1.rmse) / warp.peak.gof.poly1.rmse;
    warp.peak.nonLinIdx6 = (warp.peak.gof.exp1.rmse - warp.peak.gof.poly1.rmse) / (warp.peak.gof.poly1.rmse + warp.peak.gof.exp1.rmse);
    warp.peak.nonLinIdx7 = warp.peak.gof.exp1.sse / warp.peak.gof.poly1.sse;
    warp.peak.nonLinIdx8 = (warp.peak.gof.exp1.sse - warp.peak.gof.poly1.sse) / warp.peak.gof.poly1.sse;
    warp.peak.nonLinIdx9 = (warp.peak.gof.exp1.sse - warp.peak.gof.poly1.sse) / (warp.peak.gof.poly1.sse + warp.peak.gof.exp1.sse);
    this_data_lin = warp.peak.gof.poly1.rmse.^2; this_data_nonlin = warp.peak.gof.exp1.rmse.^2;
    warp.peak.nonLinIdx10 = this_data_nonlin / this_data_lin;
    warp.peak.nonLinIdx11 = (this_data_nonlin - this_data_lin) / this_data_lin;
    warp.peak.nonLinIdx12 = (this_data_nonlin - this_data_lin) / (this_data_lin + this_data_nonlin);
        
    % get sequence cell count
    warp.input.numSequenceCells = length(warp.input.peakTimes);
    
    % get com locations
    if strcmp(input_type,'all_stimVersion')
        warp.input.comTimes_Aonly = tng.firingField.Acatch_AW.comLocation_s(warp.input.idcs_Aonly);
        warp.input.comTimes_Xonly = tng.firingField.Xcatch_AW.comLocation_s(warp.input.idcs_Xonly);
    else
        warp.input.comTimes_Aonly = tng.firingField.A_AW.comLocation_s(warp.input.idcs_Aonly);
        warp.input.comTimes_Xonly = tng.firingField.X_AW.comLocation_s(warp.input.idcs_Xonly);
    end
    warp.input.comTimes = [warp.input.comTimes_Aonly; warp.input.comTimes_Xonly];

    % get com histogram
    [warp.com.hist.prob_Aonly,~] = histcounts(warp.input.comTimes_Aonly,p.warp.histBins,'Normalization','probability');
    [warp.com.hist.prob_Xonly,~] = histcounts(warp.input.comTimes_Xonly,p.warp.histBins,'Normalization','probability');
    [warp.com.hist.prob,warp.com.hist.edges] = histcounts(warp.input.comTimes,p.warp.histBins,'Normalization','probability');
    
    % fit com histograms
    [warp.com.poly1_Aonly,warp.com.gof.poly1_Aonly] = fit(warp.com.hist.edges(1:end-1)'+(warp.com.hist.edges(2)-warp.com.hist.edges(1))/2, warp.com.hist.prob_Aonly','poly1');
    [warp.com.exp1_Aonly,warp.com.gof.exp1_Aonly] = fit(warp.com.hist.edges(1:end-1)'+(warp.com.hist.edges(2)-warp.com.hist.edges(1))/2, warp.com.hist.prob_Aonly','exp1');
    [warp.com.power2_Aonly,warp.com.gof.power2_Aonly] = fit(warp.com.hist.edges(1:end-1)'+(warp.com.hist.edges(2)-warp.com.hist.edges(1))/2, warp.com.hist.prob_Aonly','power2');
    [warp.com.poly1_Xonly,warp.com.gof.poly1_Xonly] = fit(warp.com.hist.edges(1:end-1)'+(warp.com.hist.edges(2)-warp.com.hist.edges(1))/2, warp.com.hist.prob_Xonly','poly1');
    [warp.com.exp1_Xonly,warp.com.gof.exp1_Xonly] = fit(warp.com.hist.edges(1:end-1)'+(warp.com.hist.edges(2)-warp.com.hist.edges(1))/2, warp.com.hist.prob_Xonly','exp1');
    [warp.com.power2_Xonly,warp.com.gof.power2_Xonly] = fit(warp.com.hist.edges(1:end-1)'+(warp.com.hist.edges(2)-warp.com.hist.edges(1))/2, warp.com.hist.prob_Xonly','power2');
    [warp.com.poly1,warp.com.gof.poly1] = fit(warp.com.hist.edges(1:end-1)'+(warp.com.hist.edges(2)-warp.com.hist.edges(1))/2, warp.com.hist.prob','poly1');
    [warp.com.exp1,warp.com.gof.exp1] = fit(warp.com.hist.edges(1:end-1)'+(warp.com.hist.edges(2)-warp.com.hist.edges(1))/2, warp.com.hist.prob','exp1');
    [warp.com.power2,warp.com.gof.power2] = fit(warp.com.hist.edges(1:end-1)'+(warp.com.hist.edges(2)-warp.com.hist.edges(1))/2, warp.com.hist.prob','power2');
    
    % calculate linearity index
    warp.com.linIdx1 = warp.com.gof.poly1.rsquare / warp.com.gof.exp1.rsquare;
    warp.com.linIdx2 = (warp.com.gof.poly1.rsquare - warp.com.gof.exp1.rsquare) / warp.com.gof.exp1.rsquare;
    warp.com.linIdx3 = (warp.com.gof.poly1.rsquare - warp.com.gof.exp1.rsquare) / (warp.com.gof.poly1.rsquare + warp.com.gof.exp1.rsquare);
    warp.com.linIdx4 = warp.com.gof.poly1.rmse / warp.com.gof.exp1.rmse;
    warp.com.linIdx5 = (warp.com.gof.poly1.rmse - warp.com.gof.exp1.rmse) / warp.com.gof.exp1.rmse;
    warp.com.linIdx6 = (warp.com.gof.poly1.rmse - warp.com.gof.exp1.rmse) / (warp.com.gof.poly1.rmse + warp.com.gof.exp1.rmse);
    warp.com.linIdx7 = warp.com.gof.poly1.sse / warp.com.gof.exp1.sse;
    warp.com.linIdx8 = (warp.com.gof.poly1.sse - warp.com.gof.exp1.sse) / warp.com.gof.exp1.sse;
    warp.com.linIdx9 = (warp.com.gof.poly1.sse - warp.com.gof.exp1.sse) / (warp.com.gof.poly1.sse + warp.com.gof.exp1.sse);
    this_data_lin = warp.com.gof.poly1.rmse.^2; this_data_nonlin = warp.com.gof.exp1.rmse.^2;
    warp.com.linIdx10 = this_data_lin / this_data_nonlin;
    warp.com.linIdx11 = (this_data_lin - this_data_nonlin) / this_data_nonlin;
    warp.com.linIdx12 = (this_data_lin - this_data_nonlin) / (this_data_lin + this_data_nonlin);
    warp.com.nonLinIdx1 = warp.com.gof.exp1.rsquare / warp.com.gof.poly1.rsquare;
    warp.com.nonLinIdx2 = (warp.com.gof.exp1.rsquare - warp.com.gof.poly1.rsquare) / warp.com.gof.poly1.rsquare;
    warp.com.nonLinIdx3 = (warp.com.gof.exp1.rsquare - warp.com.gof.poly1.rsquare) / (warp.com.gof.poly1.rsquare + warp.com.gof.exp1.rsquare);
    warp.com.nonLinIdx4 = warp.com.gof.exp1.rmse / warp.com.gof.poly1.rmse;
    warp.com.nonLinIdx5 = (warp.com.gof.exp1.rmse - warp.com.gof.poly1.rmse) / warp.com.gof.poly1.rmse;
    warp.com.nonLinIdx6 = (warp.com.gof.exp1.rmse - warp.com.gof.poly1.rmse) / (warp.com.gof.poly1.rmse + warp.com.gof.exp1.rmse);
    warp.com.nonLinIdx7 = warp.com.gof.exp1.sse / warp.com.gof.poly1.sse;
    warp.com.nonLinIdx8 = (warp.com.gof.exp1.sse - warp.com.gof.poly1.sse) / warp.com.gof.poly1.sse;
    warp.com.nonLinIdx9 = (warp.com.gof.exp1.sse - warp.com.gof.poly1.sse) / (warp.com.gof.poly1.sse + warp.com.gof.exp1.sse);
    this_data_lin = warp.com.gof.poly1.rmse.^2; this_data_nonlin = warp.com.gof.exp1.rmse.^2;
    warp.com.nonLinIdx10 = this_data_nonlin / this_data_lin;
    warp.com.nonLinIdx11 = (this_data_nonlin - this_data_lin) / this_data_lin;
    warp.com.nonLinIdx12 = (this_data_nonlin - this_data_lin) / (this_data_lin + this_data_nonlin);
    
%     % get better absolute histogram
%     warp.input.peakTimes_unique = unique(warp.input.peakTimes);
%     warp.peak.hist.count = nan(1,length(warp.input.peakTimes_unique));
%     for i=1:length(warp.input.peakTimes_unique)
%         warp.peak.hist.count(i) = nansum(warp.input.peakTimes==warp.input.peakTimes_unique(i));
%     end
    
end


%% Core - warp_100t and warp_100t_stimVersion

if (strcmp(input_type,'100t') || strcmp(input_type,'100t_stimVersion'))
    warp = {};
    for k=1:length(tng)
        this_tng = tng{k};
        
        % get data
        warp{k}.p = p;
        warp{k}.input.numCells = length(find(tng{k}.prop.iscell==1));
        if strcmp(input_type,'100t_stimVersion')
            warp{k}.input.idcs_Aonly = find(tng{k}.passed_catch.AW.Acatchonly==1);
            warp{k}.input.idcs_Xonly = find(tng{k}.passed_catch.AW.Xcatchonly==1);
        else
            warp{k}.input.idcs_Aonly = find(tng{k}.passed.AW.Aonly==1);
            warp{k}.input.idcs_Xonly = find(tng{k}.passed.AW.Xonly==1);
        end

        % get peak locations
        if strcmp(input_type,'100t_stimVersion')
            warp{k}.input.peakTimes_Aonly = tng{k}.firingField.Acatch_AW.peakLocation_s(warp{k}.input.idcs_Aonly);
            warp{k}.input.peakTimes_Xonly = tng{k}.firingField.Xcatch_AW.peakLocation_s(warp{k}.input.idcs_Xonly);
        else
            warp{k}.input.peakTimes_Aonly = tng{k}.firingField.A_AW.peakLocation_s(warp{k}.input.idcs_Aonly);
            warp{k}.input.peakTimes_Xonly = tng{k}.firingField.X_AW.peakLocation_s(warp{k}.input.idcs_Xonly);
        end
        warp{k}.input.peakTimes = [warp{k}.input.peakTimes_Aonly; warp{k}.input.peakTimes_Xonly];
        
        [warp{k}.peak.hist.prob_Aonly,~] = histcounts(warp{k}.input.peakTimes_Aonly,p.warp.histBins,'Normalization','probability');
        [warp{k}.peak.hist.prob_Xonly,~] = histcounts(warp{k}.input.peakTimes_Xonly,p.warp.histBins,'Normalization','probability');
        [warp{k}.peak.hist.prob,warp{k}.peak.hist.edges] = histcounts(warp{k}.input.peakTimes,p.warp.histBins,'Normalization','probability');

        [warp{k}.peak.poly1_Aonly,warp{k}.peak.gof.poly1_Aonly] = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob_Aonly','poly1');
        [warp{k}.peak.exp1_Aonly,warp{k}.peak.gof.exp1_Aonly] = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob_Aonly','exp1');
        [warp{k}.peak.power2_Aonly,warp{k}.peak.gof.power2_Aonly] = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob_Aonly','power2');
%         warp{k}.peak.exp1const_Aonly = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob_Aonly',fittype(@(a,b,c,x) a.*exp(b.*x) + c),'StartPoint',[0,-0.5,0]);
        [warp{k}.peak.poly1_Xonly,warp{k}.peak.gof.poly1_Xonly] = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob_Xonly','poly1');
        [warp{k}.peak.exp1_Xonly,warp{k}.peak.gof.exp1_Xonly] = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob_Xonly','exp1');
        [warp{k}.peak.power2_Xonly,warp{k}.peak.gof.power2_Xonly] = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob_Xonly','power2');
%         warp{k}.peak.exp1const_Xonly = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob_Xonly',fittype(@(a,b,c,x) a.*exp(b.*x) + c),'StartPoint',[0,-0.5,0]);
        [warp{k}.peak.poly1,warp{k}.peak.gof.poly1] = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob','poly1');  
        [warp{k}.peak.exp1,warp{k}.peak.gof.exp1] = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob','exp1');  
        [warp{k}.peak.power2,warp{k}.peak.gof.power2] = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob','power2');
%         warp{k}.peak.exp1const = fit(warp{k}.peak.hist.edges(1:end-1)'+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob',fittype(@(a,b,c,x) a.*exp(b.*x) + c),'StartPoint',[0,-0.5,0]);

        % calculate linearity index
        warp{k}.peak.linIdx1 = warp{k}.peak.gof.poly1.rsquare / warp{k}.peak.gof.exp1.rsquare;
        warp{k}.peak.linIdx2 = (warp{k}.peak.gof.poly1.rsquare - warp{k}.peak.gof.exp1.rsquare) / warp{k}.peak.gof.exp1.rsquare;
        warp{k}.peak.linIdx3 = (warp{k}.peak.gof.poly1.rsquare - warp{k}.peak.gof.exp1.rsquare) / (warp{k}.peak.gof.poly1.rsquare + warp{k}.peak.gof.exp1.rsquare);
        warp{k}.peak.linIdx4 = warp{k}.peak.gof.poly1.rmse / warp{k}.peak.gof.exp1.rmse;
        warp{k}.peak.linIdx5 = (warp{k}.peak.gof.poly1.rmse - warp{k}.peak.gof.exp1.rmse) / warp{k}.peak.gof.exp1.rmse;
        warp{k}.peak.linIdx6 = (warp{k}.peak.gof.poly1.rmse - warp{k}.peak.gof.exp1.rmse) / (warp{k}.peak.gof.poly1.rmse + warp{k}.peak.gof.exp1.rmse);
        warp{k}.peak.linIdx7 = warp{k}.peak.gof.poly1.sse / warp{k}.peak.gof.exp1.sse;
        warp{k}.peak.linIdx8 = (warp{k}.peak.gof.poly1.sse - warp{k}.peak.gof.exp1.sse) / warp{k}.peak.gof.exp1.sse;
        warp{k}.peak.linIdx9 = (warp{k}.peak.gof.poly1.sse - warp{k}.peak.gof.exp1.sse) / (warp{k}.peak.gof.poly1.sse + warp{k}.peak.gof.exp1.sse);
        this_data_lin = warp{k}.peak.gof.poly1.rmse.^2; this_data_nonlin = warp{k}.peak.gof.exp1.rmse.^2;
        warp{k}.peak.linIdx10 = this_data_lin / this_data_nonlin;
        warp{k}.peak.linIdx11 = (this_data_lin - this_data_nonlin) / this_data_nonlin;
        warp{k}.peak.linIdx12 = (this_data_lin - this_data_nonlin) / (this_data_lin + this_data_nonlin);
        warp{k}.peak.nonLinIdx1 = warp{k}.peak.gof.exp1.rsquare / warp{k}.peak.gof.poly1.rsquare;
        warp{k}.peak.nonLinIdx2 = (warp{k}.peak.gof.exp1.rsquare - warp{k}.peak.gof.poly1.rsquare) / warp{k}.peak.gof.poly1.rsquare;
        warp{k}.peak.nonLinIdx3 = (warp{k}.peak.gof.exp1.rsquare - warp{k}.peak.gof.poly1.rsquare) / (warp{k}.peak.gof.poly1.rsquare + warp{k}.peak.gof.exp1.rsquare);
        warp{k}.peak.nonLinIdx4 = warp{k}.peak.gof.exp1.rmse / warp{k}.peak.gof.poly1.rmse;
        warp{k}.peak.nonLinIdx5 = (warp{k}.peak.gof.exp1.rmse - warp{k}.peak.gof.poly1.rmse) / warp{k}.peak.gof.poly1.rmse;
        warp{k}.peak.nonLinIdx6 = (warp{k}.peak.gof.exp1.rmse - warp{k}.peak.gof.poly1.rmse) / (warp{k}.peak.gof.poly1.rmse + warp{k}.peak.gof.exp1.rmse);
        warp{k}.peak.nonLinIdx7 = warp{k}.peak.gof.exp1.sse / warp{k}.peak.gof.poly1.sse;
        warp{k}.peak.nonLinIdx8 = (warp{k}.peak.gof.exp1.sse - warp{k}.peak.gof.poly1.sse) / warp{k}.peak.gof.poly1.sse;
        warp{k}.peak.nonLinIdx9 = (warp{k}.peak.gof.exp1.sse - warp{k}.peak.gof.poly1.sse) / (warp{k}.peak.gof.poly1.sse + warp{k}.peak.gof.exp1.sse);
        this_data_lin = warp{k}.peak.gof.poly1.rmse.^2; this_data_nonlin = warp{k}.peak.gof.exp1.rmse.^2;
        warp{k}.peak.nonLinIdx10 = this_data_nonlin / this_data_lin;
        warp{k}.peak.nonLinIdx11 = (this_data_nonlin - this_data_lin) / this_data_lin;
        warp{k}.peak.nonLinIdx12 = (this_data_nonlin - this_data_lin) / (this_data_lin + this_data_nonlin);
        
        % get sequence cell count
        warp{k}.input.numSequenceCells = length(warp{k}.input.peakTimes);
        
        % get com locations
        if strcmp(input_type,'100t_stimVersion')
            warp{k}.input.comTimes_Aonly = tng{k}.firingField.Acatch_AW.comLocation_s(warp{k}.input.idcs_Aonly);
            warp{k}.input.comTimes_Xonly = tng{k}.firingField.Xcatch_AW.comLocation_s(warp{k}.input.idcs_Xonly);
        else
            warp{k}.input.comTimes_Aonly = tng{k}.firingField.A_AW.comLocation_s(warp{k}.input.idcs_Aonly);
            warp{k}.input.comTimes_Xonly = tng{k}.firingField.X_AW.comLocation_s(warp{k}.input.idcs_Xonly);
        end
        warp{k}.input.comTimes = [warp{k}.input.comTimes_Aonly; warp{k}.input.comTimes_Xonly];
                
        [warp{k}.com.hist.prob_Aonly,~] = histcounts(warp{k}.input.comTimes_Aonly,p.warp.histBins,'Normalization','probability');
        [warp{k}.com.hist.prob_Xonly,~] = histcounts(warp{k}.input.comTimes_Xonly,p.warp.histBins,'Normalization','probability');
        [warp{k}.com.hist.prob,warp{k}.com.hist.edges] = histcounts(warp{k}.input.comTimes,p.warp.histBins,'Normalization','probability');

        [warp{k}.com.poly1_Aonly,warp{k}.com.gof.poly1_Aonly] = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob_Aonly','poly1');
        [warp{k}.com.exp1_Aonly,warp{k}.com.gof.exp1_Aonly] = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob_Aonly','exp1');
        [warp{k}.com.power2_Aonly,warp{k}.com.gof.power2_Aonly] = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob_Aonly','power2');
%         warp{k}.com.exp1const_Aonly = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob_Aonly',fittype(@(a,b,c,x) a.*exp(b.*x) + c),'StartPoint',[0,-0.5,0]);
        [warp{k}.com.poly1_Xonly,warp{k}.com.gof.poly1_Xonly] = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob_Xonly','poly1');
        [warp{k}.com.exp1_Xonly,warp{k}.com.gof.exp1_Xonly] = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob_Xonly','exp1');
        [warp{k}.com.power2_Xonly,warp{k}.com.gof.power2_Xonly] = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob_Xonly','power2');
%         warp{k}.com.exp1const_Xonly = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob_Xonly',fittype(@(a,b,c,x) a.*exp(b.*x) + c),'StartPoint',[0,-0.5,0]);
        [warp{k}.com.poly1,warp{k}.com.gof.poly1] = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob','poly1');  
        [warp{k}.com.exp1,warp{k}.com.gof.exp1] = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob','exp1');  
        [warp{k}.com.power2,warp{k}.com.gof.power2] = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob','power2');
%         warp{k}.com.exp1const = fit(warp{k}.com.hist.edges(1:end-1)'+(warp{k}.com.hist.edges(2)-warp{k}.com.hist.edges(1))/2, warp{k}.com.hist.prob',fittype(@(a,b,c,x) a.*exp(b.*x) + c),'StartPoint',[0,-0.5,0]);

        % calculate linearity index
        warp{k}.com.linIdx1 = warp{k}.com.gof.poly1.rsquare / warp{k}.com.gof.exp1.rsquare;
        warp{k}.com.linIdx2 = (warp{k}.com.gof.poly1.rsquare - warp{k}.com.gof.exp1.rsquare) / warp{k}.com.gof.exp1.rsquare;
        warp{k}.com.linIdx3 = (warp{k}.com.gof.poly1.rsquare - warp{k}.com.gof.exp1.rsquare) / (warp{k}.com.gof.poly1.rsquare + warp{k}.com.gof.exp1.rsquare);
        warp{k}.com.linIdx4 = warp{k}.com.gof.poly1.rmse / warp{k}.com.gof.exp1.rmse;
        warp{k}.com.linIdx5 = (warp{k}.com.gof.poly1.rmse - warp{k}.com.gof.exp1.rmse) / warp{k}.com.gof.exp1.rmse;
        warp{k}.com.linIdx6 = (warp{k}.com.gof.poly1.rmse - warp{k}.com.gof.exp1.rmse) / (warp{k}.com.gof.poly1.rmse + warp{k}.com.gof.exp1.rmse);
        warp{k}.com.linIdx7 = warp{k}.com.gof.poly1.sse / warp{k}.com.gof.exp1.sse;
        warp{k}.com.linIdx8 = (warp{k}.com.gof.poly1.sse - warp{k}.com.gof.exp1.sse) / warp{k}.com.gof.exp1.sse;
        warp{k}.com.linIdx9 = (warp{k}.com.gof.poly1.sse - warp{k}.com.gof.exp1.sse) / (warp{k}.com.gof.poly1.sse + warp{k}.com.gof.exp1.sse);
        this_data_lin = warp{k}.com.gof.poly1.rmse.^2; this_data_nonlin = warp{k}.com.gof.exp1.rmse.^2;
        warp{k}.com.linIdx10 = this_data_lin / this_data_nonlin;
        warp{k}.com.linIdx11 = (this_data_lin - this_data_nonlin) / this_data_nonlin;
        warp{k}.com.linIdx12 = (this_data_lin - this_data_nonlin) / (this_data_lin + this_data_nonlin);
        warp{k}.com.nonLinIdx1 = warp{k}.com.gof.exp1.rsquare / warp{k}.com.gof.poly1.rsquare;
        warp{k}.com.nonLinIdx2 = (warp{k}.com.gof.exp1.rsquare - warp{k}.com.gof.poly1.rsquare) / warp{k}.com.gof.poly1.rsquare;
        warp{k}.com.nonLinIdx3 = (warp{k}.com.gof.exp1.rsquare - warp{k}.com.gof.poly1.rsquare) / (warp{k}.com.gof.poly1.rsquare + warp{k}.com.gof.exp1.rsquare);
        warp{k}.com.nonLinIdx4 = warp{k}.com.gof.exp1.rmse / warp{k}.com.gof.poly1.rmse;
        warp{k}.com.nonLinIdx5 = (warp{k}.com.gof.exp1.rmse - warp{k}.com.gof.poly1.rmse) / warp{k}.com.gof.poly1.rmse;
        warp{k}.com.nonLinIdx6 = (warp{k}.com.gof.exp1.rmse - warp{k}.com.gof.poly1.rmse) / (warp{k}.com.gof.poly1.rmse + warp{k}.com.gof.exp1.rmse);
        warp{k}.com.nonLinIdx7 = warp{k}.com.gof.exp1.sse / warp{k}.com.gof.poly1.sse;
        warp{k}.com.nonLinIdx8 = (warp{k}.com.gof.exp1.sse - warp{k}.com.gof.poly1.sse) / warp{k}.com.gof.poly1.sse;
        warp{k}.com.nonLinIdx9 = (warp{k}.com.gof.exp1.sse - warp{k}.com.gof.poly1.sse) / (warp{k}.com.gof.poly1.sse + warp{k}.com.gof.exp1.sse);
        this_data_lin = warp{k}.com.gof.poly1.rmse.^2; this_data_nonlin = warp{k}.com.gof.exp1.rmse.^2;
        warp{k}.com.nonLinIdx10 = this_data_nonlin / this_data_lin;
        warp{k}.com.nonLinIdx11 = (this_data_nonlin - this_data_lin) / this_data_lin;
        warp{k}.com.nonLinIdx12 = (this_data_nonlin - this_data_lin) / (this_data_lin + this_data_nonlin);
        
    end
end
    

%% Figure

% if strcmp(input_type,'all_stimVersion') || strcmp(input_type,'all')
%     
%     nrows = 1; ncols = 1; m=0;
%     F = default_figure([0,3.5,3,3]); hold on;
% 
%     temp = [];
% 
% %     m = m+1; subplot(nrows,ncols,m); hold on;
% %     plot(warp.peak.power2,warp.peak.hist.edges(1:end-1)+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob);
% %     legend('data','power2')
% %     title(['a=',num2str(warp.peak.power2.a,2),', b=',num2str(warp.peak.power2.b,2),', c=',num2str(warp.peak.power2.c,2)],'FontSize',9)
% %     temp = [temp; warp.peak.hist.prob];
% 
%     m = m+1; subplot(nrows,ncols,m); hold on;
%     h=plot(warp.peak.exp1,warp.peak.hist.edges(1:end-1)+(warp.peak.hist.edges(2)-warp.peak.hist.edges(1))/2, warp.peak.hist.prob);
%     h(1).Color = p.col.black; h(2).Color = p.col.black; legend('hide');
%     %legend('data','exp1')
%     %title(['a=',num2str(warp.peak.exp1.a,2),', b=',num2str(warp.peak.exp1.b,2)],'FontSize',9)
%     temp = [temp; warp.peak.hist.prob];
% 
%     for k=1:1
%         %subplot(nrows,ncols,k);
%         xlim([0,5.3])
%         xticks([0,5.3])
%         xticklabels({'0','5.3'})
%         ylim([0,0.2])
%         yticks([0,0.2])
%         xlabel('Time (s)')
%         ylabel('Proportion of sequence cells')
%     end
% %    savefig(F,'C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Examples\Frida (expert imaging, runner)\new\fit.fig');
% %    saveas(F,'C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Examples\Frida (expert imaging, runner)\new\fit.png');
% 
%     
% %     if strcmp(input_type,'all_stimVersion') 
% %         suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', sequence development (catch trials)'])
% %         drawnow;
% %         savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_all_stimVersion.fig']);
% %         saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_all_stimVersion.png']);
% %         disp(['--- Saved sequence warp analysis figure to ',path.filepart_outX,'plots.'])
% %     elseif strcmp(input_type,'all')
% %         suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', sequence development'])
% %         drawnow;
% %         savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_all.fig']);
% %         saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_all.png']);
% %         disp(['--- Saved sequence warp analysis figure to ',path.filepart_outX,'plots.'])
% %     end
% % end
% % 
% % if strcmp(input_type,'100t_stimVersion') || strcmp(input_type,'100t')
% %     
% %     nrows = 4; ncols = length(tng); m=0;
% %     F = default_figure();
% %     
% %     temp = [];
% %     for k=1:length(tng)
% %         m = m+1; subplot(nrows,ncols,m); hold on;
% %         plot(warp{k}.peak.poly1,warp{k}.peak.hist.edges(1:end-1)+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob);
% %         legend('data','poly1')
% %         title(['m=',num2str(warp{k}.peak.poly1.p1,2),', b=',num2str(warp{k}.peak.poly1.p2,2)],'FontSize',9)
% %         temp = [temp; warp{k}.peak.hist.prob];
% %     end
% %     for k=1:length(tng)
% %         m = m+1; subplot(nrows,ncols,m); hold on;
% %         plot(warp{k}.peak.power2,warp{k}.peak.hist.edges(1:end-1)+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob);
% %         legend('data','power2')
% %         title(['a=',num2str(warp{k}.peak.power2.a,2),', b=',num2str(warp{k}.peak.power2.b,2),', c=',num2str(warp{k}.peak.power2.c,2)],'FontSize',9)
% %         temp = [temp; warp{k}.peak.hist.prob];
% %     end
% %     for k=1:length(tng)
% %         m = m+1; subplot(nrows,ncols,m); hold on;
% %         plot(warp{k}.peak.exp1,warp{k}.peak.hist.edges(1:end-1)+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob);
% %         legend('data','exp1')
% %         title(['a=',num2str(warp{k}.peak.exp1.a,2),', b=',num2str(warp{k}.peak.exp1.b,2)],'FontSize',9)
% %         temp = [temp; warp{k}.peak.hist.prob];
% %     end
% %     for k=1:length(tng)
% %         m = m+1; subplot(nrows,ncols,m); hold on;
% %         plot(warp{k}.peak.exp1const,warp{k}.peak.hist.edges(1:end-1)+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob);
% %         legend('data','exp1const')
% %         title(['a=',num2str(warp{k}.peak.exp1const.a,2),', b=',num2str(warp{k}.peak.exp1const.b,2),', c=',num2str(warp{k}.peak.exp1const.c,2)],'FontSize',9)
% %         temp = [temp; warp{k}.peak.hist.prob];
% %     end
% %     for k=1:length(tng)*4
% %         subplot(nrows,ncols,k);
% %         xlim([0,5.3])
% %         ylim([0,nanmax(temp(:))])
% %         xlabel('Time (s)')
% %         ylabel('Proportion of sequence cells')
% %     end
% %     
% %     if strcmp(input_type,'100t_stimVersion') 
% %         suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', sequence development (catch trials) across 100t blocks'])
% %         drawnow;
% %         savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_100t_stimVersion.fig']);
% %         saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_100t_stimVersion.png']);
% %         disp(['--- Saved sequence warp analysis figure to ',path.filepart_outX,'plots.'])
% %     elseif strcmp(input_type,'100t')
% %         suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', sequence development across 100t blocks'])
% %         drawnow;
% %         savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_100t.fig']);
% %         saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_100t.png']);
% %         disp(['--- Saved sequence warp analysis figure to ',path.filepart_outX,'plots.'])
% %     end
% end
% 
% 
% 
% if strcmp(input_type,'100t_stimVersion') || strcmp(input_type,'100t')
%     
%     disp('SfN poster figure')
%     
%     nrows = 1; ncols = length(tng); m=0;
%     F = default_figure([-20,0.5,11.5,1.5]); %default_figure([-20,0.5,18,1.5]);
%     these_rgbs = discretisedColourMap('winter',false,length(tng));
% 
%     temp = [];
%     for k=1:length(tng)
%         m = m+1; subplot(nrows,ncols,m); hold on;
%         h=plot(warp{k}.peak.exp1,warp{k}.peak.hist.edges(1:end-1)+(warp{k}.peak.hist.edges(2)-warp{k}.peak.hist.edges(1))/2, warp{k}.peak.hist.prob);
%         h(1).Color = p.col.darkGray; h(2).Color = p.col.black; legend('hide'); % these_rgbs(k,:);
%         %legend('data','exp1')
%         title(['a=',num2str(warp{k}.peak.exp1.a,2),', b=',num2str(warp{k}.peak.exp1.b,2)],'FontSize',9)
%         temp = [temp; warp{k}.peak.hist.prob];
%     end
%     for k=1:length(tng)
%         subplot(nrows,ncols,k);
%         xlim([0,5.3])
%         xticks([0,5.3])
%         xticklabels({'0 s','5.3 s'})
%         xlabel('')
%         ylim([0,nanmax(temp(:))])
%         yticks([0,0.1])
%         yticklabels({'',''})
%         ylabel('')
%         if k==1
%             yticklabels({'0','0.1'})
%             ylabel('Proportion of sequence cells')
%         end
%     end
%     
%     if strcmp(input_type,'100t_stimVersion') 
%         %suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', sequence development (catch trials) across 100t blocks'])
%         drawnow;
%         savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_100t_stimVersion_cmpr.fig']);
%         saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_100t_stimVersion_cmpr.png']);
%         disp(['--- Saved sequence warp analysis figure to ',path.filepart_outX,'plots.'])
%     elseif strcmp(input_type,'100t')
%         %suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', sequence development across 100t blocks'])
%         drawnow;
%         savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_100t_cmpr.fig']);
%         saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','sequenceWarpAnalysis_100t_cmpr.png']);
%         disp(['--- Saved sequence warp analysis figure to ',path.filepart_outX,'plots.'])
%     end
% end


%% Save and return

if strcmp(input_type,'all')
    if strcmp(ops.warp.input,"tng")
        warp_all = warp;
        save([path.filepart_outX,'warp_all.mat'],'warp_all','-v7.3');
        disp(['--- Saved warp_all file as ',[path.filepart_outX,'warp_all.mat'],'.'])
    elseif strcmp(ops.warp.input,"tngn")
        warpn_all = warp;
        save([path.filepart_outX,'warpn_all.mat'],'warpn_all','-v7.3');
        disp(['--- Saved warpn_all file as ',[path.filepart_outX,'warpn_all.mat'],'.'])
    elseif strcmp(ops.warp.input,"tngnn")
        warpnn_all = warp;
        save([path.filepart_outX,'warpnn_all.mat'],'warpnn_all','-v7.3');
        disp(['--- Saved warpnn_all file as ',[path.filepart_outX,'warpnn_all.mat'],'.'])
    end
end
if strcmp(input_type,'all_stimVersion')
    if strcmp(ops.warp.input,"tng")
        warp_all_stimVersion = warp;
        save([path.filepart_outX,'warp_all_stimVersion.mat'],'warp_all_stimVersion','-v7.3');
        disp(['--- Saved warp_all_stimVersion file as ',[path.filepart_outX,'warp_all_stimVersion.mat'],'.'])
    elseif strcmp(ops.warp.input,"tngn")
        warpn_all_stimVersion = warp;
        save([path.filepart_outX,'warpn_all_stimVersion.mat'],'warpn_all_stimVersion','-v7.3');
        disp(['--- Saved warpn_all_stimVersion file as ',[path.filepart_outX,'warpn_all_stimVersion.mat'],'.'])
    elseif strcmp(ops.warp.input,"tngnn")
        warpnn_all_stimVersion = warp;
        save([path.filepart_outX,'warpnn_all_stimVersion.mat'],'warpnn_all_stimVersion','-v7.3');
        disp(['--- Saved warpnn_all_stimVersion file as ',[path.filepart_outX,'warpnn_all_stimVersion.mat'],'.'])
    end
end
if strcmp(input_type,'100t')
    if strcmp(ops.warp.input,"tng")
        warp_100t = warp;
        save([path.filepart_outX,'warp_100t.mat'],'warp_100t','-v7.3');
        disp(['--- Saved warp_100t file as ',[path.filepart_outX,'warp_100t.mat'],'.'])
    elseif strcmp(ops.warp.input,"tngn")
        warpn_100t = warp;
        save([path.filepart_outX,'warpn_100t.mat'],'warpn_100t','-v7.3');
        disp(['--- Saved warpn_100t file as ',[path.filepart_outX,'warpn_100t.mat'],'.'])
    elseif strcmp(ops.warp.input,"tngnn")
        warpnn_100t = warp;
        save([path.filepart_outX,'warpnn_100t.mat'],'warpnn_100t','-v7.3');
        disp(['--- Saved warpnn_100t file as ',[path.filepart_outX,'warpnn_100t.mat'],'.'])
    end
end
if strcmp(input_type,'100t_stimVersion')
    if strcmp(ops.warp.input,"tng")
        warp_100t_stimVersion = warp;
        save([path.filepart_outX,'warp_100t_stimVersion.mat'],'warp_100t_stimVersion','-v7.3');
        disp(['--- Saved warp_100t_stimVersion file as ',[path.filepart_outX,'warp_100t_stimVersion.mat'],'.'])
    elseif strcmp(ops.warp.input,"tngn")
        warpn_100t_stimVersion = warp;
        save([path.filepart_outX,'warpn_100t_stimVersion.mat'],'warpn_100t_stimVersion','-v7.3');
        disp(['--- Saved warpn_100t_stimVersion file as ',[path.filepart_outX,'warpn_100t_stimVersion.mat'],'.'])
    elseif strcmp(ops.warp.input,"tngnn")
        warpnn_100t_stimVersion = warp;
        save([path.filepart_outX,'warpnn_100t_stimVersion.mat'],'warpnn_100t_stimVersion','-v7.3');
        disp(['--- Saved warpnn_100t_stimVersion file as ',[path.filepart_outX,'warpnn_100t_stimVersion.mat'],'.'])
    end
end

if ops.close_figures
    close all;
end



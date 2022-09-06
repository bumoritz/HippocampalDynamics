%% Preprocessing_Stanage

% load data
load('D:\SniffinHippo\Repo\Stanage\Stanage_20210924\Stanage_20210924_paq_beh.mat')
load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_F_beh.mat')
load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_Fneu_beh.mat')
load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_background.mat');
load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_Fhalo_all_2.mat') % 2 or 3 doesn't really matter so use two because of faster processing
sync_beh = paq_beh;

% save_path
save_path = 'C:\SniffinHippo\Snaps\Preprocessing\Stanage\';

% crop background and halo data
Fbkgd_beh = background(:,108001:end-108000);
Fhalo_beh = Fhalo_all_2(:,108001:end-108000);
Fhbkgd_beh = nanmean(Fhalo_beh,1);

% subtract background from data
F_bs_beh = F_beh - Fbkgd_beh;
Fneu_bs_beh = Fneu_beh - Fbkgd_beh;
Fhalo_bs_beh = Fhalo_beh - Fbkgd_beh;

% subtract halo-background from data
F_hbs_beh = F_beh - Fhbkgd_beh;
Fneu_hbs_beh = Fneu_beh - Fhbkgd_beh;
Fhalo_hbs_beh = Fhalo_beh - Fhbkgd_beh;

% 4-tiled halo
% load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_Fhalo_all_2_4.mat')
% load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_Fhalo_all_2_1.mat')
% load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_Fhalo_all_2_2.mat')
% load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_Fhalo_all_2_3.mat')
% Fhalo_1_beh = Fhalo_all_2_1(:,108001:end-108000);
% Fhalo_2_beh = Fhalo_all_2_2(:,108001:end-108000);
% Fhalo_3_beh = Fhalo_all_2_3(:,108001:end-108000);
% Fhalo_4_beh = Fhalo_all_2_4(:,108001:end-108000);


%% Make trial-averaged traces

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[~,~,Fbkgd_beh_avg] = preprocessActivityMeasure(repmat(Fbkgd_beh,size(F_beh,1),1),p.inh,p,sync_beh,iscell);
Fbkgd_beh_avg = nanmean(Fbkgd_beh_avg,3);
[~,~,Fhbkgd_beh_avg] = preprocessActivityMeasure(repmat(Fhbkgd_beh,size(F_beh,1),1),p.inh,p,sync_beh,iscell);
Fhbkgd_beh_avg = nanmean(Fhbkgd_beh_avg,3);
[~,~,F_beh_avg] = preprocessActivityMeasure(F_beh,p.inh,p,sync_beh,iscell);
F_beh_avg = nanmean(F_beh_avg,3);
[~,~,F_bs_beh_avg] = preprocessActivityMeasure(F_bs_beh,p.inh,p,sync_beh,iscell);
F_bs_beh_avg = nanmean(F_bs_beh_avg,3);
[~,~,F_hbs_beh_avg] = preprocessActivityMeasure(F_hbs_beh,p.inh,p,sync_beh,iscell);
F_hbs_beh_avg = nanmean(F_hbs_beh_avg,3);
[~,~,Fhalo_beh_avg] = preprocessActivityMeasure(Fhalo_beh,p.inh,p,sync_beh,iscell);
Fhalo_beh_avg = nanmean(Fhalo_beh_avg,3);
[~,~,Fhalo_bs_beh_avg] = preprocessActivityMeasure(Fhalo_bs_beh,p.inh,p,sync_beh,iscell);
Fhalo_bs_beh_avg = nanmean(Fhalo_bs_beh_avg,3);
[~,~,Fhalo_hbs_beh_avg] = preprocessActivityMeasure(Fhalo_hbs_beh,p.inh,p,sync_beh,iscell);
Fhalo_hbs_beh_avg = nanmean(Fhalo_hbs_beh_avg,3);
[~,~,Fneu_beh_avg] = preprocessActivityMeasure(Fneu_beh,p.inh,p,sync_beh,iscell);
Fneu_beh_avg = nanmean(Fneu_beh_avg,3);
[~,~,Fneu_bs_beh_avg] = preprocessActivityMeasure(Fneu_bs_beh,p.inh,p,sync_beh,iscell);
Fneu_bs_beh_avg = nanmean(Fneu_bs_beh_avg,3);
[~,~,Fneu_hbs_beh_avg] = preprocessActivityMeasure(Fneu_hbs_beh,p.inh,p,sync_beh,iscell);
Fneu_hbs_beh_avg = nanmean(Fneu_hbs_beh_avg,3);


%% Calculate factors for halo-subtraction from background-subtracted traces based on averages

% % rho method
% subtractionFactors = -100:0.1:100;
% rho = nan(size(F_bs_beh_avg,1),length(subtractionFactors));
% for i=1:size(F_bs_beh_avg,1)
%     for j=1:length(subtractionFactors)
%         temp = subtractionFactors(j);
%         rho(i,j) = corr((F_bs_beh_avg(i,:)-temp*Fhalo_bs_beh_avg(i,:))',Fhalo_bs_beh_avg(i,:)','Type','Pearson','Rows','Complete');
%     end
%     i
% end
% optimalSubtractionFactor = nan(size(F_bs_beh_avg,1),1);
% [temp1,temp2] = nanmin(abs(rho),[],2);
% for i=1:size(F_bs_beh_avg,1)
%     if ~isnan(temp1(i)) && ~isnan(temp2(i))
%         optimalSubtractionFactor(i) = subtractionFactors(temp2(i));
%     end
% end

% robust regression method
optimalSubtractionFactor = nan(size(F_bs_beh_avg,1),1);
for i=1:size(F_bs_beh_avg,1)
    if iscell(i)==1
        temp = robustfit(Fhalo_bs_beh_avg(i,:),F_bs_beh_avg(i,:));
        optimalSubtractionFactor(i) = temp(2);
    end
end

% robust regression method
optimalSubtractionFactor_fromTraces = nan(size(F_bs_beh,1),1);
for i=1:size(F_bs_beh,1)
    if iscell(i)==1
        temp = robustfit(Fhalo_bs_beh(i,:),F_bs_beh(i,:));
        optimalSubtractionFactor_fromTraces(i) = temp(2);
    end
end

% robust regression method -hbs
optimalSubtractionFactor_hbs = nan(size(F_hbs_beh_avg,1),1);
for i=1:size(F_hbs_beh_avg,1)
    if iscell(i)==1
        temp = robustfit(Fhalo_hbs_beh_avg(i,:),F_hbs_beh_avg(i,:));
        optimalSubtractionFactor_hbs(i) = temp(2);
    end
end

% same for not background subtracted

% rho_raw = nan(size(F_beh_avg,1),length(subtractionFactors));
% for i=1:size(F_beh_avg,1)
%     for j=1:length(subtractionFactors)
%         temp = subtractionFactors(j);
%         rho_raw(i,j) = corr((F_beh_avg(i,:)-temp*Fhalo_beh_avg(i,:))',Fhalo_beh_avg(i,:)','Type','Pearson','Rows','Complete');
%     end
%     i
% end
% 
% optimalSubtractionFactor_raw = nan(size(F_beh_avg,1),1);
% [temp1,temp2] = nanmin(abs(rho_raw),[],2);
% for i=1:size(F_beh_avg,1)
%     if ~isnan(temp1(i)) && ~isnan(temp2(i))
%         optimalSubtractionFactor_raw(i) = subtractionFactors(temp2(i));
%     end
% end

% robust regression method - raw
optimalSubtractionFactor_raw = nan(size(F_beh_avg,1),1);
for i=1:size(F_beh_avg,1)
    if iscell(i)==1
        temp = robustfit(Fhalo_beh_avg(i,:),F_beh_avg(i,:));
        optimalSubtractionFactor_raw(i) = temp(2);
    end
end


%% Do halo-subtraction from background-subtracted traces

F_hs_beh = F_beh - optimalSubtractionFactor_raw .* Fhalo_beh;
F_hs_bs_beh_fromTraces = F_bs_beh - optimalSubtractionFactor_fromTraces .* Fhalo_bs_beh;
F_hs_bs_beh = F_bs_beh - optimalSubtractionFactor .* Fhalo_bs_beh;
F_hs_hbs_beh = F_hbs_beh - optimalSubtractionFactor_hbs .* Fhalo_hbs_beh;

[~,~,F_hs_beh_avg] = preprocessActivityMeasure(F_hs_beh,p.inh,p,sync_beh,iscell);
F_hs_beh_avg = nanmean(F_hs_beh_avg,3);
[~,~,F_hs_bs_beh_avg] = preprocessActivityMeasure(F_hs_bs_beh,p.inh,p,sync_beh,iscell);
F_hs_bs_beh_avg = nanmean(F_hs_bs_beh_avg,3);
[~,~,F_hs_hbs_beh_avg] = preprocessActivityMeasure(F_hs_hbs_beh,p.inh,p,sync_beh,iscell);
F_hs_hbs_beh_avg = nanmean(F_hs_hbs_beh_avg,3);

[~,~,F_hs_bs_beh_fromTraces_avg] = preprocessActivityMeasure(F_hs_bs_beh_fromTraces,p.inh,p,sync_beh,iscell);
F_hs_bs_beh_fromTraces_avg = nanmean(F_hs_bs_beh_fromTraces_avg,3);

% F_ns_beh = F_beh - 1 .* Fneu_beh;
% F_ns_bs_beh = F_bs_beh - 1 .* Fneu_bs_beh;
% 
% [~,~,F_ns_beh_avg] = preprocessActivityMeasure(F_ns_beh,p.inh,p,sync_beh,iscell);
% F_ns_beh_avg = nanmean(F_ns_beh_avg,3);
% [~,~,F_ns_bs_beh_avg] = preprocessActivityMeasure(F_ns_bs_beh,p.inh,p,sync_beh,iscell);
% F_ns_bs_beh_avg = nanmean(F_ns_bs_beh_avg,3);


%% --- FIGURES ---

%% Example raw traces

idx = 1;

wdw = 1:10000;
left = 1;
right = 10000;

F=default_figure([20,0.5,20,9.9]);

subplot(4,1,1); hold on;
plot(F_beh(idx,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['F (idx=',num2str(idx),')'])

subplot(4,1,2); hold on;
plot(Fneu_beh(idx,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fneu (idx=',num2str(idx),')'])

subplot(4,1,3); hold on;
plot(Fhalo_beh(idx,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fhalo (idx=',num2str(idx),')'])

subplot(4,1,4); hold on;
plot(Fbkgd_beh(:,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fbkgd'])

saveas(F,[save_path,'exampleRawTraces_idx',num2str(idx),'.png']);


%% Example background-subtracted traces

idx = 5;

wdw = 1:10000;
left = 1;
right = 10000;

F=default_figure([20,0.5,20,9.9]);

subplot(5,1,1); hold on;
plot(Fbkgd_beh(:,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fbkgd'])

subplot(5,1,2); hold on;
plot(F_beh(idx,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['F (idx=',num2str(idx),')'])

subplot(5,1,3); hold on;
plot(F_bs_beh(idx,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['F-bs (idx=',num2str(idx),')'])

subplot(5,1,4); hold on;
plot(Fhalo_beh(idx,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fhalo (idx=',num2str(idx),')'])

subplot(5,1,5); hold on;
plot(Fhalo_bs_beh(idx,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fhalo-bs (idx=',num2str(idx),')'])

saveas(F,[save_path,'exampleBkgdSubTraces_idx',num2str(idx),'.png']);


%% Example halo-subtracted background-subtracted average traces

for i=1:30
    temp=find(iscell==1);
    idx = temp(i);

    nrows = 4; ncols = 3;
    F = default_figure([20,0.5,20,9.9]);

    r=2; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(F_beh_avg(idx,:))
    set(gca,'LineWidth',3)
    title(['F'])

    r=3; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(Fhalo_beh_avg(idx,:))
    title(['Fhalo'])

    r=4; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(F_hs_beh_avg(idx,:))
    title(['F-',num2str(optimalSubtractionFactor_raw(idx),2),'*Fhalo'])

    r=1; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(Fbkgd_beh_avg(idx,:))
    title(['Fbkgd'])

    r=2; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(F_bs_beh_avg(idx,:))
    title(['F(-bkgd)'])

    r=3; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(Fhalo_bs_beh_avg(idx,:))
    title(['Fhalo(-bkgd)'])

    r=4; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(F_hs_bs_beh_avg(idx,:))
    set(gca,'LineWidth',3)
    title(['F(-bkgd)-',num2str(optimalSubtractionFactor(idx),2),'*Fhalo(-bkgd)'])

    r=1; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(Fbkgd_beh_avg(idx,:))
    title(['Fbkgd'])

    r=2; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(F_bs_beh_avg(idx,:))
    title(['F(-bkgd)'])

    r=3; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(Fhalo_bs_beh_avg(idx,:))
    title(['Fhalo(-bkgd)'])

    r=4; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
    plot(F_hs_bs_beh_fromTraces_avg(idx,:))
    set(gca,'LineWidth',3)
    title(['F(-bkgd)-',num2str(optimalSubtractionFactor_fromTraces(idx),2),'*Fhalo(-bkgd), from Traces'])

    %%%

    % r=2; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
    % plot(F_beh_avg(idx,:))
    % set(gca,'LineWidth',3)
    % title(['F'])
    % 
    % r=3; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
    % plot(Fneu_beh_avg(idx,:))
    % title(['Fneu'])
    % 
    % r=4; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
    % plot(F_ns_beh_avg(idx,:))
    % title(['F-',num2str(1),'*Fneu'])
    % 
    % r=1; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
    % plot(Fbkgd_beh_avg(idx,:))
    % title(['Fbkgd'])
    % 
    % r=2; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
    % plot(F_bs_beh_avg(idx,:))
    % title(['F(-bkgd)'])
    % 
    % r=3; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
    % plot(Fneu_bs_beh_avg(idx,:))
    % title(['Fneu(-bkgd)'])
    % 
    % r=4; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
    % plot(F_ns_bs_beh_avg(idx,:))
    % set(gca,'LineWidth',3)
    % title(['F(-bkgd)-',num2str(1),'*Fneu(-bkgd)'])

    suptitle(['idx=',num2str(idx)])

    %saveas(F,[save_path,'exampleHaloSubAvgTraces_idx',num2str(idx),'.png']);
end






%% 

idx = 6;

[b, stats] = robustfit(Fhalo_bs_beh_avg(idx,:),F_bs_beh_avg(idx,:));
sub_trace = F_bs_beh_avg(idx,:) - b(2) * Fhalo_bs_beh_avg(idx,:);


sub_trace2 = F_bs_beh_avg(idx,:) - b(2) * Fhalo_bs_beh_avg(idx,:) - b(1);

%F_bs_beh_avg(idx,:) = b(1) + b(2) * Fhalo_bs_beh_avg(idx,:);



figure; hold on;
plot(F_bs_beh_avg(idx,:),'b')
plot(Fhalo_bs_beh_avg(idx,:),'r')
plot(sub_trace,'g')
plot(sub_trace2,'g:')

legend('F-bs','Fhalo-bs','subtracted (slope only)')
suptitle(['idx=',num2str(idx)])



%%
idx = 6;
%[b, stats] = robustfit(Fhalo_bs_beh_avg(idx,:),F_bs_beh_avg(idx,:));
% coeffs = estimateNeuropilCoefficients(F_bs_beh_avg(idx,:), Fhalo_bs_beh_avg(idx,:))



%% Transform data

% how about somehow minimising diff

% order of transformations
% - shift baseline of both neuropil signal and F signal to zero
% - scale neuropil signal by ideal factor to match shape of F (by correlation or absolute deviation or LSE or fit?)
% - subtract scaled neuropil signal from F
% - potentially scale whole baseline back up to where it was

idx = 6;
bwdw = 1:10;

% define data
this_F = F_bs_beh_avg(idx,:);
this_Fhalo = Fhalo_bs_beh_avg(idx,:);

% shift traces to zero
this_zeroshift_F = nanmean(this_F(bwdw));
this_zeroshift_Fhalo = nanmean(this_Fhalo(bwdw));
this_F_shifted = this_F - this_zeroshift_F;
this_Fhalo_shifted = this_Fhalo - this_zeroshift_Fhalo;

% identify appropriate scaling factor
scalingFactors = -20:0.1:20;
loss = nan(1,length(scalingFactors));
for j=1:length(scalingFactors)
    scalingFactor = scalingFactors(j);
    loss(j) = nanmean((this_F_shifted-scalingFactor*this_Fhalo_shifted) - this_Fhalo_shifted);
end
[~,temp]=nanmin(abs(loss));
optimalSubtractionFactor = scalingFactors(temp);

% neuropil trace was scaled to lead to the 








figure; hold on;
plot(this_F_shifted,'b')
plot(this_Fhalo_shifted,'r')
plot(this_F_shifted-optimalSubtractionFactor*this_Fhalo_shifted,'g')
% plot(sub_trace,'g')
% plot(sub_trace2,'g:')

legend('F-bs','Fhalo-bs')
suptitle(['idx=',num2str(idx)])



%%

this_dFF = dFF_hs_bs_beh;
[c, s, options] = deconvolveCa(this_dFF);


%writeNPY(dFF_hs_bs_beh,['E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\dFF_hs_bs_beh.npy']);







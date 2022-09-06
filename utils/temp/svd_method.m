%% Load data

% load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_F_beh.mat');
% F = F_beh;
% load('D:\SniffinHippo\Repo\Stanage\Stanage_20210924\Stanage_20210924_s2p_meta.mat');
% iscell = s2p_meta.iscell(:,1);
% load('D:\SniffinHippo\Repo\Stanage\Stanage_20210924\Stanage_20210924_paq_beh.mat')
% sync_beh = paq_beh;

% load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210925\Stanage_20210925_F_beh.mat');
% F = F_beh;
% load('D:\SniffinHippo\Repo\Stanage\Stanage_20210925\Stanage_20210925_s2p_meta.mat');
% iscell = s2p_meta.iscell(:,1);
% load('D:\SniffinHippo\Repo\Stanage\Stanage_20210925\Stanage_20210925_paq_beh.mat')
% sync_beh = paq_beh;

% load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210927\Stanage_20210927_F_beh.mat');
% F = F_beh;
% load('D:\SniffinHippo\Repo\Stanage\Stanage_20210927\Stanage_20210927_s2p_meta.mat');
% iscell = s2p_meta.iscell(:,1);
% load('D:\SniffinHippo\Repo\Stanage\Stanage_20210927\Stanage_20210927_paq_beh.mat')
% sync_beh = paq_beh;

% load('E:\SniffinHippo\RepoX\Carlo\Carlo_20210313\Carlo_20210313_F_beh.mat');
% F = F_beh;
% load('D:\SniffinHippo\Repo\Carlo\Carlo_20210313\Carlo_20210313_s2p_meta.mat');
% iscell = s2p_meta.iscell(:,1);
% load('D:\SniffinHippo\Repo\Carlo\Carlo_20210313\Carlo_20210313_paq_beh.mat')
% sync_beh = paq_beh;

% -----


% spks
s2p = load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\suite2p\plane0\Fall.mat');
numFrames = size(s2p.F,2);
F_.unblocked_765 = s2p.F(:,1:25333); % 1, 25333
F_.blocked_765 = s2p.F(:,25333+1:25333+23676); % 2, 23676
F_.unblocked_930 = s2p.F(:,25333+23676+1:25333+23676+22739); % 3, 22739
F_.blocked_930 = s2p.F(:,25333+23676+22739+1:25333+23676+22739+22724); % 4, 22724
iscell = s2p.iscell(:,1);

% sniffinSyncs
temp = load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_765_blocked_sniffinSyncs.mat');
sniffinSyncs.blocked_765 = temp.sniffinSyncs;
temp = load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_765_unblocked_sniffinSyncs.mat');
sniffinSyncs.unblocked_765 = temp.sniffinSyncs;
temp = load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_930_blocked_sniffinSyncs.mat');
sniffinSyncs.blocked_930 = temp.sniffinSyncs;
temp = load('C:\Data\SniffinHippo\InhibitionControl\Austin_20220726_765ctrl\Behaviour\Austin_20220726_930_unblocked_sniffinSyncs.mat');
sniffinSyncs.unblocked_930 = temp.sniffinSyncs;

%%

F = F_.blocked_930; sync_beh.sync = sniffinSyncs.blocked_930; 


%% Do SVD on F data

data_in = single(F(find(iscell==1),:));
M = data_in';
[U,S,V] = svdecon(M);
M_hat = U*S*V';
data_out = M_hat';


%% Calculate zero crossing rates of temporal singular vectors

zcr = nansum(abs(diff(U>=0,[],1)),1);
zcrd = abs(diff(zcr));
zcrd_threshold = 3*nanstd(zcrd);
singularValuesToKeep = ones(1,length(zcr));
singularValuesToKeep([true,zcrd>zcrd_threshold]) = 0;

figure; hold on;
subplot(2,1,1)
plot(zcr)
xlabel('Temporal singular vectors')
ylabel('Zero crossing rate')
subplot(2,1,2)
plot(zcrd)
yline(zcrd_threshold,':');
xlabel('Temporal singular vectors')
ylabel('Absolute derivative of zero crossing rate')
title(['Discarding ',num2str(sum(1-singularValuesToKeep)),' components'])


%% Reconstruct data

S_wo1 = S;
S_wo1(:,1) = 0;
M_hat_wo1 = U*S_wo1*V';
data_out_wo1 = M_hat_wo1';

S_wo2 = S;
S_wo2(:,1:2) = 0;
M_hat_wo2 = U*S_wo2*V';
data_out_wo2 = M_hat_wo2';

S_wo3 = S;
S_wo3(:,1:3) = 0;
M_hat_wo3 = U*S_wo3*V';
data_out_wo3 = M_hat_wo3';

S_wo4 = S;
S_wo4(:,1:4) = 0;
M_hat_wo4 = U*S_wo4*V';
data_out_wo4 = M_hat_wo4';

S_wo5 = S;
S_wo5(:,1:5) = 0;
M_hat_wo5 = U*S_wo5*V';
data_out_wo5 = M_hat_wo5';

S_zcrFiltered = S;
S_zcrFiltered(:,find(singularValuesToKeep==0)) = 0;
M_hat_zcrFiltered = U*S_zcrFiltered*V';
data_out_zcrFiltered = M_hat_zcrFiltered';


%% Calculate trial-averages

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[~,~,data_in_avg] = preprocessActivityMeasure(data_in,p.inh,p,sync_beh,ones(size(data_out,1),1));
data_in_avg = nanmean(data_in_avg,3);
[~,~,data_out_avg] = preprocessActivityMeasure(data_out,p.inh,p,sync_beh,ones(size(data_out,1),1));
data_out_avg = nanmean(data_out_avg,3);
[~,~,data_out_wo1_avg] = preprocessActivityMeasure(data_out_wo1,p.inh,p,sync_beh,ones(size(data_out_wo1,1),1));
data_out_wo1_avg = nanmean(data_out_wo1_avg,3);
[~,~,data_out_wo2_avg] = preprocessActivityMeasure(data_out_wo2,p.inh,p,sync_beh,ones(size(data_out_wo2,1),1));
data_out_wo2_avg = nanmean(data_out_wo2_avg,3);
[~,~,data_out_wo3_avg] = preprocessActivityMeasure(data_out_wo3,p.inh,p,sync_beh,ones(size(data_out_wo3,1),1));
data_out_wo3_avg = nanmean(data_out_wo3_avg,3);
[~,~,data_out_wo4_avg] = preprocessActivityMeasure(data_out_wo4,p.inh,p,sync_beh,ones(size(data_out_wo4,1),1));
data_out_wo4_avg = nanmean(data_out_wo4_avg,3);
[~,~,data_out_wo5_avg] = preprocessActivityMeasure(data_out_wo5,p.inh,p,sync_beh,ones(size(data_out_wo5,1),1));
data_out_wo5_avg = nanmean(data_out_wo5_avg,3);
[~,~,data_out_zcrFiltered_avg] = preprocessActivityMeasure(data_out_zcrFiltered,p.inh,p,sync_beh,ones(size(data_out_wo5,1),1));
data_out_zcrFiltered_avg = nanmean(data_out_zcrFiltered_avg,3);
[~,~,U_avg] = preprocessActivityMeasure(U',p.inh,p,sync_beh,ones(size(data_out_wo5,1),1));
U_avg = nanmean(U_avg,3);
[~,~,U1_avg] = preprocessActivityMeasure(repmat(U(:,1)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
U1_avg = nanmean(U1_avg,3);
[~,~,U2_avg] = preprocessActivityMeasure(repmat(U(:,2)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
U2_avg = nanmean(U2_avg,3);
[~,~,U3_avg] = preprocessActivityMeasure(repmat(U(:,3)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
U3_avg = nanmean(U3_avg,3);
[~,~,U4_avg] = preprocessActivityMeasure(repmat(U(:,4)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
U4_avg = nanmean(U4_avg,3);
[~,~,U5_avg] = preprocessActivityMeasure(repmat(U(:,5)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
U5_avg = nanmean(U5_avg,3);
[~,~,U6_avg] = preprocessActivityMeasure(repmat(U(:,6)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
U6_avg = nanmean(U6_avg,3);
[~,~,U7_avg] = preprocessActivityMeasure(repmat(U(:,7)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
U7_avg = nanmean(U7_avg,3);
[~,~,U8_avg] = preprocessActivityMeasure(repmat(U(:,8)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
U8_avg = nanmean(U8_avg,3);
[~,~,U9_avg] = preprocessActivityMeasure(repmat(U(:,9)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
U9_avg = nanmean(U9_avg,3);
[~,~,U10_avg] = preprocessActivityMeasure(repmat(U(:,10)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
U10_avg = nanmean(U10_avg,3);


%% Plot data that was reconstructed based on individual trials

nrows = 4; ncols = 5;
default_figure();

r=1; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_in_avg',[],1))
xlabel('Cell (n)')
ylabel('Frame (m)')
title('M (z-scored)')

r=1; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U_avg(1:10,:)')
title('Left (temporal) singular vectors U')
xlabel('Component (m)')
ylabel('Frame (m)')

r=1; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(S(1:10,1:10)')
title('Singular values S')
xlabel('Component (n)')
ylabel('Component (m)')

r=1; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(V',[0,0.05])
title('Right ("spatial") singular vectors V^T')
xlabel('Cell or Component (n)')
ylabel('Cell or Component (n)')

r=1; c=5; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_out_avg',[],1))
xlabel('Cell (n)')
ylabel('Frame (m)')
title('M hat (z-scored)')

r=2; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_out_wo1_avg',[],1))
xlabel('Cell')
ylabel('Frame')
title('M hat (-1) (z-scored)')

r=2; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_out_wo2_avg',[],1))
xlabel('Cell')
ylabel('Frame')
title('M hat (-2) (z-scored)')

r=2; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_out_wo3_avg',[],1))
xlabel('Cell')
ylabel('Frame')
title('M hat (-3) (z-scored)')

r=2; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_out_wo4_avg',[],1))
xlabel('Cell')
ylabel('Frame')
title('M hat (-4) (z-scored)')

r=2; c=5; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(zscore(data_out_zcrFiltered_avg',[],1))
xlabel('Cell')
ylabel('Frame')
title('M hat (zcr filtered) (z-scored)')

r=3; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U1_avg(1,:)')
%xlabel('Cell')
ylabel('Frame')
title('U(1)')

r=3; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U2_avg(1,:)')
%xlabel('Cell')
ylabel('Frame')
title('U(2)')

r=3; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U3_avg(1,:)')
%xlabel('Cell')
ylabel('Frame')
title('U(3)')

r=3; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U4_avg(1,:)')
%xlabel('Cell')
ylabel('Frame')
title('U(4)')

r=3; c=5; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U5_avg(1,:)')
%xlabel('Cell')
ylabel('Frame')
title('U(5)')

r=4; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U6_avg(1,:)')
%xlabel('Cell')
ylabel('Frame')
title('U(6)')

r=4; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U7_avg(1,:)')
%xlabel('Cell')
ylabel('Frame')
title('U(7)')

r=4; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U8_avg(1,:)')
%xlabel('Cell')
ylabel('Frame')
title('U(8)')

r=4; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U9_avg(1,:)')
%xlabel('Cell')
ylabel('Frame')
title('U(9)')

r=4; c=5; subplot(nrows,ncols,(r-1)*ncols+c);
imagesc(U10_avg(1,:)')
%xlabel('Cell')
ylabel('Frame')
title('U(10)')

suptitle('Austin, 765 nm, blocked')


%%

% reinsert notiscells
data_in_full = nan(size(F));
data_in_full(find(iscell==1),:) = data_in;
data_out_zcrFiltered_full = nan(size(F));
data_out_zcrFiltered_full(find(iscell==1),:) = data_out_zcrFiltered;










%%

% temp = data_in - nanmean(data_in,1);
% [~,~,Fcar_avg] = preprocessActivityMeasure(temp,p.inh,p,sync_beh,ones(size(U,2),1));
% Fcar_avg = nanmean(Fcar_avg,3);
% 
% temp0 = zscore(data_in,[],2);
% temp = temp0 - nanmean(temp0,1);
% [~,~,Fcarz_avg] = preprocessActivityMeasure(temp,p.inh,p,sync_beh,ones(size(U,2),1));
% Fcarz_avg = nanmean(Fcarz_avg,3);
% 
% 
% %%
% 
% figure
% 
% subplot(3,2,1)
% imagesc(data_in_avg,[prctile(data_in_avg(:),5),prctile(data_in_avg(:),95)])
% title('F')
% subplot(3,2,2)
% imagesc(Fcar_avg,[prctile(Fcar_avg(:),5),prctile(Fcar_avg(:),95)])
% title('F - common average ref')
% subplot(3,2,3)
% imagesc(zscore(data_in_avg,[],2))
% title('F (z-scored for plotting)')
% subplot(3,2,4)
% imagesc(zscore(Fcar_avg,[],2))
% title('F - common average ref (z-scored for plotting)')
% subplot(3,2,6)
% imagesc(zscore(Fcar_avg,[],2))
% title('F - common average ref (on z-scored)')
% 
% 
% %% Do ICA on F data
% 
% data_in = single(F(find(iscell==1),:));
% M = data_in';
% 
% mdl = rica(M,30);
% M_ica5 = transform(mdl,M);
% 
% [~,~,M_ica_avg1] = preprocessActivityMeasure(repmat(M_ica5(:,1)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
% M_ica_avg1 = nanmean(M_ica_avg1,3);
% [~,~,M_ica_avg2] = preprocessActivityMeasure(repmat(M_ica5(:,2)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
% M_ica_avg2 = nanmean(M_ica_avg2,3);
% [~,~,M_ica_avg3] = preprocessActivityMeasure(repmat(M_ica5(:,3)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
% M_ica_avg3 = nanmean(M_ica_avg3,3);
% [~,~,M_ica_avg4] = preprocessActivityMeasure(repmat(M_ica5(:,4)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
% M_ica_avg4 = nanmean(M_ica_avg4,3);
% [~,~,M_ica_avg5] = preprocessActivityMeasure(repmat(M_ica5(:,5)',size(U,2),1),p.inh,p,sync_beh,ones(size(U,2),1));
% M_ica_avg5 = nanmean(M_ica_avg5,3);
% 
% 
% %%
% 
% figure;
% subplot(1,5,1)
% imagesc(M_ica_avg1(1,:))
% title('IC1')
% subplot(1,5,2)
% imagesc(M_ica_avg2(1,:))
% title('IC2')
% subplot(1,5,3)
% imagesc(M_ica_avg3(1,:))
% title('IC3')
% subplot(1,5,4)
% imagesc(M_ica_avg4(1,:))
% title('IC4')
% subplot(1,5,5)
% imagesc(M_ica_avg5(1,:))
% title('IC5')
% 

%%


%%







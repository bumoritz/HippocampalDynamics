seqNMFstuff

%% seqNMF preparations

rng(p.general.rgnSeed);
K = 2; % somehow related to number of sequences
L = length(p.general.bins_analysisWindow);
lambda =.008;


%% seqNMF execution

this_window = 1:5000;

shg; clf
[W,H] = seqNMF(nf_binned(find(iscell),this_window),'K',K, 'L', L,'lambda', lambda);


% [W,H,cost,loading,power] = seqNMF(nf_binned(find(iscell),:), 'K',K,'L',L,'lambda',lambda);


%% 

figure;
imagesc(nf_binned(find(iscell),5000:5500));
colormap('parula')
colorbar

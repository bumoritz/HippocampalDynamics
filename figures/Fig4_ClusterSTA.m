%% Fig4_ClusterSTA

% for cluster response maps
% import data using Analyses_Master with ops.do_responseAnalysis = true;
% use Philip_20211004 (seq) with cluster 34 (8 is nice as well)

% for STA image
% import data using Data2repo_Master with ops.do_zcorrData = true;
% use Philip_20211004 (seq) with cluster 34 (34, but 8 is nice as well)

load([path.filepart_out,'resp.mat']);

save_root_fig = [path.root_summary,'figures\Fig4_fig\'];
save_root_png = [path.root_summary,'figures\Fig4_png\'];
save_root_pdf = [path.root_summary,'figures\Fig4_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig4_txt\'];


%% Fig4_ClusterSTA_ResponseMap

clusterID = 8;
F = clusterResponseMap(s2p_meta,iscell,trg,resp,clusterID,0,0,[],1);

% save plot
savefig(F,[save_root_fig,'\Fig4_ClusterSTA_ResponseMap.fig']);
saveas(F,[save_root_png,'\Fig4_ClusterSTA_ResponseMap.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ClusterSTA_ResponseMap.pdf']); set(gcf,'Color',[1,1,1])


%% Fig4_ClusterSTA_ResponseMapWithRadius

clusterID = 8;
F = clusterResponseMap(s2p_meta,iscell,trg,resp,clusterID,0,0,[],1,1,info);

% save plot
savefig(F,[save_root_fig,'\Fig4_ClusterSTA_ResponseMapWithRadius.fig']);
saveas(F,[save_root_png,'\Fig4_ClusterSTA_ResponseMapWithRadius.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ClusterSTA_ResponseMapWithRadius.pdf']); set(gcf,'Color',[1,1,1])


%% Fig4_ClusterSTA_ResponseMap_Colorbar

clusterID = 8;
F = clusterResponseMap(s2p_meta,iscell,trg,resp,clusterID,0,0,[],0);

savefig(F,[save_root_fig,'\Fig4_ClusterSTA_ResponseMap_Colorbar.fig']);
saveas(F,[save_root_png,'\Fig4_ClusterSTA_ResponseMap_Colorbar.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ClusterSTA_ResponseMap_Colorbar.pdf']); set(gcf,'Color',[1,1,1])

%% --------------------------------

%% --- Fig4_ClusterSTA_STAImage ---

path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig4_fig\'];
save_root_png = [path.root_summary,'figures\Fig4_png\'];
save_root_pdf = [path.root_summary,'figures\Fig4_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig4_txt\'];

% load data
[path,raw] = data2repo_gatherRawDataInfo(info,path);
load([path.filepart_out,'task.mat']);
load([path.filepart_out,'s2p_meta.mat']);
load([path.filepart_out,'paq_beh.mat']);
load([path.filepart_out,'thor_beh.mat']);
load([path.filepart_out,'trg_rigid.mat']);
iscell = s2p_meta.iscell;

path.imagingFile    = ['G:\Data\2021\2021-10\2021-10-04\Philip\Imaging\registered_movie.raw'];
skipFirstXFrames = 108000;

% select data (first cluster of seq stim)
numAvgFrames = 3;
[~,these_stimTrials,~] = intersect(find(task.var==1),find(task.odour1=='X'));

% make STA
frames_base = nan(numAvgFrames,length(these_stimTrials));
frames_resp = nan(numAvgFrames,length(these_stimTrials));
for i=1:length(these_stimTrials)
    frames_base(:,i) = (thor_beh.stimPatternOnset(1,these_stimTrials(i))-numAvgFrames:thor_beh.stimPatternOnset(1,these_stimTrials(i))-1) + skipFirstXFrames;
    frames_resp(:,i) = (thor_beh.stimPatternComplete(1,these_stimTrials(i))+1:thor_beh.stimPatternComplete(1,these_stimTrials(i))+numAvgFrames) + skipFirstXFrames;
end
images_base = zeros(512,512,length(these_stimTrials),'uint16');
images_resp = zeros(512,512,length(these_stimTrials),'uint16');
images = zeros(512,512,length(these_stimTrials),'uint16');
sta_base = zeros(512,512,'uint16');
sta_resp = zeros(512,512,'uint16');
sta = zeros(512,512,'uint16');
for i=1:length(these_stimTrials)
    
    % base
    this_fid = fopen(path.imagingFile, 'r');
    fseek(this_fid,0,'bof');
    these_frames = zeros(512,512,numAvgFrames,'uint16');
    for j=1:numAvgFrames
        fseek(this_fid,(frames_base(j,i)-1) *512*512*2,'bof');
        temp = uint16(fread(this_fid,512*512,'uint16',0,'l'));
        these_frames(:,:,j) = reshape(temp,512,512);
        frewind(this_fid);
    end
    fclose(this_fid);
    images_base(:,:,i) = nanmean(these_frames,3);
    
    % resp
    this_fid = fopen(path.imagingFile, 'r');
    fseek(this_fid,0,'bof');
    these_frames = zeros(512,512,numAvgFrames,'uint16');
    for j=1:numAvgFrames
        fseek(this_fid,(frames_resp(j,i)-1) *512*512*2,'bof');
        temp = uint16(fread(this_fid,512*512,'uint16',0,'l'));
        these_frames(:,:,j) = reshape(temp,512,512);
        frewind(this_fid);
    end
    fclose(this_fid);
    images_resp(:,:,i) = nanmean(these_frames,3);

    images(:,:,i) = images_resp(:,:,i) - images_base(:,:,i);
end
sta_base(:,:) = nanmean(images_base,3);
sta_resp(:,:) = nanmean(images_resp,3); 
sta(:,:) = nanmean(images,3);

% z-score
zscore_mean = nanmean(images_base,3);
zscore_std = nanstd(double(images_base),[],3);
sta_resp_zscored = (double(sta_resp) - zscore_mean) ./ zscore_std;

F = paper_figure([0,0.5,mm2inch(3*34),mm2inch(3*34)]); hold on;

sta_resp_zscored_woRegistrationArtefact = sta_resp_zscored;
for i=1:512
    for j=1:512
        if sta_resp_zscored(i,j)>3
            sta_resp_zscored_woRegistrationArtefact(i,j) = 0;
        end
    end
end

imagesc(sta_resp_zscored_woRegistrationArtefact',[0,1])
colormap('gray')
scatter(trg.laser_y(:,clusterID),trg.laser_x(:,clusterID),'r+','LineWidth',0.5,'SizeData',5)

daspect([1,1,1]);
xlim([1,s2p_meta.ops.Lx])
ylim([1,s2p_meta.ops.Ly])
box on;
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);

savefig(F,[save_root_fig,'\Fig4_ClusterSTA_STAImage.fig']);
saveas(F,[save_root_png,'\Fig4_ClusterSTA_STAImage.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ClusterSTA_STAImage.pdf']); set(gcf,'Color',[1,1,1])


%% Fig4_ClusterSTA_ResponseMap_Colorbar

imagesc(sta_resp_zscored_woRegistrationArtefact',[0,1])
colormap('gray');
temp=colorbar;
    
savefig(F,[save_root_fig,'\Fig4_ClusterSTA_STAImage_Colorbar.fig']);
saveas(F,[save_root_png,'\Fig4_ClusterSTA_STAImage_Colorbar.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ClusterSTA_STAImage_Colorbar.pdf']); set(gcf,'Color',[1,1,1])


%% --- Fig4_ClusterSTA_ExpressionImage ---

path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig4_fig\'];
save_root_png = [path.root_summary,'figures\Fig4_png\'];
save_root_pdf = [path.root_summary,'figures\Fig4_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig4_txt\'];

% load data
[path,raw] = data2repo_gatherRawDataInfo(info,path);
load([path.filepart_out,'s2p_meta.mat']);

F = paper_figure([0,0.5,mm2inch(3*34),mm2inch(3*34)]); hold on;

this_img = s2p_meta.ops.meanImg;
imagesc(this_img',[128,170])
colormap('gray')
scatter(trg.laser_y(:,clusterID),trg.laser_x(:,clusterID),'r+','LineWidth',0.5,'SizeData',5)

daspect([1,1,1]);
xlim([1,s2p_meta.ops.Lx])
ylim([1,s2p_meta.ops.Ly])
box on;
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);

savefig(F,[save_root_fig,'\Fig4_ClusterSTA_ExpressionImage.fig']);
saveas(F,[save_root_png,'\Fig4_ClusterSTA_ExpressionImage.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ClusterSTA_ExpressionImage.pdf']); set(gcf,'Color',[1,1,1])


%% --- Fig4_ClusterSTA_CorrelationImage ---

path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig4_fig\'];
save_root_png = [path.root_summary,'figures\Fig4_png\'];
save_root_pdf = [path.root_summary,'figures\Fig4_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig4_txt\'];

% load data
[path,raw] = data2repo_gatherRawDataInfo(info,path);
load([path.filepart_out,'s2p_meta.mat']);

F = paper_figure([0,0.5,mm2inch(3*34),mm2inch(3*34)]); hold on;

this_img = zeros(512,512);
this_img(s2p_meta.ops.yrange(1)+1:s2p_meta.ops.yrange(2),s2p_meta.ops.xrange(1)+1:s2p_meta.ops.xrange(2)) = s2p_meta.ops.Vcorr;
imagesc(this_img',[0,1])
colormap('gray')
scatter(trg.laser_y(:,clusterID),trg.laser_x(:,clusterID),'r+','LineWidth',0.5,'SizeData',5)

daspect([1,1,1]);
xlim([1,s2p_meta.ops.Lx])
ylim([1,s2p_meta.ops.Ly])
box on;
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);

savefig(F,[save_root_fig,'\Fig4_ClusterSTA_CorrelationImage.fig']);
saveas(F,[save_root_png,'\Fig4_ClusterSTA_CorrelationImage.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ClusterSTA_CorrelationImage.pdf']); set(gcf,'Color',[1,1,1])


%% Fig4_ClusterSTA_Traces

% import data using Analyses_Master with ops.do_responseAnalysis = true;
% use Philip_20211004 (seq) with cluster 34 (8 is nice as well)

load([path.filepart_out,'resp.mat']);

save_root_fig = [path.root_summary,'figures\Fig4_fig\'];
save_root_png = [path.root_summary,'figures\Fig4_png\'];
save_root_pdf = [path.root_summary,'figures\Fig4_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig4_txt\'];


%% Preparations

% overkill but okay
[prop,~,~] = preprocessActivityMeasure(act,p.genPrep,p,thor_beh,iscell);
trials = createTrialsStruct(task,1:prop.numTrials_incl);

% smooth act
%act_smoothed = smoothdata(act,2,'gaussian',3*5);
act_smoothed = act;
act_smoothed = nanzscore(act_smoothed,[],2);

% split unbinned data into trials: act -> actt
these_trial_bounds(1,:) = prop.sync_all-length(p.general.frames_pre);
these_trial_bounds(2,:) = these_trial_bounds(1,:)+p.general.numBins*p.general.binSize-1;
for i=1:length(prop.sync_all_binned)
    prop.trial_frames_unbinned(:,i) = these_trial_bounds(1,i) : these_trial_bounds(2,i);
end
actt = reshape(act_smoothed(:,prop.trial_frames_unbinned),[],size(prop.trial_frames_unbinned,1),length(prop.sync_all_binned));


%% Select data and make figure

stimType = 'seq';
trialType = 'X';
these_neurons = [182,83,91,242,770,1076]; % (Philip_20211004: seq, X, trg.idcs_targetedCells(:,34)==[83,770,91,1076,182,242]

% get data
sta = actt(these_neurons,:,intersect(find(trg.seq(2,:)),trials.stimuli.(trialType)));
sta_mean = nanmean(sta,3);

F = paper_figure([0,0.5,mm2inch(0.5*34),mm2inch(34)]); hold on;

distanceScaling = 6; ymax = 41;
this_t = p.general.t_unbinned;
patch([interp1(this_t,1:length(this_t),0),interp1(this_t,1:length(this_t),0),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1)],...
    [0,ymax,ymax,0],mean([p.col.odour;p.col.white]),'EdgeColor','none');
patch([interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+0.09),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+0.09)],...
    [0,ymax,ymax,0],mean([p.col.photostim;p.col.white]),'EdgeColor','none');

for i=1:length(these_neurons)
    shadedErrorBar(1:length(p.general.t_unbinned),nanmean(sta(i,:,:),3)+(length(these_neurons)+1-i)*distanceScaling,nanstd(sta(i,:,:),[],3),'lineProps',p.col.black);
end
xlim([60,150]) % 0: 91, 0.3: 100
ylim([4,ymax+1])
box off; axis off;

% save plot
% savefig(F,[save_root_fig,'\Fig4_ClusterSTA_Traces.fig']);
% saveas(F,[save_root_png,'\Fig4_ClusterSTA_Traces.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_ClusterSTA_Traces.pdf']); set(gcf,'Color',[1,1,1])


%%








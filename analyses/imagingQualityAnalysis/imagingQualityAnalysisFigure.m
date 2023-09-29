%function F = imagingQualityAnalysisFigure(info,p,iscell,s2p_meta,thor_beh)

%F = default_figure([20,0.5,20,9.9]);
default_figure([20,0.5,20,9.9]);

nrows = 2;
ncols = 4;


%% Noise level map

subplot(nrows,ncols,1)
hold on

img = nan(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
for n=1:length(iscell)
    if iscell(n)
        for i=1:length(s2p_meta.stat{n}.xpix)
            img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = iqa.noise(n); %s2p_meta.noise.dFF.all(n);
        end
    end
end

%img = flip(img,1); % THIS WAS MISSING IN RESPONSE ANALYSIS. IS NECESSARY, OTHERWISE AP axis is mirrored
meds = zeros(length(s2p_meta.stat),2);
for i = 1:length(s2p_meta.stat)
    meds(i,:) = [double(s2p_meta.stat{i}.med(2)),double(s2p_meta.stat{i}.med(1))];
end
%meds = flip(meds,2); % THIS WAS MISSING IN RESPONSE ANALYSIS. IS NECESSARY, OTHERWISE AP axis is mirrored

%plt.clim = [0,3];
% if any(img(:)<plt.clim(1)) | any(img(:)>plt.clim(2))
%     warning('Some data points are out of colour range')
% end
plt.clabel = 'Noise level';

imAlpha=ones(size(img)); imAlpha(isnan(img))=0;
imagesc(img,'AlphaData',imAlpha);
% set(gca,'CLim',[plt.clim(1),plt.clim(2)]);
colormap(jet);
temp=colorbar; 
temp.Label.String = plt.clabel;
daspect([1,1,1]);
hold on

xlim([0,s2p_meta.ops.Lx])
ylim([0,s2p_meta.ops.Ly])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','on');
% set(gca,'xcolor','none');
% set(gca,'ycolor','none');
hold off
title('Noise level map')


%% Noise level histogram

if any(img(:)<0) | any(img(:)>10)
    warning('Some data points are out of xlim range')
end

subplot(nrows,ncols,(2-1)*ncols+1)
plt = struct(); plt.xlabel = 'Noise level'; plt.xlim = [0,10]; plt.ylabel = 'Cells'; plt.norm = 'count'; plt.edges = [0:0.5:10];
myHistogram(s2p_meta.noise.dFF.all(find(iscell)),p,plt);


%% Lateral displacement

% offset_ft = absoff(prop.trial_frames+info.scope.numFrames_pre);
% 
% % binning
% temp = movmean(offset_ft,p.sca.binSize,1,'omitnan');
% offset_ft_binned = temp(floor(p.sca.binSize/2)+1:p.sca.binSize:end,:);
% offset_ft_binned = offset_ft_binned(1:end-1,:);

subplot(nrows,ncols,(2-1)*ncols+2)
plt = struct(); plt.xlabel = 'Lateral displacement'; plt.ylabel = 'Frames'; plt.norm = 'count'; %plt.xlim = [0,10];  %plt.edges = [0:0.5:10];
myHistogram(absoff,p,plt);

subplot(nrows,ncols,(1-1)*ncols+2)
plot(offset_ft_binned','Color',p.col.gray,'LineWidth',0.5)
temp=shadedErrorBar(1:size(offset_ft_binned,2),nanmean(offset_ft_binned,1),nansem(offset_ft_binned,1),'lineProps',p.col.black); temp.mainLine.LineWidth = 2;   
plt = struct(); traces_task(sca,p,info,plt); 



%% Lateral displacement - NEW















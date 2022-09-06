%function F = imagingQualityAnalysisFigure_maps(info,p,iscell,s2p_meta,thor_beh)

%F = default_figure([20,0.5,20,9.9]);
default_figure([20,0.5,20,9.9]);

nrows = 2;
ncols = 4;


%% Map - Mean

subplot(nrows,ncols,1)
hold on

img = nan(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
for n=1:length(iscell)
    if iscell(n)
        for i=1:length(s2p_meta.stat{n}.xpix)
            img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = iqa.mean(n); %s2p_meta.noise.dFF.all(n);
        end
    end
end
img = flip(img,1); % THIS WAS MISSING IN RESPONSE ANALYSIS. IS NECESSARY, OTHERWISE AP axis is mirrored

imAlpha=ones(size(img)); imAlpha(isnan(img))=0;
imagesc(img,'AlphaData',imAlpha);
% set(gca,'CLim',[plt.clim(1),plt.clim(2)]);
colormap(jet);
temp=colorbar; 
temp.Label.String = 'Mean (dFF)';
daspect([1,1,1]);
hold on

xlim([0,s2p_meta.ops.Lx])
ylim([0,s2p_meta.ops.Ly])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','on');
hold off
title('Imaging quality map - Mean')


%% Map - Std

subplot(nrows,ncols,2)
hold on

img = nan(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
for n=1:length(iscell)
    if iscell(n)
        for i=1:length(s2p_meta.stat{n}.xpix)
            img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = iqa.std(n); %s2p_meta.noise.dFF.all(n);
        end
    end
end
img = flip(img,1); % THIS WAS MISSING IN RESPONSE ANALYSIS. IS NECESSARY, OTHERWISE AP axis is mirrored

imAlpha=ones(size(img)); imAlpha(isnan(img))=0;
imagesc(img,'AlphaData',imAlpha);
% set(gca,'CLim',[plt.clim(1),plt.clim(2)]);
colormap(jet);
temp=colorbar; 
temp.Label.String = 'Std (dFF)';
daspect([1,1,1]);
hold on

xlim([0,s2p_meta.ops.Lx])
ylim([0,s2p_meta.ops.Ly])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','on');
hold off
title('Imaging quality map - Std')


%% Map - SNR

subplot(nrows,ncols,3)
hold on

img = nan(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
for n=1:length(iscell)
    if iscell(n)
        for i=1:length(s2p_meta.stat{n}.xpix)
            img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = iqa.snr(n); %s2p_meta.noise.dFF.all(n);
        end
    end
end
img = flip(img,1); % THIS WAS MISSING IN RESPONSE ANALYSIS. IS NECESSARY, OTHERWISE AP axis is mirrored

imAlpha=ones(size(img)); imAlpha(isnan(img))=0;
imagesc(img,'AlphaData',imAlpha);
% set(gca,'CLim',[plt.clim(1),plt.clim(2)]);
colormap(jet);
temp=colorbar; 
temp.Label.String = 'SNR (dFF)';
daspect([1,1,1]);
hold on

xlim([0,s2p_meta.ops.Lx])
ylim([0,s2p_meta.ops.Ly])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','on');
hold off
title('Imaging quality map - SNR')


%% Map - Noise

subplot(nrows,ncols,4)
hold on

img = nan(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
for n=1:length(iscell)
    if iscell(n)
        for i=1:length(s2p_meta.stat{n}.xpix)
            img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = iqa.noise(n); %s2p_meta.noise.dFF.all(n);
        end
    end
end
img = flip(img,1); % THIS WAS MISSING IN RESPONSE ANALYSIS. IS NECESSARY, OTHERWISE AP axis is mirrored

imAlpha=ones(size(img)); imAlpha(isnan(img))=0;
imagesc(img,'AlphaData',imAlpha);
% set(gca,'CLim',[plt.clim(1),plt.clim(2)]);
colormap(jet);
temp=colorbar; 
temp.Label.String = 'Noise level (z-scored dFF)';
daspect([1,1,1]);
hold on

xlim([0,s2p_meta.ops.Lx])
ylim([0,s2p_meta.ops.Ly])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','on');
hold off
title('Imaging quality map - Noise level')


%% Map - Skewness

subplot(nrows,ncols,5)
hold on

img = nan(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
for n=1:length(iscell)
    if iscell(n)
        for i=1:length(s2p_meta.stat{n}.xpix)
            img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = iqa.skew(n); %s2p_meta.noise.dFF.all(n);
        end
    end
end
img = flip(img,1); % THIS WAS MISSING IN RESPONSE ANALYSIS. IS NECESSARY, OTHERWISE AP axis is mirrored

imAlpha=ones(size(img)); imAlpha(isnan(img))=0;
imagesc(img,'AlphaData',imAlpha);
% set(gca,'CLim',[plt.clim(1),plt.clim(2)]);
colormap(jet);
temp=colorbar; 
temp.Label.String = 'Skewness (dFF)';
daspect([1,1,1]);
hold on

xlim([0,s2p_meta.ops.Lx])
ylim([0,s2p_meta.ops.Ly])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','on');
hold off
title('Imaging quality map - Skewness')



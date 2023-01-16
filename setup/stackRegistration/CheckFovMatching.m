%% Check FOV Matching (and/or Targeting)

clear; clc; 

path.refFrameFile_previous  = 'F:\Data\2020\2020-11\2020-11-19\Changa\Expression\Changa_20201119_930A_before_003_REF\ChanA_001_001_001_001.tif';
path.refFrameFile_today     = 'F:\Data\2020\2020-11\2020-11-20\Changa\Expression\Changa_20201120_930A_before_003_REF\ChanA_001_001_001_001.tif'; % 'C:\Users\user\Documents\MATLAB\SniffinHippo\StackRegistration\exampleData\REF.tif';

p.plotTargets               = true;
path.targetsPath            = 'F:\Data\2020\2020-11\2020-11-19\Changa\Targeting\TargetSelection\';

p.meanImg_upperDiv          = 3;
p.meanImg_inset             = 150:512-150;
p.meanImg_50umLines         = true;


%% Load data

% Load reference frame
refFrame_previous = imread(path.refFrameFile_previous);
refFrame_today = imread(path.refFrameFile_today);

% Load targets
if p.plotTargets
    temp = dir([path.targetsPath,'Targeting*.mat']);
    load([path.targetsPath,temp.name]);
end


%% Figure

meanImg_upperDiv = 2%p.meanImg_upperDiv;

meanImg_lower = min(refFrame_previous(:));
meanImg_upper = max(refFrame_previous(:))/meanImg_upperDiv;

if p.plotTargets
    subPlotRows = 2;
else
    subPlotRows = 1;
end

G = figure%default_figure([0,3.5,20,5]);

Gp1 = subplot(subPlotRows,3,1);
imshow(refFrame_previous,[meanImg_lower,meanImg_upper]);
colormap('gray');
xlim([p.meanImg_inset(1),p.meanImg_inset(end)]);
ylim([p.meanImg_inset(1),p.meanImg_inset(end)]);
if p.meanImg_50umLines
    hold on
    for i=0:20
        xline(i*50*512/1000.78,'r');
        yline(i*50*512/1000.78,'r');
    end
    hold off
end
title('Previous REF')
xlabel('medial (X, +) <-> lateral (X, -)')  
ylabel('posterior (Y, +) <-> anterior (Y, -)')

Gp2 = subplot(subPlotRows,3,2);
imshow(refFrame_today,[meanImg_lower,meanImg_upper]);
colormap('gray');
xlim([p.meanImg_inset(1),p.meanImg_inset(end)]);
ylim([p.meanImg_inset(1),p.meanImg_inset(end)]);
if p.meanImg_50umLines
    hold on
    for i=0:20
        xline(i*50*512/1000.78,'r');
        yline(i*50*512/1000.78,'r');
    end
    hold off
end
title('REF of today')
xlabel('medial (X, +) <-> lateral (X, -)')  
ylabel('posterior (Y, +) <-> anterior (Y, -)')

Gp3 = subplot(subPlotRows,3,3);
imshowpair(refFrame_previous,refFrame_today,'falsecolor','ColorChannels','green-magenta','Scaling','independent');
xlim([p.meanImg_inset(1),p.meanImg_inset(end)]);
ylim([p.meanImg_inset(1),p.meanImg_inset(end)]);
if p.meanImg_50umLines
    hold on
    for i=0:20
        xline(i*50*512/1000.78,'r');
        yline(i*50*512/1000.78,'r');
    end
    hold off
end
title('Comparison (green: previous, magenta: today)')
xlabel('medial (X, +) <-> lateral (X, -)')  
ylabel('posterior (Y, +) <-> anterior (Y, -)')

if p.plotTargets
    Gp4 = subplot(subPlotRows,3,4);
    imshow(trg.s2p.ops.meanImg,[meanImg_lower,meanImg_upper]);
    colormap('gray');
    xlim([p.meanImg_inset(1),p.meanImg_inset(end)]);
    ylim([p.meanImg_inset(1),p.meanImg_inset(end)]);
    if p.meanImg_50umLines
        hold on
        for i=0:20
            xline(i*50*512/1000.78,'r');
            yline(i*50*512/1000.78,'r');
        end
        hold off
    end
    title('suite2p mean image for targeting')
    xlabel('medial (X, +) <-> lateral (X, -)')  
    ylabel('posterior (Y, +) <-> anterior (Y, -)')

    Gp5 = subplot(subPlotRows,3,5);
    imshow(refFrame_today,[meanImg_lower,meanImg_upper]);
    colormap('gray');
    xlim([p.meanImg_inset(1),p.meanImg_inset(end)]);
    ylim([p.meanImg_inset(1),p.meanImg_inset(end)]);
    if p.meanImg_50umLines
        hold on
        for i=0:20
            xline(i*50*512/1000.78,'r');
            yline(i*50*512/1000.78,'r');
        end
        hold off
    end
    title('REF of today')
    xlabel('medial (X, +) <-> lateral (X, -)')  
    ylabel('posterior (Y, +) <-> anterior (Y, -)')

    Gp6 = subplot(subPlotRows,3,6);
    imshowpair(trg.s2p.ops.meanImg,refFrame_today,'falsecolor','ColorChannels','green-magenta','Scaling','independent');
    xlim([p.meanImg_inset(1),p.meanImg_inset(end)]);
    ylim([p.meanImg_inset(1),p.meanImg_inset(end)]);
    if p.meanImg_50umLines
        hold on
        for i=0:20
            xline(i*50*512/1000.78,'r');
            yline(i*50*512/1000.78,'r');
        end
        hold off
    end
    title('Comparison (green: suite2p, magenta: today)')
    xlabel('medial (X, +) <-> lateral (X, -)')  
    ylabel('posterior (Y, +) <-> anterior (Y, -)')
end

if p.plotTargets
    for k=1:6
        subplot(subPlotRows,3,k)
        hold on
        for i=1:trg.p.numClustersPerGroup*2  
            for j = 1:trg.p.numCellsPerCluster
                scatter(trg.f.med2(trg.clustering(j,i)),trg.f.med1(trg.clustering(j,i)),'y.');
            end
        end
        hold off
    end
end


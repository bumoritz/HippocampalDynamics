%% Semi-automated online z-correction

clear; clc; 
restoredefaultpath; addpath(genpath('C:\Users\user\Documents\MATLAB\SniffinHippo'));

path.stackPath      = ['F:\Data\2020\2020-11\2020-11-19\Changa\Expression\Changa-20201119-930A-stack1\';...
    'F:\Data\2020\2020-11\2020-11-19\Changa\Expression\Changa-20201119-930A-stack2\';...
    'F:\Data\2020\2020-11\2020-11-19\Changa\Expression\Changa-20201119-930A-stack3\'];

path.refFrameFile   = 'F:\Data\2020\2020-11\2020-11-19\Changa\Expression\Changa_20201119_930A_before_003_REF\ChanA_001_001_001_001.tif'; % 'C:\Users\user\Documents\MATLAB\SniffinHippo\StackRegistration\exampleData\REF.tif';
path.targetsPath    = 'F:\Data\2020\2020-11\2020-11-19\Changa\Targeting\TargetSelection\';

path.imagingPath    = 'F:\Data\2020\2020-11\2020-11-19\Changa\Targeting\Changa-20201119-spont\'; % incl. backslash
path.savePath       = 'F:\Data\2020\2020-11\2020-11-19\Changa\Targeting\dev\';

p.stackTop          = 25; % [um]
p.stackBottom       = -25 ; % [um]
p.stackStepSize     = 1; % [um]
p.frameRate         = 30.0; % [Hz]

p.repeatInterval    = 15; % [s]
p.avgNumFrames      = 90;
p.poolingInterval   = 5; % [s]
p.zMethod           = 'xcorr'; % 'max' or 'xcorr'
p.tolerance         = 2; % [um]

p.do_meanImg            = true;
p.meanImg_upperDiv      = 3;
p.meanImg_inset         = 150:512-150;%1:512;
p.meanImg_50umLines     = true;
p.plotTargets           = true;

% path.stackFile      = ['F:\2021\2021-02\2021-02-12\Biontech\Expression\stack1.tif';...
%     'F:\2021\2021-02\2021-02-12\Biontech\Expression\stack2.tif';...
%     'F:\2021\2021-02\2021-02-12\Biontech\Expression\stack3.tif'];    


%% Prepare execution

diary([path.savePath,'Zcorr.log']);
disp(['Semi-automated online z-correction.'])
disp(['- Loading data and preparing execution.'])

% Load reference frame and set up registration process
refFrame = imread(path.refFrameFile);
ops = setup_registration_phasecorr(refFrame);
p.bytesPerFrame = size(refFrame,1)*size(refFrame,2)*2; % 2 bytes per value (precision) for uint16

% Load targets
if p.plotTargets
    temp = dir([path.targetsPath,'Targeting*.mat']);
    load([path.targetsPath,temp.name]);
end

% Load stack
p.numSlices = length(p.stackBottom:p.stackStepSize:p.stackTop);
% if p.numSlices ~= length(imfinfo(path.stackFile(1,:)))
%     error('Number of slices in stack not consistent between file and input.')
% end
stack = zeros(size(refFrame,1),size(refFrame,2),p.numSlices,'uint16');
for i=1:p.numSlices
    temp = zeros(size(refFrame,1),size(refFrame,2),size(path.stackPath,1),'uint16');
    for j=1:size(path.stackPath,1)        
        temp2 = dir([path.stackPath(j,:),'Chan*']);
        this_slice_path = [path.stackPath(j,:),temp2(i).name];
        temp(:,:,j) = imread(this_slice_path); %imread(path.stackFile(j,:),i);
    end
    stack(:,:,i) = mean(temp,3);
    if i==1
        imwrite(stack(:,:,i),[path.savePath,'avgStack.tif'],'Compression','none');
    else
        imwrite(stack(:,:,i),[path.savePath,'avgStack.tif'],'WriteMode','append','Compression','none');
    end
end

% Register reference against stack
out0.shift = zeros(p.numSlices,2);
out0.corr = zeros(p.numSlices,1);
for i=1:p.numSlices
    [this_regFrame, out0.shift(i,:), ~] = return_offsets_phasecorr(single(stack(:,:,i)), ops);
    out0.corr(i) = corr2(this_regFrame,refFrame);
end
out0.corr = smoothdata(out0.corr,'gaussian',5);
[~,out0.best_slice] = max(out0.corr);
out0.best_um = p.stackTop - (out0.best_slice-1)*p.stackStepSize;

% Initialise output structs
out = out0;
out.minute = NaN;
out.best_slice = NaN;
out.best_um = NaN;
out.diff_slice = NaN;
out.diff_um = NaN;
out.best_slice_max = NaN;
out.best_um_max = NaN;
out.diff_slice_max = NaN;
out.diff_um_max = NaN;
out.best_slice_xcorr = NaN;
out.best_um_xcorr = NaN;
out.diff_slice_xcorr = NaN;
out.diff_um_xcorr = NaN;
out.driftDirection = 'none';
out.compensationDirection = 'none';
out.moveCommand = false;
out.plotCol = 'k';
outPrev = out;

%% Initialise figure

disp(['- Initialising figure.'])

F = figure%default_figure([0,3.5,20,5]);

Fp1 = subplot(1,3,1);
yyaxis left; set(gca,'YColor','k');
plot(out0.corr,1:p.numSlices,'Color','k');
xlim([-0.1,1])
xlabel('Correlation: REFERENCE')
hold on
line(get(gca,'XLim'),[out0.best_slice,out0.best_slice],'Color','k')
hold off
set(gca,'ydir','reverse');
ylim([1,p.numSlices])
yticks(0:10:p.numSlices+1)
ylabel('Slice number')
yyaxis right; set(gca,'YColor','k');
temp = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(temp(1),temp(2),['MAX CORR',newline,'slice = ',num2str(out0.best_slice),newline,'z = ',num2str(out0.best_um),' \mum',newline,'vert. = ',num2str(out0.shift(out0.best_slice,1)),' \mum',newline,'hor. = ',num2str(out0.shift(out0.best_slice,2)),' \mum'],'FontSize',10,'VerticalAlignment','top','HorizontalAlignment','right');
ylim([p.stackBottom,p.stackTop])
yticks(p.stackBottom:10:p.stackTop)
ylabel('z-Coordinate (\mum)')

Fp2 = subplot(1,3,2);
yyaxis left; set(gca,'YColor','k');
plot(out0.corr,1:p.numSlices,'Color',[0.7,0.7,0.7]);
xlim([-0.1,1])
xlabel('Correlation: PREVIOUS')
hold on
line(get(gca,'XLim'),[out0.best_slice,out0.best_slice],'Color',[0.7,0.7,0.7],'LineStyle',':')
hold off
set(gca,'ydir','reverse');
ylim([1,p.numSlices])
yticks(0:10:p.numSlices+1)
ylabel('Slice number')
yyaxis right; set(gca,'YColor','k');
temp = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(temp(1),temp(2),['MAX CORR',newline,'slice = ',num2str(outPrev.best_slice_max),newline,'z = ',num2str(outPrev.best_um_max),' \mum',newline,'vert. = ',num2str(NaN),' \mum',newline,'hor. = ',num2str(NaN),' \mum',newline,'drift = ',num2str(outPrev.diff_um_max),' \mum',newline,newline,'XCORR LAG',newline,'vert. = ',num2str(NaN),' \mum',newline,'hor. = ',num2str(NaN),' \mum',newline,'drift = ',num2str(outPrev.diff_um_xcorr),' \mum'],'FontSize',10,'VerticalAlignment','top','HorizontalAlignment','right');
                
temp = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(temp(1),temp(2),['time = ',num2str(outPrev.minute,4),' min'],'FontSize',10,'VerticalAlignment','bottom','HorizontalAlignment','right');
ylim([p.stackBottom,p.stackTop])
yticks(p.stackBottom:10:p.stackTop)
ylabel('z-Coordinate (\mum)')

Fp3 = subplot(1,3,3);
yyaxis left; set(gca,'YColor','k');
plot(out0.corr,1:p.numSlices,'Color',[0.7,0.7,0.7]);
xlim([-0.1,1])
xlabel('Correlation: CURRENT')
hold on
line(get(gca,'XLim'),[out0.best_slice,out0.best_slice],'Color',[0.7,0.7,0.7],'LineStyle',':')
hold off
set(gca,'ydir','reverse');
ylim([1,p.numSlices])
yticks(0:10:p.numSlices+1)
ylabel('Slice number')
yyaxis right; set(gca,'YColor','k');
temp = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(temp(1),temp(2),['MAX CORR',newline,'slice = ',num2str(out.best_slice_max),newline,'z = ',num2str(out.best_um_max),' \mum',newline,'vert. = ',num2str(NaN),' \mum',newline,'hor. = ',num2str(NaN),' \mum',newline,'drift = ',num2str(out.diff_um_max),' \mum',newline,newline,'XCORR LAG',newline,'vert. = ',num2str(NaN),' \mum',newline,'hor. = ',num2str(NaN),' \mum',newline,'drift = ',num2str(out.diff_um_xcorr),' \mum'],'FontSize',10,'VerticalAlignment','top','HorizontalAlignment','right');
temp = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(temp(1),temp(2),['time = ',num2str(out.minute,4),' min'],'FontSize',10,'VerticalAlignment','bottom','HorizontalAlignment','right');
ylim([p.stackBottom,p.stackTop])
yticks(p.stackBottom:10:p.stackTop)
ylabel('z-Coordinate (\mum)')
drawnow;


%% Initialise image figures

if p.do_meanImg
    meanImg_lower = min(refFrame(:));
    meanImg_upper = max(refFrame(:))/p.meanImg_upperDiv;
    
    G = figure%default_figure([0,3.5,20,5]);

    Gp1 = subplot(1,3,1);
    imshow(refFrame,[meanImg_lower,meanImg_upper]);
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
    title('Reference mean image')
    xlabel('medial (X, +) <-> lateral (X, -)')  
    ylabel('posterior (Y, +) <-> anterior (Y, -)')

    Gp2 = subplot(1,3,2);
    imshow(refFrame,[meanImg_lower,meanImg_upper]);
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
    title('Current mean image')
    xlabel('medial (X, +) <-> lateral (X, -)')  
    ylabel('posterior (Y, +) <-> anterior (Y, -)')

    Gp3 = subplot(1,3,3);
    imshowpair(refFrame,refFrame,'falsecolor','ColorChannels','green-magenta','Scaling','independent');
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
    title('Comparison (green: REF, magenta: CURRENT)')
    xlabel('medial (X, +) <-> lateral (X, -)')  
    ylabel('posterior (Y, +) <-> anterior (Y, -)')
    
    if p.plotTargets
        for k=1:3
            subplot(1,3,k)
            hold on
            for i=1:trg.p.numClustersPerGroup*2  
                for j = 1:trg.p.numCellsPerCluster
                    scatter(trg.f.med2(trg.clustering(j,i)),trg.f.med1(trg.clustering(j,i)),'y.');
                end
            end
            hold off
        end
    end
end


%% Infinite loop

disp(['- Preparations completed. Start execution by pressing any key!'])
pause;
disp(['- Online registration started.'])

n=0;
t0=tic;
execution_ongoing = true;
while execution_ongoing
    t=toc(t0);
    if t > (n+1)*p.repeatInterval
        tic;
        n=n+1;
        outPrev = out;
        
        % Find and open imagaing file
        temp =  dir([path.imagingPath,'Image_*']);
        if isempty(temp)
            temp = dir([path.imagingPath,'tis*']);
        end
        path.imagingFile = [path.imagingPath,temp.name];
        this_fid = fopen(path.imagingFile, 'r');
        fseek(this_fid,0,'bof');
        
        % Load frames and close file
        this_pool =  ceil((n*p.repeatInterval-p.poolingInterval)*p.frameRate)  : floor(n*p.repeatInterval*p.frameRate);
        if max(this_pool) > temp.bytes/p.bytesPerFrame
            execution_ongoing = false;
        else
            these_frame_numbers = randsample(this_pool,p.avgNumFrames);
            these_frames = zeros(size(refFrame,1),size(refFrame,2),p.avgNumFrames,'uint16');
            for i=1:p.avgNumFrames
                fseek(this_fid,(these_frame_numbers(i)-1)*p.bytesPerFrame,'bof');
                temp = uint16(fread(this_fid,size(refFrame,1)*size(refFrame,2),'uint16',0,'l'));
                these_frames(:,:,i) = reshape(temp,size(refFrame,1),size(refFrame,2));
                frewind(this_fid);
            end
            fclose(this_fid);

            % Calculate average image and set up registration process
            this_avgImg = nanmean(these_frames,3);
            this_avgImg = this_avgImg';
            this_ops = setup_registration_phasecorr(this_avgImg);

            % Register average image against stack
            for i=1:p.numSlices
                [this_regFrame, out.shift(i,:), ~] = return_offsets_phasecorr(single(stack(:,:,i)), this_ops);
                out.corr(i) = corr2(this_regFrame,this_avgImg);
            end
            out.corr = smoothdata(out.corr,'gaussian',5);
            [~,out.best_slice_max] = max(out.corr);
            out.best_um_max = p.stackTop - (out.best_slice_max-1)*p.stackStepSize;
            out.minute = t/60;
            [temp,temp1] = xcorr(out0.corr,out.corr);
            [~,temp2] = max(temp);
            temp1(temp2);
            out.best_slice_xcorr = out0.best_slice-temp1(temp2);
            out.best_um_xcorr = p.stackTop - (out.best_slice_xcorr-1)*p.stackStepSize;
                       
            % Assess outcome
            out.diff_slice_max = out.best_slice_max - out0.best_slice;
            out.diff_um_max = out.best_um_max - out0.best_um;
            out.diff_slice_xcorr = out.best_slice_xcorr - out0.best_slice;
            out.diff_um_xcorr = out.best_um_xcorr - out0.best_um;
            if strcmp(p.zMethod,'max')
                out.best_slice = out.best_slice_max;
                out.best_um = out.best_um_max;
                out.diff_slice = out.diff_slice_max;
                out.diff_um = out.diff_um_max;
            elseif strcmp(p.zMethod,'xcorr')
                out.best_slice = out.best_slice_xcorr;
                out.best_um = out.best_um_xcorr;
                out.diff_slice = out.diff_slice_xcorr;
                out.diff_um = out.diff_um_xcorr;
            end
            if out.diff_um > 0
                out.driftDirection = 'upwards';
                out.compensationDirection = 'downwards';
            elseif out.diff_um < 0
                out.driftDirection = 'downwards';
                out.compensationDirection = 'upwards';
            else
                out.driftDirection = 'none';
                out.compensationDirection = 'none';
            end
            if abs(out.diff_um) > p.tolerance
                out.plotCol = 'r';
            elseif abs(out.diff_um)==0
                out.plotCol = 'g';
            else
                out.plotCol = 'y';
            end
            if (abs(out.diff_um) > p.tolerance) && (abs(outPrev.diff_um) > p.tolerance)
                out.moveCommand = true;
            else
                out.moveCommand = false;
            end

            % Update plots: previous
            if n~=1
                delete(Fp2);
                Fp2 = subplot(1,3,2);
                yyaxis left; set(gca,'YColor','k');
                plot(out0.corr,1:p.numSlices,'Color',[0.7,0.7,0.7]);
                xlim([-0.1,1])
                xlabel('Correlation: PREVIOUS')
                hold on
                line(get(gca,'XLim'),[out0.best_slice,out0.best_slice],'Color',[0.7,0.7,0.7],'LineStyle',':')
                plot(outPrev.corr,1:p.numSlices,'-');
                line(get(gca,'XLim'),[outPrev.best_slice,outPrev.best_slice],'Color',outPrev.plotCol)
                hold off
                set(gca,'ydir','reverse');
                ylim([1,p.numSlices])
                yticks(0:10:p.numSlices+1)
                ylabel('Slice number')
                yyaxis right; set(gca,'YColor','k');
                temp = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
                text(temp(1),temp(2),['MAX CORR',newline,'slice = ',num2str(outPrev.best_slice_max),newline,'z = ',num2str(outPrev.best_um_max),' \mum',newline,'vert. = ',num2str(outPrev.shift(outPrev.best_slice_max,1)),' \mum',newline,'hor. = ',num2str(outPrev.shift(outPrev.best_slice_max,2)),' \mum',newline,'drift = ',num2str(outPrev.diff_um_max),' \mum',newline,newline,'XCORR LAG',newline,'vert. = ',num2str(outPrev.shift(outPrev.best_slice_xcorr,1)),' \mum',newline,'hor. = ',num2str(outPrev.shift(outPrev.best_slice_xcorr,2)),' \mum',newline,'drift = ',num2str(outPrev.diff_um_xcorr),' \mum'],'FontSize',10,'VerticalAlignment','top','HorizontalAlignment','right');
                temp = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
                text(temp(1),temp(2),['time = ',num2str(outPrev.minute,4),' min'],'FontSize',10,'VerticalAlignment','bottom','HorizontalAlignment','right');
                ylim([p.stackBottom,p.stackTop])
                yticks(p.stackBottom:10:p.stackTop)
                ylabel('z-Coordinate (\mum)')
            end

            % Update plots: current
            delete(Fp3);
            Fp3 = subplot(1,3,3);
            yyaxis left; set(gca,'YColor','k');
            plot(out0.corr,1:p.numSlices,'Color',[0.7,0.7,0.7]);
            xlim([-0.1,1])
            xlabel('Correlation: CURRENT')
            hold on
            line(get(gca,'XLim'),[out0.best_slice,out0.best_slice],'Color',[0.7,0.7,0.7],'LineStyle',':')
            plot(out.corr,1:p.numSlices,'-');
            line(get(gca,'XLim'),[out.best_slice,out.best_slice],'Color',out.plotCol)
            hold off
            set(gca,'ydir','reverse');
            ylim([1,p.numSlices])
            yticks(0:10:p.numSlices+1)
            ylabel('Slice number')
            yyaxis right; set(gca,'YColor','k');
            temp = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
            text(temp(1),temp(2),['MAX CORR',newline,'slice = ',num2str(out.best_slice_max),newline,'z = ',num2str(out.best_um_max),' \mum',newline,'vert. = ',num2str(out.shift(out.best_slice_max,1)),' \mum',newline,'hor. = ',num2str(out.shift(out.best_slice_max,2)),' \mum',newline,'drift = ',num2str(out.diff_um_max),' \mum',newline,newline,'XCORR LAG',newline,'vert. = ',num2str(out.shift(out.best_slice_xcorr,1)),' \mum',newline,'hor. = ',num2str(out.shift(out.best_slice_xcorr,2)),' \mum',newline,'drift = ',num2str(out.diff_um_xcorr),' \mum'],'FontSize',10,'VerticalAlignment','top','HorizontalAlignment','right');
            temp = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
            text(temp(1),temp(2),['time = ',num2str(out.minute,4),' min'],'FontSize',10,'VerticalAlignment','bottom','HorizontalAlignment','right');
            ylim([p.stackBottom,p.stackTop])
            yticks(p.stackBottom:10:p.stackTop)
            ylabel('z-Coordinate (\mum)')
            drawnow;
            
            % Update images
            if p.do_meanImg
                
                Gp2 = subplot(1,3,2);
                imshow(this_avgImg,[meanImg_lower,meanImg_upper]);
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
                title('Current mean image')
                xlabel('medial (X, +) <-> lateral (X, -)')  
                ylabel('posterior (Y, +) <-> anterior (Y, -)')
                
                Gp3 = subplot(1,3,3);
                imshowpair(refFrame,this_avgImg,'falsecolor','ColorChannels','green-magenta','Scaling','independent');
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
                title('Comparison (green: REF, magenta: CURRENT)')
                xlabel('medial (X, +) <-> lateral (X, -)')  
                ylabel('posterior (Y, +) <-> anterior (Y, -)')
            end

            % Print output to command window
            temp = toc;
            if ~strcmp(out.driftDirection,'none')
                disp(['--- Cycle ',num2str(n),' @ ',num2str(out.minute,2),' min (executed within ',num2str(temp,1),' s): Focal plane drifted cum. ',num2str(abs(out.diff_um),2),' um ',out.driftDirection,'.'])
                if out.moveCommand
                    disp(['%%% MOVE ',num2str(abs(out.diff_um)/2,2),' UM ',upper(out.compensationDirection),'! %%%'])
                end              
            else
                disp(['--- Cycle ',num2str(n),' @ ',num2str(out.minute,2),' min (executed within ',num2str(temp,1),' s): Focal plane remained stable.'])
            end
        end
    end
end

% Write output to log file
disp(['- Acquisition finished. Writing log file.']) 
path
p
disp(['- Done.']) 
diary off


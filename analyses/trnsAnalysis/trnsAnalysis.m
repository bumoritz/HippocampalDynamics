function [trns] = trnsAnalysis(path,data_in)
% data_in = act;

%% Preparations

% standardise
this_median = nanmedian(data_in,2);
this_std = nanstd(data_in,[],2);
dFFz = (data_in - this_median) ./ this_std;


%% Core

tic;
dFFns_beh = nan(size(data_in));
parfor i=1:length(this_median)
    if ~isnan(this_median(i))
        this_data_in = dFFz(i,:);

        % get thresholds and adjacent frame criteria
        thresholds = 1:0.2:4;
        numTresholdLevels = length(thresholds);
        maxAdjacentFrames = 10;
        numFrames = length(this_data_in);

        % calculate false positive rates
        putPos = nan(numTresholdLevels,numFrames);
        putNeg = nan(numTresholdLevels,numFrames);
        fpr = nan(numTresholdLevels,maxAdjacentFrames);
        positiveTransients = cell(numTresholdLevels,maxAdjacentFrames);
        negativeTransients = cell(numTresholdLevels,maxAdjacentFrames);
        for t=1:numTresholdLevels
            putPos(t,:) = this_data_in > thresholds(t);
            putNeg(t,:) = this_data_in < -thresholds(t);
            for n=1:maxAdjacentFrames
                if n==maxAdjacentFrames
                    positiveTransients{t,n} = pattern(putPos(t,:),[0,ones(1,n)])+1;
                    negativeTransients{t,n} = pattern(putNeg(t,:),[0,ones(1,n)])+1;
                else
                    positiveTransients{t,n} = pattern(putPos(t,:),[0,ones(1,n),0])+1;
                    negativeTransients{t,n} = pattern(putNeg(t,:),[0,ones(1,n),0])+1;
                end
                fpr(t,n) = length(negativeTransients{t,n}) / length(positiveTransients{t,n});
            end
        end

        % get threshold level
        thresholdLevel = nan(1,maxAdjacentFrames);
        for n=1:maxAdjacentFrames
            try
                thresholdLevel(n) = nanmin(find(fpr(:,n) < 0.001));
            catch
                thresholdLevel(n) = 1;
            end
        end

        % create output mask
        outputMask = zeros(1,numFrames);
        for n=1:maxAdjacentFrames
            if n<maxAdjacentFrames
                temp = positiveTransients{thresholdLevel(n),n}+([1:n]'-1);
                outputMask(temp(:)) = 1;
            else
                temp = positiveTransients{thresholdLevel(n),n};
                temp = temp(:);
                for k=1:length(temp)
                    outputMask(temp(k):temp(k)+nanmin(find(putPos(thresholdLevel(n),temp(k)+1:end)==0))-1) = 1;
                end
            end
        end

        % finalising
        outputMask(pattern(outputMask,[1,0,1])+1)=1;
        outputMask(pattern(outputMask,[0,1,0])+1)=0;
%         outputMask([pattern(outputMask,[1,0,0,1])+1,pattern(outputMask,[1,0,0,1])+2])=1;
%         outputMask([pattern(outputMask,[0,1,1,0])+1,pattern(outputMask,[0,1,1,0])+2])=0;
%         outputMask([pattern(outputMask,[1,0,0,0,1])+1,pattern(outputMask,[1,0,0,0,1])+2,pattern(outputMask,[1,0,0,0,1])+3])=1;
%         outputMask([pattern(outputMask,[0,1,1,1,0])+1,pattern(outputMask,[0,1,1,1,0])+2,pattern(outputMask,[0,1,1,1,0])+3])=0;
        this_data_out = this_data_in;
        this_data_out(outputMask==0) = 0;
        dFFns_beh(i,:) = this_data_out;
    end
end
toc;

save([path.filepart_in,'dFFns_beh.mat'],'dFFns_beh','-v7.3');
disp(['--- Saved dFFns_beh file as ',[path.filepart_in,'dFFns_beh.mat'],'.'])

end
%% Figures
% figure;
% plot(fpr')
% xlabel('Transient with at least n frames')
% ylabel('False positive rate')


% wdw = 10000:10000+240*30; %1:60*30;
% 
% figure; hold on;
% 
% subplot(2,1,1)
% plot(dFFz(i,wdw));
% title('dFF')
% 
% subplot(2,1,2)
% plot(dFFz_out(i,wdw));
% title('sig dFF out')



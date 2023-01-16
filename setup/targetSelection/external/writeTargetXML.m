function writeTargetXML(trg, savePath)
% required functions: colornames
% required trg fields: p, xml, f, clustering

% check X and Y
% check this weird pixel conversion
% check if galvo stuff works correctly


%% ThorImageSLM

% create document node and root element
docNode = com.mathworks.xml.XMLUtils.createDocument('ThorImageSLM');

% identify root element 
toc = docNode.getDocumentElement;


%% SLM Patterns, Pattern, ROI

% SLM Patterns
this_element = docNode.createElement('SLMPatterns');
toc.appendChild(this_element);

% Pattern
[~,all_cols] = colornames('CSS');
available_cols = 1:length(all_cols);
for i=1:trg.p.numClustersPerGroup*2    
    this_pattern = docNode.createElement('Pattern');
       
    % Pattern info
    temp = randsample(available_cols,1);
    this_pattern.setAttribute('name',[trg.xml.name_root,num2str(i,'%03.f')]);
    this_pattern.setAttribute('patternID',num2str(i));
    this_pattern.setAttribute('shape',trg.xml.shape);
    this_pattern.setAttribute('roiWidthPx',trg.xml.roiWidthPx);
    this_pattern.setAttribute('roiHeightPx',trg.xml.roiHeightPx);
    this_pattern.setAttribute('red',num2str(round(all_cols(temp,1)*255)));
    this_pattern.setAttribute('green',num2str(round(all_cols(temp,2)*255)));
    this_pattern.setAttribute('blue',num2str(round(all_cols(temp,3)*255)));
    this_pattern.setAttribute('pxSpacing',trg.xml.pxSpacing);
    this_pattern.setAttribute('durationMS',trg.xml.durationMS);
    this_pattern.setAttribute('iterations',trg.xml.iterations);
    this_pattern.setAttribute('power',trg.xml.power);
    this_pattern.setAttribute('prePatIdleMS',trg.xml.prePatIdleMS);
    this_pattern.setAttribute('postPatIdleMS',trg.xml.postPatIdleMS);
    this_pattern.setAttribute('preIteIdleMS',trg.xml.preIteIdleMS);
    this_pattern.setAttribute('postIteIdleMS',trg.xml.postIteIdleMS);
    this_pattern.setAttribute('measurePowerMW',trg.xml.measurePowerMW);
    this_pattern.setAttribute('measurePowerMWPerUM2',trg.xml.measurePowerMWPerUM2);
    this_element.appendChild(this_pattern);
    available_cols = setdiff(available_cols,temp);
    
    % ROI (galvo)
    this_roi = docNode.createElement('ROI');
    this_roi.setAttribute('subID',num2str(0));
    this_roi.setAttribute('centerX',num2str(trg.galvo(2,i)));
    this_roi.setAttribute('centerY',num2str(trg.galvo(1,i)));
    this_pattern.appendChild(this_roi);
    
    % ROI (targets)
    for j = 1:trg.p.numCellsPerCluster
        this_roi = docNode.createElement('ROI');

        % ROI info
        this_roi.setAttribute('subID',num2str(j));
        this_roi.setAttribute('centerX',num2str(trg.f.med2(trg.clustering(j,i))));
        this_roi.setAttribute('centerY',num2str(trg.f.med1(trg.clustering(j,i))));
        this_pattern.appendChild(this_roi);
    end
end


%% SLM Sequences, SequenceEpoch

this_element = docNode.createElement('SLMSequences');
toc.appendChild(this_element);

for i=1:trg.p.numStimTrials  
    this_epoch = docNode.createElement('SequenceEpoch');
       
    % Epoch info
    temp = trg.sequenceClusters(:,trg.sequenceOrder(i));
    this_epoch.setAttribute('sequenceID',num2str(i));
    this_epoch.setAttribute('sequence',char(join(string(temp),':')));
    this_epoch.setAttribute('sequenceEpochCount',trg.xml.sequenceEpochCount);
    this_element.appendChild(this_epoch);
end

%% Write xml file

xmlwrite(savePath,docNode);
%type(savePath);

end


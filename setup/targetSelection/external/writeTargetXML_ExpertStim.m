function writeTargetXML_ExpertStim(trg, savePath)
% required trg fields: p, xml, f, clustering


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
for i=1:trg.p.numClusterSlots*2    
    this_pattern = docNode.createElement('Pattern');
       
    % Pattern info
    temp = [trg.p.col.AB_rainbow;trg.p.col.XY_rainbow];
    this_pattern.setAttribute('name',[trg.xml.name_root,num2str(i,'%03.f')]);
    this_pattern.setAttribute('patternID',num2str(i));
    this_pattern.setAttribute('shape',trg.xml.shape);
    this_pattern.setAttribute('roiWidthPx',trg.xml.roiWidthPx);
    this_pattern.setAttribute('roiHeightPx',trg.xml.roiHeightPx);
    this_pattern.setAttribute('red',num2str(round(temp(i,1)*255)));
    this_pattern.setAttribute('green',num2str(round(temp(i,2)*255)));
    this_pattern.setAttribute('blue',num2str(round(temp(i,3)*255)));
    this_pattern.setAttribute('pxSpacing',trg.xml.pxSpacing);
    this_pattern.setAttribute('durationMS',trg.xml.durationMS);
    this_pattern.setAttribute('iterations',trg.xml.iterations);
    this_pattern.setAttribute('power',trg.xml.power{i});
    this_pattern.setAttribute('prePatIdleMS',trg.xml.prePatIdleMS);
    this_pattern.setAttribute('postPatIdleMS',trg.xml.postPatIdleMS);
    this_pattern.setAttribute('preIteIdleMS',trg.xml.preIteIdleMS);
    this_pattern.setAttribute('postIteIdleMS',trg.xml.postIteIdleMS);
    this_pattern.setAttribute('measurePowerMW',trg.xml.measurePowerMW);
    this_pattern.setAttribute('measurePowerMWPerUM2',trg.xml.measurePowerMWPerUM2);
    this_element.appendChild(this_pattern);
    
    % ROI (galvo)
    this_roi = docNode.createElement('ROI');
    this_roi.setAttribute('subID',num2str(0));
    this_roi.setAttribute('centerX',num2str(trg.galvo(2,i)));
    this_roi.setAttribute('centerY',num2str(trg.galvo(1,i)));
    this_pattern.appendChild(this_roi);
    
    % ROI (targets)
    temp = trg.clustering(:,i);
    temp = temp(~isnan(temp));
    for j = 1:length(temp)
        this_roi = docNode.createElement('ROI');

        % ROI info
        this_roi.setAttribute('subID',num2str(j));
        if trg.clustering(j,i)~=0
            this_roi.setAttribute('centerX',num2str(trg.f.med2(trg.clustering(j,i))));
            this_roi.setAttribute('centerY',num2str(trg.f.med1(trg.clustering(j,i))));
        else
            this_roi.setAttribute('centerX',num2str(trg.cluster_centres(2,i)));
            this_roi.setAttribute('centerY',num2str(trg.cluster_centres(1,i)));
        end
        this_pattern.appendChild(this_roi);
    end
end


%% SLM Sequences, SequenceEpoch

this_element = docNode.createElement('SLMSequences');
toc.appendChild(this_element);

for i=1:trg.info.numStimTrials
    this_epoch = docNode.createElement('SequenceEpoch');
       
    % Epoch info
    this_epoch.setAttribute('sequenceID',num2str(i));
    this_epoch.setAttribute('sequence',char(join(string(trg.sequences(:,i)),':')));
    this_epoch.setAttribute('sequenceEpochCount',trg.xml.sequenceEpochCount);
    this_element.appendChild(this_epoch);
end

%% Write xml file

xmlwrite(savePath,docNode);
%type(savePath);

end


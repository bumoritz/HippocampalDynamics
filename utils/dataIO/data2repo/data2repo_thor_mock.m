function [path] = data2repo_thor_mock(info,path,paq_beh,task)

warning('Generating mock thor_beh file from paq_beh file as thor_beh file is missing.')

thor_beh.sync = paq_beh.sync;
thor_beh.numFrameTriggers_original = paq_beh.numFrameTriggers_original;
thor_beh.numFrameTriggers_postDropping = paq_beh.numFrameTriggers_postDropping;
thor_beh.numFrames = thor_beh.numFrameTriggers_postDropping;

if info.stimSession
    thor_beh.stimSequence = thor_beh.sync(find(task.var))+9;
    thor_beh.stimPatternOnset = nan(info.scope.patternsPerSequence,length(thor_beh.stimSequence));
    thor_beh.stimPatternOnset(1,:) = thor_beh.stimSequence;
    temp = [8;7;8;7;8;7;8;7;8;8;7;8;7;8;7;8;7;8;7];
    for i=2:info.scope.patternsPerSequence
        thor_beh.stimPatternOnset(i,:) = thor_beh.stimPatternOnset(i-1,:)+temp(i-1);
    end
    thor_beh.stimPatternComplete = nan(info.scope.patternsPerSequence,length(thor_beh.stimSequence));
    temp2= [3;2;3;2;3;3;3;3;3;3;2;3;2;3;2;3;2;3;3;3];
    for i=1:info.scope.patternsPerSequence
        thor_beh.stimPatternComplete(i,:) = thor_beh.stimPatternOnset(i,:)+temp2(i);
    end
    temp1 = thor_beh.stimPatternComplete(:);
    temp2 = thor_beh.stimPatternOnset(:);
    thor_beh.artefact = [];
    for i=1:length(temp1)
    	thor_beh.artefact = [thor_beh.artefact, temp2(i):temp1(i)];
    end
end

thor_beh = orderfields(thor_beh);
save([path.filepart_out,'thor_beh.mat'],'thor_beh','-v7.3');
disp(['--- Added thor_beh file to repo as ',[path.filepart_out,'thor_beh.mat'],'.'])

end
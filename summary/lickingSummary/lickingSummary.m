%function [licks] = lickingSummary(d_info,d,ops)

%%

latency = extractVariable(d,'paq_beh.licks_od2','cell','firstOfEach');
%latency = fillEmptyFields(latency);

rigour = extractVariable(d,'paq_beh.licks_od2','cell','number');
%rigour = fillEmptyFields(rigour);

task = extractVariable(d,'task.type','cell','all');


%%



%%

cell2mat(cellfun(@(c) c.perf.('blocks_general').H,{d{[find(d_info.group==0);find(d_info.group==1);find(d_info.group==2)],j}}','Uniform',0))
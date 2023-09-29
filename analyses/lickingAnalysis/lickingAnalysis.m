function [lick] = lickingAnalysis(info,ops,p,path,task,perf,paq_beh)

%% Lick latency

latency = extractVariable(paq_beh,'licks_od2','array','firstOfEach');

% latency
lick.latency.AB_H = latency(intersect(find(task.type==1),find(task.response=="H")));
lick.latency.XY_H = latency(intersect(find(task.type==2),find(task.response=="H")));
lick.latency.AY_FA = latency(intersect(find(task.type==3),find(task.response=="FA")));
lick.latency.XB_FA = latency(intersect(find(task.type==4),find(task.response=="FA")));
if info.stimSession
    lick.latency.AB_H_stim = latency(intersect(find(task.type==1),intersect(find(task.response=="H"),find(task.var==1))));
    lick.latency.XY_H_stim = latency(intersect(find(task.type==2),intersect(find(task.response=="H"),find(task.var==1))));
    lick.latency.AY_FA_stim = latency(intersect(find(task.type==3),intersect(find(task.response=="FA"),find(task.var==1))));
    lick.latency.XB_FA_stim = latency(intersect(find(task.type==4),intersect(find(task.response=="FA"),find(task.var==1))));
    lick.latency.AB_H_catch = latency(intersect(find(task.type==1),intersect(find(task.response=="H"),find(task.var==0))));
    lick.latency.XY_H_catch = latency(intersect(find(task.type==2),intersect(find(task.response=="H"),find(task.var==0))));
    lick.latency.AY_FA_catch = latency(intersect(find(task.type==3),intersect(find(task.response=="FA"),find(task.var==0))));
    lick.latency.XB_FA_catch = latency(intersect(find(task.type==4),intersect(find(task.response=="FA"),find(task.var==0))));
end
lick.latency = orderfields(lick.latency);

% latency_avg
temp = fields(lick.latency);
for k=1:length(temp)
    lick.latency_avg.(char(temp(k))) = nanmean(lick.latency.(char(temp(k))));
end
lick.latency_avg = orderfields(lick.latency_avg);

% latency_blocks
for b=1:info.task.numBlocks
    lick.latency_blocks.AB_H(b) = nanmean(latency(intersect(find(task.type==1),intersect(find(task.response=="H"),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
    lick.latency_blocks.XY_H(b) = nanmean(latency(intersect(find(task.type==2),intersect(find(task.response=="H"),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
    lick.latency_blocks.AY_FA(b) = nanmean(latency(intersect(find(task.type==3),intersect(find(task.response=="FA"),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
    lick.latency_blocks.XB_FA(b) = nanmean(latency(intersect(find(task.type==4),intersect(find(task.response=="FA"),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
    if info.stimSession
        lick.latency_blocks.AB_H_stim(b) = nanmean(latency(intersect(find(task.type==1),intersect(find(task.response=="H"),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.latency_blocks.XY_H_stim(b) = nanmean(latency(intersect(find(task.type==2),intersect(find(task.response=="H"),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.latency_blocks.AY_FA_stim(b) = nanmean(latency(intersect(find(task.type==3),intersect(find(task.response=="FA"),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.latency_blocks.XB_FA_stim(b) = nanmean(latency(intersect(find(task.type==4),intersect(find(task.response=="FA"),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.latency_blocks.AB_H_catch(b) = nanmean(latency(intersect(find(task.type==1),intersect(find(task.response=="H"),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.latency_blocks.XY_H_catch(b) = nanmean(latency(intersect(find(task.type==2),intersect(find(task.response=="H"),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.latency_blocks.AY_FA_catch(b) = nanmean(latency(intersect(find(task.type==3),intersect(find(task.response=="FA"),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.latency_blocks.XB_FA_catch(b) = nanmean(latency(intersect(find(task.type==4),intersect(find(task.response=="FA"),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
    end
end
lick.latency_blocks = orderfields(lick.latency_blocks);


%% Lick rigour

numLicks_od2 = extractVariable(paq_beh,'licks_od2','array','number');
numLicks_rw = extractVariable(paq_beh,'licks_rw','array','number');
rigour = numLicks_od2;
rigour(find(task.response=="H")) = numLicks_od2(find(task.response=="H")) - numLicks_rw(find(task.response=="H")) + 1;
rigour(min(find(isnan(paq_beh.sync))):end) = NaN;

% rigour
lick.rigour.AB_H = rigour(intersect(find(task.type==1),find(task.response=="H")));
lick.rigour.XY_H = rigour(intersect(find(task.type==2),find(task.response=="H")));
lick.rigour.AY_FA = rigour(intersect(find(task.type==3),find(task.response=="FA")));
lick.rigour.XB_FA = rigour(intersect(find(task.type==4),find(task.response=="FA")));
if info.stimSession
    lick.rigour.AB_H_stim = rigour(intersect(find(task.type==1),intersect(find(task.response=="H"),find(task.var==1))));
    lick.rigour.XY_H_stim = rigour(intersect(find(task.type==2),intersect(find(task.response=="H"),find(task.var==1))));
    lick.rigour.AY_FA_stim = rigour(intersect(find(task.type==3),intersect(find(task.response=="FA"),find(task.var==1))));
    lick.rigour.XB_FA_stim = rigour(intersect(find(task.type==4),intersect(find(task.response=="FA"),find(task.var==1))));
    lick.rigour.AB_H_catch = rigour(intersect(find(task.type==1),intersect(find(task.response=="H"),find(task.var==0))));
    lick.rigour.XY_H_catch = rigour(intersect(find(task.type==2),intersect(find(task.response=="H"),find(task.var==0))));
    lick.rigour.AY_FA_catch = rigour(intersect(find(task.type==3),intersect(find(task.response=="FA"),find(task.var==0))));
    lick.rigour.XB_FA_catch = rigour(intersect(find(task.type==4),intersect(find(task.response=="FA"),find(task.var==0))));
end
lick.rigour = orderfields(lick.rigour);

% rigour_avg
temp = fields(lick.rigour);
for k=1:length(temp)
    lick.rigour_avg.(char(temp(k))) = nanmean(lick.rigour.(char(temp(k))));
end
lick.rigour_avg = orderfields(lick.rigour_avg);

% rigour_blocks
for b=1:info.task.numBlocks
    lick.rigour_blocks.AB_H(b) = nanmean(rigour(intersect(find(task.type==1),intersect(find(task.response=="H"),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
    lick.rigour_blocks.XY_H(b) = nanmean(rigour(intersect(find(task.type==2),intersect(find(task.response=="H"),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
    lick.rigour_blocks.AY_FA(b) = nanmean(rigour(intersect(find(task.type==3),intersect(find(task.response=="FA"),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
    lick.rigour_blocks.XB_FA(b) = nanmean(rigour(intersect(find(task.type==4),intersect(find(task.response=="FA"),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
    if info.stimSession
        lick.rigour_blocks.AB_H_stim(b) = nanmean(rigour(intersect(find(task.type==1),intersect(find(task.response=="H"),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.rigour_blocks.XY_H_stim(b) = nanmean(rigour(intersect(find(task.type==2),intersect(find(task.response=="H"),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.rigour_blocks.AY_FA_stim(b) = nanmean(rigour(intersect(find(task.type==3),intersect(find(task.response=="FA"),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.rigour_blocks.XB_FA_stim(b) = nanmean(rigour(intersect(find(task.type==4),intersect(find(task.response=="FA"),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.rigour_blocks.AB_H_catch(b) = nanmean(rigour(intersect(find(task.type==1),intersect(find(task.response=="H"),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.rigour_blocks.XY_H_catch(b) = nanmean(rigour(intersect(find(task.type==2),intersect(find(task.response=="H"),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.rigour_blocks.AY_FA_catch(b) = nanmean(rigour(intersect(find(task.type==3),intersect(find(task.response=="FA"),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
        lick.rigour_blocks.XB_FA_catch(b) = nanmean(rigour(intersect(find(task.type==4),intersect(find(task.response=="FA"),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)))));
    end
end
% temp = fields(lick.rigour_blocks);
% for k=1:length(temp)
%     if k<=(min(find(isnan(paq_beh.sync)))-1)/info.task.trialsPerBlock
%         temp2 = lick.rigour_blocks.(char(temp(k)));
%         temp2(isnan(temp2))=0;
%         lick.rigour_blocks.(char(temp(k))) = temp2;
%     end
% end
lick.rigour_blocks = orderfields(lick.rigour_blocks);


%% Lick probability

licks_binary = numLicks_od2>0;

% lick probability
lick.lprob.AB = nanmean(licks_binary(find(task.type==1)));
lick.lprob.XY = nanmean(licks_binary(find(task.type==2)));
lick.lprob.AY = nanmean(licks_binary(find(task.type==3)));
lick.lprob.XB = nanmean(licks_binary(find(task.type==4)));
if info.stimSession
    lick.lprob.AB_stim = nanmean(licks_binary(intersect(find(task.type==1),find(task.var==1))));
    lick.lprob.XY_stim = nanmean(licks_binary(intersect(find(task.type==2),find(task.var==1))));
    lick.lprob.AY_stim = nanmean(licks_binary(intersect(find(task.type==3),find(task.var==1))));
    lick.lprob.XB_stim = nanmean(licks_binary(intersect(find(task.type==4),find(task.var==1))));
    lick.lprob.AB_catch = nanmean(licks_binary(intersect(find(task.type==1),find(task.var==0))));
    lick.lprob.XY_catch = nanmean(licks_binary(intersect(find(task.type==2),find(task.var==0))));
    lick.lprob.AY_catch = nanmean(licks_binary(intersect(find(task.type==3),find(task.var==0))));
    lick.lprob.XB_catch = nanmean(licks_binary(intersect(find(task.type==4),find(task.var==0))));
end
lick.lprob = orderfields(lick.lprob);

% lprob_blocks
for b=1:info.task.numBlocks
    lick.lprob_blocks.AB(b) = nanmean(licks_binary(intersect(find(task.type==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)));
    lick.lprob_blocks.XY(b) = nanmean(licks_binary(intersect(find(task.type==2),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)));
    lick.lprob_blocks.AY(b) = nanmean(licks_binary(intersect(find(task.type==3),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)));
    lick.lprob_blocks.XB(b) = nanmean(licks_binary(intersect(find(task.type==4),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock)));
    if info.stimSession
        lick.lprob_blocks.AB_stim(b) = nanmean(licks_binary(intersect(find(task.type==1),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
        lick.lprob_blocks.XY_stim(b) = nanmean(licks_binary(intersect(find(task.type==2),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
        lick.lprob_blocks.AY_stim(b) = nanmean(licks_binary(intersect(find(task.type==3),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
        lick.lprob_blocks.XB_stim(b) = nanmean(licks_binary(intersect(find(task.type==4),intersect(find(task.var==1),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
        lick.lprob_blocks.AB_catch(b) = nanmean(licks_binary(intersect(find(task.type==1),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
        lick.lprob_blocks.XY_catch(b) = nanmean(licks_binary(intersect(find(task.type==2),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
        lick.lprob_blocks.AY_catch(b) = nanmean(licks_binary(intersect(find(task.type==3),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
        lick.lprob_blocks.XB_catch(b) = nanmean(licks_binary(intersect(find(task.type==4),intersect(find(task.var==0),(b-1)*info.task.trialsPerBlock+1:b*info.task.trialsPerBlock))));
    end
end
lick.lprob_blocks = orderfields(lick.lprob_blocks);


%% Response bias

% engagement
lick.engagement = nanmean([perf.type1.H,perf.type2.H,perf.type3.FA,perf.type4.FA]); % = avg response probability

% response bias - general (pos. is more lick than non-lick)
temp1 = nanmean([perf.type1.H,perf.type2.H,perf.type3.FA,perf.type4.FA]);
temp2 = nanmean([perf.type1.M,perf.type2.M,perf.type3.CR,perf.type4.CR]);
lick.rbias.general = temp1-temp2/(temp1+temp2);

% response bias - A vs X (pos. is more A)
temp1 = nanmean([perf.type1.H,perf.type3.FA]);
temp2 = nanmean([perf.type2.H,perf.type4.FA]);
lick.rbias.odourA = temp1-temp2/(temp1+temp2);

% response bias - B vs Y (pos. is more B)
temp1 = nanmean([perf.type1.H,perf.type4.FA]);
temp2 = nanmean([perf.type2.H,perf.type3.FA]);
lick.rbias.odourB = temp1-temp2/(temp1+temp2);

if info.stimSession
    lick.engagement_stim = nanmean([perf.type1.H_stim,perf.type2.H_stim,perf.type3.FA_stim,perf.type4.FA_stim]);  
    temp1 = nanmean([perf.type1.H_stim,perf.type2.H_stim,perf.type3.FA_stim,perf.type4.FA_stim]);
    temp2 = nanmean([perf.type1.M_stim,perf.type2.M_stim,perf.type3.CR_stim,perf.type4.CR_stim]);
    lick.rbias.general_stim = temp1-temp2/(temp1+temp2);
    temp1 = nanmean([perf.type1.H_stim,perf.type3.FA_stim]);
    temp2 = nanmean([perf.type2.H_stim,perf.type4.FA_stim]);
    lick.rbias.odourA_stim = temp1-temp2/(temp1+temp2);
    temp1 = nanmean([perf.type1.H_stim,perf.type4.FA_stim]);
    temp2 = nanmean([perf.type2.H_stim,perf.type3.FA_stim]);
    lick.rbias.odourB_stim = temp1-temp2/(temp1+temp2);    
    lick.engagement_catch = nanmean([perf.type1.H_catch,perf.type2.H_catch,perf.type3.FA_catch,perf.type4.FA_catch]);  
    temp1 = nanmean([perf.type1.H_catch,perf.type2.H_catch,perf.type3.FA_catch,perf.type4.FA_catch]);
    temp2 = nanmean([perf.type1.M_catch,perf.type2.M_catch,perf.type3.CR_catch,perf.type4.CR_catch]);
    lick.rbias.general_catch = temp1-temp2/(temp1+temp2);
    temp1 = nanmean([perf.type1.H_catch,perf.type3.FA_catch]);
    temp2 = nanmean([perf.type2.H_catch,perf.type4.FA_catch]);
    lick.rbias.odourA_catch = temp1-temp2/(temp1+temp2);
    temp1 = nanmean([perf.type1.H_catch,perf.type4.FA_catch]);
    temp2 = nanmean([perf.type2.H_catch,perf.type3.FA_catch]);
    lick.rbias.odourB_catch = temp1-temp2/(temp1+temp2);    
end
lick.rbias = orderfields(lick.rbias);

% rbias_blocks
for b=1:info.task.numBlocks
    lick.engagement_blocks(b) = nanmean([perf.blocks_type1.H(b),perf.blocks_type2.H(b),perf.blocks_type3.FA(b),perf.blocks_type4.FA(b)]);

    temp1 = nanmean([perf.blocks_type1.H(b),perf.blocks_type2.H(b),perf.blocks_type3.FA(b),perf.blocks_type4.FA(b)]);
    temp2 = nanmean([perf.blocks_type1.M(b),perf.blocks_type2.M(b),perf.blocks_type3.CR(b),perf.blocks_type4.CR(b)]);
    lick.rbias_blocks.general(b) = temp1-temp2/(temp1+temp2);  
    temp1 = nanmean([perf.blocks_type1.H(b),perf.blocks_type3.FA(b)]);
    temp2 = nanmean([perf.blocks_type2.M(b),perf.blocks_type4.CR(b)]);
    lick.rbias_blocks.odourA(b) = temp1-temp2/(temp1+temp2);
    temp1 = nanmean([perf.blocks_type1.H(b),perf.blocks_type4.FA(b)]);
    temp2 = nanmean([perf.blocks_type2.M(b),perf.blocks_type3.CR(b)]);
    lick.rbias_blocks.odourB(b) = temp1-temp2/(temp1+temp2);
    
    if info.stimSession
        temp1 = nanmean([perf.blocks_type1.H_stim(b),perf.blocks_type2.H_stim(b),perf.blocks_type3.FA_stim(b),perf.blocks_type4.FA_stim(b)]);
        temp2 = nanmean([perf.blocks_type1.M_stim(b),perf.blocks_type2.M_stim(b),perf.blocks_type3.CR_stim(b),perf.blocks_type4.CR_stim(b)]);
        lick.rbias_blocks.general_stim(b) = temp1-temp2/(temp1+temp2);  
        temp1 = nanmean([perf.blocks_type1.H_stim(b),perf.blocks_type3.FA_stim(b)]);
        temp2 = nanmean([perf.blocks_type2.M_stim(b),perf.blocks_type4.CR_stim(b)]);
        lick.rbias_blocks.odourA_stim(b) = temp1-temp2/(temp1+temp2);
        temp1 = nanmean([perf.blocks_type1.H_stim(b),perf.blocks_type4.FA_stim(b)]);
        temp2 = nanmean([perf.blocks_type2.M_stim(b),perf.blocks_type3.CR_stim(b)]);
        lick.rbias_blocks.odourB_stim(b) = temp1-temp2/(temp1+temp2);       
        temp1 = nanmean([perf.blocks_type1.H_catch(b),perf.blocks_type2.H_catch(b),perf.blocks_type3.FA_catch(b),perf.blocks_type4.FA_catch(b)]);
        temp2 = nanmean([perf.blocks_type1.M_catch(b),perf.blocks_type2.M_catch(b),perf.blocks_type3.CR_catch(b),perf.blocks_type4.CR_catch(b)]);
        lick.rbias_blocks.general_catch(b) = temp1-temp2/(temp1+temp2);  
        temp1 = nanmean([perf.blocks_type1.H_catch(b),perf.blocks_type3.FA_catch(b)]);
        temp2 = nanmean([perf.blocks_type2.M_catch(b),perf.blocks_type4.CR_catch(b)]);
        lick.rbias_blocks.odourA_catch(b) = temp1-temp2/(temp1+temp2);
        temp1 = nanmean([perf.blocks_type1.H_catch(b),perf.blocks_type4.FA_catch(b)]);
        temp2 = nanmean([perf.blocks_type2.M_catch(b),perf.blocks_type3.CR_catch(b)]);
        lick.rbias_blocks.odourB_catch(b) = temp1-temp2/(temp1+temp2);            
    end
end
lick.rbias_blocks = orderfields(lick.rbias_blocks);


%% Save

lick = orderfields(lick);
save([path.filepart_out,'lick.mat'],'lick','-v7.3');
disp(['--- Saved lick file as ',[path.filepart_out,'lick.mat'],'.'])


%% Make figures

% licking analysis figure
F = lickingAnalysisFigure(lick,perf,p,info);
savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','lickingAnalysis.fig']);
saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','lickingAnalysis.png']);
disp(['--- Saved licking analysis figure to ',path.filepart_outX,'plots.'])

% licking analysis figure - stimMode
if info.stimSession
    F = lickingAnalysisFigure_stim(lick,perf,p,info);
    savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','lickingAnalysis_stim.fig']);
    saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','lickingAnalysis_stim.png']);
    disp(['--- Saved licking analysis figure (stim version) to ',path.filepart_outX,'plots.'])
end

if ops.close_figures
    close all;
end
end














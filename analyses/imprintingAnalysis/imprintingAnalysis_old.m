function imprintingAnalysis_old


%% temp

%%% start with running Preparations section of sequence cell analysis
% trg = trg_nonrigid;
% 
% 
% %% kind of from snake plot utility
% 
% data_A = sca.traces_A;
% data_X = sca.traces_X;
% 
% avgTraces_A_1 = nanmean(data_A(:,:, find(task.var(find(task.odour1=="A"))) ),3);
% avgTraces_A_2 = nanmean(data_A(:,:, find(~task.var(find(task.odour1=="A"))) ),3);
% avgTraces_X_1 = nanmean(data_X(:,:, find(task.var(find(task.odour1=="X"))) ),3);
% avgTraces_X_2 = nanmean(data_X(:,:, find(~task.var(find(task.odour1=="X"))) ),3);
% 
% % 1 is stim, 2 is catch
% normAvgTraces_A_1_Aseq = (avgTraces_A_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_A_2_Aseq = (avgTraces_A_2-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_X_1_Aseq = (avgTraces_X_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_X_2_Aseq = (avgTraces_X_2-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_A_1_Xseq = (avgTraces_A_1-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_A_2_Xseq = (avgTraces_A_2-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_X_1_Xseq = (avgTraces_X_1-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_X_2_Xseq = (avgTraces_X_2-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
% 
% temp = trg.idcs_targetedCells(:,trg.sequenceClusters(:,1));
% idcs_sorted_Aseq = rmmissing(temp(:));
% numSel_A = NaN;
% temp = trg.idcs_targetedCells(:,trg.sequenceClusters(:,2));
% idcs_sorted_Xseq = rmmissing(temp(:));
% numSel_X = NaN;
% 
% AstimAseq = normAvgTraces_A_1_Aseq(idcs_sorted_Aseq,:);
% AcatchAseq = normAvgTraces_A_2_Aseq(idcs_sorted_Aseq,:);
% XstimAseq = normAvgTraces_X_1_Aseq(idcs_sorted_Aseq,:);
% XcatchAseq = normAvgTraces_X_2_Aseq(idcs_sorted_Aseq,:);
% AstimXseq = normAvgTraces_A_1_Xseq(idcs_sorted_Xseq,:);
% AcatchXseq = normAvgTraces_A_2_Xseq(idcs_sorted_Xseq,:);
% XstimXseq = normAvgTraces_X_1_Xseq(idcs_sorted_Xseq,:);
% XcatchXseq = normAvgTraces_X_2_Xseq(idcs_sorted_Xseq,:);
% 
% 
% %% A seq correlation 
% 
% corr_AstimAseq_AcatchAseq = nan(length(idcs_sorted_Aseq),1);
% for i=1:length(idcs_sorted_Aseq)
%     corr_AstimAseq_AcatchAseq(i) = corr(AstimAseq(i,:)',AcatchAseq(i,:)');
% end
% 
% corr_AstimAseq_XcatchAseq = nan(length(idcs_sorted_Aseq),1);
% for i=1:length(idcs_sorted_Aseq)
%     corr_AstimAseq_XcatchAseq(i) = corr(AstimAseq(i,:)',XcatchAseq(i,:)');
% end
% 
% figure;
% v = violinplot([corr_AstimAseq_AcatchAseq,corr_AstimAseq_XcatchAseq],{'corr(A-stim,A-nostim)','corr(A-stim,X-nostim)'});
% title('SeqA cell activity correlation between stim trials and no-stim trials')
% 
% 
% %% X seq correlation 
% 
% corr_XstimXseq_XcatchXseq = nan(length(idcs_sorted_Xseq),1);
% for i=1:length(idcs_sorted_Xseq)
%     corr_XstimXseq_XcatchXseq(i) = corr(XstimXseq(i,:)',XcatchXseq(i,:)');
% end
% 
% corr_XstimXseq_AcatchXseq = nan(length(idcs_sorted_Xseq),1);
% for i=1:length(idcs_sorted_Xseq)
%     corr_XstimXseq_AcatchXseq(i) = corr(XstimXseq(i,:)',AcatchXseq(i,:)');
% end
% 
% figure;
% v = violinplot([corr_XstimXseq_XcatchXseq,corr_XstimXseq_AcatchXseq],{'corr(X-stim,X-nostim)','corr(X-stim,A-nostim)'});
% title('SeqX cell activity correlation between stim trials and no-stim trials')
% 
% 
% %% --- %%% --- %%%
% 
% %% Pop vector correlation - step size of 5
% 
% figure;
% 
% binSize = 3;
% 
% temp = movmean(avgTraces_A_1(idcs_sorted_Aseq,sca.prop.analysisWindow),binSize,2,'omitnan');
% templateA = temp(:,floor(binSize/2)+1:binSize:end,:);
% 
% catch_trials = find(~task.var(find(task.odour1=="A")));
% catchDataA = {};
% for i=1:16
%     temp = movmean(   nanmean(data_A(idcs_sorted_Aseq,sca.prop.analysisWindow, catch_trials((i-1)*5+1:i*5) ),3)   ,binSize,2,'omitnan');
%     catchDataA{i} = temp(:,floor(binSize/2)+1:binSize:end,:);
% end
% 
% corr_templateA_catchDataA = nan(size(templateA,2),16);
% for i=1:16
%     for j=1:size(templateA,2)
%         corr_templateA_catchDataA(j,i) = corr(templateA(:,j),catchDataA{i}(:,j));
%     end
% end
% 
% subplot(1,2,1)
% imagesc(corr_templateA_catchDataA)
% colorbar
% xlabel('Catch trial bin (1 bin = 5 catch trials)')
% ylabel('Time bin (1 bin = 300 ms)')
% title('PV corr for A: catch trial bin with template from all stim trials')
% 
% 
% binSize = 3;
% temp = movmean(avgTraces_X_1(idcs_sorted_Xseq,sca.prop.analysisWindow),binSize,2,'omitnan');
% templateX = temp(:,floor(binSize/2)+1:binSize:end,:);
% 
% catch_trials = find(~task.var(find(task.odour1=="X")));
% catchDataX = {};
% for i=1:16
%     temp = movmean(   nanmean(data_X(idcs_sorted_Xseq,sca.prop.analysisWindow, catch_trials((i-1)*5+1:i*5) ),3)   ,binSize,2,'omitnan');
%     catchDataX{i} = temp(:,floor(binSize/2)+1:binSize:end,:);
% end
% 
% corr_templateX_catchDataX = nan(size(templateX,2),16);
% for i=1:16
%     for j=1:size(templateX,2)
%         corr_templateX_catchDataX(j,i) = corr(templateX(:,j),catchDataX{i}(:,j));
%     end
% end
% 
% subplot(1,2,2)
% imagesc(corr_templateX_catchDataX)
% colorbar
% xlabel('Catch trial bin (1 bin = 5 catch trials)')
% ylabel('Time bin (1 bin = 300 ms)')
% title('PV corr for X: catch trial bin with template from all stim trials')
% 
% 
% %% Pop vector correlation - step size of 1
% 
% figure;
% 
% binSize = 3;
% 
% temp = movmean(avgTraces_A_1(idcs_sorted_Aseq,sca.prop.analysisWindow),binSize,2,'omitnan');
% templateA = temp(:,floor(binSize/2)+1:binSize:end,:);
% 
% catch_trials = find(~task.var(find(task.odour1=="A")));
% catchDataA = {};
% for i=1:size(data_A,3)/5
%     temp = movmean(   nanmean(data_A(idcs_sorted_Aseq,sca.prop.analysisWindow, catch_trials((i-1)*1+1:i*1) ),3)   ,binSize,2,'omitnan');
%     catchDataA{i} = temp(:,floor(binSize/2)+1:binSize:end,:);
% end
% 
% corr_templateA_catchDataA = nan(size(templateA,2),size(data_A,3)/5);
% for i=1:size(data_A,3)/5
%     for j=1:size(templateA,2)
%         corr_templateA_catchDataA(j,i) = corr(templateA(:,j),catchDataA{i}(:,j));
%     end
% end
% 
% subplot(1,2,1)
% imagesc(corr_templateA_catchDataA)
% colorbar
% xlabel('Catch trial')
% ylabel('Time bin (1 bin = 300 ms)')
% title('PV corr for A: catch trial bin with template from all stim trials')
% 
% 
% binSize = 3;
% temp = movmean(avgTraces_X_1(idcs_sorted_Xseq,sca.prop.analysisWindow),binSize,2,'omitnan');
% templateX = temp(:,floor(binSize/2)+1:binSize:end,:);
% 
% catch_trials = find(~task.var(find(task.odour1=="X")));
% catchDataX = {};
% for i=1:size(data_X,3)/5
%     temp = movmean(   nanmean(data_X(idcs_sorted_Xseq,sca.prop.analysisWindow, catch_trials((i-1)*1+1:i*1) ),3)   ,binSize,2,'omitnan');
%     catchDataX{i} = temp(:,floor(binSize/2)+1:binSize:end,:);
% end
% 
% corr_templateX_catchDataX = nan(size(templateX,2),size(data_X,3)/5);
% for i=1:size(data_X,3)/5
%     for j=1:size(templateX,2)
%         corr_templateX_catchDataX(j,i) = corr(templateX(:,j),catchDataX{i}(:,j));
%     end
% end
% 
% subplot(1,2,2)
% imagesc(corr_templateX_catchDataX)
% colorbar
% xlabel('Catch trial')
% ylabel('Time bin (1 bin = 300 ms)')
% title('PV corr for X: catch trial bin with template from all stim trials')
% %%
% 
% figure;
% hold on
% plot(nanmean(corr_templateA_catchDataA,1),'Color',p.col.odourA)
% plot(nanmean(corr_templateX_catchDataX,1),'Color',p.col.odourX)
% xlabel('Catch trial bin (1 bin = 5 catch trials)')
% ylabel('Correlation')
% title('PV corr: catch trial bin with template from corresponding stim trials')
% hold off
% 
% 
% %%
% temp = nchoosek(1:prop.numTrials/100,2)
% stability = ones(prop.numTrials/100,prop.numTrials/100);
% for i=1:length(temp)
%     temp3 = size(sca_100t{1,temp(i,1)}.avgTraces_A(find(s2p_meta.iscell(:,1)),:),2);
%     temp2 = nan(temp3,1);
%     for j=1:temp3
%         temp2(j) = corr(sca_100t{1,temp(i,1)}.avgTraces_A(find(s2p_meta.iscell(:,1)),j),sca_100t{1,temp(i,2)}.avgTraces_A(find(s2p_meta.iscell(:,1)),j),'Rows','complete');
%     end
%     stability(temp(i,1),temp(i,2)) = nanmean(temp2);
%     stability(temp(i,2),temp(i,1)) = nanmean(temp2);
%     i
% end


%%





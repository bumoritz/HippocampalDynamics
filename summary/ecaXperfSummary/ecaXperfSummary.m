%% Preparations

if ~exist([path.root_summary,'plots\eca'],'dir')
    mkdir([path.root_summary,'plots\eca']);
end


%% Select data

this_field = 'dyn';

suptitles = {'1^{st}','A','X','sens','temp','A_{sens}','A_{temp}','X_{sens}','X_{temp}','2^{nd}','B','Y','int','A_{int}','X_{int}','R+L','R','L','run','v','a'};

suptitles2 = {'A+ (nLR)','X+ (nLR)','B+ (nLR)','Y+ (nLR)','A_{sens}+ (nLR)','X_{sens}+ (nLR)',...
    'int+ (nLR)','int- (nLR)','A_{int}+ (nLR)','A_{int}- (nLR)','X_{int}+ (nLR)','X_{int}- (nLR)',...
    'R+ (nL)','R- (nL)',...
    'A+ (nX+)','X+ (nA+)','B+ (nY+)','Y+ (nB+)','A_{sens}+ (nX_{sens}+)','X_{sens}+ (nA_{sens}+)'};

conditions = {'A_H_M','A_M_M','A_CR_M','A_FA_M',...
            'X_H_M','X_M_M','X_CR_M','X_FA_M'};
titles = {'AB trials (matched)','AY trials (matched)','XY trials (matched)','XB trials (matched)'};
cols = {p.col.AB,mean([p.col.AB;p.col.white]),p.col.AY,mean([p.col.AY;p.col.white]),...
    p.col.XY,mean([p.col.XY;p.col.white]),p.col.XB,mean([p.col.XB;p.col.white])};

ref_struct = d{1,1}.eca_all.dyn.A_H_M;
these_fields = fields(ref_struct);
for t=4
    for g=1:length(ref_struct.(these_fields{t}))

        %% Extract data

        activity = {};
        for k=1:length(conditions)
            activity.(conditions{k}) = nan(d_info.numAnimals,p.general.numBins);
        end
        for i=1:d_info.numAnimals
            try
                for k=1:length(conditions)
                    activity.(conditions{k})(i,:) = d{i,1}.eca_all.(this_field).(conditions{k}).(these_fields{t}){g};
                end
            catch
            end
        end


        %% Make figure

        nrows = 2; ncols = 2;
        F = default_figure([20,0.5,20,9.9]);
        
        temp0 = nan(length(conditions)/2,2);
        for k=1:length(conditions)/2
            subplot(nrows,ncols,k)
            hold on

            temp=shadedErrorBar(1:size(activity.(conditions{k*2-1}),2),nanmean(activity.(conditions{k*2-1}),1),nansem(activity.(conditions{k*2-1}),1),'lineProps',cols{k*2-1}); temp.mainLine.LineWidth = 2;  
            temp=shadedErrorBar(1:size(activity.(conditions{k*2}),2),nanmean(activity.(conditions{k*2}),1),nansem(activity.(conditions{k*2}),1),'lineProps',cols{k*2}); temp.mainLine.LineWidth = 2;  

            %ylim([0.35,0.5])
            temp_info.task.trialStructure.tOdour1 = 0.3;
            temp_info.task.trialStructure.tGap = 5;
            temp_info.task.trialStructure.tOdour2 = 0.3;
            temp_info.task.trialStructure.tRespDelay = 0.5;
            temp_info.task.trialStructure.tRespWindow = 1;
            taskLines(p,temp_info);
            temp0(k,:) = get(gca,'ylim');
            ylabel('Mean activity')
            title([titles{k}])
        end
        for k=1:length(conditions)/2
            subplot(nrows,ncols,k)
            ylim([min(temp0(:)),max(temp0(:))])
        end
        if t==1
            suptitle(['Encoding cell activity, ',suptitles{g},', +/-'])
            filename = ['encodingCellActivity_',suptitles{g},'_+-'];
        elseif t==2
            suptitle(['Encoding cell activity, ',suptitles{g},', +'])
            filename = ['encodingCellActivity_',suptitles{g},'_+'];
        elseif t==3
            suptitle(['Encoding cell activity, ',suptitles{g},', -'])
            filename = ['encodingCellActivity_',suptitles{g},'_-'];
        elseif t==4
            suptitle(['Encoding cell activity, ',suptitles2{g}])
            filename = ['encodingCellActivity_',suptitles2{g}];
        end
        
        savefig(F,[path.root_summary,'plots\eca\',filename,'.fig']);
        saveas(F,[path.root_summary,'plots\eca\',filename,'.png']);
    end
end
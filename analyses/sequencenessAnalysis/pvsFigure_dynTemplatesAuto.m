%%

this_colormap = flipud(pink); %jet %flipud(pink)


%% 

this_metric = 'Pearson'; % 'Pearson','cosTheta'

nrows = 4; ncols = 4;
F = default_figure([20,0.5,20,9.9]);

these_trialTypes = {'AB','AY','XY','XB'};
n=0;
for r=1:nrows
    for c=1:ncols
        n=n+1;
        subplot(nrows,ncols,n);

        this_data = povSim.(['dynTemplatesAuto_',this_metric]).iscells_full.([these_trialTypes{r},'_',these_trialTypes{c}]);
        
        imshow(this_data);
        %hold on; taskLines(p,info,1,1);
        set(gca,'CLim',[0,1]);
        colormap(gca,this_colormap); 
        daspect([1,1,1])
        colorbar;
        
        title([these_trialTypes{r},' - ',these_trialTypes{c}])
    end
end


%%

this_metric = 'Pearson'; % 'Pearson','cosTheta'

nrows = 2; ncols = 2;
F = default_figure([20,0.5,20,9.9]);

these_trialTypes = {'A','X'};
n=0;
for r=1:nrows
    for c=1:ncols
        n=n+1;
        subplot(nrows,ncols,n);

        this_data = povSim.(['dynTemplatesAuto_',this_metric]).iscells_AW.([these_trialTypes{r},'_',these_trialTypes{c}]);
        
        imshow(this_data);
%         hold on; taskLines(p,info,1,1);
        set(gca,'CLim',[0,1]);
        colormap(gca,this_colormap); 
        daspect([1,1,1])
        colorbar;
        
        title([these_trialTypes{r},' - ',these_trialTypes{c}])
    end
end


%%



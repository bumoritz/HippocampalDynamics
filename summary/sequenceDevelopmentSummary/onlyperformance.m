%% First try with mixed-effect model

this_data = correct(:,blocks(1)+1:end);

tbl_in.y = this_data(:);
temp = repmat([1:13,1:13],d_info.numAnimals,1);
tbl_in.block = temp(:);
tbl_in.animal = nominal(repmat((1:d_info.numAnimals)',size(this_data,2),1));
temp = repmat([1*ones(1,13),2*ones(1,13)],d_info.numAnimals,1);
tbl_in.switch = nominal(temp(:));
temp = nan(size(this_data));
for j=2:5
    if j==2
        temp(:,1:8) = repmat(d_info.running(:,j),1,8);
    elseif j==3
        temp(:,9:13) = repmat(d_info.running(:,j),1,5);
    elseif j==4
        temp(:,14:21) = repmat(d_info.running(:,j),1,8);
    elseif j==5
        temp(:,22:26) = repmat(d_info.running(:,j),1,5);
    end
end
tbl_in.running = nominal(temp(:));
temp = nan(size(this_data));
for i=1:d_info.numAnimals
    if d_info.group(i)==2
        temp(i,:) = zeros(1,26);
    elseif d_info.group(i)==7
        temp(i,:) = [1*ones(1,13),2*ones(1,13)];
    elseif d_info.group(i)==8
        temp(i,:) = [2*ones(1,13),1*ones(1,13)];
    end
end
tbl_in.stim = nominal(temp(:));

tbl = table(tbl_in.animal,tbl_in.switch,tbl_in.block,tbl_in.running,tbl_in.stim,tbl_in.y,...
    'VariableNames',{'animal','switch','block','running','stim','y'});

lme = fitlme(tbl,'y ~ running*block + stim*block + (1|animal) + (1|switch)')
% 


%% First try with mixed-effect model

this_data = correct(:,blocks(1)+1:end);

tbl_in.y = this_data(:);
temp = repmat([1:13,1:13],d_info.numAnimals,1);
tbl_in.block = temp(:);
temp = repmat([1:8,1:5,1:8,1:5],d_info.numAnimals,1);
tbl_in.block_sess = temp(:);
temp = repmat([0:1/(8-1):1,0:1/(5-1):1,0:1/(8-1):1,0:1/(5-1):1],d_info.numAnimals,1);
tbl_in.block_sess_norm = temp(:);
tbl_in.animal = nominal(repmat((1:d_info.numAnimals)',size(this_data,2),1));
temp = repmat([1*ones(1,13),2*ones(1,13)],d_info.numAnimals,1);
tbl_in.switch = nominal(temp(:));
temp = repmat([1*ones(1,8),2*ones(1,5),1*ones(1,8),2*ones(1,5)],d_info.numAnimals,1);
tbl_in.switchday = nominal(temp(:));
temp = repmat([1*ones(1,13),2*ones(1,13)],d_info.numAnimals,1);
tbl_in.switch = nominal(temp(:));
temp = nan(size(this_data));
for j=2:5
    if j==2
        temp(:,1:8) = repmat(d_info.running(:,j),1,8);
    elseif j==3
        temp(:,9:13) = repmat(d_info.running(:,j),1,5);
    elseif j==4
        temp(:,14:21) = repmat(d_info.running(:,j),1,8);
    elseif j==5
        temp(:,22:26) = repmat(d_info.running(:,j),1,5);
    end
end
tbl_in.running = nominal(temp(:));

tbl = table(tbl_in.animal,tbl_in.switch,tbl_in.block,tbl_in.running,tbl_in.y,...
    'VariableNames',{'animal','switch','block','running','y'});

tbl = table(tbl_in.animal,tbl_in.switch,tbl_in.switchday,tbl_in.block,tbl_in.block_sess,tbl_in.block_sess_norm,tbl_in.running,tbl_in.y,...
    'VariableNames',{'animal','switch','switchday','block','block_sess','block_sess_norm','running','y'});



lme = fitlme(tbl,'y ~ running*block + (1|animal) + (1|switch)')



% % 


lme = fitlme(tbl,'y ~ running*block_sess*switchday + (1|animal) + (1|switch)')


%%



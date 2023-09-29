%% Looking more at nem 

idx = idcs_B(76)
F = singleCellFigure(info,p,tng_all,idx,paq_beh);


%% New nem_all

%
this_testGroup = 10; % A
this_data = nem_all.testGroupResidual{this_testGroup}.coefs(:,1+(1:5));

these_idcs_pos = find(this_data(:,1)>0);
these_idcs_neg = find(this_data(:,1)<0);
these_idcs_zero = find(this_data(:,1)==0);

popActSuppPanel(avgTraces_all.A,find(iscell==1),these_idcs_pos,these_idcs_zero,these_idcs_neg,p,info,{})

figure;
plot(this_data')
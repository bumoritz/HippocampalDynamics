function [corr_r,corr_p] = fitLine(data_x,data_y,col)

[corr_r,corr_p] = corr(data_x,data_y,'Type','Pearson','Rows','Complete');

these_nans = any([isnan(data_x),isnan(data_y)],2);
temp = polyfit(data_x(~these_nans),data_y(~these_nans),1);
fit_x = linspace(nanmin(data_x),nanmax(data_x),1000);
fit_y = polyval(temp,fit_x);

if corr_p<0.05
    plot(fit_x,fit_y,'-','Color',col,'LineWidth',1)
else
    plot(fit_x,fit_y,'--','Color',col,'LineWidth',1)
end

end

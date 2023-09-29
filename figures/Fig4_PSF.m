%% Fig4_PSF
% start with Summary_master.m

%data_in = readtable('C:\SniffinHippo\Summary\other\PSF.xlsx'); % for xy and old z
data_in = readtable('C:\SniffinHippo\Summary\other\PSF2.xlsx'); % for new z

save_root_fig = [path.root_summary,'figures\Fig4_fig\'];
save_root_png = [path.root_summary,'figures\Fig4_png\'];
save_root_pdf = [path.root_summary,'figures\Fig4_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig4_txt\'];


%% Extract data - lateral PSF

% lateral PSF
psf_xy = table2array(rmmissing(data_in(:,2:6)));

% norm. lateral PSF
this_min = prctile(psf_xy,30);
this_max = nanmax(psf_xy);
normPsf_xy = (psf_xy - this_min) ./ (this_max-this_min);
avgNormPsf_xy = nanmean(normPsf_xy,2);

% norm. lateral PSF - all data points
normPsf_xy_all = normPsf_xy(:);
this_x = (1:size(normPsf_xy,1))';
this_x_all = repmat((1:size(normPsf_xy,1))',size(normPsf_xy,2),1);

% norm. lateral PSF - fit
normPsf_xy_fit = fit(this_x_all,normPsf_xy_all,'gauss1'); %this_lower = min(temp.YData); %this_upper = min(temp.YData);
clear log; normPsf_xy_hwhm = sqrt(log(2))*normPsf_xy_fit.c1; % 5.8481 -> 0.58 um


%% Extract data - axial PSF

% axial PSF
try
    psf_z = rmmissing(data_in.Var10);
    psf_z(82:103) = NaN; % removing erroneous data
catch
    psf_z = rmmissing(data_in.x50_267);
end
psf_z = smoothdata(psf_z,'gaussian',10); % window size=10*0.2um= 5um
psf_z = flip(psf_z);

% norm. axial PSF
this_min = prctile(psf_z,20);
this_max = nanmax(psf_z);
normPsf_z = (psf_z - this_min) ./ (this_max-this_min);
avgNormPsf_z = nanmean(normPsf_z,2);

% norm. axial PSF - all data points
normPsf_z_all = normPsf_z(:);
this_y = (1:size(normPsf_z,1))';
this_y_all = repmat((1:size(normPsf_z,1))',size(normPsf_z,2),1);

% norm. axial PSF - fit
temp = ~isnan(normPsf_z_all);
normPsf_z_fit = fit(this_y_all(temp),normPsf_z_all(temp),'gauss1'); %this_lower = min(temp.YData); %this_upper = min(temp.YData);
clear log; normPsf_z_hwhm = sqrt(log(2))*normPsf_z_fit.c1; % 38.0926 -> 7.6185  um


%% Fig4_PSF_xyResolutionCurve

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

shadedErrorBar(1:size(normPsf_xy,1),nanmean(normPsf_xy,2),nanstd(normPsf_xy,[],2),'lineProps',p.col.darkGray);
hold on; temp=plot(normPsf_xy_fit); temp.Color = p.col.photostim; legend off;
% [~,temp2] = max(temp.YData) % 484/1001 * 59 = 28.5 (10 steps are 1 um)
% double check: FOV size was 58um (on 512 pixels)

xticks(28.5+[-30:10:30])
xticklabels({'-3','','','0','','','3'})
xlim([28.5-30,28.5+30])
yticks([0,0.5,1])
yticklabels({'0%','50%','100%'})
ylim([-0.2,1.2])
xlabel({'Lateral position (um)'})
ylabel({'Normalized intensity'})

savefig(F,[save_root_fig,'\Fig4_PSF_xyResolutionCurve.fig']);
saveas(F,[save_root_png,'\Fig4_PSF_xyResolutionCurve.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_PSF_xyResolutionCurve.pdf']); set(gcf,'Color',[1,1,1])


%% Fig4_PSF_zResolutionCurve

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

temp = ~isnan(normPsf_z_all);
plot(this_y_all(temp),normPsf_z_all(temp),'Color',p.col.darkGray);
hold on; temp=plot(normPsf_z_fit); temp.Color = p.col.photostim; legend off;
% peak at 169.5, each step are 0.2 um -> 50 steps are 10 um, 75 steps are 15um

xticks(169.5+[-150:75:150])
xticklabels({'-30','-15','0','15','30'})
xlim([19.5,319.5])
yticks([0,0.5,1])
yticklabels({'0%','50%','100%'})
ylim([-0.2,1.2])
ylabel({'Normalized intensity'})
xlabel({'Axial position (um)'})

view([90 -90]);
savefig(F,[save_root_fig,'\Fig4_PSF_zResolutionCurve.fig']);
saveas(F,[save_root_png,'\Fig4_PSF_zResolutionCurve.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig4_PSF_zResolutionCurve.pdf']); set(gcf,'Color',[1,1,1])



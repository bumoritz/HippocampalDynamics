function F = propertyMap(s2p_meta,iscell,property,issubplot,plt)

img = zeros(s2p_meta.ops.Ly,s2p_meta.ops.Lx);
for n=1:length(iscell)
    if iscell(n)
        for i=1:length(s2p_meta.stat{n}.xpix)
            img(double(s2p_meta.stat{n}.ypix(i)),double(s2p_meta.stat{n}.xpix(i))) = property(n);
        end
    end
end

meds = zeros(length(s2p_meta.stat),2);
for i = 1:length(s2p_meta.stat)
    meds(i,:) = [double(s2p_meta.stat{i}.med(2)),double(s2p_meta.stat{i}.med(1))];
end

if nargin < 5
    c_abs = 1;
    c_label = 'z-score';
    titleText = ['Property map'];
else
    c_abs = plt.c_abs;
    c_label = plt.c_label;
    titleText = plt.titleText;
end

if ~issubplot
    F = default_figure([20,0.5,20,9.9]);
else
    F = NaN;
end

imagesc(img);
set(gca,'CLim',[-c_abs,c_abs]);
colormap(redblue);
temp=colorbar; 
temp.Label.String = c_label;
daspect([1,1,1]);
hold on

xlim([0,s2p_meta.ops.Lx])
ylim([0,s2p_meta.ops.Ly])
title(titleText)

set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);

if ~issubplot
    drawnow;
end
end
% file information
ver1 = 'v2.3'; % version
ver2 = 'v1.0-HR';
var1 = 'O2'; % var11iable
path1 = ['/raid/Data/GOBAI-' var1 '/' ver1 '/']; % file path
path2 = ['/raid/Data/GOBAI-' var1 '/' ver2 '/']; % file path
var2 = 'NO3'; % var11iable
path3 = ['/raid/Data/GOBAI-' var2 '/' ver2 '/']; % file path

% establish figure
figure('Position',[100 100 1200 800]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

% time and depth
pres = 100;
time = datenum(2017,6,15);
o2_lims = [0 300];

% download dimensions
gobai1.lon = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'lon'));
gobai1.lon = convert_lon(gobai1.lon,'format','0-360');
gobai1.lat = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'lat'));
gobai1.pres = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'pres'));
gobai1.time = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'time'));
gobai2.lon = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'lon'));
gobai2.lat = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'lat'));
gobai2.pres = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'pres'));
gobai2.time = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'time'));
gobai3.lon = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'lon'));
gobai3.lat = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'lat'));
gobai3.pres = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'pres'));
gobai3.time = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'time'));

% index
[~,idx_pres_1] = min(abs(pres-gobai1.pres));
[~,idx_pres_2] = min(abs(pres-gobai2.pres));
[~,idx_pres_3] = min(abs(pres-gobai3.pres));
[~,idx_time_1] = min(abs(time-(datenum(1950,1,1)+gobai1.time)));
[~,idx_time_2] = min(abs(time-(datenum(1950,1,1)+gobai2.time)));
[~,idx_time_3] = min(abs(time-(datenum(1950,1,1)+gobai3.time))); 

% plot v2.3 map (O2)
axis1 = nexttile;
worldmap([-90 90],[20 380]); setm(gca,'FontSize',12);
set(findobj(axis1.Children,'Tag','MLabel'),'FontSize',6);
set(findobj(axis1.Children,'Tag','PLabel'),'FontSize',6);
gobai1.oxy = squeeze(double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'oxy',...
    [1 1 idx_pres_1 1],[Inf Inf 1 Inf])));
pcolorm(gobai1.lat-0.25,gobai1.lon-0.25,mean(gobai1.oxy,3,'omitnan')');
plot_land('map');
c=colorbar(axis1);
c.Label.String = '[O_{2}] (\mumol kg^{-1})';
c.Label.FontSize = 12;
colormap(axis1,cmocean('ice'));
clim(axis1,o2_lims);
title(axis1,['GOBAI-O_{2}-' ver1 ' at ' num2str(gobai1.pres(idx_pres_1)) ' dbar']);
clear gobai1

% plot v1.0-HR map (O2)
axis2 = nexttile;
worldmap([-90 90],[20 380]); setm(gca,'FontSize',12);
set(findobj(axis2.Children,'Tag','MLabel'),'FontSize',6);
set(findobj(axis2.Children,'Tag','PLabel'),'FontSize',6);
gobai2.oxy = squeeze(double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'o2',...
    [1 1  idx_pres_2 1],[Inf Inf 1 Inf])));
pcolorm(gobai2.lat-0.25,gobai2.lon-0.25,mean(gobai2.oxy,3,'omitnan')');
plot_land('map');
c=colorbar(axis2);
c.Label.String = '[O_{2}] (\mumol kg^{-1})';
c.Label.FontSize = 12;
colormap(axis2,cmocean('ice'));
clim(axis2,o2_lims);
title(axis2,['GOBAI-O_{2}-' ver2 ' at ' num2str(gobai2.pres(idx_pres_2)) ' dbar']);
clear gobai2

% plot v1.0-HR map (NO3)
axis3 = nexttile;
worldmap([-90 90],[20 380]); setm(gca,'FontSize',12);
set(findobj(axis3.Children,'Tag','MLabel'),'FontSize',6);
set(findobj(axis3.Children,'Tag','PLabel'),'FontSize',6);
gobai3.nit = squeeze(double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'no3',...
    [1 1  idx_pres_2 1],[Inf Inf 1 Inf])));
pcolorm(gobai3.lat-0.25,gobai3.lon-0.25,mean(gobai3.nit,3,'omitnan')');
plot_land('map');
c=colorbar(axis3);
c.Label.String = '[NO_{3}] (\mumol kg^{-1})';
c.Label.FontSize = 12;
colormap(axis3,cmocean('ice'));
clim(axis3,o2_lims);
title(axis3,['GOBAI-NO_{3}-' ver2 ' at ' num2str(gobai3.pres(idx_pres_2)) ' dbar']);
clear gobai3


% export_fig(gcf,'Paper_Figs/Fig2.png','-transparent');
exportgraphics(gcf,'Paper_Figs/Fig2.png');

%% plot O2 observational data in 3D
% load data
load('Data/processed_all_o2_data_Feb-2024_D.mat')
%load('Data/processed_all_no3_data_Feb-2024_D.mat')
%load('O2/Data/processed_float_o2_data_Feb-2024_D.mat')

% create plot
figure('visible','on');
scatter3(all_data.longitude,all_data.latitude,-all_data.pressure,5,all_data.oxygen,'.');
%scatter3(all_data.longitude,all_data.latitude,-all_data.pressure,5,all_data.nitrate,'.');
%scatter3(float_data.LON,float_data.LAT,-float_data.PRES,5,float_data.OXY,'.');
hold on; set(gcf,'Position',[20 400 1000 600]);
ax = gca; ax.FontSize = 16;
ax.CameraPosition = [-269.2745 -351.3955 3.9736e+03];
% plot land
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
% modify axis limits
zlim([-500 0]);
xlim([-180 180]);
ylim([-90 90]);
clim([0 400]);
%clim([0 40]);
% edit colorbar
c=colorbar;
c.Label.String = '[O_{2}] (\mumol kg^{-1})';
%c.Label.String = '[NO_{3}] (\mumol kg^{-1})';
zlabel('Depth (dbar)')
colormap(cmocean('ice'));
%colormap(cmocean('solar'));
zticks([-500 -250 0]);
% edit background
set(gca,'color','none');
% export
exportgraphics(gcf,'Figures/3d_obs_oxygen.png','BackgroundColor','none')
%exportgraphics(gcf,'Figures/3d_obs_nitrate.png','BackgroundColor','none')

%% plot NO3 observational data in 3D
% load data
%load('Data/processed_all_o2_data_Feb-2024_D.mat')
load('NO3/Data/processed_float_no3_data_Feb-2024_D.mat')
% create plot
figure('visible','on');
%scatter3(all_data.longitude,all_data.latitude,-all_data.pressure,1,all_data.oxygen);
scatter3(float_data.LON,float_data.LAT,-float_data.PRES,5,float_data.NIT,'.');
hold on; set(gcf,'Position',[20 400 1000 600]);
ax = gca; ax.FontSize = 16;
ax.CameraPosition = [-269.2745 -351.3955 3.9736e+03];
% plot land
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
% modify axis limits
zlim([-500 0]);
xlim([-180 180]);
ylim([-90 90]);
clim([0 40]);
% edit colorbar
c=colorbar;
c.Label.String = '[NO_{3}] (\mumol kg^{-1})';
zlabel('Depth (dbar)')
colormap(cmocean('solar'));
zticks([-500 -250 0]);
% edit background
set(gca,'color','none');
% export
exportgraphics(gcf,'Figures/3d_obs_nitrate.gif','BackgroundColor','none')

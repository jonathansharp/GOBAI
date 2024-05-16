%% Plot mean [O2] with depth over time

% initialize figure
figure;
set(gcf,'units','inches','position',[0 14 20 10]); hold on;
set(gca,'fontsize',26,'YDir','reverse','linewidth',3,'layer','top','box','on');
title('GOBAI-O_{2} Global Mean [O_{2}] Anomaly Over Time (\mumol kg^{-1})','fontsize',32);

% file information
ver = 'v2.2'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate 2d mean
gobai_oxy_3d = reshape(GOBAI.oxy,[length(GOBAI.lon)*length(GOBAI.lat) length(GOBAI.pres) length(GOBAI.time)]);
gobai_vol_2d = reshape(GOBAI.vol,[length(GOBAI.lon)*length(GOBAI.lat) length(GOBAI.pres)]);
global_mean_2d = sum(gobai_oxy_3d.*gobai_vol_2d,'omitnan')./sum(gobai_vol_2d,'omitnan');
global_mean_2d = squeeze(global_mean_2d);

% plot data
if strcmp(ver,'v1.0') || strcmp(ver,'v2.0')
    time_axis = double(GOBAI.time);
elseif strcmp(ver,'v2.1')
    time_axis = datenum(2004,0,0) + double(GOBAI.time);
elseif strcmp(ver,'v2.2')
    time_axis = datenum(1950,0,0) + double(GOBAI.time);
end
contourf(time_axis,double(GOBAI.pres),...
    double(global_mean_2d-mean(global_mean_2d,2)),-3.75:0.5:3.75);

% modify axes
xlim([min(time_axis)-15 max(time_axis)+15]);
datetick('x','keeplimits');
% ylim([0 1500]);
ylabel('Pressure (dbars)');
clim([-3.75 3.75]);

% colorbar
c=colorbar;
c.TickLength = 0;
c.LineWidth = 3;
c.Position(3) = 2*c.Position(3);
c.Position(1) = c.Position(1)+0.04;
colormap(cmocean('balance',15,'pivot',0));

% save figure
exportgraphics(gcf,['Figures/global_O2_hovmoeller_GOBAI_' ver '.png']);
close
clear

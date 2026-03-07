%% Plot mean [O2] uncertainty with depth over time

%% download data
% file information
ver = 'v2.3'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen and uncertainty
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
GOBAI.uncer_D = ncread([path 'gobai-o2-uncer.nc'],'u_alg_o2');
GOBAI.uncer_DA = ncread([path 'gobai-o2-uncer-alg-new-DA.nc'],'u_alg_o2');
% GOBAI.uncer_GLODAP = ncread([path 'gobai-o2-uncer-alg-new.nc'],'u_alg_o2');
% calculate weights
[GOBAI.vol,~,GOBAI.h] = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.uncer_D,4,'omitnan'))) = NaN;
% calculate 2d mean
gobai_oxy_3d = reshape(GOBAI.oxy,[length(GOBAI.lon)*length(GOBAI.lat) length(GOBAI.pres) length(GOBAI.time)]);
gobai_uncer_D_3d = reshape(GOBAI.uncer_D,[length(GOBAI.lon)*length(GOBAI.lat) length(GOBAI.pres) length(GOBAI.time)]);
gobai_uncer_DA_3d = reshape(GOBAI.uncer_DA,[length(GOBAI.lon)*length(GOBAI.lat) length(GOBAI.pres) length(GOBAI.time)]);
gobai_vol_2d = reshape(GOBAI.vol,[length(GOBAI.lon)*length(GOBAI.lat) length(GOBAI.pres)]);
global_mean_2d = sum(gobai_oxy_3d.*gobai_vol_2d,'omitnan')./sum(gobai_vol_2d,'omitnan');
global_mean_2d = squeeze(global_mean_2d);
global_mean_uncer_D_2d = sum(gobai_uncer_D_3d.*gobai_vol_2d,'omitnan')./sum(gobai_vol_2d,'omitnan');
global_mean_uncer_D_2d = squeeze(global_mean_uncer_D_2d);
global_mean_uncer_DA_2d = sum(gobai_uncer_DA_3d.*gobai_vol_2d,'omitnan')./sum(gobai_vol_2d,'omitnan');
global_mean_uncer_DA_2d = squeeze(global_mean_uncer_DA_2d);

% process time axis
if strcmp(ver,'v1.0') || strcmp(ver,'v2.0')
    time_axis = double(GOBAI.time);
elseif strcmp(ver,'v2.1')
    time_axis = datenum(2004,0,0) + double(GOBAI.time);
elseif strcmp(ver,'v2.2')
    time_axis = datenum(1950,0,0) + double(GOBAI.time);
elseif strcmp(ver,'v2.3')
    time_axis = datenum(1950,0,0) + double(GOBAI.time);
end

%% plot absolute uncertainty
% initialize figure
figure;
set(gcf,'units','inches','position',[0 0 20 10]); hold on;
set(gca,'fontsize',26,'YDir','reverse','linewidth',3,...
    'layer','top','box','on','CLim',[5 12]);
title('GOBAI-O_{2} Global Mean [O_{2}] Algorithm Uncertainty Over Time (\mumol kg^{-1})','fontsize',24);

% plot data
contourf(time_axis,double(GOBAI.pres),...
    double(global_mean_uncer_D_2d),0:0.5:9);

% modify axes
xlim([min(time_axis)-15 max(time_axis)+15]);
xticks(datenum(2004:2:2026,1,1))
xticklabels({'2004' '2006' '2008' '2010' '2012' '2014' '2016' '2018' ...
    '2020' '2022' '2024' '2026'});
% datetick('x','yyyy','keeplimits');
% ylim([0 1500]);
ylabel('Pressure (dbars)');
clim([0 9]);

% colorbar
c=colorbar;
c.TickLength = 0;
c.LineWidth = 3;
c.Position(3) = 2*c.Position(3);
c.Position(1) = c.Position(1)+0.04;
colormap(cmocean('amp',14));

% save figure
exportgraphics(gcf,['Figures/global_O2_uncertainty_alg_hovmoeller_GOBAI_' ver '.png']);
close

%% plot uncertainty anomaly
% initialize figure
figure;
set(gcf,'units','inches','position',[0 0 20 10]); hold on;
set(gca,'fontsize',26,'YDir','reverse','linewidth',3,'layer','top','box','on','CLim',[-0.95 0.95]);
title('GOBAI-O_{2} Global Mean [O_{2}] Algorithm Uncertainty Anomaly Over Time (\mumol kg^{-1})','fontsize',24);

% plot data
contourf(time_axis,double(GOBAI.pres),...
    double(global_mean_uncer_D_2d-mean(global_mean_uncer_D_2d,2)),-0.95:0.1:0.95);

% modify axes
xlim([min(time_axis)-15 max(time_axis)+15]);
xticks(datenum(2004:2:2026,1,1))
xticklabels({'2004' '2006' '2008' '2010' '2012' '2014' '2016' '2018' ...
    '2020' '2022' '2024' '2026'});
% datetick('x','yyyy','keeplimits');
% ylim([0 1500]);
ylabel('Pressure (dbars)');
clim([-0.95 0.95]);

% colorbar
c=colorbar;
c.TickLength = 0;
c.LineWidth = 3;
c.Position(3) = 2*c.Position(3);
c.Position(1) = c.Position(1)+0.04;
colormap(cmocean('balance',19,'pivot',0));

% save figure
exportgraphics(gcf,['Figures/global_O2_uncertainty_alg_anom_hovmoeller_GOBAI_' ver '.png']);
close

%% Plot global mean [O2] uncertainty and number of floats over time
% initialize line plot
figure;
set(gcf,'units','inches','position',[0 0 20 3]); hold on;
set(gca,'fontsize',22,'linewidth',3,'layer','top','box','on');

% % determine active floats in each month
% load('O2/Data/processed_float_o2_data_adjusted_Oct-2025_D.mat')
% time_bins = [time_axis-14;time_axis(end)+16];
% active_floats =nan(size(time_axis));
% for t = 1:length(time_axis)
%     t_idx = float_data_adjusted.TIME >= time_bins(t) & ...
%         float_data_adjusted.TIME < time_bins(t+1);
%     active_floats(t) = length(unique(float_data_adjusted.FLOAT(t_idx)));
% end

% determine active floats in each month
load('O2/Data/processed_float_o2_data_Nov-2025_D_A.mat')
time_bins = [time_axis-14;time_axis(end)+16];
active_floats =nan(size(time_axis));
for t = 1:length(time_axis)
    t_idx = float_data.TIME >= time_bins(t) & ...
        float_data.TIME < time_bins(t+1);
    active_floats(t) = length(unique(float_data.FLOAT(t_idx)));
end

% plot data
global_mean_uncer_DA = sum(global_mean_uncer_DA_2d.*...
    squeeze(GOBAI.h(1,1,:)))./sum(squeeze(GOBAI.h(1,1,:)));
plot(time_axis,global_mean_uncer_DA,'k','linewidth',3);
ylabel('{\itu}[O_{2}]_{alg.} (\mumol kg^{-1})','FontSize',21);
yyaxis right
plot(time_axis,active_floats,'linewidth',3);
ylabel('Active Floats');

% determine active floats in each month
load('O2/Data/processed_float_o2_data_Oct-2025_D.mat')
time_bins = [time_axis-14;time_axis(end)+16];
active_floats =nan(size(time_axis));
for t = 1:length(time_axis)
    t_idx = float_data.TIME >= time_bins(t) & ...
        float_data.TIME < time_bins(t+1);
    active_floats(t) = length(unique(float_data.FLOAT(t_idx)));
end
plot(time_axis,active_floats,'linewidth',3);

% modify axes
xlim([min(time_axis)-15 max(time_axis)+15]);
xticks(datenum(2004:2:2026,1,1))
xticklabels({'2004' '2006' '2008' '2010' '2012' '2014' '2016' '2018' ...
    '2020' '2022' '2024' '2026'});

legend({'D & A Floats' 'D Floats Only'},'Location','north');

% save figure
exportgraphics(gcf,['Figures/global_O2_uncertainty_alg_wfloats_D_vs_DA_GOBAI_' ver '.png']);
close

%% Plot global mean [O2] uncertainty and number of cruises over time
% initialize line plot
figure;
set(gcf,'units','inches','position',[0 0 20 3]); hold on;
set(gca,'fontsize',22,'linewidth',3,'layer','top','box','on');

% determine cruises in each month
load('O2/Data/processed_glodap_o2_data_2023.mat');
time_bins = [time_axis-14;time_axis(end)+16];
cruises =nan(size(time_axis));
for t = 1:length(time_axis)
    t_idx = glodap_data.TIME >= time_bins(t) & ...
        glodap_data.TIME < time_bins(t+1);
    cruises(t) = length(unique(glodap_data.CRU(t_idx)));
end

% plot data
global_mean_uncer_D = sum(global_mean_uncer_D_2d.*...
    squeeze(GOBAI.h(1,1,:)))./sum(squeeze(GOBAI.h(1,1,:)));
plot(time_axis,global_mean_uncer_D,'k','linewidth',3);
ylabel('{\itu}[O_{2}]_{alg.} (\mumol kg^{-1})','FontSize',21);
yyaxis right
scatter(time_axis,cruises,'o','filled');
ylabel('Cruisess');

% modify axes
xlim([min(time_axis)-15 max(time_axis)+15]);
xticks(datenum(2004:2:2026,1,1))
xticklabels({'2004' '2006' '2008' '2010' '2012' '2014' '2016' '2018' ...
    '2020' '2022' '2024' '2026'});

% save figure
exportgraphics(gcf,['Figures/global_O2_uncertainty_alg_wcruises_GOBAI_' ver '.png']);
close

%% plot difference between two uncertainties
% initialize figure
figure;
set(gcf,'units','inches','position',[0 0 20 10]); hold on;
set(gca,'fontsize',26,'YDir','reverse','linewidth',3,'layer','top','box','on','CLim',[-0.95 0.95]);
title('GOBAI-O_{2} Global Mean [O_{2}] Algorithm Uncertainty Anomaly Over Time (\mumol kg^{-1})','fontsize',24);

% plot data
contourf(time_axis,double(GOBAI.pres),...
    double(global_mean_uncer_D_2d-global_mean_uncer_DA_2d),-0.95:0.1:0.95);

% modify axes
xlim([min(time_axis)-15 max(time_axis)+15]);
xticks(datenum(2004:2:2026,1,1))
xticklabels({'2004' '2006' '2008' '2010' '2012' '2014' '2016' '2018' ...
    '2020' '2022' '2024' '2026'});
% datetick('x','yyyy','keeplimits');
% ylim([0 1500]);
ylabel('Pressure (dbars)');
clim([-0.95 0.95]);

% colorbar
c=colorbar;
c.TickLength = 0;
c.LineWidth = 3;
c.Position(3) = 2*c.Position(3);
c.Position(1) = c.Position(1)+0.04;
colormap(cmocean('balance',19,'pivot',0));

% save figure
exportgraphics(gcf,['Figures/global_O2_uncertainty_alg_anom_D_vs_DA_hovmoeller_GOBAI_' ver '.png']);
close

%% clean up
clear

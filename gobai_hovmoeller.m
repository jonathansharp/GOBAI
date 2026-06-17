%% Plot mean GOBAI variable with depth over time

% file information
ver = 'v1.1-HR'; % version
var = 'O2'; % variable
varname = 'o2';
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);

% process time axis
if strcmp(ver,'v1.0') || strcmp(ver,'v2.0')
    time_axis = double(GOBAI.time);
elseif strcmp(ver,'v2.1')
    time_axis = datenum(2004,0,0) + double(GOBAI.time);
elseif strcmp(ver,'v2.2') || strcmp(ver,'v2.3')
    time_axis = datenum(1950,0,0) + double(GOBAI.time);
elseif strcmp(ver,'v1.0-HR') || strcmp(ver,'v1.1-HR')
    time_axis = datenum(1950,1,1) + double(GOBAI.time);
end

% % load mask
% if strcmp(ver,'v1.0-HR')
%     load([path 'mask-v1.0.mat']);
% elseif
% end

% define mask
if strcmp(ver,'v1.0-HR') || strcmp(ver,'v1.1-HR')
    GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
    if ~isfile([path 'mask-' ver '.mat'])
        mask = true(size(GOBAI.vol));
        for t = 1:length(GOBAI.time)
            gobai_tmp = ncread([path 'GOBAI-' var '-' ver '.nc'],...
                varname,[1 1 1 t],[Inf Inf Inf 1]);
            mask(isnan(gobai_tmp)) = false;
        end
        save([path 'mask-' ver '.mat'],'mask');
    else
        load([path 'mask-' ver '.mat'],'mask');
    end
end

% pre-allocate
global_mean_2d = nan(length(GOBAI.pres),length(GOBAI.time));
global_inv_2d = nan(length(GOBAI.pres),length(GOBAI.time));
% determine global means over depth
if ~isfile([path '2d_global_mean-' ver '.mat'])
    for t = 1:length(time_axis)
        % download variable
        GOBAI.var = ncread([path 'GOBAI-' var '-' ver '.nc'],varname,...
            [1 1 1 t],[Inf Inf Inf 1]);
        % adjust volume
        volume = GOBAI.vol;
        if strcmp(ver,'v1.0-HR') || strcmp(ver,'v1.1-HR')
            volume(isnan(GOBAI.var) | ~mask) = NaN;
        else
            volume(isnan(GOBAI.var)) = NaN;
        end
        % calculate 2d mean
        gobai_var_2d = reshape(GOBAI.var,[length(GOBAI.lon)*length(GOBAI.lat) length(GOBAI.pres)]);
        gobai_vol_2d = reshape(volume,[length(GOBAI.lon)*length(GOBAI.lat) length(GOBAI.pres)]);
        global_mean_2d(:,t) = sum(gobai_var_2d.*gobai_vol_2d,'omitnan')./sum(gobai_vol_2d,'omitnan');
        % umol/kg * m3 * kg/m3 * mol/umol = Pmol
        global_inv_2d(:,t) = sum(gobai_var_2d.*(gobai_vol_2d.*1027),'omitnan')./(10^21); 
    end
    save([path '2d_global_mean-' ver '.mat'],'global_mean_2d');
else
    load([path 'mask-' ver '.mat'],'mask');
end

% limits
if strcmp(var,'O2')
    var_lims_mon = -3.75:0.5:3.75;
    var_lims_ann = -2.75:0.25:2.75;
elseif strcmp(var,'NO3')
    var_lims_mon = -0.5:0.05:0.5;
    var_lims_ann = -0.3:0.03:0.3;
end

%% plot with seasonal cycle
% initialize figure
figure;
set(gcf,'units','inches','position',[0 0 20 10]); hold on;
set(gca,'fontsize',26,'YDir','reverse','linewidth',3,'layer','top','box','on','CLim',[-3.75 3.75]);
if strcmp(var,'O2')
    title('GOBAI-O_{2} Global Mean [O_{2}] Anomaly Over Time (\mumol kg^{-1})','fontsize',24);
elseif strcmp(var,'NO3')
    title('GOBAI-NO_{3} Global Mean [NO_{3}] Anomaly Over Time (\mumol kg^{-1})','fontsize',24);
end
% plot data
contourf(time_axis,double(GOBAI.pres),...
    double(global_mean_2d-mean(global_mean_2d,2)),var_lims_mon);
% modify axes
xlim([min(time_axis)-15 max(time_axis)+15]);
date = datevec(time_axis); % get years
xticks(datenum(date(2,1):2:date(end,1),1,1));
xticklabels(compose('%d',(date(2,1):2:date(end,1))'));
% datetick('x','yyyy','keeplimits');
% ylim([0 1500]);
ylabel('Pressure (dbars)');
clim([min(var_lims_mon) max(var_lims_mon)]);

% colorbar
c=colorbar;
c.TickLength = 0;
c.LineWidth = 3;
c.Position(3) = 2*c.Position(3);
c.Position(1) = c.Position(1)+0.04;
colormap(cmocean('balance',15,'pivot',0));

% save figure
exportgraphics(gcf,['Figures/global_' var '_hovmoeller_GOBAI_' ver '.png']);
% export_fig(['Figures/global_' var '_hovmoeller_GOBAI_' ver '.png'],'-transparent');
close

%% plot anomaly line plot
% load pressure bounds
GOBAI.prs_bnds = ncread([path 'GOBAI-' var '-' ver '.nc'],'prs_bnds');
GOBAI.prs_height = GOBAI.prs_bnds(:,2)-GOBAI.prs_bnds(:,1);
global_mean_global = sum(global_mean_2d.*GOBAI.prs_height)./sum(GOBAI.prs_height);
[yf,yr,x,err] = leastsq2(time_axis,global_mean_global,time_axis(1),2,[365.2425 365.2425/2]);
% initialize figure
figure;
set(gcf,'units','inches','position',[0 0 10 4]); hold on;
plot(time_axis,global_mean_global,time_axis,yf)
% modify axes
xlim([min(time_axis)-15 max(time_axis)+15]);
date = datevec(time_axis); % get years
xticks(datenum(date(2,1):2:date(end,1),1,1));
xticklabels(compose('%d',(date(2,1):2:date(end,1))'));
% datetick('x','yyyy','keeplimits');
% ylim([0 1500]);
ylabel('[O_{2}] Anomaly');
% save figure
exportgraphics(gcf,['Figures/global_' var '_anomaly_GOBAI_' ver '.png']);
close

%% plot annual
% initialize figure
figure;
set(gcf,'units','inches','position',[0 0 20 10]); hold on;
set(gca,'fontsize',26,'YDir','reverse','linewidth',3,'layer','top','box','on','CLim',[-3.75 3.75]);
if strcmp(var,'O2')
    title('GOBAI-O_{2} Global Mean [O_{2}] Anomaly Over Time (\mumol kg^{-1})','fontsize',24);
elseif strcmp(var,'NO3')
    title('GOBAI-NO_{3} Global Mean [NO_{3}] Anomaly Over Time (\mumol kg^{-1})','fontsize',24);
end
% plot data
anom_ann = movmean(double(global_mean_2d-mean(global_mean_2d,2)),52,2);
contourf(time_axis,double(GOBAI.pres),anom_ann,var_lims_ann);

% modify axes
xlim([min(time_axis)-15 max(time_axis)+15]);
date = datevec(time_axis); % get years
xticks(datenum(date(2,1):2:date(end,1),1,1));
xticklabels(compose('%d',(date(2,1):2:date(end,1))'));
% datetick('x','yyyy','keeplimits');
% ylim([0 1500]);
ylabel('Pressure (dbars)');
clim([min(var_lims_ann) max(var_lims_ann)]);

% colorbar
c=colorbar;
c.TickLength = 0;
c.LineWidth = 3;
c.Position(3) = 2*c.Position(3);
c.Position(1) = c.Position(1)+0.04;
colormap(cmocean('balance',23,'pivot',0));

% save figure
exportgraphics(gcf,['Figures/global_' var '_hovmoeller_ann_GOBAI_' ver '.png']);
% export_fig(['Figures/global_' var '_hovmoeller_GOBAI_' ver '.png'],'-transparent');
close

%% plot de-seasonalized
% initialize figure
figure;
set(gcf,'units','inches','position',[0 0 20 10]); hold on;
set(gca,'fontsize',26,'YDir','reverse','linewidth',3,'layer','top','box','on','CLim',[-3.75 3.75]);
if strcmp(var,'O2')
    title('GOBAI-O_{2} Global Mean [O_{2}] Anomaly Over Time (\mumol kg^{-1})','fontsize',24);
elseif strcmp(var,'NO3')
    title('GOBAI-NO_{3} Global Mean [NO_{3}] Anomaly Over Time (\mumol kg^{-1})','fontsize',24);
end
% calculate trends
global_mean_2d_deseas = nan(size(global_mean_2d));
for d = 1:size(global_mean_2d,1)
    [yf,yr,x] = leastsq2(time_axis,global_mean_2d(d,:),time_axis(1),...
        2,[365.2425 365.2425/2]);
    global_mean_2d_deseas(d,:) = (x(1) + x(2).*(time_axis-time_axis(1)) + yr) - mean(yf);
end
% plot data
contourf(time_axis,double(GOBAI.pres),movmean(global_mean_2d_deseas,4,2),var_lims_mon);

% modify axes
xlim([min(time_axis)-15 max(time_axis)+15]);
date = datevec(time_axis); % get years
xticks(datenum(date(2,1):2:date(end,1),1,1));
xticklabels(compose('%d',(date(2,1):2:date(end,1))'));
% datetick('x','yyyy','keeplimits');
% ylim([0 1500]);
ylabel('Pressure (dbars)');
clim([min(var_lims_ann) max(var_lims_ann)]);

% colorbar
c=colorbar;
c.TickLength = 0;
c.LineWidth = 3;
c.Position(3) = 2*c.Position(3);
c.Position(1) = c.Position(1)+0.04;
colormap(cmocean('balance',23,'pivot',0));

% save figure
exportgraphics(gcf,['Figures/global_' var '_hovmoeller_deseas_smoothed_GOBAI_' ver '.png']);
% export_fig(['Figures/global_' var '_hovmoeller_GOBAI_' ver '.png'],'-transparent');
close

%% clean up
clear

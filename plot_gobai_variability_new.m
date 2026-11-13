% plot gobai variability over years
pressures = [10 100 200 400 1250];
longitudes = [64 185 310];

% gobai and rfrom data
ver = 'HR-v202608'; % product version
modifier = 'fg-yes_pH_corr-no_dic_adj';
path = ['/fast6/dic/GOBAI/RFROM/FFNN/c20_Aug-2026_D/train80_val10_test10/' modifier '/']; % file path
gobai_f = [path 'GOBAI-DIC-' ver '.nc']; % file name

% download dimensions
gobai.lon = ncread(gobai_f,'lon');
gobai.lat = ncread(gobai_f,'lat');
gobai.pres = ncread(gobai_f,'pres');
gobai.time = datenum(1950,1,1) + ncread(gobai_f,'time');

% calculate mass per grid cell
gobai.vol = weights3d(gobai.lon,gobai.lat,gobai.pres);

% determine mask (time varying coverage)
if ~isfile([path 'mask_' ver '.mat'])
    mask = true(size(gobai.vol));
    year = year(gobai.time);
    year_curr = min(year);
    for t = 1:length(gobai.time)
        dic = ncread(gobai_f,'dic',[1 1 1 t],[Inf Inf Inf 1]);
        mask(isnan(dic)) = false;
    if year(t) - year_curr > 0
        fprintf('Year change detected: %d -> %d\n', year_curr, year(t));
        year_curr = year(t);
    end
    end
    save([path 'mask_' ver '.mat'],'mask'); % save mask
else
    load([path 'mask_' ver '.mat'],'mask');
end

% plot gobai variability on depth levels
for z = 1:length(pressures)

    % load on the current pressure level
    pres_idx = find(gobai.pres == pressures(z));
    gobai_temp = squeeze(ncread(gobai_f,'dic',[1,1,pres_idx,1],[Inf Inf 1 Inf]));

    % reshape for matrix vectorization (space x time)
    gobai_mat = reshape(gobai_temp, [], length(gobai.time))';
    omega1 = 2 * pi / 365.2425; % Annual frequency (assuming time in days)
    omega2 = 2 * pi / (365.2425 / 2); % Semi-annual frequency
    % [intercept, trend, sin(p1), cos(p1), sin(p2), cos(p2)]
    x = [ones(size(gobai.time)), gobai.time, sin(omega1 * gobai.time), ...
        cos(omega1 * gobai.time), sin(omega2 * gobai.time), cos(omega2 * gobai.time)];
    beta = x \ gobai_mat;
    % compute modeled trend + seasonality
    gobai_modeled = x * beta;
    % calculate anomalies (difference between observed and modeled values)
    gobai_anoms = gobai_mat - gobai_modeled;
    % calculate total variability (standard deviation of the observed values over time)
    gobai_total_variability = std(gobai_mat,[],1,'omitnan'); % standard deviation over time
    gobai_total_variability_mean = mean(gobai_total_variability,2,'omitnan'); % mean
    gobai_total_variability = reshape(gobai_total_variability,length(gobai.lon),length(gobai.lat));
    % calculate seasonality (standard deviation of modeled seasonal component)
    gobai_seasonal = gobai_modeled - x(:,1:2) * beta(1:2,:); % modeled seasonal component
    gobai_seasonality = std(gobai_seasonal,[],1,'omitnan'); % standard deviation over time
    gobai_seasonality_mean = mean(gobai_seasonality,2,'omitnan'); % average
    gobai_seasonality = reshape(gobai_seasonality,length(gobai.lon),length(gobai.lat));
    % calculate residual variability (standard deviation of detrended anomalies)
    gobai_residual = gobai_anoms - (x(:,3:end) * beta(3:end,:)); % detrended anomalies
    gobai_residual_variability = std(gobai_residual,[],1,'omitnan'); % standard deviation over time
    gobai_residual_variability_mean = mean(gobai_residual_variability,2,'omitnan'); % mean
    gobai_residual_variability = reshape(gobai_residual_variability,length(gobai.lon),length(gobai.lat));
    % calculate trend (standard deviation of modeled trend component)
    gobai_trend = x(:,1:2) * beta(1:2,:); % modeled trend component
    gobai_trend = beta(2,:) * 365.2425 * 10; % umol/kg / dec
    gobai_trend_mean = mean(gobai_trend,2,'omitnan'); % standard deviation over time
    gobai_trend = reshape(gobai_trend,length(gobai.lon),length(gobai.lat));

    % print results
    fprintf('Pressure level %d trend: %f umol/kg/dec.\n', pressures(z), gobai_trend_mean);
    fprintf('Pressure level %d seasonality: %f umol/kg\n', pressures(z), gobai_seasonality_mean);
    fprintf('Pressure level %d residual variability: %f umol/kg\n', pressures(z), gobai_residual_variability_mean);
    fprintf('Pressure level %d total variability: %f umol/kg\n', pressures(z), gobai_total_variability_mean);

    % create plot for residual variability
    figure('color','w','visible','off','position',[100 100 800 500]);
    worldmap([-90 90],[20 380]);
    pcolorm(double(gobai.lat),double(gobai.lon),gobai_residual_variability');
    title(['GOBAI-DIC Residual Variability at ' num2str(pressures(z)) ' dbar'],'fontsize',14)
    c = colorbar(gca,'southoutside');
    c.Label.String = 'DIC Residual Variability (\mumol kg^{-1})';
    c.Label.FontSize = 12; 
    clim([0 60]);
    colormap(cmocean('dense'));
    plot_land('map');
    plabel off; mlabel off;
    export_fig(gcf,['DIC/Figures/gobai_residual_variability_' ver '_' num2str(pressures(z)) 'dbar.png']);
    close(gcf);

    % create plot for seasonality
    figure('color','w','visible','off','position',[100 100 800 500]);
    worldmap([-90 90],[20 380]);
    pcolorm(double(gobai.lat),double(gobai.lon),gobai_seasonality');
    title(['GOBAI-DIC Seasonality at ' num2str(pressures(z)) ' dbar'],'fontsize',14)
    c = colorbar(gca,'southoutside');
    c.Label.String = 'DIC Seasonality (\mumol kg^{-1})';
    c.Label.FontSize = 12;
    clim([0 60]);
    colormap(cmocean('speed'));
    plot_land('map');
    plabel off; mlabel off;
    export_fig(gcf,['DIC/Figures/gobai_seasonality_' ver '_' num2str(pressures(z)) 'dbar.png']);
    close(gcf);

    % create plot for trend
    figure('color','w','visible','off','position',[100 100 800 500]);
    worldmap([-90 90],[20 380]);
    pcolorm(double(gobai.lat),double(gobai.lon),gobai_trend');
    title(['GOBAI-DIC Trend at ' num2str(pressures(z)) ' dbar'],'fontsize',14)
    c = colorbar(gca,'southoutside');
    c.Label.String = 'DIC Trend (\mumol kg^{-1} dec.^{-1})';
    c.Label.FontSize = 12;
    clim([-10 20]);
    colormap(cmocean('balance','pivot', 0));
    plot_land('map');
    plabel off; mlabel off;
    export_fig(gcf,['DIC/Figures/gobai_trend_' ver '_' num2str(pressures(z)) 'dbar.png']);
    close(gcf);

    % create plot for total variability
    figure('color','w','visible','off','position',[100 100 800 500]);
    worldmap([-90 90],[20 380]);
    pcolorm(double(gobai.lat),double(gobai.lon),gobai_total_variability');
    title(['GOBAI-DIC Total Variability at ' num2str(pressures(z)) ' dbar'],'fontsize',14)
    c = colorbar(gca,'southoutside');
    c.Label.String = 'DIC Total Variability (\mumol kg^{-1})';
    c.Label.FontSize = 12;
    clim([0 60]);
    colormap(cmocean('matter'));
    plot_land('map');
    plabel off; mlabel off;
    export_fig(gcf,['DIC/Figures/gobai_total_variability_' ver '_' num2str(pressures(z)) 'dbar.png']);
    close(gcf);

    clearvars -except gobai gobai_f ver pressures z longitudes mask modifier

end

% plot gobai variability on longitudinal transects
for x = 1:length(longitudes)

    % load the nearest longitude
    lon_idx = find(min(abs(gobai.lon-longitudes(x))) == abs(gobai.lon-longitudes(x)),1);
    gobai_temp = squeeze(ncread(gobai_f,'dic',[lon_idx,1,1,1],[1 Inf Inf Inf]));
    transect_mask = squeeze(mask(lon_idx,:,:));

    % reshape for matrix vectorization (latitude-pressure x time)
    gobai_mat = reshape(gobai_temp, [], length(gobai.time))';
    omega1 = 2 * pi / 365.2425; % Annual frequency (assuming time in days)
    omega2 = 2 * pi / (365.2425 / 2); % Semi-annual frequency
    % [intercept, trend, sin(p1), cos(p1), sin(p2), cos(p2)]
    design = [ones(size(gobai.time)), gobai.time, sin(omega1 * gobai.time), ...
        cos(omega1 * gobai.time), sin(omega2 * gobai.time), cos(omega2 * gobai.time)];
    beta = design \ gobai_mat;
    % compute modeled trend + seasonality
    gobai_modeled = design * beta;
    % calculate anomalies (difference between observed and modeled values)
    gobai_anoms = gobai_mat - gobai_modeled;
    % calculate total variability (standard deviation of the observed values over time)
    gobai_total_variability = std(gobai_mat,[],1,'omitnan');
    gobai_total_variability(~transect_mask(:)') = NaN;
    gobai_total_variability_mean = mean(gobai_total_variability,2,'omitnan');
    gobai_total_variability = reshape(gobai_total_variability,length(gobai.lat),length(gobai.pres));
    % calculate seasonality (standard deviation of modeled seasonal component)
    gobai_seasonal = gobai_modeled - design(:,1:2) * beta(1:2,:);
    gobai_seasonality = std(gobai_seasonal,[],1,'omitnan');
    gobai_seasonality(~transect_mask(:)') = NaN;
    gobai_seasonality_mean = mean(gobai_seasonality,2,'omitnan');
    gobai_seasonality = reshape(gobai_seasonality,length(gobai.lat),length(gobai.pres));
    % calculate residual variability (standard deviation of detrended anomalies)
    gobai_residual = gobai_anoms - (design(:,3:end) * beta(3:end,:));
    gobai_residual_variability = std(gobai_residual,[],1,'omitnan');
    gobai_residual_variability(~transect_mask(:)') = NaN;
    gobai_residual_variability_mean = mean(gobai_residual_variability,2,'omitnan');
    gobai_residual_variability = reshape(gobai_residual_variability,length(gobai.lat),length(gobai.pres));
    % calculate trend
    gobai_trend = beta(2,:) * 365.2425 * 10; % umol/kg/decade
    gobai_trend(~transect_mask(:)') = NaN;
    gobai_trend_mean = mean(gobai_trend,2,'omitnan');
    gobai_trend = reshape(gobai_trend,length(gobai.lat),length(gobai.pres));

    % print results
    fprintf('Longitude %d trend: %f umol/kg/dec.\n', longitudes(x), gobai_trend_mean);
    fprintf('Longitude %d seasonality: %f umol/kg\n', longitudes(x), gobai_seasonality_mean);
    fprintf('Longitude %d residual variability: %f umol/kg\n', longitudes(x), gobai_residual_variability_mean);
    fprintf('Longitude %d total variability: %f umol/kg\n', longitudes(x), gobai_total_variability_mean);

    % create plot for residual variability
    figure('color','w','visible','off','position',[100 100 800 500]);
    pcolor(double(gobai.lat),double(gobai.pres),gobai_residual_variability');
    shading flat; set(gca,'YDir','reverse');
    title(['GOBAI-DIC Residual Variability at ' num2str(longitudes(x)) ' deg'],'fontsize',14);
    c = colorbar(gca,'southoutside');
    c.Label.String = 'DIC Residual Variability (\mumol kg^{-1})';
    c.Label.FontSize = 12;
    clim([0 60]);
    colormap(cmocean('dense'));
    export_fig(gcf,['DIC/Figures/gobai_residual_variability_' ver '_' num2str(longitudes(x)) 'deg.png']);
    close(gcf);

    % create plot for seasonality
    figure('color','w','visible','off','position',[100 100 800 500]);
    pcolor(double(gobai.lat),double(gobai.pres),gobai_seasonality');
    shading flat; set(gca,'YDir','reverse');
    title(['GOBAI-DIC Seasonality at ' num2str(longitudes(x)) ' deg'],'fontsize',14);
    c = colorbar(gca,'southoutside');
    c.Label.String = 'DIC Seasonality (\mumol kg^{-1})';
    c.Label.FontSize = 12;
    clim([0 60]);
    colormap(cmocean('speed'));
    export_fig(gcf,['DIC/Figures/gobai_seasonality_' ver '_' num2str(longitudes(x)) 'deg.png']);
    close(gcf);

    % create plot for trend
    figure('color','w','visible','off','position',[100 100 800 500]);
    pcolor(double(gobai.lat),double(gobai.pres),gobai_trend');
    shading flat; set(gca,'YDir','reverse');
    title(['GOBAI-DIC Trend at ' num2str(longitudes(x)) ' deg'],'fontsize',14);
    c = colorbar(gca,'southoutside');
    c.Label.String = 'DIC Trend (\mumol kg^{-1} dec.^{-1})';
    c.Label.FontSize = 12;
    clim([-10 20]);
    colormap(cmocean('balance','pivot', 0));
    export_fig(gcf,['DIC/Figures/gobai_trend_' ver '_' num2str(longitudes(x)) 'deg.png']);
    close(gcf);

    % create plot for total variability
    figure('color','w','visible','off','position',[100 100 800 500]);
    pcolor(double(gobai.lat),double(gobai.pres),gobai_total_variability');
    shading flat; set(gca,'YDir','reverse');
    title(['GOBAI-DIC Total Variability at ' num2str(longitudes(x)) ' deg'],'fontsize',14);
    c = colorbar(gca,'southoutside');
    c.Label.String = 'DIC Total Variability (\mumol kg^{-1})';
    c.Label.FontSize = 12;
    clim([0 60]);
    colormap(cmocean('matter'));
    export_fig(gcf,['DIC/Figures/gobai_total_variability_' ver '_' num2str(longitudes(x)) 'deg.png']);
    close(gcf);
    clear gobai_temp gobai_mat design beta gobai_modeled gobai_anoms
end
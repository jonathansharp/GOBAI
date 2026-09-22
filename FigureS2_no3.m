% file information
path_start = '/raid/Data/';
% GOBAI properties
file_date = '202606';
ver_rfrom = 'v2.2-2025';
ver_gobai = 'HR-v1.3';
clrs = jet(20);

% establish figure
figure('Position',[100 100 1300 800],'Visible','on');
height = 0.44; width = 0.46;
pos = [0.03 0.54 width height;...
    0.51 0.54 width height;...
    0.03 0.1 width height;...
    0.51 0.1 width height];

% time and depth
time = datenum(2024,1,15);
press = [10;200;400;1250];

% load GMM model
load('NO3/Models/GMM/c20_Jun-2026_D_A/model_fgo.mat');

for d = 1:4

    % RFROM path
    fpath_temp = [path_start 'RFROM/RFROM_TEMP_' ver_rfrom '/RFROMV22_TEMP_STABLE_' ...
        num2str(year(time)) '_' sprintf('%02d',month(time)) '.nc'];
    fpath_sal = [path_start 'RFROM/RFROM_SAL_' ver_rfrom '/RFROMV22_SAL_STABLE_' ...
        num2str(year(time)) '_' sprintf('%02d',month(time)) '.nc'];

    % download dimensions
    ds.lon = double(ncread(fpath_temp,'longitude'));
    ds.lon = convert_lon(ds.lon,'format','0-360');
    ds.lat = double(ncread(fpath_temp,'latitude'));
    ds.pres = double(ncread(fpath_temp,'mean_pressure'));
    ds.time = double(ncread(fpath_temp,'time'));

    % index
    [~,idx_pres] = min(abs(press(d)-ds.pres));
    [~,idx_time] = min(abs(time-(datenum(1950,1,1)+ds.time)));

    %
    t = squeeze(double(ncread(fpath_temp,'ocean_temperature',...
        [1 1 idx_pres idx_time],[Inf Inf 1 1])));
    s = squeeze(double(ncread(fpath_sal,'ocean_salinity',...
        [1 1 idx_pres idx_time],[Inf Inf 1 1])));
    pres_2d = repmat(press(d),size(t));
    idx = ~isnan(t) & ~isnan(s);
    X_norm = normalize([t(idx) s(idx) pres_2d(idx)],'Center',C,'Scale',S);

    % assign to clusters and obtain probabilities
    [clusters,~,p] = cluster(gmm,X_norm);
    clusters_2d = nan(size(t));
    clusters_2d(idx) = clusters;

    % plot map
    ax.(['ax' num2str(d)]) = axes('Position', pos(d,:));
    worldmap([-90 90],[20 380]);
    set(ax.(['ax' num2str(d)]),'FontSize',12);
    mlabel off; plabel off;
    pcolorm(ds.lat-0.25,ds.lon-0.25,clusters_2d');
    plot_land('map');
    clim([0.5 20.5]); colormap(clrs);
    title([num2str(press(d)) ' dbar'],'FontSize',16);

end

% colorbar
ax.ax5 = axes('Position',[0.08 0 0.9 1],'Visible','off');
c=colorbar(ax.ax5,'Location','southoutside');
clim([0.5 20.5]); colormap(clrs);
c.TickLength = 0;
c.Ticks = 1:num_clusters;
c.FontSize = 12;
annotation('textbox', [0 0 1 0.15],'String', 'Cluster Number', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 16, 'FontWeight', 'bold');

% save figure
export_fig(gcf,'Figures/FigureS2_no3_cluster_maps.png','-transparent');
close;

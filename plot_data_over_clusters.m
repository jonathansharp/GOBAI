% plot_data_over_clusters
%
% DESCRIPTION:
% This function plots data points overlying
% clusters on different depth levels.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 6/19/2024

function plot_data_over_clusters(param,base_grid,file_date,float_file_ext,num_clusters,numWorkers_train)

%% process parameter name
param1 = param_name(param);

%% plot data points by cluster
% load combined data
load([param1 '/Data/processed_all_' param '_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');
% define pressure axis
pressures = sort(unique(all_data.pressure));
% open parallel pool
tic; parpool(numWorkers_train); fprintf('Pool initiation:'); toc;
% make plots
parfor p = 1:length(pressures)
    % plot data by cluster
    figure('visible','off','Position',[100 100 800 400]); hold on;
    idx = all_data.pressure == pressures(p);
    m_proj('robinson','lon',[20 380]);
    m_coast('patch',rgb('grey'));
    m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
    lon_temp = convert_lon(convert_lon(all_data.longitude));
    lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
    % load dimensions
    if strcmp(base_grid,'RG')
        % file path
        fpath = [pwd '/Data/RG_CLIM/'];
        % load
        Longitude = ncread([fpath 'RG_Climatology_Temp.nc'],'Longitude');
        Latitude = ncread([fpath 'RG_Climatology_Temp.nc'],'Latitude');
        Pressure = ncread([fpath 'RG_Climatology_Temp.nc'],'Pressure');
    elseif strcmp(base_grid,'RFROM')
        % file path
        fpath = [pwd '/Data/RFROM/'];
        % load
        Longitude = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'longitude');
        Latitude = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'latitude');
        Pressure = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'mean_pressure');
    end
    % process longitude
    idx_20 = Longitude<20;
    Longitude(idx_20) = Longitude(idx_20)+360;
    Longitude = [Longitude(~idx_20);Longitude(idx_20)];
    % depth index
    depth_idx = find(Pressure == pressures(p));
    % plot monthly background clusters (**** should be long-term mean, but currently monthly ****)
    GMM_clusters = load(['Data/GMM_' base_grid '_' num2str(num_clusters) ...
        '/m1_w1'],'GMM_clusters');
    z = [GMM_clusters.GMM_clusters(~idx_20,:,depth_idx);...
        GMM_clusters.GMM_clusters(idx_20,:,depth_idx)];
    m_pcolor(double(Longitude),double(Latitude),double(z)');
    % plot data by cluster
    clrs = cmocean('balance',7);
    p1 =m_scatter(lon_temp(idx & all_data.id > 10^9),...
        all_data.latitude(idx & all_data.id > 10^9),1,'.','markerfacecolor',clrs(2,:));
    p2 = m_scatter(lon_temp(idx & all_data.id < 10^9),...
        all_data.latitude(idx & all_data.id < 10^9),1,'.','markerfacecolor',clrs(6,:));
    [~,icons]=legend([p1 p2],{'Argo Float Profiles' 'GLODAP Cruise Profiles'},...
        'location','northoutside','Orientation','horizontal','fontsize',16);
    icons(3).Children.MarkerSize = 10;
    icons(4).Children.MarkerSize = 10;
    colormap([1,1,1;cbrewer('qual','Pastel2',num_clusters)]);
    clim([-0.5 num_clusters+0.5]);
    c=colorbar;
    c.Limits = [0.5 num_clusters+0.5];
    c.Label.String = 'Cluster';
    c.TickLength = 0;
    % save figure
    dname = [param1 '/Figures/Clusters/' base_grid '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    exportgraphics(gcf,[dname '/data_over_clusters_' num2str(pressures(p)) '.png']);
    close
end
% end parallel session
delete(gcp('nocreate'));

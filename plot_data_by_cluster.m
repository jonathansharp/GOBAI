% plot_data_by_cluster
%
% DESCRIPTION:
% This function plots data points according to their most fitting cluster.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 3/20/2024

function plot_data_by_cluster(param_props,base_grid,file_date,float_file_ext,num_clusters,numWorkers_predict)

%% plot data points by cluster
% load combined data
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    load([param_props.dir_name '/Data/processed_all_' param_props.file_name '_data_' file_date float_file_ext '.mat'],'all_data');
else
    load([param_props.dir_name '/Data/' base_grid '_' param_props.file_name '_data_' file_date float_file_ext '.mat'],'all_data');
end
% load cluster data
load([param_props.dir_name '/Data/all_data_clusters_' base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');
% define pressure axis
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    pressures = sort(unique(all_data.pressure));
else
    pressures = sort(unique(all_data.depth));
end
% open parallel pool
tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;
% make plots
parfor p = 1:length(pressures)
    % plot data by cluster
    figure('visible','off','Position',[100 100 800 400]); hold on;
    set(gca,'fontsize',12);
    % use depth for CMIP models
    if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
        idx = all_data.pressure == pressures(p);
    else
        idx = all_data.depth == pressures(p);
    end
    m_proj('robinson','lon',[20 380]);
    m_coast('patch',rgb('grey'));
    m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
    lon_temp = convert_lon(convert_lon(all_data.longitude));
    lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
    % use depth for CMIP models
    if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
        title(['Data by Cluster at ' num2str(pressures(p)) ' dbars (in situ data)'],'fontsize',16);
    else
        title(['Data by Cluster at ' num2str(pressures(p)) ' m (' base_grid ')'],'fontsize',16);
    end
    m_scatter(lon_temp(idx),all_data.latitude(idx),3,all_data_clusters.clusters(idx),'filled');
    mycolormap = [1,1,1;flipud(jet(num_clusters))]; % white then jet
    colormap(mycolormap); % white then jet
    clim([-0.5 num_clusters+0.5]);
    c=colorbar;
    c.Limits = [0.5 num_clusters+0.5];
    c.Label.String = 'Cluster';
    c.FontSize = 12;
    c.TickLength = 0;
    % save figure
    dname = [param_props.dir_name '/Figures/Clusters/' base_grid '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    export_fig(gcf,[dname '/clustered_data_' num2str(pressures(p)) '.png'],...
        '-transparent','-silent');
    close
    % plot data by cluster probability
    for clst = 1:num_clusters
        figure('visible','off','Position',[100 100 800 400]); hold on;
        set(gca,'fontsize',12);
        % use depth for CMIP models
        if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
            idx = all_data.pressure == pressures(p);
        else
            idx = all_data.depth == pressures(p);
        end
        m_proj('robinson','lon',[20 380]);
        m_coast('patch',rgb('grey'));
        m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
        lon_temp = convert_lon(convert_lon(all_data.longitude));
        lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
        % use depth for CMIP models
        if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
            title(['Data by Cluster at ' num2str(pressures(p)) ' dbars (in situ data)'],'fontsize',16);
        else
            title(['Data by Cluster at ' num2str(pressures(p)) ' m (' base_grid ')'],'fontsize',16);
        end
        m_scatter(lon_temp(idx),all_data.latitude(idx),3,...
            all_data_clusters.(['c' num2str(clst)])(idx),'filled');
        colormap(customcolormap([0 1],[mycolormap(clst+1,:); 1 1 1]));
        clim([0 1]);
        c=colorbar;
        c.Limits = [0 1];
        c.Label.String = ['Cluster #' num2str(clst) ' Probability'];
        c.FontSize = 12;
        c.TickLength = 0;
        % save figure
        dname = [param_props.dir_name '/Figures/Clusters/' base_grid '_c' num2str(num_clusters)];
        if ~isfolder([pwd '/' dname]); mkdir(dname); end
        export_fig([dname '/clustered_data_probability_c' ...
            num2str(clst) '_' num2str(pressures(p)) '.png'],...
            '-transparent','-silent');
        close
    end
end

% end parallel session
delete(gcp('nocreate'));

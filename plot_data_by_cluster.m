%% plot data points by cluster
% load combined data
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');
% load cluster data
load(['Data/all_data_clusters_' base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');
% define pressure axis
pressures = sort(unique(all_data.pressure));
% open parallel pool
parpool;
% make plots
parfor p = 1:length(pressures)
    h=figure('visible','off','Position',[100 100 800 400]); hold on;
    idx = all_data.pressure == pressures(p);
    m_proj('robinson','lon',[20 380]);
    m_coast('patch',rgb('grey'));
    m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
    lon_temp = convert_lon(convert_lon(all_data.longitude));
    title(['Data by Cluster (' num2str(pressures(p)) ' dbars)']);
    m_scatter(lon_temp(idx),all_data.latitude(idx),3,all_data_clusters.clusters(idx),'filled');
    colormap([1,1,1;flipud(jet(num_clusters))]); % white then jet
    clim([-0.5 num_clusters+0.5]);
    c=colorbar;
    c.Limits = [0.5 num_clusters+0.5];
    c.Label.String = 'Cluster';
    c.TickLength = 0;
    % save figure
    dname = ['Figures/Clusters/' base_grid '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    exportgraphics(gcf,[dname '/clustered_data_' num2str(pressures(p)) '.png']);
    close
end
% end parallel session
delete(gcp('nocreate'));

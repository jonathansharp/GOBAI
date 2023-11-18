%% plot data points by cluster
% load combined data
load_interpolated_combined_data_to_workspace
% load cluster data
load(['Data/all_data_clusters_'  num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');
pressures = sort(unique(all_data.pressure));
% make plots
for p = 1:length(pressures)
    h=figure('visible','on','Position',[100 100 800 400]);
    idx = all_data.pressure == 10;
    worldmap([-90 90],[20 380]);
    title(['Data by Cluster (' num2str(pressures(p)) ' dbars)']);
    scatterm(all_data.latitude(idx),all_data.longitude(idx),1,all_data_clusters.clusters(idx));
    colormap([1,1,1;flipud(jet(10))]); % white then jet
    plot_land('map');
    clim([-0.5 10.5]);
    c=colorbar;
    c.Limits = [0.5 10.5];
    c.Label.String = 'Cluster';
    c.TickLength = 0;
    mlabel off; plabel off;
    exportgraphics(gcf,['Figures/Clusters/clustered data_' num2str(pressures(p)) '.png']);
    close
    % clean up
    clear h idx c 
end
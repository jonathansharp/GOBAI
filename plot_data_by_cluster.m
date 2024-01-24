%% plot data points by cluster
% load combined data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');
% load cluster data
load(['Data/all_data_clusters_' base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');
% define pressure axis
pressures = sort(unique(all_data.pressure));
% open parallel pool
parpool
% make plots
parfor p = 1:length(pressures)
    h=figure('visible','off','Position',[100 100 800 400]);
    idx = all_data.pressure == pressures(p);
    worldmap([-90 90],[20 380]);
    title(['Data by Cluster (' num2str(pressures(p)) ' dbars)']);
    scatterm(all_data.latitude(idx),all_data.longitude(idx),1,all_data_clusters.clusters(idx));
    colormap([1,1,1;flipud(jet(num_clusters))]); % white then jet
    plot_land('map');
    clim([-0.5 num_clusters+0.5]);
    c=colorbar;
    c.Limits = [0.5 num_clusters+0.5];
    c.Label.String = 'Cluster';
    c.TickLength = 0;
    mlabel off; plabel off;
    % save figure
    dname = ['Figures/Clusters/' base_grid '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    exportgraphics(gcf,[dname '/clustered_data_' num2str(pressures(p)) '.png']);
    close
end
% end parallel session
delete(gcp('nocreate'));

%% Plot cluster averages
pressures = [2.5 100 500 1000 1975];
parfor d = 1:length(pressures)
    % establish figure
    h=figure('visible','off','Position',[100 100 800 400]);
    % create folder
    dname = ['Figures/Clusters/' base_grid '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    % establish file name
    fname = ['clusters_' num2str(pressures(d)) 'dbar.png'];
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
    % depth index
    depth_idx = find(Pressure == pressures(d));
    % plot clusters each month/week
    for m = 1:timesteps
        if strcmp(base_grid,'RG')
            % load monthly clusters
            GMM_clusters = load(['Data/GMM_' base_grid '_' num2str(num_clusters) ...
                '/t' num2str(m)],'GMM_clusters');
            % make plot
            worldmap([-90 90],[20 380]);
            title(extractAfter(datestr(datenum(2004,m,1)),'-'));
            pcolorm(double(Latitude),double(Longitude),...
                double(GMM_clusters.GMM_clusters(:,:,depth_idx))');
            colormap([1,1,1;flipud(jet(num_clusters))]); % white then jet
            plot_land('map');
            clim([-0.5 num_clusters+0.5]);
            c=colorbar;
            c.Limits = [0.5 num_clusters+0.5];
            c.Label.String = 'Cluster';
            c.TickLength = 0;
            mlabel off; plabel off;
            % capture frame
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % write to file
            if m == 1
                imwrite(imind,cm,[dname '/' fname],'gif','Loopcount',inf,'DelayTime',0.1);
            else
                imwrite(imind,cm,[dname '/' fname],'gif','WriteMode','append','DelayTime',0.1);
            end
        elseif strcmp(base_grid,'RFROM')
            % determine number of weeks in file
            weeks = length(dir(['Data/GMM_' base_grid '_' num2str(num_clusters) '/m' num2str(m) '_*.mat']));
            for w = 1:weeks
                % load monthly clusters
                GMM_clusters = load(['Data/GMM_' base_grid '_' num2str(num_clusters) ...
                    '/m' num2str(m) '_w' num2str(w)],'GMM_clusters');
                % make plot
                worldmap([-90 90],[20 380]);
                title(extractAfter(datestr(datenum(2004,m,1)),'-'));
                pcolorm(double(Latitude),double(Longitude),...
                    double(GMM_clusters.GMM_clusters(:,:,depth_idx))');
                colormap([1,1,1;flipud(jet(num_clusters))]); % white then jet
                plot_land('map');
                clim([-0.5 num_clusters+0.5]);
                c=colorbar;
                c.Limits = [0.5 num_clusters+0.5];
                c.Label.String = 'Cluster';
                c.TickLength = 0;
                mlabel off; plabel off;
                % capture frame
                frame = getframe(h);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                % write to file
                if m == 1
                    imwrite(imind,cm,[dname '/' fname],'gif','Loopcount',inf,'DelayTime',0.1);
                else
                    imwrite(imind,cm,[dname '/' fname],'gif','WriteMode','append','DelayTime',0.1);
                end
            end
        end
    end
    close
end

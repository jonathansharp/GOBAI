%% Plot cluster probabilities over time

function plot_probability_animation(base_grid,num_clusters,numWorkers_train)

% set pressures
pressures = [2.5 10 50 100 200 300 500 1000 1500 1975];

% set up parallel pool
% p = setup_pool(numWorkers_train);

for d = 2%1:length(pressures)
    for cl = 13%1:length(num_clusters)
        % create folder
        dname = ['Figures/Clusters/' base_grid '_c' num2str(num_clusters) '/c' num2str(cl)];
        if ~isfolder([pwd '/' dname]); mkdir(dname); end
        % establish file name
        fname = ['probability_animation_c' num2str(cl) '_' num2str(pressures(d)) 'dbar.gif'];
        % determine number of monthly timesteps
        ds = dir(['Data/GMM_' base_grid '_' num2str(num_clusters) '/*.mat']);
        timesteps = length(ds);
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
        depth_idx = find(Pressure == pressures(d));
        % set counter
        cnt = 1;
        % establish fiugre
        h = figure('color','w','visible','off');
        axis tight manual
        % plot clusters each month/week
        for m = 121%1:length(timesteps)
            if strcmp(base_grid,'RG')
                % clear frame
                clf
                % load monthly clusters
                GMM_cluster_probs = load(['Data/GMM_' base_grid '_' num2str(num_clusters) ...
                    '/c' num2str(num_clusters(cl)) '/m' num2str(m) '_w1'],'GMM_clusters');
                % establish figure
                figure(h);
                % make plot
                m_proj('robinson','lon',[20 380]);
                z = [GMM_cluster_probs.GMM_cluster_probs(~idx_20,:,depth_idx);...
                    GMM_cluster_probs.GMM_cluster_probs(idx_20,:,depth_idx)];
                m_pcolor(double(Longitude),double(Latitude),double(z)');
                title(gca,extractAfter(datestr(datenum(2004,m,1)),'-'));
                colormap(cmocean('amp'));
                m_coast('patch',rgb('grey'));
                m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
                clim([0 100]);
                c=colorbar;
                c.Limits = [0 100];
                c.Label.String = ['Cluster ' num2str(num_clusters(cl)) ' Probability'];
                c.TickLength = 0;
                % save frame
                if ~isfolder([dname '/' num2str(pressures(d)) 'dbars/probs'])
                    mkdir([dname '/' num2str(pressures(d)) 'dbars/probs']);
                end
                exportgraphics(h,[dname '/' num2str(pressures(d)) ...
                    'dbars/probs/m' num2str(m) '.png']);
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
                    % clear frame
                    clf
                    % load monthly clusters
                    GMM_cluster_probs = load(['Data/GMM_' base_grid '_' num2str(num_clusters) ...
                        '/c' num2str(cl) '/m' num2str(m) '_w' num2str(w)],'GMM_cluster_probs');
                    % establish figure
                    figure(h);
                    % make plot
                    m_proj('robinson','lon',[20 380]);
                    z = 100.*[GMM_cluster_probs.GMM_cluster_probs(~idx_20,:,depth_idx);...
                        GMM_cluster_probs.GMM_cluster_probs(idx_20,:,depth_idx)];
                    m_pcolor(double(Longitude),double(Latitude),double(z)');
                    title(gca,extractAfter(datestr(datenum(2004,m,1)),'-'));
                    clrs = flipud(jet(num_clusters)); % jet
                    colormap(customcolormap([0 1],[clrs(cl,:); 1 1 1]));
                    m_coast('patch',rgb('grey'));
                    m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
                    clim([0 100]);
                    c=colorbar;
                    c.Limits = [0 100];
                    c.Label.String = ['Cluster ' num2str(cl) ' Probability'];
                    c.TickLength = 0;
                    % save frame
                    if ~isfolder([dname '/' num2str(pressures(d)) 'dbars/probs'])
                        mkdir([dname '/' num2str(pressures(d)) 'dbars/probs']);
                    end
                    exportgraphics(h,[dname '/' num2str(pressures(d)) ...
                        'dbars/probs/m' num2str(m) '_w' num2str(w) '.png']);
                    % capture frame
                    frame = getframe(h);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    % write to file
                    if cnt == 1
                        imwrite(imind,cm,[dname '/' fname],'gif','Loopcount',inf,'DelayTime',0.1);
                    else
                        imwrite(imind,cm,[dname '/' fname],'gif','WriteMode','append','DelayTime',0.1);
                    end
                    % increase counter
                    cnt = cnt + 1;
                end
            end
        end
        close
    end
end

% end parallel session
delete(gcp('nocreate'));

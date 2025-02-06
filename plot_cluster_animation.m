%% Plot clusters over time

function plot_cluster_animation(param_props,fpath,base_grid,num_clusters,start_year,snap_date,numWorkers_train)

% process date
date_str = num2str(snap_date);

% set pressures
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    pressures = [2.5 10 50 100 200 300 500 1000 1500 1975];
else
    % define paths
    path2 = ['_Omon_' base_grid '_'];
    path3 = '_r1i1p1f1_gr';
    % define filename
    filename = [fpath 'combined/regridded/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    % load pressures
    pressures = ncread(filename,'depth');
end

% set up parallel pool
tic; parpool(numWorkers_train); fprintf('Pool initiation: '); toc;

parfor d = 1:length(pressures)
    % create folder
    dname = [param_props.p1 '/Figures/Clusters/' base_grid '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    % establish file name
    fname = ['cluster_animation_' num2str(pressures(d)) 'dbar.gif'];
    % determine number of monthly timesteps
    folder_name = [param_props.p1 '/Data/GMM_' base_grid '_' num2str(num_clusters)];
    cluster_inf = ncinfo([folder_name '/clusters.nc']);
    for dims = 1:length(cluster_inf.Dimensions)
        if strcmp(cluster_inf.Dimensions(dims).Name,'time')
            t_idx = dims;
        end
    end
    timesteps = cluster_inf.Dimensions(t_idx).Length;
    % load dimensions
    if strcmp(base_grid,'RG')
        % load
        Longitude = ncread([fpath '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Longitude');
        Latitude = ncread([fpath '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Latitude');
        Pressure = ncread([fpath '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Pressure');
    elseif strcmp(base_grid,'RFROM')
        % load
        Longitude = ncread([fpath '/Data/RFROM/RFROM_TEMP_STABLE_CLIM.nc'],'longitude');
        Latitude = ncread([fpath '/Data/RFROM/RFROM_TEMP_STABLE_CLIM.nc'],'latitude');
        Pressure = ncread([fpath '/Data/RFROM/RFROM_TEMP_STABLE_CLIM.nc'],'mean_pressure');
    else
        % define paths
        path2 = ['_Omon_' base_grid '_'];
        path3 = '_r1i1p1f1_gr';
        % define filename
        filename = [fpath 'combined/regridded/abs_sal' path2 ...
            'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
        % load
        Longitude = ncread(filename,'lon');
        Latitude = ncread(filename,'lat');
        Pressure = ncread(filename,'depth');
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
    h = figure('color','w','visible','on');
    axis tight manual
    % plot clusters each month/week
    for m = 1:timesteps
        if strcmp(base_grid,'RFROM')
            % determine number of weeks in file
            weeks = length(dir(['Data/GMM_' base_grid '_' num2str(num_clusters) '/m' num2str(m) '_*.mat']));
            for w = 1:weeks
                % clear frame
                clf
                % load monthly clusters
                GMM_clusters = ncread([folder_name '/clusters.nc'],...
                'clusters',[1 1 1 cnt],[Inf Inf Inf 1]);
                % make plot
                m_proj('robinson','lon',[20 380]);
                z = [GMM_clusters.GMM_clusters(~idx_20,:,depth_idx);...
                    GMM_clusters.GMM_clusters(idx_20,:,depth_idx)];
                m_pcolor(double(Longitude),double(Latitude),double(z)');
                title(gca,extractAfter(datestr(datenum(2004,m,1)),'-'));
                colormap([1,1,1;flipud(jet(num_clusters))]); % white then jet
                m_coast('patch',rgb('grey'));
                m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
                clim([-0.5 num_clusters+0.5]);
                c=colorbar;
                c.Limits = [0.5 num_clusters+0.5];
                c.Label.String = 'Cluster';
                c.TickLength = 0;
                % save frame
                if ~isfolder([dname '/' num2str(pressures(d)) 'dbars'])
                    mkdir([dname '/' num2str(pressures(d)) 'dbars']);
                end
                export_fig([dname '/' num2str(pressures(d)) 'dbars/m' num2str(m) '_w' num2str(w) '.png']);
                % capture frame
                frame = getframe(gcf);
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
        else % for RG and CMIP models
            % clear frame
            clf
            % load monthly clusters
            GMM_clusters = ncread([folder_name '/clusters.nc'],...
                'clusters',[1 1 1 cnt],[Inf Inf Inf 1]);
            % make plot
            m_proj('robinson','lon',[20 380]);
            z = [GMM_clusters(~idx_20,:,depth_idx);GMM_clusters(idx_20,:,depth_idx)];
            m_pcolor(double(Longitude),double(Latitude),double(z)');
            title(gca,extractAfter(datestr(datenum(2004,m,1)),'-'));
            colormap([1,1,1;flipud(jet(num_clusters))]); % white then jet
            m_coast('patch',rgb('grey'));
            m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
            clim([-0.5 num_clusters+0.5]);
            c=colorbar;
            c.Limits = [0.5 num_clusters+0.5];
            c.Label.String = 'Cluster';
            c.TickLength = 0;
            % save frame
            if ~isfolder([dname '/' num2str(pressures(d)) 'dbars'])
                mkdir([dname '/' num2str(pressures(d)) 'dbars']);
            end
            export_fig([dname '/' num2str(pressures(d)) 'dbars/m' num2str(m) '.png']);
            % capture frame
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % write to file
            if m == 1
                imwrite(imind,cm,[dname '/' fname],'gif','Loopcount',inf,'DelayTime',0.1);
            else
                imwrite(imind,cm,[dname '/' fname],'gif','WriteMode','append','DelayTime',0.1);
            end
            % increase counter
            cnt = cnt + 1;
        end
    end
    close
end

% end parallel session
delete(gcp('nocreate'));

%% Plot clusters over time

function plot_cluster_animation(param_props,fpaths,base_grid,num_clusters,...
    start_year,snap_date,numWorkers_train,flt,gld,ctd)

% define dataset extensions
if flt == 1; float_ext = 'f'; else float_ext = ''; end
if gld == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if ctd == 1; ctd_ext = 'w'; else ctd_ext = ''; end

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
    filename = [fpaths.param_path 'combined/regridded/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    % load pressures
    pressures = ncread(filename,'depth');
end

% set up parallel pool
tic; parpool(numWorkers_train); fprintf('Pool initiation: '); toc;

parfor d = 1:length(pressures)
    % create folder
    dname = [param_props.dir_name '/Figures/Clusters/' base_grid '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    % establish file name
    fname = ['cluster_animation_' float_ext glodap_ext ctd_ext '_' ...
        num2str(pressures(d)) 'dbar.gif'];
    % determine number of monthly timesteps
    folder_name = [fpaths.param_path 'GMM_' base_grid '_' num2str(num_clusters)];
    cluster_inf = ncinfo([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc']);
    for dims = 1:length(cluster_inf.Dimensions)
        if strcmp(cluster_inf.Dimensions(dims).Name,'time')
            t_idx = dims;
        end
    end
    timesteps = cluster_inf.Dimensions(t_idx).Length;
    % load dimensions
    Longitude = ncread([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc'],'lon');
    Latitude = ncread([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc'],'lat');
    Pressure = ncread([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc'],'pres');
    % process longitude
    idx_20 = Longitude<20;
    Longitude(idx_20) = Longitude(idx_20)+360;
    Longitude = [Longitude(~idx_20);Longitude(idx_20)];
    % depth index
    depth_idx = find(Pressure == pressures(d));
    % set counter
    cnt = 1;
    % establish figure
    h = figure('color','w','visible','off');
    axis tight manual
    % plot clusters each month/week
    for t = 1:timesteps
        % clear frame
        clf
        % load monthly clusters
        GMM_clusters = ncread([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc'],...
            'clusters',[1 1 1 t],[Inf Inf Inf 1]);
        time = ncread([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc'],'time',t,1);
        % make plot
        m_proj('robinson','lon',[20 380]);
        z = [GMM_clusters(~idx_20,:,depth_idx);GMM_clusters(idx_20,:,depth_idx)];
        m_pcolor(double(Longitude),double(Latitude),double(z)');
        if strcmp(base_grid,'RFROM')
            title(gca,datestr(time));
        else
            title(gca,extractAfter(datestr(datenum(2004,t,1)),'-'));
        end
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
        export_fig([dname '/' num2str(pressures(d)) 'dbars/t' num2str(t) '.png'],...
            '-transparent','-silent');
        % capture frame
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % write to file
        if t == 1
            imwrite(imind,cm,[dname '/' fname],'gif','Loopcount',inf,'DelayTime',0.1);
        else
            imwrite(imind,cm,[dname '/' fname],'gif','WriteMode','append','DelayTime',0.1);
        end
    end
    close
end

% end parallel session
delete(gcp('nocreate'));

%% Plot GOBAI over time

% set pressures
pressures = [2.5 10 50 100 200 300 500 1000 1500 1975];

% set up parallel pool
tic; parpool; fprintf('Pool initiation:'); toc;

parfor d = 1:length(pressures)
    % establish figure
    h=figure('visible','off','Position',[100 100 800 400]);
    axis tight manual
    % create folder
    dname = ['Figures/GOBAI/' base_grid '_' mod_type '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    % establish file name
    fname = ['gobai_animation_' num2str(pressures(d)) 'dbar.gif'];
    % determine number of monthly timesteps
    ds = dir(['Data/GOBAI/' dir_base '/*.nc']);
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
    % depth index
    depth_idx = find(Pressure == pressures(d));
    % plot gobai each month/week
    for m = 1:timesteps
        if strcmp(base_grid,'RG')
            % load monthly gobai
            gobai = ncread(['Data/GOBAI/' dir_base '/m' num2str(m) '_w1.nc'],'o2');
            % make plot
            worldmap([-90 90],[20 380]);
            title(extractAfter(datestr(datenum(2004,m,1)),'-'),'fontsize',16);
            pcolorm(double(Latitude),double(Longitude),...
                double(gobai(:,:,depth_idx))');
            colormap(cmocean('ice')); % white then jet
            plot_land('map');
            clim([0 350]);
            c=colorbar;
            c.Limits = [0 350];
            c.Label.String = '[O_{2}] (\mumol kg^{-1})';
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
                % load weekly gobai
                gobai = ncread(['Data/GOBAI/' dir_base '/m' num2str(m) '_w' num2str(w) '.nc'],'o2');
                % make plot
                worldmap([-90 90],[20 380]);
                title(extractAfter(datestr(datenum(2004,m,1)),'-'));
                pcolorm(double(Latitude),double(Longitude),...
                    double(gobai(:,:,depth_idx))');
                colormap(cmocean('ice')); % white then jet
                plot_land('map');
                clim([0 350]);
                c=colorbar;
                c.Limits = [0 350];
                c.Label.String = '[O_{2}] (\mumol kg^{-1})';
                c.TickLength = 0;
                mlabel off; plabel off;
                % capture frame
                frame = getframe(h);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                % write to file
                if m == 1 && w == 1
                    imwrite(imind,cm,[dname '/' fname],'gif','Loopcount',inf,'DelayTime',0.1);
                else
                    imwrite(imind,cm,[dname '/' fname],'gif','WriteMode','append','DelayTime',0.1);
                end
            end
        end
    end
    close
end

% end parallel session
delete(gcp('nocreate'));

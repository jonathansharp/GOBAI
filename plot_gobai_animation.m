%% Plot GOBAI over time

function plot_gobai_animation(param,dir_base,base_grid,num_clusters,mod_type)

% process parameter name
param1 = param_name(param);

% set pressures
pressures = [2.5 10 50 100 200 300 500 1000 1500 1975];

% set up parallel pool
tic; parpool(length(pressures)); fprintf('Pool initiation:'); toc;

parfor d = 1:length(pressures)
    % create folder
    dname = [param1 '/Figures/GOBAI/' base_grid '_' mod_type '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    % establish file name
    fname = ['gobai_animation_' num2str(pressures(d)) 'dbar.gif'];
    % determine number of monthly timesteps
    ds = dir([param1 '/Data/GOBAI/' dir_base '/*.nc']);
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
    for m = 1:timesteps
        if strcmp(base_grid,'RG')
            % clear frame
            clf
            % load monthly gobai
            gobai = ncread([param1 '/Data/GOBAI/' dir_base '/m' num2str(m) '_w1.nc'],param);
            % establish figure
            figure(h);
            % make plot
            m_proj('robinson','lon',[20 380]);
            z = [gobai(~idx_20,:,depth_idx);gobai(idx_20,:,depth_idx)];
            m_pcolor(double(Longitude),double(Latitude),double(z)');
            title(gca,extractAfter(datestr(datenum(2004,m,1)),'-'));
            colormap(cmocean('ice'));
            m_coast('patch',rgb('grey'));
            m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
            clim([0 350]);
            c=colorbar;
            c.Limits = [0 350];
            c.Label.String = '[O_{2}] (\mumol kg^{-1})';
            c.TickLength = 0;
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
                % load weekly gobai
                gobai = ncread([param1 '/Data/GOBAI/' dir_base '/m' num2str(m) '_w' num2str(w) '.nc'],param);
                % establish figure
                figure(h);
                % make plot
                m_proj('robinson','lon',[20 380]);
                z = [gobai(~idx_20,:,depth_idx);gobai(idx_20,:,depth_idx)];
                m_pcolor(double(Longitude),double(Latitude),double(z)');
                title(gca,extractAfter(datestr(datenum(2004,m,1)),'-'));
                colormap(cmocean('ice'));
                m_coast('patch',rgb('grey'));
                m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
                clim([0 350]);
                c=colorbar;
                c.Limits = [0 350];
                c.Label.String = '[O_{2}] (\mumol kg^{-1})';
                c.TickLength = 0;
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
    % display information
    dname
    disp(['GOBAI-' param1 ' (' mod_type ') animation at ' num2str(pressures(d)) ' dbar plotted'])
end

% end parallel session
delete(gcp('nocreate'));

end

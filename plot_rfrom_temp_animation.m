
%% Plot RFROM over time

function plot_rfrom_temp_animation(fpaths,ver,base_grid,start_year,end_year)

%% set pressures
pressures = [2.5 10 50 100 200 300 500 1000 1500 1975];

%% set delay time
if strcmp(base_grid,'RFROM')
    delay_time = 0.1;
else
    delay_time = (52/12)*0.1;
end

%% set up parallel pool
tic; parpool(length(pressures)); fprintf('Pool initiation: '); toc;

%% plot frames
parfor d = 1:length(pressures)
    % create folder for figures
    dname = ['Figures/RFROM/' ver '/'];
    % load dimensions
    TS = load_RFROM_dim(fpaths.temp_path,start_year,end_year);
    % depth index
    depth_idx = find(TS.Pressure == pressures(d));
    % establish fiugre
    h = figure('color','w','visible','off','Position',[616 474 1200 800]);
    axis tight manual
    % plot clusters each month/week
    for m = 1:length(TS.months)
        % establish file name
        fname = ['rfrom_temp_animation_' num2str(pressures(d)) 'dbar.gif'];
        % load dimensions
        TS = load_RFROM_dim(fpaths.temp_path,start_year,end_year);
        % process longitude
        idx_20 = TS.Longitude<20;
        TS.Longitude(idx_20) = TS.Longitude(idx_20)+360;
        TS.Longitude = [TS.Longitude(~idx_20);TS.Longitude(idx_20)];
        % TS = replicate_dims(base_grid,TS,1);
        % determine number of weeks in file
        nc_atts = ncinfo([fpaths.temp_path 'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
            % clear frame
            clf
            % counter
            cnt = TS.cnt{m}(w);
            % get RFROM T
            rfrom = ncread([fpaths.temp_path 'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_' ...
                num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                    'ocean_temperature',[1 1 depth_idx w],[Inf Inf 1 1]);
            time = datenum(1950,1,1) + TS.Time(cnt);
            % make plot
            m_proj('robinson','lon',[20 380]);
            z = [rfrom(~idx_20,:);rfrom(idx_20,:)];
            m_pcolor(double(TS.Longitude),double(TS.Latitude),double(z)');
            title(gca,datestr(time,'mmm-YYYY'),'FontSize',20);
            m_coast('patch',rgb('grey'));
            m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
            clim([-5 30]); c=colorbar;
            c.Label.String = ['Temperature (' char(176) 'C)'];
            c.TickLength = 0;
            colormap(cmocean('thermal'));
            set(gca,'FontSize',20);
            % save frame
            export_fig(h,[dname '/' num2str(pressures(d)) ...
                'dbars/temp_t' num2str(cnt) '.png'],'-transparent','-silent');
            % create folder
            if ~isfolder([dname '/' num2str(pressures(d)) 'dbars'])
                mkdir([dname '/' num2str(pressures(d)) 'dbars']);
            end
            if cnt == 1
                exportgraphics(h,[dname '/' fname],'Append',false);
            else
                exportgraphics(h,[dname '/' fname],'Append',true);
            end
        end 
    end
    close
    % display information
    disp(['RFROM Temperature animation at ' num2str(pressures(d)) ' dbar plotted'])
end

% end parallel session
delete(gcp('nocreate'));

end

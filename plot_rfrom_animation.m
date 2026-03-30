
%% Plot RFROM over time

function plot_rfrom_animation(fpaths,param,ver,base_grid,start_year,...
    end_year,min_val,max_val,step_val,varargin)

%% process optional input arguments
% pre-allocate
proj = 'normal';
plot_year = '';
% process inputs
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'projection')
        proj = varargin{i+1};
    elseif strcmpi(varargin{i}, 'plot_year')
        plot_year = varargin{i+1};
    end
end

%% set pressures
pressures = [2.5 10 50 100 200 300 500 1000 1500 1975];

%% set delay time
if strcmp(base_grid,'RFROM')
    delay_time = 0.1;
else
    delay_time = (52/12)*0.1;
end

%% set up parallel pool
% tic; parpool(length(pressures)); fprintf('Pool initiation: '); toc;

%% plot frames
for d = 4%1:length(pressures)
    % create folder for figures
    dname = ['Figures/RFROM/' ver];
    % load dimensions
    TS = load_RFROM_dim(fpaths.([param '_path']),ver,start_year,end_year);
    if isnumeric(plot_year)
        time_lim = TS.years == plot_year;
        TS.years = TS.years(time_lim);
        TS.months = TS.months(time_lim);
        TS.cnt = TS.cnt(time_lim);
        year_tag = ['_' num2str(plot_year)];
    else
        year_tag = '';
    end
    % depth index
    depth_idx = find(TS.Pressure == pressures(d));
    % establish fiugre
    h = figure('color','w','visible','off','Position',[616 474 1200 800]);
    axis tight manual
    % plot clusters each month/week
    for m = 1:length(TS.months)
        % establish file name
        fname = ['rfrom_' param '_animation_' num2str(pressures(d)) 'dbar' ...
            '_' proj year_tag '.gif'];
        % load dimensions
        TS = load_RFROM_dim(fpaths.([param '_path']),ver,start_year,end_year);
        % process longitude
        if strcmp(proj,'polar'); else
            idx_20 = TS.Longitude<20;
            TS.Longitude(idx_20) = TS.Longitude(idx_20)+360;
            TS.Longitude = [TS.Longitude(~idx_20);TS.Longitude(idx_20)];
        end
        % TS = replicate_dims(base_grid,TS,1);
        % determine number of weeks in file
        if strcmp(param,'temp')
            nc_atts = ncinfo([fpaths.([param '_path']) 'RFROM_TEMP_' ver '/RFROMV' ...
                ver(2) ver(4) '_TEMP_STABLE_' ...
                num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        elseif strcmp(param,'sal')
            nc_atts = ncinfo([fpaths.([param '_sal']) 'RFROM_SAL_' ver '/RFROMV' ...
                ver(2) ver(4) '_SAL_STABLE_' ...
                num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        end
        for w = 1:nc_atts.Dimensions(3).Length
            % clear frame
            clf
            % counter
            cnt = TS.cnt{m}(w);
            % get RFROM T
            if strcmp(param,'temp')
                rfrom = ncread([fpaths.([param '_path']) 'RFROM_TEMP_' ver '/RFROMV' ...
                    ver(2) ver(4) '_TEMP_STABLE_' ...
                    num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                        'ocean_temperature',[1 1 depth_idx w],[Inf Inf 1 1]);
            elseif strcmp(param,'sal')
                rfrom = ncread([fpaths.([param '_path']) 'RFROM_SAL_' ver '/RFROMV' ...
                    ver(2) ver(4) '_SAL_STABLE_' ...
                    num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                        'ocean_salinity',[1 1 depth_idx w],[Inf Inf 1 1]);
            end
            time = datenum(1950,1,1) + TS.Time(cnt);
            % make plot
            if strcmp(proj,'polar')
                m_proj('stereographic','lon',0,'lat',-90,'radius',60);
                m_pcolor(double([TS.Longitude;TS.Longitude(end)+1]),...
                    double(TS.Latitude),double([rfrom;rfrom(end,:)])');
                m_grid('linestyle','-','ytick',-90:20:90,...
                    'xtick',-180:30:180,'xaxislocation','top');
                set(gca,'XAxisLocation','bottom')
                m_coast('patch',[0.9 0.9 0.9]);
                m_text(0,-87,datestr(time,'dd-mmm-yyyy'),...
                    'HorizontalAlignment','center','FontSize',12);
            else
                m_proj('robinson','lon',[20 380]);
                z = [rfrom(~idx_20,:);rfrom(idx_20,:)];
                m_pcolor(double(TS.Longitude),double(TS.Latitude),double(z)');
                title(gca,datestr(time,'mmm-YYYY'),'FontSize',20);
                m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
                m_coast('patch',[0.9 0.9 0.9]);
            end
            set(gca,'FontSize',12);
            clim([min_val max_val]);
            c=colorbar('Location','southoutside','Ticks',min_val:step_val:max_val,...
                'TickLength',0.035,'TickLabelsMode','auto');
            if strcmp(param,'temp')
                c.Label.String = ['Temperature (' char(176) 'C)'];
            elseif strcmp(param,'sal')
                c.Label.String = 'Salinity';
            end
            colormap(cmocean('thermal',21));
            % save frame
            dirname = [dname '/' num2str(pressures(d)) 'dbars'];
            if ~exist(dirname,'dir'); mkdir(dirname); end
            export_fig(h,[dirname 'temp_t' num2str(cnt) '_' ...
                proj year_tag '.png'],'-transparent','-silent');
            frame=getframe(gcf);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % create folder
            if cnt == 1
                % imwrite(imind,cm,[dname '/' fname],'gif','Loopcount',inf,'DelayTime',.25);
                exportgraphics(h,[dname '/' fname],'Append',false);
            else
                % imwrite(imind,cm,[dname '/' fname],'gif','DelayTime',.25,'WriteMode','append');
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


%% Plot GOBAI over time

function plot_gobai_animation(param_props,base_grid,num_clusters,...
    alg_type,file_date,float_file_ext,numWorkers_predict,varargin)

%% set pressures
pressures = [2.5 10 50 100 200 300 500 1000 1500 1975];

%% process necessary input arguments for model parameters
% pre-allocate
train_ratio = NaN;
val_ratio = NaN;
test_ratio = NaN;
numtrees = NaN;
minLeafSize = NaN;
numstumps = NaN;
numbins = NaN;
% process based on algorithm type
if strcmp(alg_type,'FFNN')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'train_ratio')
            train_ratio = varargin{i+1};
        elseif strcmpi(varargin{i}, 'val_ratio')
            val_ratio = varargin{i+1};
        elseif strcmpi(varargin{i}, 'test_ratio')
            test_ratio = varargin{i+1};
        end
    end
elseif strcmp(alg_type,'RFR')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'numtrees')
            numtrees = varargin{i+1};
        elseif strcmpi(varargin{i}, 'minLeafSize')
            minLeafSize = varargin{i+1};
        end
    end
elseif strcmp(alg_type,'GBM')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'numstumps')
            numstumps = varargin{i+1};
        elseif strcmpi(varargin{i}, 'numbins')
            numbins = varargin{i+1};
        end
    end
else
    disp('"alg_type" must be "FFNN", "RFR", or "GBM"')
end

%% set up parallel pool
tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;

%% plot frames
parfor d = 1:length(pressures)
    % create folder for figures
    dname = [param_props.p1 '/Figures/GOBAI/' base_grid '_' alg_type '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    % define directory for file
    dir_base = create_dir_base(alg_type,{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
    % establish file name
    fname = ['gobai_animation_' num2str(pressures(d)) 'dbar.gif'];
    % determine number of monthly timesteps
    gobai_fname = [param_props.p1 '/Data/GOBAI/' dir_base '/gobai-' param_props.p2 '.nc'];
    gobai_inf = ncinfo(gobai_fname);
    for dims = 1:length(gobai_inf.Dimensions)
        if strcmp(gobai_inf.Dimensions(dims).Name,'time')
            t_idx = dims;
        end
    end
    timesteps = gobai_inf.Dimensions(t_idx).Length;
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
    h = figure('color','w','visible','on');
    axis tight manual
    % plot clusters each month/week
    for m = 1:timesteps
        if strcmp(base_grid,'RG')
            % clear frame
            clf
            % load monthly gobai
            gobai = ncread([param_props.p1 '/Data/GOBAI/' dir_base ...
                '/gobai-' param_props.p2 '.nc'],param_props.p2,[1 1 depth_idx m],[Inf Inf 1 1]);
            % establish figure
            figure(h);
            % make plot
            m_proj('robinson','lon',[20 380]);
            z = [gobai(~idx_20,:);gobai(idx_20,:)];
            m_pcolor(double(Longitude),double(Latitude),double(z)');
            title(gca,extractAfter(datestr(datenum(2004,m,1)),'-'));
            colormap(param_props.cmap);
            m_coast('patch',rgb('grey'));
            m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
            clim([param_props.edges(1) param_props.edges(end)]);
            c=colorbar;
            c.Limits = [param_props.edges(1) param_props.edges(end)];
            c.Label.String = [param_props.p3 ' ' param_props.units];
            c.TickLength = 0;
            % save frame
            if ~isfolder([dname '/' num2str(pressures(d)) 'dbars'])
                mkdir([dname '/' num2str(pressures(d)) 'dbars']);
            end
            export_fig(h,[dname '/' num2str(pressures(d)) ...
                'dbars/m' num2str(m) '_w1.png'],'-transparent');
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
                gobai = ncread([param_props.p1 '/Data/GOBAI/' dir_base '/m' num2str(m) '_w' num2str(w) '.nc'],param_props.p2);
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
                % save frame
                if ~isfolder([dname '/' num2str(pressures(d)) 'dbars'])
                    mkdir([dname '/' num2str(pressures(d)) 'dbars']);
                end
                export_fig(h,[dname '/' num2str(pressures(d)) ...
                    'dbars/m' num2str(m) '_w' num2str(w) '.png'],'-transparent');
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
    disp(['GOBAI-' param_props.p1 ' (' alg_type ') animation at ' num2str(pressures(d)) ' dbar plotted'])
end

% end parallel session
delete(gcp('nocreate'));

end

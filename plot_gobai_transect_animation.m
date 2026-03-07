%% Plot GOBAI transect over time

function plot_gobai_transect_animation(param_props,param_path,base_grid,num_clusters,...
    alg_type,file_date,float_file_ext,numWorkers_predict,varargin)

%% set longitudes
longitudes = 220;

%% set delay time
if strcmp(base_grid,'RFROM')
    delay_time = 0.1;
else
    delay_time = (52/12)*0.1;
end
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
elseif strcmp(alg_type,'AVG')
    % do nothing
else
    disp('"alg_type" must be "FFNN", "RFR", "GBM", or "AVG"')
end

%% set up parallel pool
% tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;

%% plot frames
for l = 1:length(longitudes)
    % create folder for figures
    dname = [param_props.dir_name '/Figures/GOBAI/transects/' base_grid '_' alg_type '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    % define directory for file
    if strcmp(alg_type,'FFNN')
        dir_base = [param_path 'GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) ...
            '_' file_date float_file_ext '/train' num2str(100*train_ratio) ...
            '_val' num2str(100*test_ratio) '_test' num2str(100*val_ratio)];
    elseif strcmp(alg_type,'RFR')
        dir_base = [param_path 'GOBAI/' base_grid '/RFR/c' num2str(num_clusters) ...
            '_' file_date float_file_ext '/tr' ...
            num2str(numtrees) '_lf' num2str(minLeafSize)];
    elseif strcmp(alg_type,'GBM')
        dir_base = [param_path 'GOBAI/' base_grid '/GBM/c' num2str(num_clusters) ...
            '_' file_date float_file_ext '/tr' num2str(numstumps) ...
            '_bin' num2str(numbins)];
    elseif strcmp(alg_type,'AVG')
        dir_base = [param_path 'GOBAI/' base_grid '/AVG/c' num2str(num_clusters) ...
            '_' file_date float_file_ext];
    end
    % establish file name
    fname = ['gobai_transect_animation_' num2str(longitudes(l)) 'deg.gif'];
    % determine number of timesteps
    %gobai_fname = [dir_base '/gobai-' param_props.file_name '.nc'];
    gobai_fname = [param_path '/GOBAI-' param_props.dir_name '-v1.0-HR.nc'];
    gobai_inf = ncinfo(gobai_fname);
    for dims = 1:length(gobai_inf.Dimensions)
        if strcmp(gobai_inf.Dimensions(dims).Name,'time')
            t_idx = dims;
        end
    end
    timesteps = gobai_inf.Dimensions(t_idx).Length;
    % load dimensions
    Longitude = ncread(gobai_fname,'lon');
    Latitude = ncread(gobai_fname,'lat',281,200);
    Pressure = ncread(gobai_fname,'pres',1,34);
    % longitude index
    [~,lon_idx] = min(abs(Longitude-longitudes(l)));
    % establish fiugre
    h = figure('color','w','visible','off','Position',[616 474 500 600]);
    axis tight manual
    % plot clusters each month/week
    for t = 1:timesteps
        % clear frame
        clf
        % load monthly gobai
        gobai = squeeze(ncread(gobai_fname,param_props.file_name,...
            [lon_idx 281 1 t],[1 200 34 1]));
        time = ncread(gobai_fname,'time',t,1)+datenum(1950,0,0);
        % establish figure
        figure(h);
        % make plot
        %pcolor(double(Latitude),double(Pressure),double(gobai)'); shading flat;
        contourf(double(Latitude),double(Pressure),double(gobai)');
        if strcmp(base_grid,'RFROM')
            title(gca,[param_props.label ' ' param_props.units ', ' datestr(time,'mmm-YYYY')],'FontSize',14);
        else
            title(gca,extractAfter(datestr(datenum(2004,t,1)),'-'),'FontSize',20);
        end
        set(gca,'YDir','reverse','FontSize',16);
        colormap(param_props.cmap);
        clim([param_props.edges(1) param_props.edges(end)]);
        c=colorbar;
        c.Limits = [param_props.edges(1) param_props.edges(end)];
        % c.Label.String = [param_props.label ' ' param_props.units];
        c.TickLength = 0;
        % save frame
        if ~isfolder([dname '/' num2str(longitudes(l)) 'deg'])
            mkdir([dname '/' num2str(longitudes(l)) 'deg']);
        end
        export_fig(h,[dname '/' num2str(longitudes(l)) ...
            'deg/t' num2str(t) '.png'],'-transparent','-silent');
        % capture frame
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % write to file
        if t == 1
            imwrite(imind,cm,[dname '/' fname],'gif','Loopcount',inf,'DelayTime',delay_time);
        else
            imwrite(imind,cm,[dname '/' fname],'gif','WriteMode','append','DelayTime',delay_time);
        end
    end
    close
    % display information
    disp(['GOBAI-' param_props.dir_name ' (' alg_type ') animation at ' num2str(longitudes(l)) ' deg plotted'])
end

% end parallel session
delete(gcp('nocreate'));

end

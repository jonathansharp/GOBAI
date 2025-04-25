
%% Plot GOBAI over time

function plot_gobai_animation(param_props,param_path,base_grid,num_clusters,...
    alg_type,file_date,float_file_ext,numWorkers_predict,varargin)

%% set pressures
pressures = [2.5 10 50 100 200 300 500 1000 1500 1975];

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
for d = 1:length(pressures)
    % create folder for figures
    dname = [param_props.dir_name '/Figures/GOBAI/' base_grid '_' alg_type '_c' num2str(num_clusters)];
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
    fname = ['gobai_animation_' num2str(pressures(d)) 'dbar.gif'];
    % determine number of timesteps
    gobai_fname = [dir_base '/gobai-' param_props.file_name '.nc'];
    gobai_inf = ncinfo(gobai_fname);
    for dims = 1:length(gobai_inf.Dimensions)
        if strcmp(gobai_inf.Dimensions(dims).Name,'time')
            t_idx = dims;
        end
    end
    timesteps = gobai_inf.Dimensions(t_idx).Length;
    % load dimensions
    Longitude = ncread(gobai_fname,'lon');
    Latitude = ncread(gobai_fname,'lat');
    Pressure = ncread(gobai_fname,'pres');
    % process longitude
    idx_20 = Longitude<20;
    Longitude(idx_20) = Longitude(idx_20)+360;
    Longitude = [Longitude(~idx_20);Longitude(idx_20)];
    % depth index
    depth_idx = find(Pressure == pressures(d));
    % establish fiugre
    h = figure('color','w','visible','on','Position',[616 474 1200 800]);
    axis tight manual
    % plot clusters each month/week
    for t = 1:timesteps
        % clear frame
        clf
        % load monthly gobai
        gobai = ncread(gobai_fname,param_props.file_name,...
            [1 1 depth_idx t],[Inf Inf 1 1]);
        time = ncread(gobai_fname,'time',t,1);
        % establish figure
        figure(h);
        set(gca,'FontSize',20);
        % make plot
        m_proj('robinson','lon',[20 380]);
        z = [gobai(~idx_20,:);gobai(idx_20,:)];
        m_pcolor(double(Longitude),double(Latitude),double(z)');
        if strcmp(base_grid,'RFROM')
            title(gca,datestr(datenum(1950,0,0)+time),'FontSize',20);
        else
            title(gca,extractAfter(datestr(datenum(2004,t,1)),'-'),'FontSize',20);
        end
        colormap(param_props.cmap);
        m_coast('patch',rgb('grey'));
        m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
        clim([param_props.edges(1) param_props.edges(end)]);
        c=colorbar;
        c.Limits = [param_props.edges(1) param_props.edges(end)];
        c.Label.String = [param_props.label ' ' param_props.units];
        c.TickLength = 0;
        % save frame
        if ~isfolder([dname '/' num2str(pressures(d)) 'dbars'])
            mkdir([dname '/' num2str(pressures(d)) 'dbars']);
        end
        export_fig(h,[dname '/' num2str(pressures(d)) ...
            'dbars/t' num2str(t) '.png'],'-transparent','-silent');
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
    disp(['GOBAI-' param_props.dir_name ' (' alg_type ') animation at ' num2str(pressures(d)) ' dbar plotted'])
end

% end parallel session
delete(gcp('nocreate'));

end

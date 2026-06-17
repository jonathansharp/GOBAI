
%% Plot GOBAI over time

param = 'o2';
param_props = param_config(param);
vrs = 'v1.1-HR';
pressure = 100;


folder = ['/raid/Data/GOBAI-' param_props.dir_name '/' vrs '/'];
file = ['GOBAI-' param_props.dir_name '-' vrs '.nc'];
fig_folder = [param_props.dir_name '/Figures/GOBAI/animations/'];
fig_file = ['GOBAI-' param_props.dir_name '-' vrs '.gif'];
if ~exist(fig_folder,"dir"); mkdir(fig_folder); end

% load dimensions
Longitude = ncread([folder file],'lon');
Latitude = ncread([folder file],'lat');
Pressure = ncread([folder file],'pres');
Time = ncread([folder file],'time');
% process longitude
idx_20 = Longitude<20;
Longitude(idx_20) = Longitude(idx_20)+360;
Longitude = [Longitude(~idx_20);Longitude(idx_20)];
% depth index
depth_idx = find(Pressure == pressure);
% establish fiugre
h = figure('color','w','visible','off','Position',[616 474 1200 800]);
axis tight manual
% establish video file
% v = VideoWriter([dname '/' v_fname]);
% plot clusters each month/week
for t = 1:length(Time)
    % clear frame
    clf
    % load monthly gobai
    gobai = ncread([folder file],param_props.file_name,...
        [1 1 depth_idx t],[Inf Inf 1 1]);
    time = datenum(1950,1,1) + ncread([folder file],'time',t,1);
    % make plot
    m_proj('robinson','lon',[20 380]);
    z = [gobai(~idx_20,:);gobai(idx_20,:)];
    m_pcolor(double(Longitude),double(Latitude),double(z)');
    title(gca,datestr(time,'mmm-YYYY'),'FontSize',20);
    m_coast('patch',rgb('grey'));
    m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
    clim([param_props.edges(1) param_props.edges(end)]);
    c=colorbar;
    c.Label.String = [param_props.label ' ' param_props.units];
    c.TickLength = 0;
    % define colormap
    colormap(param_props.cmap);
    set(gca,'FontSize',20);
    % write to file
    if t == 1
        exportgraphics(h,[fig_folder fig_file],Append=false);
    else
        exportgraphics(h,[fig_folder fig_file],Append=true);
    end
end

close


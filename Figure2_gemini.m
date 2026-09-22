% file information
path_start = '/raid/Data/';
var = {'TEMP';'SAL';'O2';'NO3';'DIC'};
varname = {'ocean_temperature';'ocean_salinity';'o2';'no3';'dic'};
varlabel = {['Temperature (' char(176) 'C)'];'Salinity';'[O_{2}] (\mumol kg^{-1})';...
    '[NO_{3}] (\mumol kg^{-1})';'DIC (\mumol kg^{-1})'};
cmap = {'thermal';'haline';'ice';'speed';'matter'};
% GOBAI properties
file_date = '202606';
ver_rfrom = 'v2.2-2025';
ver_gobai = 'HR-v1.3';

% establish figure
figure('Position',[100 100 1500 600]);
height = 0.4; width = 0.28; % slightly reduced height so bottom colorbars fit nicely
pos = [0.03 0.52 width height;...
       0.36 0.52 width height;...
       0.03 0.05 width height;...
       0.36 0.05 width height;...
       0.69 0.05 width height];

% time and depth
time = datenum(2017,6,15);
pres = 100;
lims = [0 25; 34 38; 0 300; 0 30; 1800 2200];

for p = 1:5
    if strcmp(var{p},'TEMP') || strcmp(var{p},'SAL')
        fpath = [path_start 'RFROM/RFROM_' var{p} '_' ver_rfrom '/RFROMV22_' ...
            var{p} '_STABLE_' num2str(year(time)) '_' sprintf('%02d',month(time)) '.nc'];
    else
        fpath = [path_start 'GOBAI-' var{p} '/' ver_gobai '/GOBAI-' var{p} ...
            '-HR-v' file_date '.nc'];
    end
    
    % download dimensions
    if strcmp(var{p},'TEMP') || strcmp(var{p},'SAL')
        ds.lon = double(ncread(fpath,'longitude'));
        ds.lon = convert_lon(ds.lon,'format','0-360');
        ds.lat = double(ncread(fpath,'latitude'));
        ds.pres = double(ncread(fpath,'mean_pressure'));
        ds.time = double(ncread(fpath,'time'));
    else
        ds.lon = double(ncread(fpath,'lon'));
        ds.lon = convert_lon(ds.lon,'format','0-360');
        ds.lat = double(ncread(fpath,'lat'));
        ds.pres = double(ncread(fpath,'pres'));
        ds.time = double(ncread(fpath,'time'));
    end
    
    % index
    [~,idx_pres] = min(abs(pres-ds.pres));
    [~,idx_time] = min(abs(time-(datenum(1950,1,1)+ds.time)));
    
    % --- MAP CREATION & LOCKING GCA ---
    ax_handle = worldmap([-90 90],[20 380]);
    ax.(['ax' num2str(p)]) = ax_handle;
    
    % Force current axes to this specific map handle
    axes(ax_handle);
    set(ax_handle, 'FontSize', 8, 'Position', pos(p,:));
    
    % Keep existing map features when adding pcolorm/land
    hold(ax_handle, 'on');
    
    mlabel off; plabel off;
    
    plot_var = squeeze(double(ncread(fpath,varname{p},...
        [1 1 idx_pres idx_time],[Inf Inf 1 1])));
    
    % Draw data and land (these use gca, which is locked to ax_handle now)
    pcolorm(ds.lat-0.25, ds.lon-0.25, plot_var');
    plot_land('map');
    
    % Colorbar & formatting
    c = colorbar(ax_handle, 'southoutside');
    c.Label.String = varlabel{p};
    c.Label.FontSize = 10;
    
    colormap(ax_handle, cmocean(cmap{p}));
    clim(ax_handle, lims(p,:));
    
    if strcmp(var{p},'TEMP') || strcmp(var{p},'SAL')
        title(ax_handle, ['RFROM-' var{p} ' at ' num2str(ds.pres(idx_pres)) ' dbar']);
    else
        title(ax_handle, ['GOBAI-' var{p} ' at ' num2str(ds.pres(idx_pres)) ' dbar']);
    end
    
    hold(ax_handle, 'off');
end

exportgraphics(gcf, 'Paper_Figs/Fig2.png', 'Resolution', 300);
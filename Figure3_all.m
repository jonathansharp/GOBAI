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
figure('Position',[100 100 816 1056],'Visible','on');
height = 0.3; width = 0.46;
pos = [0.02 0.67 width height;...
       0.52 0.67 width height;...
       0.02 0.33 width height;...
       0.52 0.33 width height;...
       0.26 0.01 width height];

% time and depth
time = datenum(2017,6,15);
pres = 100;
lims = [0 25; 34 38; 0 300; 0 30; 1850 2250];

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
    
    % plot map
    ax.(['ax' num2str(p)]) = axes('Position', pos(p,:));
    worldmap([-90 90],[20 380]);
    set(ax.(['ax' num2str(p)]),'FontSize',10);
    mlabel off; plabel off;
    plot_var = squeeze(double(ncread(fpath,varname{p},...
        [1 1 idx_pres idx_time],[Inf Inf 1 1])));
    pcolorm(ds.lat-0.25,ds.lon-0.25,plot_var');
    plot_land('map');
    c=colorbar(ax.(['ax' num2str(p)]),'southoutside');
    c.Label.String = varlabel{p};
    %c.Label.FontSize = 10;
    colormap(ax.(['ax' num2str(p)]),cmocean(cmap{p}));
    clim(ax.(['ax' num2str(p)]),lims(p,:));
    if strcmp(var{p},'TEMP') || strcmp(var{p},'SAL')
        title(ax.(['ax' num2str(p)]),['RFROM-' var{p} ' at ' num2str(ds.pres(idx_pres)) ' dbar']);
    else
        title(ax.(['ax' num2str(p)]),['GOBAI-' var{p} ' at ' num2str(ds.pres(idx_pres)) ' dbar']);
    end
    % reset position because matlab is dumb, or maybe I am
    %set(ax.(['ax' num2str(p)]),'Position', pos(p,:));

end

export_fig(gcf,'Figures/Figure3_global_maps.png','-transparent');
close;

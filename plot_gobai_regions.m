% file information
ver1 = 'v2.3'; % version
ver2 = 'v1.0-HR';
var1 = 'O2'; % var11iable
path1 = ['/raid/Data/GOBAI-' var1 '/' ver1 '/']; % file path
path2 = ['/raid/Data/GOBAI-' var1 '/' ver2 '/']; % file path
var2 = 'NO3'; % var11iable
path3 = ['/raid/Data/GOBAI-' var2 '/' ver2 '/']; % file path

% 
regions = {'NEP','GLOBAL','CCS','HUM','CAN','BENG'};
lat = [39 52.5; -90 90; 25 45.5; ...
    -30 0.5; 5 27.5; -35 -14.5];
lon = [219 240.5; 0 360; 220 250.5; ...
    270 290.5; 330 350.5; 5 19.5];
o2_lims = [120 300; 120 300; 120 300; ...
    0 300; 0 300; 0 300];
nit_lims = [0 30; 0 30; 0 25; ...
    0 25; 0 25; 0 25];
pres = [100; 100; 100; 100; 100; 100];
time = [datenum(2017,6,15);datenum(2017,6,15);...
    datenum(2017,6,15);datenum(2017,6,15);...
    datenum(2017,6,15);datenum(2017,6,15)];

for r = 1:length(regions)

    % download dimensions
    gobai1.lon = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'lon'));
    gobai1.lon = convert_lon(gobai1.lon,'format','0-360');
    gobai1.lat = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'lat'));
    gobai1.pres = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'pres'));
    gobai1.time = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'time'));
    gobai2.lon = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'lon'));
    gobai2.lat = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'lat'));
    gobai2.pres = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'pres'));
    gobai2.time = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'time'));
    gobai3.lon = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'lon'));
    gobai3.lat = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'lat'));
    gobai3.pres = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'pres'));
    gobai3.time = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'time'));
    
    % index
    idx_lon_1 = find(gobai1.lon > lon(r,1) & gobai1.lon <= lon(r,2));
    idx_lat_1 = find(gobai1.lat > lat(r,1) & gobai1.lat <= lat(r,2));
    idx_lon_2 = find(gobai2.lon > lon(r,1) & gobai2.lon <= lon(r,2));
    idx_lat_2 = find(gobai2.lat > lat(r,1) & gobai2.lat <= lat(r,2));
    [~,idx_pres_1] = min(abs(pres(r)-gobai1.pres));
    [~,idx_pres_2] = min(abs(pres(r)-gobai2.pres));
    [~,idx_pres_3] = min(abs(pres(r)-gobai3.pres));
    [~,idx_time_1] = min(abs(time(r)-(datenum(1950,1,1)+gobai1.time)));
    [~,idx_time_2] = min(abs(time(r)-(datenum(1950,1,1)+gobai2.time)));
    [~,idx_time_3] = min(abs(time(r)-(datenum(1950,1,1)+gobai3.time)));
    
    % download oxygen
    gobai1.oxy = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'oxy',...
        [min(idx_lon_1) min(idx_lat_1) idx_pres_1 idx_time_1],...
        [length(idx_lon_1) length(idx_lat_1) 1 1]));
    gobai2.oxy = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'o2',...
        [min(idx_lon_2) min(idx_lat_2) idx_pres_2 idx_time_2],...
        [length(idx_lon_2) length(idx_lat_2) 1 1]));
    gobai2.nit = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'no3',...
        [min(idx_lon_2) min(idx_lat_2) idx_pres_3 idx_time_3],...
        [length(idx_lon_2) length(idx_lat_2) 1 1]));
    
    % plot v2.3 map (O2)
    h1 = figure; worldmap(lat(r,:),lon(r,:));
    pcolorm(gobai1.lat(idx_lat_1)-0.25,gobai1.lon(idx_lon_1)-0.25,gobai1.oxy');
    plot_land('map');
    c=colorbar;
    colormap(cmocean('ice'));
    clim([o2_lims(r,1) o2_lims(r,2)]);
    title(['GOBAI-O_{2}-' ver1 ' at ' num2str(gobai1.pres(idx_pres_1)) ' dbar']);
    c.Label.String = '[O_{2}] (\mumol kg^{-1})';
    set(gca,'FontSize',12);
    mlabel('south');
    export_fig(['O2/Figures/Surface_Plots/' regions{r} '_' var1 '_' ...
        num2str(gobai1.pres(idx_pres_1)) 'dbar_' ver1 '.png'],'-transparent');
    close
    % plot v1.0-HR map (O2)
    h2 = figure; worldmap(lat(r,:),lon(r,:));
    pcolorm(gobai2.lat(idx_lat_2)-0.0625,gobai2.lon(idx_lon_2)-0.0625,gobai2.oxy');
    plot_land('map');
    c=colorbar;
    colormap(cmocean('ice'));
    clim([o2_lims(r,1) o2_lims(r,2)]);
    title(['GOBAI-O_{2}-' ver2 ' at ' num2str(gobai2.pres(idx_pres_2)) ' dbar']);
    c.Label.String = '[O_{2}] (\mumol kg^{-1})';
    set(gca,'FontSize',12);
    mlabel('south');
    export_fig(['O2/Figures/Surface_Plots/' regions{r} '_' var1 '_' ...
        num2str(gobai2.pres(idx_pres_2)) 'dbar_' ver2 '.png'],'-transparent');
    close
    % plot v1.0-HR map (NO3)
    h3 = figure; worldmap(lat(r,:),lon(r,:));
    pcolorm(gobai2.lat(idx_lat_2)-0.0625,gobai2.lon(idx_lon_2)-0.0625,gobai2.nit');
    plot_land('map');
    c=colorbar;
    colormap(cmocean('speed'));
    clim([nit_lims(r,1) nit_lims(r,2)]);
    title(['GOBAI-NO_{3}-' ver2 ' at ' num2str(gobai3.pres(idx_pres_3)) ' dbar']);
    c.Label.String = '[NO_{3}] (\mumol kg^{-1})';
    set(gca,'FontSize',12);
    mlabel('south');
    export_fig(['NO3/Figures/Surface_Plots/' regions{r} '_' var2 '_' ...
        num2str(gobai3.pres(idx_pres_3)) 'dbar_' ver2 '.png'],'-transparent');
    close

end
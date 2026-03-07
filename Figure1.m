% file information
ver1 = 'v2.3'; % version
ver2 = 'v1.0-HR';
var1 = 'O2'; % var11iable
path1 = ['/raid/Data/GOBAI-' var1 '/' ver1 '/']; % file path
path2 = ['/raid/Data/GOBAI-' var1 '/' ver2 '/']; % file path
var2 = 'NO3'; % var11iable
path3 = ['/raid/Data/GOBAI-' var2 '/' ver2 '/']; % file path

% 
regions = {'CCS','HUM','CAN','BENG'};
lat = [5 25.5; -25 -4.5; 5 25.5; -25 -4.5];
lon = [225 255.5; 265 295.5; 325 355.5; -10 19.5];
o2_lims = [0 300; 0 300; 0 300; 0 300];
nit_lims = [0 30; 0 30; 0 30; 0 30];
pres = [100; 100; 100; 100];
time = [datenum(2017,6,15);datenum(2017,6,15);...
    datenum(2017,6,15);datenum(2017,6,15)];

% establish figure
figure('Position',[100 100 1200 800]);
tiledlayout(4,3,'TileSpacing','compact','Padding','compact');

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
    lon = convert_lon(lon,'format','0-360');
    idx_lat_1 = find(gobai1.lat > lat(r,1) & gobai1.lat <= lat(r,2));
    idx_lat_2 = find(gobai2.lat > lat(r,1) & gobai2.lat <= lat(r,2));
    if r == 4
        idx_lon_1_1 = find(gobai1.lon > lon(r,1) & gobai1.lon <= 360);
        idx_lon_1_2 = find(gobai1.lon > 0 & gobai1.lon <= lon(r,2));
        idx_lon_2_1 = find(gobai2.lon > lon(r,1) & gobai2.lon <= 360);
        idx_lon_2_2 = find(gobai2.lon > 0 & gobai2.lon <= lon(r,2));
    else
        idx_lon_1 = find(gobai1.lon > lon(r,1) & gobai1.lon <= lon(r,2));
        idx_lon_2 = find(gobai2.lon > lon(r,1) & gobai2.lon <= lon(r,2));
    end
    [~,idx_pres_1] = min(abs(pres(r)-gobai1.pres));
    [~,idx_pres_2] = min(abs(pres(r)-gobai2.pres));
    [~,idx_pres_3] = min(abs(pres(r)-gobai3.pres));
    [~,idx_time_1] = min(abs(time(r)-(datenum(1950,1,1)+gobai1.time)));
    [~,idx_time_2] = min(abs(time(r)-(datenum(1950,1,1)+gobai2.time)));
    [~,idx_time_3] = min(abs(time(r)-(datenum(1950,1,1)+gobai3.time))); 
    
    % 
    if r == 4
        oxy_temp1 = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'oxy',...
            [min(idx_lon_1_2) min(idx_lat_1) idx_pres_1 idx_time_1],...
            [length(idx_lon_1_2) length(idx_lat_1) 1 1]));
        oxy_temp2 = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'oxy',...
            [min(idx_lon_1_1) min(idx_lat_1) idx_pres_1 idx_time_1],...
            [length(idx_lon_1_1) length(idx_lat_1) 1 1]));
        gobai1.oxy = [oxy_temp2;oxy_temp1];
        oxy_temp1 = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'o2',...
            [min(idx_lon_2_1) min(idx_lat_2) idx_pres_1 idx_time_1],...
            [length(idx_lon_2_1) length(idx_lat_2) 1 1]));
        oxy_temp2 = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'o2',...
            [min(idx_lon_2_2) min(idx_lat_2) idx_pres_1 idx_time_1],...
            [length(idx_lon_2_2) length(idx_lat_2) 1 1]));
        gobai2.oxy = [oxy_temp1;oxy_temp2];
        nit_temp1 = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'no3',...
            [min(idx_lon_2_1) min(idx_lat_2) idx_pres_1 idx_time_1],...
            [length(idx_lon_2_1) length(idx_lat_2) 1 1]));
        nit_temp2 = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'no3',...
            [min(idx_lon_2_2) min(idx_lat_2) idx_pres_1 idx_time_1],...
            [length(idx_lon_2_2) length(idx_lat_2) 1 1]));
        gobai2.nit = [nit_temp1;nit_temp2];
    else
        % download variables
        gobai1.oxy = double(ncread([path1 'GOBAI-' var1 '-' ver1 '.nc'],'oxy',...
            [min(idx_lon_1) min(idx_lat_1) idx_pres_1 idx_time_1],...
            [length(idx_lon_1) length(idx_lat_1) 1 1]));
        gobai2.oxy = double(ncread([path2 'GOBAI-' var1 '-' ver2 '.nc'],'o2',...
            [min(idx_lon_2) min(idx_lat_2) idx_pres_2 idx_time_2],...
            [length(idx_lon_2) length(idx_lat_2) 1 1]));
        gobai2.nit = double(ncread([path3 'GOBAI-' var2 '-' ver2 '.nc'],'no3',...
            [min(idx_lon_2) min(idx_lat_2) idx_pres_3 idx_time_3],...
            [length(idx_lon_2) length(idx_lat_2) 1 1]));
    end

    % plot v2.3 map (O2)
    axis.(['ax1_' num2str(r)]) = nexttile;
    worldmap(lat(r,:),lon(r,:));
    setm(axis.(['ax1_' num2str(r)]),'MapProjection','mercator','FontSize',12);
    mlabel('south');
    set(findobj(axis.(['ax1_' num2str(r)]).Children,'Tag','MLabel'),'FontSize',6)
    set(findobj(axis.(['ax1_' num2str(r)]).Children,'Tag','PLabel'),'FontSize',6)
    if r == 4    
        pcolorm(gobai1.lat(idx_lat_1)-0.25,...
            [gobai1.lon(idx_lon_1_1);gobai1.lon(idx_lon_1_2)]-0.25,gobai1.oxy');
        textm(mean(gobai1.lat(idx_lat_1(end)-2)),gobai1.lon(idx_lon_1_1(1)),...
                regions{r},'FontWeight','bold','FontSize',16);
    else
        pcolorm(gobai1.lat(idx_lat_1)-0.25,gobai1.lon(idx_lon_1)-0.25,gobai1.oxy');
        if r == 2 
            textm(mean(gobai1.lat(idx_lat_1(end)-2)),gobai1.lon(idx_lon_1(1)),...
                regions{r},'FontWeight','bold','FontSize',16,'Color','w');
        else
            textm(mean(gobai1.lat(idx_lat_1(end)-2)),gobai1.lon(idx_lon_1(1)),...
                regions{r},'FontWeight','bold','FontSize',16);
        end
    end
    plot_land('map');
    c=colorbar(axis.(['ax1_' num2str(r)]));
    c.Label.String = '[O_{2}] (\mumol kg^{-1})';
    c.Label.FontSize = 12;
    colormap(axis.(['ax1_' num2str(r)]),cmocean('ice'));
    clim(axis.(['ax1_' num2str(r)]),[o2_lims(r,1) o2_lims(r,2)]);
    if r == 1
        title(axis.(['ax1_' num2str(r)]),['GOBAI-O_{2}-' ver1 ' at ' num2str(gobai1.pres(idx_pres_1)) ' dbar']);
    end

    % plot v1.0-HR map (O2)
    axis.(['ax2_' num2str(r)]) = nexttile;
    worldmap(lat(r,:),lon(r,:));
    setm(axis.(['ax2_' num2str(r)]),'MapProjection','mercator','FontSize',12);
    mlabel('south');
    set(findobj(axis.(['ax2_' num2str(r)]).Children,'Tag','MLabel'),'FontSize',6)
    set(findobj(axis.(['ax2_' num2str(r)]).Children,'Tag','PLabel'),'FontSize',6);
    if r == 4    
        pcolorm(gobai2.lat(idx_lat_2)-0.25,...
            [gobai2.lon(idx_lon_2_1);gobai2.lon(idx_lon_2_2)]-0.25,gobai2.oxy');
    else
        pcolorm(gobai2.lat(idx_lat_2)-0.25,gobai2.lon(idx_lon_2)-0.25,gobai2.oxy');
    end
    plot_land('map');
    c=colorbar(axis.(['ax2_' num2str(r)]));
    c.Label.String = '[O_{2}] (\mumol kg^{-1})';
    c.Label.FontSize = 12;
    colormap(axis.(['ax2_' num2str(r)]),cmocean('ice'));
    clim(axis.(['ax2_' num2str(r)]),[o2_lims(r,1) o2_lims(r,2)]);
    if r == 1
        title(axis.(['ax2_' num2str(r)]),['GOBAI-O_{2}-' ver2 ' at ' num2str(gobai2.pres(idx_pres_2)) ' dbar']);
    end

    % plot v1.0-HR map (NO3)
    axis.(['ax3_' num2str(r)]) = nexttile;
    worldmap(lat(r,:),lon(r,:));
    setm(axis.(['ax3_' num2str(r)]),'MapProjection','mercator','FontSize',12);
    mlabel('south');
    set(findobj(axis.(['ax3_' num2str(r)]).Children,'Tag','MLabel'),'FontSize',6)
    set(findobj(axis.(['ax3_' num2str(r)]).Children,'Tag','PLabel'),'FontSize',6)
    if r == 4    
        pcolorm(gobai2.lat(idx_lat_2)-0.25,...
            [gobai2.lon(idx_lon_2_1);gobai2.lon(idx_lon_2_2)]-0.25,gobai2.nit');
    else
        pcolorm(gobai2.lat(idx_lat_2)-0.25,gobai2.lon(idx_lon_2)-0.25,gobai2.nit');
    end
    plot_land('map');
    c=colorbar(axis.(['ax3_' num2str(r)]));
    c.Label.String = '[NO_{3}] (\mumol kg^{-1})';
    c.Label.FontSize = 12;
    colormap(axis.(['ax3_' num2str(r)]),cmocean('speed'));
    clim(axis.(['ax3_' num2str(r)]),[nit_lims(r,1) nit_lims(r,2)]);
    if r ==1
        title(axis.(['ax3_' num2str(r)]),['GOBAI-NO_{3}-' ver2 ' at ' num2str(gobai3.pres(idx_pres_3)) ' dbar']);
    end

end

% export_fig(gcf,'Paper_Figs/Fig1.png','-transparent');
exportgraphics(gcf,'Paper_Figs/Fig1.png');

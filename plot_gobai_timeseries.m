

% latitudes and longitudes to test
lats = [39 36 20 -40 45 -50 -33 13];
lons = [5 16 320 30 190 200 348 147];
press = [2 2 2 2 2 2 2 2];

% load dimensions
gobaiv23 = load_gobai('O2','v2.3','/raid/Data/GOBAI-O2/v2.3/',{});
gobaiv23.lon = convert_lon(gobaiv23.lon,'format','0-360');
gobaiv22 = load_gobai('O2','v2.2','/raid/Data/GOBAI-O2/v2.2/',{});
gobaiv22.lon = convert_lon(gobaiv22.lon,'format','0-360');
gobaivhr.lon = ncread('/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-v1.0-HR.nc','lon');
gobaivhr.lat = ncread('/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-v1.0-HR.nc','lat');
gobaivhr.pres = ncread('/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-v1.0-HR.nc','pres');
gobaivhr.time = ncread('/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-v1.0-HR.nc','time');
% map figure
map_fig = figure; hold on;
worldmap('world'); plot_land('map');
clrs = repmat(colororder,2,1);

for n = 1:length(lats)
    % plot location
    figure(map_fig);
    scatterm(lats(n),lons(n),100,clrs(n,:),'.');
    % plot timeseries
    f.(['timeseries_' num2str(n)]) = figure('position',[100 300 1600 400]);
    hold on
    % v2.2
    [~,idx_lat] = min(abs(gobaiv22.lat-lats(n)));
    [~,idx_lon] = min(abs(gobaiv22.lon-lons(n)));
    [~,idx_pres] = min(abs(gobaiv22.pres-press(n)));
    o2 = ncread('/raid/Data/GOBAI-O2/v2.2/GOBAI-O2-v2.2.nc','oxy',...
        [idx_lon(1),idx_lat(1),idx_pres(1),1],[1 1 1 Inf]);
    plot(double(datenum(1950,1,1)+gobaiv22.time),double(squeeze(o2))); 
    % v2.3
    [~,idx_lat] = min(abs(gobaiv23.lat-lats(n)));
    [~,idx_lon] = min(abs(gobaiv23.lon-lons(n)));
    [~,idx_pres] = min(abs(gobaiv23.pres-press(n)));
    o2 = ncread('/raid/Data/GOBAI-O2/v2.3/GOBAI-O2-v2.3.nc','oxy',...
        [idx_lon(1),idx_lat(1),idx_pres(1),1],[1 1 1 Inf]);
    plot(double(datenum(1950,1,1)+gobaiv23.time),double(squeeze(o2)));
    % v1.0-HR
%     [~,idx_lat] = min(abs(gobaivhr.lat-lats(n)));
%     [~,idx_lon] = min(abs(gobaivhr.lon-lons(n)));
%     [~,idx_pres] = min(abs(gobaivhr.pres-press(n)));
%     o2 = ncread('/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-v1.0-HR.nc','o2',...
%         [idx_lon(1),idx_lat(1),idx_pres(1),1],[1 1 1 Inf]);
%     plot(double(datenum(1950,1,1)+gobaivhr.time),double(squeeze(o2)));
    % figure properties
    title(['Lat = ' num2str(lats(n)) ', Lon = ' num2str(lons(n)) ...
        ', Pres = ' num2str(press(n))]);
    box on;
    ylabel('[O_{2}] (\mumol kg^{-1})');
    legend({'v2.2' 'v2.3'});
    % legend({'v2.2' 'v2.3' 'v1.0-HR'});
    datetick('x');
    % save figure
    hold off
    export_fig(['Figures/Lat' num2str(lats(n)) '_Lon' num2str(lons(n)) ...
        '_Pres' num2str(press(n)) '.png'],'-transparent');
end

figure(map_fig);
set(map_fig,'color','w');
export_fig('Figures/timeseries_map.png','');
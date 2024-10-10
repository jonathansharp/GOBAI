%% Plot O2 over months (gif)

% file information
ver = 'v2.2'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
lon = -150.5; % longitude
% file names
filename = ['Figures/O2_section_by_month_lon' num2str(lon) '_' ver '.gif'];
% set loop parameters
n=1; year = 2004:2023;
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% set longitude index
lon_idx = find(GOBAI.lon==lon);
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
h1=figure;

% loop over time dimensions
for y = 1:length(year)
    for m = 1:12

        % download oxygen
        GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy',[1 1 1 (y-1)*12+m],[Inf Inf Inf 1]);
        lat_idx = any(GOBAI.oxy(lon_idx,:,:),3);

        % O2 map
        figure(h1); set(h1,'Position',[100 100 800 400]); clf;
        % m_proj('orthographic','lon',260,'lat',0)
        contourf(double(GOBAI.lat(lat_idx)),double(GOBAI.pres),...
            double(squeeze(GOBAI.oxy(lon_idx,lat_idx,:)))',[0:25:350],'LineStyle','-','ShowText','on');
        set(gca,'YDir','reverse');      
        title([num2str(m) '/' num2str(year(y)) ' - Mean [O_{2}] at ' ...
            num2str(lon) ' E (\mumol kg^{-1})'],'fontsize',12);
        c=colorbar('location','southoutside');
        caxis([0 350]);
        colormap(cmocean('ice',28));
        c.FontSize = 12;
        c.TickLength = 0;
        if n == 1
            exportgraphics(h1,filename,'Append',false);
        else
            exportgraphics(h1,filename,'Append',true);
        end
        % increase counter
        n=n+1;

    end
end

close all
clear

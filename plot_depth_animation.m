% Plot O2 over depth (gif)

% file information
ver = 'v2.1'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% file names
filename1 = ['Figures/O2_by_depth_' ver '.gif'];
% set loop parameters
n=1; year = 2004:2023;
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% process longitude
idx_20 = GOBAI.lon<20;
GOBAI.lon(idx_20) = GOBAI.lon(idx_20)+360;
GOBAI.lon = [GOBAI.lon(~idx_20);GOBAI.lon(idx_20)];
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
GOBAI.oxy = mean(GOBAI.oxy,4,'omitnan');
% establish figure
h1=figure;

% loop over time dimensions
for z = 1:length(GOBAI.pres)

    % O2 map
    figure(h1); clf;
    m_proj('orthographic','lon',260,'lat',0)
    m_pcolor(double(GOBAI.lon)-0.5,double(GOBAI.lat)-0.5,GOBAI.oxy(:,:,z)');
    m_coast('patch',rgb('grey'));
    m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
    title('Mean [O_{2}] (\mumol kg^{-1})','fontsize',12)
    c=colorbar('location','southoutside');
    caxis([0 400]);
    colormap(cmocean('ice',28));
    c.FontSize = 12;
    c.TickLength = 0;
    m_text(283,-7,[num2str(GOBAI.pres(z)) ' dbar'],'fontsize',8,...
        'fontweight','bold','color','k');
    if n == 1
        exportgraphics(h1,filename1,'Append',false);
    else
        exportgraphics(h1,filename1,'Append',true);
    end

    % increase counter
    n=n+1;

end

close all
clear
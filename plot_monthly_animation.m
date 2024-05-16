%% Plot O2 over months (gif)

% file information
ver = 'v2.2'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
pres = 200; % pressure level
% file names
filename1 = ['Figures/O2_by_month_' num2str(pres) 'dbar_' ver '.gif'];
filename2 = ['Figures/O2_inv_by_month_' num2str(pres) 'dbar_' ver '.gif'];
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
% set pressure index
pres_idx = find(GOBAI.pres==pres);
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
h1=figure; h2=figure;

% loop over time dimensions
for y = 1:length(year)
    for m = 1:12

        % download oxygen
        GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy',[1 1 1 (y-1)*12+m],[Inf Inf Inf 1]);
        
        % O2 map
        figure(h1); clf;
        m_proj('orthographic','lon',260,'lat',0)
        m_pcolor(double(GOBAI.lon)-0.5,double(GOBAI.lat)-0.5,GOBAI.oxy(:,:,pres_idx)');
        m_coast('patch',rgb('grey'));
        m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
        title(['Mean [O_{2}] at ' num2str(pres) ' dbar (\mumol kg^{-1})'],'fontsize',12)
        c=colorbar('location','southoutside');
        caxis([0 350]);
        colormap(cmocean('ice',28));
        c.FontSize = 12;
        c.TickLength = 0;
        m_text(50,55,[num2str(m) '/' num2str(year(y))],'fontsize',12,...
            'fontweight','bold','color','w');
        if n == 1
            exportgraphics(h1,filename1,'Append',false);
        else
            exportgraphics(h1,filename1,'Append',true);
        end

%         % O2 inventory
%         figure(h2); clf;
%         mean_O2 = squeeze((sum(sum(temp_O2.*area_weights_temp,1,'omitnan'),2,'omitnan'))./...
%             (sum(sum(area_weights,1,'omitnan'),2,'omitnan')));
%         set(h2,'color','white');
%         set(gcf,'position',[200 200 675 400]);
%         plot(GOBAI.time(1:(y-1)*12+m),mean_O2(1:(y-1)*12+m),'-','linewidth',3,'color',GOBAIb('blue'));
%         ylim([min(mean_O2)-0.5 max(mean_O2)+0.5]);
%         ylabel(['Mean Global [O_{2}] at ' num2str(pres) ' dbars (\mumol kg^{-1})']);
%         xlim([min(GOBAI.time) max(GOBAI.time)]);
%         datetick('x','keeplimits');
%         set(gca,'fontsize',18);
%         if n == 1
%             exportgraphics(h2,filename2,'Append',false);
%         else
%             exportgraphics(h2,filename2,'Append',true);
%         end

        % increase counter
        n=n+1;

    end
end

close all
clear

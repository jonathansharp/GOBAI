%% Plot O2 trend at specified depths
plot_O2_trend([-90 90],[20 380],20,RG,-2,2,21,xdim,ydim);
plot_O2_trend([-90 90],[20 380],300,RG,-2,2,21,xdim,ydim);
plot_O2_trend([-90 90],[20 380],1000,RG,-2,2,21,xdim,ydim);

%% Plot integrated O2 trend
plot_O2_trend_int([-90 90],[20 380],[0 2000],RG,-1,1,21,xdim,ydim,volume_weights,area_weights);

%% Embedded function
function plot_O2_trend(latlim,lonlim,depth,RG,min_trend,max_trend,levels,xdim,ydim)
figure; worldmap(latlim,lonlim);
title(['Oxygen trend at ' num2str(depth) ' dbar (\mumol kg^{-1} yr^{-1})'],'fontsize',16);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
idx = find(abs(RG.pressure-depth)==min(abs(RG.pressure-depth)));
trend = nan(xdim,ydim);
for a = 1:xdim
    for b = 1:ydim
        coeffs = polyfit(RG.month,squeeze(RG.oxy_ENS(a,b,idx,:)),1);
        trend(a,b) = coeffs(1).*12;
    end
end
pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
    double(repmat(RG.longitude,1,size(RG.latitude,1))),trend);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
caxis([min_trend-0.1 max_trend+0.1]);
colormap(cmocean('balance',levels));
c.FontSize = 12;
c.TickLength = 0;
mlabel off; plabel off;
exportgraphics(gcf,['Figures/Surface Plots/oxy_trend_' num2str(depth) '_dbar_ENS.png']);
end

%% Embedded function
function plot_O2_trend_int(latlim,lonlim,depths,RG,min_trend,max_trend,levels,xdim,ydim,volume_weights,area_weights)
figure; worldmap(latlim,lonlim);
title(['Average Oxygen Trend (\mumol kg^{-1} yr^{-1})'],'fontsize',16);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
idx = find(RG.pressure>min(depths) & RG.pressure<max(depths));
oxy_tmp = squeeze(sum(RG.oxy_ENS(:,:,idx,:).*volume_weights(:,:,idx,:),...
    3,'omitnan')./sum(volume_weights(:,:,idx,:),3,'omitnan'));
trend = nan(xdim,ydim);
for a = 1:xdim
    for b = 1:ydim
        coeffs = polyfit(RG.month,squeeze(oxy_tmp(a,b,:)),1);
        trend(a,b) = coeffs(1).*12;
    end
end
avg_trend = sum(trend(:).*area_weights(:),'omitnan')./sum(area_weights(:),'omitnan');
pcolorm(double(RG.latitude),double([RG.longitude;RG.longitude(end)]),...
    [trend;trend(end,:)]');
textm(58,60,num2str(round(10*avg_trend,1)),'color','w','fontweight','b');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
caxis([min_trend-0.05 max_trend+0.05]);
colormap(cmocean('balance',levels));
c.FontSize = 12;
c.TickLength = 0;
mlabel off; plabel off;
exportgraphics(gcf,['Figures/Surface Plots/oxy_trend_' num2str(depths(1)) '_' num2str(depths(2)) '_dbar_ENS.png']);
end

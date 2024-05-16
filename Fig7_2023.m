%% Figure 7
% Interannual global means on different depth levels

% file information
ver = 'v2.2'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-annual-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-annual-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-annual-' ver '.nc'],'pres');
GOBAI.year = ncread([path 'GOBAI-' var '-annual-' ver '.nc'],'year');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-annual-' ver '.nc'],'oxy');
% set dimension
xdim = length(GOBAI.lon);
ydim = length(GOBAI.lat);
zdim = length(GOBAI.pres);
tdim = length(GOBAI.year);

%% initialize figure
figure; set(gcf,'units','inches','position',[0 14 24 12]);

%% set intervals
int1 = GOBAI.pres >= 0 & GOBAI.pres <= 100;
int2 = GOBAI.pres >= 100 & GOBAI.pres <= 600;
int3 = GOBAI.pres >= 600 & GOBAI.pres <= 2000;

%% set colors
clrs_o = cmocean('ice',10); % oxygen
clrs_t = cmocean('solar',10); % temp
clrs_s = flipud(cmocean('algae',10)); % saturation

%% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;

%% calculate weighted means of oxygen
mean_upper = squeeze(sum(sum(sum(GOBAI.oxy(:,:,int1,:).*...
    GOBAI.vol(:,:,int1),3,'omitnan'),2,'omitnan'),1,'omitnan')./...
    sum(sum(sum(GOBAI.vol(:,:,int1),3,'omitnan'),2,'omitnan'),1,'omitnan'));
mean_middle = squeeze(sum(sum(sum(GOBAI.oxy(:,:,int2,:).*...
    GOBAI.vol(:,:,int2),3,'omitnan'),2,'omitnan'),1,'omitnan')./...
    sum(sum(sum(GOBAI.vol(:,:,int2),3,'omitnan'),2,'omitnan'),1,'omitnan'));
mean_lower = squeeze(sum(sum(sum(GOBAI.oxy(:,:,int3,:).*...
    GOBAI.vol(:,:,int3),3,'omitnan'),2,'omitnan'),1,'omitnan')./...
    sum(sum(sum(GOBAI.vol(:,:,int3),3,'omitnan'),2,'omitnan'),1,'omitnan'));
mean_depth = squeeze(sum(sum(GOBAI.oxy.*...
    GOBAI.vol,2,'omitnan'),1,'omitnan')./...
    sum(sum(GOBAI.vol,2,'omitnan'),1,'omitnan'));

%% plot weighted means of oxygen
mean_upper_uncer = 0.420;
mean_middle_uncer = 0.177;
mean_lower_uncer = 0.104;
ax1 = axes; hold on; box on; grid off;
set(ax1,'units','normalized','position',[0.05 0.56 0.25 0.42]);
set(ax1,'fontsize',18,'xticklabels',[],'linewidth',2);
ylabel('[O_{2}] Anomaly (\mumol kg^{-1})');
ylim([-2.5 2.5]);
xlim([0 20]);
xticks([2 7 12 17]);
xticklabels({'2005' '2010' '2015' '2020'});
p1=plot(1:tdim,mean_upper-mean(mean_upper),'linewidth',5,'color',clrs_o(8,:));
scatter(1:tdim,mean_upper-mean(mean_upper),100,'o','MarkerFaceColor',clrs_o(8,:),'MarkerEdgeColor','none');
fill([1:tdim fliplr(1:tdim)],[mean_upper+mean_upper_uncer-mean(mean_upper);...
    flipud(mean_upper-mean_upper_uncer-mean(mean_upper))],clrs_o(8,:),...
    'FaceAlpha',0.2,'LineStyle','none');
p2=plot(1:tdim,mean_middle-mean(mean_middle),'linewidth',5,'color',clrs_o(6,:));
scatter(1:tdim,mean_middle-mean(mean_middle),100,'o','MarkerFaceColor',clrs_o(6,:),'MarkerEdgeColor','none');
fill([1:tdim fliplr(1:tdim)],[mean_middle+mean_middle_uncer-mean(mean_middle);...
    flipud(mean_middle-mean_middle_uncer-mean(mean_middle))],clrs_o(6,:),...
    'FaceAlpha',0.2,'LineStyle','none');
p3=plot(1:tdim,mean_lower-mean(mean_lower),'linewidth',5,'color',clrs_o(4,:));
scatter(1:tdim,mean_lower-mean(mean_lower),100,'o','MarkerFaceColor',clrs_o(4,:),'MarkerEdgeColor','none');
fill([1:tdim fliplr(1:tdim)],[mean_lower+mean_lower_uncer-mean(mean_lower);...
    flipud(mean_lower-mean_lower_uncer-mean(mean_lower))],clrs_o(4,:),...
    'FaceAlpha',0.2,'LineStyle','none');
text(1,2.25,'a','fontsize',40);

%% add legend
legend([p1 p2 p3],{'0-100 dbars' '100-600 dbars' '600-2000 dbars'},...
    'location','southwest','NumColumns',1);

%% plot weighted means of oxygen with depth
ax2 = axes; hold on; box on;
set(ax2,'units','normalized','position',[0.06 0.05 0.23 0.42]);
set(ax2,'fontsize',18,'xticklabels',[],'YDir','reverse');
set(ax2,'layer','top','Box','on','linewidth',2);
contourf(1:tdim,GOBAI.pres,mean_depth-mean(mean_depth,2),...
    -2.625:0.25:2.625,'LineStyle','none');
ylabel('Pressure (dbar)');
ylim([0 1975]);
xlim([1 19]);
xticks([2 7 12 17]);
xticklabels({'2005' '2010' '2015' '2020'});
caxis([-2.625 2.625]);
mycolormap = cmocean('balance',21,'pivot',0);
mycolormap(11,:) = 1;
ax2.Colormap = mycolormap;
c=colorbar;
c.Label.String = '[O_{2}] Anomaly (\mumol kg^{-1})';
c.TickLength = 0;
c.Label.FontSize = 18;
c.Position(3) = c.Position(3).*2;
c.Position(1) = c.Position(1)+0.03;
text(-1.2,0,'d','fontsize',40);

%% calculate weighted means of temperature
temp_mean_upper = squeeze(sum(sum(sum(RG_ann.temp(:,:,int1,:).*...
    GOBAI.vol(:,:,int1),3,'omitnan'),2,'omitnan'),1,'omitnan')./...
    sum(sum(sum(GOBAI.vol(:,:,int1),3,'omitnan'),2,'omitnan'),1,'omitnan'));
temp_mean_middle = squeeze(sum(sum(sum(RG_ann.temp(:,:,int2,:).*...
    GOBAI.vol(:,:,int2),3,'omitnan'),2,'omitnan'),1,'omitnan')./...
    sum(sum(sum(GOBAI.vol(:,:,int2),3,'omitnan'),2,'omitnan'),1,'omitnan'));
temp_mean_lower = squeeze(sum(sum(sum(RG_ann.temp(:,:,int3,:).*...
    GOBAI.vol(:,:,int3),3,'omitnan'),2,'omitnan'),1,'omitnan')./...
    sum(sum(sum(GOBAI.vol(:,:,int3),3,'omitnan'),2,'omitnan'),1,'omitnan'));
temp_mean_depth = squeeze(sum(sum(RG_ann.temp.*...
    GOBAI.vol,2,'omitnan'),1,'omitnan')./...
    sum(sum(GOBAI.vol,2,'omitnan'),1,'omitnan'));

%% plot weighted means of temperature
ax3 = axes; hold on; box on; grid off;
set(ax3,'units','normalized','position',[0.37 0.56 0.25 0.42]);
set(ax3,'fontsize',18,'xticklabels',[],'linewidth',2);
ylabel(['Temperature Anomaly (' char(176) 'C)']);
ylim([-0.2 0.2]);
xlim([0 20]);
xticks([2 7 12 17]);
xticklabels({'2005' '2010' '2015' '2020'});
p1=plot(1:tdim,temp_mean_upper-mean(temp_mean_upper),'linewidth',5,'color',clrs_t(8,:));
scatter(1:tdim,temp_mean_upper-mean(temp_mean_upper),100,'o','MarkerFaceColor',clrs_t(8,:),'MarkerEdgeColor','none');
p2=plot(1:tdim,temp_mean_middle-mean(temp_mean_middle),'linewidth',5,'color',clrs_t(6,:));
scatter(1:tdim,temp_mean_middle-mean(temp_mean_middle),100,'o','MarkerFaceColor',clrs_t(6,:),'MarkerEdgeColor','none');
p3=plot(1:tdim,temp_mean_lower-mean(temp_mean_lower),'linewidth',5,'color',clrs_t(4,:));
scatter(1:tdim,temp_mean_lower-mean(temp_mean_lower),100,'o','MarkerFaceColor',clrs_t(4,:),'MarkerEdgeColor','none');
text(1,0.18,'b','fontsize',40);

%% add legend
legend([p1 p2 p3],{'0-100 dbars' '100-600 dbars' '600-2000 dbars'},...
    'location','southeast','NumColumns',1);

%% plot weighted means of temperature with depth
ax4 = axes; hold on; box on;
set(ax4,'units','normalized','position',[0.38 0.05 0.23 0.42]);
set(ax4,'fontsize',18,'xticklabels',[],'yticklabels',[],'YDir','reverse');
set(ax4,'layer','top','Box','on','linewidth',2);
contourf(1:tdim,GOBAI.pres,temp_mean_depth-mean(temp_mean_depth,2),...
    -0.2125:0.025:0.2125,'LineStyle','none');
ylim([0 1975]);
xlim([1 19]);
xticks([2 7 12 17]);
xticklabels({'2005' '2010' '2015' '2020'});
caxis([-0.2125 0.2125]);
mycolormap = cmocean('balance',17,'pivot',0);
mycolormap(9,:) = 1;
ax4.Colormap = mycolormap;
c=colorbar;
c.Label.String = ['Temperature Anomaly (' char(176) 'C)'];
c.Label.FontSize = 18;
c.TickLength = 0;
c.Position(3) = c.Position(3).*2;
c.Position(1) = c.Position(1)+0.03;
text(-0.5,0,'e','fontsize',40);

%% calculate weighted means of saturation
mean_upper_sat = squeeze(sum(sum(sum(RG_ann.oxy_sat(:,:,int1,:).*...
    GOBAI.vol(:,:,int1),3,'omitnan'),2,'omitnan'),1,'omitnan')./...
    sum(sum(sum(GOBAI.vol(:,:,int1),3,'omitnan'),2,'omitnan'),1,'omitnan'));
mean_middle_sat = squeeze(sum(sum(sum(RG_ann.oxy_sat(:,:,int2,:).*...
    GOBAI.vol(:,:,int2),3,'omitnan'),2,'omitnan'),1,'omitnan')./...
    sum(sum(sum(GOBAI.vol(:,:,int2),3,'omitnan'),2,'omitnan'),1,'omitnan'));
mean_lower_sat = squeeze(sum(sum(sum(RG_ann.oxy_sat(:,:,int3,:).*...
    GOBAI.vol(:,:,int3),3,'omitnan'),2,'omitnan'),1,'omitnan')./...
    sum(sum(sum(GOBAI.vol(:,:,int3),3,'omitnan'),2,'omitnan'),1,'omitnan'));
mean_depth_sat = squeeze(sum(sum(RG_ann.oxy_sat.*...
    GOBAI.vol,2,'omitnan'),1,'omitnan')./...
    sum(sum(GOBAI.vol,2,'omitnan'),1,'omitnan'));

%% plot weighted means of saturation
ax5 = axes; hold on; box on; grid off;
set(ax5,'units','normalized','position',[0.69 0.56 0.25 0.42]);
set(ax5,'fontsize',18,'xticklabels',[],'linewidth',2);
ylabel('[O_{2}] Saturation Anomaly (\mumol kg^{-1})');
ylim([-1 1]);
xlim([0 20]);
xticks([2 7 12 17]);
xticklabels({'2005' '2010' '2015' '2020'});
p1=plot(1:tdim,mean_upper_sat-mean(mean_upper_sat),'linewidth',5,'color',clrs_s(8,:));
scatter(1:tdim,mean_upper_sat-mean(mean_upper_sat),100,'o','MarkerFaceColor',clrs_s(8,:),'MarkerEdgeColor','none');
p2=plot(1:tdim,mean_middle_sat-mean(mean_middle_sat),'linewidth',5,'color',clrs_s(6,:));
scatter(1:tdim,mean_middle_sat-mean(mean_middle_sat),100,'o','MarkerFaceColor',clrs_s(6,:),'MarkerEdgeColor','none');
p3=plot(1:tdim,mean_lower_sat-mean(mean_lower_sat),'linewidth',5,'color',clrs_s(4,:));
scatter(1:tdim,mean_lower_sat-mean(mean_lower_sat),100,'o','MarkerFaceColor',clrs_s(4,:),'MarkerEdgeColor','none');
text(1,0.9,'c','fontsize',40);

%% add legend
legend([p1 p2 p3],{'0-100 dbars' '100-600 dbars' '600-2000 dbars'},...
    'location','southwest','NumColumns',1);

%% plot weighted means of saturation with depth
ax6 = axes; hold on; box on;
set(ax6,'units','normalized','position',[0.70 0.05 0.23 0.42]);
set(ax6,'fontsize',18,'xticklabels',[],'yticklabels',[],'YDir','reverse');
set(ax6,'layer','top','Box','on','linewidth',2);
contourf(1:tdim,GOBAI.pres,mean_depth_sat-mean(mean_depth_sat,2),...
    -1.0625:0.125:1.0625,'LineStyle','none');
ylim([0 1975]);
xlim([1 19]);
xticks([2 7 12 17]);
xticklabels({'2005' '2010' '2015' '2020'});
caxis([-1.0625 1.0625]);
mycolormap = cmocean('balance',17,'pivot',0);
mycolormap(9,:) = 1;
ax6.Colormap = mycolormap;
c=colorbar;
c.Label.String = '[O_{2}] Saturation Anomaly (\mumol kg^{-1})';
c.Label.FontSize = 18;
c.TickLength = 0;
c.Position(3) = c.Position(3).*2;
c.Position(1) = c.Position(1)+0.03;
text(-0.5,0,'f','fontsize',40);

%% export figure
exportgraphics(gcf,'/Users/sharp/Desktop/Figure7.png');
close

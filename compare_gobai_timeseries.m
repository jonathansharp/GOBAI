%% v1.0
% file information
ver = 'v1.0'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
figure; plot(double(datenum(0,0,0)+GOBAI.time),global_mean,'LineWidth',2); hold on

%% v2.0
% file information
ver = 'v2.0'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(datenum(0,0,0)+GOBAI.time),global_mean,'LineWidth',2);

%% v2.1
% file information
ver = 'v2.1'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(datenum(2004,0,0)+GOBAI.time),global_mean,'LineWidth',2);

%% v 2.2
% file information
ver = 'v2.2'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(datenum(1950,0,0)+GOBAI.time),global_mean,'LineWidth',2);

%% v 3.0
% file information
ver = 'v3.0'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxygen');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(GOBAI.time),global_mean,'LineWidth',2);

%% v 3.0-FFNN
% file information
ver = 'v3.0'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-FFNN-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-FFNN-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-FFNN-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-FFNN-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-FFNN-' ver '.nc'],'oxygen');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(GOBAI.time),global_mean,'LineWidth',2);

%% figure information
f = gcf;
f.Position(3) = f.Position(3)*2;
datetick('x','yyyy');
title('Global Mean GOBAI-O_{2} Versions');
ylabel('Weighted Average [O_{2}]');
legend({'v1.0' 'v2.0' 'v2.1' 'v2.2' 'v3.0' 'v3.0-FFNN'},'Location','northeast');

%% export figure
exportgraphics(gcf,'/raid/Data/GOBAI-O2/global_gobai_comparison.png');
close
clear

%% v 1.0
% file information
ver = 'v1.0'; % version
var = 'NO3'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download nitrate
GOBAI.nit = ncread([path 'GOBAI-' var '-' ver '.nc'],'nitrate');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.nit,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.nit,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(GOBAI.time),global_mean,'LineWidth',2);

%% figure information
f = gcf;
f.Position(3) = f.Position(3)*2;
datetick('x','yyyy');
title('Global Mean GOBAI-NO_{3} Versions');
ylabel('Weighted Average [NO_{3}^{2-}]');
legend({'v1.0'},'Location','northeast');

%% export figure
exportgraphics(gcf,'/raid/Data/GOBAI-NO3/global_gobai_comparison.png');
close
clear
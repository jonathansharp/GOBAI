%% v 2.3
% file information
ver = 'v2.3'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen uncertainty
GOBAI.uncer = ncread([path 'GOBAI-' var '-' ver '.nc'],'uncer');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.uncer,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.uncer,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(datenum(1950,0,0)+double(GOBAI.time),global_mean,'LineWidth',4);

%% v 1.0-HR-FFNN
% % file information
% ver = 'v1.0'; % version
% var = 'O2'; % variable
% path = ['/raid/Data/GOBAI-' var '/' ver '-HR/']; % file path
% % download dimensions
% GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],'lon');
% GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],'lat');
% GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],'pres');
% GOBAI.time = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],'time');
% % calculate weights
% GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
% mask = true(size(GOBAI.vol));
% for t = 1:length(GOBAI.time)
%     gobai_tmp = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],...
%         'o2',[1 1 1 t],[Inf Inf Inf 1]);
%     mask(isnan(gobai_tmp)) = false;
% end
% GOBAI.vol(~mask) = NaN;
% % download oxygen
% for t = 1:length(GOBAI.time)
%     gobai_tmp = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],...
%         'o2',[1 1 1 t],[Inf Inf Inf 1]);
%     gobai_tmp(~mask) = NaN;
%     global_mean(t) = sum(gobai_tmp(:).*...
%         GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
% end
% % calculate global mean
% plot(datenum(1950,0,0)+double(GOBAI.time),global_mean,'LineWidth',2);
% save([path 'mask'],'mask');
% save([path 'global_mean'],'global_mean');

%% figure information
f = gcf;
f.Position(3) = f.Position(3)*2;
set(gca,'fontsize',20);
datetick('x','yyyy');
title('Global');
ylabel('Weighted Average [O_{2}] Uncertainty');

%% export figure
exportgraphics(gcf,'/raid/Data/GOBAI-O2/global_gobai_uncer_comparison.png');
close
clear

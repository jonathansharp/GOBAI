%% create old mask
% file information
ver = 'v2.3'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
mask_original = any(GOBAI.oxy,4);

%% create new mask
% file information
ver = 'v1.2'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/HR-' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'longitude');
GOBAI.lat = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'latitude');
GOBAI.pres = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'mean_pressure');
GOBAI.time = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'time');
% calculate weights
mask_new = true(length(GOBAI.lon),length(GOBAI.lat),length(GOBAI.pres));
for t = 1:length(GOBAI.time)
    gobai_tmp = ...
        ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],...
        'o2',[1 1 1 t],[Inf Inf Inf 1]);
    mask_new(isnan(gobai_tmp)) = false;
end

%% reformat masks
plt_mask_original = mask_original & [mask_new(21:end,26:170,:);mask_new(1:20,26:170,:)];
plt_mask_new = mask_new & [false(360,25,58),[mask_original(341:end,:,:);mask_original(1:340,:,:)],false(360,10,58)];
clear GOBAI gobai_tmp mask_original mask_new t

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
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(~plt_mask_original) = NaN;
% calculate global mean
GOBAI.(['global_mean_' replace(ver,'.','_')]) = [];
for t = 1:length(GOBAI.time)
     gobai_tmp = ...
        ncread([path 'GOBAI-' var '-' ver '.nc'],...
        'oxy',[1 1 1 t],[Inf Inf Inf 1]);
    gobai_tmp(~plt_mask_original) = NaN;
    GOBAI.(['global_mean_' replace(ver,'.','_')])(t) = ...
        sum(gobai_tmp(:).*GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
end
figure;
% fill([datenum(0,0,0)+double(GOBAI.time);flipud(datenum(0,0,0)+double(GOBAI.time))],...
%     [global_mean+global_uncer,fliplr(global_mean-global_uncer)],'k','FaceAlpha',0.2);
plot(double(datenum(0,0,0)+GOBAI.time),GOBAI.(['global_mean_' replace(ver,'.','_')]),'LineWidth',2); hold on

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
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(~plt_mask_original) = NaN;
% calculate global mean
GOBAI.(['global_mean_' replace(ver,'.','_')]) = [];
for t = 1:length(GOBAI.time)
     gobai_tmp = ...
        ncread([path 'GOBAI-' var '-' ver '.nc'],...
        'oxy',[1 1 1 t],[Inf Inf Inf 1]);
    gobai_tmp(~plt_mask_original) = NaN;
    GOBAI.(['global_mean_' replace(ver,'.','_')])(t) = ...
        sum(gobai_tmp(:).*GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
end
% fill([datenum(0,0,0)+double(GOBAI.time);flipud(datenum(0,0,0)+double(GOBAI.time))],...
%     [global_mean+global_uncer,fliplr(global_mean-global_uncer)],'k','FaceAlpha',0.2);
plot(double(datenum(0,0,0)+GOBAI.time),GOBAI.(['global_mean_' replace(ver,'.','_')]),'LineWidth',2);

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
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(~plt_mask_original) = NaN;
% calculate global mean
GOBAI.(['global_mean_' replace(ver,'.','_')]) = [];
for t = 1:length(GOBAI.time)
     gobai_tmp = ...
        ncread([path 'GOBAI-' var '-' ver '.nc'],...
        'oxy',[1 1 1 t],[Inf Inf Inf 1]);
    gobai_tmp(~plt_mask_original) = NaN;
    GOBAI.(['global_mean_' replace(ver,'.','_')])(t) = ...
        sum(gobai_tmp(:).*GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
end
% fill([datenum(2004,0,0)+double(GOBAI.time);flipud(datenum(2004,0,0)+double(GOBAI.time))],...
%     [global_mean+global_uncer,fliplr(global_mean-global_uncer)],'k','FaceAlpha',0.2);
plot(double(datenum(2004,0,0)+GOBAI.time),GOBAI.(['global_mean_' replace(ver,'.','_')]),'LineWidth',2);

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
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(~plt_mask_original) = NaN;
% calculate global mean
GOBAI.(['global_mean_' replace(ver,'.','_')]) = [];
for t = 1:length(GOBAI.time)
     gobai_tmp = ...
        ncread([path 'GOBAI-' var '-' ver '.nc'],...
        'oxy',[1 1 1 t],[Inf Inf Inf 1]);
    gobai_tmp(~plt_mask_original) = NaN;
    GOBAI.(['global_mean_' replace(ver,'.','_')])(t) = ...
        sum(gobai_tmp(:).*GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
end
% fill([datenum(1950,0,0)+double(GOBAI.time);flipud(datenum(1950,0,0)+double(GOBAI.time))],...
%     [global_mean+global_uncer,fliplr(global_mean-global_uncer)],'k','FaceAlpha',0.2);
plot(double(datenum(1950,0,0)+GOBAI.time),GOBAI.(['global_mean_' replace(ver,'.','_')]),'LineWidth',2);

%% v 2.3-prelim
% % file information
% ver = 'v2.3'; % version
% var = 'O2'; % variable
% path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% % download dimensions
% GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '-prelim.nc'],'lon');
% GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '-prelim.nc'],'lat');
% GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '-prelim.nc'],'pres');
% GOBAI.time = ncread([path 'GOBAI-' var '-' ver '-prelim.nc'],'time');
% % download oxygen
% GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '-prelim.nc'],'oxy');
% GOBAI.uncer = ncread([path 'GOBAI-' var '-' ver '.nc'],'uncer');
% % calculate weights
% GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
% GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% % calculate global mean
% global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
%     length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
%     GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
% global_uncer = sqrt(sum(reshape(GOBAI.uncer.^2,[length(GOBAI.lon)*...
%     length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
%     GOBAI.vol(:).^2,'omitnan')./(sum(GOBAI.vol(:).^2,'omitnan')));
% fill([datenum(1950,0,0)+double(GOBAI.time);flipud(datenum(1950,0,0)+double(GOBAI.time))],...
%     [global_mean+global_uncer,fliplr(global_mean-global_uncer)],'k','FaceAlpha',0.2);
% plot(datenum(1950,0,0)+double(GOBAI.time),global_mean,'LineWidth',2);

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
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(~plt_mask_original) = NaN;
% calculate global mean
GOBAI.(['global_mean_' replace(ver,'.','_')]) = [];
for t = 1:length(GOBAI.time)
     gobai_tmp = ...
        ncread([path 'GOBAI-' var '-' ver '.nc'],...
        'oxy',[1 1 1 t],[Inf Inf Inf 1]);
    gobai_tmp(~plt_mask_original) = NaN;
    GOBAI.(['global_mean_' replace(ver,'.','_')])(t) = ...
        sum(gobai_tmp(:).*GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
end
% fill([datenum(1950,0,0)+double(GOBAI.time);flipud(datenum(1950,0,0)+double(GOBAI.time))],...
%     [global_mean+global_uncer,fliplr(global_mean-global_uncer)],'k','FaceAlpha',0.2);
plot(datenum(1950,0,0)+double(GOBAI.time),GOBAI.(['global_mean_' replace(ver,'.','_')]),'LineWidth',2);

%% v 2.4
% % file information
% ver = 'v2.4'; % version
% var = 'O2'; % variable
% path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% % download dimensions
% GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
% GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
% GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
% GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% % download oxygen
% GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'o2');
% % GOBAI.uncer = ncread([path 'GOBAI-' var '-' ver '.nc'],'uncer');
% % calculate weights
% GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
% GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% % calculate global mean
% global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
%     length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
%     GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
% % global_uncer = sqrt(sum(reshape(GOBAI.uncer.^2,[length(GOBAI.lon)*...
% %     length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
% %     GOBAI.vol(:).^2,'omitnan')./(sum(GOBAI.vol(:).^2,'omitnan')));
% % fill([datenum(1950,0,0)+double(GOBAI.time);flipud(datenum(1950,0,0)+double(GOBAI.time))],...
% %     [global_mean+global_uncer,fliplr(global_mean-global_uncer)],'k','FaceAlpha',0.2);
% plot(datenum(1950,0,0)+double(GOBAI.time),global_mean,'LineWidth',2);

%% v1.1-HR
% file information
ver = 'v1.1'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/HR-' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'longitude');
GOBAI.lat = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'latitude');
GOBAI.pres = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'mean_pressure');
GOBAI.time = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'time');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(~plt_mask_new) = NaN;
% download oxygen
GOBAI.(['global_mean_HR_' replace(ver,'.','_')]) = [];
for t = 1:length(GOBAI.time)
    gobai_tmp = ...
        ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],...
        'o2',[1 1 1 t],[Inf Inf Inf 1]);
    gobai_tmp(~plt_mask_new) = NaN;
    GOBAI.(['global_mean_HR_' replace(ver,'.','_')])(t) = ...
        sum(gobai_tmp(:).*GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
end
% calculate global mean
plot(datenum(1950,0,0)+double(GOBAI.time),GOBAI.(['global_mean_HR_' replace(ver,'.','_')]),'LineWidth',2);
% save([path 'global_mean_HR_' ver],['GOBAI.global_mean_HR_' replace(ver,'.','_')]);

%% v1.2-HR
% file information
ver = 'v1.2'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/HR-' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'longitude');
GOBAI.lat = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'latitude');
GOBAI.pres = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'mean_pressure');
GOBAI.time = ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],'time');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(~plt_mask_new) = NaN;
% download oxygen
GOBAI.(['global_mean_HR_' replace(ver,'.','_')]) = [];
for t = 1:length(GOBAI.time)
    gobai_tmp = ...
        ncread([path 'GOBAI-' var '-HR-' ver '-Monthly-Mean-1x1.nc'],...
        'o2',[1 1 1 t],[Inf Inf Inf 1]);
    gobai_tmp(~plt_mask_new) = NaN;
    GOBAI.(['global_mean_HR_' replace(ver,'.','_')])(t) = ...
        sum(gobai_tmp(:).*GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
end
% calculate global mean
plot(datenum(1950,0,0)+double(GOBAI.time),GOBAI.(['global_mean_HR_' replace(ver,'.','_')]),'LineWidth',2);
% save([path 'global_mean_HR_' ver],['GOBAI.global_mean_HR_' replace(ver,'.','_')]);

%% v1.2-HR
% file information
ver = 'v1.3'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/HR-' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-Monthly-Mean-1x1-' var '-HR-v202606.nc'],'longitude');
GOBAI.lat = ncread([path 'GOBAI-Monthly-Mean-1x1-' var '-HR-v202606.nc'],'latitude');
GOBAI.pres = ncread([path 'GOBAI-Monthly-Mean-1x1-' var '-HR-v202606.nc'],'mean_pressure');
GOBAI.time = ncread([path 'GOBAI-Monthly-Mean-1x1-' var '-HR-v202606.nc'],'time');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(~plt_mask_new) = NaN;
% download oxygen
GOBAI.(['global_mean_HR_' replace(ver,'.','_')]) = [];
for t = 1:length(GOBAI.time)
    gobai_tmp = ...
        ncread([path 'GOBAI-Monthly-Mean-1x1-' var '-HR-v202606.nc'],...
        'o2',[1 1 1 t],[Inf Inf Inf 1]);
    gobai_tmp(~plt_mask_new) = NaN;
    GOBAI.(['global_mean_HR_' replace(ver,'.','_')])(t) = ...
        sum(gobai_tmp(:).*GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
end
% calculate global mean
plot(datenum(1950,0,0)+double(GOBAI.time),GOBAI.(['global_mean_HR_' replace(ver,'.','_')]),'LineWidth',2);
% save([path 'global_mean_HR_' ver],['GOBAI.global_mean_HR_' replace(ver,'.','_')]);

%% figure information
f = gcf;
f.Position(3) = f.Position(3)*2;
datetick('x','yyyy');
title('Global Mean GOBAI-O_{2} Versions');
ylabel('Weighted Average [O_{2}]');
% legend({'v1.0' 'v2.0' 'v2.1' 'v2.2' 'v2.3-prelim' 'v2.3'},'Location','northeast');
% legend({'v1.0' 'v2.0' 'v2.1' 'v2.2' 'v2.3' 'v3.0' 'v3.0-FFNN'},'Location','northeast');
legend({'v1.0' 'v2.0' 'v2.1' 'v2.2' 'v2.3' 'HR-v1.1' 'HR-v1.2' 'HR-v1.3'},'Location','northeast');

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
GOBAI.nit = ncread([path 'GOBAI-' var '-' ver '.nc'],'no3');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.nit,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.nit,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(GOBAI.time),global_mean,'LineWidth',2);

%% v1.1-HR
% file information
ver = 'v1.1-HR'; % version
var = 'NO3'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '-Monthly-Mean-1x1.nc'],'longitude');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '-Monthly-Mean-1x1.nc'],'latitude');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '-Monthly-Mean-1x1.nc'],'mean_pressure');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '-Monthly-Mean-1x1.nc'],'time');
% download nitrate
GOBAI.nit = ncread([path 'GOBAI-' var '-' ver '-Monthly-Mean-1x1.nc'],'no3');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.nit,4,'omitnan'))) = NaN;
% determine mask (time varying coverage) and blank out
mask = true(size(GOBAI.vol));
for t = 1:length(GOBAI.time); mask(isnan(GOBAI.nit(:,:,:,t))) = false; end
for t = 1:length(GOBAI.time)
    gobai_tmp = GOBAI.nit(:,:,:,t);
    gobai_tmp(~mask) = NaN;
    GOBAI.nit(:,:,:,t) = gobai_tmp;
end
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

%% v1.1-HR
% file information
ver = 'HR-v1.1'; % version
var = 'DIC'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '-Monthly-Mean-1x1.nc'],'longitude');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '-Monthly-Mean-1x1.nc'],'latitude');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '-Monthly-Mean-1x1.nc'],'mean_pressure');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '-Monthly-Mean-1x1.nc'],'time');
% download nitrate
GOBAI.dic = ncread([path 'GOBAI-' var '-' ver '-Monthly-Mean-1x1.nc'],'dic');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.dic,4,'omitnan'))) = NaN;
% determine mask (time varying coverage) and blank out
mask = true(size(GOBAI.vol));
for t = 1:length(GOBAI.time); mask(isnan(GOBAI.dic(:,:,:,t))) = false; end
for t = 1:length(GOBAI.time)
    gobai_tmp = GOBAI.dic(:,:,:,t);
    gobai_tmp(~mask) = NaN;
    GOBAI.dic(:,:,:,t) = gobai_tmp;
end
% calculate global mean
global_mean = sum(reshape(GOBAI.dic,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(GOBAI.time),global_mean,'LineWidth',2);

%% figure information
f = gcf;
f.Position(3) = f.Position(3)*2;
datetick('x','yyyy');
title('Global Mean GOBAI-DIC Versions');
ylabel('Weighted Average [NO_{3}^{2-}]');
legend({'v1.0'},'Location','northeast');

%% export figure
exportgraphics(gcf,'/raid/Data/GOBAI-NO3/global_gobai_comparison.png');
close
clear
%% load and process basin masks

% define RECCAP2 mask file
filen = 'RECCAP2_region_masks_all_v20221025.nc';
% load and process RECCAP2 mask latitude
mask.lat = ncread(filen,'lat');
mask.lat = mask.lat(26:170);
% load and process RECCAP2 mask longitude
mask.lon = ncread(filen,'lon');
mask.lon = convert_lon(mask.lon);
mask.lon = [mask.lon(21:end);mask.lon(1:20)];
% load RECCAP2 basin masks
basins = {'atlantic' 'pacific' 'indian' 'arctic' 'southern'};
basin_names = {'Atlantic' 'Pacific' 'Indian' 'Arctic' 'Southern'};
for b = 1:length(basins)
    mask.(basins{b}) = ncread(filen,basins{b});
    mask.(basins{b}) = [mask.(basins{b})(21:end,26:170);mask.(basins{b})(1:20,26:170)];
end
% test plot
% figure;
% worldmap([-90 90],[0 360]);
% pcolorm(lat,lon,mask.(basins{2})');
% colorbar;
% clean up
clear filen

%%
% for b = 1:length(basins)
% 
% %% v1.0
% % file information
% ver = 'v1.0'; % version
% var = 'O2'; % variable
% path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% % download dimensions
% GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
% GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
% GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
% GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% % download oxygen
% GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% % calculate weights
% GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
% GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% % process basin
% mask_3d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres));
% mask_4d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres),length(GOBAI.time));
% GOBAI.oxy(~mask_4d) = NaN;
% GOBAI.vol(~mask_3d) = NaN;
% % calculate global mean
% basin_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
%     length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
%     GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
% figure; plot(double(datenum(0,0,0)+GOBAI.time),basin_mean,'LineWidth',2); hold on
% 
% %% v2.0
% % file information
% ver = 'v2.0'; % version
% var = 'O2'; % variable
% path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% % download dimensions
% GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
% GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
% GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
% GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% % download oxygen
% GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% % calculate weights
% GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
% GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% % process basin
% mask_3d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres));
% mask_4d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres),length(GOBAI.time));
% GOBAI.oxy(~mask_4d) = NaN;
% GOBAI.vol(~mask_3d) = NaN;
% % calculate global mean
% basin_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
%     length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
%     GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
% plot(double(datenum(0,0,0)+GOBAI.time),basin_mean,'LineWidth',2);
% 
% %% v2.1
% % file information
% ver = 'v2.1'; % version
% var = 'O2'; % variable
% path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% % download dimensions
% GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
% GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
% GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
% GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% % download oxygen
% GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% % calculate weights
% GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
% GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% % process basin
% mask_3d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres));
% mask_4d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres),length(GOBAI.time));
% GOBAI.oxy(~mask_4d) = NaN;
% GOBAI.vol(~mask_3d) = NaN;
% % calculate global mean
% basin_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
%     length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
%     GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
% plot(double(datenum(2004,0,0)+GOBAI.time),basin_mean,'LineWidth',2);
% 
% %% v 2.2
% % file information
% ver = 'v2.2'; % version
% var = 'O2'; % variable
% path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% % download dimensions
% GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
% GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
% GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
% GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% % download oxygen
% GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% % calculate weights
% GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
% GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% % process basin
% mask_3d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres));
% mask_4d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres),length(GOBAI.time));
% GOBAI.oxy(~mask_4d) = NaN;
% GOBAI.vol(~mask_3d) = NaN;
% % calculate global mean
% basin_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
%     length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
%     GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
% plot(double(datenum(1950,0,0)+GOBAI.time),basin_mean,'LineWidth',2);
% 
% %% v 3.0
% % file information
% ver = 'v3.0'; % version
% var = 'O2'; % variable
% path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% % download dimensions
% GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
% GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
% GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
% GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% % download oxygen
% GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxygen');
% % calculate weights
% GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
% GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% % process basin
% mask_3d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres));
% mask_4d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres),length(GOBAI.time));
% GOBAI.oxy(~mask_4d) = NaN;
% GOBAI.vol(~mask_3d) = NaN;
% % calculate global mean
% basin_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
%     length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
%     GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
% plot(double(GOBAI.time),basin_mean,'LineWidth',2);
% 
% %% figure information
% f = gcf;
% f.Position(3) = f.Position(3)*2;
% datetick('x','yyyy');
% title([basin_names{b} ' Mean GOBAI-O_{2} Versions']);
% ylabel('Weighted Average [O_{2}]');
% legend({'v1.0' 'v2.0' 'v2.1' 'v2.2' 'v3.0'},'Location','northeast');
% 
% %% export figure
% exportgraphics(gcf,['/raid/Data/GOBAI-O2/' basins{b} '_gobai_comparison.png']);
% close
% 
% % clean up
% clear f GOBAI mask_3d mask_4d var ver path
% 
% end



for b = 1:length(basins)

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
% process basin
mask_3d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres));
mask_4d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres),length(GOBAI.time));
GOBAI.nit(~mask_4d) = NaN;
GOBAI.vol(~mask_3d) = NaN;
% calculate global mean
basin_mean = sum(reshape(GOBAI.nit,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(GOBAI.time),basin_mean,'LineWidth',2);

%% figure information
f = gcf;
f.Position(3) = f.Position(3)*2;
datetick('x','yyyy');
title([basin_names{b} ' Mean GOBAI-NO_{3} Versions']);
ylabel('Weighted Average [NO_{3}^{2-}]');
legend({'v1.0'},'Location','northeast');

%% export figure
exportgraphics(gcf,['/raid/Data/GOBAI-NO3/' basins{b} '_gobai_comparison.png']);
close

% clean up
clear f GOBAI mask_3d mask_4d var ver path

end
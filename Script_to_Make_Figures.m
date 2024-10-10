
%% load v2.2
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
GOBAI.temp = ncread([path 'GOBAI-' var '-' ver '.nc'],'temp');

%% determine dimension lengths
xdim = length(GOBAI.lon);
ydim = length(GOBAI.lat);
zdim = length(GOBAI.pres);
tdim = length(GOBAI.time);

%% set limits
latlim=[-90 90];
lonlim=[-180 180];
t1=[2004 1 1 0 0 0];
t2=datevec(date);
% add months
RG.month = (1:tdim)';

%% determine ocean mask
ocean_mask = sum(~isnan(GOBAI.temp),4) == tdim;
% figure; worldmap(latlim,[20 380]);
% pcolorm(double(RG.latitude),double(RG.longitude),double(ocean_mask(:,:,58)'))

%% deternine area, depth, and volume weights
Weights

%% calculate oxygen saturation, AOU, density, depth, mass
RG.oxy_sat = o2satv2b(RG.sal,RG.temp);
RG.AOU = RG.oxy_sat-RG.oxy_ENS;
if ~isfield(RG,'sal_abs') && ~isfield(RG,'temp_cns')
    RG.sal_abs = nan(size(RG.oxy_ENS));
    RG.temp_cns = nan(size(RG.oxy_ENS));
    for t = 1:length(RG.time)
        RG.sal_abs(:,:,:,t) = gsw_SA_from_SP(RG.sal(:,:,:,t),...
            repmat(permute(RG.pressure,[3 2 1]),xdim,ydim,1),...
            repmat(RG.longitude,1,ydim,zdim),...
            repmat(RG.latitude',xdim,1,zdim));
        RG.temp_cns(:,:,:,t) = gsw_CT_from_t(RG.sal_abs(:,:,:,t),...
            RG.temp(:,:,:,t),repmat(permute(RG.pressure,[3 2 1]),xdim,ydim,1));
    end
end
RG.dens = gsw_rho(RG.sal_abs,RG.temp_cns,... % [kg/m3]
    repmat(permute(RG.pressure,[3 2 1]),xdim,ydim,1,tdim));
RG.depth = ... % meters
    -gsw_z_from_p(repmat(RG.pressure,1,ydim),repmat(RG.latitude',zdim,1));
RG.kg = RG.dens.*repmat(volume_weights.*1e9,1,1,1,tdim);

%% determine annual and climatological means on pressure levels
annual
annual_uncer
climatology
climatology_uncer

%% calculate and display global trends

%% interpolate variables onto a consistent density grid
interpolateSigma
interpolateSigma_clim

%% determine annual and climatological means on potential density levels
annual_sig
% climatology_sig

%% import World Ocean Atlas data
importWOA
WOA.OXY_SAT = o2satv2b(WOA.SAL,WOA.TMP);

%% plot O2 on depth levels and difference from WOA18
plotO2;
plotDiff_WOA;

%% calculate anomaly from O2 seasonal cycle on depth levels
O2_anom

%% calculate total oxygen utilization
preformed_O2
RG.TOU = RG.AOU - RG.oxy_dis;

%% sub-sample and scale according to Ito et al. (2017) procedure
Ito_Scaling

%% interpolate variables onto a consistent depth grid
% interpolateDepth

%% save annual anomalies as netcdf


%% import mixed layer depth
importMLD
clear MLD
% annual mean MLD
RG.mld_mean = mean(RG.mld,3,'omitnan');
% annual winter max mld
RG.mld_max = nan(xdim,ydim,tdim/12);
for a = 1:tdim/12
    RG.mld_max(:,:,a) = max(RG.mld(:,:,(a-1)*12+1:(a-1)*12+12),[],3,'omitnan');
end

%% import surface chlorophyll and calculate euphotic depth
chl_and_zeu
clear CHL

%% calculate full oxygen column inventory
RG.O2_col = nan(xdim,ydim,tdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:tdim
            if sum(~isnan(RG.oxy_ENS(a,b,:,c))) > 0
                % calculate oxygen column inventory
                RG.O2_col(a,b,c) = sum(RG.oxy_ENS(a,b,:,c).*... % umol/kg
                            RG.kg(a,b,:,c)./... % umol
                            1e6,3,'omitnan')./ ... % mol
                            (area_weights(a,b).*1e6); % mol/m2
            else
                RG.O2_col(a,b,c) = NaN;
            end
        end
    end
end

%% calculate trends in oxygen column inventory
RG.O2_col_trend = nan(xdim,ydim);
for a = 1:xdim
    for b = 1:ydim
        if sum(~isnan(RG.O2_col(a,b,:))) > 0
            fit = polyfit(RG.month,squeeze(RG.O2_col(a,b,:)),1);
            RG.O2_col_trend(a,b) = fit(1).*12;
        else
            RG.O2_col_trend(a,b) = NaN;
        end
    end
end

%% calculate full oxygen column saturation inventory
RG.O2_col_sat = nan(xdim,ydim,tdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:tdim
            if sum(~isnan(RG.oxy_sat(a,b,:,c))) > 0
                % calculate oxygen column inventory
                RG.O2_col_sat(a,b,c) = sum(RG.oxy_sat(a,b,:,c).*... % umol/kg
                            RG.kg(a,b,:,c)./... % umol
                            1e6,3,'omitnan')./ ... % mol
                            (area_weights(a,b).*1e6); % mol/m2
            else
                RG.O2_col_sat(a,b,c) = NaN;
            end
        end
    end
end

%% calculate trends in oxygen column saturation inventory
RG.O2_col_sat_trend = nan(xdim,ydim);
for a = 1:xdim
    for b = 1:ydim
        if sum(~isnan(RG.O2_col_sat(a,b,:))) > 0
            fit = polyfit(RG.month,squeeze(RG.O2_col_sat(a,b,:)),1);
            RG.O2_col_sat_trend(a,b) = fit(1).*12;
        else
            RG.O2_col_sat_trend(a,b) = NaN;
        end
    end
end

%% calculate trends in oxygen concentration on depth levels
RG.oxy_ENS_trend = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim
            if sum(~isnan(RG.oxy_ENS(a,b,c,:))) > 0
                [~,~,x] = leastsq2(RG.month,...
                    squeeze(RG.oxy_ENS(a,b,c,:)),0,2,[6 12]);
                RG.oxy_ENS_trend(a,b,c) = x(2).*12;
            else
                RG.oxy_ENS_trend(a,b,c) = NaN;
            end
        end
    end
end

% Print trends
disp(['Global Avg. Trend = ' num2str(round(sum(sum(sum(RG.oxy_ENS_trend.*...
    volume_weights,1,'omitnan'),2,'omitnan'),3,'omitnan')./...
    sum(sum(sum(volume_weights,1,'omitnan'),2,'omitnan'),3,'omitnan'),2))...
    'umol/kg'])
disp(['100-1000 dbars Trend = ' num2str(round(sum(sum(sum(RG.oxy_ENS_trend(:,:,11:44).*...
    volume_weights(:,:,11:44),1,'omitnan'),2,'omitnan'),3,'omitnan')./...
    sum(sum(sum(volume_weights(:,:,11:44),1,'omitnan'),2,'omitnan'),3,'omitnan'),2))...
    'umol/kg'])

%% calculate trends in oxygen saturation on depth levels
RG.oxy_sat_trend = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim
            if sum(~isnan(RG.oxy_sat(a,b,c,:))) > 0
                fit = polyfit(RG.month,squeeze(RG.oxy_sat(a,b,c,:)),1);
                RG.oxy_sat_trend(a,b,c) = fit(1).*12;
            else
                RG.oxy_sat_trend(a,b,c) = NaN;
            end
        end
    end
end

%% calculate trends in oxygen concentration on sigma levels
RG_sigma.oxy_ENS_trend = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim_sigma
            if sum(~isnan(RG_sigma.oxy_ENS(a,b,c,:))) >= 20
                idx = ~isnan(RG_sigma.oxy_ENS(a,b,c,:));
                fit = polyfit(RG.month(idx),squeeze(RG_sigma.oxy_ENS(a,b,c,idx)),1);
                RG_sigma.oxy_ENS_trend(a,b,c) = fit(1).*12;
            else
                RG_sigma.oxy_ENS_trend(a,b,c) = NaN;
            end
        end
    end
end

%% calculate trends in AOU on depth levels
RG.aou_trend = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim
            if sum(~isnan(RG.AOU(a,b,c,:))) > 0
                fit = polyfit(RG.month,squeeze(RG.AOU(a,b,c,:)),1);
                RG.aou_trend(a,b,c) = fit(1).*12;
            else
                RG.aou_trend(a,b,c) = NaN;
            end
        end
    end
end

%% calculate trends in AOU on sigma levels
RG_sigma.aou_trend = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim_sigma
            if sum(~isnan(RG_sigma.aou(a,b,c,:))) >= 20
                idx = ~isnan(RG_sigma.aou(a,b,c,:));
                fit = polyfit(RG.month(idx),squeeze(RG_sigma.aou(a,b,c,idx)),1);
                RG_sigma.aou_trend(a,b,c) = fit(1).*12;
            else
                RG_sigma.aou_trend(a,b,c) = NaN;
            end
        end
    end
end

%% calculate trends in density anomaly on depth levels
RG.sigma0_trend = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim
            if sum(~isnan(RG.sigma0(a,b,c,:))) > 0
                fit = polyfit(RG.month,squeeze(RG.sigma0(a,b,c,:)),1);
                RG.sigma0_trend(a,b,c) = fit(1).*12;
            else
                RG.sigma0_trend(a,b,c) = NaN;
            end
        end
    end
end

%% calculate trends in temperature on depth levels
RG.temp_trend = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim
            if sum(~isnan(RG.temp(a,b,c,:))) > 0
                fit = polyfit(RG.month,squeeze(RG.temp(a,b,c,:)),1);
                RG.temp_trend(a,b,c) = fit(1).*12;
            else
                RG.temp_trend(a,b,c) = NaN;
            end
        end
    end
end

%% calculate trends in salinity on depth levels
RG.sal_trend = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim
            if sum(~isnan(RG.sal(a,b,c,:))) > 0
                fit = polyfit(RG.month,squeeze(RG.sal(a,b,c,:)),1);
                RG.sal_trend(a,b,c) = fit(1).*12;
            else
                RG.sal_trend(a,b,c) = NaN;
            end
        end
    end
end

%% calculate correlations with oxygen on depth levels
RG.oxy_temp_corr = nan(xdim,ydim,zdim);
RG.oxy_sal_corr = nan(xdim,ydim,zdim);
RG.oxy_sigma0_corr = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim
            if sum(~isnan(RG.oxy_ENS(a,b,c,:))) > 0
                RG.oxy_temp_corr(a,b,c) = ...
                    corr(squeeze(RG.oxy_ENS(a,b,c,:)),squeeze(RG.temp(a,b,c,:)));
                RG.oxy_sal_corr(a,b,c) = ...
                    corr(squeeze(RG.oxy_ENS(a,b,c,:)),squeeze(RG.sal(a,b,c,:)));
                RG.oxy_sigma0_corr(a,b,c) = ...
                    corr(squeeze(RG.oxy_ENS(a,b,c,:)),squeeze(RG.sigma0(a,b,c,:)));
            else
                RG.oxy_temp_corr(a,b,c) = NaN;
                RG.oxy_sal_corr(a,b,c) = NaN;
                RG.oxy_sigma0_corr(a,b,c) = NaN;
            end
        end
    end
end

%% calculate correlations with AOU on depth levels
RG.aou_temp_corr = nan(xdim,ydim,zdim);
RG.aou_sal_corr = nan(xdim,ydim,zdim);
RG.aou_sigma0_corr = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim
            if sum(~isnan(RG.oxy_ENS(a,b,c,:))) > 0
                RG.aou_temp_corr(a,b,c) = ...
                    corr(squeeze(RG.AOU(a,b,c,:)),squeeze(RG.temp(a,b,c,:)));
                RG.aou_sal_corr(a,b,c) = ...
                    corr(squeeze(RG.AOU(a,b,c,:)),squeeze(RG.sal(a,b,c,:)));
                RG.aou_sigma0_corr(a,b,c) = ...
                    corr(squeeze(RG.AOU(a,b,c,:)),squeeze(RG.sigma0(a,b,c,:)));
            else
                RG.aou_temp_corr(a,b,c) = NaN;
                RG.aou_sal_corr(a,b,c) = NaN;
                RG.aou_sigma0_corr(a,b,c) = NaN;
            end
        end
    end
end

%% calculate correlations between oxygen contributors on depth levels
RG.oxy_aou_corr = nan(xdim,ydim,zdim);
RG.oxy_sol_corr = nan(xdim,ydim,zdim);
RG.aou_sol_corr = nan(xdim,ydim,zdim);
for a = 1:xdim
    for b = 1:ydim
        for c = 1:zdim
            if sum(~isnan(RG.oxy_ENS(a,b,c,:))) > 0
                RG.oxy_aou_corr(a,b,c) = ...
                    corr(squeeze(RG.oxy_ENS(a,b,c,:)),squeeze(RG.AOU(a,b,c,:)));
                RG.oxy_sol_corr(a,b,c) = ...
                    corr(squeeze(RG.oxy_ENS(a,b,c,:)),squeeze(RG.oxy_sat(a,b,c,:)));
                RG.aou_sol_corr(a,b,c) = ...
                    corr(squeeze(RG.AOU(a,b,c,:)),squeeze(RG.oxy_sat(a,b,c,:)));
            else
                RG.oxy_aou_corr(a,b,c) = NaN;
                RG.oxy_sol_corr(a,b,c) = NaN;
                RG.aou_sol_corr(a,b,c) = NaN;
            end
        end
    end
end

%% determine maximum mixed layer differential between years
RG.mld_diff = nan(xdim,ydim,tdim/12-1);
for a = 1:xdim
    for b = 1:ydim
        for t = 1:tdim/12-1
            RG.mld_diff(a,b,t) = ...
                max(RG.mld(a,b,(t-1)*12+1:(t-1)*12+12))-...
                max(RG.mld(a,b,t*12+1:t*12+12));
        end
    end
end

%% plot MLD change (GIF)
h = figure;
set(h,'color','white');
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'Figures/GIFs/MLD_increase.gif';
year = 2005:2020; n=1;
for d = 1:tdim/12-1
    worldmap(latlim,lonlim);
    title(['Maximum Mixed Layer Depth Increase (' num2str(year(d)-1) ...
        ' to ' num2str(year(d)) ')'],'fontsize',14);
    setm(gca,'ffacecolor','w');
    setm(gca,'fontsize',12);
    pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
        double(repmat(RG.longitude,1,size(RG.latitude,1))),...
        RG.mld_diff(:,:,d));
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(land,'FaceColor',rgb('grey'));
    c=colorbar;
    caxis([-25 25]);
    cm = flipud(cmocean('balance',10,'pivot',0));
    cm(5:6,:) = 1; colormap(cm);
    c.FontSize = 12;
    c.Label.String = 'Max MLD Increase (m)';
    mlabel off; plabel off;
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1); 
    end
    n=n+1;
end

%% determine thickness of hypoxia
% define grid cells as hypoxic or not
RG.hypoxic_cells = false(size(RG.oxy_ENS));
RG.hypoxic_cells = RG.oxy_ENS < 62;

% determine thickness of each grid cell
RG.hypoxic_thickness = nan(size(RG.oxy_ENS));
RG.hypoxic_thickness(~isnan(RG.oxy_ENS)) = 0;
depth_weights_4D = repmat(depth_weights,1,1,1,tdim);
RG.hypoxic_thickness(RG.hypoxic_cells) = depth_weights_4D(RG.hypoxic_cells);
clear depth_weights_4D

% sum grid cells through depth
RG.hypoxic_thickness = squeeze(sum(RG.hypoxic_thickness,3));
RG.hypoxic_thickness(RG.hypoxic_thickness==0) = NaN;

%% plot annual mean hypoxic thickness
figure; worldmap(latlim,[20 380]);
title('Thickness (m) of hypoxia ([O_{2}] < 20 \mumol kg^{-1})','fontsize',16);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
    double(repmat(RG.longitude,1,size(RG.latitude,1))),...
    mean(RG.hypoxic_thickness,3,'omitnan'));
% Indian
lon_idx_Ind = find(RG.longitude == 70.5);
lat_idx_Ind = RG.latitude > -65 & RG.latitude < 15;
lat_plot_Ind = double(RG.latitude(lat_idx_Ind));
lon_plot_Ind = double(repmat(RG.longitude(lon_idx_Ind),...
    length(RG.latitude(lat_idx_Ind)),1));
plotm(lat_plot_Ind,lon_plot_Ind,'linewidth',3,'color','k');
% Pacific
lon_idx_Pac = find(RG.longitude == -163.5);
lat_idx_Pac = RG.latitude > -65 & RG.latitude < 55;
lat_plot_Pac = double(RG.latitude(lat_idx_Pac));
lon_plot_Pac = double(repmat(RG.longitude(lon_idx_Pac),...
    length(RG.latitude(lat_idx_Pac)),1));
plotm(lat_plot_Pac,lon_plot_Pac,'linewidth',3,'color','k');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
caxis([0 2000]);
colormap(cmocean('deep',20));
c.FontSize = 12;
mlabel off; plabel off;
exportgraphics(gcf,'/Users/sharp/Desktop/hypoxic_thickness_20uM.jpg');

%% calculate trends in hypoxia
RG.hypoxic_trend = nan(xdim,ydim);
RG.hypoxic_trend_err = nan(xdim,ydim);
for a = 1:xdim
    for b = 1:ydim
        if sum(~isnan(RG.hypoxic_thickness(a,b,:))) > 0
            [~,yr,x,err] = leastsq2(RG.month,...
                squeeze(RG.hypoxic_thickness(a,b,:)),...
                0,0,0);
            % scale from monthly to annual
            RG.hypoxic_trend(a,b) = x(2).*12;
            [acov,acor,lag,dof] = autocov2(RG.month,yr,36);
            RG.hypoxic_trend_err(a,b) = err(2).*12.*... % annual standard error
                (sqrt(length(RG.month))./sqrt(dof)).*... % scaled by degrees of freedom
                2; % scaled to 95% confidence (estimated scaling as 2)
        else
            RG.hypoxic_trend(a,b) = NaN;
            RG.hypoxic_trend_err(a,b) = NaN;
        end
    end
end

% plot trend in hypoxic thickness
figure; worldmap(latlim,[20 380]); hold on;
%set(gcf,'units','normalized','position',[0 0 1 1]);
title('Thickness of suboxia trend (m yr^{-1})','fontsize',16);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
    double(repmat(RG.longitude,1,size(RG.latitude,1))),...
    RG.hypoxic_trend);
idx = (RG.hypoxic_trend+RG.hypoxic_trend_err > 0 & ...
       RG.hypoxic_trend-RG.hypoxic_trend_err > 0) | ...
      (RG.hypoxic_trend+RG.hypoxic_trend_err < 0 & ...
       RG.hypoxic_trend-RG.hypoxic_trend_err < 0);
idx = ~idx; idx(RG.hypoxic_trend == 0) = 0; idx = idx(:);
lat_tmp = double(repmat(RG.latitude',size(RG.longitude,1),1));
lon_tmp = double(repmat(RG.longitude,1,size(RG.latitude,1)));
%scatterm(lat_tmp(idx),lon_tmp(idx),1,'ok');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
caxis([-10 10]);
cm = cmocean('balance',200,'pivot',0); cm(100:101,:) = 1; colormap(cm);
% colormap(cmocean('balance',200,'pivot',0));
c.FontSize = 12;
mlabel off; plabel off;
exportgraphics(gcf,'/Users/sharp/Desktop/hypoxic_thickness_trend_10uM.jpg');


%% determine top and bottom of hypoxia

% determine top and bottom of hypoxia in each month for each lat/lon
RG.hypoxic_top = nan(xdim,ydim,tdim);
RG.hypoxic_bot = nan(xdim,ydim,tdim);
for a = 1:xdim
    for b = 1:ydim
        for d = 1:tdim
            if sum(RG.hypoxic_cells(a,b,:,d)) > 0
                idx_top = ...
                    find(squeeze(RG.hypoxic_cells(a,b,:,d)),1,'first');
                RG.hypoxic_top(a,b,d) = RG.pressure(idx_top);
                idx_bot = ...
                    find(squeeze(RG.hypoxic_cells(a,b,:,d)),1,'last');
                RG.hypoxic_bot(a,b,d) = RG.pressure(idx_bot);
            else
                RG.hypoxic_top(a,b,d) = NaN;
                RG.hypoxic_bot(a,b,d) = NaN;
            end
        end
    end
end

% plot average top of hypoxia
figure; worldmap(latlim,[20 380]);
title('Annual mean top depth (m) of hypoxia ([O_{2}] < 62 \mumol kg^{-1})','fontsize',16);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
    double(repmat(RG.longitude,1,size(RG.latitude,1))),...
    mean(RG.hypoxic_top,3));
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
caxis([0 1000]);
colormap(cmocean('thermal',10));
set(c,'YDir','reverse');
c.FontSize = 12;
mlabel off; plabel off;
exportgraphics(gcf,'/Users/sharp/Desktop/hypoxic_top.jpg');

% plot average bottom of hypoxia
figure; worldmap(latlim,[20 380]);
title('Annual mean bottom depth (m) of hypoxia ([O_{2}] < 62 \mumol kg^{-1})','fontsize',16);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
    double(repmat(RG.longitude,1,size(RG.latitude,1))),...
    mean(RG.hypoxic_bot,3));
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
caxis([500 2000]);
colormap(cmocean('thermal',15));
set(c,'YDir','reverse');
c.FontSize = 12;
mlabel off; plabel off;
exportgraphics(gcf,'/Users/sharp/Desktop/hypoxic_bot.jpg');

%% determine trend in top and bottom of hypoxia

% determine trend in top and bottom of hypoxia in each month for each lat/lon
for a = 1:xdim
    for b = 1:ydim
        if sum(isnan(RG.hypoxic_top(a,b,:))) < 102
            fit = polyfit(RG.month,squeeze(RG.hypoxic_top(a,b,:)),1);
            RG.hypoxic_top_trend(a,b) = fit(1).*12;
            fit = polyfit(RG.month,squeeze(RG.hypoxic_bot(a,b,:)),1);
            RG.hypoxic_bot_trend(a,b) = fit(1).*12;
        else
           RG.hypoxic_top_trend(a,b) = NaN;
           RG.hypoxic_bot_trend(a,b) = NaN;
        end
    end
end

% plot trend in top of hypoxia
figure; worldmap(latlim,[20 380]);
title('Trend in top depth of hypoxia (m yr^{-1})','fontsize',16);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
    double(repmat(RG.longitude,1,size(RG.latitude,1))),...
    RG.hypoxic_top_trend);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
caxis([-10 10]);
cm = flipud(cmocean('balance',200,'pivot',0)); cm(95:105,:) = 1; colormap(cm);
set(c,'YDir','reverse');
c.FontSize = 12;
mlabel off; plabel off;
exportgraphics(gcf,'/Users/sharp/Desktop/hypoxic_top_trend.jpg');

% plot trend in bottom of hypoxia
figure; worldmap(latlim,[20 380]);
title('Trend in bottom depth of hypoxia (m yr^{-1})','fontsize',16);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
    double(repmat(RG.longitude,1,size(RG.latitude,1))),...
    RG.hypoxic_bot_trend);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
caxis([-10 10]);
cm = cmocean('balance',200,'pivot',0); cm(95:105,:) = 1; colormap(cm);
set(c,'YDir','reverse');
c.FontSize = 12;
mlabel off; plabel off;
exportgraphics(gcf,'/Users/sharp/Desktop/hypoxic_bot_trend.jpg');

%% determine O2 beneath seasonal maximum MLD
RG.O2_beneath_MLD = nan(xdim,ydim,tdim);
for a = 1:xdim
    for b = 1:ydim
        for t = 1:tdim/12
            if sum(~isnan(RG.oxy_ENS(a,b,:,t))) > 0 && ...
                    sum(~isnan(RG.mld(a,b,(t-1)*12+1:(t-1)*12+12))) > 0
                maxMLD = max(RG.mld(a,b,(t-1)*12+1:(t-1)*12+12));
                index = find(abs(RG.pressure-maxMLD+40) == ...
                    min(abs(RG.pressure-maxMLD+40)));
                RG.O2_beneath_MLD(a,b,(t-1)*12+1:(t-1)*12+12) = ...
                    RG.oxy_ENS(a,b,index,(t-1)*12+1:(t-1)*12+12);
            else
                RG.O2_beneath_MLD(a,b,(t-1)*12+1:(t-1)*12+12) = NaN;
            end
        end
    end
end

RG.O2_beneath_MLD_clim = nan(xdim,ydim,12);
RG.O2_beneath_MLD_std = nan(xdim,ydim,12);
for a = 1:xdim
    for b = 1:ydim
        for m = 1:12
            RG.O2_beneath_MLD_clim(a,b,m) = ...
                mean(RG.O2_beneath_MLD(a,b,m:12:end));
            RG.O2_beneath_MLD_std(a,b,m) = ...
                std(RG.O2_beneath_MLD(a,b,m:12:end));
        end
    end
end

figure; set(gcf,'units','normalized','position',[0 0 0.5 1]);
lat = -70:20:70;
for b = 1:8
    subplot(8,1,9-b); hold on;
    latmin = find(abs(RG.latitude-lat(b)) == ...
        min(abs(RG.latitude-lat(b))),2);
    latmax = find(abs(RG.latitude-(lat(b)+20)) == ...
        min(abs(RG.latitude-(lat(b)+20))),1);
    o2 = squeeze(mean(mean(RG.O2_beneath_MLD_clim(:,latmin:latmax,:),...
    1,'omitnan'),2,'omitnan'));
    o2_std = squeeze(mean(mean(RG.O2_beneath_MLD_std(:,latmin:latmax,:),...
    1,'omitnan'),2,'omitnan'));
    errorbar([1:12]',o2,o2_std,'linewidth',2);
    scatter([1:12]',o2,'ok','filled');
    ylabel([num2str(lat(b)) ' N to ' num2str(lat(b)+20) ' N']);
    ylim([round(mean(o2)-20) round(mean(o2)+20)]);
    xlim([0 13]);
    if b > 1
        set(gca,'XTick',[]);
    else
        xticks(1:12);
        xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
    end
end
title('Climatological oxygen 40m beneath mixed layer in latitudinal segmets [\mumol kg^{-1}]','fontsize',16);
exportgraphics(gcf,'/Users/sharp/Desktop/o2_beneath_MLD.jpg');

%% determine monthly O2 inventory (above 2000m)
% convert umol/kg to mol/km3
O2_temp = RG.oxy_ENS./1e6.*RG.dens.*1e9;
O2_inv_2000 = nan(length(RG.time),1);
for m = 1:length(RG.time)
% convert mol/km3 to moles each year
O2_inv_2000(m) = ...
    double(mean(squeeze(sum(sum(sum(...
    (O2_temp(:,:,:,m).*volume_weights),...
    1,'omitnan'),2,'omitnan'),3,'omitnan')),'omitnan'));
end
O2_inv_2000_rm = smoothdata(O2_inv_2000,'movmean',12);
clear O2_temp

%% determine monthly hypoxic volume (above 2000m)


%% determine monthly ocean heat content (above 2000m)
RG.Cp = sw_cp(RG.sal,RG.temp,repmat(permute(RG.pressure,[3 2 1]),xdim,ydim,1,tdim));
RG.OHC = RG.dens.*RG.Cp.*RG.temp.*repmat(volume_weights.*1e9,1,1,1,tdim);
OHC_2000 = squeeze(sum(sum(sum(RG.OHC,1,'omitnan'),2,'omitnan'),3,'omitnan'));
OHC_2000_anom = squeeze(sum(sum(sum(RG.OHC,1,'omitnan'),2,'omitnan'),3,'omitnan')) - ...
    mean(sum(sum(sum(RG.OHC,1,'omitnan'),2,'omitnan'),3,'omitnan'));
OHC_2000_rm = smoothdata(OHC_2000,'movmean',12);
OHC_2000_anom_rm = smoothdata(OHC_2000_anom,'movmean',12);

%% determine annual mean O2 inventory (0-1000m)
% convert umol/kg to mol/km3
O2_temp = RG.oxy_ENS./1e6.*RG.dens.*1e9;
annual_O2_inv_1000 = nan(length(RG.time)/12,1);
for y = 1:length(RG.time)/12
% convert mol/km3 to moles each year
annual_O2_inv_1000(y) = ...
    double(mean(squeeze(sum(sum(sum(...
    (O2_temp(:,:,1:44,(y-1)*12+1:(y-1)*12+12).*...
    volume_weights(:,:,1:44)),...
    1,'omitnan'),2,'omitnan'),3,'omitnan')),'omitnan'));
end
clear O2_temp
annual_O2_inv_anom_1000 = (annual_O2_inv_1000 - mean(annual_O2_inv_1000))./10^14;

%% determine annual mean O2 saturation inventory (200-1000m)
% convert umol/kg to mol/km3
O2_temp = RG.oxy_sat./1e6.*RG.dens.*1e9;
annual_O2_sat_inv_1000 = nan(length(RG.time)/12,1);
for y = 1:length(RG.time)/12
% convert mol/km3 to moles each year
annual_O2_sat_inv_1000(y) = ...
    double(mean(squeeze(sum(sum(sum(...
    (O2_temp(:,:,44:58,(y-1)*12+1:(y-1)*12+12).*...
    volume_weights(:,:,44:58)),...
    1,'omitnan'),2,'omitnan'),3,'omitnan')),'omitnan'));
end
clear O2_temp
annual_O2_sat_inv_anom_1000 = (annual_O2_sat_inv_1000 - mean(annual_O2_sat_inv_1000))./10^14;

%% plot annual mean O2 inventory and o2 saturation inventory (200-1000m)
figure; hold on;
set(gcf,'units','normalized','position',[1 1 0.5 0.3]);
set(gca,'fontsize',16);
o2_slp = polyfit(2004:2021,annual_O2_inv_anom_1000,1);
o2_slp = o2_slp(1);
p1=plot(2004:2021,annual_O2_inv_anom_1000,'-r','linewidth',4);
scatter(2004:2021,annual_O2_inv_anom_1000,50,'ok','filled');
sat_slp = polyfit(2004:2021,annual_O2_sat_inv_anom_1000,1);
sat_slp = sat_slp(1);
p2=plot(2004:2021,annual_O2_sat_inv_anom_1000,'-b','linewidth',4);
scatter(2004:2021,annual_O2_sat_inv_anom_1000,50,'ok','filled');
plot([2000 2025],[0 0],'--k');
xlim([2003 2022]);
ylim([-6 4]);
legend([p1 p2],...
    {['[O_{2}] Anom. (' num2str(round(o2_slp,2)) '\mumol kg^{-1} yr^{-1})'] ...
    ['[O_{2}]_{sol.} Anom. (' num2str(round(sat_slp,2)) '\mumol kg^{-1} yr^{-1})']},...
    'location','southwest');
exportgraphics(gcf,'/Users/sharp/Desktop/o2_vs_sol_1000_2000.jpg');

%% determine annual mean O2 inventory (above 200m)
% convert umol/kg to mol/km3
O2_temp = RG.oxy_ENS./1e6.*RG.dens.*1e9;
annual_O2_inv = nan(tdim/12,1);
% index to depth interval
depth1 = find(RG.pressure == 200);
depth2 = find(RG.pressure == 1000);
for y = 1:tdim/12
% convert mol/km3 to moles each year
annual_O2_inv(y) = ...
    double(mean(squeeze(sum(sum(sum(...
    (O2_temp(:,:,depth1:depth2,(y-1)*12+1:(y-1)*12+12).*...
    volume_weights(:,:,depth1:depth2)),...
    1,'omitnan'),2,'omitnan'),3,'omitnan')),'omitnan'));
end
clear O2_temp
annual_O2_inv_anom = (annual_O2_inv - mean(annual_O2_inv))./10^14;

%% determine annual average O2 concentration anomaly (above 1000m)
annual_O2_mean = nan(length(RG.time)/12,1);
for y = 1:length(RG.time)/12
% convert mol/km3 to moles each year
weights_temp = volume_weights(:,:,1:44);
annual_O2_mean(y) = ...
    double(mean(squeeze(sum(sum(sum(...
    (RG.oxy_ENS(:,:,1:44,(y-1)*12+1:(y-1)*12+12).*...
    weights_temp)./...
    sum(sum(sum(weights_temp,...
    1,'omitnan'),2,'omitnan'),3,'omitnan'),...
    1,'omitnan'),2,'omitnan'),3,'omitnan'))));
end
clear weights_temp
annual_O2_anom = annual_O2_mean - mean(annual_O2_mean);

%% determine Ito et al. (2017) annual average O2 concentration anomaly (above 1000m)
% load Ito data
Ito = load_Ito;
% define weights of Ito dataset
[area_weights_Ito,depth_weights_Ito,volume_weights_Ito] = ...
    Weights_Ito(Ito);
annual_O2_anom_Ito = nan(length(Ito.time.data),1);
for y = 1:length(Ito.time.data)
% convert mol/km3 to moles each year
annual_O2_anom_Ito(y) = ...
    squeeze(sum(sum(sum(...
    (Ito.o2.data(:,:,:,y).*...
    volume_weights_Ito)./...
    sum(sum(sum(volume_weights_Ito,...
    1,'omitnan'),2,'omitnan'),3,'omitnan'),...
    1,'omitnan'),2,'omitnan'),3,'omitnan'));
end
year = datevec(datenum([repmat(1950,67,1) ones(67,1) ones(67,1) ...
    zeros(67,1) Ito.time.data zeros(67,1)])); year = year(:,1);


figure; plot(year,annual_O2_anom_Ito,unique(squeeze(RG.year(1,1,1,:))),annual_O2_anom)

%% plot O2 and AOU trends (2004-2020) spatially

type = 'depth';
minimum = 100;
maximum = 500;
tail = [type '_' num2str(minimum) '_' num2str(maximum)];

% calculate oxygen column inventory and AOU inventory
[O2_col,AOU_col] = ...
    col_O2_AOU(RG,type,minimum,maximum,xdim,ydim,tdim,area_weights);

% calculate trends in oxygen column inventory and AOU inventory
[O2_col_trend,AOU_col_trend] = ...
    col_trend_O2_AOU(RG.month,O2_col,AOU_col,xdim,ydim);

% plot oxygen column inventory and AOU inventory
plot_col_O2(RG,O2_col,type,minimum,maximum,latlim)
exportgraphics(gcf,['/Users/sharp/Desktop/o2_column_inv_' tail '.jpg']);
close
plot_col_AOU(RG,AOU_col,type,minimum,maximum,latlim)
exportgraphics(gcf,['/Users/sharp/Desktop/AOU_column_inv_' tail '.jpg']);
close

% plot trends in oxygen column inventory and AOU inventory
plot_col_trend_O2(RG,O2_col_trend,type,minimum,maximum,latlim)
exportgraphics(gcf,['/Users/sharp/Desktop/o2_column_inv_trend_' tail '.jpg']);
close
plot_col_trend_AOU(RG,AOU_col_trend,type,minimum,maximum,latlim)
exportgraphics(gcf,['/Users/sharp/Desktop/AOU_column_trend_' tail '.jpg']);
close

% plot depths of sigma level or sigma at depth level
plot_mean_depth_level(RG,type,minimum,maximum,latlim,xdim,ydim,tdim)
exportgraphics(gcf,['/Users/sharp/Desktop/mean_depth_levels_' tail '.jpg']);
close

% plot GIF of oxygen column inventory and AOU inventory

%% plot O2 on isopycnals at HOT
figure; hold on; title('Output at HOT');
scatter(RG.time,squeeze(RG_sigma.oxy_ENS(182,88,56,:)));
datetick('x','keeplimits');
ylabel('Oxygen on Sigma 26.5 (\mumol/kg)');
ylim([0 180]);
set(gca,'fontsize',15);
exportgraphics(gcf,'/Users/sharp/Desktop/HOT_oxy_sig_26.5.jpg');

figure; hold on; title('Output at HOT');
scatter(RG.time,squeeze(RG_sigma.oxy_ENS(182,88,61,:)));
datetick('x','keeplimits');
ylabel('Oxygen on Sigma 27.0 (\mumol/kg)');
ylim([0 70]);
set(gca,'fontsize',15);
exportgraphics(gcf,'/Users/sharp/Desktop/HOT_oxy_sig_27.0.jpg');

%% plot O2 on sigma levels
sigma = 27.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_idx = RG_sigma.sigma0 == sigma;
figure; worldmap(latlim,[20 380]);
title(['Oxygen on \sigma = ' num2str(round(sigma,1))],'fontsize',16);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
    double(repmat(RG.longitude,1,size(RG.latitude,1))),...
    mean(RG_sigma.oxy_ENS(:,:,sig_idx,:),4));
land = shaperead('landareas', 'UseGeoCoords',true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
caxis([0 350]);
colormap(cmocean('haline',14));
c.FontSize = 14;
c.Label.String = '[O_{2}] (\mumol kg^{-1})';
mlabel off; plabel off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Trend in O2 on sigma levels
RG_sigma.oxy_trend = nan(xdim,ydim,zdim_sigma);
for n = 1:zdim_sigma
    sig_idx = find(abs(RG_sigma.sigma0 - RG_sigma.sigma0(n)) == ...
        min(abs(RG_sigma.sigma0 - RG_sigma.sigma0(n))));
    % calculate trends
    for a = 1:xdim
        for b = 1:ydim
            if sum(~isnan(RG_sigma.oxy_ENS(a,b,sig_idx,:))) > 102
                [~,~,x,~,~,~,~] = ...
                    leastsq2(RG.month,squeeze(RG_sigma.oxy_ENS(a,b,sig_idx,:)),0,2,[6 12]);
                RG_sigma.oxy_trend(a,b,n) = x(2)*12;
            else
                RG_sigma.oxy_trend(a,b,n) = NaN;
            end
        end
    end
    % plot trends
    figure; worldmap(latlim,[20 380]);
    title(['Oxygen trend on \sigma = ' num2str(round(RG_sigma.sigma0(n),1))],'fontsize',16);
    setm(gca,'ffacecolor','w');
    setm(gca,'fontsize',12);
    pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
        double(repmat(RG.longitude,1,size(RG.latitude,1))),...
        RG_sigma.oxy_trend(:,:,n));
    land = shaperead('landareas', 'UseGeoCoords',true);
    geoshow(land,'FaceColor',rgb('grey'));
    c=colorbar;
    caxis([-2 2]);
    colormap(cmocean('balance',16,'pivot',0));
    c.FontSize = 14;
    c.Label.String = '\mumol kg^{-1} yr^{-1}';
    mlabel off; plabel off;
    exportgraphics(gcf,['/Users/sharp/Desktop/O2_sigma/sig' num2str(10*RG_sigma.sigma0(n)) '.jpg'])
    close
end

%% Trend in depth of sigma levels
RG_sigma.pres_trend = nan(xdim,ydim,zdim_sigma);
for n = 1:zdim_sigma
    sig_idx = find(abs(RG_sigma.sigma0 - RG_sigma.sigma0(n)) == ...
        min(abs(RG_sigma.sigma0 - RG_sigma.sigma0(n))));
    % calculate trends
    for a = 1:xdim
        for b = 1:ydim
            if sum(~isnan(RG_sigma.pressure(a,b,sig_idx,:))) > 102
                [~,~,x,~,~,~,~] = ...
                    leastsq2(RG.month,squeeze(RG_sigma.pressure(a,b,sig_idx,:)),0,2,[6 12]);
                RG_sigma.pres_trend(a,b,n) = x(2)*12;
            else
                RG_sigma.pres_trend(a,b,n) = NaN;
            end
        end
    end
    % plot trends
    figure; worldmap(latlim,[20 380]);
    title(['Trend in depth of \sigma = ' num2str(round(RG_sigma.sigma0(n),1))],'fontsize',16);
    setm(gca,'ffacecolor','w');
    setm(gca,'fontsize',12);
    pcolorm(double(repmat(RG.latitude',size(RG.longitude,1),1)),...
        double(repmat(RG.longitude,1,size(RG.latitude,1))),...
        RG_sigma.pres_trend(:,:,n));
    land = shaperead('landareas', 'UseGeoCoords',true);
    geoshow(land,'FaceColor',rgb('grey'));
    c=colorbar;
    caxis([-5 5]);
    colormap(cmocean('balance',16,'pivot',0));
    c.FontSize = 14;
    c.Label.String = 'dbars yr^{-1}';
    mlabel off; plabel off;
    exportgraphics(gcf,['/Users/sharp/Desktop/Depth_sigma/sig' num2str(10*RG_sigma.sigma0(n)) '.jpg'])
    close
end

%% plot trend in O2 on sigma levels


% import and average surface currents
% currents

% save models and variables
% vars = who;
% save(['CC-SSBE_' strrep(strrep(datestr(clock),':','-'),' ','-') '.mat'],'-v7.3');

% load data
%load('CC-SSBE_26-Aug-2021-11-59-54.mat');

%% correct for diapycnal diffusivity and horizontal advection
correctPhysics

% determine deepest of MLD/Ez
%RG.deep_layer = max(cat(4,squeeze(RG.zeu),squeeze(RG.mld)),[],4,'omitnan');
%RG.deep_layer = RG.zeu;
%RG.deep_layer = nan(size(RG.mld)); RG.deep_layer(:,:,:) = 100;

%% define spatial box and plot its extent along with surface O2 and currents
minlat = 16;
maxlat = 17;
minlon = 88;
maxlon = 89;
spatial_box

% plot O2 section for box of interest
plotSection([minlat maxlat],[minlon maxlon],'oxy_ENS',1000,RG,area_weights,0,150)
plotSection_climatology([minlat maxlat],[minlon maxlon],'oxy_ENS',300,RG,area_weights,100,250)
exportgraphics(gcf,['/Users/sharp/Desktop/Box_' num2str(minlat) ...
    '_' num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) '.jpg']);

% calculate respiration on presseure levels beneath Ez and maximum MLD
Resp_pres_levels

% calculate respiration on density levels beneath Ez and maximum MLD
Resp_dens_levels

% calculate and plot mixed layer and euphotic
% depth for box of interest
plotDepths

% calculate and plot seasonal respiration beneath euphotic
% depth for box of interest
calculateResp

% 
% figure; pcolor(RG.longitude(:,:,1,1),RG.latitude(:,:,1,1),mean(RG.ZLEE(:,:,1,:),4)); colorbar
% figure; pcolor(RG.longitude(:,:,1,1),RG.latitude(:,:,1,1),mean(RG.mld(:,:,1,:),4)); colorbar
figure; plot(RG.time,squeeze(RG.zeu(9,31,:)))
hold on; plot(RG.time,squeeze(RG.mld(9,31,:)))

% fit trend/harmonics to each lat/lon/depth grid cell
fitTrends
fitTrends_sigma

% plot hypoxic volume
plotHypox

% plot O2 inventory
plotO2

%% plot mean O2 profiles by latitude (uncorrected)
figure; set(gcf,'units','normalized','position',[0 0 0.4 1]); hold on;
set(gca,'XAxisLocation','top','YDir','reverse','fontsize',20);
ylabel('Pressure (dbar)');
xlabel('Oxygen (\mumol kg^{-1})');
% clrs = [0,0.4470,0.7410;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250]; n=1;
for l=1:15:145
    area_weights_temp = area_weights;
    area_weights_temp(:,RG.latitude < RG.latitude(l) | ...
        RG.latitude > RG.latitude(l+14),:,:) = NaN;
    mean_prof = squeeze(wmean(RG.oxy_ENS,area_weights_temp,'dim',[1 2 4],'omitnan'));
    %std_prof = squeeze(std(RG.oxy_ENS,area_weights_temp,[1 2 4]));
    plot(mean_prof,squeeze(RG.pressure),'linewidth',4);
    n=n+1;
end
legend({'15 to 30 N','30 to 45 N','45 to 60 N'},'location','southeast','fontsize',28)
exportgraphics(gcf,'/Users/sharp/Desktop/mean_oxy_profiles_ENS.jpg');

%% plot variables on depth levels
idxlat = find(RG.latitude > 10 & RG.latitude < 20);
idxlon = find(RG.longitude > -180 & RG.longitude < 180);
idxdepth = find(RG.pressure(1,1,:,1)==100);
figure; xlim([0 13]); hold on;
for y = 1:tdim/12
    plot(RG.month((y-1)*12+1:(y-1)*12+12)-...
        RG.month((y-1)*12+1)+1,...
         squeeze(RG.oxy_ENS(idxlon,idxlat,idxdepth,(y-1)*12+1:(y-1)*12+12)),...
         'linewidth',2);
end
legend(arrayfun(@num2str, 2004:2020, 'UniformOutput', 0))

%% plot variables on sigma levels
idxlat = find(RG_sigma.latitude(1,:,1,1)==40.5);
idxlon = find(RG_sigma.longitude(:,1,1,1)==225.5);
idxsig = find(RG_sigma.sigma0(1,1,:,1)==26.2);
figure; xlim([0 13]); hold on;
for y = 1:tdim/12
    plot(squeeze(RG_sigma.month(idxlon,idxlat,idxsig,(y-1)*12+1:(y-1)*12+12))-...
        RG_sigma.month(idxlon,idxlat,idxsig,(y-1)*12+1)+1,...
         squeeze(RG_sigma.oxy_ENS(idxlon,idxlat,idxsig,(y-1)*12+1:(y-1)*12+12)),...
         'linewidth',2);
end
legend(arrayfun(@num2str, 2004:2020, 'UniformOutput', 0))

%% plot oxygen climatology by latitude (uncorrected)
clrs = [0,0.4470,0.7410;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250]; n=1;
figure; set(gcf,'units','normalized','position',[0 0 1 0.5]);
set(gca,'fontsize',20); hold on;
plot_climatology_depth([-30 -40],[-180 180],'oxy_ENS',100,500,RG,'k',area_weights)
plot_climatology_depth([30 45],[-150 -105],'oxy_ENS',0,50,RG,clrs(2,:),area_weights)
plot_climatology_depth([45 60],[-150 -105],'oxy_ENS',0,50,RG,clrs(3,:),area_weights)
lgd = legend({'15 to 30 N','30 to 45 N','45 to 60 N'},'location',...
    'eastoutside','fontsize',28);
title(lgd,'0 to 50 dbars');
exportgraphics(gcf,'/Users/sharp/Desktop/integrated_oxy_climatologies_0_50_ENS.jpg');

%% plot how spice correction is done (39.5 N, 140.5 W, 40 dbars)
lon_idx = 1; lat_idx = 1; prs_idx = 25;
figure; set(gcf,'units','normalized','position',[0,1,0.5,0.5]);
set(gca,'fontsize',18); hold on;
scatter(squeeze(RG.spice_anom(lon_idx,lat_idx,prs_idx,:)),squeeze(RG.oxy_ENS_per_sat_anom(lon_idx,lat_idx,prs_idx,:)),...
    'or','filled');
plot([min(RG.spice_anom(lon_idx,lat_idx,prs_idx,:)) max(RG.spice_anom(lon_idx,lat_idx,prs_idx,:))],...
     [min(RG.spice_anom(lon_idx,lat_idx,prs_idx,:)).*min(RG.oxy_ENS_per_sat_slope(lon_idx,lat_idx,prs_idx,:))+...
          min(RG.oxy_ENS_per_sat_intercept(10,25,11,:)) ...
      max(RG.spice_anom(lon_idx,lat_idx,prs_idx,:)).*min(RG.oxy_ENS_per_sat_slope(lon_idx,lat_idx,prs_idx,:))+...
          min(RG.oxy_ENS_per_sat_intercept(10,25,11,:))],...
      '-k','linewidth',2);
text(min(RG.spice_anom(lon_idx,lat_idx,prs_idx,:))+0.05,min(RG.oxy_ENS_per_sat_anom(lon_idx,lat_idx,prs_idx,:)+2),...
    ['p-value = ' num2str(round(RG.oxy_ENS_per_sat_pval(lon_idx,lat_idx,prs_idx,1),2))]);
xlabel('Spice Anomaly (kg m^{-3})');
ylabel('Percent Oxygen Saturation Anomaly (\mumol kg^{-1})');
exportgraphics(gcf,'/Users/sharp/Desktop/spice_correction_example.jpg');

%% plot fits in example grid cell (39.5 N, 140.5 W, 40 dbars)
lon_idx = find(RG.longitude(:,1,1,1) == 360-49.5);
lat_idx = find(RG.latitude(1,:,1,1) == 39.5);
prs_idx = find(RG.pressure(1,1,:,1) == 40);
figure; set(gcf,'units','normalized','position',[0,1,1,0.5]);
set(gca,'fontsize',20); hold on;
scatter(squeeze(RG.time(lon_idx,lat_idx,prs_idx,:)),squeeze(RG.oxy_ENS(lon_idx,lat_idx,prs_idx,:)),100,'o',...
    'filled','MarkerEdgeColor',[0.8500,0.3250,0.0980],'MarkerFaceColor',[0.8500,0.3250,0.0980]);
plot(squeeze(RG.time(lon_idx,lat_idx,prs_idx,:)),cell2mat(RG.oxy_ENS_fit(lon_idx,lat_idx,prs_idx)),'k','linewidth',4);
datetick('x'); ylabel('Oxygen (\mumol kg^{-1})');
exportgraphics(gcf,'/Users/sharp/Desktop/oxy_fit_example.jpg');

figure; set(gcf,'units','normalized','position',[0,1,1,0.5]);
set(gca,'fontsize',20); hold on;
plot([RG.time(lon_idx,lat_idx,prs_idx,1),RG.time(lon_idx,lat_idx,prs_idx,end)],[0,0],'--k','linewidth',1);
plot(squeeze(RG.time(lon_idx,lat_idx,prs_idx,:)),cell2mat(RG.oxy_ENS_resid(lon_idx,lat_idx,prs_idx)),'k','linewidth',4);
plot(squeeze(RG.time(lon_idx,lat_idx,prs_idx,:)),smoothdata(cell2mat(RG.oxy_ENS_resid(lon_idx,lat_idx,prs_idx)),...
    'movmean',12),'color',[0.8500,0.3250,0.0980],'linewidth',2);
datetick('x'); ylabel('Oxygen Anomaly (\mumol kg^{-1})');
exportgraphics(gcf,'/Users/sharp/Desktop/oxy_anom_example.jpg');

%% plot integrated temperature in entire domain over time (depth)
plot_integrated_anom_depth(15:60,-155:-105,'temp',1000,2000,RG,area_weights)
exportgraphics(gcf,'/Users/sharp/Desktop/temp_all_1000_2000.jpg');
plot_integrated_anom_depth(15:60,-155:-105,'temp_resid_notrend',1000,2000,RG,area_weights)
exportgraphics(gcf,'/Users/sharp/Desktop/temp_anom_all_1000_2000.jpg');

%% plot integrated O2 in entire domain over time (depth)
plot_integrated_anom_depth(39:41,-136:-134,'temp',280,320,RG,area_weights)
exportgraphics(gcf,'/Users/sharp/Desktop/temp_40_225_300.jpg');
plot_integrated_anom_depth(15:60,-155:-105,'oxy_ENS_resid_notrend',0,2000,RG,area_weights)
exportgraphics(gcf,'/Users/sharp/Desktop/oxy_anom_all.jpg');

%% plot integrated temperature in entire domain over time (sigma)
plot_integrated_anom_sigma(15:60,-155:-105,'temp',26.5,27,RG,area_weights)
%exportgraphics(gcf,'/Users/sharp/Desktop/temp_all_1000_2000.jpg');
plot_integrated_anom_sigma(15:60,-155:-105,'temp_resid_notrend',26.5,27,RG,area_weights)
%exportgraphics(gcf,'/Users/sharp/Desktop/temp_anom_all_1000_2000.jpg');

%% plot integrated O2 in entire domain over time (sigma)
plot_integrated_anom_depth(15:60,-155:-105,'oxy_ENS',26.5,27,RG,area_weights)
text(datenum([2005 1 1]),132,'26.5 to 27.0')
exportgraphics(gcf,'/Users/sharp/Desktop/oxy_lon_150_lat_40_60_sigma_265_270.jpg');
plot_integrated_anom_sigma(15:35,-151:-149,'oxy_ENS_resid_notrend',26.5,27,RG,area_weights)
text(datenum([2005 1 1]),132,'26.5 to 27.0')
exportgraphics(gcf,'/Users/sharp/Desktop/oxy_lon_150_lat_15_35_sigma_265_270.jpg');
plot_integrated_anom_sigma(39:41,-155:-125,'oxy_ENS_resid_notrend',26.5,27,RG,area_weights)
exportgraphics(gcf,'/Users/sharp/Desktop/oxy_lat_40_lon_205_235_sigma_265_270.jpg');


%% plot trends in O2 as longitudinal segment

for l = 150.5:-10:110.5
    lon_idx = find(RG.longitude(:,1,1,1) == 360-l);
    trend_oxy_ENS = cell2mat(RG.oxy_ENS(lon_idx,:,:));
    trend_oxy_ENS = permute(trend_oxy_ENS(2,:,:),[2 3 1]);
    figure; hold on;
    contourf(squeeze(RG.latitude(1,:,1,1)),squeeze(RG.pressure(1,1,:,1)),...
        12.*trend_oxy_ENS');
    contour(squeeze(RG.latitude(1,:,1,1)),squeeze(RG.pressure(1,1,:,1)),...
        squeeze(mean(RG.sigma0(lon_idx,:,:,:),4,'omitnan'))','-k',...
        'ShowText','on','linewidth',1);
    set(gca,'YDir','reverse','fontsize',14);
    ylim([0 1000]);
    caxis([-1 1]);
    text(RG.latitude(1,2,1,1),900,['Lon = ' num2str(floor(l))],...
        'fontsize',20,'fontweight','bold');
    c=colorbar; c.Label.String = 'Corrected O_{2} Saturation Trend (% yr^{-1})';
    colormap(cmocean('balance','pivot',0));
    xlabel('Latitude');
    ylabel('Depth (dbar)');
    exportgraphics(gcf,['/Users/sharp/Desktop/oxy_lon_' num2str(floor(l)) '.jpg']);
end

%% plot trends in O2 as latitudinal segment

for l = 20.5:5:55.5
    lat_idx = find(RG.latitude(1,:,1,1) == l);
    trend_oxy_ENS = cell2mat(permute(RG.oxy_ENS_per_sat_corr_x(:,lat_idx,:),[2 1 3]));
    trend_oxy_ENS = permute(trend_oxy_ENS(2,:,:),[2 3 1]);
    figure; hold on;
    contourf(squeeze(RG.longitude(:,1,1,1)),squeeze(RG.pressure(1,1,:,1)),...
        12.*trend_oxy_ENS');
    contour(squeeze(RG.longitude(:,1,1,1)),squeeze(RG.pressure(1,1,:,1)),...
        squeeze(mean(RG.sigma0(:,lat_idx,:,:),4,'omitnan'))','-k',...
        'ShowText','on','linewidth',1);
    set(gca,'YDir','reverse','fontsize',14);
    ylim([0 1000]);
    caxis([-1 1]);
    text(RG.longitude(2,1,1,1),900,['Lat = ' num2str(floor(l))],...
        'fontsize',20,'fontweight','bold');    
    c=colorbar; c.Label.String = 'Corrected O_{2} Saturation Trend (% yr^{-1})';
    colormap(cmocean('balance','pivot',0));
    xlabel('Longitude');
    ylabel('Depth (dbar)');
    exportgraphics(gcf,['/Users/sharp/Desktop/oxy_sat_corr_trend_lat_' num2str(floor(l)) '.jpg']);
end

%% calculate respiration rate


%% plot anomaly sections

% plot anomalies fit along pressure levels with respect to pressure
plot_anom_pres([40 45],[-122 -117],'temp',0,400,GOBAI);
exportgraphics(gcf,['Figures/Anomaly Plots/temp_PresFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_pres.jpg']);
plot_anom_pres([43 44],[-146 -145],'sal',0,400,RG);
exportgraphics(gcf,['Figures/Anomaly Plots/sal_PresFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_pres.jpg']);
plot_anom_pres([43 44],[-146 -145],'oxy_ENS',0,400,RG);
exportgraphics(gcf,['Figures/Anomaly Plots/oxy_PresFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_pres.jpg']);
plot_anom_pres([minlat maxlat],[minlon maxlon],'oxy_ENS_per_sat',0,400,RG);
exportgraphics(gcf,['Figures/Anomaly Plots/oxy_sat_PresFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_pres.jpg']);

% plot anomalies fit along sigma levels with respect to pressure
sig_plot_anom_pres([minlat maxlat],[minlon maxlon],'temp',0,400,RG_sigma);
exportgraphics(gcf,['Figures/Anomaly Plots/temp_SigFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_pres.jpg']);
sig_plot_anom_pres([minlat maxlat],[minlon maxlon],'sal',0,400,RG_sigma);
exportgraphics(gcf,['Figures/Anomaly Plots/sal_SigFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_pres.jpg']);
sig_plot_anom_pres([minlat maxlat],[minlon maxlon],'oxy_ENS',0,400,RG_sigma);
exportgraphics(gcf,['Figures/Anomaly Plots/oxy_SigFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_pres.jpg']);
sig_plot_anom_pres([minlat maxlat],[minlon maxlon],'oxy_ENS_per_sat',0,400,RG_sigma);
exportgraphics(gcf,['Figures/Anomaly Plots/oxy_sat_SigFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_pres.jpg']);

% plot anomalies fit along sigma levels with respect to pressure
plot_anom_sigma([minlat maxlat],[minlon maxlon],'temp',27,RG);
exportgraphics(gcf,['Figures/Anomaly Plots/temp_PresFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_sigma.jpg']);
plot_anom_sigma([minlat maxlat],[minlon maxlon],'sal',27,RG);
exportgraphics(gcf,['Figures/Anomaly Plots/sal_PresFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_sigma.jpg']);
plot_anom_sigma([minlat maxlat],[minlon maxlon],'oxy_ENS',27,RG);
exportgraphics(gcf,['Figures/Anomaly Plots/oxy_PresFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_sigma.jpg']);
plot_anom_sigma([minlat maxlat],[minlon maxlon],'oxy_ENS_per_sat',27,RG);
exportgraphics(gcf,['Figures/Anomaly Plots/oxy_sat_PresFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_sigma.jpg']);

% plot anomalies fit along sigma levels with respect to sigma
sig_plot_anom_sigma([minlat maxlat],[minlon maxlon],'temp',27,RG_sigma);
exportgraphics(gcf,['Figures/Anomaly Plots/temp_SigFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_sigma.jpg']);
sig_plot_anom_sigma([minlat maxlat],[minlon maxlon],'sal',27,RG_sigma);
exportgraphics(gcf,['Figures/Anomaly Plots/sal_SigFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_sigma.jpg']);
sig_plot_anom_sigma([minlat maxlat],[minlon maxlon],'oxy_ENS',27,RG_sigma);
exportgraphics(gcf,['Figures/Anomaly Plots/oxy_SigFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_sigma.jpg']);
sig_plot_anom_sigma([minlat maxlat],[minlon maxlon],'oxy_ENS_per_sat',27,RG_sigma);
exportgraphics(gcf,['Figures/Anomaly Plots/oxy_sat_SigFit_' num2str(minlat) '_' ...
    num2str(maxlat) '_' num2str(360+minlon) '_' num2str(360+maxlon) ...
    '_sigma.jpg']);

%% plot integrated anomalies

plot_integrated_anom_depth(15,60,'oxy_ENS_resid_notrend',0,2000,RG,area_weights)
plot_integrated_anom_depth(45:50,-150:-145,'spice',0,10,RG,area_weights)

plot_integrated_anom_sigma(43.5,145.5,'temp',26,26.5,RG,area_weights)
plot_integrated_anom_sigma(45:50,-150:-145,'oxy_ENS_corr',26,26.5,RG,area_weights)

%% plot correlation between temperature anomaly and corrected O2 anomaly
RG.corr_temp_oxy = nan(size(RG.temp_mean));
RG.corr_temp_oxy_pval = nan(size(RG.temp_mean));
for a = 1:size(RG.longitude,1)
    for b = 1:size(RG.latitude,2)
        for c = 1:size(RG.pressure,3)
            [RG.corr_temp_oxy(a,b,c),RG.corr_temp_oxy_pval(a,b,c)] = ...
            corr(squeeze(RG.temp(a,b,1,:)),...
                squeeze(RG.oxy_ENS_per_sat(a,b,c,:)));
        end
    end
end
%% plot correlation spatially
depthidx = 1:20;
figure; worldmap(latlim,lonlim);
title('Correlation','fontsize',16)
%set(gcf,'Position',[617, 599, 820, 820])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
lat = RG.latitude(:,:,1,1); lon = RG.longitude(:,:,1,1);
pcolorm(lat,lon,mean(RG.corr_temp_oxy(:,:,depthidx),3));
scatterm(lat(mean(RG.corr_temp_oxy_pval(:,:,depthidx),3)>0.05)+0.5,...
    lon(mean(RG.corr_temp_oxy_pval(:,:,depthidx),3)>0.05)+0.5,'.k');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
colormap(cmocean('balance'));
caxis([-1 1]);
c.Label.String = 'Correlation';
c.FontSize = 12;
exportgraphics(gcf,['/Users/sharp/Desktop/corr_temp_O2.jpg']);

%% plot correlation between salinity anomaly and O2 anomaly
RG.corr_sal_oxy = nan(size(RG.temp_mean));
RG.corr_sal_oxy_pval = nan(size(RG.temp_mean));
for a = 1:size(RG.longitude,1)
    for b = 1:size(RG.latitude,2)
        for c = 1:size(RG.pressure,3)
            [RG.corr_sal_oxy(a,b,c),RG.corr_sal_oxy_pval(a,b,c)] = ...
            corr(squeeze(RG.sal(a,b,1,:)),...
                squeeze(RG.oxy_ENS(a,b,c,:)));
        end
    end
end
%% plot correlation spatially
depthidx = 20;
figure; worldmap(latlim,lonlim);
title('Correlation','fontsize',16)
%set(gcf,'Position',[617, 599, 820, 820])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
lat = RG.latitude(:,:,1,1); lon = RG.longitude(:,:,1,1);
pcolorm(lat,lon,mean(RG.corr_sal_oxy(:,:,depthidx),3));
scatterm(lat(mean(RG.corr_sal_oxy_pval(:,:,depthidx),3)>0.05)+0.5,...
    lon(mean(RG.corr_sal_oxy_pval(:,:,depthidx),3)>0.05)+0.5,'.k');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar;
colormap(cmocean('balance'));
caxis([-1 1]);
c.Label.String = 'Correlation';
c.FontSize = 12;
exportgraphics(gcf,['/Users/sharp/Desktop/corr_sal_O2_200m.jpg']);

%% plot all anomalies in depth layers

% % 0 to 150 dbar
% temp_anom_0_150 = cell2mat(RG.temp_resid(:,:,1:10));
% temp_anom_0_150 = temp_anom_0_150(~isnan(temp_anom_0_150));
% oxy_anom_0_150 = cell2mat(RG.oxy_resid(:,:,1:10));
% oxy_anom_0_150 = oxy_anom_0_150(~isnan(oxy_anom_0_150));
% figure; scatter(temp_anom_0_150,oxy_anom_0_150);
% [r_0_150,p_0_150] = corr(temp_anom_0_150,oxy_anom_0_150)
% % 160 to 500 dbar
% temp_anom_160_500 = cell2mat(RG.temp_resid(:,:,17:34));
% temp_anom_160_500 = temp_anom_160_500(~isnan(temp_anom_160_500));
% oxy_anom_160_500 = cell2mat(RG.oxy_resid(:,:,17:34));
% oxy_anom_160_500 = oxy_anom_160_500(~isnan(oxy_anom_160_500));
% figure; scatter(temp_anom_160_500,oxy_anom_160_500);
% [r_160_500,p_160_500] = corr(temp_anom_160_500,oxy_anom_160_500)
% % 550 to 1000 dbar
% temp_anom_550_1000 = cell2mat(RG.temp_resid(:,:,35:44));
% temp_anom_550_1000 = temp_anom_550_1000(~isnan(temp_anom_550_1000));
% oxy_anom_550_1000 = cell2mat(RG.oxy_resid(:,:,35:44));
% oxy_anom_550_1000 = oxy_anom_550_1000(~isnan(oxy_anom_550_1000));
% figure; scatter(temp_anom_550_1000,oxy_anom_550_1000);
% [r_550_1000,p_550_1000] = corr(temp_anom_550_1000,oxy_anom_550_1000)
% % 1050 to 1975 dbar
% temp_anom_1050_1975 = cell2mat(RG.temp_resid(:,:,45:58));
% temp_anom_1050_1975 = temp_anom_1050_1975(~isnan(temp_anom_1050_1975));
% oxy_anom_1050_1975 = cell2mat(RG.oxy_resid(:,:,45:58));
% oxy_anom_1050_1975 = oxy_anom_1050_1975(~isnan(oxy_anom_1050_1975));
% figure; scatter(temp_anom_1050_1975,oxy_anom_1050_1975);
% [r_1050_1975,p_1050_1975] = corr(temp_anom_1050_1975,oxy_anom_1050_1975)
% 

% % Compare floats to CC-SSBE, WOA, and ESPER
% compareFloat

%% interpolate oxygen for each "profile" and calculate column inventory and hypoxic thickness
% RG.hypoxic_thickness = nan(xdim,ydim,tdim);
% RG.O2_col = nan(xdim,ydim,tdim);
% for a = 1:xdim
%     for b = 1:ydim
%         for t = 1:tdim
%             if sum(~isnan(RG.oxy_ENS(a,b,:,t))) ~= 0
%                 % interpolated profile (oxygen)
%                 oxy_temp = double(interp1(RG.pressure,...
%                     squeeze(RG.oxy_ENS(a,b,:,t)),[1:2000]','linear','extrap'));
%                 % interpolated profile (mass)
%                 kg_temp = double(interp1(RG.pressure,...
%                     squeeze(RG.kg(a,b,:,t)),[1:2000]','linear','extrap'));
%                 % inventory
%                 RG.O2_col(a,b,t) = sum(oxy_temp.*... % umol/kg
%                         kg_temp./... % umol
%                         1e6,'omitnan')./ ... % mol
%                         (area_weights(a,b).*1e6); % mol/m2
%                 % hypoxia
%                 is_hypoxic = oxy_temp < 10;
%                 RG.hypoxic_thickness(a,b,t) = sum(is_hypoxic);
%             else
%             end
%         end
%     end
% end
% RG.hypoxic_thickness(isnan(RG.oxy_ENS(:,:,:,1))) = NaN;
% clear oxy_temp kg_temp is_hypoxic
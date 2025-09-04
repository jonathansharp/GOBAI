% adjust_o2_float_data
%
% DESCRIPTION:
% This function bins float and glodap data, co-locates corresponding bins,
% compares the two datasets, and calculates a bulk correction factor based
% on the discrepancies between the two.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/22/2023

function adjust_o2_float_data(float_file_ext,glodap_year,snap_date)

%% load interpolated float and glodap data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['O2/Data/processed_float_o2_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date');
load(['O2/Data/processed_glodap_o2_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');
load(['O2/Data/processed_wod_ctd_o2_data_' num2str(glodap_year) '.mat'],...
    'wod_data');

% idx_f = (float_data.LON > 0 & float_data.LON < 35) & ...
%     (float_data.LAT > 18 & float_data.LAT < 43) & ...
%     float_data.YEAR < 2009;
% idx_g = (glodap_data.LON > 0 & glodap_data.LON < 35) & ...
%     (glodap_data.LAT > 18 & glodap_data.LAT < 43) & ...
%     glodap_data.YEAR < 2009;
% figure; scatter(float_data.OXY(idx_f),float_data.PRES(idx_f));
% figure; scatter(glodap_data.OXY(idx_g),glodap_data.PRES(idx_g));
% 
% prof_ids = unique(float_data.PROF_ID(idx_f));
% figure; hold on;
% set(gca,'ydir','reverse');
% for p = 1:length(prof_ids)
%     idx = prof_ids(p) == float_data.PROF_ID;
%     plot(float_data.OXY(idx_f),float_data.PRES(idx_f)); 
% end

%% import WOA climatologies
temp_path = [pwd '/Data/WOA/TEMPERATURE/'];
sal_path = [pwd '/Data/WOA/SALINITY/'];
oxy_path = [pwd '/Data/WOA/OXYGEN/'];
for m = 1:12
    WOA.TMP(:,:,:,m) = ncread([temp_path 'woa18_decav_t' sprintf('%02d',m) '_01.nc'],'t_an');
    WOA.SAL(:,:,:,m) = ncread([sal_path 'woa18_decav_s' sprintf('%02d',m) '_01.nc'],'s_an');
    WOA.OXY(:,:,:,m) = ncread([oxy_path 'woa18_all_o',sprintf('%02d',m),'_01.nc'],'o_an');
end
% Import and format WOA latitude and longitude
WOA.LAT   = ncread([oxy_path 'woa18_all_o01_01.nc'],'lat');
WOA.LON   = ncread([oxy_path 'woa18_all_o01_01.nc'],'lon');
WOA.DEPTH = ncread([oxy_path 'woa18_all_o01_01.nc'],'depth');
[~,lat_3d,depth_3d] = ndgrid(WOA.LON,WOA.LAT,WOA.DEPTH);
WOA.PRES = -gsw_p_from_z(depth_3d,lat_3d);
% clean up
clear m temp_path sal_path oxy_path

%% compare float data to WOA
% index to valid data points
idx = find(~isnan(float_data.OXY) & ~isnan(float_data.LAT) & ...
    ~isnan(float_data.LON) & ~isnan(float_data.PRES) & ...
    ~isnan(float_data.TIME));
% obtain month from date
date = datevec(datenum(float_data.YEAR,0,float_data.DAY));
float_data.MONTH = date(:,2);
clear date
% pre-allocate matches
WOA_match = nan(size(float_data.OXY));
% match float data with WOA
for i=1:length(idx)
    idx_lon = find(min(abs(WOA.LON - float_data.LON(idx(i)))) == ...
        abs(WOA.LON - float_data.LON(idx(i))));
    idx_lat = find(min(abs(WOA.LAT - float_data.LAT(idx(i)))) == ...
        abs(WOA.LAT - float_data.LAT(idx(i))));
    idx_pres = find(min(abs(squeeze(WOA.PRES(idx_lon,idx_lat,:)) - float_data.PRES(idx(i)))) == ...
        abs(squeeze(WOA.PRES(idx_lon,idx_lat,:)) - float_data.PRES(idx(i))));
    idx_mnth = find(min(abs((1:12)' - float_data.MONTH(idx(i)))) == ...
        abs((1:12)' - float_data.MONTH(idx(i))));
    WOA_match(idx(i)) = WOA.OXY(idx_lon(1),idx_lat(1),idx_pres(1),idx_mnth(1));
end
% calculate differences
WOA_delta = float_data.OXY-WOA_match;
WOA_delta_per = WOA_delta./float_data.OXY;
WOA_delta_per(isinf(WOA_delta_per))=NaN;
st_dev_delta = std(WOA_delta,[],'omitnan');
mean_delta = mean(WOA_delta,'omitnan');
% clean up
clear idx idx_lon idx_lat idx_pres idx_mnth WOA i

%% plot differences between float data and WOA
% histogram of differences
figure; hold on
histogram(WOA_delta);
set(gca,'fontsize',16);
ax = gca;
plot([mean_delta+3*st_dev_delta mean_delta+3*st_dev_delta],ax.YLim,'r','linewidth',2);
plot([mean_delta-3*st_dev_delta mean_delta-3*st_dev_delta],ax.YLim,'r','linewidth',2);
xlabel('Float [O_{2}] - WOA [O_{2}]');
exportgraphics(gcf,[pwd '/O2/Figures/Data/WOA_comp_histogram_' file_date float_file_ext '.png']);
close
% scatter of delta values
figure; hold on
set(gca,'fontsize',16);
xlabel('WOA [O_{2}]');
ylabel('Float [O_{2}] - WOA [O_{2}]');
[counts,bin_centers] = ...
    hist3([WOA_match,WOA_delta],'Edges',{0:10:460 -345:15:345});
h=pcolor(bin_centers{1}-mean(diff(bin_centers{1}))/2,...
    bin_centers{2}-mean(diff(bin_centers{2}))/2,counts');
plot([0 450],[0 0],'k--');
set(h,'EdgeColor','none');
xlim([-0.5 460]); ylim([-345.5 345]);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log')
clim([1e0 1e5]);
c=colorbar;
c.Label.String = 'log10(Bin Counts)';
exportgraphics(gcf,[pwd '/O2/Figures/Data/WOA_comp_scatter_' file_date float_file_ext '.png']);
close
% clean up
clear counts bin_centers c h myColorMap

%% remove data points based on global range test
% 479 umol/kg is max for WOD range checks
idx_rem = float_data.OXY < 0 | float_data.OXY > 480;
disp([num2str(sum(idx_rem)) ' data points removed by global range test (' ...
    num2str(100*(sum(idx_rem)/length(float_data.PROF_ID))) ' % of data)']);
vars = fieldnames(float_data);
for v = 1:length(vars)
    float_data.(vars{v})(idx_rem) = [];
end
WOA_delta(idx_rem) = [];
WOA_delta_per(idx_rem) = [];
WOA_match(idx_rem) = [];

%% remove data points more than 3 sigmas from WOA value for each depth
pres_levels = unique(float_data.PRES);
idx_rem = false(length(WOA_delta),1);
st_dev_delta_pres = nan(size(pres_levels));
mean_delta_pres = nan(size(pres_levels));
for z = 1:length(pres_levels)
    pres_idx = float_data.PRES == pres_levels(z);
    st_dev_delta_pres(z) = std(WOA_delta(pres_idx),[],'omitnan');
    mean_delta_pres(z) = mean(WOA_delta(pres_idx),'omitnan');
    idx_rem(pres_idx) = WOA_delta(pres_idx) > mean_delta_pres(z)+3.*st_dev_delta_pres(z) | ...
        WOA_delta(pres_idx) < mean_delta_pres(z)-3.*st_dev_delta_pres(z);
end
disp([num2str(sum(idx_rem)) ' data points removed by WOA comparison test (' ...
    num2str(100*(sum(idx_rem)/length(float_data.PROF_ID))) ' % of data)']);
vars = fieldnames(float_data);
% plot(mean_delta_pres,pres_levels); xlabel('Float - WOA'); ylabel('pres.');
for v = 1:length(vars)
    float_data.(vars{v})(idx_rem) = [];
end
WOA_delta(idx_rem) = [];
WOA_delta_per(idx_rem) = [];
WOA_match(idx_rem) = [];

%% histogram of differences (corrected)
figure; hold on
histogram(WOA_delta);
set(gca,'fontsize',16);
ax = gca;
plot([mean_delta+3*st_dev_delta mean_delta+3*st_dev_delta],ax.YLim,'r','linewidth',2);
plot([mean_delta-3*st_dev_delta mean_delta-3*st_dev_delta],ax.YLim,'r','linewidth',2);
xlabel('Float [O_{2}] - WOA [O_{2}]');
exportgraphics(gcf,[pwd '/O2/Figures/Data/WOA_comp_histogram_' file_date float_file_ext '.png']);
close
% clean up
clear WOA_match WOA_delta WOA_delta_per mean_delta st_dev_delta idx_rem v vars

% prof_ids = unique(float_data.PROF_ID);
% for p = 1:length(prof_ids)
%     idx = prof_ids(p) == float_data.PROF_ID;
%     figure;
%     plot(float_data.OXY(idx),float_data.PRES(idx));
%     set(gca,'ydir','reverse');
%     keyboard
%     close
% end

%% Co-locate and compare float and glodap data (via profile crossovers)
pres_levels = [2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975]';
float_profile_IDs = unique(float_data.PROF_ID);
% pre-allocate
crossover.id = [];
crossover.lon = [];
crossover.lat = [];
crossover.pres = [];
crossover.sigma_float = [];
crossover.oxy_float = [];
crossover.oxy_float_gradient = [];
crossover.temp_float = [];
crossover.sal_float = [];
crossover.oxy_ship = [];
crossover.oxy_ship_gradient = [];
crossover.sigma_ship = [];
crossover.oxy_delta = [];
crossover.oxy_delta_per = [];
% cycle through float profiles 
ship_lon = [glodap_data.LON;wod_data.LON];
ship_lat = [glodap_data.LAT;wod_data.LAT];
ship_time = [glodap_data.TIME;wod_data.TIME];
ship_sigma = [glodap_data.SIGMA;wod_data.SIGMA];
ship_oxy = [glodap_data.OXY;wod_data.OXY];
for p = 1:length(float_profile_IDs)
    % float position
    lon = mean(float_data.LON(float_data.PROF_ID==float_profile_IDs(p)));
    lat = mean(float_data.LAT(float_data.PROF_ID==float_profile_IDs(p)));
    time = mean(float_data.TIME(float_data.PROF_ID==float_profile_IDs(p)));
    % match index
    idx_lon = abs(ship_lon - lon) < 0.5;
    idx_lat = abs(ship_lat - lat) < 0.5;
    idx_time = abs(ship_time - time) < 30;
    idx_all = idx_lon & idx_lat & idx_time;
    % if there is a matching ship profile
    if sum(idx_all) > 0
        % float oxygen and depth
        oxy_float = float_data.OXY(float_data.PROF_ID==float_profile_IDs(p));
        pres_float = float_data.PRES(float_data.PROF_ID==float_profile_IDs(p));
        sigma_float = float_data.SIGMA(float_data.PROF_ID==float_profile_IDs(p));
        lon_float = float_data.LON(float_data.PROF_ID==float_profile_IDs(p));
        lat_float = float_data.LAT(float_data.PROF_ID==float_profile_IDs(p));
        temp_float = float_data.TEMP(float_data.PROF_ID==float_profile_IDs(p));
        sal_float = float_data.SAL(float_data.PROF_ID==float_profile_IDs(p));
        % combine float profiles if there are more than one
        if length(oxy_float) > 58
            oxy_temp = nan(length(pres_levels),1);
            pres_temp = nan(length(pres_levels),1);
            sigma_temp = nan(length(pres_levels),1);
            lon_temp = nan(length(pres_levels),1);
            lat_temp = nan(length(pres_levels),1);
            temp_temp = nan(length(pres_levels),1);
            sal_temp = nan(length(pres_levels),1);
            for z = 1:length(pres_levels)
                idx = pres_float == pres_levels(z);
                oxy_temp(z) = mean(oxy_float(idx),'omitnan');
                pres_temp(z) = mean(pres_float(idx),'omitnan');
                sigma_temp(z) = mean(sigma_float(idx),'omitnan');
                lon_temp(z) = mean(lon_float(idx),'omitnan');
                lat_temp(z) = mean(lat_float(idx),'omitnan');
                temp_temp(z) = mean(temp_float(idx),'omitnan');
                sal_temp(z) = mean(sal_float(idx),'omitnan');
            end
            oxy_float = oxy_temp(~isnan(oxy_temp));
            pres_float = pres_temp(~isnan(pres_temp));
            sigma_float = sigma_temp(~isnan(sigma_temp));
            lon_float = lon_temp(~isnan(lon_temp));
            lat_float = lat_temp(~isnan(lat_temp));
            temp_float = temp_temp(~isnan(temp_temp));
            sal_float = sal_temp(~isnan(sal_temp));
        end
        % get oxygen gradient
        oxy_float_gradient = [0;diff(oxy_float)];
        % pressure index
        idx_pres = ismember(pres_levels,pres_float);
        % extract ship data
        oxy_ship_temp = ship_oxy(idx_all);
        sigma_ship_temp = ship_sigma(idx_all);
        % take average of all matching ship profiles (helps if there are multiple)
        oxy_ship = nan(58,1);
        sigma_ship = nan(58,1);
        for gld_z = 1:58
            oxy_ship(gld_z) = mean(oxy_ship_temp(gld_z:58:end),'omitnan');
            sigma_ship(gld_z) = mean(sigma_ship_temp(gld_z:58:end),'omitnan');
        end
        % get oxygen gradient
        oxy_ship_gradient = [0;diff(oxy_ship)];
        % log matched data
        crossover.id = [crossover.id;repmat(float_profile_IDs(p),length(pres_float),1)];
        crossover.lon = [crossover.lon;lon_float];
        crossover.lat = [crossover.lat;lat_float];
        crossover.pres = [crossover.pres;pres_float];
        crossover.oxy_float = [crossover.oxy_float;oxy_float];
        crossover.oxy_float_gradient = [crossover.oxy_float_gradient;oxy_float_gradient];
        crossover.temp_float = [crossover.temp_float;temp_float];
        crossover.sal_float = [crossover.sal_float;sal_float];
        crossover.oxy_ship = [crossover.oxy_ship;oxy_ship(idx_pres)];
        crossover.oxy_ship_gradient = [crossover.oxy_ship_gradient;oxy_ship_gradient(idx_pres)];
        crossover.sigma_float = [crossover.sigma_float;sigma_float];
        crossover.sigma_ship = [crossover.sigma_ship;sigma_ship];
        crossover.oxy_delta = [crossover.oxy_delta;oxy_float-oxy_ship(idx_pres)];
        crossover.oxy_delta_per = [crossover.oxy_delta_per;...
            (oxy_float-oxy_ship(idx_pres))./oxy_float];
    end
end

% calculate float oxygen saturation
Ts = log((298.15 - crossover.temp_float)./(273.15 + crossover.temp_float)); 
% The coefficents below are from the second column of Table 1 of Garcia and
% Gordon (1992)
a0 =  5.80871; 
a1 =  3.20291;
a2 =  4.17887;
a3 =  5.10006;
a4 = -9.86643e-2;
a5 =  3.80369;
b0 = -7.01577e-3;
b1 = -7.70028e-3;
b2 = -1.13864e-2;
b3 = -9.51519e-3;
c0 = -2.75915e-7;
crossover.oxy_sat_float = exp(a0 + Ts.*(a1 + Ts.*(a2 + Ts.*(a3 + Ts.*(a4 + a5*Ts)))) ...
          + crossover.sal_float.*(b0 + Ts.*(b1 + Ts.*(b2 + b3*Ts)) + c0*crossover.sal_float));
crossover.oxy_float_sat_per = 100*(crossover.oxy_float./crossover.oxy_sat_float);
crossover.oxy_delta_sat_per = 100*(crossover.oxy_delta./crossover.oxy_sat_float);

% index to below 300 dbars
idx = crossover.pres > 300 & ~isnan(crossover.oxy_float) & ~isnan(crossover.oxy_ship);
% clean up
clear pres_levels float_profile_IDs lon lat time oxy_float pres_float
clear idx_lon idx_lat idx_time oxy_temp pres_temp idx_pres oxy_ship
% fit delta against [O2]
mdl = fitlm(crossover.oxy_sat_float(idx),crossover.oxy_delta_sat_per(idx),'Intercept',true);
slp = mdl.Coefficients.Estimate(2);
int = mdl.Coefficients.Estimate(1);
% apply linear correction to float O2
corr_fac = slp .* crossover.oxy_sat_float + int; % modelled saturation delta
crossover.oxy_float_sat_per_corr = crossover.oxy_float_sat_per - corr_fac;
crossover.oxy_float_corr = (crossover.oxy_float_sat_per_corr./100).*crossover.oxy_sat_float;
% re-calculate delta
crossover.oxy_delta_corr = crossover.oxy_float_corr-crossover.oxy_ship;
crossover.oxy_delta_per_corr = (crossover.oxy_float_corr-crossover.oxy_ship)./crossover.oxy_float_corr;
% save crossover data
if ~exist([pwd '/O2/Data'],'dir'); mkdir('Data'); end
save(['O2/Data/crossover_data_' file_date float_file_ext],'crossover','idx','-v7.3')
% save correction factors
if ~exist([pwd '/O2/Data'],'dir'); mkdir('O2/Data'); end
save(['O2/Data/float_corr_' file_date float_file_ext],'slp','int');
clear slp int

%% plot uncorrected float vs. ship residuals
load(['O2/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
scatter(crossover.oxy_ship(idx),crossover.oxy_float(idx),20,'.');
plot([0,450],[0 450],'k--');
ylabel('Binned BGC Argo Oxygen Data (\mumol kg^{-1})');
xlabel('Binned GLODAP/WOD Oxygen Data (\mumol kg^{-1})');
crossover.err_md = median(crossover.oxy_delta(idx));
crossover.std = std(crossover.oxy_delta(idx));
crossover.rmse = sqrt(mean((crossover.oxy_delta(idx)).^2));
crossover.r2 = corr(crossover.oxy_ship(idx),crossover.oxy_float(idx));
crossover.slope = polyfit(crossover.oxy_ship(idx),crossover.oxy_float(idx),1);
crossover.slope = crossover.slope(1);
text(300,120,['Med. Err. = ' num2str(round(crossover.err_md,2))],'fontsize',12);
text(300,90,['RMSE = ' num2str(round(crossover.rmse,1))],'fontsize',12);
text(300,60,['R^{2} = ' num2str(round(crossover.r2,2))],'fontsize',12);
exportgraphics(gcf,['O2/Figures/Data/delta_vs_float_300_uncorr_' file_date float_file_ext '.png']);
close
% histogram
figure; hold on
histogram(crossover.oxy_delta(idx),'Normalization','probability');
ylim('manual');
xlim([-30 30]);
set(gca,'fontsize',16);
sig2_min = double(mean(crossover.oxy_delta(idx))-2*std(crossover.oxy_delta(idx)));
sig2_max = double(mean(crossover.oxy_delta(idx))+2*std(crossover.oxy_delta(idx)));
plot([sig2_min sig2_min],[0 1],'r','linewidth',1);
plot([sig2_max sig2_max],[0 1],'r','linewidth',1);
xlabel('Float [O_{2}] - Ship [O_{2}]');
xL = xlim; yL = ylim;
text(0.9*xL(1),0.9*yL(2),['Mean \Delta[O_{2}] = ' ...
    num2str(round(mean(crossover.oxy_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
text(0.9*xL(1),0.8*yL(2),['Std. \Delta[O_{2}] = ' ...
    num2str(round(std(crossover.oxy_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
exportgraphics(gcf,[pwd '/O2/Figures/Data/Ship_comp_histogram_uncorr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot corrected float vs. ship residuals
load(['O2/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
scatter(crossover.oxy_ship(idx),crossover.oxy_float_corr(idx),20,'.');
plot([0,450],[0 450],'k--');
ylabel('Corrected, Binned BGC Argo Oxygen Data (\mumol kg^{-1})');
xlabel('Binned GLODAP/WOD Oxygen Data (\mumol kg^{-1})');
crossover.err_md_corr = median(crossover.oxy_delta_corr(idx));
crossover.std_corr = std(crossover.oxy_delta_corr(idx));
crossover.rmse_corr = sqrt(mean((crossover.oxy_delta_corr(idx)).^2));
crossover.r2_corr = corr(crossover.oxy_ship(idx),crossover.oxy_float_corr(idx));
crossover.slope_corr = polyfit(crossover.oxy_ship(idx),crossover.oxy_float_corr(idx),1);
crossover.slope_corr = crossover.slope_corr(1);
text(300,120,['Med. Err. = ' num2str(round(crossover.err_md_corr,2))],'fontsize',12);
text(300,90,['RMSE = ' num2str(round(crossover.rmse_corr,1))],'fontsize',12);
text(300,60,['R^{2} = ' num2str(round(crossover.r2_corr,2))],'fontsize',12);
exportgraphics(gcf,['O2/Figures/Data/delta_vs_float_300_corr_' file_date float_file_ext '.png']);
close
% histogram
figure; hold on
histogram(crossover.oxy_delta_corr(idx),'Normalization','probability');
ylim('manual');
xlim([-30 30]);
set(gca,'fontsize',16);
sig2_min = double(mean(crossover.oxy_delta_corr(idx))-2*std(crossover.oxy_delta_corr(idx)));
sig2_max = double(mean(crossover.oxy_delta_corr(idx))+2*std(crossover.oxy_delta_corr(idx)));
plot([sig2_min sig2_min],[0 1],'r','linewidth',1);
plot([sig2_max sig2_max],[0 1],'r','linewidth',1);
xlabel('Float [O_{2}] - Ship [O_{2}]');
xL = xlim; yL = ylim;
text(0.9*xL(1),0.9*yL(2),['Mean \Delta[O_{2}] = ' ...
    num2str(round(mean(crossover.oxy_delta_corr(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
text(0.9*xL(1),0.8*yL(2),['Std. \Delta[O_{2}] = ' ...
    num2str(round(std(crossover.oxy_delta_corr(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
exportgraphics(gcf,[pwd '/O2/Figures/Data/Ship_comp_histogram_corr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot uncorrected float vs. ship residuals
load(['O2/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
pressures = unique(crossover.pres);
avg_delta_by_p = nan(size(pressures));
for p = 1:length(pressures)
    idx_p = crossover.pres == pressures(p);
    avg_delta_by_p(p) = mean(crossover.oxy_delta(idx_p),'omitnan');
end
plot(avg_delta_by_p,pressures);
set(gca,'YDir','reverse');

ylabel('Pressure (dbar)');
xlabel('Binned GLODAP/WOD Oxygen Data (\mumol kg^{-1})');
crossover.err_md = median(crossover.oxy_delta(idx));
crossover.std = std(crossover.oxy_delta(idx));
crossover.rmse = sqrt(mean((crossover.oxy_delta(idx)).^2));
crossover.r2 = corr(crossover.oxy_ship(idx),crossover.oxy_float(idx));
crossover.slope = polyfit(crossover.oxy_ship(idx),crossover.oxy_float(idx),1);
crossover.slope = crossover.slope(1);
text(300,120,['Med. Err. = ' num2str(round(crossover.err_md,2))],'fontsize',12);
text(300,90,['RMSE = ' num2str(round(crossover.rmse,1))],'fontsize',12);
text(300,60,['R^{2} = ' num2str(round(crossover.r2,2))],'fontsize',12);
exportgraphics(gcf,['O2/Figures/Data/delta_vs_float_300_uncorr_' file_date float_file_ext '.png']);
close
% histogram
figure; hold on
histogram(crossover.oxy_delta(idx),'Normalization','probability');
ylim('manual');
xlim([-30 30]);
set(gca,'fontsize',16);
sig2_min = double(mean(crossover.oxy_delta(idx))-2*std(crossover.oxy_delta(idx)));
sig2_max = double(mean(crossover.oxy_delta(idx))+2*std(crossover.oxy_delta(idx)));
plot([sig2_min sig2_min],[0 1],'r','linewidth',1);
plot([sig2_max sig2_max],[0 1],'r','linewidth',1);
xlabel('Float [O_{2}] - GLODAP/WOD [O_{2}]');
xL = xlim; yL = ylim;
text(0.9*xL(1),0.9*yL(2),['Mean \Delta[O_{2}] = ' ...
    num2str(round(mean(crossover.oxy_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
text(0.9*xL(1),0.8*yL(2),['Std. \Delta[O_{2}] = ' ...
    num2str(round(std(crossover.oxy_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
exportgraphics(gcf,[pwd '/O2/Figures/Data/Ship_comp_histogram_uncorr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot mapped crossover data by profile mean delta
load(['O2/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
lon_temp = convert_lon(crossover.lon,'0-360');
lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
% set up figure for uncorrected data
figure(1); hold on;
set(gcf,'visible','on','position',[100 100 1600 800]);
m_proj('robinson','lon',[20 380]);
m_coast('patch',rgb('gray'));
m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
clim([-30 30]);
colormap(cmocean('balance','pivot',0));
% for each crossover profile
float_profs = unique(crossover.id);
for p = 1:length(float_profs)
    idx = crossover.id == float_profs(p) & crossover.pres >= 300;
    m_scatter(mean(lon_temp(idx),'omitnan'),...
        mean(crossover.lat(idx),'omitnan'),50,...
        mean(crossover.oxy_delta(idx),'omitnan'),'filled',...
        'markeredgecolor','k');
end
colorbar;
if ~exist(['O2/Figures/Data'],'dir'); mkdir('pH/Figures/Data'); end
export_fig(gcf,'O2/Figures/Data/mapped_comparison_uncorr.png','-transparent');
close
% set up figure for corrected data
figure(1); hold on;
set(gcf,'visible','on','position',[100 100 1600 800]);
m_proj('robinson','lon',[20 380]);
m_coast('patch',rgb('gray'));
m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
clim([-30 30]);
colormap(cmocean('balance','pivot',0));
% for each crossover profile
float_profs = unique(crossover.id);
for p = 1:length(float_profs)
    idx = crossover.id == float_profs(p) & crossover.pres >= 300;
    m_scatter(mean(lon_temp(idx),'omitnan'),...
        mean(crossover.lat(idx),'omitnan'),50,...
        mean(crossover.oxy_delta_corr(idx),'omitnan'),'filled',...
        'markeredgecolor','k');
end
colorbar;
if ~exist(['O2/Figures/Data'],'dir'); mkdir('pH/Figures/Data'); end
export_fig(gcf,['O2/Figures/Data/mapped_comparison_corr.png'],'-transparent');
close
% clean up
clear crossover lon_temp

%% adjust and save float data
load(['O2/Data/float_corr_' file_date float_file_ext],'slp','int');
vars = fieldnames(float_data);
for v = 1:length(vars)
    if strcmp(vars{v},'OXY')
        % calculate float oxygen saturation
        Ts = log((298.15 - float_data.TEMP)./(273.15 + float_data.TEMP)); 
        % The coefficents below are from the second column of Table 1 of Garcia and
        % Gordon (1992)
        a0 =  5.80871; 
        a1 =  3.20291;
        a2 =  4.17887;
        a3 =  5.10006;
        a4 = -9.86643e-2;
        a5 =  3.80369;
        b0 = -7.01577e-3;
        b1 = -7.70028e-3;
        b2 = -1.13864e-2;
        b3 = -9.51519e-3;
        c0 = -2.75915e-7;
        float_data.OXY_SAT = exp(a0 + Ts.*(a1 + Ts.*(a2 + Ts.*(a3 + Ts.*(a4 + a5*Ts)))) ...
                  + float_data.SAL.*(b0 + Ts.*(b1 + Ts.*(b2 + b3*Ts)) + c0*float_data.SAL));
        float_data.OXY_SAT_PER = 100*(float_data.OXY./float_data.OXY_SAT);
        % apply linear correction to float O2
        float_data_adjusted.OXY = double(((float_data.OXY_SAT_PER - ...
            (slp .* float_data.OXY_SAT + int))./100) .* float_data.OXY_SAT);
    else
        float_data_adjusted.(vars{v}) = float_data.(vars{v});
    end
end

% save
if ~exist([pwd '/O2/Data'],'dir'); mkdir('O2/Data'); end
save(['O2/Data/processed_float_o2_data_adjusted_' file_date float_file_ext '.mat'],...
    'float_data_adjusted','file_date','-v7.3');
clear slp int float_data float_data_adjusted v vars

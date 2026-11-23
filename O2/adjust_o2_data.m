% adjust_o2_float_data
%
% DESCRIPTION:
% This function bins float and glodap data, co-locates corresponding bins,
% compares the two datasets, and calculates a bulk correction factor based
% on the discrepancies between the two.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 07/24/2026

function adjust_o2_data(param_props,float_file_ext,glodap_vrs,snap_date,...
    include_float,include_glodap,include_osd,include_ctd)

%% define dataset extensions
if include_float == 1; float_ext = 'f'; else float_ext = ''; end
if include_glodap == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if include_osd == 1; osd_ext = 'o'; else osd_ext = ''; end
if include_ctd == 1; ctd_ext = 'c'; else ctd_ext = ''; end

%% define year
todays_date = datevec(datetime('today'));
end_year = todays_date(1);

%% load interpolated float and glodap data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
if include_float; load([param_props.data_path 'Data/processed_float_o2_data_' ...
        file_date float_file_ext '.mat'],'float_data','file_date'); end
if include_glodap; load([param_props.data_path 'Data/processed_glodap_o2_data_' ...
        glodap_vrs '.mat'],'glodap_data'); end
if include_osd || include_ctd; load([param_props.data_path 'Data/processed_wod_o2_data_' ...
        num2str(end_year) '.mat'],'wod_data'); end

%% only do if file doesn't exist
if exist([param_props.data_path 'Data/crossover_data_' file_date float_file_ext '.mat'],'file') ~= 2

%% remove data points based on global range test
% 480 umol/kg used as max for WOD range checks
if include_float
    idx_rem = float_data.OXY < 0 | float_data.OXY > 480;
    disp([num2str(sum(idx_rem)) ' float data points removed by global range test (' ...
        num2str(100*(sum(idx_rem)/length(float_data.PROF_ID))) ' % of data)']);
    vars = fieldnames(float_data);
    for v = 1:length(vars)
        float_data.(vars{v})(idx_rem) = [];
    end
end

if include_glodap
    idx_rem = glodap_data.OXY < 0 | glodap_data.OXY > 480;
    disp([num2str(sum(idx_rem)) ' glodap data points removed by global range test (' ...
        num2str(100*(sum(idx_rem)/length(glodap_data.ID))) ' % of data)']);
    vars = fieldnames(glodap_data);
    for v = 1:length(vars)
        glodap_data.(vars{v})(idx_rem) = [];
    end
end

if include_osd || include_ctd
    idx_rem = wod_data.OXY < 0 | wod_data.OXY > 480;
    disp([num2str(sum(idx_rem)) ' wod data points removed by global range test (' ...
        num2str(100*(sum(idx_rem)/length(wod_data.ID))) ' % of data)']);
    vars = fieldnames(wod_data);
    for v = 1:length(vars)
        wod_data.(vars{v})(idx_rem) = [];
    end
end

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
temp_path = [param_props.ancillary_data_path 'WOA/TEMPERATURE/'];
sal_path = [param_props.ancillary_data_path 'WOA/SALINITY/'];
oxy_path = [param_props.ancillary_data_path 'WOA/OXYGEN/'];
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
if include_float

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
    clear idx idx_lon idx_lat idx_pres idx_mnth i
    
    %% plot differences between float data and WOA
    % histogram of differences
    figure; hold on
    histogram(WOA_delta);
    set(gca,'fontsize',16);
    ax = gca;
    plot([mean_delta+3*st_dev_delta mean_delta+3*st_dev_delta],ax.YLim,'r','linewidth',2);
    plot([mean_delta-3*st_dev_delta mean_delta-3*st_dev_delta],ax.YLim,'r','linewidth',2);
    xlabel('Float [O_{2}] - WOA [O_{2}]');
    exportgraphics(gcf,[param_props.data_path 'Figures/Data/WOA_comp_histogram_' file_date float_file_ext '.png']);
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
    exportgraphics(gcf,[param_props.data_path 'Figures/Data/WOA_Float_comp_scatter_' file_date float_file_ext '.png']);
    close
    % clean up
    clear counts bin_centers c h myColorMap
    
    %% remove float data points more than 3 sigmas from WOA value for each depth
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
    disp([num2str(sum(idx_rem)) ' float data points removed by WOA comparison test (' ...
        num2str(100*(sum(idx_rem)/length(float_data.PROF_ID))) ' % of data)']);
    vars = fieldnames(float_data);
    % plot(mean_delta_pres,pres_levels); xlabel('Float - WOA'); ylabel('pres.');
    for v = 1:length(vars)
        float_data.(vars{v})(idx_rem) = [];
    end
    WOA_delta(idx_rem) = [];
    WOA_delta_per(idx_rem) = [];
    WOA_match(idx_rem) = [];
    
    %% histogram of float differences (corrected)
    figure; hold on
    histogram(WOA_delta);
    set(gca,'fontsize',16);
    ax = gca;
    plot([mean_delta+3*st_dev_delta mean_delta+3*st_dev_delta],ax.YLim,'r','linewidth',2);
    plot([mean_delta-3*st_dev_delta mean_delta-3*st_dev_delta],ax.YLim,'r','linewidth',2);
    xlabel('Float [O_{2}] - WOA [O_{2}]');
    exportgraphics(gcf,[param_props.data_path 'Figures/Data/WOA_Float_comp_histogram_' file_date float_file_ext '.png']);
    close
    % clean up
    clear WOA_match WOA_delta WOA_delta_per mean_delta st_dev_delta idx_rem v vars

end

% prof_ids = unique(float_data.PROF_ID);
% for p = 1:length(prof_ids)
%     idx = prof_ids(p) == float_data.PROF_ID;
%     figure;
%     plot(float_data.OXY(idx),float_data.PRES(idx));
%     set(gca,'ydir','reverse');
%     keyboard
%     close
% end

%% compare WOD data to WOA
if include_osd || include_ctd
% index to valid data points
idx = find(~isnan(wod_data.OXY) & ~isnan(wod_data.LAT) & ...
    ~isnan(wod_data.LON) & ~isnan(wod_data.PRES) & ...
    ~isnan(wod_data.TIME));
% obtain month from date
date = datevec(datenum(wod_data.YEAR,0,wod_data.DAY));
wod_data.MONTH = date(:,2);
clear date
% pre-allocate matches
WOA_match = nan(size(wod_data.OXY));
% match float data with WOA
for i=1:length(idx)
    idx_lon = find(min(abs(WOA.LON - wod_data.LON(idx(i)))) == ...
        abs(WOA.LON - wod_data.LON(idx(i))));
    idx_lat = find(min(abs(WOA.LAT - wod_data.LAT(idx(i)))) == ...
        abs(WOA.LAT - wod_data.LAT(idx(i))));
    idx_pres = find(min(abs(squeeze(WOA.PRES(idx_lon,idx_lat,:)) - wod_data.PRES(idx(i)))) == ...
        abs(squeeze(WOA.PRES(idx_lon,idx_lat,:)) - wod_data.PRES(idx(i))));
    idx_mnth = find(min(abs((1:12)' - wod_data.MONTH(idx(i)))) == ...
        abs((1:12)' - wod_data.MONTH(idx(i))));
    WOA_match(idx(i)) = WOA.OXY(idx_lon(1),idx_lat(1),idx_pres(1),idx_mnth(1));
end
% calculate differences
WOA_delta = wod_data.OXY-WOA_match;
WOA_delta_per = WOA_delta./wod_data.OXY;
WOA_delta_per(isinf(WOA_delta_per))=NaN;
st_dev_delta = std(WOA_delta,[],'omitnan');
mean_delta = mean(WOA_delta,'omitnan');
% clean up
clear idx idx_lon idx_lat idx_pres idx_mnth WOA i

%% plot differences between WOD data and WOA
% histogram of differences
figure; hold on
histogram(WOA_delta);
set(gca,'fontsize',16);
ax = gca;
plot([mean_delta+3*st_dev_delta mean_delta+3*st_dev_delta],ax.YLim,'r','linewidth',2);
plot([mean_delta-3*st_dev_delta mean_delta-3*st_dev_delta],ax.YLim,'r','linewidth',2);
xlabel('WOD [O_{2}] - WOA [O_{2}]');
exportgraphics(gcf,[param_props.data_path 'Figures/Data/WOA_vs_WOD_comp_histogram_' file_date float_file_ext '.png']);
close
% scatter of delta values
figure; hold on
set(gca,'fontsize',16);
xlabel('WOA [O_{2}]');
ylabel('WOD [O_{2}] - WOA [O_{2}]');
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
exportgraphics(gcf,[param_props.data_path 'Figures/Data/WOA_vs_WOD_comp_scatter_' file_date float_file_ext '.png']);
close
% clean up
clear counts bin_centers c h myColorMap

%% remove WOD data points more than 3 sigmas from WOA value for each depth
pres_levels = unique(wod_data.PRES);
idx_rem = false(length(WOA_delta),1);
st_dev_delta_pres = nan(size(pres_levels));
mean_delta_pres = nan(size(pres_levels));
for z = 1:length(pres_levels)
    pres_idx = wod_data.PRES == pres_levels(z);
    st_dev_delta_pres(z) = std(WOA_delta(pres_idx),[],'omitnan');
    mean_delta_pres(z) = mean(WOA_delta(pres_idx),'omitnan');
    idx_rem(pres_idx) = WOA_delta(pres_idx) > mean_delta_pres(z)+3.*st_dev_delta_pres(z) | ...
        WOA_delta(pres_idx) < mean_delta_pres(z)-3.*st_dev_delta_pres(z);
end
disp([num2str(sum(idx_rem)) ' WOD data points removed by WOA comparison test (' ...
    num2str(100*(sum(idx_rem)/length(wod_data.ID))) ' % of data)']);
vars = fieldnames(wod_data);
% plot(mean_delta_pres,pres_levels); xlabel('Float - WOA'); ylabel('pres.');
for v = 1:length(vars)
    wod_data.(vars{v})(idx_rem) = [];
end
WOA_delta(idx_rem) = [];
WOA_delta_per(idx_rem) = [];
WOA_match(idx_rem) = [];

%% histogram of WOD differences (corrected)
figure; hold on
histogram(WOA_delta);
set(gca,'fontsize',16);
ax = gca;
plot([mean_delta+3*st_dev_delta mean_delta+3*st_dev_delta],ax.YLim,'r','linewidth',2);
plot([mean_delta-3*st_dev_delta mean_delta-3*st_dev_delta],ax.YLim,'r','linewidth',2);
xlabel('WOD [O_{2}] - WOA [O_{2}]');
exportgraphics(gcf,[param_props.data_path 'Figures/Data/WOA_WOD_comp_histogram_' file_date float_file_ext '.png']);
close
% clean up
clear WOA_match WOA_delta WOA_delta_per mean_delta st_dev_delta idx_rem v vars
end

% prof_ids = unique(wod_data.PROF_ID);
% for p = 1:length(prof_ids)
%     idx = prof_ids(p) == wod_data.PROF_ID;
%     figure;
%     plot(wod_data.OXY(idx),wod_data.PRES(idx));
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
% include WOD data or not
if include_ctd || include_osd
    ship_lon = [glodap_data.LON;wod_data.LON];
    ship_lat = [glodap_data.LAT;wod_data.LAT];
    ship_time = [glodap_data.TIME;wod_data.TIME];
    ship_sigma = [glodap_data.SIGMA;wod_data.SIGMA];
    ship_pres = [glodap_data.PRES;wod_data.PRES];
    ship_oxy = [glodap_data.OXY;wod_data.OXY];
else
    ship_lon = [glodap_data.LON];
    ship_lat = [glodap_data.LAT];
    ship_time = [glodap_data.TIME];
    ship_sigma = [glodap_data.SIGMA];
    ship_pres = [glodap_data.PRES];
    ship_oxy = [glodap_data.OXY];
end
% cycle through float profiles 
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
        % extract float data
        oxy_temp = float_data.OXY(float_data.PROF_ID==float_profile_IDs(p));
        pres_temp = float_data.PRES(float_data.PROF_ID==float_profile_IDs(p));
        sigma_temp = float_data.SIGMA(float_data.PROF_ID==float_profile_IDs(p));
        lon_temp = float_data.LON(float_data.PROF_ID==float_profile_IDs(p));
        lat_temp = float_data.LAT(float_data.PROF_ID==float_profile_IDs(p));
        temp_temp = float_data.TEMP(float_data.PROF_ID==float_profile_IDs(p));
        sal_temp = float_data.SAL(float_data.PROF_ID==float_profile_IDs(p));
        % combine float profiles (helpful if there are more than one or
        % multiple measurements at the same depths)
        oxy_float = nan(length(pres_levels),1);
        pres_float = nan(length(pres_levels),1);
        sigma_float = nan(length(pres_levels),1);
        lon_float = nan(length(pres_levels),1);
        lat_float = nan(length(pres_levels),1);
        temp_float = nan(length(pres_levels),1);
        sal_float = nan(length(pres_levels),1);
        for z = 1:length(pres_levels)
            idx = pres_temp == pres_levels(z);
            oxy_float(z) = mean(oxy_temp(idx),'omitnan');
            pres_float(z) = mean(pres_temp(idx),'omitnan');
            sigma_float(z) = mean(sigma_temp(idx),'omitnan');
            lon_float(z) = mean(lon_temp(idx),'omitnan');
            lat_float(z) = mean(lat_temp(idx),'omitnan');
            temp_float(z) = mean(temp_temp(idx),'omitnan');
            sal_float(z) = mean(sal_temp(idx),'omitnan');
        end
        % get oxygen gradient
        oxy_float_gradient = [0;diff(oxy_float)];
        % extract ship data
        pres_temp = ship_pres(idx_all);
        oxy_temp = ship_oxy(idx_all);
        sigma_temp = ship_sigma(idx_all);
        % combine ship profiles (helpful if there are more than one or
        % multiple measurements at the same depths)
        oxy_ship = nan(length(pres_levels),1);
        pres_ship = nan(length(pres_levels),1);
        sigma_ship = nan(length(pres_levels),1);
        for z = 1:length(pres_levels)
            idx = pres_temp == pres_levels(z);
            pres_ship(z) = mean(pres_temp(idx),'omitnan');
            oxy_ship(z) = mean(oxy_temp(idx),'omitnan');
            sigma_ship(z) = mean(sigma_temp(idx),'omitnan');
        end
        % oxy_ship = oxy_temp(~isnan(oxy_temp));
        % pres_ship = pres_temp(~isnan(pres_temp));
        % sigma_ship = sigma_temp(~isnan(sigma_temp));
        % get oxygen gradient
        oxy_ship_gradient = [0;diff(oxy_ship)];
        % plot matched profiles
        % figure; plot(oxy_float,pres_levels(idx_pres),oxy_ship,pres_levels(idx_pres)); set(gca,'YDir','reverse'); close;
        % log matched data
        crossover.id = [crossover.id;repmat(float_profile_IDs(p),length(pres_float),1)];
        crossover.lon = [crossover.lon;lon_float];
        crossover.lat = [crossover.lat;lat_float];
        crossover.pres = [crossover.pres;pres_float];
        crossover.oxy_float = [crossover.oxy_float;oxy_float];
        crossover.oxy_float_gradient = [crossover.oxy_float_gradient;oxy_float_gradient];
        crossover.temp_float = [crossover.temp_float;temp_float];
        crossover.sal_float = [crossover.sal_float;sal_float];
        crossover.oxy_ship = [crossover.oxy_ship;oxy_ship];
        crossover.oxy_ship_gradient = [crossover.oxy_ship_gradient;oxy_ship_gradient];
        crossover.sigma_float = [crossover.sigma_float;sigma_float];
        crossover.sigma_ship = [crossover.sigma_ship;sigma_ship];
        try
        crossover.oxy_delta = [crossover.oxy_delta;oxy_float-oxy_ship];
        catch
            keyboard
        end
        crossover.oxy_delta_per = [crossover.oxy_delta_per;...
            (oxy_float-oxy_ship)./oxy_float];
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

% clean up
clear pres_levels float_profile_IDs lon lat time oxy_float pres_float
clear idx_lon idx_lat idx_time oxy_temp pres_temp idx_pres oxy_ship

% save crossover data
if ~exist([param_props.data_path 'Data'],'dir'); mkdir([param_props.data_path 'Data']); end
save([param_props.data_path 'Data/crossover_data_' file_date float_file_ext],'crossover','idx','-v7.3')


else

disp(['Float data crossovers (' erase(float_file_ext,'_') ') for ' file_date ' already determined.']);

end

%% develop O2 correction
% load crossover data
load([param_props.data_path 'Data/crossover_data_' file_date float_file_ext],'crossover','idx')

% % plot deltas spatially
% pressures = unique(crossover.pres);
% ids = unique(crossover.id);
% for p = 40%:length(pressures)
%     figure('Position',[100 100 1200 800])
%     worldmap([-90 90],[-180 180]); hold on;
%     mlabel off; plabel off;
%     lon_to_plot = NaN(size(ids));
%     lat_to_plot = NaN(size(ids));
%     delta_to_plot = NaN(size(ids));
%     for i = 1:length(ids)
%         idx = find(crossover.pres == pressures(p) & crossover.id == ids(i));
%         if ~isempty(idx)
%         lon_to_plot(i) = crossover.lon(idx);
%         lat_to_plot(i) = crossover.lat(idx);
%         delta_to_plot(i) = crossover.oxy_delta(idx);
%         end
%     end
%     scatterm(lat_to_plot,lon_to_plot,20,...
%         delta_to_plot,'filled');
%     plot_land('map');
%     colorbar;
%     caxis([-30 30]);
%     colormap(cmocean('balance'));
% end

%% use constant offset
% if strcmp(float_file_ext,'_D')
%     slp = 0;
%     int = -1.69;
% elseif strcmp(float_file_ext,'_D_A')
%     % load([param_props.data_path 'Data/processed_float_o2_data_' ...
%     %     file_date float_file_ext '.mat'],'float_data','file_date');
%     slp = 0;
%     int_d = -1.69;
%     int_a = -4.69;
%     int = int_d*0.5 + int_a*0.5;
% end
% % apply linear correction to float O2 (with absolute O2)
% corr_fac = slp .* crossover.oxy_float + int; % modelled saturation delta
% crossover.oxy_float_corr = crossover.oxy_float - corr_fac;

%% inspect the distriibution of the oxygen deltas
% Keep finite float/ship pairs with a finite pressure and offset
valid = isfinite(crossover.oxy_float) & ...
    isfinite(crossover.oxy_ship) & ...
    isfinite(crossover.pres) & ...
    isfinite(crossover.oxy_delta);

delta = double(crossover.oxy_delta(valid));
pressure = double(crossover.pres(valid));

if isempty(delta)
    error('No valid crossover offsets are available.');
end

% Robust summary
q = prctile(delta, [1 5 25 50 75 95 99]);
fprintf('N = %d\n', numel(delta));
fprintf('Mean = %.2f, SD = %.2f umol/kg\n', mean(delta), std(delta));
fprintf('Median = %.2f, IQR = %.2f umol/kg\n', q(4), q(5)-q(3));
fprintf('Percentiles [1 5 25 50 75 95 99]:\n');
disp(q);

% Summaries by pressure level
levels = unique(pressure);
n = zeros(size(levels));
med = nan(size(levels));
lo = nan(size(levels));
hi = nan(size(levels));

for k = 1:numel(levels)
    d = delta(pressure == levels(k));
    n(k) = numel(d);
    med(k) = median(d);
    pct = prctile(d, [25 75]);
    lo(k) = pct(1);
    hi(k) = pct(2);
end

% Plot global distribution and pressure dependence
figure;
tiledlayout(2,2);

nexttile;
histogram(delta, 60);
xline(q(4), 'r-', 'Median');
xlabel('Float - ship O_2 (\mumol kg^{-1})');
ylabel('Count');
title('All valid crossovers');

nexttile;
histogram(delta, 60);
xlim(q([1 7]));  % zoom to the 1st–99th percentile range
xline(q(4), 'r-', 'Median');
xlabel('Float - ship O_2 (\mumol kg^{-1})');
ylabel('Count');
title('Central 98% of offsets');

nexttile;
errorbar(levels, med, med-lo, hi-med, 'o-');
yline(0, 'k:');
xlabel('Pressure (dbar)');
ylabel('Median offset and IQR (\mumol kg^{-1})');
title('Offset by pressure level');

nexttile;
semilogy(levels, n, 'o-');
xlabel('Pressure (dbar)');
ylabel('Number of crossovers');
title('Sample count by pressure');

% save figure
exportgraphics(gcf,[param_props.data_path 'Figures/Data/corr_summary_' file_date float_file_ext '.png']);
close

%% develop a depth-dependent correction (using unified equation with a delta maximum)
% index
idx = ~isnan(crossover.oxy_float) & ~isnan(crossover.oxy_ship);
% determine average delta by pressure level
pressures = unique(crossover.pres(idx));
avg_delta_by_p = nan(size(pressures));
avg_by_p = nan(size(pressures));
for p = 1:length(pressures)
    idx_p = crossover.pres == pressures(p) & idx;
    avg_delta_by_p(p) = mean(crossover.oxy_delta(idx_p),'omitnan');
    avg_by_p(p) = mean(crossover.oxy_ship(idx_p),'omitnan');
end
% find maximum crossover offset
valid_delta = isfinite(pressures) & isfinite(avg_delta_by_p);
if ~any(valid_delta)
    error('adjust_o2_data:NoValidPressureDeltas', ...
        'Cannot fit a pressure correction without finite oxygen deltas.');
end
candidate_pressures = pressures(valid_delta);
candidate_deltas = avg_delta_by_p(valid_delta);
[~,idx_max_delta] = max(abs(candidate_deltas));
max_delta_pressure = candidate_pressures(idx_max_delta);
disp(['Maximum crossover offset at ' num2str(max_delta_pressure) ' dbar.'])
correction_model = fit_pressure_correction(pressures, avg_delta_by_p, ...
    max_delta_pressure);
fit_est = apply_pressure_correction(correction_model, pressures);
crossover.corr_fac = apply_pressure_correction(correction_model, crossover.pres);
crossover.oxy_float_corr = crossover.oxy_float - crossover.corr_fac;
% plot offset as a function of pressure
figure; hold on;
plot(avg_delta_by_p,pressures,'k','linewidth',2);
plot(fit_est,pressures,'g-','linewidth',2);
set(gca,'ydir','reverse');
xlabel('Float [O_{2}] - Ship [O_{2}]');
ylabel('Pressure (dbar)');
legend('Average \Delta[O_{2}]','Breakpoint fit');
exportgraphics(gcf,[param_props.data_path 'Figures/Data/delta_vs_pressure_' file_date float_file_ext '.png']);
close

%% fit delta against [O2] solubility
% mdl = fitlm(crossover.oxy_sat_float(idx),crossover.oxy_delta_sat_per(idx),'Intercept',true);
% slp = mdl.Coefficients.Estimate(2);
% int = mdl.Coefficients.Estimate(1);

% apply linear correction to float O2 (with per sat)
% corr_fac = slp .* crossover.oxy_sat_float + int; % modelled saturation delta
% crossover.oxy_float_sat_per_corr = crossover.oxy_float_sat_per - corr_fac;
% crossover.oxy_float_corr = (crossover.oxy_float_sat_per_corr./100).*crossover.oxy_sat_float;

%% re-calculate deltas
crossover.oxy_delta_corr = crossover.oxy_float_corr-crossover.oxy_ship;
crossover.oxy_delta_per_corr = (crossover.oxy_float_corr-crossover.oxy_ship)./crossover.oxy_float_corr;

%% save crossover data and correction factors
if ~exist([param_props.data_path 'Data'],'dir'); mkdir([param_props.data_path 'Data']); end
% save(['O2/Data/float_corr_' file_date float_file_ext],'slp','int');
% clear slp int
save([param_props.data_path 'Data/crossover_data_' file_date float_file_ext],'crossover','idx','-v7.3')
slp_shal = correction_model.coefficients(2);
slp_deep = correction_model.coefficients(3);
int_shal = correction_model.coefficients(1) - slp_shal*max_delta_pressure;
int_deep = correction_model.coefficients(1) - slp_deep*max_delta_pressure;
save([param_props.data_path 'Data/float_corr_' file_date float_file_ext], ...
    'correction_model','slp_deep','int_deep','slp_shal','int_shal', ...
    'max_delta_pressure');

%% plot uncorrected float vs. ship residuals
load([param_props.data_path 'Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
ylabel('Binned BGC Argo Oxygen Data (\mumol kg^{-1})');
xlabel('Binned GLODAP/WOD Oxygen Data (\mumol kg^{-1})');
[counts,bin_centers] = hist3([crossover.oxy_ship(idx),...
    crossover.oxy_float(idx)],'Edges',{0:5:460 0:5:460});
h=pcolor(bin_centers{1}-mean(diff(bin_centers{1}))/2,...
    bin_centers{2}-mean(diff(bin_centers{2}))/2,counts');
plot([0,450],[0 450],'k--');
set(h,'EdgeColor','none');
xlim([-0.5 460]); ylim([-0.5 460]);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log')
clim([1e0 5e3]);
c=colorbar;
c.Label.String = 'log10(Bin Counts)';
crossover.err_md = median(crossover.oxy_delta(idx));
crossover.std = std(crossover.oxy_delta(idx));
crossover.rmse = sqrt(mean((crossover.oxy_delta(idx)).^2));
crossover.r2 = corr(crossover.oxy_ship(idx),crossover.oxy_float(idx));
crossover.slope = polyfit(crossover.oxy_ship(idx),crossover.oxy_float(idx),1);
crossover.slope = crossover.slope(1);
text(20,420,['Med. Err. = ' num2str(round(crossover.err_md,2))],'fontsize',12);
text(20,390,['RMSE = ' num2str(round(crossover.rmse,1))],'fontsize',12);
text(20,360,['R^{2} = ' num2str(round(crossover.r2,2))],'fontsize',12);
exportgraphics(gcf,[param_props.data_path 'Figures/Data/delta_vs_float_300_uncorr_' file_date float_file_ext '.png']);
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
exportgraphics(gcf,[param_props.data_path 'Figures/Data/Ship_comp_histogram_uncorr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot corrected float vs. ship residuals
load([param_props.data_path 'Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scattered 2D histogram
figure; hold on;
ylabel('Adjusted, Binned BGC Argo Oxygen Data (\mumol kg^{-1})');
xlabel('Binned GLODAP/WOD Oxygen Data (\mumol kg^{-1})');
[counts,bin_centers] = hist3([crossover.oxy_ship(idx),...
    crossover.oxy_float_corr(idx)],'Edges',{0:5:460 0:5:460});
h=pcolor(bin_centers{1}-mean(diff(bin_centers{1}))/2,...
    bin_centers{2}-mean(diff(bin_centers{2}))/2,counts');
plot([0,450],[0 450],'k--');
set(h,'EdgeColor','none');
xlim([-0.5 460]); ylim([-0.5 460]);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log')
clim([1e0 5e3]);
c=colorbar;
c.Label.String = 'log10(Bin Counts)';
crossover.err_md_corr = median(crossover.oxy_delta_corr(idx));
crossover.std_corr = std(crossover.oxy_delta_corr(idx));
crossover.rmse_corr = sqrt(mean((crossover.oxy_delta_corr(idx)).^2));
crossover.r2_corr = corr(crossover.oxy_ship(idx),crossover.oxy_float_corr(idx));
crossover.slope_corr = polyfit(crossover.oxy_ship(idx),crossover.oxy_float_corr(idx),1);
crossover.slope_corr = crossover.slope_corr(1);
text(20,420,['Med. Err. = ' num2str(round(crossover.err_md_corr,2))],'fontsize',12);
text(20,390,['RMSE = ' num2str(round(crossover.rmse_corr,1))],'fontsize',12);
text(20,360,['R^{2} = ' num2str(round(crossover.r2_corr,2))],'fontsize',12);
exportgraphics(gcf,[param_props.data_path 'Figures/Data/delta_vs_float_300_corr_' file_date float_file_ext '.png']);
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
exportgraphics(gcf,[param_props.data_path 'Figures/Data/Ship_comp_histogram_corr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot float vs. ship residual profiles
load([param_props.data_path 'Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% uncorrected
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
exportgraphics(gcf,[param_props.data_path 'Figures/Data/delta_vs_float_by_pressure_uncorr_' file_date float_file_ext '.png']);
close
% corrected
figure; hold on;
pressures = unique(crossover.pres);
avg_delta_by_p = nan(size(pressures));
for p = 1:length(pressures)
    idx_p = crossover.pres == pressures(p);
    avg_delta_by_p(p) = mean(crossover.oxy_delta_corr(idx_p),'omitnan');
end
plot(avg_delta_by_p,pressures);
set(gca,'YDir','reverse');
ylabel('Pressure (dbar)');
xlabel('Binned GLODAP/WOD Oxygen Data (\mumol kg^{-1})');
crossover.err_md = median(crossover.oxy_delta_corr(idx));
crossover.std = std(crossover.oxy_delta_corr(idx));
crossover.rmse = sqrt(mean((crossover.oxy_delta_corr(idx)).^2));
crossover.r2 = corr(crossover.oxy_ship(idx),crossover.oxy_float(idx));
crossover.slope = polyfit(crossover.oxy_ship(idx),crossover.oxy_float(idx),1);
crossover.slope = crossover.slope(1);
text(300,120,['Med. Err. = ' num2str(round(crossover.err_md,2))],'fontsize',12);
text(300,90,['RMSE = ' num2str(round(crossover.rmse,1))],'fontsize',12);
text(300,60,['R^{2} = ' num2str(round(crossover.r2,2))],'fontsize',12);
exportgraphics(gcf,[param_props.data_path 'Figures/Data/delta_vs_float_by_pressure_corr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot mapped crossover data by profile mean delta
load([param_props.data_path 'Data/crossover_data_' file_date float_file_ext],'crossover','idx');
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
if ~exist([param_props.data_path 'Figures/Data'],'dir'); mkdir([param_props.data_path 'Figures/Data']); end
export_fig(gcf,[param_props.data_path 'Figures/Data/mapped_comparison_uncorr.png'],'-transparent');
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
if ~exist([param_props.data_path 'Figures/Data'],'dir'); mkdir([param_props.data_path 'Figures/Data']); end
export_fig(gcf,[param_props.data_path 'Figures/Data/mapped_comparison_corr.png'],'-transparent');
close
% clean up
clear crossover lon_temp

%% adjust and save float data
load([param_props.data_path 'Data/float_corr_' file_date float_file_ext], ...
    'correction_model');
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

        % apply linear correction to float O2 (with per sat)
        % float_data_adjusted.OXY = double(((float_data.OXY_SAT_PER - ...
         %    (slp .* float_data.OXY_SAT + int))./100) .* float_data.OXY_SAT);

        % apply linear correction to float O2 (with absolute O2)
        % corr_fac = slp .* float_data.OXY + int; % modelled [O2] delta
        % float_data_adjusted.OXY = float_data.OXY - corr_fac;

        % apply the fitted pressure-dependent correction
        corr_fac = apply_pressure_correction(correction_model, float_data.PRES);
        float_data_adjusted.OXY = float_data.OXY - corr_fac;

    else
        float_data_adjusted.(vars{v}) = float_data.(vars{v});
    end
end

if include_float
    % save adjusted float data
    if ~exist([param_props.data_path 'Data'],'dir'); mkdir([param_props.data_path 'Data']); end
    save([param_props.data_path 'Data/processed_float_o2_data_adjusted_' ...
        file_date float_file_ext '.mat'],...
        'float_data_adjusted','file_date','-v7.3');
    clear slp int float_data float_data_adjusted v vars
end

if include_osd || include_ctd
    % save adjusted WOD data
    if ~exist([param_props.data_path 'Data'],'dir'); mkdir([param_props.data_path 'Data']); end
    save([param_props.data_path 'Data/processed_wod_o2_data_adjusted_' ...
        num2str(end_year) '.mat'],...
        'wod_data','file_date','-v7.3');
end

if include_glodap
    % save adjusted GLODAP data
    if ~exist([param_props.data_path 'Data'],'dir'); mkdir([param_props.data_path 'Data']); end
    save([param_props.data_path 'Data/processed_glodap_o2_data_adjusted_' ...
        glodap_vrs '.mat'],'glodap_data','file_date','-v7.3');
end

end

function model = fit_pressure_correction(pressures, oxy_deltas, breakpoint)
if ~isvector(pressures) || ~isvector(oxy_deltas) || ...
        numel(pressures) ~= numel(oxy_deltas)
    error('adjust_o2_data:InvalidFitData', ...
        'Pressures and oxygen deltas must be vectors of equal length.');
end
if ~isscalar(breakpoint) || ~isfinite(breakpoint)
    error('adjust_o2_data:InvalidBreakpoint', ...
        'The pressure breakpoint must be a finite scalar.');
end

pressures = pressures(:);
oxy_deltas = oxy_deltas(:);
valid = isfinite(pressures) & isfinite(oxy_deltas);
pressures = pressures(valid);
oxy_deltas = oxy_deltas(valid);
if ~any(pressures < breakpoint) || ~any(pressures > breakpoint)
    error('adjust_o2_data:BreakpointNeedsBothSides', ...
        'The breakpoint must have valid pressure data on both sides.');
end

d = pressures - breakpoint;
design = [ones(size(d)), min(d,0), max(d,0)];
if rank(design) < 3
    error('adjust_o2_data:InsufficientFitData', ...
        'The pressure data do not support fitting both segment slopes.');
end

% Constrain the profile to have a minimum at the breakpoint.
constraints = [0 1 0; 0 0 -1];
coefficients = lsqlin(design, oxy_deltas, constraints, [0; 0]);
model = struct('breakpoint', breakpoint, ...
    'coefficients', coefficients(:));
end

function correction = apply_pressure_correction(model, input_pressures)
if ~isstruct(model) || ~isfield(model,'breakpoint') || ...
        ~isfield(model,'coefficients') || numel(model.coefficients) ~= 3
    error('adjust_o2_data:InvalidCorrectionModel', ...
        'Expected a model returned by fit_pressure_correction.');
end

d = input_pressures - model.breakpoint;
coefficients = model.coefficients(:);
correction = coefficients(1) + coefficients(2).*min(d,0) + ...
    coefficients(3).*max(d,0);
end
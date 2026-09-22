% adjust_no3_float_data
%
% DESCRIPTION:
% This function bins float and glodap data, co-locates corresponding bins,
% compares the two datasets, and calculates a bulk correction factor based
% on the discrepancies between the two.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 07/15/2026

function adjust_no3_data(float_file_ext,glodap_vrs,snap_date,...
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
if include_float; load(['NO3/Data/processed_float_no3_data_' ...
        file_date float_file_ext '.mat'],'float_data','file_date'); end
if include_glodap; load(['NO3/Data/processed_glodap_no3_data_' ...
        glodap_vrs '.mat'],'glodap_data'); end
if include_osd || include_ctd; load(['NO3/Data/processed_wod_no3_data_' ...
        num2str(end_year) '.mat'],'wod_data'); end

%% remove data points based on global range test
% 100 umol/kg used as max for WOD range checks
if include_float
    idx_rem = float_data.NIT < 0 | float_data.NIT > 100;
    disp([num2str(sum(idx_rem)) ' float data points removed by global range test (' ...
        num2str(100*(sum(idx_rem)/length(float_data.PROF_ID))) ' % of data)']);
    vars = fieldnames(float_data);
    for v = 1:length(vars)
        float_data.(vars{v})(idx_rem) = [];
    end
end

if include_glodap
    idx_rem = glodap_data.NIT < 0 | glodap_data.NIT > 100;
    disp([num2str(sum(idx_rem)) ' glodap data points removed by global range test (' ...
        num2str(100*(sum(idx_rem)/length(glodap_data.ID))) ' % of data)']);
    vars = fieldnames(glodap_data);
    for v = 1:length(vars)
        glodap_data.(vars{v})(idx_rem) = [];
    end
end

if include_osd || include_ctd
    idx_rem = wod_data.NIT < 0 | wod_data.NIT > 100;
    disp([num2str(sum(idx_rem)) ' wod data points removed by global range test (' ...
        num2str(100*(sum(idx_rem)/length(wod_data.ID))) ' % of data)']);
    vars = fieldnames(wod_data);
    for v = 1:length(vars)
        wod_data.(vars{v})(idx_rem) = [];
    end
end

%% import WOA climatologies
temp_path = [pwd '/Data/WOA/TEMPERATURE/'];
sal_path = [pwd '/Data/WOA/SALINITY/'];
oxy_path = [pwd '/Data/WOA/OXYGEN/'];
nit_path = [pwd '/Data/WOA/NITRATE/'];
for m = 1:12
    WOA.TMP(:,:,:,m) = ncread([temp_path 'woa18_decav_t' sprintf('%02d',m) '_01.nc'],'t_an');
    WOA.SAL(:,:,:,m) = ncread([sal_path 'woa18_decav_s' sprintf('%02d',m) '_01.nc'],'s_an');
    WOA.OXY(:,:,:,m) = ncread([oxy_path 'woa18_all_o',sprintf('%02d',m),'_01.nc'],'o_an');
    WOA.NIT(:,:,:,m) = ncread([nit_path 'woa18_all_n',sprintf('%02d',m),'_01.nc'],'n_an');
end
% Import and format WOA latitude and longitude
WOA.LAT   = ncread([nit_path 'woa18_all_n01_01.nc'],'lat');
WOA.LON   = ncread([nit_path 'woa18_all_n01_01.nc'],'lon');
WOA.DEPTH = ncread([nit_path 'woa18_all_n01_01.nc'],'depth');
[~,lat_3d,depth_3d] = ndgrid(WOA.LON,WOA.LAT,WOA.DEPTH);
WOA.PRES = -gsw_p_from_z(depth_3d,lat_3d);
% clean up
clear m temp_path sal_path nit_path

%% compare float data to WOA
if include_float

    % index to valid data points
    idx = find(~isnan(float_data.OXY) & ~isnan(float_data.NIT) & ...
        ~isnan(float_data.LAT) & ~isnan(float_data.LON) & ...
        ~isnan(float_data.PRES) & ~isnan(float_data.TIME));
    % obtain month from date
    date = datevec(datenum(float_data.YEAR,0,float_data.DAY));
    float_data.MONTH = date(:,2);
    clear date
    % pre-allocate matches
    WOA_NIT_match = nan(size(float_data.NIT));
    WOA_OXY_match = nan(size(float_data.OXY));
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
        WOA_OXY_match(idx(i)) = WOA.OXY(idx_lon(1),idx_lat(1),idx_pres(1),idx_mnth(1));
        WOA_NIT_match(idx(i)) = WOA.NIT(idx_lon(1),idx_lat(1),idx_pres(1),idx_mnth(1));
    end
    % calculate differences for oxygen
    WOA_oxy_delta = float_data.OXY-WOA_OXY_match;
    WOA_oxy_delta_per = WOA_oxy_delta./float_data.OXY;
    WOA_oxy_delta_per(isinf(WOA_oxy_delta_per))=NaN;
    st_dev_oxy_delta = std(WOA_oxy_delta,[],'omitnan');
    mean_oxy_delta = mean(WOA_oxy_delta,'omitnan');
    % calculate differences for nitrate
    WOA_nit_delta = float_data.NIT-WOA_NIT_match;
    WOA_nit_delta_per = WOA_nit_delta./float_data.NIT;
    WOA_nit_delta_per(isinf(WOA_nit_delta_per))=NaN;
    st_dev_nit_delta = std(WOA_nit_delta,[],'omitnan');
    mean_nit_delta = mean(WOA_nit_delta,'omitnan');
    % index data above 800dbar
    idx_woa_pres = float_data.PRES <= 800;

    %% plot differences between float data and WOA
    % histogram of differences
    figure; hold on
    histogram(WOA_nit_delta(idx_woa_pres));
    set(gca,'fontsize',16);
    ax = gca;
    plot([mean_nit_delta+3*st_dev_nit_delta mean_nit_delta+3*st_dev_nit_delta],ax.YLim,'r','linewidth',2);
    plot([mean_nit_delta-3*st_dev_nit_delta mean_nit_delta-3*st_dev_nit_delta],ax.YLim,'r','linewidth',2);
    xlabel('Float [NO_{3}] - WOA [NO_{3}]');
    if ~exist([pwd '/NO3/Figures/Data'],'dir'); mkdir('NO3/Figures/Data'); end
    exportgraphics(gcf,[pwd '/NO3/Figures/Data/WOA_comp_histogram_' file_date float_file_ext '.png']);
    close
    % scatter of delta values
    figure; hold on
    set(gca,'fontsize',16);
    xlabel('WOA [NO_{3}]');
    ylabel('Float [NO_{3}] - WOA [NO_{3}]');
    [counts,bin_centers] = ...
        hist3([WOA_NIT_match(idx_woa_pres),WOA_nit_delta(idx_woa_pres)],'Edges',{0:1:46 -34.5:1.5:34.5});
    h=pcolor(bin_centers{1}-mean(diff(bin_centers{1}))/2,...
        bin_centers{2}-mean(diff(bin_centers{2}))/2,counts');
    plot([0 45],[0 0],'k--');
    set(h,'EdgeColor','none');
    xlim([-0.5 46]); ylim([-34.5 34]);
    myColorMap = flipud(hot(256.*32));
    myColorMap(1,:) = 1;
    colormap(myColorMap);
    set(gca,'ColorScale','log')
    caxis([1e0 1e5]);
    c=colorbar;
    c.Label.String = 'log10(Bin Counts)';
    if ~exist([pwd '/NO3/Figures/Data'],'dir'); mkdir('NO3/Figures/Data'); end
    exportgraphics(gcf,[pwd '/NO3/Figures/Data/WOA_comp_scatter_' file_date float_file_ext '.png']);
    close
    % clean up
    clear counts bin_centers c h myColorMap

    %% remove profiles with data points that average more than 3 sigmas from WOA values
    % calculate standard deviation at each depth level
    pres_levels = unique(float_data.PRES);
    st_dev_delta_pres = nan(size(pres_levels));
    mean_delta_pres = nan(size(pres_levels));
    for z = 1:length(pres_levels)
        pres_idx = float_data.PRES == pres_levels(z);
        st_dev_delta_pres(z) = std(WOA_nit_delta(pres_idx),[],'omitnan');
        mean_delta_pres(z) = mean(WOA_nit_delta(pres_idx),'omitnan');
    end
    % loop through each profile and remove if 
    profiles = unique(float_data.PROF_ID);
    rem_idx = [];
    for p = 1:length(profiles)
        try
        idx = float_data.PROF_ID == profiles(p);
        pres_temp = float_data.PRES(idx & idx_woa_pres);
        nit_delta_temp = WOA_nit_delta(idx & idx_woa_pres);
        % average multiple profiles if necessary
        if length(pres_temp) > length(unique(pres_temp))
            for d = 1:length(unique(pres_temp))
                pres_temp_temp(d,1) = mean(pres_temp(pres_temp==pres_levels(d)),'omitnan');
                nit_delta_temp_temp(d,1) = mean(nit_delta_temp(pres_temp==pres_levels(d)),'omitnan');
            end
            idx_temp_temp = ~isnan(pres_temp_temp) & ~isnan(nit_delta_temp_temp);
            pres_temp = pres_temp_temp(idx_temp_temp);
            nit_delta_temp = nit_delta_temp_temp(idx_temp_temp);
            clear pres_temp_temp nit_delta_temp_temp idx_temp_temp
        end
        nit_stdev_temp = st_dev_delta_pres(ismember(pres_levels,pres_temp));
        if mean(nit_delta_temp./nit_stdev_temp) > 3
            disp(['Average difference for profile ' num2str(profiles(p)) ' is ' ...
                num2str(mean(nit_delta_temp./nit_stdev_temp)) ' sigma']);
            rem_idx = [rem_idx;find(idx)];
        end
        catch
            keyboard
        end
    end

end

%% compare WOD data to WOA
if include_osd || include_ctd

    % index to valid data points
    idx = find(~isnan(wod_data.OXY) & ~isnan(wod_data.NIT) & ...
        ~isnan(wod_data.LAT) & ~isnan(wod_data.LON) & ...
        ~isnan(wod_data.PRES) & ~isnan(wod_data.TIME));
    % obtain month from date
    date = datevec(datenum(wod_data.YEAR,0,wod_data.DAY));
    wod_data.MONTH = date(:,2);
    clear date
    % pre-allocate matches
    WOA_NIT_match = nan(size(wod_data.NIT));
    WOA_OXY_match = nan(size(wod_data.OXY));
    % match WOD data with WOA
    for i=1:length(idx)
        idx_lon = find(min(abs(WOA.LON - wod_data.LON(idx(i)))) == ...
            abs(WOA.LON - wod_data.LON(idx(i))));
        idx_lat = find(min(abs(WOA.LAT - wod_data.LAT(idx(i)))) == ...
            abs(WOA.LAT - wod_data.LAT(idx(i))));
        idx_pres = find(min(abs(squeeze(WOA.PRES(idx_lon,idx_lat,:)) - wod_data.PRES(idx(i)))) == ...
            abs(squeeze(WOA.PRES(idx_lon,idx_lat,:)) - wod_data.PRES(idx(i))));
        idx_mnth = find(min(abs((1:12)' - wod_data.MONTH(idx(i)))) == ...
            abs((1:12)' - wod_data.MONTH(idx(i))));
        WOA_OXY_match(idx(i)) = WOA.OXY(idx_lon(1),idx_lat(1),idx_pres(1),idx_mnth(1));
        WOA_NIT_match(idx(i)) = WOA.NIT(idx_lon(1),idx_lat(1),idx_pres(1),idx_mnth(1));
    end
    % calculate differences for oxygen
    WOA_oxy_delta = wod_data.OXY-WOA_OXY_match;
    WOA_oxy_delta_per = WOA_oxy_delta./wod_data.OXY;
    WOA_oxy_delta_per(isinf(WOA_oxy_delta_per))=NaN;
    st_dev_oxy_delta = std(WOA_oxy_delta,[],'omitnan');
    mean_oxy_delta = mean(WOA_oxy_delta,'omitnan');
    % calculate differences for nitrate
    WOA_nit_delta = wod_data.NIT-WOA_NIT_match;
    WOA_nit_delta_per = WOA_nit_delta./wod_data.NIT;
    WOA_nit_delta_per(isinf(WOA_nit_delta_per))=NaN;
    st_dev_nit_delta = std(WOA_nit_delta,[],'omitnan');
    mean_nit_delta = mean(WOA_nit_delta,'omitnan');
    % index data above 800dbar
    idx_woa_pres = wod_data.PRES <= 800;
    % clean up
    clear idx idx_lon idx_lat idx_pres idx_mnth WOA i

    %% plot differences between WOD data and WOA
    % histogram of differences
    figure; hold on
    histogram(WOA_nit_delta(idx_woa_pres));
    set(gca,'fontsize',16);
    ax = gca;
    plot([mean_nit_delta+3*st_dev_nit_delta mean_nit_delta+3*st_dev_nit_delta],ax.YLim,'r','linewidth',2);
    plot([mean_nit_delta-3*st_dev_nit_delta mean_nit_delta-3*st_dev_nit_delta],ax.YLim,'r','linewidth',2);
    xlabel('Float [NO_{3}] - WOA [NO_{3}]');
    if ~exist([pwd '/NO3/Figures/Data'],'dir'); mkdir('NO3/Figures/Data'); end
    exportgraphics(gcf,[pwd '/NO3/Figures/Data/WOA_comp_histogram_' file_date float_file_ext '.png']);
    close
    % scatter of delta values
    figure; hold on
    set(gca,'fontsize',16);
    xlabel('WOA [NO_{3}]');
    ylabel('Float [NO_{3}] - WOA [NO_{3}]');
    [counts,bin_centers] = ...
        hist3([WOA_NIT_match(idx_woa_pres),WOA_nit_delta(idx_woa_pres)],'Edges',{0:1:46 -34.5:1.5:34.5});
    h=pcolor(bin_centers{1}-mean(diff(bin_centers{1}))/2,...
        bin_centers{2}-mean(diff(bin_centers{2}))/2,counts');
    plot([0 45],[0 0],'k--');
    set(h,'EdgeColor','none');
    xlim([-0.5 46]); ylim([-34.5 34]);
    myColorMap = flipud(hot(256.*32));
    myColorMap(1,:) = 1;
    colormap(myColorMap);
    set(gca,'ColorScale','log')
    caxis([1e0 1e5]);
    c=colorbar;
    c.Label.String = 'log10(Bin Counts)';
    if ~exist([pwd '/NO3/Figures/Data'],'dir'); mkdir('NO3/Figures/Data'); end
    exportgraphics(gcf,[pwd '/NO3/Figures/Data/WOA_comp_scatter_' file_date float_file_ext '.png']);
    close
    % clean up
    clear counts bin_centers c h myColorMap

    %% remove profiles with data points that average more than 3 sigmas from WOA values
    % calculate standard deviation at each depth level
    pres_levels = unique(wod_data.PRES);
    st_dev_delta_pres = nan(size(pres_levels));
    mean_delta_pres = nan(size(pres_levels));
    for z = 1:length(pres_levels)
        pres_idx = wod_data.PRES == pres_levels(z);
        st_dev_delta_pres(z) = std(WOA_nit_delta(pres_idx),[],'omitnan');
        mean_delta_pres(z) = mean(WOA_nit_delta(pres_idx),'omitnan');
    end
    % loop through each profile and remove if 
    profiles = unique(wod_data.ID);
    rem_idx = [];
    for p = 1:length(profiles)
        try
        idx = wod_data.ID == profiles(p);
        pres_temp = wod_data.PRES(idx & idx_woa_pres);
        nit_delta_temp = WOA_nit_delta(idx & idx_woa_pres);
        % average multiple profiles if necessary
        if length(pres_temp) > length(unique(pres_temp))
            for d = 1:length(unique(pres_temp))
                pres_temp_temp(d,1) = mean(pres_temp(pres_temp==pres_levels(d)),'omitnan');
                nit_delta_temp_temp(d,1) = mean(nit_delta_temp(pres_temp==pres_levels(d)),'omitnan');
            end
            idx_temp_temp = ~isnan(pres_temp_temp) & ~isnan(nit_delta_temp_temp);
            pres_temp = pres_temp_temp(idx_temp_temp);
            nit_delta_temp = nit_delta_temp_temp(idx_temp_temp);
            clear pres_temp_temp nit_delta_temp_temp idx_temp_temp
        end
        nit_stdev_temp = st_dev_delta_pres(ismember(pres_levels,pres_temp));
        if mean(nit_delta_temp./nit_stdev_temp) > 3
            disp(['Average difference for profile ' num2str(profiles(p)) ' is ' ...
                num2str(mean(nit_delta_temp./nit_stdev_temp)) ' sigma']);
            rem_idx = [rem_idx;find(idx)];
        end
        catch
            keyboard
        end
    end

end

% clean up
clear idx idx_lon idx_lat idx_pres idx_mnth WOA i

%% Co-locate and compare float and ship data (via profile crossovers)
pres_levels = [2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975]';
float_profile_IDs = unique(float_data.PROF_ID);
% pre-allocate
crossover.id = [];
crossover.lon = [];
crossover.lat = [];
crossover.pres = [];
crossover.sigma_float = [];
crossover.nit_float = [];
crossover.nit_float_gradient = [];
crossover.temp_float = [];
crossover.sal_float = [];
crossover.nit_ship = [];
crossover.nit_ship_gradient = [];
crossover.sigma_ship = [];
crossover.nit_delta = [];
crossover.nit_delta_per = [];
% include WOD data or not
if include_ctd || include_osd
    ship_lon = [glodap_data.LON;wod_data.LON];
    ship_lat = [glodap_data.LAT;wod_data.LAT];
    ship_time = [glodap_data.TIME;wod_data.TIME];
    ship_sigma = [glodap_data.SIGMA;wod_data.SIGMA];
    ship_pres = [glodap_data.PRES;wod_data.PRES];
    ship_nit = [glodap_data.NIT;wod_data.NIT];
else
    ship_lon = [glodap_data.LON];
    ship_lat = [glodap_data.LAT];
    ship_time = [glodap_data.TIME];
    ship_sigma = [glodap_data.SIGMA];
    ship_pres = [glodap_data.PRES];
    ship_nit = [glodap_data.NIT];
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
        nit_temp = float_data.NIT(float_data.PROF_ID==float_profile_IDs(p));
        pres_temp = float_data.PRES(float_data.PROF_ID==float_profile_IDs(p));
        sigma_temp = float_data.SIGMA(float_data.PROF_ID==float_profile_IDs(p));
        lon_temp = float_data.LON(float_data.PROF_ID==float_profile_IDs(p));
        lat_temp = float_data.LAT(float_data.PROF_ID==float_profile_IDs(p));
        temp_temp = float_data.TEMP(float_data.PROF_ID==float_profile_IDs(p));
        sal_temp = float_data.SAL(float_data.PROF_ID==float_profile_IDs(p));
        % combine float profiles (helpful if there are more than one or
        % multiple measurements at the same depths)
        nit_float = nan(length(pres_levels),1);
        pres_float = nan(length(pres_levels),1);
        sigma_float = nan(length(pres_levels),1);
        lon_float = nan(length(pres_levels),1);
        lat_float = nan(length(pres_levels),1);
        temp_float = nan(length(pres_levels),1);
        sal_float = nan(length(pres_levels),1);
        for z = 1:length(pres_levels)
            idx = pres_temp == pres_levels(z);
            nit_float(z) = mean(nit_temp(idx),'omitnan');
            pres_float(z) = mean(pres_temp(idx),'omitnan');
            sigma_float(z) = mean(sigma_temp(idx),'omitnan');
            lon_float(z) = mean(lon_temp(idx),'omitnan');
            lat_float(z) = mean(lat_temp(idx),'omitnan');
            temp_float(z) = mean(temp_temp(idx),'omitnan');
            sal_float(z) = mean(sal_temp(idx),'omitnan');
        end
        % get nitgen gradient
        nit_float_gradient = [0;diff(nit_float)];
        % extract ship data
        pres_temp = ship_pres(idx_all);
        nit_temp = ship_nit(idx_all);
        sigma_temp = ship_sigma(idx_all);
        % combine ship profiles (helpful if there are more than one or
        % multiple measurements at the same depths)
        nit_ship = nan(length(pres_levels),1);
        pres_ship = nan(length(pres_levels),1);
        sigma_ship = nan(length(pres_levels),1);
        for z = 1:length(pres_levels)
            idx = pres_temp == pres_levels(z);
            pres_ship(z) = mean(pres_temp(idx),'omitnan');
            nit_ship(z) = mean(nit_temp(idx),'omitnan');
            sigma_ship(z) = mean(sigma_temp(idx),'omitnan');
        end
        % nit_ship = nit_temp(~isnan(nit_temp));
        % pres_ship = pres_temp(~isnan(pres_temp));
        % sigma_ship = sigma_temp(~isnan(sigma_temp));
        % get nitrogen gradient
        nit_ship_gradient = [0;diff(nit_ship)];
        % plot matched profiles
        % figure; plot(nit_float,pres_levels(idx_pres),nit_ship,pres_levels(idx_pres)); set(gca,'YDir','reverse'); close;
        % log matched data
        crossover.id = [crossover.id;repmat(float_profile_IDs(p),length(pres_float),1)];
        crossover.lon = [crossover.lon;lon_float];
        crossover.lat = [crossover.lat;lat_float];
        crossover.pres = [crossover.pres;pres_float];
        crossover.nit_float = [crossover.nit_float;nit_float];
        crossover.nit_float_gradient = [crossover.nit_float_gradient;nit_float_gradient];
        crossover.temp_float = [crossover.temp_float;temp_float];
        crossover.sal_float = [crossover.sal_float;sal_float];
        crossover.nit_ship = [crossover.nit_ship;nit_ship];
        crossover.nit_ship_gradient = [crossover.nit_ship_gradient;nit_ship_gradient];
        crossover.sigma_float = [crossover.sigma_float;sigma_float];
        crossover.sigma_ship = [crossover.sigma_ship;sigma_ship];
        try
        crossover.nit_delta = [crossover.nit_delta;nit_float-nit_ship];
        catch
            keyboard
        end
        crossover.nit_delta_per = [crossover.nit_delta_per;...
            (nit_float-nit_ship)./nit_float];
    end
end

% index to below 300 dbars
idx = crossover.pres > 300 & ~isnan(crossover.nit_float) & ~isnan(crossover.nit_ship);
% clean up
clear pres_levels float_profile_IDs lon lat time nit_float pres_float
clear idx_lon idx_lat idx_time nit_temp pres_temp idx_pres nit_ship

% save crossover data
if ~exist([pwd '/NO3/Data'],'dir'); mkdir('Data'); end
save(['NO3/Data/crossover_data_' file_date float_file_ext],'crossover','idx','-v7.3')

%% depth-dependent correction (using two depths)
pressures = unique(crossover.pres);
avg_delta_by_p = nan(size(pressures));
avg_by_p = nan(size(pressures));
for p = 1:length(pressures)
    idx_p = crossover.pres == pressures(p);
    avg_delta_by_p(p) = mean(crossover.nit_delta(idx_p),'omitnan');
    avg_by_p(p) = mean(crossover.nit_ship(idx_p),'omitnan');
end
[~,idx_min_delta] = min(abs(avg_delta_by_p)); % find maximum crossover offset
min_delta_pressure = pressures(idx_min_delta);
disp(['Minimum crossover offset at ' num2str(min_delta_pressure) ' dbar.'])
% deep adjustment equation
mdl = fitlm(pressures(idx_min_delta:end),...
    avg_delta_by_p(idx_min_delta:end),'Intercept',true);
slp_deep = mdl.Coefficients.Estimate(2);
int_deep = mdl.Coefficients.Estimate(1);
% shallow adjustment equation
mdl = fitlm(pressures(1:idx_min_delta-1),...
    avg_delta_by_p(1:idx_min_delta-1),'Intercept',true);
slp_shal = mdl.Coefficients.Estimate(2);
int_shal = mdl.Coefficients.Estimate(1);
% apply linear correction to float NO3 (with pressure)
corr_fac = nan(size(crossover.nit_float));
idx_deep = crossover.pres >= pressures(idx_min_delta);
corr_fac(idx_deep) = slp_deep .* crossover.pres(idx_deep) + int_deep;
idx_shal = crossover.pres < pressures(idx_min_delta);
corr_fac(idx_shal) = slp_shal .* crossover.pres(idx_shal) + int_shal;
crossover.nit_float_corr = crossover.nit_float - corr_fac;

%% re-calculate deltas
crossover.nit_delta_corr = crossover.nit_float_corr-crossover.nit_ship;
crossover.nit_delta_per_corr = (crossover.nit_float_corr-crossover.nit_ship)./crossover.nit_float_corr;

%% save crossover data and correction factors
if ~exist([pwd '/NO3/Data'],'dir'); mkdir('O2/Data'); end
% save(['O2/Data/float_corr_' file_date float_file_ext],'slp','int');
% clear slp int
save(['NO3/Data/crossover_data_' file_date float_file_ext],'crossover','idx','-v7.3')
save(['NO3/Data/float_corr_' file_date float_file_ext],'slp_deep',...
    'int_deep','slp_shal','int_shal','min_delta_pressure');
clear slp_deep int_deep slp_shal int_shal

%% plot uncorrected float vs. glodap residuals
load(['NO3/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
scatter(crossover.nit_ship(idx),crossover.nit_float(idx),20,'.');
plot([0,45],[0 45],'k--');
ylabel('Binned BGC Argo Nitrate Data (\mumol kg^{-1})');
xlabel('Binned GLODAP Nitrate Data (\mumol kg^{-1})');
crossover.err_md = median(crossover.nit_delta(idx));
crossover.std = std(crossover.nit_delta(idx));
crossover.rmse = sqrt(mean((crossover.nit_delta(idx)).^2));
crossover.r2 = corr(crossover.nit_ship(idx),crossover.nit_float(idx));
crossover.slope = polyfit(crossover.nit_ship(idx),crossover.nit_float(idx),1);
crossover.slope = crossover.slope(1);
text(30,12,['Med. Err. = ' num2str(round(crossover.err_md,2))],'fontsize',12);
text(30,9,['RMSE = ' num2str(round(crossover.rmse,1))],'fontsize',12);
text(30,6,['R^{2} = ' num2str(round(crossover.r2,2))],'fontsize',12);
exportgraphics(gcf,['NO3/Figures/Data/delta_vs_float_300_uncorr_' file_date float_file_ext '.png']);
close
% histogram
figure; hold on
histogram(crossover.nit_delta(idx),'Normalization','probability');
ylim('manual');
xlim([-3 3]);
set(gca,'fontsize',16);
sig2_min = double(mean(crossover.nit_delta(idx))-2*std(crossover.nit_delta(idx)));
sig2_max = double(mean(crossover.nit_delta(idx))+2*std(crossover.nit_delta(idx)));
plot([sig2_min sig2_min],[0 1],'r','linewidth',1);
plot([sig2_max sig2_max],[0 1],'r','linewidth',1);
xlabel('Float [NO_{3}] - GLODAP [NO_{3}]');
xL = xlim; yL = ylim;
text(0.9*xL(1),0.9*yL(2),['Mean \Delta[NO_{3}] = ' ...
    num2str(round(mean(crossover.nit_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
text(0.9*xL(1),0.8*yL(2),['Std. \Delta[NO_{3}] = ' ...
    num2str(round(std(crossover.nit_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
exportgraphics(gcf,[pwd '/NO3/Figures/Data/GLODAP_comp_histogram_uncorr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot corrected float vs. glodap residuals
load(['NO3/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
scatter(crossover.nit_ship(idx),crossover.nit_float_corr(idx),20,'.');
plot([0,45],[0 45],'k--');
ylabel('Corrected, Binned BGC Argo Nitrate Data (\mumol kg^{-1})');
xlabel('Binned GLODAP Nitrate Data (\mumol kg^{-1})');
crossover.err_md_corr = median(crossover.nit_delta_corr(idx));
crossover.std_corr = std(crossover.nit_delta_corr(idx));
crossover.rmse_corr = sqrt(mean((crossover.nit_delta_corr(idx)).^2));
crossover.r2_corr = corr(crossover.nit_ship(idx),crossover.nit_float_corr(idx));
crossover.slope_corr = polyfit(crossover.nit_ship(idx),crossover.nit_float_corr(idx),1);
crossover.slope_corr = crossover.slope_corr(1);
text(30,12,['Med. Err. = ' num2str(round(crossover.err_md_corr,2))],'fontsize',12);
text(30,9,['RMSE = ' num2str(round(crossover.rmse_corr,1))],'fontsize',12);
text(30,6,['R^{2} = ' num2str(round(crossover.r2_corr,2))],'fontsize',12);
exportgraphics(gcf,['NO3/Figures/Data/delta_vs_float_300_corr_' file_date float_file_ext '.png']);
close
% histogram
figure; hold on
histogram(crossover.nit_delta_corr(idx),'Normalization','probability');
ylim('manual');
xlim([-3 3]);
set(gca,'fontsize',16);
sig2_min = double(mean(crossover.nit_delta_corr(idx))-2*std(crossover.nit_delta_corr(idx)));
sig2_max = double(mean(crossover.nit_delta_corr(idx))+2*std(crossover.nit_delta_corr(idx)));
plot([sig2_min sig2_min],[0 1],'r','linewidth',1);
plot([sig2_max sig2_max],[0 1],'r','linewidth',1);
xlabel('Float [NO_{3}] - GLODAP [NO_{3}]');
xL = xlim; yL = ylim;
text(0.9*xL(1),0.9*yL(2),['Mean \Delta[NO_{3}] = ' ...
    num2str(round(mean(crossover.nit_delta_corr(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
text(0.9*xL(1),0.8*yL(2),['Std. \Delta[NO_{3}] = ' ...
    num2str(round(std(crossover.nit_delta_corr(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
exportgraphics(gcf,[pwd '/NO3/Figures/Data/GLODAP_comp_histogram_corr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% adjust and save float data
load(['NO3/Data/float_corr_' file_date float_file_ext]);
vars = fieldnames(float_data);
for v = 1:length(vars)
    if strcmp(vars{v},'NIT')

        % apply linear correction to float O2 (with pressure)
        corr_fac = nan(size(float_data.NIT));
        idx_deep = float_data.PRES >= min_delta_pressure;
        corr_fac(idx_deep) = float_data.PRES(idx_deep).*slp_deep + int_deep;
        idx_shal = float_data.PRES < min_delta_pressure;
        corr_fac(idx_shal) = float_data.PRES(idx_shal).*slp_shal + int_shal;
        float_data_adjusted.NIT = float_data.NIT - corr_fac;
        disp(['Average adjustment = ' ...
            num2str(round(mean(corr_fac,'omitnan'),2)) ' umol/kg']);
    else
        float_data_adjusted.(vars{v}) = float_data.(vars{v});
    end
end

if include_float
    % save adjusted float data
    if ~exist([pwd '/NO3/Data'],'dir'); mkdir('NO3/Data'); end
    save(['NO3/Data/processed_float_no3_data_adjusted_' ...
        file_date float_file_ext '.mat'],...
        'float_data_adjusted','file_date','-v7.3');
    clear slp int float_data float_data_adjusted v vars
end

if include_osd || include_ctd
    % save adjusted WOD data
    if ~exist([pwd '/NO3/Data'],'dir'); mkdir('NO3/Data'); end
    save(['NO3/Data/processed_wod_no3_data_adjusted_' ...
        num2str(end_year) '.mat'],...
        'wod_data','file_date','-v7.3');
end

if include_glodap
    % save adjusted GLODAP data
    if ~exist([pwd '/NO3/Data'],'dir'); mkdir('NO3/Data'); end
    save(['NO3/Data/processed_glodap_no3_data_adjusted_' ...
        glodap_vrs '.mat'],...
        'glodap_data','file_date','-v7.3');
end

% 
% else
% 
% % display information
% disp('Float data already adjusted.')
% 
% end

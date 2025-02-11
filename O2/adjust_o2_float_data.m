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

function adjust_o2_float_data(float_file_ext,file_date,glodap_year)

%% load interpolated float and glodap data
load(['O2/Data/processed_float_oxygen_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date');
load(['O2/Data/processed_glodap_oxygen_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');
% load(['O2/Data/processed_wod_ctd_oxygen_data_' num2str(glodap_year) '.mat'],...
%     'wod_data');

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
    idx_pres = find(min(abs(WOA.DEPTH - float_data.PRES(idx(i)))) == ...
        abs(WOA.DEPTH - float_data.PRES(idx(i))));
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
plot([mean_delta+5*st_dev_delta mean_delta+5*st_dev_delta],ax.YLim,'r','linewidth',2);
plot([mean_delta-5*st_dev_delta mean_delta-5*st_dev_delta],ax.YLim,'r','linewidth',2);
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
caxis([1e0 1e5]);
c=colorbar;
c.Label.String = 'log10(Bin Counts)';
exportgraphics(gcf,[pwd '/O2/Figures/Data/WOA_comp_scatter_' file_date float_file_ext '.png']);
close
% clean up
clear counts bin_centers c h myColorMap

%% remove data points more than 5 sigmas from WOA value
idx_rem = WOA_delta > mean_delta+5.*st_dev_delta | WOA_delta < mean_delta-5.*st_dev_delta;
% sum(idx_rem)
vars = fieldnames(float_data);
for v = 1:length(vars)
    float_data.(vars{v})(idx_rem) = [];
end
% clean up
clear WOA_match WOA_delta WOA_delta_per mean_delta st_dev_delta idx_rem v vars

%% determine histogram counts and indices
% % establish edges of bins
% x_edges = -180:180;
% x_bins = -179.5:179.5;
% y_edges = -85:85;
% y_bins = -84.5:84.5;
% z_edges = [0 5:10:175 190:20:450 475:50:1375 1450:100:1950 2000];
% z_bins = [2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975];
% mn_edges = 1:12;
% yr_edges = 2004:glodap_year;
% t_edges = datenum([[repelem(yr_edges,1,length(mn_edges)) yr_edges(end)+1]', ...
%                   [repmat(mn_edges,1,length(yr_edges)) 1]', ...
%                   [zeros(1,length(mn_edges)*length(yr_edges)+1)]']);
% t_bins = datenum([[repelem(yr_edges,1,length(mn_edges))]', ...
%                   [repmat(mn_edges,1,length(yr_edges))]', ...
%                   [repmat(15,1,length(mn_edges)*length(yr_edges))]']);
% % get histogram counts in each bin
% [~,~,Xnum_float] = histcounts(float_data.LON,x_edges);
% [~,~,Ynum_float] = histcounts(float_data.LAT,y_edges);
% [~,~,Znum_float] = histcounts(float_data.PRES,z_edges);
% [~,~,Tnum_float] = histcounts(float_data.TIME,t_edges);
% [~,~,Xnum_glodap] = histcounts(glodap_data.LON,x_edges);
% [~,~,Ynum_glodap] = histcounts(glodap_data.LAT,y_edges);
% [~,~,Znum_glodap] = histcounts(glodap_data.PRES,z_edges);
% [~,~,Tnum_glodap] = histcounts(glodap_data.TIME,t_edges);
% % accumulate index of counts
% subs_float = [Xnum_float,Ynum_float,Znum_float,Tnum_float];
% idx_float = ~any(subs_float==0,2);
% subs_glodap = [Xnum_glodap,Ynum_glodap,Znum_glodap,Tnum_glodap];
% idx_glodap = ~any(subs_glodap==0,2);
% clear Xnum_float Ynum_float Znum_float Tnum_float
% clear Xnum_glodap Ynum_glodap Znum_glodap Tnum_glodap
% % determine size of 4D grid
% sz = [length(x_bins),length(y_bins),length(z_bins),length(t_bins)];
% % clean up
% clear x_edges y_edges z_edges mn_edges yr_edges t_edges
% 
% %% Bin float and glodap data
% binned_data.oxy_float = single(nan(sz));
% binned_data.oxy_glodap = single(nan(sz));
% for m = 1:length(t_bins)
%     % month-specific float index
%     idx_float_tmp = idx_float;
%     idx_float_tmp(subs_float(:,4)~=m) = false;
%     % month-specific glodap index
%     idx_glodap_tmp = idx_glodap;
%     idx_glodap_tmp(subs_glodap(:,4)~=m) = false;
%     % bin oxygen data
%     binned_data.oxy_float(:,:,:,m) = single(accumarray(subs_float(idx_float_tmp,1:3),float_data.OXY(idx_float_tmp),sz(1:3),@nanmean,nan));
%     binned_data.oxy_glodap(:,:,:,m) = single(accumarray(subs_glodap(idx_glodap_tmp,1:3),glodap_data.OXY(idx_glodap_tmp),sz(1:3),@nanmean,nan));
% end
% % add pressure bins
% binned_data.pres = repmat(permute(single(z_bins),[3 1 2]),length(x_bins),length(y_bins),1,length(t_bins));
% % save binned data
% if ~exist([pwd '/O2/Data'],'dir'); mkdir('O2/Data'); end
% save(['O2/Data/binned_data_' file_date float_file_ext],'binned_data','-v7.3')
% 
% % clean up
% clear binned_data idx_float idx_glodap
% clear idx_float_tmp idx_glodap_tmp m sz subs_float subs_glodap
% clear x_edges x_bins y_bins z_edges z_bins t_bins

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
crossover.oxy_glodap = [];
crossover.oxy_glodap_gradient = [];
crossover.sigma_glodap = [];
crossover.oxy_delta = [];
crossover.oxy_delta_per = [];
% cycle through float profiles 
for p = 1:length(float_profile_IDs)
    % float position
    lon = mean(float_data.LON(float_data.PROF_ID==float_profile_IDs(p)));
    lat = mean(float_data.LAT(float_data.PROF_ID==float_profile_IDs(p)));
    time = mean(float_data.TIME(float_data.PROF_ID==float_profile_IDs(p)));
    % match index
    idx_lon = abs(glodap_data.LON - lon) < 0.5;
    idx_lat = abs(glodap_data.LAT - lat) < 0.5;
    idx_time = abs(glodap_data.TIME - time) < 30;
    idx_all = idx_lon & idx_lat & idx_time;
    % if there is a matching glodap profile
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
        % extract glodap data
        oxy_glodap_temp = glodap_data.OXY(idx_all);
        sigma_glodap_temp = glodap_data.SIGMA(idx_all);
        % take average of all matching glodap profiles (helps if there are multiple)
        oxy_glodap = nan(58,1);
        sigma_glodap = nan(58,1);
        for gld_z = 1:58
            oxy_glodap(gld_z) = mean(oxy_glodap_temp(gld_z:58:end),'omitnan');
            sigma_glodap(gld_z) = mean(sigma_glodap_temp(gld_z:58:end),'omitnan');
        end
        % get oxygen gradient
        oxy_glodap_gradient = [0;diff(oxy_glodap)];
        % log matched data
        crossover.id = [crossover.id;repmat(float_profile_IDs(p),length(pres_float),1)];
        crossover.lon = [crossover.lon;lon_float];
        crossover.lat = [crossover.lat;lat_float];
        crossover.pres = [crossover.pres;pres_float];
        crossover.oxy_float = [crossover.oxy_float;oxy_float];
        crossover.oxy_float_gradient = [crossover.oxy_float_gradient;oxy_float_gradient];
        crossover.temp_float = [crossover.temp_float;temp_float];
        crossover.sal_float = [crossover.sal_float;sal_float];
        crossover.oxy_glodap = [crossover.oxy_glodap;oxy_glodap(idx_pres)];
        crossover.oxy_glodap_gradient = [crossover.oxy_glodap_gradient;oxy_glodap_gradient(idx_pres)];
        crossover.sigma_float = [crossover.sigma_float;sigma_float];
        crossover.sigma_glodap = [crossover.sigma_glodap;sigma_glodap];
        crossover.oxy_delta = [crossover.oxy_delta;oxy_float-oxy_glodap(idx_pres)];
        crossover.oxy_delta_per = [crossover.oxy_delta_per;...
            (oxy_float-oxy_glodap(idx_pres))./oxy_float];
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
idx = crossover.pres > 300 & ~isnan(crossover.oxy_float) & ~isnan(crossover.oxy_glodap);
% clean up
clear pres_levels float_profile_IDs lon lat time oxy_float pres_float
clear idx_lon idx_lat idx_time oxy_temp pres_temp idx_pres oxy_glodap
% fit delta against [O2]
mdl = fitlm(crossover.oxy_sat_float(idx),crossover.oxy_delta_sat_per(idx),'Intercept',true);
slp = mdl.Coefficients.Estimate(2);
int = mdl.Coefficients.Estimate(1);
% apply linear correction to float O2
corr_fac = slp .* crossover.oxy_sat_float + int; % modelled saturation delta
crossover.oxy_float_sat_per_corr = crossover.oxy_float_sat_per - corr_fac;
crossover.oxy_float_corr = (crossover.oxy_float_sat_per_corr./100).*crossover.oxy_sat_float;
% re-calculate delta
crossover.oxy_delta_corr = crossover.oxy_float_corr-crossover.oxy_glodap;
crossover.oxy_delta_per_corr = (crossover.oxy_float_corr-crossover.oxy_glodap)./crossover.oxy_float_corr;
% save crossover data
if ~exist([pwd '/O2/Data'],'dir'); mkdir('Data'); end
save(['O2/Data/crossover_data_' file_date float_file_ext],'crossover','idx','-v7.3')
% save correction factors
if ~exist([pwd '/O2/Data'],'dir'); mkdir('O2/Data'); end
save(['O2/Data/float_corr_' file_date float_file_ext],'slp','int');
clear slp int

%% Co-locate and compare float and glodap data (via bins)
% % load binned data
% load(['O2/Data/binned_data_' file_date float_file_ext],'binned_data')
% % index where binned float and glodap data overlap under 300 dbars
% idx = ~isnan(binned_data.oxy_float) & ~isnan(binned_data.oxy_glodap) & ...
%     binned_data.pres > 300;
% % add matched data to structure
% matched_data.oxy_float = binned_data.oxy_float(idx);
% matched_data.oxy_glodap = binned_data.oxy_glodap(idx);
% matched_data.oxy_delta = matched_data.oxy_glodap-matched_data.oxy_float;
% matched_data.oxy_delta_per = (matched_data.oxy_glodap-matched_data.oxy_float)./matched_data.oxy_glodap;
% matched_data.pres = binned_data.pres(idx);
% % fit delta against [O2]
% mdl = fitlm(matched_data.oxy_float,matched_data.oxy_delta_per,'Intercept',true);
% slp = mdl.Coefficients.Estimate(2);
% int = mdl.Coefficients.Estimate(1);
% % apply linear correction to float O2
% matched_data.oxy_float_corr = matched_data.oxy_float + ...
%     (slp.*matched_data.oxy_float + int).*matched_data.oxy_float;
% % re-calculate delta
% matched_data.oxy_delta_corr = matched_data.oxy_glodap-matched_data.oxy_float_corr;
% matched_data.oxy_delta_per_corr = matched_data.oxy_glodap-matched_data.oxy_float_corr;
% % clean up
% clear idx mdl binned_data
% 
% % save matched float and glodap data
% if ~exist([pwd '/O2/Data'],'dir'); mkdir('O2/Data'); end
% save(['O2/Data/matched_data_' file_date float_file_ext],'matched_data','-v7.3');
% clear matched_data
% % save correction factors
% if ~exist([pwd '/O2/Data'],'dir'); mkdir('O2/Data'); end
% save(['O2/Data/float_corr_' file_date float_file_ext],'slp','int');
% clear slp int

%% plot uncorrected float vs. glodap residuals
load(['O2/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
scatter(crossover.oxy_glodap(idx),crossover.oxy_float(idx),20,'.');
plot([0,450],[0 450],'k--');
ylabel('Binned BGC Argo Oxygen Data (\mumol kg^{-1})');
xlabel('Binned GLODAP Oxygen Data (\mumol kg^{-1})');
crossover.err_md = median(crossover.oxy_delta(idx));
crossover.std = std(crossover.oxy_delta(idx));
crossover.rmse = sqrt(mean((crossover.oxy_delta(idx)).^2));
crossover.r2 = corr(crossover.oxy_glodap(idx),crossover.oxy_float(idx));
crossover.slope = polyfit(crossover.oxy_glodap(idx),crossover.oxy_float(idx),1);
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
xlabel('GLODAP [O_{2}] - Float [O_{2}]');
xL = xlim; yL = ylim;
text(0.9*xL(1),0.9*yL(2),['Mean \Delta[O_{2}] = ' ...
    num2str(round(mean(crossover.oxy_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
text(0.9*xL(1),0.8*yL(2),['Std. \Delta[O_{2}] = ' ...
    num2str(round(std(crossover.oxy_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
exportgraphics(gcf,[pwd '/O2/Figures/Data/GLODAP_comp_histogram_uncorr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot corrected float vs. glodap residuals
load(['O2/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
scatter(crossover.oxy_glodap(idx),crossover.oxy_float_corr(idx),20,'.');
plot([0,450],[0 450],'k--');
ylabel('Corrected, Binned BGC Argo Oxygen Data (\mumol kg^{-1})');
xlabel('Binned GLODAP Oxygen Data (\mumol kg^{-1})');
crossover.err_md_corr = median(crossover.oxy_delta_corr(idx));
crossover.std_corr = std(crossover.oxy_delta_corr(idx));
crossover.rmse_corr = sqrt(mean((crossover.oxy_delta_corr(idx)).^2));
crossover.r2_corr = corr(crossover.oxy_glodap(idx),crossover.oxy_float_corr(idx));
crossover.slope_corr = polyfit(crossover.oxy_glodap(idx),crossover.oxy_float_corr(idx),1);
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
xlabel('GLODAP [O_{2}] - Float [O_{2}]');
xL = xlim; yL = ylim;
text(0.9*xL(1),0.9*yL(2),['Mean \Delta[O_{2}] = ' ...
    num2str(round(mean(crossover.oxy_delta_corr(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
text(0.9*xL(1),0.8*yL(2),['Std. \Delta[O_{2}] = ' ...
    num2str(round(std(crossover.oxy_delta_corr(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
exportgraphics(gcf,[pwd '/O2/Figures/Data/GLODAP_comp_histogram_corr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot uncorrected float vs. glodap residuals
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
xlabel('Binned GLODAP Oxygen Data (\mumol kg^{-1})');
crossover.err_md = median(crossover.oxy_delta(idx));
crossover.std = std(crossover.oxy_delta(idx));
crossover.rmse = sqrt(mean((crossover.oxy_delta(idx)).^2));
crossover.r2 = corr(crossover.oxy_glodap(idx),crossover.oxy_float(idx));
crossover.slope = polyfit(crossover.oxy_glodap(idx),crossover.oxy_float(idx),1);
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
xlabel('GLODAP [O_{2}] - Float [O_{2}]');
xL = xlim; yL = ylim;
text(0.9*xL(1),0.9*yL(2),['Mean \Delta[O_{2}] = ' ...
    num2str(round(mean(crossover.oxy_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
text(0.9*xL(1),0.8*yL(2),['Std. \Delta[O_{2}] = ' ...
    num2str(round(std(crossover.oxy_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
exportgraphics(gcf,[pwd '/O2/Figures/Data/GLODAP_comp_histogram_uncorr_' file_date float_file_ext '.png']);
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

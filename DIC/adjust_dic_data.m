% adjust_dic_float_data
%
% DESCRIPTION:
% This function bins float and glodap data, co-locates corresponding bins,
% compares the two datasets, and calculates a bulk correction factor based
% on the discrepancies between the two.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 07/23/2026

function adjust_dic_data(float_file_ext,glodap_vrs,snap_date,...
    include_float,include_glodap,include_osd,include_ctd)

%% load interpolated float and glodap data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['DIC/Data/processed_float_dic_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date');
load(['DIC/Data/processed_glodap_dic_data_' glodap_vrs '.mat'],...
    'glodap_data');

%% pH range test
idx = float_data.PH < 7.3 | float_data.PH > 8.5;
vars = fieldnames(float_data);
for v = 1:length(vars)
    float_data.(vars{v})(idx) = [];
end

%% calculate DIC from float data
% convert pressure to depth
float_data.DEPTH = -gsw_z_from_p(float_data.PRES,float_data.LAT);
% make ESPER estimates
esper_preds = ESPER_Mixed([1 2 3 4 6],[float_data.LON float_data.LAT float_data.DEPTH],...
    [float_data.TEMP float_data.SAL float_data.OXY],[2 1 6],'Equations',7,...
    'EstDates',float_data.YEAR + (float_data.DAY/365.25));
float_data.TA_est = esper_preds.TA;
float_data.DIC_est = esper_preds.DIC;
float_data.pH_est = esper_preds.pH;
float_data.PHOS_est = esper_preds.phosphate;
float_data.SIL_est = esper_preds.silicate;
idx = ~isnan(float_data.PH) & ~isnan(float_data.pH_est);

%% plot measured vs. estimated pH
figure('visible','off','position',[100 100 1000 400]);
set(gca,'fontsize',15,'Box','on','LineWidth',2); hold on;
histogram2(gca,float_data.PH,float_data.pH_est,'DisplayStyle', 'tile');
plot([7.4 8.3],[7.4 8.3],'b');
xlim([7.4 8.3]); ylim([7.4 8.3]);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
clim([1e0 1e4]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text(7.5,8.2,['Average Diff. = ' ...
    num2str(round(mean(float_data.PH(idx)-float_data.pH_est(idx)),4)) ...
    ' +/- ' num2str(round(std(float_data.PH(idx)-float_data.pH_est(idx)),4))],...
    'FontSize',12);
xlabel('Float measured pH'); ylabel('Estimated pH (ESPER w/ T/S/O2)');
export_fig(gcf,['DIC/Figures/Data/pH_meas_vs_esper_uncorr_' ...
    file_date float_file_ext '.png'],'-transparent'); close;

%% adjust float pH via compariaon with ESPER estimates (profile by profile)
unique_profs = unique(float_data.PROF_ID);
float_data.PH_ESPER_ADJUSTED = nan(size(float_data.PH));
adjustment = nan(size(unique_profs));
latitude = nan(size(unique_profs));
time = nan(size(unique_profs));
% % adjustment by full-column pH v. delta_pH fit
% for p = 1:length(unique_profs)
%     idx_p = find(float_data.PROF_ID == unique_profs(p) & ~isnan(float_data.PH) & ~isnan(float_data.pH_est));
%     fit_params = polyfit(float_data.PH(idx_p),float_data.PH(idx_p)-float_data.pH_est(idx_p),1);
%     float_data.PH_ESPER_ADJUSTED(idx_p) = float_data.PH(idx_p) - (float_data.PH(idx_p).*fit_params(1) + fit_params(2));
%     adjustment(p) = mean((float_data.PH(idx_p).*fit_params(1) + fit_params(2)));
%     latitude(p) = mean(float_data.LAT(idx_p));
%     time(p) = mean(float_data.TIME(idx_p));
% end
% adjustment by deep mean offset (w/ temp ratio)
for p = 1:length(unique_profs)
    idx_p = find(float_data.PROF_ID == unique_profs(p) & ~isnan(float_data.PH) & ~isnan(float_data.pH_est));
    idx_p_adj = find(float_data.PROF_ID == unique_profs(p) & ~isnan(float_data.PH) & ~isnan(float_data.pH_est) & float_data.PRES > 1000);
    adjustment(p) = mean(float_data.PH(idx_p_adj)-float_data.pH_est(idx_p_adj));
    % figure; plot(float_data.PH(idx_p),float_data.PRES(idx_p),float_data.pH_est(idx_p),float_data.PRES(idx_p));
    t_cor = (mean(float_data.TEMP(idx_p_adj))+273.15)./(float_data.TEMP(idx_p)+273.15);
    float_data.PH_ESPER_ADJUSTED(idx_p) = float_data.PH(idx_p) - adjustment(p).*t_cor;
    % figure; plot(float_data.PH_ESPER_ADJUSTED(idx_p),float_data.PRES(idx_p),float_data.pH_est(idx_p),float_data.PRES(idx_p));
    latitude(p) = mean(float_data.LAT(idx_p));
    time(p) = mean(float_data.TIME(idx_p));
end

%% plot adjustments over time and latitude
% set bins
xEdges = datenum(2012:2026',0,0);
yEdges = (-90:10:90)';
numXBins = length(xEdges) - 1;
numYBins = length(yEdges) - 1;
% get the bin index for each data point
[~, ~, ~, binX, binY] = histcounts2(time,latitude,xEdges,yEdges);
% filter out any points that fell outside the specified edges (where index is 0)
validMask = (binX > 0) & (binY > 0);
% accumulate Z values directly into a 2D grid using their mean
% Note: y-indices match rows (1st dim) and x-indices match columns (2nd dim)
avgMatrix = accumarray([binY(validMask), binX(validMask)], adjustment(validMask), ...
    [numYBins, numXBins], @mean, NaN);
% plot the results using pcolor
figure('color','w','visible','on');
pcolor(xEdges, yEdges, [avgMatrix, nan(numYBins, 1); nan(1, numXBins + 1)]);
shading flat; clim([-0.015 0.015]);
colormap(cmocean('balance','pivot',0));
datetick('x','keeplimits');
colorbar;
ylabel('Latitude');
title('Average ESPER-based pH Adjustment');
text(xEdges(2),yEdges(end-1),['Avg. pH adj. = ' num2str(round(mean(adjustment,'omitnan'),3))]);
export_fig(gcf,['DIC/Figures/Data/adjustment_histogram_' ...
    file_date float_file_ext '.png']); close;

%% determine average profiles
pres = unique(float_data.PRES);
pH_by_depth = nan(length(pres),2);
pH_est_by_depth = nan(length(pres),2);
pH_adj_by_depth = nan(length(pres),2);
for p = 1:length(pres)
    idx_z = float_data.PRES == pres(p);
    pH_by_depth(p,1) = mean(float_data.PH(idx_z),'omitnan');
    pH_by_depth(p,2) = std(float_data.PH(idx_z),[],'omitnan');
    pH_est_by_depth(p,1) = mean(float_data.pH_est(idx_z),'omitnan');
    pH_est_by_depth(p,2) = std(float_data.pH_est(idx_z),[],'omitnan');
    pH_adj_by_depth(p,1) = mean(float_data.PH_ESPER_ADJUSTED(idx_z),'omitnan');
    pH_adj_by_depth(p,2) = std(float_data.PH_ESPER_ADJUSTED(idx_z),[],'omitnan');
end

%% plot average float pH profile vs. ESPER Estimate
clrs = orderedcolors('gem'); fntsz = 14;
figure('visible','on','position',[100 100 1000 1000]);
ax1 = axes('position',[0.05 0.05 0.4 0.9],'Box','on','LineWidth',2);
set(ax1,'YDir','reverse','fontsize',fntsz,'XAxisLocation','top'); hold on;
% fill(ax1,[pH_by_depth(:,1)-pH_by_depth(:,2);flipud(sum(pH_by_depth,2))],...
%     [pres;flipud(pres)],clrs(1,:),'FaceAlpha',0.25,'LineStyle','none');
% fill(ax1,[pH_est_by_depth(:,1)-pH_est_by_depth(:,2);flipud(sum(pH_est_by_depth,2))],...
%     [pres;flipud(pres)],clrs(2,:),'FaceAlpha',0.25,'LineStyle','none');
% fill(ax1,[pH_adj_by_depth(:,1)-pH_adj_by_depth(:,2);flipud(sum(pH_adj_by_depth,2))],...
%     [pres;flipud(pres)],clrs(3,:),'FaceAlpha',0.25,'LineStyle','none');
%     plot(pH_by_depth(:,1),pres,pH_est_by_depth(:,1),pres,pH_adj_by_depth(:,1),pres,'linewidth',2);
plot(pH_by_depth(:,1),pres,pH_est_by_depth(:,1),pres,pH_adj_by_depth(:,1),pres,'linewidth',2);
legend({'Measured pH' 'ESPER pH' 'Adjusted pH'},'location','southeast');
ylabel(ax1,'Depth (dbar)'); xlabel('pH'); xlim([7.75 8.05]);
ax2 = axes('position',[0.55 0.05 0.4 0.9],'Box','on','LineWidth',2);
set(ax2,'YDir','reverse','fontsize',fntsz,'XAxisLocation','top'); hold on;
% fill(ax2,[pH_by_depth(:,1)-pH_by_depth(:,2);flipud(sum(pH_by_depth,2))],...
%     [pres;flipud(pres)],clrs(1,:),'FaceAlpha',0.25,'LineStyle','none');
% fill(ax2,[pH_est_by_depth(:,1)-pH_est_by_depth(:,2);flipud(sum(pH_est_by_depth,2))],...
%     [pres;flipud(pres)],clrs(2,:),'FaceAlpha',0.25,'LineStyle','none');
plot(pH_by_depth(:,1)-pH_est_by_depth(:,1),pres,...
    pH_adj_by_depth(:,1)-pH_est_by_depth(:,1),pres,'linewidth',2);
plot([0 0],[0 2000],'--k');
xlim([-0.015 0.015]);
legend({'Measured pH - ESPER pH' 'Adjusted pH - ESPER pH' ''},'location','southwest');
ylabel(ax1,'Depth (dbar)'); xlabel('\DeltapH');
export_fig(gcf,['DIC/Figures/Data/pH_meas_vs_esper_profile_' ...
    file_date float_file_ext '.png'],'-transparent'); close;

%%  calculate DIC from ADJUSTED float data
float_data.PH(float_data.PH_ESPER_ADJUSTED > 8.5) = NaN; % eliminate obviously bad pH values
carb = CO2SYS(float_data.TA_est,float_data.PH_ESPER_ADJUSTED,1,3,float_data.SAL,...
    float_data.TEMP,NaN,float_data.PRES,NaN,float_data.SIL_est,float_data.PHOS_est,...
    0,0,1,10,1,2,2);
float_data.DIC = carb(:,2);
float_data.DIC(float_data.DIC==-999) = NaN;
idx = ~isnan(float_data.PH_ESPER_ADJUSTED) & ~isnan(float_data.pH_est);

%% plot measured (adjusted) vs. estimated pH
figure('visible','on','position',[100 100 1000 400]);
set(gca,'fontsize',15,'Box','on','LineWidth',2); hold on;
histogram2(gca,float_data.PH_ESPER_ADJUSTED,float_data.pH_est,'DisplayStyle', 'tile');
plot([7.4 8.3],[7.4 8.3],'b');
xlim([7.4 8.3]); ylim([7.4 8.3]);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
clim([1e0 1e4]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text(7.5,8.2,['Average Diff. = ' ...
    num2str(round(mean(float_data.PH_ESPER_ADJUSTED(idx)-float_data.pH_est(idx)),4)) ...
    ' +/- ' num2str(round(std(float_data.PH_ESPER_ADJUSTED(idx)-float_data.pH_est(idx)),4))],...
    'FontSize',12);
xlabel('Float measured pH (Adjusted to ESPER)'); ylabel('Estimated pH (ESPER w/ T/S/O2)');
export_fig(gcf,['DIC/Figures/Data/pH_meas_vs_esper_corr_' ...
    file_date float_file_ext '.png'],'-transparent'); close;

%% Co-locate and compare float and glodap data (via profile crossovers)
pres_levels = [2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975]';
float_profile_IDs = unique(float_data.PROF_ID);
% pre-allocate
crossover.id = [];
crossover.lon = [];
crossover.lat = [];
crossover.pres = [];
crossover.sigma_float = [];
crossover.dic_float = [];
crossover.dic_float_gradient = [];
crossover.temp_float = [];
crossover.sal_float = [];
crossover.dic_glodap = [];
crossover.dic_glodap_gradient = [];
crossover.sigma_glodap = [];
crossover.dic_delta = [];
crossover.dic_delta_per = [];
% cycle through float profiles 
for p = 1:length(float_profile_IDs)
    % float position
    lon = mean(float_data.LON(float_data.PROF_ID==float_profile_IDs(p)));
    lat = mean(float_data.LAT(float_data.PROF_ID==float_profile_IDs(p)));
    time = mean(float_data.TIME(float_data.PROF_ID==float_profile_IDs(p)));
    % match index
    idx_lon = abs(glodap_data.LON - lon) <= 0.5;
    idx_lat = abs(glodap_data.LAT - lat) <= 0.5;
    idx_time = abs(glodap_data.TIME - time) <= 365/2;
    % idx_time = true(size(glodap_data.TIME));
    idx_all = idx_lon & idx_lat & idx_time;
    % if there is a matching glodap profile
    if sum(idx_all) > 0
        % float dic and depth
        dic_float = float_data.DIC(float_data.PROF_ID==float_profile_IDs(p));
        pres_float = float_data.PRES(float_data.PROF_ID==float_profile_IDs(p));
        sigma_float = float_data.SIGMA(float_data.PROF_ID==float_profile_IDs(p));
        lon_float = float_data.LON(float_data.PROF_ID==float_profile_IDs(p));
        lat_float = float_data.LAT(float_data.PROF_ID==float_profile_IDs(p));
        temp_float = float_data.TEMP(float_data.PROF_ID==float_profile_IDs(p));
        sal_float = float_data.SAL(float_data.PROF_ID==float_profile_IDs(p));
        % combine float profiles if there are more than one
        if length(dic_float) > 58
            dic_temp = nan(length(pres_levels),1);
            pres_temp = nan(length(pres_levels),1);
            sigma_temp = nan(length(pres_levels),1);
            lon_temp = nan(length(pres_levels),1);
            lat_temp = nan(length(pres_levels),1);
            temp_temp = nan(length(pres_levels),1);
            sal_temp = nan(length(pres_levels),1);
            for z = 1:length(pres_levels)
                idx = pres_float == pres_levels(z);
                dic_temp(z) = mean(dic_float(idx),'omitnan');
                pres_temp(z) = mean(pres_float(idx),'omitnan');
                sigma_temp(z) = mean(sigma_float(idx),'omitnan');
                lon_temp(z) = mean(lon_float(idx),'omitnan');
                lat_temp(z) = mean(lat_float(idx),'omitnan');
                temp_temp(z) = mean(temp_float(idx),'omitnan');
                sal_temp(z) = mean(sal_float(idx),'omitnan');
            end
            dic_float = dic_temp(~isnan(dic_temp));
            pres_float = pres_temp(~isnan(pres_temp));
            sigma_float = sigma_temp(~isnan(sigma_temp));
            lon_float = lon_temp(~isnan(lon_temp));
            lat_float = lat_temp(~isnan(lat_temp));
            temp_float = temp_temp(~isnan(temp_temp));
            sal_float = sal_temp(~isnan(sal_temp));
        end
        % get dic gradient
        dic_float_gradient = [0;diff(dic_float)];
        % pressure index
        idx_pres = ismember(pres_levels,pres_float);
        % extract glodap data
        dic_glodap_temp = glodap_data.DIC(idx_all);
        sigma_glodap_temp = glodap_data.SIGMA(idx_all);
        % take average of all matching glodap profiles (helps if there are multiple)
        dic_glodap = nan(58,1);
        sigma_glodap = nan(58,1);
        for gld_z = 1:58
            dic_glodap(gld_z) = mean(dic_glodap_temp(gld_z:58:end),'omitnan');
            sigma_glodap(gld_z) = mean(sigma_glodap_temp(gld_z:58:end),'omitnan');
        end
        % get dic gradient
        dic_glodap_gradient = [0;diff(dic_glodap)];
        % log matched data
        crossover.id = [crossover.id;repmat(float_profile_IDs(p),length(pres_float),1)];
        crossover.lon = [crossover.lon;lon_float];
        crossover.lat = [crossover.lat;lat_float];
        crossover.pres = [crossover.pres;pres_float];
        crossover.dic_float = [crossover.dic_float;dic_float];
        crossover.dic_float_gradient = [crossover.dic_float_gradient;dic_float_gradient];
        crossover.temp_float = [crossover.temp_float;temp_float];
        crossover.sal_float = [crossover.sal_float;sal_float];
        crossover.dic_glodap = [crossover.dic_glodap;dic_glodap(idx_pres)];
        crossover.dic_glodap_gradient = [crossover.dic_glodap_gradient;dic_glodap_gradient(idx_pres)];
        crossover.sigma_float = [crossover.sigma_float;sigma_float];
        crossover.sigma_glodap = [crossover.sigma_glodap;sigma_glodap];
        crossover.dic_delta = [crossover.dic_delta;dic_float-dic_glodap(idx_pres)];
        crossover.dic_delta_per = [crossover.dic_delta_per;...
            (dic_float-dic_glodap(idx_pres))./dic_float];
    end
end

% index to below 300 dbars
idx = crossover.pres > 300 & ~isnan(crossover.dic_float) & ~isnan(crossover.dic_glodap);
% clean up
clear pres_levels float_profile_IDs lon lat time dic_float pres_float
clear idx_lon idx_lat idx_time dic_temp pres_temp idx_pres dic_glodap

% fit delta against DIC
mdl = fitlm(crossover.dic_float(idx),crossover.dic_delta(idx),'Intercept',true);
slp = mdl.Coefficients.Estimate(2);
int = mdl.Coefficients.Estimate(1);
% apply linear correction to float DIC
corr_fac = slp .* crossover.dic_float + int;
crossover.dic_float_corr = crossover.dic_float - corr_fac;
% re-calculate delta
crossover.dic_delta_corr = crossover.dic_float_corr-crossover.dic_glodap;
crossover.dic_delta_per_corr = (crossover.dic_float_corr-crossover.dic_glodap)./crossover.dic_float_corr;
% save crossover data
if ~exist([pwd '/DIC/Data'],'dir'); mkdir('Data'); end
save(['DIC/Data/crossover_data_' file_date float_file_ext],'crossover','idx','-v7.3')
% save correction factors
if ~exist([pwd '/DIC/Data'],'dir'); mkdir('DIC/Data'); end
save(['DIC/Data/float_corr_' file_date float_file_ext],'slp','int');
clear slp int

%% Co-locate and compare float and glodap data (via bins)
% % load binned data
% load(['DIC/Data/binned_data_' file_date float_file_ext],'binned_data')
% % index where binned float and glodap data overlap under 300 dbars
% idx = ~isnan(binned_data.dic_float) & ~isnan(binned_data.dic_glodap) & ...
%     binned_data.pres > 300;
% % add matched data to structure
% matched_data.dic_float = binned_data.dic_float(idx);
% matched_data.dic_glodap = binned_data.dic_glodap(idx);
% matched_data.dic_delta = matched_data.dic_glodap-matched_data.dic_float;
% matched_data.dic_delta_per = (matched_data.dic_glodap-matched_data.dic_float)./matched_data.dic_glodap;
% matched_data.pres = binned_data.pres(idx);
% % fit delta against DIC
% mdl = fitlm(matched_data.dic_float,matched_data.dic_delta_per,'Intercept',true);
% slp = mdl.Coefficients.Estimate(2);
% int = mdl.Coefficients.Estimate(1);
% % apply linear correction to float dic
% matched_data.dic_float_corr = matched_data.dic_float + ...
%     (slp.*matched_data.dic_float + int).*matched_data.dic_float;
% % re-calculate delta
% matched_data.dic_delta_corr = matched_data.dic_glodap-matched_data.dic_float_corr;
% matched_data.dic_delta_per_corr = matched_data.dic_glodap-matched_data.dic_float_corr;
% % clean up
% clear idx mdl binned_data
% 
% % save matched float and glodap data
% if ~exist([pwd '/DIC/Data'],'dir'); mkdir('DIC/Data'); end
% save(['DIC/Data/matched_data_' file_date float_file_ext],'matched_data','-v7.3');
% clear matched_data
% % save correction factors
% if ~exist([pwd '/DIC/Data'],'dir'); mkdir('DIC/Data'); end
% save(['DIC/Data/float_corr_' file_date float_file_ext],'slp','int');
% clear slp int

%% plot uncorrected float vs. glodap residuals
load(['DIC/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
scatter(crossover.dic_glodap(idx),crossover.dic_float(idx),20,'.');
plot([2000,2400],[2000 2400],'k--');
ylabel('Binned BGC Argo DIC Data (\mumol kg^{-1})');
xlabel('Binned GLODAP DIC Data (\mumol kg^{-1})');
crossover.err_md = median(crossover.dic_delta(idx));
crossover.std = std(crossover.dic_delta(idx));
crossover.rmse = sqrt(mean((crossover.dic_delta(idx)).^2));
crossover.r2 = corr(crossover.dic_glodap(idx),crossover.dic_float(idx));
crossover.slope = polyfit(crossover.dic_glodap(idx),crossover.dic_float(idx),1);
crossover.slope = crossover.slope(1);
text(2050,2350,['Med. Err. = ' num2str(round(crossover.err_md,2))],'fontsize',12);
text(2050,2325,['RMSE = ' num2str(round(crossover.rmse,1))],'fontsize',12);
text(2050,2300,['R^{2} = ' num2str(round(crossover.r2,2))],'fontsize',12);
exportgraphics(gcf,['DIC/Figures/Data/delta_vs_float_300_uncorr_' file_date float_file_ext '.png']);
close
% histogram
figure; hold on
histogram(crossover.dic_delta(idx),'Normalization','probability');
ylim('manual');
xlim([-40 40]);
set(gca,'fontsize',16);
sig2_min = double(mean(crossover.dic_delta(idx))-2*std(crossover.dic_delta(idx)));
sig2_max = double(mean(crossover.dic_delta(idx))+2*std(crossover.dic_delta(idx)));
plot([sig2_min sig2_min],[0 1],'r','linewidth',1);
plot([sig2_max sig2_max],[0 1],'r','linewidth',1);
xlabel('Float DIC - GLODAP DIC');
xL = xlim; yL = ylim;
text(0.9*xL(1),0.9*yL(2),['Mean \DeltaDIC = ' ...
    num2str(round(mean(crossover.dic_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
text(0.9*xL(1),0.8*yL(2),['Std. \DeltaDIC = ' ...
    num2str(round(std(crossover.dic_delta(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
exportgraphics(gcf,[pwd '/DIC/Figures/Data/GLODAP_comp_histogram_uncorr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot corrected float vs. glodap residuals
load(['DIC/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
scatter(crossover.dic_glodap(idx),crossover.dic_float_corr(idx),20,'.');
plot([2000,2400],[2000 2400],'k--');
ylabel('Corrected, Binned BGC Argo DIC Data (\mumol kg^{-1})');
xlabel('Binned GLODAP DIC Data (\mumol kg^{-1})');
crossover.err_md_corr = median(crossover.dic_delta_corr(idx));
crossover.std_corr = std(crossover.dic_delta_corr(idx));
crossover.rmse_corr = sqrt(mean((crossover.dic_delta_corr(idx)).^2));
crossover.r2_corr = corr(crossover.dic_glodap(idx),crossover.dic_float_corr(idx));
crossover.slope_corr = polyfit(crossover.dic_glodap(idx),crossover.dic_float_corr(idx),1);
crossover.slope_corr = crossover.slope_corr(1);
text(2050,2350,['Med. Err. = ' num2str(round(crossover.err_md_corr,2))],'fontsize',12);
text(2050,2325,['RMSE = ' num2str(round(crossover.rmse_corr,1))],'fontsize',12);
text(2050,2300,['R^{2} = ' num2str(round(crossover.r2_corr,2))],'fontsize',12);
exportgraphics(gcf,['DIC/Figures/Data/delta_vs_float_300_corr_' file_date float_file_ext '.png']);
close
% histogram
figure; hold on
histogram(crossover.dic_delta_corr(idx),'Normalization','probability');
ylim('manual');
xlim([-40 40]);
set(gca,'fontsize',16);
sig2_min = double(mean(crossover.dic_delta_corr(idx))-2*std(crossover.dic_delta_corr(idx)));
sig2_max = double(mean(crossover.dic_delta_corr(idx))+2*std(crossover.dic_delta_corr(idx)));
plot([sig2_min sig2_min],[0 1],'r','linewidth',1);
plot([sig2_max sig2_max],[0 1],'r','linewidth',1);
xlabel('Float DIC - GLODAP DIC');
xL = xlim; yL = ylim;
text(0.9*xL(1),0.9*yL(2),['Mean \DeltaDIC = ' ...
    num2str(round(mean(crossover.dic_delta_corr(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
text(0.9*xL(1),0.8*yL(2),['Std. \DeltaDIC = ' ...
    num2str(round(std(crossover.dic_delta_corr(idx)),2))],...
    'HorizontalAlignment','left','VerticalAlignment','top');
exportgraphics(gcf,[pwd '/DIC/Figures/Data/GLODAP_comp_histogram_corr_' file_date float_file_ext '.png']);
close
% clean up
clear crossover

%% plot mapped crossover data by profile mean delta
load(['DIC/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
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
        mean(crossover.dic_delta(idx),'omitnan'),'filled',...
        'markeredgecolor','k');
end
colorbar;
if ~exist('DIC/Figures/Data','dir'); mkdir('pH/Figures/Data'); end
export_fig(gcf,'DIC/Figures/Data/mapped_comparison_uncorr.png','-transparent');
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
        mean(crossover.dic_delta_corr(idx),'omitnan'),'filled',...
        'markeredgecolor','k');
end
colorbar;
if ~exist('DIC/Figures/Data','dir'); mkdir('pH/Figures/Data'); end
export_fig(gcf,'DIC/Figures/Data/mapped_comparison_corr.png','-transparent');
close
% clean up
clear crossover lon_temp

%% adjust and save float data
load(['DIC/Data/float_corr_' file_date float_file_ext],'slp','int');
vars = fieldnames(float_data);
for v = 1:length(vars)      
    float_data_adjusted.(vars{v}) = float_data.(vars{v});
end

if include_float
    % save adjusted float data
    if ~exist([pwd '/DIC/Data'],'dir'); mkdir('DIC/Data'); end
    save(['DIC/Data/processed_float_dic_data_adjusted_' ...
        file_date float_file_ext '.mat'],...
        'float_data_adjusted','file_date','-v7.3');
    clear slp int float_data float_data_adjusted v vars
end

if include_glodap
    % save adjusted GLODAP data
    if ~exist([pwd '/DIC/Data'],'dir'); mkdir('DIC/Data'); end
    save(['DIC/Data/processed_glodap_dic_data_adjusted_' ...
        glodap_vrs '.mat'],...
        'glodap_data','file_date','-v7.3');
end

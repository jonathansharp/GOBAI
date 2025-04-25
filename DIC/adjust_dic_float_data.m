% adjust_dic_float_data
%
% DESCRIPTION:
% This function bins float and glodap data, co-locates corresponding bins,
% compares the two datasets, and calculates a bulk correction factor based
% on the discrepancies between the two.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 04/21/2025

function adjust_dic_float_data(float_file_ext,glodap_year,snap_date)

%% load interpolated float and glodap data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['DIC/Data/processed_float_dic_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date');
load(['DIC/Data/processed_glodap_dic_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');

%% calculate DIC from float data
% convert pressure to depth
float_data.DEPTH = -gsw_z_from_p(float_data.PRES,float_data.LAT);
% make ESPER estimates
esper_preds = ESPER_Mixed([1 4 6],[float_data.LON float_data.LAT float_data.DEPTH],...
    [float_data.TEMP float_data.SAL float_data.OXY],[2 1 6],'Equations',7);
float_data.TA = esper_preds.TA;
float_data.PHOS = esper_preds.phosphate;
float_data.SIL = esper_preds.silicate;
% carbonate system calculations for DIC
float_data.PH(float_data.PH > 8.5) = NaN; % eliminate obviously bad pH values
carb = CO2SYS(float_data.TA,float_data.PH,1,3,float_data.SAL,...
    float_data.TEMP,NaN,float_data.PRES,NaN,float_data.SIL,float_data.PHOS,...
    0,0,1,10,1,2,2);
float_data.DIC = carb(:,2);
float_data.DIC(float_data.DIC==-999) = NaN;

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
% binned_data.dic_float = single(nan(sz));
% binned_data.dic_glodap = single(nan(sz));
% for m = 1:length(t_bins)
%     % month-specific float index
%     idx_float_tmp = idx_float;
%     idx_float_tmp(subs_float(:,4)~=m) = false;
%     % month-specific glodap index
%     idx_glodap_tmp = idx_glodap;
%     idx_glodap_tmp(subs_glodap(:,4)~=m) = false;
%     % bin dic data
%     binned_data.dic_float(:,:,:,m) = single(accumarray(subs_float(idx_float_tmp,1:3),float_data.DIC(idx_float_tmp),sz(1:3),@nanmean,nan));
%     binned_data.dic_glodap(:,:,:,m) = single(accumarray(subs_glodap(idx_glodap_tmp,1:3),glodap_data.DIC(idx_glodap_tmp),sz(1:3),@nanmean,nan));
% end
% % add pressure bins
% binned_data.pres = repmat(permute(single(z_bins),[3 1 2]),length(x_bins),length(y_bins),1,length(t_bins));
% % save binned data
% if ~exist([pwd '/DIC/Data'],'dir'); mkdir('DIC/Data'); end
% save(['DIC/Data/binned_data_' file_date float_file_ext],'binned_data','-v7.3')
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
    idx_time = abs(glodap_data.TIME - time) <= 30;
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

% save
if ~exist('/DIC/Data','dir'); mkdir('DIC/Data'); end
save(['DIC/Data/processed_float_dic_data_adjusted_' file_date float_file_ext '.mat'],...
    'float_data_adjusted','file_date','-v7.3');
clear slp int float_data float_data_adjusted v vars

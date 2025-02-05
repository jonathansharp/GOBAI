% adjust_no3_float_data
%
% DESCRIPTION:
% This function bins float and glodap data, co-locates corresponding bins,
% compares the two datasets, and calculates a bulk correction factor based
% on the discrepancies between the two.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 3/28/2024

function adjust_no3_float_data(float_file_ext,file_date,glodap_year)

%% load interpolated float and glodap data
load(['NO3/Data/processed_float_no3_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date');
load(['NO3/Data/processed_glodap_no3_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');

%% import WOA climatologies
temp_path = [pwd '/Data/WOA/TEMPERATURE/'];
sal_path = [pwd '/Data/WOA/SALINITY/'];
nit_path = [pwd '/Data/WOA/NITRATE/'];
for m = 1:12
    WOA.TMP(:,:,:,m) = ncread([temp_path 'woa18_decav_t' sprintf('%02d',m) '_01.nc'],'t_an');
    WOA.SAL(:,:,:,m) = ncread([sal_path 'woa18_decav_s' sprintf('%02d',m) '_01.nc'],'s_an');
    WOA.NIT(:,:,:,m) = ncread([nit_path 'woa18_all_n',sprintf('%02d',m),'_01.nc'],'n_an');
end
% Import and format WOA latitude and longitude
WOA.LAT   = ncread([nit_path 'woa18_all_n01_01.nc'],'lat');
WOA.LON   = ncread([nit_path 'woa18_all_n01_01.nc'],'lon');
WOA.DEPTH = ncread([nit_path 'woa18_all_n01_01.nc'],'depth');
% clean up
clear m temp_path sal_path nit_path

%% compare float data to WOA
% index to valid data points
idx = find(~isnan(float_data.NIT) & ~isnan(float_data.LAT) & ...
    ~isnan(float_data.LON) & ~isnan(float_data.PRES) & ...
    ~isnan(float_data.TIME));
% obtain month from date
date = datevec(datenum(float_data.YEAR,0,float_data.DAY));
float_data.MONTH = date(:,2);
clear date
% pre-allocate matches
WOA_match = nan(size(float_data.NIT));
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
    WOA_match(idx(i)) = WOA.NIT(idx_lon(1),idx_lat(1),idx_pres(1),idx_mnth(1));
end
% calculate differences
WOA_delta = float_data.NIT-WOA_match;
WOA_delta_per = WOA_delta./float_data.NIT;
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
    hist3([WOA_match,WOA_delta],'Edges',{0:1:46 -34.5:1.5:34.5});
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
% establish edges of bins
x_edges = -180:180;
x_bins = -179.5:179.5;
y_edges = -85:85;
y_bins = -84.5:84.5;
z_edges = [0 5:10:175 190:20:450 475:50:1375 1450:100:1950 2000];
z_bins = [2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975];
mn_edges = 1:12;
yr_edges = 2004:glodap_year;
t_edges = datenum([[repelem(yr_edges,1,length(mn_edges)) yr_edges(end)+1]', ...
                  [repmat(mn_edges,1,length(yr_edges)) 1]', ...
                  [zeros(1,length(mn_edges)*length(yr_edges)+1)]']);
t_bins = datenum([[repelem(yr_edges,1,length(mn_edges))]', ...
                  [repmat(mn_edges,1,length(yr_edges))]', ...
                  [repmat(15,1,length(mn_edges)*length(yr_edges))]']);
% get histogram counts in each bin
[~,~,Xnum_float] = histcounts(float_data.LON,x_edges);
[~,~,Ynum_float] = histcounts(float_data.LAT,y_edges);
[~,~,Znum_float] = histcounts(float_data.PRES,z_edges);
[~,~,Tnum_float] = histcounts(float_data.TIME,t_edges);
[~,~,Xnum_glodap] = histcounts(glodap_data.LON,x_edges);
[~,~,Ynum_glodap] = histcounts(glodap_data.LAT,y_edges);
[~,~,Znum_glodap] = histcounts(glodap_data.PRES,z_edges);
[~,~,Tnum_glodap] = histcounts(glodap_data.TIME,t_edges);
% accumulate index of counts
subs_float = [Xnum_float,Ynum_float,Znum_float,Tnum_float];
idx_float = ~any(subs_float==0,2);
subs_glodap = [Xnum_glodap,Ynum_glodap,Znum_glodap,Tnum_glodap];
idx_glodap = ~any(subs_glodap==0,2);
clear Xnum_float Ynum_float Znum_float Tnum_float
clear Xnum_glodap Ynum_glodap Znum_glodap Tnum_glodap
% determine size of 4D grid
sz = [length(x_bins),length(y_bins),length(z_bins),length(t_bins)];
% clean up
clear x_edges y_edges z_edges mn_edges yr_edges t_edges

%% Bin float and glodap data
binned_data.nit_float = single(nan(sz));
binned_data.nit_glodap = single(nan(sz));
for m = 1:length(t_bins)
    % month-specific float index
    idx_float_tmp = idx_float;
    idx_float_tmp(subs_float(:,4)~=m) = false;
    % month-specific glodap index
    idx_glodap_tmp = idx_glodap;
    idx_glodap_tmp(subs_glodap(:,4)~=m) = false;
    % bin nitrate data
    binned_data.nit_float(:,:,:,m) = single(accumarray(subs_float(idx_float_tmp,1:3),float_data.NIT(idx_float_tmp),sz(1:3),@nanmean,nan));
    binned_data.nit_glodap(:,:,:,m) = single(accumarray(subs_glodap(idx_glodap_tmp,1:3),glodap_data.NIT(idx_glodap_tmp),sz(1:3),@nanmean,nan));
end
% add pressure bins
binned_data.pres = repmat(permute(single(z_bins),[3 1 2]),length(x_bins),length(y_bins),1,length(t_bins));
% save binned data
if ~exist([pwd '/NO3/Data'],'dir'); mkdir('NO3/Data'); end
save(['NO3/Data/binned_data_' file_date float_file_ext],'binned_data','-v7.3')

% clean up
clear binned_data idx_float idx_glodap
clear idx_float_tmp idx_glodap_tmp m sz subs_float subs_glodap
clear x_edges x_bins y_bins z_edges z_bins t_bins

%% Co-locate and compare float and glodap data (via profile crossovers)
pres_levels = [2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975]';
float_profile_IDs = unique(float_data.PROF_ID);
% pre-allocate
crossover.pres = [];
crossover.nit_float = [];
crossover.temp_float = [];
crossover.sal_float = [];
crossover.nit_glodap = [];
crossover.nit_delta = [];
crossover.nit_delta_per = [];
% cycle through float profiles 
for p = 1:length(float_profile_IDs)
    % float position
    lon = mean(float_data.LON(float_data.PROF_ID==float_profile_IDs(p)));
    lat = mean(float_data.LAT(float_data.PROF_ID==float_profile_IDs(p)));
    time = mean(float_data.TIME(float_data.PROF_ID==float_profile_IDs(p)));
    % match index
    idx_lon = abs(glodap_data.LON - lon) < 0.25;
    idx_lat = abs(glodap_data.LAT - lat) < 0.25;
    idx_time = abs(glodap_data.TIME - time) < 30;
    idx_all = idx_lon & idx_lat & idx_time;
    % if there is a matching glodap profile
    if sum(idx_all) > 0
        % float nitrate and depth
        nit_float = float_data.NIT(float_data.PROF_ID==float_profile_IDs(p));
        pres_float = float_data.PRES(float_data.PROF_ID==float_profile_IDs(p));
        temp_float = float_data.TEMP(float_data.PROF_ID==float_profile_IDs(p));
        sal_float = float_data.SAL(float_data.PROF_ID==float_profile_IDs(p));
        % combine float profiles if there are more than one
        if length(nit_float) > 58
            nit_temp = nan(length(pres_levels),1);
            pres_temp = nan(length(pres_levels),1);
            temp_temp = nan(length(pres_levels),1);
            sal_temp = nan(length(pres_levels),1);
            for z = 1:length(pres_levels)
                idx = pres_float == pres_levels(z);
                nit_temp(z) = mean(nit_float(idx),'omitnan');
                pres_temp(z) = mean(pres_float(idx),'omitnan');
                temp_temp(z) = mean(temp_float(idx),'omitnan');
                sal_temp(z) = mean(sal_float(idx),'omitnan');
            end
            nit_float = nit_temp(~isnan(nit_temp));
            pres_float = pres_temp(~isnan(pres_temp));
            temp_float = temp_temp(~isnan(temp_temp));
            sal_float = sal_temp(~isnan(sal_temp));
        end
        % pressure index
        idx_pres = ismember(pres_levels,pres_float);
        % extract glodap data
        nit_glodap = glodap_data.NIT(idx_all);
        % log matched data
        crossover.pres = [crossover.pres;pres_float];
        crossover.nit_float = [crossover.nit_float;nit_float];
        crossover.temp_float = [crossover.temp_float;temp_float];
        crossover.sal_float = [crossover.sal_float;sal_float];
        crossover.nit_glodap = [crossover.nit_glodap;nit_glodap(idx_pres)];
        crossover.nit_delta = [crossover.nit_delta;nit_float-nit_glodap(idx_pres)];
        crossover.nit_delta_per = [crossover.nit_delta_per;...
            (nit_float-nit_glodap(idx_pres))./nit_float];
    end
end

% index to below 300 dbars
idx = crossover.pres > 300 & ~isnan(crossover.nit_float) & ~isnan(crossover.nit_glodap);
% clean up
clear pres_levels float_profile_IDs lon lat time nit_float pres_float
clear idx_lon idx_lat idx_time nit_temp pres_temp idx_pres nit_glodap
% fit delta against [NO3]
mdl = fitlm(crossover.nit_float(idx),crossover.nit_delta_per(idx),'Intercept',true);
slp = mdl.Coefficients.Estimate(2);
int = mdl.Coefficients.Estimate(1);
% apply linear correction to float NO3
corr_fac = slp .* crossover.nit_float + int; % modelled percent delta
crossover.nit_float_corr = crossover.nit_float - (corr_fac.*crossover.nit_float);
% re-calculate delta
crossover.nit_delta_corr = crossover.nit_float_corr-crossover.nit_glodap;
crossover.nit_delta_per_corr = (crossover.nit_float_corr-crossover.nit_glodap)./crossover.nit_float_corr;
% save crossover data
if ~exist([pwd '/NO3/Data'],'dir'); mkdir('Data'); end
save(['NO3/Data/crossover_data_' file_date float_file_ext],'crossover','idx','-v7.3')
% save correction factors
if ~exist([pwd '/NO3/Data'],'dir'); mkdir('NO3/Data'); end
save(['NO3/Data/float_corr_' file_date float_file_ext],'slp','int');
clear slp int
% clean up

%% Co-locate and compare float and glodap data (via bins)
% % load binned data
% load(['NO3/Data/binned_data_' file_date float_file_ext],'binned_data')
% % index where binned float and glodap data overlap under 300 dbars
% idx = ~isnan(binned_data.nit_float) & ~isnan(binned_data.nit_glodap) & ...
%     binned_data.pres > 300;
% % add matched data to structure
% matched_data.nit_float = binned_data.nit_float(idx);
% matched_data.nit_glodap = binned_data.nit_glodap(idx);
% matched_data.nit_delta = matched_data.nit_glodap-matched_data.nit_float;
% matched_data.nit_delta_per = (matched_data.nit_glodap-matched_data.nit_float)./matched_data.nit_glodap;
% matched_data.pres = binned_data.pres(idx);
% % fit delta against [NO3]
% mdl = fitlm(matched_data.nit_float,matched_data.nit_delta_per,'Intercept',true);
% slp = mdl.Coefficients.Estimate(2);
% int = mdl.Coefficients.Estimate(1);
% % apply linear correction to float NO3
% matched_data.nit_float_corr = matched_data.nit_float + ...
%     (slp.*matched_data.nit_float + int).*matched_data.nit_float;
% % re-calculate delta
% matched_data.nit_delta_corr = matched_data.nit_glodap-matched_data.nit_float_corr;
% matched_data.nit_delta_per_corr = matched_data.nit_glodap-matched_data.nit_float_corr;
% % clean up
% clear idx mdl binned_data
% 
% save matched float and glodap data
% if ~exist([pwd '/NO3/Data'],'dir'); mkdir('NO3/Data'); end
% save(['NO3/Data/matched_data_' file_date float_file_ext],'matched_data','-v7.3');
% clear matched_data
% % save correction factors
% if ~exist([pwd '/NO3/Data'],'dir'); mkdir('NO3/Data'); end
% save(['NO3/Data/float_corr_' file_date float_file_ext],'slp','int');
% clear slp int

%% plot uncorrected float vs. glodap residuals
load(['NO3/Data/crossover_data_' file_date float_file_ext],'crossover','idx');
% scatter
figure; hold on;
scatter(crossover.nit_glodap(idx),crossover.nit_float(idx),20,'.');
plot([0,45],[0 45],'k--');
ylabel('Binned BGC Argo Nitrate Data (\mumol kg^{-1})');
xlabel('Binned GLODAP Nitrate Data (\mumol kg^{-1})');
crossover.err_md = median(crossover.nit_delta(idx));
crossover.std = std(crossover.nit_delta(idx));
crossover.rmse = sqrt(mean((crossover.nit_delta(idx)).^2));
crossover.r2 = corr(crossover.nit_glodap(idx),crossover.nit_float(idx));
crossover.slope = polyfit(crossover.nit_glodap(idx),crossover.nit_float(idx),1);
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
xlabel('GLODAP [NO_{3}] - Float [NO_{3}]');
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
scatter(crossover.nit_glodap(idx),crossover.nit_float_corr(idx),20,'.');
plot([0,45],[0 45],'k--');
ylabel('Corrected, Binned BGC Argo Nitrate Data (\mumol kg^{-1})');
xlabel('Binned GLODAP Nitrate Data (\mumol kg^{-1})');
crossover.err_md_corr = median(crossover.nit_delta_corr(idx));
crossover.std_corr = std(crossover.nit_delta_corr(idx));
crossover.rmse_corr = sqrt(mean((crossover.nit_delta_corr(idx)).^2));
crossover.r2_corr = corr(crossover.nit_glodap(idx),crossover.nit_float_corr(idx));
crossover.slope_corr = polyfit(crossover.nit_glodap(idx),crossover.nit_float_corr(idx),1);
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
xlabel('GLODAP [NO_{3}] - Float [NO_{3}]');
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
load(['NO3/Data/float_corr_' file_date float_file_ext],'slp','int');
vars = fieldnames(float_data);
for v = 1:length(vars)
    if strcmp(vars{v},'NIT')
        float_data_adjusted.NIT = double(float_data.NIT - ...
            (slp.*float_data.NIT + int).*float_data.NIT);
    else
        float_data_adjusted.(vars{v}) = float_data.(vars{v});
    end
end
% save
if ~exist([pwd '/NO3/Data'],'dir'); mkdir('NO3/Data'); end
save(['NO3/Data/processed_float_no3_data_adjusted_' file_date float_file_ext '.mat'],...
    'float_data_adjusted','file_date','-v7.3');
clear slp int float_data float_data_adjusted v vars

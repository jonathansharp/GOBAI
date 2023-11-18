% adjust_float_data
%
% DESCRIPTION:
% This function bins float and glodap data, co-locates corresponding bins,
% compares the two datasets, and calculates a bulk correction factor based
% on the discrepancies between the two.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/22/2023

%% load interpolated float and glodap data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_float_o2_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date');
load(['Data/processed_glodap_o2_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');

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
idx = find(~isnan(float_data.OXY) & ~isnan(float_data.OXY_LAT) & ...
    ~isnan(float_data.OXY_LON) & ~isnan(float_data.OXY_PRES) & ...
    ~isnan(float_data.OXY_TIME));
% obtain month from date
date = datevec(datenum(float_data.OXY_YEAR,0,float_data.OXY_DAY));
float_data.OXY_MONTH = date(:,2);
clear date
% pre-allocate matches
WOA_match = nan(size(float_data.OXY));
% match float data with WOA
for i=1:length(idx)
    idx_lon = find(min(abs(WOA.LON - float_data.OXY_LON(idx(i)))) == ...
        abs(WOA.LON - float_data.OXY_LON(idx(i))));
    idx_lat = find(min(abs(WOA.LAT - float_data.OXY_LAT(idx(i)))) == ...
        abs(WOA.LAT - float_data.OXY_LAT(idx(i))));
    idx_pres = find(min(abs(WOA.DEPTH - float_data.OXY_PRES(idx(i)))) == ...
        abs(WOA.DEPTH - float_data.OXY_PRES(idx(i))));
    idx_mnth = find(min(abs((1:12)' - float_data.OXY_MONTH(idx(i)))) == ...
        abs((1:12)' - float_data.OXY_MONTH(idx(i))));
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
plot([mean_delta+5*st_dev_delta mean_delta+5*st_dev_delta],[0 80000],'r','linewidth',2);
plot([mean_delta-5*st_dev_delta mean_delta-5*st_dev_delta],[0 80000],'r','linewidth',2);
xlabel('Float [O_{2}] - WOA [O_{2}]');
exportgraphics(gcf,[pwd '/Figures/WOA_comp_histogram_' file_date float_file_ext '.png']);
close
% scatter of delta values
figure; hold on
set(gca,'fontsize',16);
xlabel('WOA [O_{2}]');
ylabel('Float [O_{2}] - WOA [O_{2}]');
[counts,bin_centers] = ...
    hist3([WOA_match,WOA_delta],'Edges',{0:10:460 -345:15:345});
h=pcolor(bin_centers{1},bin_centers{2},counts');
plot([0 450],[0 0],'k--');
set(h,'EdgeColor','none');
xlim([0 450]); ylim([-340 340]);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log')
caxis([1e0 1e5]);
c=colorbar;
c.Label.String = 'log10(Bin Counts)';
exportgraphics(gcf,[pwd '/Figures/WOA_comp_scatter_' file_date float_file_ext '.png']);
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
yr_edges = 2004:2022;
t_edges = datenum([[repelem(yr_edges,1,length(mn_edges)) yr_edges(end)+1]', ...
                  [repmat(mn_edges,1,length(yr_edges)) 1]', ...
                  [zeros(1,length(mn_edges)*length(yr_edges)+1)]']);
t_bins = datenum([[repelem(yr_edges,1,length(mn_edges))]', ...
                  [repmat(mn_edges,1,length(yr_edges))]', ...
                  [repmat(15,1,length(mn_edges)*length(yr_edges))]']);
% get histogram counts in each bin
[~,~,Xnum_float] = histcounts(float_data.OXY_LON,x_edges);
[~,~,Ynum_float] = histcounts(float_data.OXY_LAT,y_edges);
[~,~,Znum_float] = histcounts(float_data.OXY_PRES,z_edges);
[~,~,Tnum_float] = histcounts(float_data.OXY_TIME,t_edges);
[~,~,Xnum_glodap] = histcounts(glodap_data.OXY_LON,x_edges);
[~,~,Ynum_glodap] = histcounts(glodap_data.OXY_LAT,y_edges);
[~,~,Znum_glodap] = histcounts(glodap_data.OXY_PRES,z_edges);
[~,~,Tnum_glodap] = histcounts(glodap_data.OXY_TIME,t_edges);
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
binned_data.oxy_float = single(nan(sz));
binned_data.oxy_glodap = single(nan(sz));
for m = 1:length(t_bins)
    % month-specific float index
    idx_float_tmp = idx_float;
    idx_float_tmp(subs_float(:,4)~=m) = false;
    % month-specific glodap index
    idx_glodap_tmp = idx_glodap;
    idx_glodap_tmp(subs_glodap(:,4)~=m) = false;
    % bin oxygen data
    binned_data.oxy_float(:,:,:,m) = single(accumarray(subs_float(idx_float_tmp,1:3),float_data.OXY(idx_float_tmp),sz(1:3),@nanmean,nan));
    binned_data.oxy_glodap(:,:,:,m) = single(accumarray(subs_glodap(idx_glodap_tmp,1:3),glodap_data.OXY(idx_glodap_tmp),sz(1:3),@nanmean,nan));
end
% add pressure bins
binned_data.pres = repmat(permute(single(z_bins),[3 1 2]),length(x_bins),length(y_bins),1,length(t_bins));
if ~exist([pwd '/Data'],'dir'); mkdir('Data'); end
save(['Data/binned_data_' file_date float_file_ext],'binned_data','-v7.3')

% clean up
clear glodap_data binned_data idx_float idx_glodap
clear idx_float_tmp idx_glodap_tmp m sz subs_float subs_glodap
clear x_edges x_bins y_bins z_edges z_bins t_bins

%% Co-locate and compare float and glodap data
% load binned data
load(['Data/binned_data_' file_date float_file_ext],'binned_data')
% index where binned float and glodap data overlap under 300 dbars
idx = ~isnan(binned_data.oxy_float) & ~isnan(binned_data.oxy_glodap) & ...
    binned_data.pres > 300;
% add matched data to structure
matched_data.oxy_float = binned_data.oxy_float(idx);
matched_data.oxy_glodap = binned_data.oxy_glodap(idx);
matched_data.oxy_delta = matched_data.oxy_glodap-matched_data.oxy_float;
matched_data.oxy_delta_per = (matched_data.oxy_glodap-matched_data.oxy_float)./matched_data.oxy_glodap;
matched_data.pres = binned_data.pres(idx);
% fit delta against [O2]
mod = fitlm(matched_data.oxy_float,matched_data.oxy_delta_per,'Intercept',true);
slp = mod.Coefficients.Estimate(2);
int = mod.Coefficients.Estimate(1);
% apply linear correction to float O2
matched_data.oxy_float_corr = matched_data.oxy_float + ...
    (slp.*matched_data.oxy_float + int).*matched_data.oxy_float;
% re-calculate delta
matched_data.oxy_delta_corr = matched_data.oxy_glodap-matched_data.oxy_float_corr;
% clean up
clear idx mod binned_data

%% save matched float and glodap data
if ~exist([pwd '/Data'],'dir'); mkdir('Data'); end
save(['Data/matched_data_' file_date float_file_ext],'matched_data','-v7.3');
clear matched_data

%% save correction factors
if ~exist([pwd '/Data'],'dir'); mkdir('Data'); end
save(['Data/float_corr_' file_date float_file_ext],'slp','int');
clear slp int

%% plot uncorrected float vs. glodap residuals
load(['Data/matched_data_' file_date float_file_ext],'matched_data');
figure; hold on;
scatter(matched_data.oxy_glodap,matched_data.oxy_float,20,'.');
plot([0,450],[0 450],'k--');
ylabel('Binned BGC Argo Oxygen Data (\mumol kg^{-1})');
xlabel('Binned GLODAP Oxygen Data (\mumol kg^{-1})');
matched_data.err_md = median(matched_data.oxy_delta);
matched_data.std = std(matched_data.oxy_delta);
matched_data.rmse = sqrt(mean((matched_data.oxy_delta).^2));
matched_data.r2 = corr(matched_data.oxy_glodap,matched_data.oxy_float);
matched_data.slope = polyfit(matched_data.oxy_glodap,matched_data.oxy_float,1);
matched_data.slope = matched_data.slope(1);
text(300,120,['Med. Err. = ' num2str(round(matched_data.err_md,2))],'fontsize',12);
text(300,90,['RMSE = ' num2str(round(matched_data.rmse,1))],'fontsize',12);
text(300,60,['R^{2} = ' num2str(round(matched_data.r2,2))],'fontsize',12);
% save figure
exportgraphics(gcf,['Figures/delta_vs_float_300_uncorr' file_date float_file_ext '.png']);
close
% clean up
clear matched_data

%% plot corrected float vs. glodap residuals
load(['Data/matched_data_' file_date float_file_ext],'matched_data');
figure; hold on;
scatter(matched_data.oxy_glodap,matched_data.oxy_float_corr,20,'.');
plot([0,450],[0 450],'k--');
ylabel('Corrected, Binned BGC Argo Oxygen Data (\mumol kg^{-1})');
xlabel('Binned GLODAP Oxygen Data (\mumol kg^{-1})');
matched_data.err_md_corr = median(matched_data.oxy_delta_corr);
matched_data.std_corr = std(matched_data.oxy_delta_corr);
matched_data.rmse_corr = sqrt(mean((matched_data.oxy_delta_corr).^2));
matched_data.r2_corr = corr(matched_data.oxy_glodap,matched_data.oxy_float_corr);
matched_data.slope_corr = polyfit(matched_data.oxy_glodap,matched_data.oxy_float_corr,1);
matched_data.slope_corr = matched_data.slope_corr(1);
text(300,120,['Med. Err. = ' num2str(round(matched_data.err_md_corr,2))],'fontsize',12);
text(300,90,['RMSE = ' num2str(round(matched_data.rmse_corr,1))],'fontsize',12);
text(300,60,['R^{2} = ' num2str(round(matched_data.r2_corr,2))],'fontsize',12);
% save figure
exportgraphics(gcf,['Figures/delta_vs_float_300_corr_' file_date float_file_ext '.png']);
close
% clean up
clear matched_data

%% adjust and save float data
load(['Data/float_corr_' file_date float_file_ext],'slp','int');
vars = fieldnames(float_data);
for v = 1:length(vars)
    if strcmp(vars{v},'OXY')
        float_data_adjusted.OXY = double(float_data.OXY + ...
            (slp.*float_data.OXY + int).*float_data.OXY);
    else
        float_data_adjusted.(vars{v}) = float_data.(vars{v});
    end
end
% save
if ~exist([pwd '/Data'],'dir'); mkdir('Data'); end
save(['Data/processed_float_o2_data_adjusted_' file_date float_file_ext '.mat'],...
    'float_data_adjusted','file_date','-v7.3');
clear slp int float_data float_data_adjusted v vars

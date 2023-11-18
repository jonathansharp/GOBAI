% calculate_gridding_uncertainty
%
% DESCRIPTION:
% This function bins float and glodap data and uses the standard deviations
% of observations within bins to estimate uncertainty associated with
% binning to spatiotemporal grid cells.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/12/2023

%% load float and glodap data
load_interpolated_data_to_workspace

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

%% Calculate average sigma in grid cells
grid_uncer.oxy_std = single(nan(sz));
grid_uncer.sig = single(nan(sz));
for m = 1:length(t_bins)
    % month-specific float index
    idx_float_tmp = idx_float;
    idx_float_tmp(subs_float(:,4)~=m) = false;
    % month-specific glodap index
    idx_glodap_tmp = idx_glodap;
    idx_glodap_tmp(subs_glodap(:,4)~=m) = false;
    % calculate stadard deviations of float and glodap oxy and average them
    float_oxy_std = single(accumarray(subs_float(idx_float_tmp,1:3),float_data.OXY(idx_float_tmp),sz(1:3),@nanstd,nan));
    glodap_oxy_std = single(accumarray(subs_glodap(idx_glodap_tmp,1:3),glodap_data.OXY(idx_glodap_tmp),sz(1:3),@nanstd,nan));    
    grid_uncer.oxy_std(:,:,:,m) = mean(cat(4,float_oxy_std,glodap_oxy_std),4,'omitnan');
    % calculate mean of float and glodap sigma and average them
    float_sig = single(accumarray(subs_float(idx_float_tmp,1:3),float_data.OXY_SIGMA(idx_float_tmp),sz(1:3),@nanmean,nan));
    glodap_sig = single(accumarray(subs_glodap(idx_glodap_tmp,1:3),glodap_data.OXY_SIGMA(idx_glodap_tmp),sz(1:3),@nanmean,nan));
    grid_uncer.sig(:,:,:,m) = mean(cat(4,float_sig,glodap_sig),4,'omitnan');
end
save('Data/grid_uncer','grid_uncer','-v7.3')
clear grid_uncer

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
% add dimensional bins
%binned_data.LON = single(repmat(x_bins',1,length(y_bins),length(z_bins),length(t_bins)));
%binned_data.LAT = single(repmat(y_bins,length(x_bins),1,length(z_bins),length(t_bins)));
binned_data.pres = repmat(permute(single(z_bins),[3 1 2]),length(x_bins),length(y_bins),1,length(t_bins));
%binned_data.TIME = single(repmat(permute(t_bins,[4 3 2 1]),length(x_bins),length(y_bins),length(z_bins)));
save('Data/binned_data','binned_data','-v7.3')
clear binned_data

%% clean up
clear float_data glodap_data x_edges x_bins y_edges y_bins
clear z_edges z_bins mn_edges yr_edges t_edges t_bins
clear subs_float idx_float subs_glodap idx_glodap sz 

%% Co-locate and compare float and glodap data
% load binned data
load('Data/binned_data','binned_data')
% index where binned float and glodap data overlap under 300 dbars
idx = ~isnan(binned_data.oxy_float) & ~isnan(binned_data.oxy_glodap) & ...
    binned_data.pres > 300;
% add matched data to structure
matched_data.oxy_float = binned_data.oxy_float(idx);
matched_data.oxy_glodap = binned_data.oxy_glodap(idx);
matched_data.oxy_delta = matched_data.oxy_glodap-matched_data.oxy_float;
% matched_data.LON = binned_data.LON(idx);
% matched_data.LAT = binned_data.LAT(idx);
matched_data.pres = binned_data.pres(idx);
% matched_data.TEMP_FLOAT = binned_data.TEMP_FLOAT(idx);
% matched_data.TEMP_GLODAP = binned_data.TEMP_GLODAP(idx);
% matched_data.TEMP = mean([matched_data.TEMP_FLOAT matched_data.TEMP_GLODAP],2);
% fit delta against [O2]
mod = fitlm(matched_data.oxy_float,matched_data.oxy_delta,'Intercept',false);
slp = mod.Coefficients.Estimate(1);
int = 0; % mod.Coefficients.Estimate(1);
% apply linear correction to float O2
matched_data.oxy_floar_corr = matched_data.oxy_float + slp.*matched_data.oxy_float + int;
% re-calculate delta
matched_data.oxy_delta_corr = matched_data.oxy_glodap-matched_data.oxy_floar_corr;

%% using standard deviations on grid to determine gridding uncertainty
% add dimensional bins
grid_uncer.lon = single(repmat(x_bins',1,length(y_bins),length(z_bins),length(t_bins)));
grid_uncer.lon = double(grid_uncer.lon(grid_uncer.idx));
grid_uncer.lat = single(repmat(y_bins,length(x_bins),1,length(z_bins),length(t_bins)));
grid_uncer.lat = double(grid_uncer.lat(grid_uncer.idx));
grid_uncer.pres = single(repmat(permute(z_bins,[3 1 2]),length(x_bins),length(y_bins),1,length(t_bins)));
grid_uncer.pres = double(grid_uncer.pres(grid_uncer.idx));
% calculate distance from shore
grid_uncer.dist = dist2coast(grid_uncer.lat,grid_uncer.lon);
% add bottom depth
grid_uncer.bot = bottom_depth(grid_uncer.lat,grid_uncer.lon);
% remove standard deviations calculated from nine or fewer measurements
oxygen_count = accumarray(subs_float(idx_float,:),1,sz) + accumarray(subs_glodap(idx_glodap,:),1,sz);
oxygen_count = oxygen_count(grid_uncer.idx);
oxygen_count_idx = oxygen_count < 10;
grid_uncer.oxy_std(oxygen_count_idx) = nan;
grid_uncer.pres(oxygen_count_idx) = nan;
grid_uncer.sig(oxygen_count_idx) = nan;
grid_uncer.bot(oxygen_count_idx) = nan;
% scatter different parameters against std
% figure; scatter(grid_uncer.pres(:),grid_uncer.oxy_std(:));
% fit model of variability vs depth, sigma, and bottom depth
[b,~,~,~,stats] = ...
    regress(grid_uncer.oxy_std,[ones(size(grid_uncer.pres)) ...
            grid_uncer.pres grid_uncer.pres.^2 grid_uncer.sig ...
            grid_uncer.sig.^2 grid_uncer.bot grid_uncer.bot.^2]);
stats(1) % R^2 statistic
save(['std_coefs-' date '-P_S_B'],'b')
% plot figure of gridding uncertainty versus n
grid_vs_n_binned
% clean up
clear oxygen_count oxygen_count_idx grid_uncer b stats

%% Bin float and glodap data

% bin oxygen
binned_data.OXY_FLOAT = single(accumarray(subs_float(idx_float,:),float_data.OXY(idx_float),sz,@nanmean,nan));
binned_data.OXY_GLODAP = single(accumarray(subs_glodap(idx_glodap,:),glodap_data.OXY(idx_glodap),sz,@nanmean,nan));
% bin temperature
binned_data.TEMP_FLOAT = single(accumarray(subs_float(idx_float,:),float_data.OXY_TEMP(idx_float),sz,@nanmean,nan));
binned_data.TEMP_GLODAP = single(accumarray(subs_glodap(idx_glodap,:),glodap_data.OXY_TEMP(idx_glodap),sz,@nanmean,nan));
% bin salinity
binned_data.SAL_FLOAT = single(accumarray(subs_float(idx_float,:),float_data.OXY_SAL(idx_float),sz,@nanmean,nan));
binned_data.SAL_GLODAP = single(accumarray(subs_glodap(idx_glodap,:),glodap_data.OXY_SAL(idx_glodap),sz,@nanmean,nan));
% add dimensional bins
binned_data.LON = single(repmat(x_bins',1,length(y_bins),length(z_bins),length(t_bins)));
binned_data.LAT = single(repmat(y_bins,length(x_bins),1,length(z_bins),length(t_bins)));
binned_data.PRES = single(repmat(permute(z_bins,[3 1 2]),length(x_bins),length(y_bins),1,length(t_bins)));
binned_data.TIME = single(repmat(permute(t_bins,[4 3 2 1]),length(x_bins),length(y_bins),length(z_bins)));

%% clean up
clear float_data glodap_data x_edges x_bins y_edges y_bins
clear z_edges z_bins mn_edges yr_edges t_edges t_bins
clear subs_float idx_float subs_glodap idx_glodap sz 

%% Co-locate and compare float and glodap data
% index where binned float and glodap data overlap
idx = ~isnan(binned_data.OXY_FLOAT) & ~isnan(binned_data.OXY_GLODAP) & ...
    binned_data.PRES > 300;
% add matched data to structure
matched_data.OXY_FLOAT = binned_data.OXY_FLOAT(idx);
matched_data.OXY_GLODAP = binned_data.OXY_GLODAP(idx);
matched_data.OXY_DELTA = matched_data.OXY_GLODAP-matched_data.OXY_FLOAT;
matched_data.LON = binned_data.LON(idx);
matched_data.LAT = binned_data.LAT(idx);
matched_data.PRES = binned_data.PRES(idx);
matched_data.TEMP_FLOAT = binned_data.TEMP_FLOAT(idx);
matched_data.TEMP_GLODAP = binned_data.TEMP_GLODAP(idx);
matched_data.TEMP = mean([matched_data.TEMP_FLOAT matched_data.TEMP_GLODAP],2);
% fit delta against [O2]
mod = fitlm(matched_data.OXY_FLOAT,matched_data.OXY_DELTA,'Intercept',false);
slp = mod.Coefficients.Estimate(1);
int = 0; % mod.Coefficients.Estimate(1);
% apply linear correction to float O2
matched_data.OXY_FLOAT_CORR = matched_data.OXY_FLOAT + slp.*matched_data.OXY_FLOAT + int;
% re-calculate delta
matched_data.OXY_DELTA_CORR = matched_data.OXY_GLODAP-matched_data.OXY_FLOAT_CORR;

%% plot uncorrected float vs. glodap residuals
figure; hold on;
scatter(matched_data.OXY_GLODAP,matched_data.OXY_FLOAT,20,'.');
plot([0,450],[0 450],'k--');
ylabel('Binned BGC Argo Oxygen Data (\mumol kg^{-1})');
xlabel('Binned GLODAP Oxygen Data (\mumol kg^{-1})');
matched_data.err_md = median(matched_data.OXY_DELTA);
matched_data.std = std(matched_data.OXY_DELTA);
matched_data.rmse = sqrt(mean((matched_data.OXY_DELTA).^2));
matched_data.r2 = corr(matched_data.OXY_GLODAP,matched_data.OXY_FLOAT);
matched_data.slope = polyfit(matched_data.OXY_GLODAP,matched_data.OXY_FLOAT,1);
matched_data.slope = matched_data.slope(1);
text(300,120,['Med. Err. = ' num2str(round(matched_data.err_md,2))],'fontsize',12);
text(300,90,['RMSE = ' num2str(round(matched_data.rmse,1))],'fontsize',12);
text(300,60,['R^{2} = ' num2str(round(matched_data.r2,2))],'fontsize',12);
% save figure
exportgraphics(gcf,'Figures/Data/delta_vs_float_300_uncorr.jpg');
close

%% plot corrected float vs. glodap residuals
figure; hold on;
scatter(matched_data.OXY_GLODAP,matched_data.OXY_FLOAT_CORR,20,'.');
plot([0,450],[0 450],'k--');
ylabel('Corrected, Binned BGC Argo Oxygen Data (\mumol kg^{-1})');
xlabel('Binned GLODAP Oxygen Data (\mumol kg^{-1})');
matched_data.err_md_corr = median(matched_data.OXY_DELTA_CORR);
matched_data.std_corr = std(matched_data.OXY_DELTA_CORR);
matched_data.rmse_corr = sqrt(mean((matched_data.OXY_DELTA_CORR).^2));
matched_data.r2_corr = corr(matched_data.OXY_GLODAP,matched_data.OXY_FLOAT_CORR);
matched_data.slope_corr = polyfit(matched_data.OXY_GLODAP,matched_data.OXY_FLOAT_CORR,1);
matched_data.slope_corr = matched_data.slope_corr(1);
text(300,120,['Med. Err. = ' num2str(round(matched_data.err_md_corr,2))],'fontsize',12);
text(300,90,['RMSE = ' num2str(round(matched_data.rmse_corr,1))],'fontsize',12);
text(300,60,['R^{2} = ' num2str(round(matched_data.r2_corr,2))],'fontsize',12);
% save figure
exportgraphics(gcf,'Figures/Data/delta_vs_float_300_corr.jpg');
close

%% save matched float and glodap data
save(['matched_data_' datestr(date) '.mat'],'matched_data','-v7.3');
% clean up
clear idx matched_data

%% save correction factors
save(['float_corr_' datestr(date) '.mat'],'slp','int');

%% average float and glodap bins that overlap (and add correction to float oxy)
% add correction to float data
binned_data.OXY_FLOAT = ...
    binned_data.OXY_FLOAT + slp.*binned_data.OXY_FLOAT + int;
% average float and glodap bins
binned_data.OXY = ...
    mean(cat(5,binned_data.OXY_FLOAT,binned_data.OXY_GLODAP),5,'omitnan');
binned_data.TEMP = ...
    mean(cat(5,binned_data.TEMP_FLOAT,binned_data.TEMP_GLODAP),5,'omitnan');
binned_data.SAL = ...
    mean(cat(5,binned_data.SAL_FLOAT,binned_data.SAL_GLODAP),5,'omitnan');

%% save binned data
save(['binned_data_' datestr(date) '.mat'],'binned_data','-v7.3');

%% Display binned data distribution

% Index gridded data
f_idx = ~isnan(binned_data.OXY_FLOAT);
f_idx = squeeze(any(f_idx,3));
g_idx = ~isnan(binned_data.OXY_GLODAP);
g_idx = squeeze(any(g_idx,3));
time_temp = squeeze(binned_data.TIME(:,:,1,:));
lat_temp = squeeze(binned_data.LAT(:,:,1,:));
lon_temp = squeeze(binned_data.LON(:,:,1,:));

% By time
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
year_temp = 2004:2023;
edges=datenum([year_temp' ones(length(year_temp),1) ones(length(year_temp),1)]);
histogram(time_temp(f_idx),edges,'FaceColor','r');
histogram(time_temp(g_idx),edges,'FaceColor','b');
legend({'Floats' 'GLODAP'},'location','northwest');
datetick('x'); xlim([datenum([2003 1 1]) datenum([2024 1 1])]);
ylabel('Profiles within each year');
if ~exist('Figures','dir'); mkdir('Figures'); end
if ~exist('Figures/Data','dir'); mkdir('Figures/Data'); end
exportgraphics(gcf,'Figures/Data/data_by_year_binned.jpg');
clear edges
close

% By latitude
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
edges=-90:5:90;
histogram(lat_temp(f_idx),edges,'FaceColor','r');
histogram(lat_temp(g_idx),edges,'FaceColor','b');
legend({'Floats' 'GLODAP'})
ylabel('Profiles within each latitude range');
exportgraphics(gcf,'Figures/Data/data_by_latitude_binned.jpg');
close

% By longitude
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
edges=-180:10:180;
histogram(lon_temp(f_idx),edges,'FaceColor','r');
histogram(lon_temp(g_idx),edges,'FaceColor','b');
legend({'Floats' 'GLODAP'})
ylabel('Profiles within each longitude range');
exportgraphics(gcf,'Figures/Data/data_by_longitude_binned.jpg');
close

% clean up
clear binned_data f_idx g_idx edges lat_temp lon_temp time_temp year_temp slp c ans

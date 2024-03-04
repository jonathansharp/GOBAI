% display_o2_data
%
% DESCRIPTION:
% This function is used to create and save some descriptive figures of the
% data distribution for GOBAI-O2 algorithm training.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/12/2023

function display_o2_data(float_file_ext,snap_date,glodap_year)

%% load interpolated float and glodap data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_float_o2_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date');
load(['Data/processed_glodap_o2_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');

%% plot for sanity
clrs=colororder;
% float (depth average)
zi = ([2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';
oxy_float_mean = nan(size(zi));
oxy_float_std = nan(size(zi));
for z = 1:length(zi)
    idx_float = float_data.OXY_PRES == zi(z);
    oxy_float_mean(z) = mean(float_data.OXY(idx_float));
    oxy_float_std(z) = std(float_data.OXY(idx_float));
end
figure; hold on;
plot(oxy_float_mean,zi);
fill([oxy_float_mean+oxy_float_std;flipud(oxy_float_mean-oxy_float_std)],...
    [zi;flipud(zi)],clrs(1,:),'FaceAlpha',0.25,'LineStyle','none');
set(gca,'YDir','reverse');
hold off;
exportgraphics(gcf,['Figures/Data/mean_float_profile_' file_date float_file_ext '.png']);
close
% glodap (depth average)
zi = ([2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';
oxy_glodap_mean = nan(size(zi));
oxy_glodap_std = nan(size(zi));
for z = 1:length(zi)
    idx_glodap = glodap_data.OXY_PRES == zi(z);
    oxy_glodap_mean(z) = mean(glodap_data.OXY(idx_glodap),'omitnan');
    oxy_glodap_std(z) = std(glodap_data.OXY(idx_glodap),'omitnan');
end
figure; hold on;
plot(oxy_glodap_mean,zi);
fill([oxy_glodap_mean+oxy_glodap_std;flipud(oxy_glodap_mean-oxy_glodap_std)],...
    [zi;flipud(zi)],clrs(1,:),'FaceAlpha',0.25,'LineStyle','none');
set(gca,'YDir','reverse');
hold off;
exportgraphics(gcf,['Figures/Data/mean_glodap_profile_' num2str(glodap_year) '.png']);
close
% glodap (individual profiles)
figure; hold on;
profs = unique(glodap_data.OXY_ID);
for z = 1:length(profs)
    idx = glodap_data.OXY_ID == profs(z);
    plot(glodap_data.OXY(idx),glodap_data.OXY_PRES(idx));
end
set(gca,'YDir','reverse');
hold off;
exportgraphics(gcf,['Figures/Data/all_float_profiles_' file_date float_file_ext '.png']);
close
% clean up
clear clrs idx idx_float idx_glodap oxy_float_mean oxy_float_std 
clear oxy_glodap_mean oxy_glodap_std profs z zi

%% Display data distribution

% float
f_profiles = unique(float_data.OXY_PROF_ID);
f_idx = nan(length(f_profiles),1);
for n = 1:length(f_profiles)
    f_idx(n) = find(f_profiles(n)==float_data.OXY_PROF_ID,1);
end
% glodap
g_profiles = unique(glodap_data.OXY_ID);
g_idx = nan(length(g_profiles),1);
for n = 1:length(g_profiles)
    g_idx(n) = find(g_profiles(n)==glodap_data.OXY_ID,1);
end
clear n

% By time
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
min_year = datevec(min(float_data.OXY_TIME));
min_year = min_year(1);
max_year = datevec(max(float_data.OXY_TIME));
max_year = max_year(1);
year_temp = (min_year:max_year)';
edges=datenum([year_temp ones(length(year_temp),1) ones(length(year_temp),1)]);
histogram(float_data.OXY_TIME(f_idx),edges,'FaceColor','r');
histogram(glodap_data.OXY_TIME(g_idx),edges,'FaceColor','b');
legend({'Floats' 'GLODAP'},'location','northwest');
datetick('x'); xlim([datenum([2003 1 1]) datenum([2024 1 1])]);
ylabel('Profiles within each year');
if ~exist([pwd '/Figures'],'dir'); mkdir('Figures'); end
if ~exist([pwd '/Figures/Data'],'dir'); mkdir('Figures/Data'); end
exportgraphics(gcf,['Figures/Data/data_by_year_' file_date float_file_ext '.png']);
clear edges year_temp
close

% By latitude
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
edges=-90:5:90;
histogram(float_data.OXY_LAT(f_idx),edges,'FaceColor','r');
histogram(glodap_data.OXY_LAT(g_idx),edges,'FaceColor','b');
legend({'Floats' 'GLODAP'})
ylabel('Profiles within each latitude range');
if ~exist([pwd '/Figures/Data'],'dir'); mkdir('Figures/Data'); end
exportgraphics(gcf,['Figures/Data/data_by_latitude_' file_date float_file_ext '.png']);
clear edges
close

% By longitude
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
edges=-180:5:180;
histogram(float_data.OXY_LON(f_idx),edges,'FaceColor','r');
histogram(glodap_data.OXY_LON(g_idx),edges,'FaceColor','b');
legend({'Floats' 'GLODAP'})
ylabel('Profiles within each longitude range');
if ~exist([pwd '/Figures/Data'],'dir'); mkdir('Figures/Data'); end
exportgraphics(gcf,['Figures/Data/data_by_longitude_' file_date float_file_ext '.png']);
clear edges
close

end
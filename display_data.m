% display_data
%
% DESCRIPTION:
% This function is used to create and save some descriptive figures of the
% data distribution for GOBAI algorithm training.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/12/2023

function display_data(param_props,float_file_ext,file_date,glodap_year)

%% load interpolated float and glodap data
load([param_props.p1 '/Data/processed_float_' param_props.p2 '_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date');
load([param_props.p1 '/Data/processed_glodap_' param_props.p2 '_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');

%% plot for sanity
clrs=colororder;
% float (depth average)
zi = ([2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';
float_mean = nan(size(zi));
float_std = nan(size(zi));
for z = 1:length(zi)
    idx_float = float_data.PRES == zi(z);
    float_mean(z) = mean(float_data.(param_props.p5)(idx_float));
    float_std(z) = std(float_data.(param_props.p5)(idx_float));
end
figure; hold on;
title('Average Float Profile');
plot(float_mean,zi,'linewidth',3);
scatter(float_mean,zi,'.k');
fill([float_mean+float_std;flipud(float_mean-float_std)],...
    [zi;flipud(zi)],clrs(1,:),'FaceAlpha',0.25,'LineStyle','none');
set(gca,'YDir','reverse');
ylabel('Depth (dbar)');
xlabel('[O_{2}] (\mumol kg^{-1})');
hold off;
exportgraphics(gcf,[param_props.p1 '/Figures/Data/mean_float_profile_' file_date float_file_ext '.png']);
close
% floats (individual profiles)
figure; hold on;
title('5% of Float Profiles');
profs = unique(float_data.PROF_ID);
for z = 1:length(profs)/20
    idx = float_data.PROF_ID == profs(z*20);
    plot(float_data.(param_props.p5)(idx),float_data.PRES(idx));
end
set(gca,'YDir','reverse');
ylabel('Depth (dbar)');
xlabel('[O_{2}] (\mumol kg^{-1})');
hold off;
exportgraphics(gcf,[param_props.p1 '/Figures/Data/all_float_profiles_' file_date float_file_ext '.png']);
close
% clean up
clear idx idx_float float_mean float_std profs z zi
% glodap (depth average)
zi = ([2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';
glodap_mean = nan(size(zi));
glodap_std = nan(size(zi));
for z = 1:length(zi)
    idx_glodap = glodap_data.PRES == zi(z);
    glodap_mean(z) = mean(glodap_data.(param_props.p5)(idx_glodap),'omitnan');
    glodap_std(z) = std(glodap_data.(param_props.p5)(idx_glodap),'omitnan');
end
figure; hold on;
plot(glodap_mean,zi,'linewidth',3);
scatter(glodap_mean,zi,'.k');
fill([glodap_mean+glodap_std;flipud(glodap_mean-glodap_std)],...
    [zi;flipud(zi)],clrs(1,:),'FaceAlpha',0.25,'LineStyle','none');
set(gca,'YDir','reverse');
hold off;
exportgraphics(gcf,[param_props.p1 '/Figures/Data/mean_glodap_profile_' num2str(glodap_year) '.png']);
close
% glodap (individual profiles)
figure; hold on;
profs = unique(glodap_data.ID);
for z = 1:length(profs)
    idx = glodap_data.ID == profs(z);
    plot(glodap_data.(param_props.p5)(idx),glodap_data.PRES(idx));
end
set(gca,'YDir','reverse');
hold off;
exportgraphics(gcf,[param_props.p1 '/Figures/Data/all_glodap_profiles_' num2str(glodap_year) '.png']);
close
% clean up
clear clrs idx idx_glodap glodap_mean glodap_std profs z zi

%% Display data distribution

% float
f_profiles = unique(float_data.PROF_ID);
f_idx = nan(length(f_profiles),1);
for n = 1:length(f_profiles)
    f_idx(n) = find(f_profiles(n)==float_data.PROF_ID,1);
end
% glodap
g_profiles = unique(glodap_data.ID);
g_idx = nan(length(g_profiles),1);
for n = 1:length(g_profiles)
    g_idx(n) = find(g_profiles(n)==glodap_data.ID,1);
end
clear n

% By time
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
min_year = datevec(min([float_data.TIME;glodap_data.TIME]));
min_year = min_year(1);
max_year = datevec(max([float_data.TIME;glodap_data.TIME]));
max_year = max_year(1);
year_temp = (min_year:max_year+1)';
edges=datenum([year_temp ones(length(year_temp),1) ones(length(year_temp),1)]);
histogram(float_data.TIME(f_idx),edges,'FaceColor','r');
histogram(glodap_data.TIME(g_idx),edges,'FaceColor','b');
legend({'Floats' 'GLODAP'},'location','northwest');
datetick('x'); xlim([datenum([2003 1 1]) datenum([2024 1 1])]);
ylabel('Profiles within each year');
if ~exist([pwd '/' param_props.p1 '/Figures'],'dir'); mkdir([param_props.p1 '/Figures']); end
if ~exist([pwd '/' param_props.p1 '/Figures/Data'],'dir'); mkdir([param_props.p1 '/Figures/Data']); end
exportgraphics(gcf,[param_props.p1 '/Figures/Data/data_by_year_' file_date float_file_ext '.png']);
clear edges year_temp
close

% By latitude
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
edges=-90:5:90;
histogram(float_data.LAT(f_idx),edges,'FaceColor','r');
histogram(glodap_data.LAT(g_idx),edges,'FaceColor','b');
legend({'Floats' 'GLODAP'})
ylabel('Profiles within each latitude range');
if ~exist([pwd '/' param_props.p1 '/Figures/Data'],'dir'); mkdir([param_props.p1 '/Figures/Data']); end
exportgraphics(gcf,[param_props.p1 '/Figures/Data/data_by_latitude_' file_date float_file_ext '.png']);
clear edges
close

% By longitude
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
edges=-180:5:180;
histogram(float_data.LON(f_idx),edges,'FaceColor','r');
histogram(glodap_data.LON(g_idx),edges,'FaceColor','b');
legend({'Floats' 'GLODAP'})
ylabel('Profiles within each longitude range');
if ~exist([pwd '/' param_props.p1 '/Figures/Data'],'dir'); mkdir([param_props.p1 '/Figures/Data']); end
exportgraphics(gcf,[param_props.p1 '/Figures/Data/data_by_longitude_' file_date float_file_ext '.png']);
clear edges
close

end
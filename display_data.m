% display_data
%
% DESCRIPTION:
% This function is used to create and save some descriptive figures of the
% data distribution for GOBAI algorithm training.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 7/28/2025

function display_data(param_props,float_file_ext,glodap_year,start_year,snap_date,flt,gld,ctd)

%% load interpolated float and glodap data
file_date = datestr(datenum(floor(snap_date/1e2),...
    mod(snap_date,1e2),1),'mmm-yyyy');
if flt == 1; load([param_props.dir_name '/Data/processed_float_' ...
    param_props.file_name '_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date'); end
if gld == 1; load([param_props.dir_name '/Data/processed_glodap_' ...
    param_props.file_name '_data_' num2str(glodap_year) '.mat'],...
    'glodap_data'); end
if ctd == 1; load([param_props.dir_name '/Data/processed_wod_ctd_' ...
    param_props.file_name '_data_' num2str(glodap_year) '.mat'],...
    'wod_data'); end

%% plot for sanity
% clrs=colororder;
% % float (depth average)
% zi = ([2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';
% float_mean = nan(size(zi));
% float_std = nan(size(zi));
% for z = 1:length(zi)
%     idx_float = float_data.PRES == zi(z);
%     float_mean(z) = mean(float_data.(param_props.temp_name)(idx_float));
%     float_std(z) = std(float_data.(param_props.temp_name)(idx_float));
% end
% figure; hold on;
% title('Average Float Profile');
% plot(float_mean,zi,'linewidth',3);
% scatter(float_mean,zi,'.k');
% fill([float_mean+float_std;flipud(float_mean-float_std)],...
%     [zi;flipud(zi)],clrs(1,:),'FaceAlpha',0.25,'LineStyle','none');
% set(gca,'YDir','reverse');
% ylabel('Depth (dbar)');
% xlabel('[O_{2}] (\mumol kg^{-1})');
% hold off;
% if ~exist([pwd '/' param_props.dir_name '/Figures/Data'],'dir')
%     mkdir([param_props.dir_name '/Figures/Data']); end
% export_fig(gcf,[param_props.dir_name '/Figures/Data/mean_float_profile_' ...
%     file_date float_file_ext '.png'],'-transparent');
% close
% % floats (individual profiles)
% figure; hold on;
% title('5% of Float Profiles');
% profs = unique(float_data.PROF_ID);
% for z = 1:length(profs)/20
%     idx = float_data.PROF_ID == profs(z*20);
%     plot(float_data.(param_props.temp_name)(idx),float_data.PRES(idx));
% end
% set(gca,'YDir','reverse');
% ylabel('Depth (dbar)');
% xlabel('[O_{2}] (\mumol kg^{-1})');
% hold off;
% if ~exist([pwd '/' param_props.dir_name '/Figures/Data'],'dir')
%     mkdir([param_props.dir_name '/Figures/Data']); end
% export_fig(gcf,[param_props.dir_name '/Figures/Data/all_float_profiles_' ...
%     file_date float_file_ext '.png'],'-transparent');
% close
% % clean up
% clear idx idx_float float_mean float_std profs z zi
% % glodap (depth average)
% zi = ([2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';
% glodap_mean = nan(size(zi));
% glodap_std = nan(size(zi));
% for z = 1:length(zi)
%     idx_glodap = glodap_data.PRES == zi(z);
%     glodap_mean(z) = mean(glodap_data.(param_props.temp_name)(idx_glodap),'omitnan');
%     glodap_std(z) = std(glodap_data.(param_props.temp_name)(idx_glodap),'omitnan');
% end
% figure; hold on;
% plot(glodap_mean,zi,'linewidth',3);
% scatter(glodap_mean,zi,'.k');
% fill([glodap_mean+glodap_std;flipud(glodap_mean-glodap_std)],...
%     [zi;flipud(zi)],clrs(1,:),'FaceAlpha',0.25,'LineStyle','none');
% set(gca,'YDir','reverse');
% hold off;
% if ~exist([pwd '/' param_props.dir_name '/Figures/Data'],'dir')
%     mkdir([param_props.dir_name '/Figures/Data']); end
% export_fig(gcf,[param_props.dir_name '/Figures/Data/mean_glodap_profile_' ...
%     num2str(glodap_year) '.png'],'-transparent');
% close
% % glodap (individual profiles)
% figure; hold on;
% profs = unique(glodap_data.ID);
% for z = 1:length(profs)
%     idx = glodap_data.ID == profs(z);
%     plot(glodap_data.(param_props.temp_name)(idx),glodap_data.PRES(idx));
% end
% set(gca,'YDir','reverse');
% hold off;
% if ~exist([pwd '/' param_props.dir_name '/Figures/Data'],'dir')
%     mkdir([param_props.dir_name '/Figures/Data']); end
% export_fig(gcf,[param_props.dir_name '/Figures/Data/all_glodap_profiles_' ...
%     num2str(glodap_year) '.png'],'-transparent');
% close
% % clean up
% clear clrs idx idx_glodap glodap_mean glodap_std profs z zi
% 
% %% Display data distribution
% 
% if ~isfile([param_props.dir_name '/Figures/Data/data_by_year_' file_date float_file_ext '.png']) && ...
%    ~isfile([param_props.dir_name '/Figures/Data/data_by_latitude_' file_date float_file_ext '.png']) && ...
%    ~isfile([param_props.dir_name '/Figures/Data/data_by_longitude_' file_date float_file_ext '.png'])

% float
if flt == 1; [f_profiles,f_idx] = unique(float_data.PROF_ID); end
% glodap
if gld == 1; [g_profiles,g_idx] = unique(glodap_data.ID); end
% wod
if ctd == 1; [w_profiles,w_idx] = unique(wod_data.ID); end

% By time
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
min_year = datevec(min([float_data.TIME;glodap_data.TIME;wod_data.TIME]));
min_year = min_year(1);
max_year = datevec(max([float_data.TIME;glodap_data.TIME;wod_data.TIME]));
max_year = max_year(1);
year_temp = (min_year:max_year+1)';
edges = datenum([year_temp ones(length(year_temp),1) ones(length(year_temp),1)]);
if gld == 1; g_counts = histcounts(glodap_data.TIME(g_idx),edges); else g_counts = zeros(size(edges)); end
if flt == 1; f_counts = histcounts(float_data.TIME(f_idx),edges); else f_counts = zeros(size(edges)); end
if ctd == 1; w_counts = histcounts(wod_data.TIME(w_idx),edges); else w_counts = zeros(size(edges)); end
bar(edges(1:end-1),[g_counts',w_counts',f_counts'],'stacked');
legend({'GLODAP' 'CTD' 'Floats'},'location','northwest');
datetick('x'); xlim([datenum([start_year-1 1 1]) max(edges)]);
ylabel('Profiles within each year');
if ~exist([pwd '/' param_props.dir_name '/Figures'],'dir'); mkdir([param_props.dir_name '/Figures']); end
if ~exist([pwd '/' param_props.dir_name '/Figures/Data'],'dir'); mkdir([param_props.dir_name '/Figures/Data']); end
export_fig(gcf,[param_props.dir_name '/Figures/Data/data_by_year_' file_date ...
    float_file_ext '.png'],'-transparent');
clear edges year_temp
close

% By latitude
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
edges=-90:5:90;
if gld == 1; g_counts = histcounts(glodap_data.LAT(g_idx),edges); else g_counts = zeros(size(edges)); end
if flt == 1; f_counts = histcounts(float_data.LAT(f_idx),edges); else f_counts = zeros(size(edges)); end
if ctd == 1; w_counts = histcounts(wod_data.LAT(w_idx),edges); else w_counts = zeros(size(edges)); end
bar(edges(1:end-1),[g_counts',w_counts',f_counts'],'stacked');
legend({'CTD' 'Floats' 'GLODAP'});
ylabel('Profiles within each latitude range');
if ~exist([pwd '/' param_props.dir_name '/Figures/Data'],'dir'); mkdir([param_props.dir_name '/Figures/Data']); end
export_fig(gcf,[param_props.dir_name '/Figures/Data/data_by_latitude_' ...
    file_date float_file_ext '.png'],'-transparent');
clear edges
close

% By longitude
figure; hold on;
set(gca,'fontsize',20);
set(gcf,'units','inches','position',[0 5 10 10]);
edges=-180:5:180;
if gld == 1; g_counts = histcounts(glodap_data.LON(g_idx),edges); else g_counts = zeros(size(edges)); end
if flt == 1; f_counts = histcounts(float_data.LON(f_idx),edges); else f_counts = zeros(size(edges)); end
if ctd == 1; w_counts = histcounts(wod_data.LON(w_idx),edges); else w_counts = zeros(size(edges)); end
bar(edges(1:end-1),[g_counts',w_counts',f_counts'],'stacked');
legend({'CTD' 'Floats' 'GLODAP'});
ylabel('Profiles within each longitude range');
if ~exist([pwd '/' param_props.dir_name '/Figures/Data'],'dir'); mkdir([param_props.dir_name '/Figures/Data']); end
export_fig(gcf,[param_props.dir_name '/Figures/Data/data_by_longitude_' ...
    file_date float_file_ext '.png'],'-transparent');
clear edges
close

%% Maps by decade
dec_start = floor(start_year/10)*10:10:2020;
dec_end = floor(start_year/10)*10+10:10:2030;
clrs = colororder;
figure('visible','on');
set(gcf,'units','inches','position',[0 5 20 10]);
tiledlayout(2,2,'Padding','none','TileSpacing','compact');
for d = 1:length(dec_start)
    nexttile; hold on;
    % index coordinates by decade
    g_idx_y = glodap_data.YEAR >= dec_start(d) & glodap_data.YEAR < dec_end(d);
    lon_g = glodap_data.LON; lon_g(~g_idx_y) = NaN;
    lat_g = glodap_data.LAT; lat_g(~g_idx_y) = NaN;
    f_idx_y = float_data.YEAR >= dec_start(d) & float_data.YEAR < dec_end(d);
    lon_f = float_data.LON; lon_f(~f_idx_y) = NaN;
    lat_f = float_data.LAT; lat_f(~f_idx_y) = NaN;
    w_idx_y = wod_data.YEAR >= dec_start(d) & wod_data.YEAR < dec_end(d);
    lon_w = wod_data.LON; lon_w(~w_idx_y) = NaN;
    lat_w = wod_data.LAT; lat_w(~w_idx_y) = NaN;
    % set map projection
    m_proj('robinson','lon',[20 380]);
    % convert longitude to proper format
    temp_lon_g = convert_lon(lon_g(g_idx),'format','0-360');
    temp_lon_g(temp_lon_g < 20) = temp_lon_g(temp_lon_g < 20) + 360;
    temp_lon_f = convert_lon(lon_f(f_idx),'format','0-360');
    temp_lon_f(temp_lon_f < 20) = temp_lon_f(temp_lon_f < 20) + 360;
    temp_lon_w = convert_lon(lon_w(w_idx),'format','0-360');
    temp_lon_w(temp_lon_w < 20) = temp_lon_w(temp_lon_w < 20) + 360;
    % scatter data
    m_scatter(temp_lon_f,lat_f(f_idx),5,clrs(1,:),'filled');
    m_scatter(temp_lon_w,lat_w(w_idx),5,clrs(2,:),'filled');
    m_scatter(temp_lon_g,lat_g(g_idx),5,clrs(3,:),'filled');
    % set properties
    title([num2str(dec_start(d)) ' to ' num2str(dec_end(d))]);
    m_coast('patch',rgb('grey'));
    m_grid('linestyle','none','xticklabels',[],...
        'yticklabels',[],'ytick',-90:30:90);
    % legend({'Float Data' 'CTD Data' 'GLODAP Data'});
end
export_fig(gcf,[param_props.dir_name '/Figures/Data/Mapped_' ...
    param_props.dir_name '_' file_date float_file_ext '.png'],'-transparent');
close

%% Map for full period
clrs = colororder;
figure('visible','on'); hold on;
set(gcf,'units','inches','position',[0 5 20 10]);
% set map projection
m_proj('robinson','lon',[20 380]);
% convert longitude to proper format
temp_lon_g = convert_lon(glodap_data.LON,'format','0-360');
temp_lon_g(temp_lon_g < 20) = temp_lon_g(temp_lon_g < 20) + 360;
temp_lon_f = convert_lon(float_data.LON,'format','0-360');
temp_lon_f(temp_lon_f < 20) = temp_lon_f(temp_lon_f < 20) + 360;
temp_lon_w = convert_lon(woa_data.LON,'format','0-360');
temp_lon_w(temp_lon_w < 20) = temp_lon_w(temp_lon_w < 20) + 360;
% scatter data
m_scatter(temp_lon_f,float_data.LAT,10,clrs(1,:),'filled');
m_scatter(temp_lon_w,woa_data.LAT,10,clrs(2,:),'filled');
m_scatter(temp_lon_g,glodap_data.LAT,10,clrs(3,:),'filled');
% set properties
title([param_props.label ' Data Distribution']);
m_coast('patch',rgb('grey'));
m_grid('linestyle','none','xticklabels',[],...
    'yticklabels',[],'ytick',-90:30:90);
legend({'Float Data' 'CTD Data' 'GLODAP Data'},'Location','northwest');
set(gca,'FontSize',24);
export_fig(gcf,[param_props.dir_name '/Figures/Data/All_Mapped_' ...
    param_props.dir_name '_' file_date float_file_ext '.png'],'-transparent');
close

end

% end
% plot Figure 1 for Sharp et al. (in prep) GOBAI-O2 and GOBAI-NO3 Paper

clrs = {'#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3'};
tiledlayout(5,2);
set(gcf,'units','inches','position',[0 5 20 10],'visible','off');
dot_size = 2;

labels = {'O_{2}' 'NO_{3}'};
dirs = {'O2' 'NO3'};
flns = {'o2' 'no3'};
tiles = [1,7;2,8];

for l = 1:length(labels)

nexttile(tiles(l,1),[3 1]); hold on;

%% load interpolated float and glodap data for o2
if strcmp(flns{l},'o2')
    load([dirs{l} '/Data/processed_float_' flns{l} ...
        '_data_adjusted_Jun-2026_D_A.mat'],'float_data_adjusted','file_date');
elseif strcmp(flns{l},'no3')
    load([dirs{l} '/Data/processed_float_' flns{l} ...
        '_data_adjusted_Jun-2026_D_A.mat'],'float_data_adjusted','file_date');
end
load([dirs{l} '/Data/processed_glodap_' flns{l} '_data_adjusted_v3.mat'],'glodap_data');
load([dirs{l} '/Data/processed_wod_' flns{l} '_data_adjusted_2026.mat'],'wod_data');

%% define indices
% float
[~,f_idx] = unique(float_data_adjusted.PROF_ID);
% glodap
[~,g_idx] = unique(glodap_data.ID);
% osd
id_temp = wod_data.ID;
id_temp(wod_data.TYPE==2) = NaN;
[~,o_idx] = unique(id_temp,'TreatMissingAsDistinct',false);
o_idx(end) = [];
% ctd
id_temp = wod_data.ID;
id_temp(wod_data.TYPE==1) = NaN;
[~,c_idx] = unique(id_temp,'TreatMissingAsDistinct',false);
c_idx(end) = [];

%% Map for full period
% set map projection
m_proj('robinson','lon',[20 380]);
% float
temp_lon_f = convert_lon(float_data_adjusted.LON(f_idx),'format','0-360');
temp_lon_f(temp_lon_f < 20) = temp_lon_f(temp_lon_f < 20) + 360;
m_scatter(temp_lon_f,float_data_adjusted.LAT(f_idx),dot_size,hex2rgb(clrs{1}),'filled');
% osd
temp_lon_o = convert_lon(wod_data.LON(o_idx),'format','0-360');
temp_lon_o(temp_lon_o < 20) = temp_lon_o(temp_lon_o < 20) + 360;
m_scatter(temp_lon_o,wod_data.LAT(o_idx),dot_size,hex2rgb(clrs{2}),'filled');
% ctd
temp_lon_c = convert_lon(wod_data.LON(c_idx),'format','0-360');
temp_lon_c(temp_lon_c < 20) = temp_lon_c(temp_lon_c < 20) + 360;
m_scatter(temp_lon_c,wod_data.LAT(c_idx),dot_size,hex2rgb(clrs{3}),'filled');
% gld
temp_lon_g = convert_lon(glodap_data.LON(g_idx),'format','0-360');
temp_lon_g(temp_lon_g < 20) = temp_lon_g(temp_lon_g < 20) + 360;
m_scatter(temp_lon_g,glodap_data.LAT(g_idx),dot_size,hex2rgb(clrs{4}),'filled');

% set properties
title([labels{l} ' Data Distribution'],'FontSize',20);
m_coast('patch',rgb('grey'));
m_grid('linestyle','none','xticklabels',[],...
    'yticklabels',[],'ytick',-90:30:90);
if strcmp(flns{l},'no3')
    map_legend = {'' 'Argo Float Profiles' 'WOD Profiles (OSD)' ...
        'WOD Profiles (CTD)' 'GLODAP Ship Profiles'};
    hl = legend(map_legend,'Location','northwest','FontSize',16);
    hl.Units = 'normalized';
    hl.Position = [0.425 0.75 0.15 0.15];
end

% next tile
nexttile(tiles(l,2),[2 1]); hold on;
set(gca,'FontSize',16);

% load data
if strcmp(flns{l},'o2')
    load([dirs{l} '/Data/processed_all_' flns{l} ...
        '_data_fgoc_Jun-2026_D_A.mat'],'all_data');
elseif strcmp(flns{l},'no3')
    load([dirs{l} '/Data/processed_all_' flns{l} ...
        '_data_fgo_Jun-2026_D.mat'],'all_data');
end

% index to unique profiles
[~,prof_idx] = unique(all_data.id);
vars = fieldnames(all_data);
for v = 1:length(vars)
    if ndims(all_data.(vars{v})) < 3
        all_data.(vars{v}) = all_data.(vars{v})(prof_idx);
    end
end

% index to each dataset type
y_osd = all_data.year(all_data.type==3);
y_ctd = all_data.year(all_data.type==4);
y_gld = all_data.year(all_data.type==2);
y_flt = all_data.year(all_data.type==1);

% count profiles in each year
counts_osd = histc(y_osd,1993:2026);
counts_ctd = histc(y_ctd,1993:2026);
counts_gld = histc(y_gld,1993:2026);
counts_flt = histc(y_flt,1993:2026);

% plot histogram
b=bar(1993:2026,[counts_osd counts_ctd counts_gld counts_flt],'stacked');
clr_indices = [2;3;4;1];
for k=1:numel(b)
    b(k).FaceColor = hex2rgb(clrs{clr_indices(k)});
end
if strcmp(flns{l},'o2')
    legend({'OSD' 'CTD' 'GLODAP' 'Argo Float'},...
        'NumColumns',4,'Location','northwest','FontSize',16);
elseif strcmp(flns{l},'no3')
    legend({'OSD' '' 'GLODAP' 'Argo Float'},...
        'NumColumns',4,'Location','northwest','FontSize',16);
end
ylim([0 30000]);
ylabel('Number of Profiles','FontSize',20);

end

export_fig(gcf,'Figures/Figure1_data_distribution.png','-transparent');
close

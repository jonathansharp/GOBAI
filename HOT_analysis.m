
inf = ncinfo('Data/HOT_Data.nc');
num_clusters = 15;

hot.lon = -157.87;
hot.lat = 21.31;

% read data
for v = 1:length(inf.Variables)
    hot.(inf.Variables(v).Name) = ...
        ncread('Data/HOT_Data.nc',inf.Variables(v).Name);
end
% index
idx = hot.mdate < 0 | hot.boxy < 0;
for v = 1:length(inf.Variables)
    if isnumeric(hot.(inf.Variables(v).Name))
        hot.(inf.Variables(v).Name)(idx) = []; 
    end
end
% process time
month = floor(hot.mdate/10000);
day = floor(mod(hot.mdate,10000)/100);
year = mod(hot.mdate,100);
year(year>50) = year(year>50)+1900;
year(year<50) = year(year<50)+2000;
time = datenum(year,month,day);
time0 = datenum(year,1,1);
doy = time-time0;
% calculate absolute salinity and conservative temperature
hot.sal_abs = gsw_SA_from_SP(hot.csal,hot.press,hot.lon,hot.lat);
hot.temp_cns = gsw_CT_from_t(hot.sal_abs,hot.temp,hot.press);
% calculate o2 with gobai algorithm
hot.gobai_o2 = gobai_nn(num_clusters,'temp_cns',hot.temp_cns,'sal_abs',hot.sal_abs,...
    'sigma',hot.sigma,'lon',-157.87,'lat',21.31,'pres',hot.press,...
    'doy',doy,'year',year);
% define gobai file name
path_ext = ['GOBAI/RFROM/FFNN/c' num2str(num_clusters) ...
        '_Jul-2025_D/train80_val10_test10/fgw/'];
gobai_alg_dir = [fpaths.param_path path_ext];
filename = [gobai_alg_dir 'gobai-' param_props.file_name '.nc'];
gobai.time = datenum(1950,1,1) + ncread(filename,'time');
gobai.lon = ncread(filename,'lon');
gobai.lat = ncread(filename,'lat');
gobai.pres = ncread(filename,'pres');
[~,idx_lon] = min(abs(gobai.lon-(hot.lon+360)));
[~,idx_lat] = min(abs(gobai.lat-hot.lat));
[~,idx_pres10] = min(abs(gobai.pres-10));
[~,idx_pres100] = min(abs(gobai.pres-100));
gobai.o2_10 = squeeze(ncread(filename,'o2',[idx_lon idx_lat idx_pres10 1],[1 1 1 Inf]));
gobai.o2_100 = squeeze(ncread(filename,'o2',[idx_lon idx_lat idx_pres100 1],[1 1 1 Inf]));
% plot HOT vs. GOBAI (10 dbar)
idx10 = hot.press >= 0 & hot.press <= 20;
h=figure; hold on;
set(gcf,'Position',[100 100 1400 400]);
scatter(time(idx10),hot.boxy(idx10),'.')
scatter(time(idx10),hot.gobai_o2(idx10),'.');
scatter(gobai.time,gobai.o2_10,'.');
datetick('x');
legend({'HOT (0-20dbar)' 'HOT_{GOBAI-Alg.}' 'GOBAI-O2 (10dbar)'});
export_fig(h,'O2/Figures/HOTvGOBAI_10dbar.png'); close
% plot HOT vs. GOBAI (100 dbar)
idx100 = hot.press >= 90 & hot.press <= 110;
h=figure; hold on;
set(gcf,'Position',[100 100 1400 400]);
scatter(time(idx100),hot.boxy(idx100),'.')
scatter(time(idx100),hot.gobai_o2(idx100),'.');
scatter(gobai.time,gobai.o2_100,'.');
datetick('x');
legend({'HOT (90-110dbar)' 'HOT_{GOBAI-Alg.}' 'GOBAI-O2 (100dbar)'});
export_fig(h,'O2/Figures/HOTvGOBAI_100dbar.png');
% plot HOT vs. GOBAI (1000 dbar)
idx1000 = hot.press >= 950 & hot.press <= 1150;
h=figure; hold on;
set(gcf,'Position',[100 100 1400 400]);
scatter(time(idx1000),hot.boxy(idx1000),'.')
scatter(time(idx1000),hot.gobai_o2(idx1000),'.');
scatter(gobai.time,gobai.o2_1000,'.');
datetick('x');
legend({'HOT (950-1150dbar)' 'HOT_{GOBAI-Alg.}' 'GOBAI-O2 (1000dbar)'});
export_fig(h,'O2/Figures/HOTvGOBAI_1000dbar.png');

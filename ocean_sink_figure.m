% GCB Sink Figure

%% gcb data
websave('Data/Global_Carbon_Budget_2025_v1.0.xlsx',...
'https://globalcarbonbudget.org/download/2341/?tmstv=1762815901');
gcb_f = 'Data/Global_Carbon_Budget_2025_v1.0.xlsx'; % filename from 2025 budget
os = readtable(gcb_f,'Sheet','Ocean Sink','NumHeaderLines',31); % fluxes in PgC

%% gobai and rfrom data
ver = 'HR-v202608'; % product version
modifier = 'fg';
path = ['/fast6/dic/GOBAI/RFROM/FFNN/c20_Aug-2026_D/train80_val10_test10/' modifier '/']; % file path

if ~isfile(['DIC/Data/dic_and_canth_inventories_' modifier '_' ver '.mat'])

    gobai_f = [path 'GOBAI-Monthly-Mean-1x1-DIC-' ver '.nc']; % file name
    % gobai_f = [path 'gobai-monthly-mean-1x1-dic.nc']; % file name
    rfrom_t_f = '/fast2/temperature/RFROM_TEMP_v2.2_2025/RFROMV22_TEMP_MONTHLY_MEAN_1x1.nc';
    rfrom_s_f = '/fast3/salinity/RFROM_SAL_v2.2_2025/RFROMV22_SAL_MONTHLY_MEAN_1x1.nc';
    % download dimensions
    gobai.lon = ncread(gobai_f,'longitude');
    gobai.lat = ncread(gobai_f,'latitude');
    gobai.pres = ncread(gobai_f,'mean_pressure');
    gobai.time = datenum(1950,1,1) + ncread(gobai_f,'time');
    % calculate mass per grid cell
    [gobai.vol,gobai.area] = weights3d(gobai.lon,gobai.lat,gobai.pres);
    gobai.area = gobai.area(:,:,1);
    % determine mask (time varying coverage)
    gobai.mask = true(size(gobai.vol));
    for t = 1:length(gobai.time)
        gobai.dic = ncread(gobai_f,'dic',[1 1 1 t],[Inf Inf Inf 1]);
        gobai.mask(isnan(gobai.dic)) = false;
    end

    %% determine inventory (PgC)
    [lon3d,lat3d,pres3d] = ndgrid(gobai.lon,gobai.lat,gobai.pres);
    yr = datevec(double(gobai.time)); yr = yr(:,1);
    dic_inv = nan(length(gobai.time),1);
    canth_inv = nan(length(gobai.time),1);
    parfor t = 1:length(gobai.time)
        dic = ncread(gobai_f,'dic',[1 1 1 t],[Inf Inf Inf 1]);
        cns_tmp = ncread(rfrom_t_f,'ocean_temperature',[1 1 1 t],[Inf Inf Inf 1]); % conservative temperature
        abs_sal = ncread(rfrom_s_f,'ocean_salinity',[1 1 1 t],[Inf Inf Inf 1]); % absolute salinity
        tmp = gsw_t_from_CT(abs_sal,cns_tmp,pres3d); % in-situ temperature
        sal = gsw_SP_from_SA(abs_sal,pres3d,lon3d,lat3d); % practical salinity
        dens = gsw_rho(abs_sal,cns_tmp,pres3d); % kg/m3
        % calculate global dic inventory 
        dic_inv(t) = ... % Pg = umol/kg * kg/m3 * m3 * g/mol * mol/umol
            sum(dic(gobai.mask).*dens(gobai.mask).*gobai.vol(gobai.mask).*12.01.*10^-21,'omitnan');
        % apply TRACE to calculate canth
        lon2d = lon3d(:,:,1); lat2d = lat3d(:,:,1);
        canth_3d = nan(size(lon3d));
        for p = 1:length(gobai.pres)
            pres2d = pres3d(:,:,p); tmp2d = tmp(:,:,p); sal2d = sal(:,:,p);
            canth_2d = nan(size(sal2d));
            canth_2d(gobai.mask(:,:,p)) = ...
                TRACEv1([lon2d(gobai.mask(:,:,p)) lat2d(gobai.mask(:,:,p)) ...
                    pres2d(gobai.mask(:,:,p))],yr(t),[tmp2d(gobai.mask(:,:,p)) ...
                    sal2d(gobai.mask(:,:,p))],[2 1],1);
            canth_3d(:,:,p) = canth_2d;
        end
        % calculate global canth inventory 
        canth_inv(t) = ... % Pg = umol/kg * kg/m3 * m3 * g/mol * mol/umol
            sum(canth_3d(gobai.mask).*dens(gobai.mask).*gobai.vol(gobai.mask).*12.01.*10^-21,'omitnan');
        % print status
        fprintf('Processed time step %d of %d\n', t, length(gobai.time));
    end
    % end parallel session
    delete(gcp('nocreate'));
    % save timeseries
    save(['DIC/Data/dic_and_canth_inventories_' modifier '_' ver '.mat'],'yr','dic_inv','canth_inv','gobai','-mat');
else
    % load timeseries
    load(['DIC/Data/dic_and_canth_inventories_' modifier '_' ver '.mat']);
end

%% load canth from Mueller et al. (2023)
m23.fname = 'Data/Muller-2023/dcant_emlr_cstar_mueller_94-04-14.nc';
m23.canth_94_04 = ncread(m23.fname,'dcant',[1 1 1 1 3 5],[Inf Inf Inf 1 1 1]);
m23.canth_04_14 = ncread(m23.fname,'dcant',[1 1 1 2 3 5],[Inf Inf Inf 1 1 1]);
m23.lat = ncread(m23.fname,'lat');
m23.lon = ncread(m23.fname,'lon');
m23.lon = convert_lon(m23.lon,'format','0-360');
m23.depth = ncread(m23.fname,'depth');
[m23.vol,m23.area] = weights3d(m23.lon,m23.lat,m23.depth);
m23.area = m23.area(:,:,1);
% umol/kg/dec * kg/m3 * m3 * g/mol * Pg/g * mol/umol = PgC/dec
m23.mask = ~isnan(m23.canth_94_04);
m23.canth_inv_94_04 = sum(m23.canth_94_04(m23.mask).*1026.*m23.vol(m23.mask).*12.01.*10^-21,'omitnan');
% m23.canth_col_inv = sum(sum(trapz(m23.depth,m23.canth_94_04.*1025,3).*m23.area,'omitnan'),'omitnan').*12.01.*10^-21;
m23.mask = ~isnan(m23.canth_04_14);
m23.canth_inv_04_14 = sum(m23.canth_04_14(m23.mask).*1027.*m23.vol(m23.mask).*12.01.*10^-21,'omitnan');

%% load DIC from Keppler et al. (2023)
if ~isfile('DIC/Data/mobo_inventories.mat')
    mobo_name = 'DIC'; dstn = '/fast6/dic/';
    mobo_local_fname = [dstn 'MPI_MOBO-DIC_2004-2019_v2.nc'];
    mobo.lon = ncread(mobo_local_fname,'lon');
    mobo.lon = convert_lon(mobo.lon,'format','0-360');
    mobo.lat = ncread(mobo_local_fname,'lat');
    mobo.depth = ncread(mobo_local_fname,'depth');
    mobo.time = datenum(2004,ncread(mobo_local_fname,'juld'),1);
    % mobo.var = ncread(mobo_local_fname,mobo_name);
    [lon3d,lat3d,depth3d] = ndgrid(mobo.lon,mobo.lat,mobo.depth);
    rg_pres = double(ncread('Data/RG_CLIM/RG_Climatology_Sal.nc','Pressure'));
    [rg_lon3d,rg_lat3d,rg_pres3d] = ndgrid(mobo.lon,mobo.lat,rg_pres);
    % determine inventory
    [mobo.vol,mobo.area] = weights3d(mobo.lon,mobo.lat,mobo.depth);
    mobo.area = mobo.area(:,:,1);
    mobo_inv = nan(length(mobo.time),1);
    parfor t = 1:length(mobo.time)
        dic = ncread(mobo_local_fname,mobo_name,[1 1 1 t],[Inf Inf Inf 1]);
        mask = dic > 0;
        % tmp = double(ncread('Data/RG_CLIM/RG_Climatology_Temp.nc','Temperature',[1 1 1 t],[Inf Inf Inf 1]));
        % sal = double(ncread('Data/RG_CLIM/RG_Climatology_Sal.nc','Salinity',[1 1 1 t],[Inf Inf Inf 1]));
        % tmp = griddata(rg_lon3d,rg_lat3d,rg_pres3d,tmp,lon3d,lat3d,depth3d);
        % sal = griddata(rg_lon3d,rg_lat3d,rg_pres3d,sal,lon3d,lat3d,depth3d);
        % SA = gsw_SA_from_SP(sal,depth3d,lon3d,lat3d);
        % CT = gsw_CT_from_t(SA,tmp,depth3d);
        % dens = gsw_rho(SA,CT,depth3d); % kg/m3
        % calculate global dic inventory
        mobo_inv(t) = ... % Pg = umol/kg * kg/m3 * m3 * g/mol * mol/umol
            sum(dic(mask).*1027.*mobo.vol(mask).*12.01.*10^-21,'omitnan');
        fprintf('Processed time step %d of %d\n', t, length(mobo.time));
    end
    % end parallel session
    delete(gcp('nocreate'));
    dic = ncread(mobo_local_fname,mobo_name,[1 1 1 1],[Inf Inf Inf 1]);
    mobo.mask = dic > 0;
    % save timeseries
    save('DIC/Data/mobo_inventories.mat','mobo','mobo_inv','-mat');
else
    % load timeseries
    load('DIC/Data/mobo_inventories.mat');
end

%% download GLODAP gridded file if needed
dstn = '/fast6/dic/';
if ~isdir([dstn '/GLODAPv2.2016b_MappedClimatologies'])
    glodap_fname = 'https://www.nodc.noaa.gov/archive/arc0107/0162565/1.1/data/0-data/mapped/GLODAPv2.2016b_MappedClimatologies.tar.gz';
    websave([dstn 'GLODAPv2.2016b_MappedClimatologies.tar.gz'],glodap_fname);
    gunzip([dstn 'GLODAPv2.2016b_MappedClimatologies.tar.gz']);
    untar([dstn 'GLODAPv2.2016b_MappedClimatologies.tar'],dstn);
    delete([dstn 'GLODAPv2.2016b_MappedClimatologies.tar.gz']);
    delete([dstn 'GLODAPv2.2016b_MappedClimatologies.tar']);
end

%% load glodap climatologies
dstn = '/fast6/dic/';
GLODAP.lon = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc'],'lon');
GLODAP.lon = convert_lon(GLODAP.lon,'format','0-360');
GLODAP.lat = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc'],'lat');
GLODAP.depth = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc'],'Depth');
GLODAP.dic = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc'],'TCO2');
GLODAP.cant = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.Cant.nc'],'Cant');
[GLODAP.lon,idx_lon] = sort(GLODAP.lon);
GLODAP.dic = GLODAP.dic(idx_lon,:,:);
GLODAP.cant = GLODAP.cant(idx_lon,:,:);
[GLODAP.vol,GLODAP.area] = weights3d(GLODAP.lon,GLODAP.lat,GLODAP.depth);
GLODAP.area = GLODAP.area(:,:,1);
% DIC inventories
GLODAP.mask = ~isnan(GLODAP.dic);
GLODAP.dic_inv = sum(GLODAP.dic(GLODAP.mask).*1027.*GLODAP.vol(GLODAP.mask).*12.01.*10^-21,'omitnan');
mask_rfrom = GLODAP.mask; mask_rfrom(:,:,27:end) = false;
mask_mobo = GLODAP.mask; mask_mobo(:,:,25:end) = false;
GLODAP.dic_inv_rfrom = sum(GLODAP.dic(GLODAP.mask).*1027.*GLODAP.vol(GLODAP.mask).*12.01.*10^-21,'omitnan')./...
sum(GLODAP.dic(mask_rfrom).*1027.*GLODAP.vol(mask_rfrom).*12.01.*10^-21,'omitnan');
GLODAP.dic_inv_mobo = sum(GLODAP.dic(mask_mobo).*1027.*GLODAP.vol(mask_mobo).*12.01.*10^-21,'omitnan');
% Cant inventories
GLODAP.mask = ~isnan(GLODAP.cant);
GLODAP.cant_inv = sum(GLODAP.cant(GLODAP.mask).*1027.*GLODAP.vol(GLODAP.mask).*12.01.*10^-21,'omitnan');
% mask_rfrom = GLODAP.mask; mask_rfrom(:,:,27:end) = false;
% mask_mobo = GLODAP.mask; mask_mobo(:,:,25:end) = false;
% GLODAP.cant_inv_rfrom = sum(GLODAP.cant(mask_rfrom).*1027.*GLODAP.vol(mask_rfrom).*12.01.*10^-21,'omitnan');
% GLODAP.cant_inv_mobo = sum(GLODAP.cant(mask_mobo).*1027.*GLODAP.vol(mask_mobo).*12.01.*10^-21,'omitnan');

%% scaling factor
% glodap_gobai_dic_ratio = GLODAP.dic_inv / mean(dic_inv); % ratio between GLODAP DIC inventory and mean GOBAI DIC inventory
% glodap_mobo_dic_ratio = GLODAP.dic_inv / mean(mobo_inv); % ratio between GLODAP DIC inventory and MOBO DIC inventory
% glodap_trace_cant_ratio = GLODAP.cant_inv / mean(canth_inv); % ratio between GLODAP Cant inventory and mean TRACE Cant inventory
glodap_gobai_ratio = sum(GLODAP.area(GLODAP.mask(:,:,1)))./...
    sum(gobai.area(gobai.mask(:,:,1)));
glodap_mobo_ratio = sum(GLODAP.area(GLODAP.mask(:,:,1)))./...
    sum(mobo.area(mobo.mask(:,:,1)));
% glodap_m23_ratio = sum(GLODAP.area(GLODAP.mask(:,:,1)))./...
%    sum(m23.area(m23.mask(:,:,1)));
% apply to scale
dic_inv = dic_inv * glodap_gobai_ratio;
mobo_inv = mobo_inv * glodap_mobo_ratio;
canth_inv = canth_inv * glodap_gobai_ratio;
canth_inv_sandborn = [88 117 130 136 149 156 169 186 219]; % PgC, from Sandborn et al. (2023)
canth_time_sandborn = datenum([1980 1994 2000 2002 2007 2010 2014 2020 2030],...
    [1 7 1 7 7 1 7 1 1],1);

%% annual mean of inventories (GOBAI, TRACE, MOBO)
dic_inv_ann = mean(reshape(dic_inv,12,[]),1); % RFROM grid coverage to 2000m
canth_inv_ann = mean(reshape(canth_inv,12,[]),1); % RFROM grid coverage to 2000m
time_ann = double(mean(reshape(gobai.time,12,[]),1));
mobo.time_ann = double(mean(reshape(mobo.time,12,[]),1));
mobo_inv_ann = mean(reshape(mobo_inv,12,[]),1); % mobo grid coverage to 1500m

%% plot annual DIC inv
clrs = orderedcolors('gem');
figure('position',[100 100 1000 600],'Color','w'); hold on;
set(gca,'FontSize',22);
gcb_inv = cumsum(os.GCB);
plot(datenum(os.year,7,0),cumsum(os.GCB),'linewidth',3,'color',clrs(1,:));
for v = [5:14,19:27]
    gcb_inv_temp = cumsum(os.(os.Properties.VariableNames{v}));
    if v < 16
        plot(datenum(os.year,7,0),gcb_inv_temp,'--k');
    else
        plot(datenum(os.year,7,0),gcb_inv_temp,':k');
    end
end
plot(mobo.time_ann+15,mobo_inv_ann-mobo_inv_ann(1)+gcb_inv(46),'linewidth',3,'color',clrs(4,:));
scatter(datenum(2004,7,0),m23.canth_inv_94_04+gcb_inv(36),100,'p','MarkerFaceColor',clrs(5,:),'MarkerEdgeColor','k');
scatter(datenum(2014,7,0),m23.canth_inv_04_14+gcb_inv(46),100,'p','MarkerFaceColor',clrs(5,:),'MarkerEdgeColor','k');
plot(time_ann(1:end-1),dic_inv_ann(1:end-1)-dic_inv_ann(1)+gcb_inv(35),'linewidth',3,'color',clrs(2,:));
% plot(gobai.time,dic_inv-dic_inv_ann(1)+gcb_inv(35),'linewidth',1,'color',clrs(2,:));
plot(time_ann(1:end-1),canth_inv_ann(1:end-1)-canth_inv_ann(1)+gcb_inv(35),'linewidth',3,'color',clrs(3,:));
% plot(canth_time_sandborn,canth_inv_sandborn-canth_inv_sandborn(1)+gcb_inv(22),'linewidth',3,'color',clrs(3,:));
% plot(gobai.time,canth_inv-canth_inv_ann(1)+gcb_inv(35),'linewidth',1,'color',clrs(3,:));
legend({'GCB','GCB (products)','','','','','','','','','','GCB (models)',...
    '','','','','','','','','MOBO-DIC','Mulller23','','GOBAI-DIC','TRACE'},...
    'location','northwest');
    datetick('x','yyyy','keeplimits');
ylabel('\DeltaDIC (PgC)');
export_fig(gcf,['DIC/Figures/GCB_vs_GOBAI_vs_TRACE_inv_' modifier '_' ver '.png']);
close;

%% plot annual ocean uptake
clrs = orderedcolors('gem');
figure('position',[100 100 1000 600],'Color','w'); hold on;
set(gca,'FontSize',22);
plot(datenum(os.year,7,0),os.GCB,'linewidth',3,'color',clrs(1,:));
for v = [5:14,19:27]
    if v < 16
        plot(datenum(os.year,7,0),os.(os.Properties.VariableNames{v}),'--k');
    else
        plot(datenum(os.year,7,0),os.(os.Properties.VariableNames{v}),':k');
    end
end
ylabel('Ocean Carbon Sink (PgC yr^{-1})');
datetick('x','yyyy','keeplimits');
% GCB only
legend({'GCB','GCB (GOBMs)','','','','','','','','','','GCB (fCO2)',...
    '','','','','','','',''},'location','northwest');
export_fig(gcf,['DIC/Figures/GCB_vs_GOBAI_vs_TRACE_flux_' modifier '_' ver '_gcb_only.png']);
% GCB and GOBAI
gb_tmp = plot(time_ann(2:end),diff(dic_inv_ann(1:end)),'linewidth',3,'color',clrs(2,:));
legend({'GCB','GCB (GOBMs)','','','','','','','','','','GCB (fCO2)',...
    '','','','','','','','','GOBAI-DIC'},...
    'location','northwest');
export_fig(gcf,['DIC/Figures/GCB_vs_GOBAI_vs_TRACE_flux_' modifier '_' ver '_gcb_gobai.png']);
delete(gb_tmp);
% GCB and MOBO
plot(mobo.time_ann(2:end),diff(mobo_inv_ann),'linewidth',3,'color',clrs(4,:));
legend({'GCB','GCB (GOBMs)','','','','','','','','','','GCB (fCO2)',...
    '','','','','','','','','MOBO-DIC'},'location','northwest');
export_fig(gcf,['DIC/Figures/GCB_vs_GOBAI_vs_TRACE_flux_' modifier '_' ver '_gcb_mobo.png']);
% GCB, MOBO and Muller
% scatter(datenum(1999,7,0),m23.canth_inv_94_04/10,200,'p','MarkerFaceColor',clrs(5,:),'MarkerEdgeColor','k');
scatter(datenum(2009,7,0),2.7,500,'p','MarkerFaceColor',clrs(5,:),'MarkerEdgeColor','k');
legend({'GCB','GCB (GOBMs)','','','','','','','','','','GCB (fCO2)',...
    '','','','','','','','','MOBO-DIC','Mulller23'},'location','northwest');
export_fig(gcf,['DIC/Figures/GCB_vs_GOBAI_vs_TRACE_flux_' modifier '_' ver '_gcb_mobo_m23.png']);
% GCB, MOBO, Muller and TRACE
plot(time_ann(2:end),diff(canth_inv_ann(1:end)),'linewidth',3,'color',clrs(3,:));
legend({'GCB','GCB (GOBMs)','','','','','','','','','','GCB (fCO2)',...
    '','','','','','','','','MOBO-DIC','Mulller23','TRACE'},...
    'location','northwest');
export_fig(gcf,['DIC/Figures/GCB_vs_GOBAI_vs_TRACE_flux_' modifier '_' ver '_gcb_mobo_m23_trace.png']);
% GCB, MOBO, Muller, TRACE and GOBAI
plot(time_ann(2:end),diff(dic_inv_ann(1:end)),'linewidth',3,'color',clrs(2,:));
legend({'GCB','GCB (GOBMs)','','','','','','','','','','GCB (fCO2)',...
    '','','','','','','','','MOBO-DIC','Mulller23','TRACE','GOBAI-DIC'},...
    'location','northwest');
export_fig(gcf,['DIC/Figures/GCB_vs_GOBAI_vs_TRACE_flux_' modifier '_' ver '_gcb_mobo_m23_trace_gobai.png']);
close;
% inventory trends
idx_gcb = os.year >= 1993 & os.year <= 2024;
gcb_trend = polyfit(datenum(os.year(idx_gcb),7,0),gcb_inv(idx_gcb),1);
gcb_trend = gcb_trend(1)*365.2425;
fprintf('GCB trend (1993-2024): %.3f PgC yr^{-1}\n', gcb_trend);
% find GOBAI-DIC trend
gobai_dic_trend = polyfit(time_ann(1:end-1),dic_inv_ann(1:end-1),1);
gobai_dic_trend = gobai_dic_trend(1)*365.2425;
fprintf('GOBAI-DIC trend (1993-2024): %.3f PgC yr^{-1}\n', gobai_dic_trend);
mobo_dic_trend = polyfit(mobo.time_ann,mobo_inv_ann,1);
mobo_dic_trend = mobo_dic_trend(1)*365.2425;
fprintf('MOBO-DIC trend (2004-2019): %.3f PgC yr^{-1}\n', mobo_dic_trend);
canth_inv_trend = polyfit(time_ann(1:end-1),canth_inv_ann(1:end-1),1);
canth_inv_trend = canth_inv_trend(1)*365.2425;
fprintf('TRACE trend (1993-2024): %.3f PgC yr^{-1}\n', canth_inv_trend);

% % pentadal means
% year_ann = datevec(time_ann); year_ann = year_ann(:,1);
% period_start = [1993, 1996:5:(max(year_ann)-4)];
% period_end = [1995, period_start(2:end)+4];
% dic_inv_pent = nan(size(period_start));
% for p = 1:length(period_start)
%     idx = year_ann >= period_start(p) & year_ann <= period_end(p);
%     dic_inv_pent(p) = mean(dic_inv_ann(idx),'omitnan');
% end
% dic_inv_pent_year = datenum((period_start+period_end)./2,7,0);
% plot(dic_inv_pent_year(2:end),diff(dic_inv_pent),'linewidth',3,'color',clrs(2,:));
% annual change
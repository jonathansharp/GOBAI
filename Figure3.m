% Pressures to plot
pressures = [10 100 500 1000];

%% Hawaii Ocean Timeseries
% HOT Dataset info
inf = ncinfo('Data/HOT_Data.nc');
num_clusters = 15;
hot.lon = -158;
hot.lat = 22.75;
% read data
for v = 1:length(inf.Variables)
    hot.(inf.Variables(v).Name) = ...
        ncread('Data/HOT_Data.nc',inf.Variables(v).Name);
    % replace -9 with NaN
    hot.(inf.Variables(v).Name)(hot.(inf.Variables(v).Name) == -9) = NaN;
end
% index
idx = hot.mdate < 0 | hot.boxy < 0;
for v = 1:length(inf.Variables)
    if isnumeric(hot.(inf.Variables(v).Name))
        hot.(inf.Variables(v).Name)(idx) = []; 
    end
end
% process time
hot.month = floor(hot.mdate/10000);
hot.day = floor(mod(hot.mdate,10000)/100);
hot.year = mod(hot.mdate,100);
hot.year(hot.year>50) = hot.year(hot.year>50)+1900;
hot.year(hot.year<50) = hot.year(hot.year<50)+2000;
hot.time = datenum(hot.year,hot.month,hot.day); % for each day
hot.time0 = datenum(hot.year,1,1);
hot.doy = hot.time-hot.time0;
% calculate absolute salinity and conservative temperature
hot.sal_abs = gsw_SA_from_SP(hot.csal,hot.press,hot.lon,hot.lat);
hot.temp_cns = gsw_CT_from_t(hot.sal_abs,hot.temp,hot.press);
% plot oxygen
plot_timeseries_comparison('v1.1-HR','','O2','o2',hot,'HOT','boxy',pressures,'[O_{2}]');
% plot nitrate
% plot_timeseries_comparison('v2.3','NO3','no3',hot,'HOT','nit',pressures,'[NO_{3}]');
% plot dic
% plot_timeseries_comparison('v2.3','DIC','dic',hot,'HOT','dic',pressures,'DIC (\mumol kg^{-1})');

clearvars -except pressures

%% Bermuda Atlantic Timeseries
web_path = 'https://datadocs.bco-dmo.org/dataset/917255/';
web_file = 'file/M7rPoqyt6vYv4V/917255_v8_bats_bval_bottle.csv';
websave('Data/BATS_data',[web_path web_file]);
data_table = readtable('BATS_data.csv');
% process variables
bats.lon = data_table.Longitude;
bats.lat = data_table.Latitude;
bats.press = gsw_p_from_z(-data_table.Depth,data_table.Latitude);
bats.o2 = data_table.O2;
% process time
bats.time = datenum(num2str(data_table.yyyymmdd),'yyyymmdd');
% plot oxygen
plot_timeseries_comparison('v2.4','-depth-dependent-adjustment','O2','o2',bats,'BATS','o2',pressures,'[O_{2}]');
% plot nitrate
%plot_timeseries_comparison('v2.3','','NO3','no3',bats,'BATS','nit',pressures,'[NO_{3}]');
% plot dic
%plot_timeseries_comparison('v2.3','','DIC','dic',bats,'BATS','dic',pressures,'DIC (\mumol kg^{-1})');

clearvars -except pressures

%% Hawaii Ocean Timeseries
% HOT Dataset info
inf = ncinfo('Data/HOT_Data.nc');
num_clusters = 15;
hot.lon = -158;
hot.lat = 22.75;
% read data
for v = 1:length(inf.Variables)
    hot.(inf.Variables(v).Name) = ...
        ncread('Data/HOT_Data.nc',inf.Variables(v).Name);
    % replace -9 with NaN
    hot.(inf.Variables(v).Name)(hot.(inf.Variables(v).Name) == -9) = NaN;
end
% index
idx = hot.mdate < 0 | hot.boxy < 0;
for v = 1:length(inf.Variables)
    if isnumeric(hot.(inf.Variables(v).Name))
        hot.(inf.Variables(v).Name)(idx) = []; 
    end
end
% process time
hot.month = floor(hot.mdate/10000);
hot.day = floor(mod(hot.mdate,10000)/100);
hot.year = mod(hot.mdate,100);
hot.year(hot.year>50) = hot.year(hot.year>50)+1900;
hot.year(hot.year<50) = hot.year(hot.year<50)+2000;
hot.time = datenum(hot.year,hot.month,hot.day); % for each day
hot.time0 = datenum(hot.year,1,1);
hot.doy = hot.time-hot.time0;
% calculate absolute salinity and conservative temperature
hot.sal_abs = gsw_SA_from_SP(hot.csal,hot.press,hot.lon,hot.lat);
hot.temp_cns = gsw_CT_from_t(hot.sal_abs,hot.temp,hot.press);
% plot oxygen
plot_timeseries_comparison('v1.1-HR','O2','o2',hot,'HOT','boxy',pressures,'[O_{2}]');
% plot nitrate
plot_timeseries_comparison('v1.1-HR','NO3','no3',hot,'HOT','nit',pressures,'[NO_{3}]');
% plot dic
plot_timeseries_comparison('v1.1-HR','DIC','dic',hot,'HOT','dic',pressures,'DIC (\mumol kg^{-1})');

clearvars -except pressures

%% plotting function
function plot_timeseries_comparison(ver,ext,var,gobai_var_name,ts_data,...
    ts_name,ts_var_name,pres,var_label)

        % calculate o2 with gobai algorithm
        % hot.(['gobai_' gobai_var_name]) = ...
        %     gobai_nn(num_clusters,'temp_cns',ts_data.temp_cns,'sal_abs',ts_data.sal_abs,...
        %         'sigma',ts_data.sigma,'lon',-157.87,'lat',21.31,'pres',ts_data.press,...
        %         'doy',doy,'year',year);
        % define gobai file name
        path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
        % download dimensions
        GOBAI.lon = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'lon');
        GOBAI.lat = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'lat');
        GOBAI.pres = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'pres');
        GOBAI.time = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'time');
        if strcmp(ver,'v2.1')
            GOBAI.time = datenum(2004,1,1)+GOBAI.time;
        else
            GOBAI.time = datenum(1950,1,1)+GOBAI.time;
        end
        % download oxygen
        [~,idx_lon] = min(abs(GOBAI.lon-mean(ts_data.lon+360)));
        [~,idx_lat] = min(abs(GOBAI.lat-mean(ts_data.lat)));
        
        GOBAI.var = squeeze(ncread([path 'GOBAI-' var '-' ver ext '.nc'],...
            gobai_var_name,[idx_lon idx_lat 1 1],[1 1 Inf Inf]));

    for p = 1:length(pres)
        [~,idx_pres] = min(abs(GOBAI.pres-pres(p)));
        % calculate o2 with gobai algorithm
        % hot.(['gobai_' gobai_var_name]) = ...
        %     gobai_nn(num_clusters,'temp_cns',ts_data.temp_cns,'sal_abs',ts_data.sal_abs,...
        %         'sigma',ts_data.sigma,'lon',-157.87,'lat',21.31,'pres',ts_data.press,...
        %         'doy',doy,'year',year);
        % define gobai file name
        % plot timeseries vs. gobai
        if pres(p) <= 100; presmin = pres(p)-10; presmax = pres(p)+10;
        elseif pres(p) <= 500; presmin = pres(p)-20; presmax = pres(p)+20;
        elseif pres(p) <= 1000; presmin = pres(p)-50; presmax = pres(p)+50;
        else; presmin = pres(p)-1000; presmax = pres(p)+1000; end
        idx_pres_ts = ts_data.press >= presmin & ts_data.press <= presmax;
        % establish figure
        h=figure; hold on;
        set(h,'Position',[100 100 1400 400]);
        set(gca,'FontSize',14);
        clrs = get(gca,'ColorOrder');
        % % average each day
        % days = unique(ts_data.time(idx_pres_ts));
        % ts_data.([var '_avg']) = nan(size(days));
        % ts_data.([var '_std']) = nan(size(days));
        % for t = 1:length(GOBAI.time)
        %     idx_t = ts_data.time == days(t);
        %     ts_data.([var '_avg'])(t) = ...
        %         mean(ts_data.(ts_var_name)(idx_pres_ts & idx_t),'omitnan');
        %     ts_data.([var '_std'])(t) = ...
        %         std(ts_data.(ts_var_name)(idx_pres_ts & idx_t),[],'omitnan');
        % end
        % average each gobai timestep
        ts_data.([var '_avg']) = nan(size(GOBAI.time));
        ts_data.([var '_std']) = nan(size(GOBAI.time));
        time_mids = (GOBAI.time(1:end-1)+GOBAI.time(2:end))./2;
        start_bound = GOBAI.time(1) - (time_mids(1) - GOBAI.time(1));
        end_bound = GOBAI.time(end) + (GOBAI.time(end) - time_mids(end));
        time_mids = [start_bound;time_mids;end_bound];
        for t = 1:length(GOBAI.time)
            idx_t = ts_data.time >= time_mids(t) & ts_data.time < time_mids(t+1);
            ts_data.([var '_avg'])(t) = ...
                mean(ts_data.(ts_var_name)(idx_pres_ts & idx_t),'omitnan');
            ts_data.([var '_std'])(t) = ...
                std(ts_data.(ts_var_name)(idx_pres_ts & idx_t),[],'omitnan');
        end
        % plot
        idx_nan = ~isnan(ts_data.([var '_avg']));
        % fill([GOBAI.time;flipud(days(idx_nan))],...
        %     [ts_data.([var '_avg'])(idx_nan)+ts_data.([var '_std'])(idx_nan);...
        %     flipud(ts_data.([var '_avg'])(idx_nan)-ts_data.([var '_std'])(idx_nan))],...
        %     clrs(1,:),'LineStyle','none','FaceAlpha',0.25);
        p1=errorbar(double(GOBAI.time(idx_nan)),ts_data.([var '_avg'])(idx_nan),...
            ts_data.([var '_std'])(idx_nan),'Color',clrs(1,:),...
            'Marker','.','MarkerSize',10,'LineStyle','none'); % timeseries data
        % scatter(time(idx),ts_data.gobai_o2(idx),'.'); % algorithm prediction
        p2=plot(double(GOBAI.time),GOBAI.var(idx_pres,:),'Color',clrs(2,:),...
            'Marker','.','MarkerSize',10); % gobai data
        %fill(GOBAI.time,GOBAI_uncer)
        datetick('x'); ylabel(var_label);
        % legend({['ts_name ' (0-20dbar)'],['HOT_{GOBAI-Alg.}'],['GOBAI-' var ' (10dbar)']});
        legend([p1 p2],{[ts_name ' (' num2str(presmin) '-' num2str(presmax) 'dbar)'],...
            ['GOBAI-' var '-' ver ' (' num2str(pres(p)) 'dbar)']});
        export_fig(h,[var '/Figures/' ts_name 'vGOBAI-' ver ext '_' ...
            num2str(pres(p)) 'dbar.png']);
        close;
    end

end
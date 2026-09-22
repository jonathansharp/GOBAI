% Pressures to plot
pressures = [10 200 400 1250];

%% Hawaii Ocean Timeseries
% % HOT Dataset info
% inf = ncinfo('Data/HOT_Data.nc');
% num_clusters = 15;
% hot.lon = -158;
% hot.lat = 22.75;
% % read data
% for v = 1:length(inf.Variables)
%     hot.(inf.Variables(v).Name) = ...
%         ncread('Data/HOT_Data.nc',inf.Variables(v).Name);
%     % replace -9 with NaN
%     hot.(inf.Variables(v).Name)(hot.(inf.Variables(v).Name) == -9) = NaN;
% end
% % index
% idx = hot.mdate < 0 | hot.boxy < 0;
% for v = 1:length(inf.Variables)
%     if isnumeric(hot.(inf.Variables(v).Name))
%         hot.(inf.Variables(v).Name)(idx) = []; 
%     end
% end
% % process time
% hot.month = floor(hot.mdate/10000);
% hot.day = floor(mod(hot.mdate,10000)/100);
% hot.year = mod(hot.mdate,100);
% hot.year(hot.year>50) = hot.year(hot.year>50)+1900;
% hot.year(hot.year<50) = hot.year(hot.year<50)+2000;
% hot.time = datenum(hot.year,hot.month,hot.day); % for each day
% hot.time0 = datenum(hot.year,1,1);
% hot.doy = hot.time-hot.time0;
% % calculate absolute salinity and conservative temperature
% hot.sal_abs = gsw_SA_from_SP(hot.csal,hot.press,hot.lon,hot.lat);
% hot.temp_cns = gsw_CT_from_t(hot.sal_abs,hot.temp,hot.press);
% % pre-allocate statistics
% hot_stats = nan(3,4);
% % plot oxygen
% [hot_stats(1,1),hot_stats(1,2),hot_stats(1,3),hot_stats(1,4)] = ...
%     plot_timeseries_comparison('HR-v1.3','','O2','o2',hot,'HOT','boxy',pressures,'[O_{2}]');
% % plot nitrate
% [hot_stats(2,1),hot_stats(2,2),hot_stats(2,3),hot_stats(2,4)] = ...
%     plot_timeseries_comparison('HR-v1.3','','NO3','no3',hot,'HOT','nit',pressures,'[NO_{3}]');
% % plot dic
% [hot_stats(3,1),hot_stats(3,2),hot_stats(3,3),hot_stats(3,4)] = ...
%     plot_timeseries_comparison('HR-v1.3','','DIC','dic',hot,'HOT','dic',pressures,'DIC (\mumol kg^{-1})');
% % clean up
% clearvars -except pressures

%% Bermuda Atlantic Timeseries
% web_path = 'https://datadocs.bco-dmo.org/dataset/917255/';
% web_file = 'file/M7rPoqyt6vYv4V/917255_v8_bats_bval_bottle.csv';
% websave('Data/BATS_data',[web_path web_file]);
% data_table = readtable('BATS_data.csv');
% % process variables
% bats.lon = data_table.Longitude;
% bats.lat = data_table.Latitude;
% bats.press = gsw_p_from_z(-data_table.Depth,data_table.Latitude);
% bats.o2 = data_table.O2;
% bats.nit = data_table.NO3_NO2;
% bats.dic = data_table.DIC;
% % process time
% bats.time = datenum(num2str(data_table.yyyymmdd),'yyyymmdd');
% % pre-allocate statistics
% bat_stats = nan(3,4);
% % plot oxygen
% [bat_stats(1,1),bat_stats(1,2),bat_stats(1,3),bat_stats(1,4)] = ...
%     plot_timeseries_comparison('HR-v1.3','','O2','o2',bats,'BATS','o2',pressures,'[O_{2}]');
% % plot nitrate
% [bat_stats(2,1),bat_stats(2,2),bat_stats(2,3),bat_stats(2,4)] = ...
%     plot_timeseries_comparison('HR-v1.3','','NO3','no3',bats,'BATS','nit',pressures,'[NO_{3}]');
% % plot dic
% [bat_stats(3,1),bat_stats(3,2),bat_stats(3,3),bat_stats(3,4)] = ...
%     plot_timeseries_comparison('HR-v1.3','','DIC','dic',bats,'BATS','dic',pressures,'DIC (\mumol kg^{-1})');
% % clean up
% clearvars -except pressures

%% Ocean Station Papa
% download OSP dataset
moorA = 'https://erddap.dataexplorer.oceanobservatories.org/erddap/tabledap/ooi-gp03flma-ris01-03-dostad000.nc?time%2Cmoles_of_oxygen_per_unit_mass_in_sea_water%2Cmoles_of_oxygen_per_unit_mass_in_sea_water_qc_agg%2Cz&time%3E%3D2013-07-21T22%3A45%3A00Z&time%3C%3D2025-10-13T16%3A00%3A00Z';
moorB = 'https://erddap.dataexplorer.oceanobservatories.org/erddap/tabledap/ooi-gp03flmb-ris01-03-dostad000.nc?time%2Cmoles_of_oxygen_per_unit_mass_in_sea_water%2Cmoles_of_oxygen_per_unit_mass_in_sea_water_qc_agg%2Cz&time%3E%3D2013-07-24T06%3A45%3A00Z&time%3C%3D2025-10-12T00%3A00%3A00Z';
websave('Data/ospA.nc',moorA);
websave('Data/ospB.nc',moorB);
% OSP Dataset info
infA = ncinfo('Data/ospA.nc');
infB = ncinfo('Data/ospB.nc');
num_clusters = 15;
ospA.lon = -144.36; ospA.lat = 50.02;
ospB.lon = -144.51; ospB.lat = 50.38;
% read data
for v = 1:length(infA.Variables)
    ospA.(infA.Variables(v).Name) = ...
        ncread('Data/ospA.nc',infA.Variables(v).Name);
end
for v = 1:length(infB.Variables)
    ospB.(infB.Variables(v).Name) = ...
        ncread('Data/ospB.nc',infB.Variables(v).Name);
end
% process time
ospA.time = ospA.time + datenum(1970,1,1);
ospB.time = ospB.time + datenum(1970,1,1);
% calculate absolute salinity and conservative temperature
% hot.sal_abs = gsw_SA_from_SP(hot.csal,hot.press,hot.lon,hot.lat);
% hot.temp_cns = gsw_CT_from_t(hot.sal_abs,hot.temp,hot.press);
% plot oxygen
plot_timeseries_comparison('HR-v202606','','O2','o2',ospA,'OSP',...
    'moles_of_oxygen_per_unit_mass_in_sea_water',pressures,'[O_{2}]');
% plot nitrate
plot_timeseries_comparison('v1.1-HR','NO3','no3',hot,'HOT','nit',pressures,'[NO_{3}]');
% plot dic
plot_timeseries_comparison('v1.1-HR','DIC','dic',hot,'HOT','dic',pressures,'DIC (\mumol kg^{-1})');

clearvars -except pressures

%% plotting function
function [bias,rms_error,r,gobai_minus_ts_ann_amp] = ...
    plot_timeseries_comparison(ver,ext,var,gobai_var_name,ts_data,...
    ts_name,ts_var_name,pres,var_label)

        % calculate o2 with gobai algorithm
        % hot.(['gobai_' gobai_var_name]) = ...
        %     gobai_nn(num_clusters,'temp_cns',ts_data.temp_cns,'sal_abs',ts_data.sal_abs,...
        %         'sigma',ts_data.sigma,'lon',-157.87,'lat',21.31,'pres',ts_data.press,...
        %         'doy',doy,'year',year);
        % define gobai file name
        path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
        % download dimensions
        if strcmp(ver,'HR-v1.4')
            files = dir([path 'monthly/']);
            GOBAI.time = [];
            for n=3:length(files)
                if n == 3
                    GOBAI.lon = ncread([path 'monthly/' files(n).name],'longitude');
                    GOBAI.lat = ncread([path 'monthly/' files(n).name],'latitude');
                    GOBAI.pres = ncread([path 'monthly/' files(n).name],'mean_pressure');
                else
                    GOBAI.time = [GOBAI.time;ncread([path 'monthly/' files(n).name],'time')];
                end
            end
        else
            GOBAI.lon = ncread([path 'GOBAI-' var '-HR-v202606' ext '.nc'],'lon');
            GOBAI.lat = ncread([path 'GOBAI-' var '-HR-v202606' ext '.nc'],'lat');
            GOBAI.pres = ncread([path 'GOBAI-' var '-HR-v202606' ext '.nc'],'pres');
            GOBAI.time = ncread([path 'GOBAI-' var '-HR-v202606' ext '.nc'],'time');
%             GOBAI.lon = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'lon');
%             GOBAI.lat = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'lat');
%             GOBAI.pres = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'pres');
%             GOBAI.time = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'time');
        end
        if strcmp(ver,'v2.1')
            GOBAI.time = datenum(2004,1,1)+GOBAI.time;
        else
            GOBAI.time = datenum(1950,1,1)+GOBAI.time;
        end
        % download oxygen
        [~,idx_lon] = min(abs(GOBAI.lon-mean(ts_data.lon+360)));
        [~,idx_lat] = min(abs(GOBAI.lat-mean(ts_data.lat)));
        if strcmp(ver,'HR-v1.4')
            files = dir([path 'monthly/']);
            GOBAI.var = [];
            for n=3:length(files)
                GOBAI.var = [GOBAI.var,squeeze(ncread([path 'monthly/' files(n).name],...
                    gobai_var_name,[idx_lon idx_lat 1 1],[1 1 Inf Inf]))];
            end
        else
            GOBAI.var = squeeze(ncread([path 'GOBAI-' var '-HR-v202606' ext '.nc'],...
                gobai_var_name,[idx_lon idx_lat 1 1],[1 1 Inf Inf]));
%             GOBAI.var = squeeze(ncread([path 'GOBAI-' var '-' ver ext '.nc'],...
%                 gobai_var_name,[idx_lon idx_lat 1 1],[1 1 Inf Inf]));
        end
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
        else; presmin = pres(p)-100; presmax = pres(p)+100; end
        if strcmp(ts_name,'OSP')
            idx_pres_ts = ts_data.z >= presmin & ts_data.z <= presmax;
        else
            idx_pres_ts = ts_data.press >= presmin & ts_data.press <= presmax;
        end
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
        % calculate statistics for paper
        gobai_var = GOBAI.var(idx_pres,:)';
        gobai_var = gobai_var(idx_nan);
        ts_var = ts_data.([var '_avg'])(idx_nan);
        time = double(GOBAI.time(idx_nan));
        bias = mean(gobai_var-ts_var);
        rms_error = rmse(gobai_var,ts_var);
        r = corrcoef(gobai_var,ts_var); r = r(1,end);
        [gobai_fit,gobai_resid,gobai_x] = leastsq2(time,gobai_var,...
            time(1),2,[365.2425 365.2425/2 ]);
        gobai_ann_amp = sqrt(gobai_x(3)^2+gobai_x(4)^2);
        gobai_semiann_amp = sqrt(gobai_x(5)^2+gobai_x(6)^2);
        [ts_fit,ts_resid,ts_x] = leastsq2(time,ts_var,...
            time(1),2,[365.2425/2 365.2425]);
        ts_ann_amp = sqrt(ts_x(3)^2+ts_x(4)^2);
        ts_semiann_amp = sqrt(ts_x(5)^2+ts_x(6)^2);
        gobai_minus_ts_ann_amp = gobai_ann_amp - ts_ann_amp;
    end

end
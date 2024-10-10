% Glob_climatology_sla.m
%
% DESCRIPTION:
% This script 
% 
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 6/6/2024

%% establish paths, folders, and options
% file paths and names
fpath_sla = '/raid/Data/CMEMS/SEALEVEL_GLO_PHY_L4_MY_008_047/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D_202112';
fpath_rfrom = [pwd '/Data/RFROM/'];
folder_rfrom = 'RFROM_TEMP_v0.1';
fname_rfrom = '/RFROM_TEMP_STABLE_';
folder_sla_rfrom = 'RFROM_SLA_CMEMS';
fname_sla_rfrom = 'RFROM_SLA_';
% make new directory if it does not exist
if ~exist([fpath_rfrom folder_sla_rfrom],'dir')
    mkdir([fpath_rfrom folder_sla_rfrom])
end
% define timespan
y1 = 2004; y2 = 2021;
% pre-allocate variables
filenames = [];
for y = y1:y2
    for m = 1:12
        dir_temp = dir([fpath_sla '/' num2str(y) '/' sprintf('%02d',m)]);
        for n = 3:length(dir_temp)
            filenames = [filenames;dir_temp(n).name];
        end
        clear dir_temp
    end
end
sla_year = nan(length(filenames),1);
sla_month = nan(length(filenames),1);
sla_day = nan(length(filenames),1);

%% get and extract dates from sla file list
for l = 1:length(filenames)
    sla_date = cell2mat(extractBetween(filenames(l,:),'l4_','_'));
    sla_year(l) = str2double(sla_date(1:4));
    sla_month(l) = str2double(sla_date(5:6));
    sla_day(l) = str2double(sla_date(7:8));
end
file_datenum = datenum(sla_year,sla_month,sla_day);
clear filenames sla_date
% get schema from sla file
f_name_temp = '/2004/01/dt_global_allsat_phy_l4_20040101_20210726.nc';
sla_schema = ncinfo([fpath_sla f_name_temp]);
% get dimensions from sla file
sla_lon = ncread([fpath_sla f_name_temp],'longitude');
sla_lat = ncread([fpath_sla f_name_temp],'latitude');
% get dimensions from RFROM file
rfrom_lon = ncread([fpath_rfrom folder_rfrom fname_rfrom '2004_01.nc'],'longitude');
rfrom_lat = ncread([fpath_rfrom folder_rfrom fname_rfrom '2004_01.nc'],'latitude');

%% loop through each year and month and save RFROM equivalent of sea level anomaly
for y = y2
    for m = 1:12
        % determine weekly times from RFROM file
        rfrom_datenum = datenum(1950,1,1)+...
            ncread([fpath_rfrom folder_rfrom fname_rfrom num2str(y) '_' ...
                sprintf('%02s',num2str(m)) '.nc'],'time');
        % get schema from RFROM file
        rfrom_schema = ncinfo([fpath_rfrom folder_rfrom fname_rfrom num2str(y) '_' ...
                sprintf('%02s',num2str(m)) '.nc']);
        % get time from RFROM file
        rfrom_time = ncread([fpath_rfrom folder_rfrom fname_rfrom num2str(y) '_' ...
                sprintf('%02s',num2str(m)) '.nc'],'time');
        % alter schema for RFROM sla file
        rfrom_sla_schema = rfrom_schema;
        idx_rem = [];
        for v = 1:length(rfrom_sla_schema.Variables)
            if contains(rfrom_sla_schema.Variables(v).Name,'pressure')
                idx_rem = [idx_rem v];
            elseif contains(rfrom_sla_schema.Variables(v).Name,'temperature')
                idx_sla = v; 
            end
        end
        rfrom_sla_schema.Variables(idx_sla).Name = 'sla';
        rfrom_sla_schema.Variables(idx_sla).Attributes(1).Value = 'm';
        rfrom_sla_schema.Variables(idx_sla).Attributes(2).Value = ...
            'The sea level anomaly is the sea surface height above mean sea surface; it is referenced to the [1993, 2012] period; see the product user manual for details';
        rfrom_sla_schema.Variables(idx_sla).Dimensions(3) = [];
        rfrom_sla_schema.Variables(idx_rem) = [];
        rfrom_sla_schema.Attributes = [];
        rfrom_sla_schema.Dimensions(4:5) = [];
        % delete NetCDF if it exists
        fname = [fpath_rfrom folder_sla_rfrom '/' fname_sla_rfrom ...
            num2str(y) '_' sprintf('%02s',num2str(m)) '.nc'];
        if exist(fname,'file'); delete(fname); end
        % write schema to NetCDF
        ncwriteschema(fname,rfrom_sla_schema);
        % save dimensions to NetCDF
        ncwrite(fname,'longitude',rfrom_lon);
        ncwrite(fname,'latitude',rfrom_lat);
        ncwrite(fname,'time',rfrom_time);
        % log weekly average sea level anomaly values in the new file
        for w = 1:length(rfrom_datenum)
            % find seven days of the week
            index = find(file_datenum >= rfrom_datenum(w)-3 & file_datenum <= rfrom_datenum(w)+3);
            % pre-allocate daily sea level anomaly
            sla = nan(length(sla_lon),length(sla_lat),length(index));
            % log daily sea level anomaly
            for d = 1:length(index)
                if datenum(sla_year(index(d)),sla_month(index(d)),sla_day(index(d))) < 738157
                    sla(:,:,d) = ncread([fpath_sla '/' ...
                        num2str(sla_year(index(d))) '/' ...
                        sprintf('%02d',sla_month(index(d))) '/' ...
                        'dt_global_allsat_phy_l4_' ... 
                        num2str(sla_year(index(d))) ...
                        sprintf('%02d',sla_month(index(d))) ...
                        sprintf('%02d',sla_day(index(d))) ...
                        '_20210726.nc'],...
                        'sla');
                elseif datenum(sla_year(index(d)),sla_month(index(d)),sla_day(index(d))) < 738371
                    sla(:,:,d) = ncread([fpath_sla '/' ...
                        num2str(sla_year(index(d))) '/' ...
                        sprintf('%02d',sla_month(index(d))) '/' ...
                        'dt_global_allsat_phy_l4_' ... 
                        num2str(sla_year(index(d))) ...
                        sprintf('%02d',sla_month(index(d))) ...
                        sprintf('%02d',sla_day(index(d))) ...
                        '_20220120.nc'],...
                        'sla');
                else
                    sla(:,:,d) = ncread([fpath_sla '/' ...
                        num2str(sla_year(index(d))) '/' ...
                        sprintf('%02d',sla_month(index(d))) '/' ...
                        'dt_global_allsat_phy_l4_' ... 
                        num2str(sla_year(index(d))) ...
                        sprintf('%02d',sla_month(index(d))) ...
                        sprintf('%02d',sla_day(index(d))) ...
                        '_20220422.nc'],...
                        'sla');
                end

            end
            % average into weekly sea level anomaly
            sla_w = mean(sla,3,'omitnan');
            % convert RFROM longitude to -180 to 180 for interpolation
            rfrom_lon_new = convert_lon(rfrom_lon);
            % interpolate 4km weekly sea level anomaly to RFROM grid
            rfrom_sla = single(interp2(sla_lon',sla_lat,sla_w',...
                rfrom_lon_new',rfrom_lat))'; clear rfrom_lon_new
            % test plots
            % figure; pcolor(sla_lon,sla_lat,sla_w'); shading flat; colorbar; caxis([0 2]);
            % figure; pcolor(rfrom_lon,rfrom_lat,rfrom_sla'); shading flat; colorbar; caxis([0 2]);
            % interpolate over gaps
            % tbd
            % test plot
            % figure; pcolor(rfrom_lon,rfrom_lat,rfrom_sla_i'); shading flat; colorbar; caxis([0 2]);
            % save weekly sea level anomaly to NetCDF
            ncwrite(fname,'sla',rfrom_sla,[1 1 w]);
        end
    end
end

%% create monthly annual files
% make new directory if it does not exist
if ~exist([fpath_rfrom folder_sla_rfrom '_annual'],'dir')
    mkdir([fpath_rfrom folder_sla_rfrom '_annual'])
end
% loop through each year to create annual files
for y = y1:y2
    % pre-allocate annual monthly RFROMs
    RFROM_sla.(['m' num2str(y)]) = [];
    RFROM_sla.(['m' num2str(y) '_time']) = [];
    % extract weekly sea level anomaly values and average to monthly
    for m = 1:12
        fname = [fpath_rfrom folder_sla_rfrom '/' fname_sla_rfrom ...
            num2str(y) '_' sprintf('%02s',num2str(m)) '.nc'];
        RFROM_sla_temp = ncread(fname,'sla');
        RFROM_sla.(['m' num2str(y)]) = ...
            cat(3,RFROM_sla.(['m' num2str(y)]),mean(RFROM_sla_temp,3,'omitnan'));
        RFROM_sla_temp = ncread([fpath_rfrom folder_sla_rfrom '/' fname_sla_rfrom num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'time');
        RFROM_sla.(['m' num2str(y) '_time']) = ...
            [RFROM_sla.(['m' num2str(y) '_time']);mean(RFROM_sla_temp,'omitnan')];
        clear RFROM_temp
    end
    % copy schema to new annual file
    schema = ncinfo([fpath_rfrom folder_sla_rfrom '/' fname_sla_rfrom num2str(y) '_01.nc']);
    schema.Dimensions(3).Length = 12;
    schema.Variables(3).Size = 12;
    schema.Variables(3).Dimensions.Length = 12;
    schema.Variables(3).Attributes(1).Value = 'month of year';
    schema.Variables(3).Attributes(2).Value = 'Number of Month, 1 (Jan.) to 12 (Dec.)';
    schema.Variables(4).Size = [1440,720,12];
    schema.Variables(4).Dimensions(3).Length = 12;
    if isfile([fpath_rfrom folder_sla_rfrom '_annual' '/' fname_sla_rfrom num2str(y) '.nc'])
        delete([fpath_rfrom folder_sla_rfrom '_annual' '/' fname_sla_rfrom num2str(y) '.nc']);
    end
    ncwriteschema([fpath_rfrom folder_sla_rfrom '_annual' '/' fname_sla_rfrom num2str(y) '.nc'],schema);
    clear schema
    % add dimensional variables to monthly file
    vars = {'longitude' 'latitude'};
    for v = 1:length(vars)
        a=ncread([fpath_rfrom folder_sla_rfrom '/' fname_sla_rfrom num2str(y) '_01.nc'],vars{v});
        ncwrite([fpath_rfrom folder_sla_rfrom '_annual' '/' fname_sla_rfrom num2str(y) '.nc'],vars{v},a);
    end
    % add averages to monthly file
    ncwrite([fpath_rfrom folder_sla_rfrom '_annual' '/' fname_sla_rfrom num2str(y) '.nc'],...
        'time',RFROM_sla.(['m' num2str(y) '_time']));
    ncwrite([fpath_rfrom folder_sla_rfrom '_annual' '/' fname_sla_rfrom num2str(y) '.nc'],...
        'sla',RFROM_sla.(['m' num2str(y)]));
    % clean up
    clear RFROM_sla a m v y vars
end

%% create climatological file
% copy schema to new climatological file
schema = ncinfo([fpath_rfrom folder_sla_rfrom '_annual/' fname_sla_rfrom '2004.nc']);
if isfile([fpath_rfrom fname_sla_rfrom 'CLIM.nc'])
    delete([fpath_rfrom fname_sla_rfrom 'CLIM.nc']);
end
ncwriteschema([fpath_rfrom fname_sla_rfrom 'CLIM.nc'],schema);
clear schema
% add dimensional variables to monthly file
vars = {'longitude' 'latitude'};
for v = 1:length(vars)
    a=ncread([fpath_rfrom folder_sla_rfrom '_annual/' fname_sla_rfrom '2004.nc'],vars{v});
    ncwrite([fpath_rfrom fname_sla_rfrom 'CLIM.nc'],vars{v},a);
end
% clean up
clear a m v y vars

% loop through each month to add to climatological file
for m = 1:12
    RFROM_sla_clim = [];
    % extract monthly sea level anomaly values
    for y = y1:y2
        RFROM_sla_temp = ncread([fpath_rfrom folder_sla_rfrom '_annual/' fname_sla_rfrom num2str(y) '.nc'],...
            'sla',[1 1 m],[Inf Inf 1]);
        RFROM_sla_clim = ...
            cat(3,RFROM_sla_clim,RFROM_sla_temp);
        clear RFROM_sla_temp
    end
    % average monthly values to climatological means
    RFROM_sla_clim = mean(RFROM_sla_clim,3,'omitnan');
    % add to NetCDF
    ncwrite([fpath_rfrom fname_sla_rfrom 'CLIM.nc'],'time',m,m);
    ncwrite([fpath_rfrom fname_sla_rfrom 'CLIM.nc'],'sla',...
        RFROM_sla_clim,[1 1 m]);
    % clean up
    clear RFROM_sla_clim a m v y vars
end

% clean up
clear

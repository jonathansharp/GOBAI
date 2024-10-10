% Glob_climatology_chlor_a.m
%
% DESCRIPTION:
% This script 
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 6/5/2024

%% establish paths, folders, and options
% file paths and names
fpath_chl = '/raid/Data/CMEMS/cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D';
fpath_rfrom = [pwd '/Data/RFROM/'];
folder_rfrom = 'RFROM_TEMP_v0.1';
fname_rfrom = '/RFROM_TEMP_STABLE_';
folder_chl_rfrom = 'RFROM_CHL_CMEMS';
fname_chl_rfrom = 'RFROM_CHL_';
% make new directory if it does not exist
if ~exist([fpath_rfrom folder_chl_rfrom],'dir')
    mkdir([fpath_rfrom folder_chl_rfrom])
end
% define timespan
y1 = 2004; y2 = 2021;
% pre-allocate variables
filenames = dir(fpath_chl);
chl_year = nan(length(filenames)-5,1);
chl_month = nan(length(filenames)-5,1);
chl_day = nan(length(filenames)-5,1);

%% get and extract dates from chl file list
for l = 4:length(filenames)-2
    if length(filenames(l).name) > 30
        chl_date = filenames(l).name(1:8);
        chl_year(l-3) = str2double(chl_date(1:4));
        chl_month(l-3) = str2double(chl_date(5:6));
        chl_day(l-3) = str2double(chl_date(7:8));
    else
        chl_date = NaN;
        chl_year(l-3) = NaN;
        chl_month(l-3) = NaN;
        chl_day(l-3) = NaN;
    end
end
file_datenum = datenum(chl_year,chl_month,chl_day);
clear filenames chl_date
% get schema from chl file
f_name_temp = '/2004/01/20040101_cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D.nc';
chl_schema = ncinfo([fpath_chl f_name_temp]);
% get dimensions from chl file
chl_lon = ncread([fpath_chl f_name_temp],'lon');
chl_lat = ncread([fpath_chl f_name_temp],'lat');
% get dimensions from RFROM file
rfrom_lon = ncread([fpath_rfrom folder_rfrom fname_rfrom '2004_01.nc'],'longitude');
rfrom_lat = ncread([fpath_rfrom folder_rfrom fname_rfrom '2004_01.nc'],'latitude');

%% loop through each year and month and save RFROM equivalent of chlorophyll
for y = y1:y2
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
        % alter schema for RFROM chl file
        rfrom_chl_schema = rfrom_schema;
        idx_rem = [];
        for v = 1:length(rfrom_chl_schema.Variables)
            if contains(rfrom_chl_schema.Variables(v).Name,'pressure')
                idx_rem = [idx_rem v];
            elseif contains(rfrom_chl_schema.Variables(v).Name,'temperature')
                idx_chl = v; 
            end
        end
        rfrom_chl_schema.Variables(idx_chl).Name = 'chlor_a';
        rfrom_chl_schema.Variables(idx_chl).Attributes(1).Value = 'mg m^-3';
        rfrom_chl_schema.Variables(idx_chl).Attributes(2).Value = 'mass_concentration_of_chlorophyll_in_sea_water';
        rfrom_chl_schema.Variables(idx_chl).Dimensions(3) = [];
        rfrom_chl_schema.Variables(idx_rem) = [];
        rfrom_chl_schema.Attributes = [];
        rfrom_chl_schema.Dimensions(4:5) = [];
        % delete NetCDF if it exists
        fname = [fpath_rfrom folder_chl_rfrom '/' fname_chl_rfrom ...
            num2str(y) '_' sprintf('%02s',num2str(m)) '.nc'];
        if exist(fname,'file'); delete(fname); end
        % write schema to NetCDF
        ncwriteschema(fname,rfrom_chl_schema);
        % save dimensions to NetCDF
        ncwrite(fname,'longitude',rfrom_lon);
        ncwrite(fname,'latitude',rfrom_lat);
        ncwrite(fname,'time',rfrom_time);
        % log weekly average chlorophyll values in the new file
        for w = 1:length(rfrom_datenum)
            % find seven days of the week
            index = find(file_datenum >= rfrom_datenum(w)-3 & file_datenum <= rfrom_datenum(w)+3);
            % pre-allocate daily chlorophyll
            chlor_a = nan(length(chl_lon),length(chl_lat),length(index));
            % log daily chlorophyll
            for d = 1:length(index)
                chlor_a(:,:,d) = ncread([fpath_chl '/' ...
                    num2str(chl_year(index(d))) '/' ...
                    sprintf('%02d',chl_month(index(d))) '/' ...
                    num2str(chl_year(index(d))) ...
                    sprintf('%02d',chl_month(index(d))) ...
                    sprintf('%02d',chl_day(index(d))) ...
                    '_cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D.nc'],...
                    'CHL');
            end
            % average into weekly chlorophyll
            chlor_a_w = mean(chlor_a,3,'omitnan');
            % convert RFROM longitude to -180 to 180 for interpolation
            rfrom_lon_new = convert_lon(rfrom_lon);
            % interpolate 4km weekly chlorophyll to RFROM grid
            rfrom_chlor_a = single(interp2(chl_lon',chl_lat,chlor_a_w',...
                rfrom_lon_new',rfrom_lat))'; clear rfrom_lon_new
            % test plots
            % figure; pcolor(chl_lon,chl_lat,chlor_a_w'); shading flat; colorbar; caxis([0 2]);
            % figure; pcolor(rfrom_lon,rfrom_lat,rfrom_chlor_a'); shading flat; colorbar; caxis([0 2]);
            % interpolate over gaps
            % tbd
            % test plot
            % figure; pcolor(rfrom_lon,rfrom_lat,rfrom_chlor_a_i'); shading flat; colorbar; caxis([0 2]);
            % save weekly chlorophyll to NetCDF
            ncwrite(fname,'chlor_a',rfrom_chlor_a,[1 1 w]);
        end
    end
end

%% create monthly annual files
% make new directory if it does not exist
if ~exist([fpath_rfrom folder_chl_rfrom '_annual'],'dir')
    mkdir([fpath_rfrom folder_chl_rfrom '_annual'])
end
% loop through each year to create annual files
for y = y1:y2
    % pre-allocate annual monthly RFROMs
    RFROM_chl.(['m' num2str(y)]) = [];
    RFROM_chl.(['m' num2str(y) '_time']) = [];
    % extract weekly chlorophyll values and average to monthly
    for m = 1:12
        fname = [fpath_rfrom folder_chl_rfrom '/' fname_chl_rfrom ...
            num2str(y) '_' sprintf('%02s',num2str(m)) '.nc'];
        RFROM_chl_temp = ncread(fname,'chlor_a');
        RFROM_chl.(['m' num2str(y)]) = ...
            cat(3,RFROM_chl.(['m' num2str(y)]),mean(RFROM_chl_temp,3,'omitnan'));
        RFROM_chl_temp = ncread([fpath_rfrom folder_chl_rfrom '/' fname_chl_rfrom num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'time');
        RFROM_chl.(['m' num2str(y) '_time']) = ...
            [RFROM_chl.(['m' num2str(y) '_time']);mean(RFROM_chl_temp,'omitnan')];
        clear RFROM_temp
    end
    % copy schema to new annual file
    schema = ncinfo([fpath_rfrom folder_chl_rfrom '/' fname_chl_rfrom num2str(y) '_01.nc']);
    schema.Dimensions(3).Length = 12;
    schema.Variables(3).Size = 12;
    schema.Variables(3).Dimensions.Length = 12;
    schema.Variables(3).Attributes(1).Value = 'month of year';
    schema.Variables(3).Attributes(2).Value = 'Number of Month, 1 (Jan.) to 12 (Dec.)';
    schema.Variables(4).Size = [1440,720,12];
    schema.Variables(4).Dimensions(3).Length = 12;
    if isfile([fpath_rfrom folder_chl_rfrom '_annual' '/' fname_chl_rfrom num2str(y) '.nc'])
        delete([fpath_rfrom folder_chl_rfrom '_annual' '/' fname_chl_rfrom num2str(y) '.nc']);
    end
    ncwriteschema([fpath_rfrom folder_chl_rfrom '_annual' '/' fname_chl_rfrom num2str(y) '.nc'],schema);
    clear schema
    % add dimensional variables to monthly file
    vars = {'longitude' 'latitude'};
    for v = 1:length(vars)
        a=ncread([fpath_rfrom folder_chl_rfrom '/' fname_chl_rfrom num2str(y) '_01.nc'],vars{v});
        ncwrite([fpath_rfrom folder_chl_rfrom '_annual' '/' fname_chl_rfrom num2str(y) '.nc'],vars{v},a);
    end
    % add averages to monthly file
    ncwrite([fpath_rfrom folder_chl_rfrom '_annual' '/' fname_chl_rfrom num2str(y) '.nc'],...
        'time',RFROM_chl.(['m' num2str(y) '_time']));
    ncwrite([fpath_rfrom folder_chl_rfrom '_annual' '/' fname_chl_rfrom num2str(y) '.nc'],...
        'chlor_a',RFROM_chl.(['m' num2str(y)]));
    % clean up
    clear RFROM_chl a m v y vars
end

%% create climatological file
% copy schema to new climatological file
schema = ncinfo([fpath_rfrom folder_chl_rfrom '_annual/' fname_chl_rfrom '2004.nc']);
if isfile([fpath_rfrom fname_chl_rfrom 'CLIM.nc'])
    delete([fpath_rfrom fname_chl_rfrom 'CLIM.nc']);
end
ncwriteschema([fpath_rfrom fname_chl_rfrom 'CLIM.nc'],schema);
clear schema
% add dimensional variables to monthly file
vars = {'longitude' 'latitude'};
for v = 1:length(vars)
    a=ncread([fpath_rfrom folder_chl_rfrom '_annual/' fname_chl_rfrom '2004.nc'],vars{v});
    ncwrite([fpath_rfrom fname_chl_rfrom 'CLIM.nc'],vars{v},a);
end
% clean up
clear a m v y vars

% loop through each month to add to climatological file
for m = 1:12
    RFROM_chl_clim = [];
    % extract monthly chlorophyll values
    for y = y1:y2
        RFROM_chl_temp = ncread([fpath_rfrom folder_chl_rfrom '_annual/' fname_chl_rfrom num2str(y) '.nc'],...
            'chlor_a',[1 1 m],[Inf Inf 1]);
        RFROM_chl_clim = ...
            cat(3,RFROM_chl_clim,RFROM_chl_temp);
        clear RFROM_chl_temp
    end
    % average monthly values to climatological means
    RFROM_chl_clim = mean(RFROM_chl_clim,3,'omitnan');
    % add to NetCDF
    ncwrite([fpath_rfrom fname_chl_rfrom 'CLIM.nc'],'time',m,m);
    ncwrite([fpath_rfrom fname_chl_rfrom 'CLIM.nc'],'chlor_a',...
        RFROM_chl_clim,[1 1 m]);
    % clean up
    clear RFROM_chl_clim a m v y vars
end

% clean up
clear

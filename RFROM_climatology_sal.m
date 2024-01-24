% RFROM_climatology_sal.m
%
% DESCRIPTION:
% This function calculates annual monthly files and a
% monthly climatology from weekly RFROM salinity files
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/30/2023

%% establish paths, folders, and options
% file paths and names
fpath = [pwd '/Data/RFROM/'];
folder = 'RFROM_SAL_v0.1';
fname = '/RFROM_SAL_STABLE_';
% create annual file folder
if ~isfolder([fpath folder '_annual']); mkdir([fpath folder '_annual']); end
% timespan
y1 = 2004;
y2 = 2022;

%% create monthly annual files
% loop through each year to create annual files
for y = y1:y2
    % pre-allocate annual monthly RFROMs
    RFROM.(['m' num2str(y)]) = [];
    RFROM.(['m' num2str(y) '_time']) = [];
    % extract weekly salinity values and average to monthly
    for m = 1:12
        RFROM_temp = ncread([fpath folder fname num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'ocean_salinity');
        RFROM.(['m' num2str(y)]) = ...
            cat(4,RFROM.(['m' num2str(y)]),mean(RFROM_temp,4,'omitnan'));
        RFROM_temp = ncread([fpath folder fname num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'time');
        RFROM.(['m' num2str(y) '_time']) = ...
            [RFROM.(['m' num2str(y) '_time']);mean(RFROM_temp,'omitnan')];
        clear RFROM_temp
    end
    % copy schema to new annual file
    schema = ncinfo([fpath folder fname num2str(y) '_01.nc']);
    schema.Dimensions(3).Length = 12;
    schema.Variables(3).Size = 12;
    schema.Variables(3).Dimensions.Length = 12;
    schema.Variables(5).Size = [1440,720,58,12];
    schema.Variables(5).Dimensions(4).Length = 12;
    if isfile([fpath folder '_annual' fname num2str(y) '.nc'])
        delete([fpath folder '_annual' fname num2str(y) '.nc']);
    end
    ncwriteschema([fpath folder '_annual' fname num2str(y) '.nc'],schema);
    clear schema
    % add dimensional variables to monthly file
    vars = {'longitude' 'latitude' 'mean_pressure' 'mean_pressure_bnds'};
    for v = 1:length(vars)
        a=ncread([fpath folder fname num2str(y) '_01.nc'],vars{v});
        ncwrite([fpath folder '_annual' fname num2str(y) '.nc'],vars{v},a);
    end
    % add averages to monthly file
    ncwrite([fpath folder '_annual' fname num2str(y) '.nc'],...
        'time',RFROM.(['m' num2str(y) '_time']));
    ncwrite([fpath folder '_annual' fname num2str(y) '.nc'],...
        'ocean_salinity',RFROM.(['m' num2str(y)]));
    % clean up
    clear RFROM a m v y vars
end

%% create climatological file
% copy schema to new climatological file
schema = ncinfo([fpath folder '_annual' fname '2004.nc']);
schema.Variables(3).Attributes(1).Value = 'month of year';
schema.Variables(3).Attributes(2).Value = 'Number of Month, 1 (Jan.) to 12 (Dec.)';
if isfile([fpath fname 'CLIM.nc'])
    delete([fpath fname 'CLIM.nc']);
end
ncwriteschema([fpath fname 'CLIM.nc'],schema);
clear schema
% add dimensional variables to monthly file
vars = {'longitude' 'latitude' 'mean_pressure' 'mean_pressure_bnds'};
for v = 1:length(vars)
    a=ncread([fpath folder '_annual' fname '2004.nc'],vars{v});
    ncwrite([fpath fname 'CLIM.nc'],vars{v},a);
end
% clean up
clear a m v y vars

% loop through each month to add to climatological file
for m = 1:12
    RFROM_clim = [];
    % extract monthly salinity values
    for y = y1:y2
        RFROM_temp = ncread([fpath folder '_annual' fname num2str(y) '.nc'],...
            'ocean_salinity',[1 1 1 m],[Inf Inf Inf 1]);
        RFROM_clim = ...
            cat(4,RFROM_clim,RFROM_temp);
        clear RFROM_temp
    end
    % average monthly values to climatological means
    RFROM_clim = mean(RFROM_clim,4,'omitnan');
    % add to NetCDF
    ncwrite([fpath fname 'CLIM.nc'],'time',m,m);
    ncwrite([fpath fname 'CLIM.nc'],'ocean_salinity',...
        RFROM_clim,[1 1 1 m]);
    % clean up
    clear RFROM_clim a m v y vars
end

clear

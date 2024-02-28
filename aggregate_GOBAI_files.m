% aggregate_GOBAI_files.m
%
% DESCRIPTION:
% This function 
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 2/27/2023

%% establish paths, folders, and options
% file paths and names
fpath = [pwd '/Data/GOBAI/' dir_base];
% create annual and monthly file folders
if ~isfolder([fpath '/monthly']); mkdir([fpath '/monthly']); end
if ~isfolder([fpath '/annual']); mkdir([fpath '/annual']); end
% timespan
y1 = 2004;
y2 = 2022;

%% create monthly files for each year
% loop through each year to create annual files
for y = y1:y2
    % pre-allocate annual monthly RFROMs
    GOBAI.(['m' num2str(y)]) = [];
    GOBAI.(['m' num2str(y) '_time']) = [];
    % extract weekly temperature values and average to monthly
    for m = 1%:12
        % determine number of weeks in file
        weeks = length(dir([fpath '/m' num2str(m) '_*.nc']));
        % pre-allocalte for weekly grids
        GOBAI_temp = [];
        GOBAI_time_temp = [];
        % load weekly data
        for w = 1:weeks
            GOBAI_temp = cat(4,GOBAI_temp,...
                ncread([fpath '/m' num2str(m) '_w' num2str(w) '.nc'],'o2'));   
            GOBAI_time_temp = [GOBAI_time_temp;...
                datenum(y,m,(w-1)*7+3.5)-datenum(2004,1,1)];
        end
        % assemble monthly means 
        GOBAI.(['m' num2str(y)]) = ...
            cat(4,GOBAI.(['m' num2str(y)]),mean(GOBAI_temp,4,'omitnan'));
        GOBAI.(['m' num2str(y) '_time']) = ...
            [GOBAI.(['m' num2str(y) '_time']);mean(GOBAI_time_temp)];
        clear GOBAI_temp GOBAI_time_temp
    end
    % copy schema to new monthly file
    schema = ncinfo([fpath '/m' num2str(1) '_w' num2str(1) '.nc']);
    if isfile([fpath '/annual/y' num2str(y) '.nc'])
        delete([fpath '/annual/y' num2str(y) '.nc']);
    end
    ncwriteschema([fpath '/annual/y' num2str(y) '.nc'],schema);
    clear schema
    % add dimensional variables to monthly file
    vars = {'lon' 'lat' 'pres'};
    for v = 1:length(vars)
        a=ncread([fpath '/m' num2str(1) '_w' num2str(1) '.nc'],vars{v});
        ncwrite([fpath '/annual/y' num2str(y) '.nc'],vars{v},a);
    end
    % add averages to monthly file
    nccreate([fpath '/annual/y'  num2str(y) '.nc'],'month','Datatype',...
        'single','Dimensions',{'time' 12},'FillValue',NaN);
    ncwrite([fpath '/annual/y'  num2str(y) '.nc'],...
        'month',GOBAI.(['m' num2str(y) '_time']));
    ncwrite([fpath '/annual/y'  num2str(y) '.nc'],...
        'o2',GOBAI.(['m' num2str(y)]));
    % clean up
    clear GOBAI a m v y vars
end

%% create annual means file

%% create climatological file
% % copy schema to new climatological file
% schema = ncinfo([fpath folder '_annual' fname '2004.nc']);
% schema.Variables(3).Attributes(1).Value = 'month of year';
% schema.Variables(3).Attributes(2).Value = 'Number of Month, 1 (Jan.) to 12 (Dec.)';
% if isfile([fpath fname 'CLIM.nc'])
%     delete([fpath fname 'CLIM.nc']);
% end
% ncwriteschema([fpath fname 'CLIM.nc'],schema);
% clear schema
% % add dimensional variables to monthly file
% vars = {'longitude' 'latitude' 'mean_pressure' 'mean_pressure_bnds'};
% for v = 1:length(vars)
%     a=ncread([fpath folder '_annual' fname '2004.nc'],vars{v});
%     ncwrite([fpath fname 'CLIM.nc'],vars{v},a);
% end
% % clean up
% clear a m v y vars
% 
% % loop through each month to add to climatological file
% for m = 1:12
%     RFROM_clim = [];
%     % extract monthly temperature values
%     for y = y1:y2
%         RFROM_temp = ncread([fpath folder '_annual' fname num2str(y) '.nc'],...
%             'ocean_temperature',[1 1 1 m],[Inf Inf Inf 1]);
%         RFROM_clim = ...
%             cat(4,RFROM_clim,RFROM_temp);
%         clear RFROM_temp
%     end
%     % average monthly values to climatological means
%     RFROM_clim = mean(RFROM_clim,4,'omitnan');
%     % add to NetCDF
%     ncwrite([fpath fname 'CLIM.nc'],'time',m,m);
%     ncwrite([fpath fname 'CLIM.nc'],'ocean_temperature',...
%         RFROM_clim,[1 1 1 m]);
%     % clean up
%     clear RFROM_clim a m v y vars
% end

clear

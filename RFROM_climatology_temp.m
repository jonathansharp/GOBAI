% RFROM_climatology_temp.m
%
% DESCRIPTION:
% This function calculates annual monthly files and a
% monthly climatology from weekly RFROM temperature files
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 5/29/2025

%% establish paths, folders, and options
% file paths and names
fpath = '/fast2/temperature/';
folder = 'RFROM_TEMP_v2.2';
fname = '/RFROMV22_TEMP_STABLE_';
% create annual file folder
if ~isfolder([fpath folder '_annual']); mkdir([fpath folder '_annual']); end
% timespan
y1 = 1993;
y2 = 2024;

%% create monthly annual files
% loop through each year to create annual files
for y = y1:y2
    % pre-allocate annual monthly RFROMs
    RFROM.(['m' num2str(y)]) = [];
    RFROM.(['m' num2str(y) '_time']) = [];
    % extract weekly temperature values and average to monthly
    for m = 1:12
        RFROM_temp = ncread([fpath folder fname num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'ocean_temperature');
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
        'ocean_temperature',RFROM.(['m' num2str(y)]));
    % clean up
    clear RFROM a m v y vars
end

%% create climatological file
% copy schema to new climatological file
schema = ncinfo([fpath folder '_annual' fname '2004.nc']);
schema.Variables(3).Attributes(1).Value = 'month of year';
schema.Variables(3).Attributes(2).Value = 'Number of Month, 1 (Jan.) to 12 (Dec.)';
if isfile([fpath fname 'CLIM_' num2str(y1) '_' num2str(y2) '.nc'])
    delete([fpath fname 'CLIM_' num2str(y1) '_' num2str(y2) '.nc']);
end
ncwriteschema([fpath fname 'CLIM_' num2str(y1) '_' num2str(y2) '.nc'],schema);
clear schema
% add dimensional variables to monthly file
vars = {'longitude' 'latitude' 'mean_pressure' 'mean_pressure_bnds'};
for v = 1:length(vars)
    a=ncread([fpath folder '_annual' fname '2004.nc'],vars{v});
    ncwrite([fpath fname 'CLIM_' num2str(y1) '_' num2str(y2) '.nc'],vars{v},a);
end
% clean up
clear a m v y vars

% loop through each month to add to climatological file
for m = 1:12
    RFROM_clim = [];
    % extract monthly temperature values
    for y = y1:y2
        RFROM_temp = ncread([fpath folder '_annual' fname num2str(y) '.nc'],...
            'ocean_temperature',[1 1 1 m],[Inf Inf Inf 1]);
        RFROM_clim = ...
            cat(4,RFROM_clim,RFROM_temp);
        clear RFROM_temp
    end
    % average monthly values to climatological means
    RFROM_clim = mean(RFROM_clim,4,'omitnan');
    % add to NetCDF
    ncwrite([fpath fname 'CLIM_' num2str(y1) '_' num2str(y2) '.nc'],'time',m,m);
    ncwrite([fpath fname 'CLIM_' num2str(y1) '_' num2str(y2) '.nc'],'ocean_temperature',...
        RFROM_clim,[1 1 1 m]);
    % clean up
    clear RFROM_clim a m v y vars
end

clear





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old code
% 
%     % extract weekly temperature values and average to monthly
%     for m = 1:12
%         RFROM_temp = ncread([fpath fname num2str(y) '_' ...
%             sprintf('%02s',num2str(m)) '.nc'],'ocean_temperature');
%         RFROM.(['y' num2str(y)]) = ...
%             cat(4,RFROM.(['y' num2str(y)]),RFROM_temp);
%         clear RFROM_temp
%     end
%     % interpolate to 
%     
%     
%     
%     
%     % extract weekly temperature values and average to monthly
%     for m = 1:12
%         RFROM_temp = ncread([fpath fname num2str(y) '_' ...
%             sprintf('%02s',num2str(m)) '.nc'],'ocean_temperature');
%         RFROM.(['m' num2str(y)]) = ...
%             cat(4,RFROM.(['m' num2str(y)]),mean(RFROM_temp,4,'omitnan'));
%         RFROM_temp = ncread([fpath fname num2str(y) '_' ...
%             sprintf('%02s',num2str(m)) '.nc'],'time');
%         RFROM.(['m' num2str(y) '_time']) = ...
%             [RFROM.(['m' num2str(y) '_time']);mean(RFROM_temp,'omitnan')];
%         clear RFROM_temp
%     end
% 
%     %% add monthly values to copied file
%     ncwrite([fpath '_annual' fname num2str(y) '.nc'],...
%         'ocean_temperature',RFROM.(['m' num2str(y)]));
%     ncwrite([fpath '_annual' fname num2str(y) '.nc'],...
%         'time',RFROM.(['m' num2str(y) '_time']),'Dimensions');
% 
% 
% 
%     % extract dimensions
%     RFROM.longitude = ncread([fpath fname num2str(y) '_' ...
%             sprintf('%02s',num2str(m)) '.nc'],'longitude');
%     RFROM.latitude = ncread([fpath fname num2str(y) '_' ...
%             sprintf('%02s',num2str(m)) '.nc'],'latitude');
%     RFROM.mean_pressure = ncread([fpath fname num2str(y) '_' ...
%             sprintf('%02s',num2str(m)) '.nc'],'mean_pressure');
%     % export monthly means as NetCDFs
%     info = ncinfo([fpath fname num2str(y) '_' sprintf('%02s',num2str(m)) '.nc']);
%     sz = length(info.Variables);
%     % log latitude, longitude, and pressure
%     for v = 1:sz
%         var = info.Variables(v).Name;
%         if matches(var,{'longitude' 'latitude' 'mean_pressure'})
%             lgth = info.Variables(v).Dimensions.Length;
%             nccreate([fpath fname num2str(y) '.nc'],var,...
%                 'Dimensions',{var,lgth},'Datatype','single');
%             ncwrite([fpath fname num2str(y) '.nc'],var,...
%                 RFROM.(var));
%             for a = 1:length(info.Variables(v).Attributes)
%                 ncwriteatt([fpath fname num2str(y) '.nc'],var,...
%                     info.Variables(v).Attributes(a).Name,...
%                     info.Variables(v).Attributes(a).Value);
%             end
%         end
%     end
%     clear v a m var lgth sz info
%     % log month
%     nccreate([fpath fname num2str(y) '.nc'],'month',...
%         'Dimensions',{'month',12},'Datatype','single');
%     ncwrite([fpath fname num2str(y) '.nc'],'month',(1:12)');
%     % log monthly mean temperature
%     nccreate([fpath fname num2str(y) '.nc'],'ocean_temperature',...
%         'Dimensions',{'longitude','latitude','mean_pressure','month'},...
%         'Datatype','single');
%     ncwrite([fpath fname num2str(y) '.nc'],'ocean_temperature',...
%         RFROM.(['m' num2str(y)]));
%     clear RFROM
% end
% 
% clear fpath fname y

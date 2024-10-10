% MODIS_climatology_chlor_a.m
%
% DESCRIPTION:
% This script calculates annual monthly files and a
% monthly climatology from 8-day MODIS chlorophyll files
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 6/4/2024

%% establish paths, folders, and options
% file paths and names
fpath_chl = '/raid/Data/Sat/MODIS/Aqua/Mapped/Daily/4km/chlor_a/';
fpath_rfrom = [pwd '/Data/RFROM/'];
folder_rfrom = 'RFROM_TEMP_v0.1';
fname_rfrom = '/RFROM_TEMP_STABLE_';
folder_chl_rfrom = 'RFROM_CHL/';
fname_chl_rfrom = 'RFROM_CHL_';
% make new directory if it does not exist
if ~exist([fpath_rfrom folder_chl_rfrom],'dir')
    mkdir([fpath_rfrom folder_chl_rfrom])
end
% define timespan
y1 = 2004; y2 = 2023;
% pre-allocate variables
filenames = dir(fpath_chl);
chl_date = cell(length(filenames)-2,1);
chl_year = nan(length(filenames)-2,1);
chl_month = nan(length(filenames)-2,1);
chl_day = nan(length(filenames)-2,1);

%% get and extract dates from chl file list
for l = 3:length(filenames)
    chl_date(l-2) = extractBetween(filenames(l).name,'AQUA_MODIS.','.L3');
    chl_year(l-2) = str2double(chl_date{l-2}(1:4));
    chl_month(l-2) = str2double(chl_date{l-2}(5:6));
    chl_day(l-2) = str2double(chl_date{l-2}(7:8));
end
file_datenum = datenum(chl_year,chl_month,chl_day);
clear filenames chl_date
% get schema from chl file
chl_schema = ncinfo([fpath_chl 'AQUA_MODIS.20040101.L3m.DAY.CHL.chlor_a.4km.nc']);
% get dimensions from chl file
chl_lon = ncread([fpath_chl 'AQUA_MODIS.20040101.L3m.DAY.CHL.chlor_a.4km.nc'],'lon');
chl_lat = ncread([fpath_chl 'AQUA_MODIS.20040101.L3m.DAY.CHL.chlor_a.4km.nc'],'lat');
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
        fname = [fpath_rfrom folder_chl_rfrom fname_chl_rfrom ...
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
                chlor_a(:,:,d) = ncread([fpath_chl 'AQUA_MODIS.' ...
                    num2str(chl_year(index(d))) ...
                    sprintf('%02d',chl_month(index(d))) ...
                    sprintf('%02d',chl_day(index(d))) ...
                    '.L3m.DAY.CHL.chlor_a.4km.nc'],'chlor_a');
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

% %% create monthly annual files
% % loop through each year to create annual files
% for y = y1:y2
%     % pre-allocate annual monthly RFROMs
%     MODIS.(['m' num2str(y)]) = [];
%     MODIS.(['m' num2str(y) '_time']) = [];
%     % extract 8-day chlorophyll values and average to monthly
%     for m = 1:12
%         RFROM_temp = ncread(,'ocean_salinity');
%         RFROM.(['m' num2str(y)]) = ...
%             cat(4,RFROM.(['m' num2str(y)]),mean(RFROM_temp,4,'omitnan'));
%         RFROM_temp = ncread([fpath folder fname num2str(y) '_' ...
%             sprintf('%02s',num2str(m)) '.nc'],'time');
%         RFROM.(['m' num2str(y) '_time']) = ...
%             [RFROM.(['m' num2str(y) '_time']);mean(RFROM_temp,'omitnan')];
%         clear RFROM_temp
%     end
%     % copy schema to new annual file
%     schema = ncinfo([fpath folder fname num2str(y) '_01.nc']);
%     schema.Dimensions(3).Length = 12;
%     schema.Variables(3).Size = 12;
%     schema.Variables(3).Dimensions.Length = 12;
%     schema.Variables(5).Size = [1440,720,58,12];
%     schema.Variables(5).Dimensions(4).Length = 12;
%     if isfile([fpath folder '_annual' fname num2str(y) '.nc'])
%         delete([fpath folder '_annual' fname num2str(y) '.nc']);
%     end
%     ncwriteschema([fpath folder '_annual' fname num2str(y) '.nc'],schema);
%     clear schema
%     % add dimensional variables to monthly file
%     vars = {'longitude' 'latitude' 'mean_pressure' 'mean_pressure_bnds'};
%     for v = 1:length(vars)
%         a=ncread([fpath folder fname num2str(y) '_01.nc'],vars{v});
%         ncwrite([fpath folder '_annual' fname num2str(y) '.nc'],vars{v},a);
%     end
%     % add averages to monthly file
%     ncwrite([fpath folder '_annual' fname num2str(y) '.nc'],...
%         'time',RFROM.(['m' num2str(y) '_time']));
%     ncwrite([fpath folder '_annual' fname num2str(y) '.nc'],...
%         'ocean_salinity',RFROM.(['m' num2str(y)]));
%     % clean up
%     clear RFROM a m v y vars
% end
% 
% %% create climatological file
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
%     % extract monthly salinity values
%     for y = y1:y2
%         RFROM_temp = ncread([fpath folder '_annual' fname num2str(y) '.nc'],...
%             'ocean_salinity',[1 1 1 m],[Inf Inf Inf 1]);
%         RFROM_clim = ...
%             cat(4,RFROM_clim,RFROM_temp);
%         clear RFROM_temp
%     end
%     % average monthly values to climatological means
%     RFROM_clim = mean(RFROM_clim,4,'omitnan');
%     % add to NetCDF
%     ncwrite([fpath fname 'CLIM.nc'],'time',m,m);
%     ncwrite([fpath fname 'CLIM.nc'],'ocean_salinity',...
%         RFROM_clim,[1 1 1 m]);
%     % clean up
%     clear RFROM_clim a m v y vars
% end

% clean up
clear

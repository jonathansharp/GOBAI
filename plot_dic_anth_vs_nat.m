%% file information
ver = 'HR-v1.1'; % version
var = 'DIC'; % variable
var2 = 'dic';
fpath = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
name1 = ['GOBAI-' var '-HR-v202602-1993-01.nc'];
fpath_temp = '/raid/Data/RFROM/RFROM_TEMP_v2.2-2025/';
fpath_sal = '/raid/Data/RFROM/RFROM_SAL_v2.2-2025/';
% download dimensions
GOBAI.lon = ncread([fpath 'monthly/' name1],'longitude');
GOBAI.lat = ncread([fpath 'monthly/' name1],'latitude');
GOBAI.pres = ncread([fpath 'monthly/' name1],'mean_pressure');
[GOBAI.lon3d,GOBAI.lat3d,GOBAI.pres3d] = ...
    ndgrid(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.depth3d = -gsw_z_from_p(GOBAI.pres3d,GOBAI.lat3d);

%% define files to use
% gobai
files = dir([fpath 'monthly/']); idx_files = false(length(files),1);
for f = 1:length(files)
    idx_files(f) = contains(files(f).name,'GOBAI');
end
files(~idx_files) = [];
% rfrom temp
files_temp = dir(fpath_temp); idx_files_temp = false(length(files_temp),1);
for f = 1:length(files_temp)
    idx_files_temp(f) = contains(files_temp(f).name,'_STABLE_');
end
files_temp(~idx_files_temp) = [];
% rfrom sal
files_sal = dir(fpath_sal); idx_files_sal = false(length(files_sal),1);
for f = 1:length(files_sal)
    idx_files_sal(f) = contains(files_sal(f).name,'_STABLE_');
end
files_sal(~idx_files_sal) = [];
% extract version date
vrs_date_str = char(extractBetween(files(1).name,'HR-','-1993'));

%% define coverage mask and calculate weights
if ~isfile([fpath 'GOBAI-' var '-HR-' vrs_date_str '-mask.nc'])
    GOBAI.mask = true(length(GOBAI.lon),length(GOBAI.lat),length(GOBAI.pres));
    for t = 1:length(files)
        inf = ncinfo([fpath 'monthly/' files(t).name]);
        for w = 1:inf.Dimensions(find(strcmp({inf.Dimensions.Name},'time'))).Length
            gobai_tmp = ncread([fpath 'monthly/' files(t).name],...
                var2,[1 1 1 w],[Inf Inf Inf 1]);
            GOBAI.mask(isnan(gobai_tmp)) = false;
        end
    end
    % calculate volumetric weights
    GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
    GOBAI.vol(~GOBAI.mask) = NaN;
    % save basin mask
    nccreate([fpath 'GOBAI-' var '-HR-' vrs_date_str '-mask.nc'],'mask','Datatype','int8',...
        'Dimensions',{'lon' length(GOBAI.lon) 'lat' length(GOBAI.lat) 'pres' length(GOBAI.pres)});
    ncwrite([fpath 'GOBAI-' var '-HR-' vrs_date_str '-mask.nc'],'mask',int8(GOBAI.mask));
    % save basin volume
    nccreate([fpath 'GOBAI-' var '-HR-' vrs_date_str '-vol.nc'],'vol','Datatype','double',...
        'Dimensions',{'lon' length(GOBAI.lon) 'lat' length(GOBAI.lat) 'pres' length(GOBAI.pres)});
    ncwrite([fpath 'GOBAI-' var '-HR-' vrs_date_str '-vol.nc'],'vol',GOBAI.vol);
else
    GOBAI.mask = ncread([fpath 'GOBAI-' var '-HR-' vrs_date_str '-mask.nc'],'mask');
    GOBAI.vol = ncread([fpath 'GOBAI-' var '-HR-' vrs_date_str '-vol.nc'],'vol');
end

%% plot

%% plot 
for t = 200%1:length(files)
    time_temp = ncread([fpath 'monthly/' files(t).name],'time');
    inf = ncinfo([fpath 'monthly/' files(t).name]);
    for w = 3%1:inf.Dimensions(find(strcmp({inf.Dimensions.Name},'time'))).Length
        gobai_dic = ncread([fpath 'monthly/' files(t).name],...
            var2,[1 1 1 w],[Inf Inf Inf 1]);
        gobai_time = ncread([fpath 'monthly/' files(t).name],'time',w,1);
        gobai_dic(~GOBAI.mask) = NaN;
        figure; worldmap([-90 90],[20 380]);
        set(gcf,'Position',[100 100 1000 600])
        pcolorm(double(GOBAI.lat),double(GOBAI.lon),...
            double(gobai_dic(:,:,1))');
        colorbar; clim([1800 2300]);
        % access t and s
        rfrom_temp = ncread([fpath_temp files_temp(t).name],...
            'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
        rfrom_sal = ncread([fpath_sal files_sal(t).name],...
            'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
        % apply trace
        Canth = nan(size(rfrom_sal));
        for d = 1:length(GOBAI.pres)
            % extract each layer
            lon_temp = GOBAI.lon3d(:,:,d);
            lat_temp = GOBAI.lat3d(:,:,d);
            depth_temp = GOBAI.depth3d(:,:,d);
            sal_temp = rfrom_sal(:,:,d);
            temp_temp = rfrom_temp(:,:,d);
            % define index
            idx_temp = ~isnan(gobai_dic(:,:,d)) & ~isnan(sal_temp) & ...
                ~isnan(temp_temp);
            % calculate antropogenic carbon
            Canth_temp = nan(size(sal_temp));
            Canth_temp(idx_temp) = ...
                TRACEv1([lon_temp(idx_temp) lat_temp(idx_temp) depth_temp(idx_temp)],...
                year(datenum(1950,1,1)+double(gobai_time)),...
                [sal_temp(idx_temp) temp_temp(idx_temp)],[1 2],1);
            Canth(:,:,d) = Canth_temp;
            % print results
            if d == 1
                fprintf(['Anthropogenic carbon calculated for ' ...
                    num2str(GOBAI.pres(d)) 'dbar, '])
            else
                fprintf([num2str(GOBAI.pres(d)) 'dbar, '])
            end
        end
        % plot Canth
        figure; worldmap([-90 90],[20 380]);
        set(gcf,'Position',[100 100 1000 600])
        pcolorm(double(GOBAI.lat),double(GOBAI.lon),...
            double(Canth(:,:,1))');
        colorbar; clim([1800 2300]);


    end
end
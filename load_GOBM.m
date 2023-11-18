function GOBM = load_GOBM(model,ext,vars)

% define path and filename
filepath = ['/raid/Data/RECCAP2/Ocean/Models/3D_All/' model '/'];
filename = [filepath vars{1} '_' model '_A_1_gr_' ext '.nc'];

% longitude
lon_names = {'lon';'LONGITUDE';'Lon'};
found_lon = 0;
i = 1;
while ~found_lon && i <= length(lon_names)
    try
        lon = ncread(filename, lon_names{i});
        found_lon = 1;
    catch
        i = i + 1;
    end
end
if ~found_lon
    error('lon variable not found in %s', filename)
end
lon(lon > 180) = lon(lon > 180) - 360;
GOBM.lon = lon;

% latitude
lat_names = {'lat';'LATITUDE';'Lat'};
found_lat = 0;
i = 1;
while ~found_lat && i <= length(lat_names)
    try
        GOBM.lat = ncread(filename, lat_names{i});
        found_lat = 1;
    catch
        i = i + 1;
    end
end
if ~found_lat
    error('lat variable not found in %s', filename)
end

% depth
depth_names = {'depth';'z_t';'z_l';'DEPTH';'Depth';'lev'};
found_depth = 0;
i = 1;
while ~found_depth && i <= length(depth_names)
    try
        GOBM.depth = ncread(filename, depth_names{i});
        name_depth = depth_names{i};
        found_depth = 1;
    catch
        i = i + 1;
    end
end
if ~found_depth
    % for CNRM-ESM1-1 model, which has depth data in separate file
    % (created from CMIP6 output)
    filename_depth = [filepath, 'lev.nc'];
    GOBM.depth = ncread(filename_depth, 'lev');
    unit_depth = ncreadatt(filename_depth, 'lev', 'units');
else
    try
        unit_depth = ncreadatt(filename, name_depth, 'units');
    catch
        unit_depth = ncreadatt(filename, name_depth, 'unit');
    end
end
if strncmpi(unit_depth, 'centimeter', 10) || ...
    strncmpi(unit_depth, 'cm', 2)
    GOBM.depth = 0.01 * GOBM.depth; % convert to m
end

% time
time_names = {'time'; 'time_ann';'TIME';'Time'};
found_time = 0;
i = 1;
while ~found_time && i <= length(time_names)
    try
        GOBM.time = ncread(filename, time_names{i});
        GOBM.unit_time = ncreadatt(filename, time_names{i}, 'units');
        found_time = 1;
    catch
        i = i + 1; % try the next name
    end
end
if ~found_time
    error('time variable not found in %s', filename)
end
nt = length(GOBM.time);
% FIXME assume that there are between ~38 and ~60 years in each run
if nt > 100
    ts_per_year = 12; % monthly data
else
    ts_per_year = 1; % annual data
end
% start and count values can be used for monthly and annual data
if GOBM.time(end)<700000
    mtime = length(GOBM.time)-1;
    yr    = 1980:1:1980+mtime;
else % CCSM WHOI is monthly
    yr    = str2num(datestr(GOBM.time,'yyyy'));
end

% temporal index
if strcmp(model,'CESM-ETHZ') || strcmp(model,'MOM6-Princeton')
    idx_t = find(datenum(1980,0,GOBM.time) > datenum(2004,0,0) & ...
        datenum(1980,0,GOBM.time) < datenum(2018,0,0));
elseif strcmp(model,'CCSM-WHOI')
    idx_t = 553:720;
else
    idx_t = find(datenum(1980,0,GOBM.time) > datenum(2004,0,0) & ...
        datenum(1980,0,GOBM.time) < datenum(2018,0,0));
end

% load variables
for v = 1:length(vars)
    % define path and filename
    filename = [filepath vars{v} '_' model '_A_1_gr_' ext '.nc'];
    % load variable
    GOBM.(vars{v}) = ncread(filename,vars{v});
    % calculate trends for oxygen
    if strcmp(vars{v},'o2')
        GOBM.([vars{v} '_tr']) = nan(length(GOBM.lon),length(GOBM.lat),length(GOBM.depth));
        for a = 1:length(GOBM.lon)
            for b = 1:length(GOBM.lat)
                for c = 1:length(GOBM.depth)
                    if any(isnan(squeeze(GOBM.(vars{v})(a,b,c,idx_t))))
                        GOBM.([vars{v} '_tr'])(a,b,c) = NaN;
                    else
                        [~,~,x] = leastsq2(GOBM.time(idx_t),squeeze(GOBM.(vars{v})(a,b,c,idx_t)),GOBM.time(1),0,0);
                        GOBM.([vars{v} '_tr'])(a,b,c) = x(2)*365.25;
                    end
                end
            end
        end
    end
    % IAV
    GOBM.([vars{v} '_iav']) = std(GOBM.(vars{v})(:,:,:,idx_t),[],4,'omitnan'); % IAV
    % mean
    GOBM.(vars{v}) = mean(GOBM.(vars{v})(:,:,:,idx_t),4,'omitnan'); % compute mean
end

% calculate density from GOBM to convert mol m-3 to umol kg-1
GOBM.pres = -gsw_p_from_z(repmat(permute(GOBM.depth,[3 2 1]),length(GOBM.lon),length(GOBM.lat)),...
    repmat(GOBM.lat',length(GOBM.lon),1,length(GOBM.depth)));
GOBM.sal_abs = single(gsw_SA_from_SP(GOBM.so,GOBM.pres,...
    repmat(GOBM.lon,1,length(GOBM.lat),length(GOBM.depth)),...
    repmat(GOBM.lat',length(GOBM.lon),1,length(GOBM.depth))));
GOBM.temp_cns = gsw_CT_from_pt(GOBM.sal_abs,GOBM.thetao);
GOBM.sigma = gsw_sigma0(GOBM.sal_abs,GOBM.temp_cns);
GOBM.dens = gsw_rho(GOBM.sal_abs,GOBM.temp_cns,GOBM.pres);
GOBM.temp = gsw_t_from_CT(GOBM.sal_abs,GOBM.temp_cns,GOBM.pres);

% adjusted concentrations if needed
for v = 1:length(vars)
    % define path and filename
    filename = [filepath vars{v} '_' model '_A_1_gr_' ext '.nc'];
    try
        this_unit = ncreadatt(filename,vars{v},'units');
    catch
        this_unit = ncreadatt(filename,vars{v},'unit');
    end
    % FIXME this should be more general
    if strcmp(this_unit,'mol m-3') || strcmp(this_unit,'mol/m3')
        GOBM.(vars{v}) = (GOBM.(vars{v})./GOBM.dens).*10^6;
        GOBM.([vars{v} '_iav']) = (GOBM.([vars{v} '_iav'])./GOBM.dens).*10^6;
        try
        GOBM.([vars{v} '_tr']) = (GOBM.([vars{v} '_tr'])./GOBM.dens).*10^6;
        catch
        end
    end
end


end

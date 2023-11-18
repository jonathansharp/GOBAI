% Process RFROM
fpath = '/raid/Data/RFROM/RFROM_TEMP_v0.1';
fname1 = '/RFROM_TEMP_STABLE_';

for y = 2004:2022
    % pre-allocate monthly RFROM
    RFROM.(['m' num2str(y)]) = [];
    RFROM.(['m' num2str(y) '_time']) = [];
    % create annual file
    copyfile([fpath fname1 num2str(y) '_01.nc'],...
        [fpath '_annual' fname1 num2str(y) '.nc']);
    % extract weekly temperature values and average to monthly
    for m = 1:12
        RFROM_temp = ncread([fpath fname1 num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'ocean_temperature');
        RFROM.(['m' num2str(y)]) = ...
            cat(4,RFROM.(['m' num2str(y)]),mean(RFROM_temp,4,'omitnan'));
        RFROM_temp = ncread([fpath fname1 num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'time');
        RFROM.(['m' num2str(y) '_time']) = ...
            [RFROM.(['m' num2str(y) '_time']);mean(RFROM_temp,'omitnan')];
        clear RFROM_temp
    end
    % add monthly values to copied file
    ncwrite([fpath '_annual' fname1 num2str(y) '.nc'],...
        'ocean_temperature',RFROM.(['m' num2str(y)]));
    ncwrite([fpath '_annual' fname1 num2str(y) '.nc'],...
        'time',RFROM.(['m' num2str(y) '_time']),'Dimensions',);



    % extract dimensions
    RFROM.longitude = ncread([fpath fname1 num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'longitude');
    RFROM.latitude = ncread([fpath fname1 num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'latitude');
    RFROM.mean_pressure = ncread([fpath fname1 num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'mean_pressure');
    % export monthly means as NetCDFs
    info = ncinfo([fpath fname1 num2str(y) '_' sprintf('%02s',num2str(m)) '.nc']);
    sz = length(info.Variables);
    % log latitude, longitude, and pressure
    for v = 1:sz
        var = info.Variables(v).Name;
        if matches(var,{'longitude' 'latitude' 'mean_pressure'})
            lgth = info.Variables(v).Dimensions.Length;
            nccreate([fpath fname1 num2str(y) '.nc'],var,...
                'Dimensions',{var,lgth},'Datatype','single');
            ncwrite([fpath fname1 num2str(y) '.nc'],var,...
                RFROM.(var));
            for a = 1:length(info.Variables(v).Attributes)
                ncwriteatt([fpath fname1 num2str(y) '.nc'],var,...
                    info.Variables(v).Attributes(a).Name,...
                    info.Variables(v).Attributes(a).Value);
            end
        end
    end
    clear v a m var lgth sz info
    % log month
    nccreate([fpath fname1 num2str(y) '.nc'],'month',...
        'Dimensions',{'month',12},'Datatype','single');
    ncwrite([fpath fname1 num2str(y) '.nc'],'month',(1:12)');
    % log monthly mean temperature
    nccreate([fpath fname1 num2str(y) '.nc'],'ocean_temperature',...
        'Dimensions',{'longitude','latitude','mean_pressure','month'},...
        'Datatype','single');
    ncwrite([fpath fname1 num2str(y) '.nc'],'ocean_temperature',...
        RFROM.(['m' num2str(y)]));
    clear RFROM
end

clear fpath fname1 y

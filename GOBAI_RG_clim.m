% define filenames
fpath = '/raid/Data/GOBAI-O2/v2.3/';
fname = 'GOBAI-O2-v2.3.nc';
% write schema to new file
info = ncinfo([fpath fname]);
info.Dimensions(4).Length = 12;
info.Variables(1).Size = [360 145 58 12];
info.Variables(1).Dimensions(4).Length = 12;
info.Variables(5).Size = 12;
info.Variables(5).Dimensions.Length = 12;
info.Variables(6).Size = [360 145 58 12];
info.Variables(6).Dimensions(4).Length = 12;
info.Variables(7).Size = [360 145 58 12];
info.Variables(7).Dimensions(4).Length = 12;
info.Variables(8).Size = [360 145 58 12];
info.Variables(8).Dimensions(4).Length = 12;
if isfile([fpath 'GOBAI-O2-clim-v2.3.nc']); delete([fpath 'GOBAI-O2-clim-v2.3.nc']); end
ncwriteschema([fpath 'GOBAI-O2-clim-v2.3.nc'],info);
% read o2
o2 = ncread([fpath fname],'oxy');
uncer = ncread([fpath fname],'uncer');
% calculate climatology
o2_clim = nan(360,145,58,12);
uncer_clim = nan(360,145,58,12);
for m = 1:12
    o2_clim(:,:,:,m) = mean(o2(:,:,:,m:12:end),4,'omitnan');
    uncer_clim(:,:,:,m) = mean(uncer(:,:,:,m:12:end),4,'omitnan');
end
ncwrite([fpath 'GOBAI-O2-clim-v2.3.nc'],'oxy',o2_clim);
ncwrite([fpath 'GOBAI-O2-clim-v2.3.nc'],'uncer',uncer_clim);
clear o2 uncer
% read temp
temp = ncread([fpath fname],'temp');
% calculate climatology
temp_clim = nan(360,145,58,12);
for m = 1:12
    temp_clim(:,:,:,m) = mean(temp(:,:,:,m:12:end),4,'omitnan');
end
ncwrite([fpath 'GOBAI-O2-clim-v2.3.nc'],'temp',temp_clim);
clear temp
% read sal
sal = ncread([fpath fname],'sal');
% calculate climatology
sal_clim = nan(360,145,58,12);
for m = 1:12
    sal_clim(:,:,:,m) = mean(sal(:,:,:,m:12:end),4,'omitnan');
end
ncwrite([fpath 'GOBAI-O2-clim-v2.3.nc'],'sal',sal_clim);
clear sal
% write dimensions
lat = ncread([fpath fname],'lat');
ncwrite([fpath 'GOBAI-O2-clim-v2.3.nc'],'lat',lat);
lon = ncread([fpath fname],'lon');
ncwrite([fpath 'GOBAI-O2-clim-v2.3.nc'],'lon',lon);
pres = ncread([fpath fname],'pres');
ncwrite([fpath 'GOBAI-O2-clim-v2.3.nc'],'pres',pres);
ncwrite([fpath 'GOBAI-O2-clim-v2.3.nc'],'time',(1:12)');
ncwriteatt([fpath 'GOBAI-O2-clim-v2.3.nc'],'time','units','month of year');

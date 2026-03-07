% define filenames
fpath = '/raid/Data/GOBAI-O2/v2.3/';
fname = 'GOBAI-O2-v2.3.nc';
% write schema to new file
info = ncinfo([fpath fname]);
years = info.Dimensions(4).Length/12;
info.Dimensions(4).Length = years;
info.Variables(1).Size = [360 145 58 years];
info.Variables(1).Dimensions(4).Length = years;
info.Variables(5).Size = years;
info.Variables(5).Dimensions.Length = years;
info.Variables(6).Size = [360 145 58 years];
info.Variables(6).Dimensions(4).Length = years;
info.Variables(7).Size = [360 145 58 years];
info.Variables(7).Dimensions(4).Length = years;
info.Variables(8).Size = [360 145 58 years];
info.Variables(8).Dimensions(4).Length = years;
if isfile([fpath 'GOBAI-O2-ann-v2.3.nc']); delete([fpath 'GOBAI-O2-ann-v2.3.nc']); end
ncwriteschema([fpath 'GOBAI-O2-ann-v2.3.nc'],info);
% read o2
o2 = ncread([fpath fname],'oxy');
uncer = ncread([fpath fname],'uncer');
% calculate annual means
o2_ann = nan(360,145,58,years);
uncer_ann = nan(360,145,58,years);
for y = 1:years
    o2_ann(:,:,:,y) = mean(o2(:,:,:,(y-1)*12+1:(y-1)*12+12),4,'omitnan');
    uncer_ann(:,:,:,y) = mean(uncer(:,:,:,(y-1)*12+1:(y-1)*12+12),4,'omitnan');
end
ncwrite([fpath 'GOBAI-O2-ann-v2.3.nc'],'oxy',o2_ann);
ncwrite([fpath 'GOBAI-O2-ann-v2.3.nc'],'uncer',uncer_ann);
clear o2 uncer
% read temp
temp = ncread([fpath fname],'temp');
% calculate annual means
temp_ann = nan(360,145,58,years);
for y = 1:years
    temp_ann(:,:,:,y) = mean(temp(:,:,:,(y-1)*12+1:(y-1)*12+12),4,'omitnan');
end
ncwrite([fpath 'GOBAI-O2-ann-v2.3.nc'],'temp',temp_ann);
clear temp
% read sal
sal = ncread([fpath fname],'sal');
% calculate annual means
sal_ann = nan(360,145,58,years);
for y = 1:years
    sal_ann(:,:,:,y) = mean(sal(:,:,:,(y-1)*12+1:(y-1)*12+12),4,'omitnan');
end
ncwrite([fpath 'GOBAI-O2-ann-v2.3.nc'],'sal',sal_ann);
clear sal
% write dimensions
lat = ncread([fpath fname],'lat');
ncwrite([fpath 'GOBAI-O2-ann-v2.3.nc'],'lat',lat);
lon = ncread([fpath fname],'lon');
ncwrite([fpath 'GOBAI-O2-ann-v2.3.nc'],'lon',lon);
pres = ncread([fpath fname],'pres');
ncwrite([fpath 'GOBAI-O2-ann-v2.3.nc'],'pres',pres);
ncwrite([fpath 'GOBAI-O2-ann-v2.3.nc'],'time',(2004:2004+years-1)');
ncwriteatt([fpath 'GOBAI-O2-ann-v2.3.nc'],'time','units','year');

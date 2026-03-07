% create GOBAI-O2-v2.3

% filename
fpath = '/raid/Data/GOBAI-O2/v2.3/';
new_fname = 'GOBAI-O2-v2.3.nc';
if isfile([fpath new_fname]); delete([fpath new_fname]); end
lon = ncread([fpath 'gobai-o2.nc'],'lon');
lat = ncread([fpath 'gobai-o2.nc'],'lat');
pres = ncread([fpath 'gobai-o2.nc'],'pres');
time = ncread([fpath 'gobai-o2.nc'],'time');
% o2
o2 = ncread([fpath 'gobai-o2.nc'],'o2');
nccreate([fpath new_fname],'oxy','Dimensions',...
    {'lon' length(lon) 'lat' length(lat) 'pres' length(pres) 'time' length(time)},...
    'Datatype','single');
ncwrite([fpath new_fname],'oxy',o2);
ncwriteatt([fpath new_fname],'oxy','long_name','Dissolved Oxygen Concentration');
ncwriteatt([fpath new_fname],'oxy','units','micromoles per kilogram');
ncwriteatt([fpath new_fname],'oxy','source','machine learning scheme described in Sharp et al. (2023)');
% lon
nccreate([fpath new_fname],'lon','Dimensions',{'lon'},'Datatype','single');
ncwrite([fpath new_fname],'lon',lon);
ncwriteatt([fpath new_fname],'lon','long_name','Longitude');
ncwriteatt([fpath new_fname],'lon','units','degrees_east');
ncwriteatt([fpath new_fname],'lon','source','grid');
% lat
nccreate([fpath new_fname],'lat','Dimensions',{'lat'},'Datatype','single');
ncwrite([fpath new_fname],'lat',lat);
ncwriteatt([fpath new_fname],'lat','long_name','Latitude');
ncwriteatt([fpath new_fname],'lat','units','degrees_north');
ncwriteatt([fpath new_fname],'lat','source','grid');
% pres
nccreate([fpath new_fname],'pres','Dimensions',{'pres'},'Datatype','single');
ncwrite([fpath new_fname],'pres',pres);
ncwriteatt([fpath new_fname],'pres','long_name','Pressure');
ncwriteatt([fpath new_fname],'pres','units','decibars');
ncwriteatt([fpath new_fname],'pres','source','grid');
% time
nccreate([fpath new_fname],'time','Dimensions',{'time'},'Datatype','single');
ncwrite([fpath new_fname],'time',time);
ncwriteatt([fpath new_fname],'time','long_name','Time');
ncwriteatt([fpath new_fname],'time','units','days since 1950-01-01');
ncwriteatt([fpath new_fname],'time','source','grid');
% uncer
uncer = ncread([fpath 'gobai-o2-uncer.nc'],'u_tot_o2');
nccreate([fpath new_fname],'uncer','Dimensions',...
    {'lon' length(lon) 'lat' length(lat) 'pres' length(pres) 'time' length(time)},...
    'Datatype','single');
ncwrite([fpath new_fname],'uncer',uncer);
ncwriteatt([fpath new_fname],'uncer','long_name','Total Uncertainty');
ncwriteatt([fpath new_fname],'uncer','units','micromoles per kilogram');
ncwriteatt([fpath new_fname],'uncer','source','produced by combining three separate uncertainty sources in quadrature, as described in Sharp et al. (2023)');
% temp
temp = ncread([fpath 'RG_Climatology_Temp.nc'],'Temperature');
nccreate([fpath new_fname],'temp','Dimensions',...
    {'lon' length(lon) 'lat' length(lat) 'pres' length(pres) 'time' length(time)},...
    'Datatype','single');
ncwrite([fpath new_fname],'temp',temp);
ncwriteatt([fpath new_fname],'temp','long_name','Temperature');
ncwriteatt([fpath new_fname],'temp','units','degrees Celcius');
ncwriteatt([fpath new_fname],'temp','source','Roemmich and Gilson, 2009 (https://sio-argo.ucsd.edu/RG_Climatology.html)');
% sal
sal = ncread([fpath 'RG_Climatology_Sal.nc'],'Salinity');
nccreate([fpath new_fname],'sal','Dimensions',...
    {'lon' length(lon) 'lat' length(lat) 'pres' length(pres) 'time' length(time)},...
    'Datatype','single');
ncwrite([fpath new_fname],'sal',sal);
ncwriteatt([fpath new_fname],'sal','long_name','Salinity');
ncwriteatt([fpath new_fname],'sal','units','N/A');
ncwriteatt([fpath new_fname],'sal','source','Roemmich and Gilson, 2009 (https://sio-argo.ucsd.edu/RG_Climatology.html)');
clear

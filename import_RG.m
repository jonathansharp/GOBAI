% import_RG
%
% DESCRIPTION:
% This function imports the Roemmich and Gilson climatology of temperature
% and salinity, calculates monthly values using the mean field and monthly
% anomalies, appends additional monthly extension values to the dataset,
% and saves the data file.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/12/2023

%% Import Roemmich-Gilson Climatology
fpath = 'Data/RG_CLIM/';
temp_file = 'RG_ArgoClim_Temperature_2019.nc';
sal_file = 'RG_ArgoClim_Salinity_2019.nc';
RG.latitude = ncread([fpath temp_file],'LATITUDE');
RG.longitude = ncread([fpath temp_file],'LONGITUDE');
RG.pressure = ncread([fpath temp_file],'PRESSURE');
RG.time = ncread([fpath temp_file],'TIME');
RG.temp_anom = single(ncread([fpath temp_file],'ARGO_TEMPERATURE_ANOMALY'));
RG.temp_mean = single(ncread([fpath temp_file],'ARGO_TEMPERATURE_MEAN'));
RG.sal_anom  = single(ncread([fpath sal_file],'ARGO_SALINITY_ANOMALY'));
RG.sal_mean  = single(ncread([fpath sal_file],'ARGO_SALINITY_MEAN'));
clear temp_file sal_file

%% Add temperature and salinity anomaly extensions
current_date = datevec(datetime('now'));
for y = 2019:current_date(1)
    for m = 1:12
        try
            % temperature
            temporary_temp = ...
                ncread([fpath 'RG_ArgoClim_' num2str(y) sprintf('%02d',m) ...
                    '_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
            RG.temp_anom = cat(4,RG.temp_anom,single(temporary_temp));
            clear temporary_temp
            % salinity
            temporary_sal = ...
                ncread([fpath 'RG_ArgoClim_' num2str(y) sprintf('%02d',m) ...
                    '_2019.nc'],'ARGO_SALINITY_ANOMALY');
            RG.sal_anom = cat(4,RG.sal_anom,single(temporary_sal));
            clear temporary_sal
            % time
            temporary_time = ...
                ncread([fpath 'RG_ArgoClim_' num2str(y) sprintf('%02d',m) ...
                    '_2019.nc'],'TIME');
            RG.time = cat(1,RG.time,temporary_time);
            clear temporary_time
        catch
        end
    end
end
clear y m current_date

%% Determine absolute temperatures and salinities from mean and anomalies
RG.temp = RG.temp_anom + RG.temp_mean;
RG = rmfield(RG,{'temp_anom' 'temp_mean'});
RG.sal = RG.sal_anom + RG.sal_mean;
RG = rmfield(RG,{'sal_anom' 'sal_mean'});

%% save RG file

% establish file names
fname_temp = [fpath 'RG_Climatology_Temp.nc'];
fname_sal = [fpath 'RG_Climatology_Sal.nc'];

% delete files if applicable
if isfile(fname_temp); delete(fname_temp); end
if isfile(fname_sal); delete(fname_sal); end

% convert time
time = datevec(datetime(2004,RG.time+0.5,15));
RG.time = single(datenum(time)-datenum(2004,1,1));

% create temperature file
nccreate(fname_temp,'Temperature','Dimensions',...
    {'Longitude' length(RG.longitude) 'Latitude' length(RG.latitude)...
    'Pressure' length(RG.pressure) 'Time' length(RG.time)},'Datatype','single');

% write variables to temperature file
ncwrite(fname_temp,'Temperature',RG.temp);
nccreate(fname_temp,'Longitude','Dimensions',{'Longitude'},'Datatype','single');
ncwrite(fname_temp,'Longitude',RG.longitude);
nccreate(fname_temp,'Latitude','Dimensions',{'Latitude'},'Datatype','single');
ncwrite(fname_temp,'Latitude',RG.latitude);
nccreate(fname_temp,'Pressure','Dimensions',{'Pressure'},'Datatype','single');
ncwrite(fname_temp,'Pressure',RG.pressure);
nccreate(fname_temp,'Time','Dimensions',{'Time'},'Datatype','single');
ncwrite(fname_temp,'Time',RG.time);

% create salinity file
nccreate(fname_sal,'Salinity','Dimensions',...
    {'Longitude' length(RG.longitude) 'Latitude' length(RG.latitude)...
    'Pressure' length(RG.pressure) 'Time' length(RG.time)},'Datatype','single');

% write variables to temperature file
ncwrite(fname_sal,'Salinity',RG.sal);
nccreate(fname_sal,'Longitude','Dimensions',{'Longitude'},'Datatype','single');
ncwrite(fname_sal,'Longitude',RG.longitude);
nccreate(fname_sal,'Latitude','Dimensions',{'Latitude'},'Datatype','single');
ncwrite(fname_sal,'Latitude',RG.latitude);
nccreate(fname_sal,'Pressure','Dimensions',{'Pressure'},'Datatype','single');
ncwrite(fname_sal,'Pressure',RG.pressure);
nccreate(fname_sal,'Time','Dimensions',{'Time'},'Datatype','single');
ncwrite(fname_sal,'Time',RG.time);

% Matlab
% if ~isfolder('Data'); mkdir('Data'); end
% save(['Data/RG_200401_' num2str(time(end,1)) sprintf('%02d',time(end,2)) '.mat'],'RG','-v7.3');

%% Plot T at 20 dbar
figure; worldmap([-90 90],[20 380]);
title('Annual mean at 20 dbar (RG09)','fontsize',16)
set(gcf,'Position',[617, 599, 820, 420])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(RG.latitude),double(RG.longitude),mean(RG.temp(:,:,3,:),4)');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar; caxis([0 30]);
colormap(cmocean('thermal',12));
c.Label.String = ['Temperature ' char(176) 'C'];
c.FontSize = 12;
c.TickLength = 0;
mlabel off; plabel off;
clear c land
if ~isfolder('Figures'); mkdir('Figures'); end
if ~isfolder('Figures/Surface_Plots'); mkdir('Figures/Surface_Plots'); end
exportgraphics(gcf,'Figures/Surface_Plots/temp_20_dbar_RG.png');
close

%% Plot S at 20 dbar
figure; worldmap([-90 90],[20 380]);
title('Annual mean at 20 dbar (RG09)','fontsize',16)
set(gcf,'Position',[617, 599, 820, 420])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(RG.latitude),double(RG.longitude),mean(RG.sal(:,:,3,:),4)');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar; caxis([32 38]);
colormap(cmocean('haline',12));
c.Label.String = 'Salinity';
c.FontSize = 12;
mlabel off; plabel off;
clear c land
if ~isfolder('Figures'); mkdir('Figures'); end
if ~isfolder('Figures/Surface_Plots'); mkdir('Figures/Surface_Plots'); end
exportgraphics(gcf,'Figures/Surface_Plots/sal_20_dbar_RG.png');
close

%% clean up
clear RG fpath time ans fname_sal fname_temp tst

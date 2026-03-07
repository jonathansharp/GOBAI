%% plot GOBAI-O2 data in 3D

% file paths
fname = 'Figures/3d_gobai_oxygen_animation_hr.gif';
gpath = '/raid/Data/GOBAI-O2/v1.0-HR/';

% create plot
f = figure('visible','off'); hold on; grid on;
set(f,'Position',[20 400 1000 600],'color','w');
ax = gca; ax.FontSize = 16; ax.Color = 'w'; ax.TickLength = [0.005 0.005];
%text(-235,-105,'J. Sharp, 2024','color',rgb('gray'));

% plot land
land = shaperead('landareas', 'UseGeoCoords', true);
l=geoshow(land,'FaceColor',rgb('grey'));

% modify axis limits
zlim([-500 0]); zticks([-500 -250 0]); zlabel({'0' '250' '500'});
xlim([-180 180]); xticks([-150 -90 -30 30 90 150]);
ylim([-90 90]); yticks([-80 -40 0 40 80]);
clim([0 400]);

% edit colorbar
c=colorbar;
c.Label.String = '[O_{2}] (\mumol kg^{-1})';
zlabel('Depth (dbar)')
colormap(cmocean('ice'));

% load dimensions
lon = ncread([gpath 'GOBAI-O2-v1.0-HR.nc'],'lon');
lat = ncread([gpath 'GOBAI-O2-v1.0-HR.nc'],'lat');
pres = ncread([gpath 'GOBAI-O2-v1.0-HR.nc'],'pres');
time = ncread([gpath 'GOBAI-O2-v1.0-HR.nc'],'time');

% process dimensions
lon = convert_lon(lon,'format','-180:180');
[lon,lon_idx] = sort(lon);
[lon3d,lat3d,pres3d] = ndgrid(lon,lat,pres);

% view coordinates
%xpos = 0;
%ypos = 0;
xpos = [repmat(0,1,104),-30/260:-30/260:-30,-30+30/260:30/260:30];
ypos = [repmat(90,1,104),90-30/260:-30/260:60,repmat(60,1,520)];
lat_idx = 50;

for t = 1:884
        % load gobai
        gobai = ncread([gpath 'GOBAI-O2-v1.0-HR.nc'],'o2',[1,1,1,t],[Inf,Inf,Inf,1]);
        gobai = gobai(lon_idx,:,:);
        % plot data
        view(xpos(t),ypos(t));
        % plot surface
        if t <= 540
            h1=surf(lon3d(:,:,1),lat3d(:,:,1),-pres3d(:,:,1),gobai(:,:,1),...
                'EdgeColor','none');
            h2=surf(squeeze(lon3d(1,:,:)),squeeze(lat3d(1,:,:)),...
                -squeeze(pres3d(1,:,:)),squeeze(gobai(1,:,:)),'EdgeColor','none');
            h3=surf(squeeze(lon3d(end,:,:)),squeeze(lat3d(end,:,:)),...
                -squeeze(pres3d(end,:,:)),squeeze(gobai(end,:,:)),'EdgeColor','none');
            h4=surf(squeeze(lon3d(:,1,:)),squeeze(lat3d(:,1,:)),...
                -squeeze(pres3d(:,1,:)),squeeze(gobai(:,1,:)),'EdgeColor','none');
            h5=surf(squeeze(lon3d(:,2,:)),squeeze(lat3d(:,2,:)),...
                -squeeze(pres3d(:,2,:)),squeeze(gobai(:,2,:)),'EdgeColor','none');
            h6=surf(squeeze(lon3d(:,3,:)),squeeze(lat3d(:,3,:)),...
                -squeeze(pres3d(:,3,:)),squeeze(gobai(:,3,:)),'EdgeColor','none');
        % plot with equatorial cutout
        elseif t <= 550
            h1=surf(lon3d(:,lat_idx:end,1),lat3d(:,lat_idx:end,1),-pres3d(:,lat_idx:end,1),gobai(:,lat_idx:end,1),...
                'EdgeColor','none');
            h2=surf(squeeze(lon3d(1,lat_idx:end,:)),squeeze(lat3d(1,lat_idx:end,:)),...
                -squeeze(pres3d(1,lat_idx:end,:)),squeeze(gobai(1,lat_idx:end,:)),'EdgeColor','none');
            h3=surf(squeeze(lon3d(end,lat_idx:end,:)),squeeze(lat3d(end,lat_idx:end,:)),...
                -squeeze(pres3d(end,lat_idx:end,:)),squeeze(gobai(end,lat_idx:end,:)),'EdgeColor','none');
            h4=surf(squeeze(lon3d(:,lat_idx,:)),squeeze(lat3d(:,lat_idx,:)),...
                -squeeze(pres3d(:,lat_idx,:)),squeeze(gobai(:,lat_idx,:)),'EdgeColor','none');
            h5=surf(squeeze(lon3d(:,lat_idx+1,:)),squeeze(lat3d(:,lat_idx+1,:)),...
                -squeeze(pres3d(:,lat_idx+1,:)),squeeze(gobai(:,lat_idx+1,:)),'EdgeColor','none');
        else
        end
        title(datestr(time(t)+datenum(1950,0,1)));
        % capture frame
        frame = getframe(f);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % write to file
        if t == 1
            imwrite(imind,cm,fname,'gif','Loopcount',inf,'DelayTime',0.075);
        else
            imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime',0.075);
        end
        % increase counter
        if t <= 540
            delete(h1); delete(h2); delete(h3); delete(h4); delete(h5); delete(h6);
        else
            delete(h1); delete(h2); delete(h3); delete(h4); delete(h5);
        end
        if t == 270
            delete(l);
            land = shaperead('landareas', 'UseGeoCoords', true);
            ant = land(1);
            land(1) = [];
            l = geoshow(land,'FaceColor',rgb('grey'));
            a = geoshow(ant,'FaceColor',rgb('grey'),'FaceAlpha',0.2);
        end
end
clear
close all

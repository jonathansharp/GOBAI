%% plot GOBAI-O2 data in 3D

% file paths
fname = 'Figures/3d_gobai_oxygen_animation.gif';

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
lon = ncread('/raid/Data/GOBAI-O2/v2.1/GOBAI-O2-v2.1.nc','lon');
lat = ncread('/raid/Data/GOBAI-O2/v2.1/GOBAI-O2-v2.1.nc','lat');
pres = ncread('/raid/Data/GOBAI-O2/v2.1/GOBAI-O2-v2.1.nc','pres');

% process dimensions
lon = [lon(161:end);lon(1:160)];
[lon3d,lat3d,pres3d] = ndgrid(lon,lat,pres);

% view coordinates
xpos = [repmat(-30/96,1,23),-30/96:-30/96:-30,-30:30/54:30];
ypos = [repmat(90-30/96,1,23),90-30/96:-30/96:60,repmat(60,1,109)];
lat_idx = 50;

cnt = 1;
for y = 2004:2022
    for m = 1:12
        % load gobai
        gobai = ncread('/raid/Data/GOBAI-O2/v2.1/GOBAI-O2-v2.1.nc','oxy',...
            [1,1,1,(y-2004)*12+m],[Inf,Inf,Inf,1]);
        gobai = [gobai(161:end,:,:);gobai(1:160,:,:)];
        % plot data
        view(xpos(cnt),ypos(cnt));
        if cnt <= 120
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
        else
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
        end
        title([num2str(m) '/' num2str(y)]);
        % capture frame
        frame = getframe(f);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % write to file
        if cnt == 1
            imwrite(imind,cm,fname,'gif','Loopcount',inf,'DelayTime',0.2);
        else
            imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime',0.2);
        end
        % increase counter
        if cnt <= 120
            delete(h1); delete(h2); delete(h3); delete(h4); delete(h5); delete(h6);
        else
            delete(h1); delete(h2); delete(h3); delete(h4); delete(h5);
        end
        if cnt == 60
            delete(l);
            land = shaperead('landareas', 'UseGeoCoords', true);
            ant = land(1);
            land(1) = [];
            l = geoshow(land,'FaceColor',rgb('grey'));
            a = geoshow(ant,'FaceColor',rgb('grey'),'FaceAlpha',0.2);
        end
        cnt = cnt+1;
    end
end
clear
close all

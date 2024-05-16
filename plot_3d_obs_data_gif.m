%% plot O2 observational data in 3D
float_only = 0;
accum = 1;
param = 'o2';

% load data
if float_only == 1
    if ~exist('float_data','var')
        if strcmp(param,'o2')
            load('O2/Data/processed_float_o2_data_Feb-2024_D.mat');
        elseif strcmp(param,'no3')
            load('NO3/Data/processed_float_no3_data_Feb-2024_D.mat');
        end
        date = datenum(float_data.YEAR,0,float_data.DAY);
    end
else
    if ~exist('all_data','var')
        if strcmp(param,'o2')
            load('O2/Data/processed_all_o2_data_Feb-2024_D.mat');
        elseif strcmp(param,'no3')
            load('NO3/Data/processed_all_no3_data_Feb-2024_D.mat');
        end
        date = datenum(all_data.year,0,all_data.day);
    end
end

% file paths
if float_only == 1
    if strcmp(param,'o2')
        if accum == 1
            fname = 'Figures/3d_float_obs_oxygen_accum_animation.gif';
        else
            fname = 'Figures/3d_float_obs_oxygen_animation.gif';
        end
    elseif strcmp(param,'no3')
        if accum == 1
            fname = 'Figures/3d_float_obs_nitrate_accum_animation.gif';
        else
            fname = 'Figures/3d_float_obs_nitrate_animation.gif';
        end
    end
else
    if strcmp(param,'o2')
        if accum == 1
            fname = 'Figures/3d_obs_oxygen_accum_animation.gif';
        else
            fname = 'Figures/3d_obs_oxygen_animation.gif';
        end
    elseif strcmp(param,'no3')
        if accum == 1
            fname = 'Figures/3d_obs_nitrate_accum_animation.gif';
        else
            fname = 'Figures/3d_obs_nitrate_animation.gif';
        end
    end
end

% create plot
f = figure('visible','off'); hold on; grid on;
set(f,'Position',[20 400 1000 600],'color','w');
ax = gca; ax.FontSize = 16; ax.Color = 'w'; ax.TickLength = [0.005 0.005];
%text(-235,-105,'J. Sharp, 2024','color',rgb('gray'));

% plot land
land = shaperead('landareas', 'UseGeoCoords', true);
l = geoshow(land,'FaceColor',rgb('grey'));

% modify axis limits
zlim([-500 0]); zticks([-500 -250 0]); zlabel({'0' '250' '500'}); 
xlim([-180 180]); xticks([-150 -90 -30 30 90 150]); 
ylim([-90 90]); yticks([-80 -40 0 40 80]);
if strcmp(param,'o2'); clim([0 400]); elseif strcmp(param,'no3'); clim([0 40]); end

% edit colorbar
c=colorbar;
if strcmp(param,'o2')
    c.Label.String = '[O_{2}] (\mumol kg^{-1})';
    colormap(cmocean('ice'));
elseif strcmp(param,'no3')
    c.Label.String = '[NO_{3}] (\mumol kg^{-1})';
    colormap(cmocean('solar'));
end
zlabel('Depth (dbar)')

% view coordinates
xpos = [repmat(-30/96,1,23),-30/96:-30/96:-30,-30:30/54:30];
ypos = [repmat(90-30/96,1,23),90-30/96:-30/96:60,repmat(60,1,109)];

cnt = 1;
for y = 2004:2022
    for m = 1:12
        % plot data
        idx = date >= datenum(y,m,1) & date < datenum(y,m+1,1);
        view(xpos(cnt),ypos(cnt));
        if float_only == 1
            if strcmp(param,'o2')
                h=scatter3(float_data.LON(idx),float_data.LAT(idx),...
                    -float_data.PRES(idx),5,float_data.OXY(idx),'.');
            elseif strcmp(param,'no3')
                h=scatter3(float_data.LON(idx),float_data.LAT(idx),...
                    -float_data.PRES(idx),5,float_data.NIT(idx),'.');
            end
        else
            if strcmp(param,'o2')
                h=scatter3(all_data.longitude(idx),all_data.latitude(idx),...
                -all_data.pressure(idx),5,all_data.oxygen(idx),'.');
            elseif strcmp(param,'no3')
                h=scatter3(all_data.longitude(idx),all_data.latitude(idx),...
                -all_data.pressure(idx),5,all_data.nitrate(idx),'.');
            end
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
        cnt = cnt+1;
        if accum == 1
            % allow to accumulate
        else
            delete(h);
        end
        if cnt == 60
            delete(l);
            land = shaperead('landareas', 'UseGeoCoords', true);
            ant = land(1);
            land(1) = [];
            l = geoshow(land,'FaceColor',rgb('grey'));
            a = geoshow(ant,'FaceColor',rgb('grey'),'FaceAlpha',0.2);
            clear land ant
        end
    end
end
clear
close all

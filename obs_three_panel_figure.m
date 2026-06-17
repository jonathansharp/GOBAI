%% plot O2 observational data in 3D
float_only = 0;
accum = 0;
param = 'no3';
param_dir = 'NO3';
param_name = 'nitrate';
data_date = 'Mar-2026';
sz1 = 10;
sz2 = 3;

% load data
if float_only == 1
    if ~exist('float_data','var')
        load([param_dir '/Data/processed_float_' param '_data_fg_' data_date '_D_A.mat']);
        date = datenum(float_data.YEAR,0,float_data.DAY);
    end
else
    if ~exist('all_data','var')
        load([param_dir '/Data/processed_all_' param '_data_fg_' data_date '_D_A.mat']);
        date = datenum(all_data.year,0,all_data.day);
    end
end

% process profiles
vars = {'type' 'id' 'latitude' 'longitude' 'time' 'year'};
for v = 1:length(vars)
    profiles_gld.(vars{v}) = [];
    profiles_flt.(vars{v}) = [];
end
profiles_gld.date = [];
profiles_flt.date = [];
profiles_gld.count = [];
profiles_flt.count = [];
count_date = [];
[~,c] = unique(all_data.id,'first');
unique_idx = false(size(all_data.id));
unique_idx(c) = true;
for y = 2000:2025
    for m = 1:12
        % plot data
        idx_gld = date >= datenum(y,m,1) & date < datenum(y,m+1,1) & ...
            all_data.type == 2 & unique_idx;
        idx_flt = date >= datenum(y,m,1) & date < datenum(y,m+1,1) & ...
            all_data.type == 1 & unique_idx;
        vars = {'type' 'id' 'latitude' 'longitude' 'time' 'year'};
        for v = 1:length(vars)
            profiles_gld.(vars{v}) = [profiles_gld.(vars{v});...
                all_data.(vars{v})(idx_gld)];
            profiles_flt.(vars{v}) = [profiles_flt.(vars{v});...
                all_data.(vars{v})(idx_flt)];
        end
        profiles_gld.date = [profiles_gld.date;date(idx_gld)];
        profiles_flt.date = [profiles_flt.date;date(idx_flt)];
        profiles_gld.count = [profiles_gld.count;sum(idx_gld)];
        profiles_flt.count = [profiles_flt.count;sum(idx_flt)];
        count_date = [count_date;datenum(y,m,15)];
    end
end
        

% animation path
if float_only == 1
    fname = ['Figures/3d_float_profiles_' param_name '_animation.gif'];
elseif float_only == 0
    fname = ['Figures/3d_profiles_' param_name '_animation.gif'];
end

%% create plot
f = figure('visible','off');
set(f,'Position',[20 400 1800 1000],'color','w');
tiledlayout(3,2,'TileSpacing','loose','Padding','compact');
clrs = get(groot,'FactoryAxesColorOrder');

% create first axis
ax1 = nexttile([2 1]);
box on; hold on;
ax1.FontSize = 16;
ax1.Color = 'w';
ax1.TickLength = [0.005 0.005];
%text(ax1,-235,-105,'J. Sharp, 2024','color',rgb('gray'));
land = shaperead('landareas', 'UseGeoCoords', true);
l1 = geoshow(ax1,land,'FaceColor',rgb('grey'));
% modify axis limits
zlim([-500 0]); zticks([-500 -250 0]);
xlim([-180 180]); xticks([-150 -90 -30 30 90 150]); 
ylim([-90 90]); yticks([-80 -40 0 40 80]);
if strcmp(param,'o2'); clim([0 400]);
elseif strcmp(param,'no3'); clim([0 40]); 
elseif strcmp(param,'ph'); clim([7.5 8.2]); end
zticklabels(ax1,{'500' '250' '0'});
zlabel('Depth (dbar)')

% create first axis
ax2 = nexttile([2 1]);
box on; hold on;
ax2.FontSize = 16;
ax2.Color = 'w';
ax2.TickLength = [0.005 0.005];
%text(ax2,-235,-105,'J. Sharp, 2024','color',rgb('gray'));
land = shaperead('landareas', 'UseGeoCoords', true);
l2 = geoshow(ax2,land,'FaceColor',rgb('grey'));
% modify axis limits
zlim([-500 0]); zticks([-500 -250 0]);
xlim([-180 180]); xticks([-150 -90 -30 30 90 150]); 
ylim([-90 90]); yticks([-80 -40 0 40 80]);
if strcmp(param,'o2'); clim([0 400]);
elseif strcmp(param,'no3'); clim([0 40]); 
elseif strcmp(param,'ph'); clim([7.5 8.2]); end
zlabel('Depth (dbar)');
zticklabels(ax2,{'0' '250' '500'});
if strcmp(param,'o2')
    title(ax2,'[O_{2}] profiles over time');
elseif strcmp(param,'no3')
     title(ax2,'[NO_{3}] profiles over time');
elseif strcmp(param,'dic')
     title(ax2,'DIC profiles over time');
end

% create third axis
ax3 = nexttile([1 2]);
box on; hold on;
ax3.FontSize = 16;
ax3.Color = 'w';
ax3.TickLength = [0.005 0.005];

% set view coordinates
xpos = [-5:-10/48:-30,-30:15/48:30];
ypos = [85:-10/48:60,repmat(60,1,193)];

cnt = 1;
for y = 2000:2025
    for m = 1:12
        % get indices
        idx_gld_data = date >= datenum(y,m,1) & date < datenum(y,m+1,1) & all_data.type == 2;
        idx_flt_data = date >= datenum(y,m,1) & date < datenum(y,m+1,1) & all_data.type == 1;
        idx_gld_prof = profiles_gld.date >= datenum(y,m,1) & profiles_gld.date < datenum(y,m+1,1);
        idx_flt_prof = profiles_flt.date >= datenum(y,m,1) & profiles_flt.date < datenum(y,m+1,1);
        % adjust axis 1 view
        view(ax1,xpos(cnt),ypos(cnt));
        if float_only == 1
            if strcmp(param,'o2')
                h=scatter3(float_data.LON(idx),float_data.LAT(idx),...
                    -float_data.PRES(idx),sz,float_data.OXY(idx),'.');
            elseif strcmp(param,'no3')
                h=scatter3(float_data.LON(idx),float_data.LAT(idx),...
                    -float_data.PRES(idx),sz,float_data.NIT(idx),'.');
            elseif strcmp(param,'ph')
                h=scatter3(float_data.LON(idx),float_data.LAT(idx),...
                    -float_data.PRES(idx),sz,float_data.PH(idx),'.');
            end
        else
                h11=scatter3(ax1,all_data.longitude(idx_gld_data),all_data.latitude(idx_gld_data),...
                    -all_data.pressure(idx_gld_data),sz1,clrs(4,:),'o','filled');
                h12=scatter(ax2,profiles_gld.longitude(idx_gld_prof),profiles_gld.latitude(idx_gld_prof),...
                    sz2,clrs(4,:),'o','filled');
                h21=scatter3(ax1,all_data.longitude(idx_flt_data),all_data.latitude(idx_flt_data),...
                    -all_data.pressure(idx_flt_data),sz1,clrs(5,:),'o','filled');
                h22=scatter(ax2,profiles_flt.longitude(idx_flt_prof),profiles_flt.latitude(idx_flt_prof),...
                    sz2,clrs(5,:),'o','filled');
                legend(ax2,{'GLODAP' 'Argo'});
                h31 = plot(ax3,count_date(1:cnt),profiles_gld.count(1:cnt),...
                    'LineWidth',2,'Color',clrs(4,:));
                h32 = plot(ax3,count_date(1:cnt),profiles_flt.count(1:cnt),...
                    'LineWidth',2,'Color',clrs(5,:));
                ylabel('Monthly Profiles');
                xlim([min(count_date) max(count_date)]);
                if strcmp(param,'o2'); ylim([0 2000]); end
                if strcmp(param,'no3'); ylim([0 1250]); end
                datetick('x','yyyy','keeplimits');
                legend(ax3,{'GLODAP Profiles' 'Argo Profiles (A and D mode)'},'Location','northwest');
        end
        % add titles
        title(ax1,datestr(datenum(y,m,15),'mmm-yyyy'));
        % capture frame
%         frame = getframe(f);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
        % write to file
        if cnt == 1
            % imwrite(imind,cm,fname,'gif','Loopcount',inf,'DelayTime',0.2);
            exportgraphics(f,fname);
        else
            % imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime',0.2);
            exportgraphics(f,fname,Append=true);
        end
        % increase counter
        cnt = cnt+1;
        if accum == 1
            % allow to accumulate
        else
            delete(h11);
            delete(h21);
            delete(h31);
            delete(h32);
        end
        if cnt == 60
            delete(l1);
            land = shaperead('landareas', 'UseGeoCoords', true);
            ant = land(1);
            land(1) = [];
            l1 = geoshow(ax1,land,'FaceColor',rgb('grey'));
            a1 = geoshow(ax1,ant,'FaceColor',rgb('grey'),'FaceAlpha',0.2);
            clear land ant
        end
        delete(h11);
        delete(h21);
    end
end
clear
close all

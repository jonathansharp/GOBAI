% This plots an animation for a specified year in the southern ocean,
% comparing between low-resolution and high-resolution products
function plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,show_fig,txt_align,stop_opt)

% define properties
param1 = {'TEMP' 'SAL' 'O2' 'NO3' 'DIC'};
param2 = {'Temp' 'Sal' 'o2' 'no3' 'dic'};
param3 = {'Temperature' 'Salinity' 'oxy' 'no3' 'dic'};
param4 = {'ocean_temperature' 'ocean_salinity' '' '' ''};

% parameter processing
if strcmp(param_to_plot,'temp')
    vrs1 = 'RG';
    vrs2 = 'v2.2-2025';
    vrs1_disp = 'RG09';
    vrs2_disp = 'RFROM-v2.3';
    type = 'phys';
    param_idx = 1;
elseif strcmp(param_to_plot,'sal')
    vrs1 = 'RG';
    vrs2 = 'v2.2-2025';
    vrs1_disp = 'RG09';
    vrs2_disp = 'RFROM-v2.3';
    type = 'phys';
    param_idx = 2;
elseif strcmp(param_to_plot,'o2')
    vrs1 = 'v2.3';
    vrs2 = 'v1.1-HR';
    vrs1_disp = 'GOBAI-v2.3';
    vrs2_disp = 'GOBAI-HR-v1.1';
    type = 'bgc';
    param_idx = 3;
elseif strcmp(param_to_plot,'no3')
    vrs1 = 'v1.0';
    vrs2 = 'v1.1-HR';
    vrs1_disp = 'GOBAI-v1.0';
    vrs2_disp = 'GOBAI-HR-v1.1';
    type = 'bgc';
    param_idx = 4;
elseif strcmp(param_to_plot,'dic')
    vrs1 = 'v1.0';
    vrs2 = 'v1.1-HR';
    vrs1_disp = 'GOBAI-v1.0';
    vrs2_disp = 'GOBAI-HR-v1.1';
    type = 'bgc';
    param_idx = 5;
end

% define folder names
if strcmp(type,'bgc')
    folder_old = ['/raid/Data/GOBAI-' param1{param_idx} '/' vrs1 '/'];
    folder_new = ['/raid/Data/GOBAI-' param1{param_idx} '/' vrs2 '/'];
elseif strcmp(type,'phys')
    folder_old = '/raid/sharp/matlab/GOBAI/Data/RG_CLIM/';
    folder_new = ['/raid/Data/RFROM/RFROM_' param1{param_idx} '_' vrs2 '/'];
end

% define file names
if strcmp(type,'bgc')
    file_old = ['GOBAI-' param1{param_idx} '-' vrs1 '.nc'];
    file_new = ['GOBAI-' param1{param_idx} '-' vrs2 '.nc'];
elseif strcmp(type,'phys')
    file_old = ['RG_Climatology_' param2{param_idx} '.nc'];
    file_new = ['RFROMV' vrs2(2) vrs2(4) '_' param1{param_idx} '_STABLE_'];
end

% define gif names
f_name = ['gobai_' param2{param_idx} '_animation_' vrs1 '_vs_' vrs2 '_' ...
    num2str(plot_depth) 'dbar_' num2str(lat_lims(1)) '_' num2str(lat_lims(2)) ...
    'N_' num2str(lon_lims(1)) '_' num2str(lon_lims(2)) 'S'];
g_name = [f_name '.gif']; v_name = [f_name '.avi'];

% process time for 1x1 file
if strcmp(type,'bgc')
    time_old = ncread([folder_old file_old],'time');
    time_old = datenum(1950,1,1+double(time_old));
    date_old = datevec(time_old);
    year_old = date_old(:,1);
    time_idx_old = find(year_old >= min(plot_year) & year_old <= max(plot_year));
    time_old = time_old(time_idx_old);
    date_old = date_old(time_idx_old,:);
    year_old = year_old(time_idx_old);
elseif strcmp(type,'phys')
    time_old = ncread([folder_old file_old],'Time');
    time_old = datenum(2004,1,1+double(time_old));
    date_old = datevec(time_old);
    year_old = date_old(:,1);
    time_idx_old = find(year_old >= min(plot_year) & year_old <= max(plot_year));
    time_old = time_old(time_idx_old);
    date_old = date_old(time_idx_old,:);
    year_old = year_old(time_idx_old);
end

% process time for high-res file
if strcmp(type,'bgc')
    time_new = ncread([folder_new file_new],'time');
    time_new = datenum(1950,1,1+double(time_new));
    date_new = datevec(time_new);
    year_new = date_new(:,1);
    time_idx_new = find(year_new >= min(plot_year) & year_new <= max(plot_year));
    time_new = time_new(time_idx_new);
    date_new = date_new(time_idx_new,:);
    year_new = year_new(time_idx_new);
elseif strcmp(type,'phys')
    files = dir(folder_new);
    time_new = [];
    for f = 1:length(files)
        if contains(files(f).name,['_' param1{param_idx} '_STABLE'])
            time_new = [time_new;ncread([folder_new files(f).name],'time')];
        end
    end
    time_new = datenum(1950,1,1+double(time_new));
    date_new = datevec(time_new);
    year_new = date_new(:,1);
    time_idx_new = find(year_new >= min(plot_year) & year_new <= max(plot_year));
    time_new = time_new(time_idx_new);
    date_new = date_new(time_idx_new,:);
    year_new = year_new(time_idx_new);
end

% load dimensions
if strcmp(type,'bgc')
    Longitude_old = ncread([folder_old file_old],'lon');
    Latitude_old = ncread([folder_old file_old],'lat');
    Pressure_old = ncread([folder_old file_old],'pres');
    Longitude_new = ncread([folder_new file_new],'lon');
    Latitude_new = ncread([folder_new file_new],'lat');
    Pressure_new = ncread([folder_new file_new],'pres');
elseif strcmp(type,'phys')
    Longitude_old = ncread([folder_old file_old],'Longitude');
    Latitude_old = ncread([folder_old file_old],'Latitude');
    Pressure_old = ncread([folder_old file_old],'Pressure');
    Longitude_new = ncread([folder_new file_new '1993_01.nc'],'longitude');
    Latitude_new = ncread([folder_new file_new '1993_01.nc'],'latitude');
    Pressure_new = ncread([folder_new file_new '1993_01.nc'],'mean_pressure');
end

% depth index
depth_idx_old = find(Pressure_old == plot_depth);
depth_idx_new = find(Pressure_new == plot_depth);

% establish video file
% v = VideoWriter(['Figures/' v_name]);

% plot data each month/week
t_new = 1; c_o = 1;
for t_old = 1:length(time_idx_old)
    
    % establish figure
    h = figure('color','w','visible',show_fig,'Position',[616 474 1600 800]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % load monthly gobai (old)
    nexttile;

    % load old parameter
    z_old = ncread([folder_old file_old],param3{param_idx},...
        [1 1 depth_idx_old time_idx_old(t_old)],[Inf Inf 1 1]);
   
    % make plot for 1x1 file
    c_o = make_plot(Latitude_old,Longitude_old,lat_lims,lon_lims,...
        z_old,time_old,val_lims,step_val,...
        t_old,param1{param_idx},reg,plot_depth,vrs1_disp,type,txt_align);

    % add text for version
%     if strcmp(type,'bgc')
%         text(-0.8,.95,['GOBAI-' vrs1],'FontSize',14,'HorizontalAlignment',txt_align);
%     elseif strcmp(type,'phys')
%         text(-0.8,.95,['RFROM-' vrs1],'FontSize',14,'HorizontalAlignment',txt_align);
%     end

    % make plot for high-resolution file
    nexttile;
    mnth_idx = find(date_new(:,1) == date_old(t_old,1) & ... % year
        date_new(:,2) == date_old(t_old,2)); % month
    for t_new = 1:length(mnth_idx)

        % load new parameter
        if strcmp(type,'bgc')
            z_new = ncread([folder_new file_new],param2{param_idx},...
                [1 1 depth_idx_new time_idx_new(mnth_idx(t_new))],[Inf Inf 1 1]);
        else
            z_new = ncread([folder_new file_new num2str(year_new(mnth_idx(t_new))) ...
                '_' sprintf('%02d',date_new(mnth_idx(t_new),2)) '.nc'],param4{param_idx},...
                [1 1 depth_idx_new t_new],[Inf Inf 1 1]);
        end

        % make plot for high-res file
        c_n = make_plot(Latitude_new,Longitude_new,lat_lims,lon_lims,...
            z_new,time_new,val_lims,step_val,...
            mnth_idx(t_new),param1{param_idx},...
            reg,plot_depth,vrs2_disp,type,txt_align);

        % capture frame and write video
        % open(v);
        % frame = getframe(h);
        % writeVideo(v,frame);
        % close(v);
        
        % stop if applicable
        if strcmp(stop_opt,'yes')
            keyboard
        end
    
        % write to gif file
        if t_old == 1 && t_new == 1
            exportgraphics(h,['Figures/' g_name]);
        else
            exportgraphics(h,['Figures/' g_name],Append=true);
        end
    
        % clear frame
        cla; delete(c_n);

    end

    % close figure
    close(h);

end

function c = make_plot(lat,lon,lat_lims,lon_lims,z,time,...
    val_lims,step_val,t,param,reg,depth,vrs_disp,type,txt_align)

    % make plot for 1x1 file
    if strcmp(reg,'so')
        m_proj('stereographic','lon',0,'lat',-90,'radius',60);
    elseif strcmp(reg,'other')
        m_proj('miller','lon',lon_lims,'lat',lat_lims);
    end
    % z = [gobai(~idx_20,:);gobai(idx_20,:)];
    z(z<val_lims(1)) = val_lims(1); z(z>val_lims(end)) = val_lims(end);
    m_contourf(double([lon;lon(end)+1]),double(lat),...
        double([z;z(end,:)])',val_lims(1):step_val:val_lims(end),...
        'linestyle','none');
    if strcmp(reg,'so')
        m_grid('linestyle','-','ytick',-90:20:90,...
            'xtick',-180:30:180,'xaxislocation','top');
    elseif strcmp(reg,'other')
        m_grid('linestyle','none','ytick',-90:10:90,...
            'xtick',0:10:360);
    end
    set(gca,'XAxisLocation','bottom')
    m_coast('patch',[0.9 0.9 0.9]);
    % define display date
    if mean(diff(time)) < 10 % check for HR vs LR
        disp_date = datestr(time(t),'dd-mmm-yyyy');
    else
        disp_date = datestr(time(t),'mmm-yyyy');
    end
    % add text to plot
    if strcmp(txt_align,'left'); lon_pos = min(lon_lims);
    elseif strcmp(txt_align,'right'); lon_pos = max(lon_lims);
    else; disp('"txt_align" must be "left" or "right"'); end
    if strcmp(reg,'so')
        m_text(0,-87,disp_date,'HorizontalAlignment','center','FontSize',12);
        text(-0.8,.95,vrs_disp,'FontSize',14,'HorizontalAlignment','right');
    elseif strcmp(reg,'other')
        m_text(lon_pos,max(lat_lims),{vrs_disp;[num2str(depth) ' dbars'];disp_date},...
            'HorizontalAlignment',txt_align,'VerticalAlignment','top',...
            'BackgroundColor','w','FontWeight','Bold','FontSize',16);
    end
    clim(gca,val_lims);
    c=colorbar('location','southoutside');
    if strcmp(param,'O2')
        c.Label.String = 'Oxygen (\mumol kg^{-1})';
        colormap(gca,cmocean('ice',length(val_lims(1):step_val:val_lims(end))-1));
    elseif strcmp(param,'NO3')
        c.Label.String = 'Nitrate (\mumol kg^{-1})';
        colormap(gca,cmocean('speed',length(val_lims(1):step_val:val_lims(end))-1));
    elseif strcmp(param,'DIC')
        c.Label.String = 'Dissolved Inorganic Carbon (\mumol kg^{-1})';
        colormap(gca,cmocean('matter',length(val_lims(1):step_val:val_lims(end))-1));
    elseif strcmp(param,'TEMP')
        c.Label.String = ['Conservative Temperature (' char(176) 'C)'];
        colormap(gca,cmocean('thermal',length(val_lims(1):step_val:val_lims(end))-1));
    elseif strcmp(param,'SAL')
        c.Label.String = 'Absolute Salinity';
        colormap(gca,cmocean('haline',length(val_lims(1):step_val:val_lims(end))-1));
    end
    c.TickLength = 0;
    set(gca,'FontSize',20);

end

end
% This plots an animation for a specified year in the southern ocean,
% comparing between physical and biogeochemical products
function plot_animation_phys_vs_bgc(plot_year,plot_depth,param_to_plot_phys,...
    param_to_plot_bgc,type,reg,quant,lat_lims,lon_lims,val_lims_1,step_val_1,...
    val_lims_2,step_val_2,show_fig,txt_align,stop_opt)

% define properties [T S O2 NO3 DIC]
param1 = {'TEMP' 'SAL' 'O2' 'NO3' 'DIC'};
param2 = {'Temp' 'Sal' 'o2' 'no3' 'dic'};
param3 = {'Temperature' 'Salinity' 'o2' 'no3' 'dic'};
param4 = {'ocean_temperature' 'ocean_salinity' '' '' ''};

% parameter processing
if strcmp(param_to_plot_phys,'temp')
    if strcmp(type,'lr')
        vrs_phys = 'RG';
    elseif strcmp(type,'hr')
        vrs_phys = 'v2.2-2025';
    end
    param_idx_phys = 1;
elseif strcmp(param_to_plot_phys,'sal')
    if strcmp(type,'lr')
        vrs_phys = 'RG';
    elseif strcmp(type,'hr')
        vrs_phys = 'v2.2-2025';
    end
    param_idx_phys = 2;
else
    disp('param_to_plot_phys must be "temp" or "sal"')
end

% parameter processing
if strcmp(param_to_plot_bgc,'o2')
    if strcmp(type,'lr')
        vrs_bgc = 'v2.3';
    elseif strcmp(type,'hr')
        vrs_bgc = 'v1.1-HR';
    end
    param_idx_bgc = 3;
elseif strcmp(param_to_plot_bgc,'no3')
    if strcmp(type,'lr')
        vrs_bgc = 'v1.0';
    elseif strcmp(type,'hr')
        vrs_bgc = 'v1.1-HR';
    end
    param_idx_bgc = 4;
elseif strcmp(param_to_plot_bgc,'dic')
    if strcmp(type,'lr')
        vrs_bgc = 'v1.0';
    elseif strcmp(type,'hr')
        vrs_bgc = 'v1.1-HR';
    end
    param_idx_bgc = 5;
else
    disp('param_to_plot_bgc must be "o2", "no3", or "dic"')
end

% define folder names
if strcmp(type,'lr')
    folder_phys = '/raid/sharp/matlab/GOBAI/Data/RG_CLIM/';
    folder_bgc = ['/raid/Data/GOBAI-' param1{param_idx_bgc} '/' vrs_bgc '/'];
elseif strcmp(type,'hr')
    folder_phys = ['/raid/Data/RFROM/RFROM_' param1{param_idx_phys} '_' vrs_phys '/'];
    folder_bgc = ['/raid/Data/GOBAI-' param1{param_idx_bgc} '/' vrs_bgc '/'];
end

% define file names
if strcmp(type,'lr')
    file_phys = ['RG_Climatology_' param2{param_idx_phys} '.nc'];
    file_bgc = ['GOBAI-' param1{param_idx_bgc} '-' vrs_bgc '.nc'];
elseif strcmp(type,'hr')
    file_phys = ['RFROMV' vrs_phys(2) vrs_phys(4) '_' param1{param_idx_phys} '_STABLE_'];
    file_bgc = ['GOBAI-' param1{param_idx_bgc} '-' vrs_bgc '.nc'];
end

% define gif names
f_name = ['gobai_' param2{param_idx_phys} '_' param2{param_idx_bgc} '_animation_' vrs_phys ...
    '_vs_' vrs_bgc '_' num2str(plot_depth) 'dbar_' num2str(lat_lims(1)) '_' num2str(lat_lims(2)) ...
    'N_' num2str(lon_lims(1)) '_' num2str(lon_lims(2)) 'S'];
g_name = [f_name '.gif']; v_name = [f_name '.avi'];

% process time for 1x1 file
if strcmp(type,'lr')
    % bgc
    time_bgc = ncread([folder_bgc file_bgc],'time');
    time_bgc = datenum(1950,1,1+double(time_bgc));
    date_bgc = datevec(time_bgc);
    year_bgc = date_bgc(:,1);
    time_idx_bgc = find(ismember(year_bgc,plot_year));
    time_bgc = time_bgc(time_idx_bgc);
    date_bgc = date_bgc(time_idx_bgc,:);
    year_bgc = year_bgc(time_idx_bgc);
    % physical
    time_phys = ncread([folder_phys file_phys],'Time');
    time_phys = datenum(2004,1,1+double(time_phys));
    date_phys = datevec(time_phys);
    year_phys = date_phys(:,1);
    time_idx_phys = find(ismember(year_phys,plot_year));
    time_phys = time_phys(time_idx_phys);
    date_phys = date_phys(time_idx_phys,:);
    year_phys = year_phys(time_idx_phys);
elseif strcmp(type,'hr')
    % bgc
    time_bgc = ncread([folder_bgc file_bgc],'time');
    time_bgc = datenum(1950,1,1+double(time_bgc));
    date_bgc = datevec(time_bgc);
    year_bgc = date_bgc(:,1);
    time_idx_bgc = find(ismember(year_bgc,plot_year));
    time_bgc = time_bgc(time_idx_bgc);
    date_bgc = date_bgc(time_idx_bgc,:);
    year_bgc = year_bgc(time_idx_bgc);
    % phys
    files = dir(folder_phys);
    time_phys = [];
    for f = 1:length(files)
        if contains(files(f).name,['_' param1{param_idx_phys} '_STABLE'])
            time_phys = [time_phys;ncread([folder_phys files(f).name],'time')];
        end
    end
    time_phys = datenum(1950,1,1+double(time_phys));
    date_phys = datevec(time_phys);
    year_phys = date_phys(:,1);
    time_idx_phys = find(ismember(year_phys,plot_year));
    time_phys = time_phys(time_idx_phys);
    date_phys = date_phys(time_idx_phys,:);
    year_phys = year_phys(time_idx_phys);
end

% load dimensions
if strcmp(type,'lr')
    Longitude_bgc = ncread([folder_bgc file_bgc],'lon');
    Latitude_bgc = ncread([folder_bgc file_bgc],'lat');
    Pressure_bgc = ncread([folder_bgc file_bgc],'pres');
    Longitude_phys = ncread([folder_phys file_phys],'Longitude');
    Latitude_phys = ncread([folder_phys file_phys],'Latitude');
    Pressure_phys = ncread([folder_phys file_phys],'Pressure');
elseif strcmp(type,'hr')
    Longitude_bgc = ncread([folder_bgc file_bgc],'lon');
    Latitude_bgc = ncread([folder_bgc file_bgc],'lat');
    Pressure_bgc = ncread([folder_bgc file_bgc],'pres');
    Longitude_phys = ncread([folder_phys file_phys '1993_01.nc'],'longitude');
    Latitude_phys = ncread([folder_phys file_phys '1993_01.nc'],'latitude');
    Pressure_phys = ncread([folder_phys file_phys '1993_01.nc'],'mean_pressure');
end

% depth index
depth_idx_phys = find(Pressure_phys == plot_depth);
depth_idx_bgc = find(Pressure_bgc == plot_depth);

% establish video file
% v = VideoWriter(['Figures/' v_name]);

% plot data each month/week
t_bgc = 1; c_o = 1;

% define month index
mnth_idx_phys = [];
for y = 1:length(plot_year)
    for m = 1:12
        weeks = sum(date_phys(:,1) == plot_year(y) & date_phys(:,2) == m);
        mnth_idx_phys = [mnth_idx_phys;(1:weeks)'];
    end
end

for t = 1:length(time_idx_phys)

    % establish figure
    h = figure('color','w','visible',show_fig,'Position',[616 474 1600 800]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % load monthly gobai (old)
    nexttile;

    % load physical parameter
    if strcmp(type,'lr')
        z_phys = ncread([folder_phys file_phys],param3{param_idx_phys},...
            [1 1 depth_idx_phys time_idx_phys(t)],[Inf Inf 1 1]);
    elseif strcmp(type,'hr')
        try
        z_phys = ncread([folder_phys file_phys num2str(year_phys(t)) ...
            '_' sprintf('%02d',date_phys(t,2)) '.nc'],param4{param_idx_phys},...
            [1 1 depth_idx_phys mnth_idx_phys(t)],[Inf Inf 1 1]);
        catch
            keyboard
        end
    end
   
    % make plot for physical file
    make_plot(Latitude_phys,Longitude_phys,lat_lims,lon_lims,z_phys,time_phys,...
        val_lims_1,step_val_1,t,param1{param_idx_phys},reg,plot_depth,...
        vrs_phys,'phys',quant,txt_align);

    % load monthly gobai (old)
    nexttile;

    % load bgc parameter
    if strcmp(quant,'anom')
        lat_idx = find(Latitude_bgc > lat_lims(1) & ...
            Latitude_bgc < lat_lims(end));
        lon_idx = find(Longitude_bgc > lon_lims(1) & ...
            Longitude_bgc < lon_lims(end));
        z_bgc = ncread([folder_bgc file_bgc],param3{param_idx_bgc},...
            [min(lon_idx) min(lat_idx) depth_idx_bgc time_idx_bgc(t)],...
            [length(lon_idx) length(lat_idx) 1 length(time_bgc)]);
        [nlon,nlat,npres,ntime] = size(z_bgc);
        z_bgc_temp = reshape(z_bgc,[],ntime)';
        ngrid = nlon*nlat*npres;
        A = [ones(ntime,1),time_bgc,...
            sin(2.*pi.*time_bgc./365.2425),cos(2.*pi.*time_bgc./365.2425),...
            sin(4.*pi.*time_bgc./365.2425),cos(4.*pi.*time_bgc./365.2425)];
        B = A \ z_bgc_temp;
        tr = reshape(B(2,:)',[nlon nlat]);
        z_bgc_fit = A * B;
        z_bgc_anom = z_bgc_temp - z_bgc_fit;
        z_bgc = reshape(z_bgc_anom,[nlon nlat ntime]);
    else
        z_bgc = ncread([folder_bgc file_bgc],param3{param_idx_bgc},...
            [1 1 depth_idx_bgc time_idx_bgc(t)],[Inf Inf 1 1]);
    end
   
    % make plot for bgc file
    make_plot(Latitude_bgc,Longitude_bgc,lat_lims,lon_lims,z_bgc,...
        time_bgc,val_lims_2,step_val_2,t,param1{param_idx_bgc},reg,plot_depth,...
        vrs_bgc,'bgc',quant,txt_align);

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
    if strcmp(type,'lr')
        if t == 1
            exportgraphics(h,['Figures/' g_name]);
            exportgraphics(h,['Figures/' g_name],Append=true);
            exportgraphics(h,['Figures/' g_name],Append=true);
            exportgraphics(h,['Figures/' g_name],Append=true);
        else
            exportgraphics(h,['Figures/' g_name],Append=true);
            exportgraphics(h,['Figures/' g_name],Append=true);
            exportgraphics(h,['Figures/' g_name],Append=true);
            exportgraphics(h,['Figures/' g_name],Append=true);
        end
    elseif strcmp(type,'hr')
        if t == 1
            exportgraphics(h,['Figures/' g_name]);
        else
            exportgraphics(h,['Figures/' g_name],Append=true);
        end
    end

    % close figure
    close(h);

end

    function c = make_plot(lat,lon,lat_lims,lon_lims,z,time,val_lims,...
    step_val,t,param,reg,depth,vrs,type,quant,txt_align)

    % make plot for 1x1 file
    if strcmp(reg,'so')
        m_proj('stereographic','lon',0,'lat',-90,'radius',60);
    elseif strcmp(reg,'other')
        m_proj('miller','lon',lon_lims,'lat',lat_lims);
    end
    % z = [gobai(~idx_20,:);gobai(idx_20,:)];
    if strcmp(quant,'anom')
        lon_idx = find(lon > lon_lims(1) & lon < lon_lims(end));
        lat_idx = find(lat > lat_lims(1) & lat < lat_lims(end));
        m_contourf(double([lon(lon_idx);lon(lon_idx(end))+1]),...
            double(lat(lat_idx)),double([z;z(end,:)])',...
            -(val_lims(end)-val_lims(1))/10:step_val/10:(val_lims(end)-val_lims(1))/10,...
            'linestyle','none');
    else
        z(z<val_lims(1)) = val_lims(1); z(z>val_lims(end)) = val_lims(end);
        m_contourf(double([lon;lon(end)+1]),double(lat),...
            double([z;z(end,:)])',val_lims(1):step_val:val_lims(end),...
            'linestyle','none');
    end
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
    % define text for version
    if strcmp(type,'bgc')
        % disp_vrs = ['GOBAI-' vrs];
        disp_vrs = 'GOBAI-HR-v1.1';
    elseif strcmp(type,'phys')
        % disp_vrs = ['RFROM-' vrs];
        disp_vrs = 'RFROM-v2.3';
    end
    % add text to plot
    if strcmp(txt_align,'left'); lon_pos = min(lon_lims);
    elseif strcmp(txt_align,'right'); lon_pos = max(lon_lims);
    else; disp('"txt_align" must be "left" or "right"'); end
    if strcmp(reg,'so')
        m_text(0,-87,disp_date,'HorizontalAlignment','center','FontSize',12);
        text(-0.8,.95,'FontSize',14,'HorizontalAlignment','right');
    elseif strcmp(reg,'other')
        m_text(lon_pos,max(lat_lims),{disp_vrs;[num2str(depth) ' dbars'];disp_date},...
            'HorizontalAlignment',txt_align,'VerticalAlignment','top',...
            'BackgroundColor','w','FontWeight','Bold','FontSize',16);
    end
    clim(val_lims);
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
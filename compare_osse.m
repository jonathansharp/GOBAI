% compare_osse_to_model
%
% DESCRIPTION:
% This function 
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/7/2024

function compare_osse(param,param1,fpath,base_grid,file_date,...
    float_file_ext,num_clusters,start_year,snap_date,rlz,grid_label,grid_type)

%% process date
date_str = num2str(snap_date);

%% create directory names
% gobai
gobai_filepath = [param1 '/Data/GOBAI/' base_grid '/AVG/c' ...
    num2str(num_clusters) '_' file_date float_file_ext '/gobai-' param '.nc'];
% cmip
path2 = ['_Omon_' base_grid '_'];
path3 = ['_' rlz '_' grid_label];
cmip_filepath = [fpath 'combined/' grid_type '/' param path2 ...
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
% delta
delta_filepath = [param1 '/Data/GOBAI/' base_grid '/DELTA/c' ...
    num2str(num_clusters) '_' file_date float_file_ext];

%% load dimensions
lat = ncread(gobai_filepath,'lat');
lon = ncread(gobai_filepath,'lon');
depth = ncread(gobai_filepath,'depth');
time = ncread(gobai_filepath,'time');

%% load and evaluate monthly reconstructed gobai in relation to cmip
gobai_mean = nan(length(time),1);
gobai_depth_mean = nan(length(time),length(depth));
cmip_mean = nan(length(time),1);
cmip_depth_mean = nan(length(time),length(depth));

%% loop through each month
for m = 1:12%length(time)

    %% load monthly output
    gobai = ncread(gobai_filepath,param,[1 1 1 m],[Inf Inf Inf 1]);
    cmip = ncread(cmip_filepath,param,[1 1 1 m],[Inf Inf Inf 1]);
    delta = gobai-cmip;

    %% calculate global means
    vol = single(calculate_volume(lat,lon,depth));
    idx_gobai = ~isnan(gobai);
    gobai_mean(m) = sum(gobai(idx_gobai).*vol(idx_gobai))./sum(vol(idx_gobai));
    idx_cmip = ~isnan(cmip);
    cmip_mean(m) = sum(cmip(idx_cmip).*vol(idx_cmip))./sum(vol(idx_cmip));

    %% calculate means by depth
    area = single(calculate_area(lat,lon));
    for p = 1:length(depth)
        gobai_temp = gobai(:,:,p);
        idx_gobai_temp = idx_gobai(:,:,p);
        gobai_depth_mean(m,p) = sum(gobai_temp(idx_gobai_temp).*...
            area(idx_gobai_temp))./sum(area(idx_gobai_temp));
        cmip_temp = cmip(:,:,p);
        idx_cmip_temp = idx_cmip(:,:,p);
        cmip_depth_mean(m,p) = sum(cmip_temp(idx_cmip_temp).*...
            area(idx_cmip_temp))./sum(area(idx_cmip_temp));
    end

    %% plot
    figure('visible','off');
    worldmap('world');
    set(gca,'fontsize',12);
    pcolorm(lat,lon,delta(:,:,1)');
    title(datestr(time(m),'mmm yyyy'));
    plot_land('map');
    c=colorbar;
    caxis([-50 50]);
    colormap(cmocean('balance'));
    c.Label.String = ['\Delta[O_{2}]_{(GOBAI - ' base_grid ')}'];
    mlabel off;
    plabel off;
    if ~isfolder([param1 '/Figures/' base_grid '/' rlz '_' grid_label])
        mkdir(['Figures/' base_grid '/' rlz '_' grid_label]);
    end
    export_fig(gcf,[param1 '/Figures/' base_grid '/' rlz '_' grid_label ...
        '/' datestr(time(m),'mmm-yyyy') '.png'],'-transparent');
    close

    %% save monthly differences
    if m == 1
        if ~exist(delta_filepath); mkdir(delta_filepath); end
        % create combined file
        ncsave_4d([delta_filepath '/delta_gobai-' param '.nc'],...
            {'lon' lon 'longitude' 'degrees east'},...
            {'lat' lat 'latitude' 'degrees north'},...
            {'depth' depth 'depth' 'meters'},...
            {'time' time(m) 'time' 'time'},...
            {['delta_' param] delta ' ' 'umol/kg'});
    else
        % append to combined file
        ncwrite([delta_filepath '/delta_gobai-' param '.nc'],'time',time(2),m);
        ncwrite([delta_filepath '/delta_gobai-' param '.nc'],['delta_' param],delta,[1 1 1 m]);
    end

end

%% save global means


%% plot timeseries
figure;
plot(time,cmip_mean(1:12),time,gobai_mean(1:12),'LineWidth',2);
legend({base_grid ['GOBAI-' param1 '_{(' base_grid ')}']});
export_fig(gcf,[param1 '/Figures/' base_grid '/' rlz '_' grid_label ...
        '/timeseries.png'],'-transparent'); close;

%% plot average profile
figure;
plot(mean(gobai_depth_mean(1:12,:)),depth,...
    mean(cmip_depth_mean(1:12,:)),depth,'linewidth',2);
set(gca,'YDir','reverse','XAxisLocation','top');
ylim([0 max(depth)]);
xlabel('[O_{2}] (\mumol kg^{-1})');
ylabel('Depth (m)');
legend({base_grid ['GOBAI-' param1 '_{(' base_grid ')}']},'Location','southeast');

%% plot average profile difference
figure; hold on;
plot(mean(gobai_depth_mean(1:12,:))-...
    mean(cmip_depth_mean(1:12,:)),depth,'linewidth',2);
plot([0 0],[min(depth) max(depth)],'k--');
set(gca,'YDir','reverse','XAxisLocation','top');
ylim([0 max(depth)]);
xlim([-3 3]);
xlabel('\Delta[O_{2}] (\mumol kg^{-1})');
ylabel('Depth (m)');


end
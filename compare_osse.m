% compare_osse
%
% DESCRIPTION:
% This function compares reconstructed fields from GOBAI OSSEs
% to full model gridded fields.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/7/2024

function compare_osse(param_props,fpath,base_grid,file_date,...
    float_file_ext,num_clusters,start_year,snap_date,train_ratio,...
    val_ratio,test_ratio,rlz)

%% process date
date_str = num2str(snap_date);

%% create directory names (multiple models)
gobai_filepath = ... % FFNN
    [fpath 'GOBAI/' base_grid '/FFNN/c' ...
    num2str(num_clusters) '_' file_date float_file_ext '/train' ...
    num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
    num2str(100*test_ratio) '/gobai-' param_props.file_name '.nc'];

%% combined filepath
% gobai
gobai_filepath = ...
    [fpath 'GOBAI/' base_grid '/FFNN/c' ...
    num2str(num_clusters) '_' file_date float_file_ext '/train' ...
        num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
        num2str(100*test_ratio) '/gobai-' param_props.file_name '.nc'];
% cmip
path2 = ['_Omon_' base_grid '_'];
path3 = ['_' rlz '_gr'];
cmip_filepath = [fpath 'combined/regridded/' param_props.file_name path2 ...
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
% delta
delta_filepath = [fpath 'GOBAI/' base_grid '/DELTA/c' ...
    num2str(num_clusters) '_' file_date float_file_ext];
% delta_filepath = [fpath 'GOBAI/' base_grid '/DELTA/' ...
%     'multiple_clusters_' file_date float_file_ext];

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
for m = 1:length(time)

    %% load monthly output
    gobai = ncread(gobai_filepath,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);
    cmip = ncread(cmip_filepath,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);

    %% apply RG mask
    % load RG grid
    lat_rg = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Latitude');
    lon_rg = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Longitude');
    temp_rg = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'],...
        'Temperature',[1 1 1 m],[Inf Inf 1 1]);
    mask_rg = (isnan(temp_rg));
    % reorder and pad mask
    mask_cmip = [mask_rg(341:end,:,:);mask_rg(1:340,:,:)];
    mask_cmip = [true(length(lon),25,length(depth)),...
        repmat(mask_cmip,1,1,length(depth)),true(length(lon),10,length(depth))];
    % apply mask to monthly output
    gobai(mask_cmip) = NaN;
    cmip(mask_cmip) = NaN;

    %% determine differences between masked grids
    delta = gobai-cmip;

    %% calculate global means
    vol = single(calculate_volume(lat,lon,depth));
    idx_gobai = ~isnan(gobai);
    % Eliminate bottom depth for MPI-ESM1-2-LR
    if strcmp(base_grid,'MPI-ESM1-2-LR'); idx_gobai(:,:,end) = false; end
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
    % figure('visible','on');
    % worldmap('world');
    % set(gca,'fontsize',12);
    % pcolorm(lat,lon,delta(:,:,1)');
    % title(datestr(datetime(1950,1,1)+time(m),'mmm yyyy'));
    % plot_land('map');
    % c=colorbar;
    % clim([-(max(param_props.edges)/10) (max(param_props.edges)/10)]);
    % colormap(cmocean('balance'));
    % c.Label.String = ['\Delta[O_{2}]_{(GOBAI - ' base_grid ')}'];
    % mlabel off;
    % plabel off;
    % if ~isfolder([param_props.dir_name '/Figures/' base_grid '/' rlz '_gr'])
    %     mkdir([param_props.dir_name '/Figures/' base_grid '/' rlz '_gr']);
    % end
    % date = datevec(time(m));
    % export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
    %     '/' num2str(date(1)) '_' num2str(date(2)) '.png'],'-transparent');
    % close

    %% save monthly differences
    if m == 1
        if ~isfolder(delta_filepath); mkdir(delta_filepath); end
        % create combined files
        ncsave_4d([delta_filepath '/delta_gobai-' param_props.file_name '.nc'],...
            {'lon' lon 'longitude' 'degrees east'},...
            {'lat' lat 'latitude' 'degrees north'},...
            {'depth' depth 'depth' 'meters'},...
            {'time' time(m) 'time' 'time'},...
            {['delta_' param_props.file_name] delta ' ' 'umol/kg'});
    else
        % append to combined files
        ncwrite([delta_filepath '/delta_gobai-' param_props.file_name '.nc'],'time',time(2),m);
        ncwrite([delta_filepath '/delta_gobai-' param_props.file_name '.nc'],['delta_' param_props.file_name],delta,[1 1 1 m]);
    end

end

%% save global means
if ~isfolder([param_props.dir_name '/Data/' base_grid]); mkdir([param_props.dir_name '/Data/' base_grid '/' rlz '_gr']); end
save([param_props.dir_name '/Data/' base_grid '/' rlz '_gr/statistics.mat'],...
    'gobai_mean','gobai_depth_mean','cmip_mean','cmip_depth_mean');

%% plot timeseries
figure;
plot(datenum(1950,1,1)+time,cmip_mean,time,gobai_mean,'LineWidth',2);
legend({base_grid ['GOBAI-' param_props.dir_name '_{(' base_grid ')}']});
datetick('x','keeplimits');
if ~isfolder([param_props.dir_name '/Figures/' base_grid '/' rlz '_gr'])
    mkdir([param_props.dir_name '/Figures/' base_grid '/' rlz '_gr']);
end
export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
        '/timeseries.png'],'-transparent'); close;

%% plot average profile
figure;
plot(mean(cmip_depth_mean,'omitnan'),depth,...
    mean(gobai_depth_mean,'omitnan'),depth,'linewidth',2);
set(gca,'YDir','reverse','XAxisLocation','top');
ylim([0 max(depth)]);
xlabel('[O_{2}] (\mumol kg^{-1})');
ylabel('Depth (m)');
legend({base_grid ['GOBAI-' param_props.dir_name '_{(' base_grid ')}']},'Location','southeast');
export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
        '/profile.png'],'-transparent'); close;

%% plot average profile difference
figure; hold on;
plot(mean(cmip_depth_mean,'omitnan')-...
    mean(gobai_depth_mean,'omitnan'),depth,'linewidth',2);
plot([0 0],[min(depth) max(depth)],'k--');
set(gca,'YDir','reverse','XAxisLocation','top');
ylim([0 max(depth)]);
xlim([-3 3]);
xlabel('\Delta[O_{2}] (\mumol kg^{-1})');
ylabel('Depth (m)');
export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
        '/profile_delta.png'],'-transparent'); close;

end
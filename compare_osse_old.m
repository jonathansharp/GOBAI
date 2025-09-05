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
    val_ratio,test_ratio,numtrees,minLeafSize,...
    numstumps,numbins,rlz)

%% process date
date_str = num2str(snap_date);

%% create directory names (multiple models)
ffnn_filepath = ... % FFNN
    [fpath 'GOBAI/' base_grid '/FFNN/c' ...
    num2str(num_clusters) '_' file_date float_file_ext '/train' ...
    num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
    num2str(100*test_ratio) '/gobai-' param_props.file_name '.nc'];
% rfr_filepath = ... % RFR
%     [fpath 'GOBAI/' base_grid '/RFR/c' ...
%     num2str(num_clusters) '_' file_date float_file_ext '/tr' ...
%     num2str(numtrees) '_lf' num2str(minLeafSize) '/gobai-' ...
%     param_props.file_name '.nc'];
% gbm_filepath = ... % GBM
%     [fpath 'GOBAI/' base_grid '/GBM/c' ...
%     num2str(num_clusters) '_' file_date float_file_ext '/tr' ...
%     num2str(numstumps) '_bin' num2str(numbins) '/gobai-' ...
%     param_props.file_name '.nc'];

%% create directory names (multiple cluster configurations)
% ffnn_filepath_1 = ... % FFNN
%     [fpath 'GOBAI/' base_grid '/FFNN/c' ...
%     num2str(num_clusters_1) '_' file_date float_file_ext '/train' ...
%     num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
%     num2str(100*test_ratio) '/gobai-' param_props.file_name '.nc'];
% ffnn_filepath_2 = ... % FFNN
%     [fpath 'GOBAI/' base_grid '/FFNN/c' ...
%     num2str(num_clusters_2) '_' file_date float_file_ext '/train' ...
%     num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
%     num2str(100*test_ratio) '/gobai-' param_props.file_name '.nc'];
% ffnn_filepath_3 = ... % FFNN
%     [fpath 'GOBAI/' base_grid '/FFNN/c' ...
%     num2str(num_clusters_3) '_' file_date float_file_ext '/train' ...
%     num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
%     num2str(100*test_ratio) '/gobai-' param_props.file_name '.nc'];

%% combined filepath
% gobai
gobai_filepath = ...
    [fpath 'GOBAI/' base_grid '/AVG/c' ...
    num2str(num_clusters) '_' file_date float_file_ext '/gobai-' ...
    param_props.file_name '.nc'];
% gobai_filepath = ...
%     [fpath 'GOBAI/' base_grid '/AVG/multiple_clusters' ...
%     '_' file_date float_file_ext '/gobai-' ...
%     param_props.file_name '.nc'];
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
% rfr_mean = nan(length(time),1);
ffnn_mean = nan(length(time),1);
% ffnn_mean_1 = nan(length(time),1);
% ffnn_mean_2 = nan(length(time),1);
% ffnn_mean_3 = nan(length(time),1);
% gbm_mean = nan(length(time),1);
gobai_mean = nan(length(time),1);
gobai_depth_mean = nan(length(time),length(depth));
cmip_mean = nan(length(time),1);
cmip_depth_mean = nan(length(time),length(depth));

%% loop through each month
for m = 1:240

    %% load monthly output
    % rfr = ncread(rfr_filepath,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);
    ffnn = ncread(ffnn_filepath,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);
    % ffnn1 = ncread(ffnn_filepath_1,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);
    % ffnn2 = ncread(ffnn_filepath_2,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);
    % ffnn3 = ncread(ffnn_filepath_3,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);
    % gbm = ncread(gbm_filepath,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);    
    % gobai = ncread(gobai_filepath,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);
    cmip = ncread(cmip_filepath,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);

    %% apply RFROM mask
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
    % rfr(mask_cmip) = NaN;
    ffnn(mask_cmip) = NaN;
    % ffnn1(mask_cmip) = NaN;
    % ffnn2(mask_cmip) = NaN;
    % ffnn3(mask_cmip) = NaN;
    % gbm(mask_cmip) = NaN;
    % gobai(mask_cmip) = NaN;
    cmip(mask_cmip) = NaN;

    %% determine differences between masked grids
    delta = gobai-cmip;
    % delta_rfr = rfr-cmip;
    % delta_ffnn = ffnn-cmip;
    delta_ffnn_1 = ffnn1-cmip;
    delta_ffnn_2 = ffnn2-cmip;
    delta_ffnn_3 = ffnn3-cmip;
    % delta_gbm = gbm-cmip;

    %% calculate global means
    vol = single(calculate_volume(lat,lon,depth));
    idx_gobai = ~isnan(gobai);
    % Eliminate bottom depth for MPI-ESM1-2-LR
    if strcmp(base_grid,'MPI-ESM1-2-LR'); idx_gobai(:,:,end) = false; end
    % rfr_mean(m) = sum(rfr(idx_gobai).*vol(idx_gobai))./sum(vol(idx_gobai));
    % ffnn_mean(m) = sum(ffnn(idx_gobai).*vol(idx_gobai))./sum(vol(idx_gobai));
    ffnn_mean_1(m) = sum(ffnn1(idx_gobai).*vol(idx_gobai))./sum(vol(idx_gobai));
    ffnn_mean_2(m) = sum(ffnn2(idx_gobai).*vol(idx_gobai))./sum(vol(idx_gobai));
    ffnn_mean_3(m) = sum(ffnn3(idx_gobai).*vol(idx_gobai))./sum(vol(idx_gobai));
    % gbm_mean(m) = sum(gbm(idx_gobai).*vol(idx_gobai))./sum(vol(idx_gobai));
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
    % figure('visible','off');
    % worldmap('world');
    % set(gca,'fontsize',12);
    % pcolorm(lat,lon,delta(:,:,1)');
    % title(datestr(time(m),'mmm yyyy'));
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
        % ncsave_4d([delta_filepath '/delta_gobai-rfr-' param_props.file_name '.nc'],...
        %     {'lon' lon 'longitude' 'degrees east'},...
        %     {'lat' lat 'latitude' 'degrees north'},...
        %     {'depth' depth 'depth' 'meters'},...
        %     {'time' time(m) 'time' 'time'},...
        %     {['delta_' param_props.file_name] delta_rfr ' ' 'umol/kg'});
        % ncsave_4d([delta_filepath '/delta_gobai-ffnn-' param_props.file_name '.nc'],...
        %     {'lon' lon 'longitude' 'degrees east'},...
        %     {'lat' lat 'latitude' 'degrees north'},...
        %     {'depth' depth 'depth' 'meters'},...
        %     {'time' time(m) 'time' 'time'},...
        %     {['delta_' param_props.file_name] delta_ffnn ' ' 'umol/kg'});
        ncsave_4d([delta_filepath '/delta_gobai-ffnn-1-' param_props.file_name '.nc'],...
            {'lon' lon 'longitude' 'degrees east'},...
            {'lat' lat 'latitude' 'degrees north'},...
            {'depth' depth 'depth' 'meters'},...
            {'time' time(m) 'time' 'time'},...
            {['delta_' param_props.file_name] delta_ffnn_1 ' ' 'umol/kg'});
        ncsave_4d([delta_filepath '/delta_gobai-ffnn-2-' param_props.file_name '.nc'],...
            {'lon' lon 'longitude' 'degrees east'},...
            {'lat' lat 'latitude' 'degrees north'},...
            {'depth' depth 'depth' 'meters'},...
            {'time' time(m) 'time' 'time'},...
            {['delta_' param_props.file_name] delta_ffnn_2 ' ' 'umol/kg'});
        ncsave_4d([delta_filepath '/delta_gobai-ffnn-3-' param_props.file_name '.nc'],...
            {'lon' lon 'longitude' 'degrees east'},...
            {'lat' lat 'latitude' 'degrees north'},...
            {'depth' depth 'depth' 'meters'},...
            {'time' time(m) 'time' 'time'},...
            {['delta_' param_props.file_name] delta_ffnn_3 ' ' 'umol/kg'});
        % ncsave_4d([delta_filepath '/delta_gobai-gbm-' param_props.file_name '.nc'],...
        %     {'lon' lon 'longitude' 'degrees east'},...
        %     {'lat' lat 'latitude' 'degrees north'},...
        %     {'depth' depth 'depth' 'meters'},...
        %     {'time' time(m) 'time' 'time'},...
        %     {['delta_' param_props.file_name] delta_gbm ' ' 'umol/kg'});
    else
        % append to combined files
        ncwrite([delta_filepath '/delta_gobai-' param_props.file_name '.nc'],'time',time(2),m);
        ncwrite([delta_filepath '/delta_gobai-' param_props.file_name '.nc'],['delta_' param_props.file_name],delta,[1 1 1 m]);
        % ncwrite([delta_filepath '/delta_gobai-rfr-' param_props.file_name '.nc'],'time',time(2),m);
        % ncwrite([delta_filepath '/delta_gobai-rfr-' param_props.file_name '.nc'],['delta_' param_props.file_name],delta_rfr,[1 1 1 m]);
        % ncwrite([delta_filepath '/delta_gobai-ffnn-' param_props.file_name '.nc'],'time',time(2),m);
        % ncwrite([delta_filepath '/delta_gobai-ffnn-' param_props.file_name '.nc'],['delta_' param_props.file_name],delta_ffnn,[1 1 1 m]);
        ncwrite([delta_filepath '/delta_gobai-ffnn-1-' param_props.file_name '.nc'],'time',time(2),m);
        ncwrite([delta_filepath '/delta_gobai-ffnn-1-' param_props.file_name '.nc'],['delta_' param_props.file_name],delta_ffnn_1,[1 1 1 m]);
        ncwrite([delta_filepath '/delta_gobai-ffnn-2-' param_props.file_name '.nc'],'time',time(2),m);
        ncwrite([delta_filepath '/delta_gobai-ffnn-2-' param_props.file_name '.nc'],['delta_' param_props.file_name],delta_ffnn_2,[1 1 1 m]);
        ncwrite([delta_filepath '/delta_gobai-ffnn-3-' param_props.file_name '.nc'],'time',time(2),m);
        ncwrite([delta_filepath '/delta_gobai-ffnn-3-' param_props.file_name '.nc'],['delta_' param_props.file_name],delta_ffnn_3,[1 1 1 m]);
        % ncwrite([delta_filepath '/delta_gobai-gbm-' param_props.file_name '.nc'],'time',time(2),m);
        % ncwrite([delta_filepath '/delta_gobai-gbm-' param_props.file_name '.nc'],['delta_' param_props.file_name],delta_gbm,[1 1 1 m]);
    end

end

%% save global means
if ~isfolder([param_props.dir_name '/Data/' base_grid]); mkdir([param_props.dir_name '/Data/' base_grid '/' rlz '_gr']); end
% save([param_props.dir_name '/Data/' base_grid '/' rlz '_gr/statistics.mat'],'rfr_mean',...
%     'ffnn_mean','gbm_mean','gobai_mean','gobai_depth_mean','cmip_mean','cmip_depth_mean');
save([param_props.dir_name '/Data/' base_grid '/' rlz '_gr/statistics.mat'],'ffnn_mean_1',...
    'ffnn_mean_2','ffnn_mean_3','gobai_mean','gobai_depth_mean','cmip_mean','cmip_depth_mean');

%% plot timeseries
figure;
plot(time,cmip_mean,time,gobai_mean,'LineWidth',2);
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
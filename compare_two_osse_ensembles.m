%% compare OSSE results
function compare_two_osse_ensembles(param_props,fpath,model_types,realizations,...
    num_clusters,file_date,float_file_ext,train_ratio,val_ratio,test_ratio,...
    start_year,snap_date,ext1,ext2)

%% loop through each model
for m = 1:length(model_types)

    %% plot summary figures

    % define filepaths
    date_str = num2str(snap_date);
    gobai_filepath_1 = [fpath.param_path 'GOBAI/' model_types{m} '/FFNN/c' ...
        num2str(num_clusters) '_' file_date float_file_ext '/train' ...
        num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
        num2str(100*test_ratio) '/' ext1 '/gobai-' param_props.file_name '.nc'];
    delta_filepath_1 = [fpath.param_path 'GOBAI/' model_types{m} '/DELTA/c' ...
        num2str(num_clusters) '_' file_date float_file_ext '/' ext1];
    gobai_filepath_2 = [fpath.param_path 'GOBAI/' model_types{m} '/FFNN/c' ...
        num2str(num_clusters) '_' file_date float_file_ext '/train' ...
        num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
        num2str(100*test_ratio) '/' ext2 '/gobai-' param_props.file_name '.nc'];
    delta_filepath_2 = [fpath.param_path 'GOBAI/' model_types{m} '/DELTA/c' ...
        num2str(num_clusters) '_' file_date float_file_ext '/' ext2];
    path2 = ['_Omon_' model_types{m} '_']; path3 = ['_' realizations{m} '_gr'];
    cmip_filepath = [fpath.model_path model_types{m} '/combined/regridded/' ...
        param_props.file_name path2 'combined' path3 '_' num2str(start_year) ...
        '01-' date_str '.nc'];

    % load time and depth
    time = ncread(gobai_filepath_1,'time')+datenum(1950,0,0);
    depth = ncread(gobai_filepath_1,'depth');

    % load global means
    stats1 = load([param_props.dir_name '/Data/' model_types{m} '/' realizations{m} ...
        '_gr/statistics_' ext1 '.mat']);
    stats2 = load([param_props.dir_name '/Data/' model_types{m} '/' realizations{m} ...
        '_gr/statistics_' ext2 '.mat']);

    % load dimensions
    lon = ncread([delta_filepath_1 '/delta_gobai-' param_props.file_name '.nc'],'lon');
    lat = ncread([delta_filepath_1 '/delta_gobai-' param_props.file_name '.nc'],'lat');
    depth = ncread([delta_filepath_1 '/delta_gobai-' param_props.file_name '.nc'],'depth');
    
    %% delta values
    % load delta values
    delta1 = ncread([delta_filepath_1 '/delta_gobai-' ...
        param_props.file_name '.nc'],['delta_' param_props.file_name]);
    delta2 = ncread([delta_filepath_2 '/delta_gobai-' ...
        param_props.file_name '.nc'],['delta_' param_props.file_name]);
    % calculate spatial means of delta values
    delta_mean_1 = double(mean(delta1,4,'omitnan'));
    delta_mean_2 = double(mean(delta2,4,'omitnan'));
    % calculate temporal rmsds of delta values
    delta_ts_1 = reshape(delta1,[length(lon)*length(lat)*length(depth) length(time)]);
    rmsd_ts_1(:,m) = sqrt(mean(delta_ts_1.^2,'omitnan'))';
    delta_ts_2 = reshape(delta2,[length(lon)*length(lat)*length(depth) length(time)]);
    rmsd_ts_2(:,m) = sqrt(mean(delta_ts_2.^2,'omitnan'))';
    % calculate spatial rmsd of delta value
    rmsd1 = double(sqrt(mean(delta1.^2,4,'omitnan')));
    rmsd2 = double(sqrt(mean(delta2.^2,4,'omitnan')));
    % calculate depth-weighted mean
    vol = weights3d(lon,lat,depth); vol(isnan(delta_mean_1)) = NaN;
    delta_wtd_mean_1(:,:,m) = sum(delta_mean_1.*vol,3,'omitnan')./sum(vol,3,'omitnan');
    vol = weights3d(lon,lat,depth); vol(isnan(delta_mean_2)) = NaN;
    delta_wtd_mean_2(:,:,m) = sum(delta_mean_2.*vol,3,'omitnan')./sum(vol,3,'omitnan');

    %% parameter values
    % load parameter from gobai
    var1 = ncread(gobai_filepath_1,param_props.file_name);
    var2 = ncread(gobai_filepath_2,param_props.file_name);
    N_idx = lat > 0; S_idx = lat <= 0;
    depth_idx = depth < 300; % index to above 300m
    var_ts_N_1 = reshape(var1(:,N_idx,depth_idx,:),[length(lon)*length(lat(N_idx))*length(depth(depth_idx)) length(time)]);
    var_mean_ts_N_1 = mean(var_ts_N_1,'omitnan')'-mean(mean(var_ts_N_1,'omitnan')');
    var_ts_S_1 = reshape(var1(:,S_idx,depth_idx,:),[length(lon)*length(lat(S_idx))*length(depth(depth_idx)) length(time)]);
    var_mean_ts_S_1 = mean(var_ts_S_1,'omitnan')'-mean(mean(var_ts_S_1,'omitnan')');
    var_ts_N_2 = reshape(var2(:,N_idx,depth_idx,:),[length(lon)*length(lat(N_idx))*length(depth(depth_idx)) length(time)]);
    var_mean_ts_N_2 = mean(var_ts_N_2,'omitnan')'-mean(mean(var_ts_N_2,'omitnan')');
    var_ts_S_2 = reshape(var2(:,S_idx,depth_idx,:),[length(lon)*length(lat(S_idx))*length(depth(depth_idx)) length(time)]);
    var_mean_ts_S_2 = mean(var_ts_S_2,'omitnan')'-mean(mean(var_ts_S_2,'omitnan')');

    % cmip
    cmip = ncread(cmip_filepath,param_props.file_name);
    cmip_ts_N = reshape(cmip(:,N_idx,depth_idx,:),[length(lon)*length(lat(N_idx))*length(depth(depth_idx)) length(time)]);
    cmip_mean_ts_N = mean(cmip_ts_N,'omitnan')'-mean(mean(cmip_ts_N,'omitnan')');
    cmip_ts_S = reshape(cmip(:,S_idx,depth_idx,:),[length(lon)*length(lat(S_idx))*length(depth(depth_idx)) length(time)]);
    cmip_mean_ts_S = mean(cmip_ts_S,'omitnan')'-mean(mean(cmip_ts_S,'omitnan')');
    
    % calculate climatology
    for mnth = 1:12
        var_clim_anom_ts_N_1(mnth,m) = mean(var_mean_ts_N_1(mnth:12:end));
        var_clim_anom_ts_S_1(mnth,m) = mean(var_mean_ts_S_1(mnth:12:end));
        var_clim_anom_ts_N_2(mnth,m) = mean(var_mean_ts_N_2(mnth:12:end));
        var_clim_anom_ts_S_2(mnth,m) = mean(var_mean_ts_S_2(mnth:12:end));
        cmip_clim_anom_ts_N(mnth,m) = mean(cmip_mean_ts_N(mnth:12:end));
        cmip_clim_anom_ts_S(mnth,m) = mean(cmip_mean_ts_S(mnth:12:end));
    end

    % var_mean = double(mean(var,4,'omitnan'));
    % vol = weights3d(lon,lat,depth);
    % vol(isnan(delta_mean)) = NaN;
    % idx_depth = find(abs(depth-500)==min(abs(depth-500)));
    % var_wtd_mean(:,:,m) = sum(var_mean(:,:,1:idx_depth).*vol(:,:,1:idx_depth),3,'omitnan')./...
    %     sum(vol(:,:,1:idx_depth),3,'omitnan');


end

% plot ensemble mean differences
z1 = [mean(delta_wtd_mean_1,3,'omitnan');mean(delta_wtd_mean_1(end,:,:),3,'omitnan')]';
z2 = [mean(delta_wtd_mean_2,3,'omitnan');mean(delta_wtd_mean_2(end,:,:),3,'omitnan')]';
figure('visible','on');
worldmap([-90 90],[20 380]);
set(gca,'fontsize',12);
pcolorm(lat,[lon;lon(end)+1],z1-z2);
title('Ensemble Average');
plot_land('map');
c=colorbar;
caxis([-15 15]);
colormap(cmocean('balance'));
c.Label.String = ['Avg. \Delta[O_{2}]_{(GOBAI - ESM)}'];
mlabel off;
plabel off;
export_fig(gcf,[param_props.dir_name '/Figures/ensemble_mean_delta_' ...
   ext1 '_versus_' ext2 '.png'],'-transparent');

% plot boxplot
model_categories = categorical(repmat(model_types,length(time)*2,1));
sub_categories = repelem({'GLODAP';'GLODAP+Float'},length(time),5);
rmsd_ts = [rmsd_ts_1;rmsd_ts_2];
% boxplot
figure('visible','on'); hold on;
title('Monthly Mean \Delta[O_{2}] RMSD (\mumol kg^{-1})');
boxchart(model_categories(:),rmsd_ts(:),'GroupByColor',sub_categories(:),...
    'LineWidth',2,'MarkerStyle','none');
legend({'GLODAP' 'GLODAP+Float'},'Location','northwest');
export_fig(gcf,[param_props.dir_name '/Figures/monthly_rmsd_' ...
   ext1 '_versus_' ext2 '.png'],'-transparent');
close

% plot climatology
figure('visible','on'); hold on;
title('0-300m Climatological [O_{2}] Anomaly (\mumol kg^{-1})');
clrs = colororder;
% for m = 1:length(model_types)
%     plot(1:12,var_clim_anom_ts_N_1(:,m),'Color',clrs(1,:),'LineWidth',1);
%     plot(1:12,var_clim_anom_ts_N_2(:,m),'Color',clrs(2,:),'LineWidth',1);
%     plot(1:12,cmip_clim_anom_ts_N(:,m),'Color',clrs(3,:),'LineWidth',1);
% end
plot(1:12,mean(var_clim_anom_ts_N_1,2),'Color',clrs(1,:),'LineWidth',3);
plot(1:12,mean(var_clim_anom_ts_N_2,2),'Color',clrs(2,:),'LineWidth',3);
plot(1:12,mean(cmip_clim_anom_ts_N,2),'Color',clrs(3,:),'LineWidth',3);
legend({'GLODAP Sampling' 'GLODAP+Float Sampling' 'Model "Truth"'},...
    'Location','southwest');
export_fig(gcf,[param_props.dir_name '/Figures/climatological_anomaly_' ...
   ext1 '_versus_' ext2 '.png'],'-transparent');
close
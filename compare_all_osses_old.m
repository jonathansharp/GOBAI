% compare_all_osses
%
% DESCRIPTION:
% This function creates plots using the results of all osses 
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/7/2024

function compare_all_osses(param_props,fpath,model_types,realizations,...
    num_clusters,file_date,float_file_ext,train_ratio,val_ratio,test_ratio)

% establish colors
clrs = colororder;

% establish timeseries of residuals figure
figure(1); fig = gcf; hold on;
fig.Position(3) = fig.Position(3).*2;
fig.Position(4) = fig.Position(4).*1.3;
set(gca,'fontsize',18);
xlim([datenum(2003,1,1) datenum(2025,1,1)]);
datetick('x');
ylim([-2 2]);
plot([datenum(2003,1,1) datenum(2025,1,1)],[0 0],'k--')
xlabel('Year');
ylabel('\Delta[O_{2}] (GOBAI - CMIP)');
hold off;

% establish detailed timeseries of residuals figure
figure(2); fig = gcf; hold on;
fig.Position(3) = fig.Position(3).*2;
fig.Position(4) = fig.Position(4).*1.3;
set(gca,'fontsize',18);
xlim([datenum(2003,1,1) datenum(2025,1,1)]);
datetick('x');
ylim([-2 2]);
plot([datenum(2003,1,1) datenum(2025,1,1)],[0 0],'k--')
xlabel('Year');
ylabel('\Delta[O_{2}] (GOBAI - CMIP)');
hold off;

% establish timeseries of global means figure
figure(3); fig = gcf; hold on;
fig.Position(3) = fig.Position(3).*2;
fig.Position(4) = fig.Position(4).*1.3;
set(gca,'fontsize',18);
xlim([datenum(2003,1,1) datenum(2025,1,1)]);
datetick('x');
plot([datenum(2003,1,1) datenum(2025,1,1)],[0 0],'k--')
xlabel('Year');
ylabel('[O_{2}] (\mumol kg^{-1})');
hold off;

% establish profile figure
figure(4); fig = gcf; hold on;
fig.Position(3) = fig.Position(3).*0.6;
fig.Position(4) = fig.Position(4).*2;
set(gca,'fontsize',16);
set(gca,'YDir','reverse');
xlim([-4 4]);
ylim([0 2]);
plot([0 0],[0 2000],'k--');
xlabel('\Delta[O_{2}] (GOBAI - CMIP)');
ylabel('Depth (km)');
hold off;

%% loop through each model
for m = 1:length(model_types)

    %% plot summary figures

    % define filepaths
    gobai_filepath = [fpath model_types{m} '/GOBAI/' model_types{m} '/FFNN/c' ...
        num2str(num_clusters) '_' file_date float_file_ext '/train' ...
        num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
        num2str(100*test_ratio) '/gobai-' param_props.file_name '.nc'];
    delta_filepath = [fpath model_types{m} '/GOBAI/' model_types{m} '/DELTA/c' ...
        num2str(num_clusters) '_' file_date float_file_ext];
    % gobai_filepath = [fpath model_types{m} '/GOBAI/' model_types{m} '/AVG/' ...
    %     'multiple_clusters_' file_date float_file_ext '/gobai-' ...
    %     param_props.file_name '.nc'];
    % delta_filepath = [fpath model_types{m} '/GOBAI/' model_types{m} '/DELTA/' ...
    %     'multiple_clusters_' file_date float_file_ext];

    % load time and depth
    time = ncread(gobai_filepath,'time')+datenum(1950,0,0);
    depth = ncread(gobai_filepath,'depth');

    % load global means
    load([param_props.dir_name '/Data/' model_types{m} '/' realizations{m} ...
        '_gr/statistics.mat']);

    % plot timeseries of residuals
    figure(1); hold on;
    idx = ~isnan(gobai_mean) & ~isnan(cmip_mean);
    % % FFNN, RFR, and GBM
    % var1 = (gobai_mean-cmip_mean)+std([rfr_mean-cmip_mean,ffnn_mean-cmip_mean,gbm_mean-cmip_mean],[],2,'omitnan');
    % var2 = (gobai_mean-cmip_mean)-std([rfr_mean-cmip_mean,ffnn_mean-cmip_mean,gbm_mean-cmip_mean],[],2,'omitnan');
    % % FFNNs with different configurations
    % var1 = (gobai_mean-cmip_mean)+std([ffnn_mean_1-cmip_mean,ffnn_mean_2-cmip_mean,ffnn_mean_3-cmip_mean],[],2,'omitnan');
    % var2 = (gobai_mean-cmip_mean)-std([ffnn_mean_1-cmip_mean,ffnn_mean_2-cmip_mean,ffnn_mean_3-cmip_mean],[],2,'omitnan');
    % just FFNN
    % var1 = (gobai_mean-cmip_mean)+std([ffnn_mean-cmip_mean],[],2,'omitnan');
    % var2 = (gobai_mean-cmip_mean)-std([ffnn_mean-cmip_mean],[],2,'omitnan');
    % pp=patch([time(idx);flipud(time(idx))],[var1(idx);flipud(var2(idx))],...
    %     clrs(m,:),'facealpha',0.3,'linestyle','none');
    % uistack(pp,'bottom');
    ts_plot(m) = plot(time,gobai_mean-cmip_mean,'color',clrs(m,:),'linewidth',3);
    hold off;

    % plot detailed timeseries of residuals
    figure(2); hold on;
    plot(time,gobai_mean-cmip_mean,'color',clrs(m,:),'linewidth',3);
    % p1=plot(time,rfr_mean-cmip_mean,'--','color',clrs(m,:),'linewidth',1);
    p2=plot(time,ffnn_mean-cmip_mean,':','color',clrs(m,:),'linewidth',1);
    % p3=plot(time,gbm_mean-cmip_mean,'-.','color',clrs(m,:),'linewidth',1);
    % p1=plot(time,ffnn_mean_1-cmip_mean,'--','color',clrs(m,:),'linewidth',1);
    % p2=plot(time,ffnn_mean_2-cmip_mean,':','color',clrs(m,:),'linewidth',1);
    % p3=plot(time,ffnn_mean_3-cmip_mean,'-.','color',clrs(m,:),'linewidth',1);
    % legend([p1 p2 p3],{'RFR' 'FFNN' 'GBM'});
    legend([p1 p2 p3],{'FFNN-15' 'FFNN-20' 'FFNN-25'});
    export_fig(gcf,[param_props.dir_name '/Figures/osse_detailed_timeseries_' ...
        model_types{m} '.png'],'-transparent');
    figure(2); cla;
    hold off;

    % plot timeseries of each global mean
    figure(3); hold on;
    ylim([floor(min([cmip_mean;gobai_mean])) ceil(max([cmip_mean;gobai_mean]))]);
    p1=plot(time,cmip_mean,':','color',clrs(m,:),'linewidth',3);
    p2=plot(time,gobai_mean,'color',clrs(m,:),'linewidth',3);
    % p3=plot(time,rfr_mean,'--','color',clrs(m,:),'linewidth',1);
    p3=plot(time,ffnn_mean,':','color',clrs(m,:),'linewidth',1);
    % p5=plot(time,gbm_mean,'-.','color',clrs(m,:),'linewidth',1);
    % p3=plot(time,ffnn_mean_1,'--','color',clrs(m,:),'linewidth',1);
    % p4=plot(time,ffnn_mean_2,':','color',clrs(m,:),'linewidth',1);
    % p5=plot(time,ffnn_mean_3,'-.','color',clrs(m,:),'linewidth',1);
    legend([p1 p2 p3],{model_types{m} ['GOBAI-O_{2(' model_types{m} ')}'] ...
        ['GOBAI-FFNN-O_{2(' model_types{m} ')}']});
    export_fig(gcf,[param_props.dir_name '/Figures/osse_global_mean_timeseries_' ...
        model_types{m} '.png'],'-transparent');
    figure(3); cla;
    hold off;

    % plot depth profile
    figure(4); hold on;
    diff = gobai_depth_mean-cmip_depth_mean;
    var1 = mean(diff,1,'omitnan')+std(diff,[],1,'omitnan');
    var2 = mean(diff,1,'omitnan')-std(diff,[],1,'omitnan');
    idx = ~isnan(var1);
    pp=patch([var1(idx),fliplr(var2(idx))],[depth(idx);flipud(depth(idx))]./1e3,...
        clrs(m,:),'facealpha',0.3,'linestyle','none');
    uistack(pp,'bottom');
    prof_plot(m) = plot(mean(diff,1,'omitnan'),depth./1e3,'color',clrs(m,:),'linewidth',3);
    hold off;

    %% plot mapped differences for each model

    % load dimensions
    lon = ncread([delta_filepath '/delta_gobai-' param_props.file_name '.nc'],'lon');
    lat = ncread([delta_filepath '/delta_gobai-' param_props.file_name '.nc'],'lat');
    depth = ncread([delta_filepath '/delta_gobai-' param_props.file_name '.nc'],'depth');
    
    % load delta values
    delta = ncread([delta_filepath '/delta_gobai-' param_props.file_name '.nc'],['delta_' param_props.file_name]);
    delta_mean = double(mean(delta,4,'omitnan'));
    % delta_rfr = ncread([delta_filepath '/delta_gobai-rfr-' param_props.file_name '.nc'],['delta_' param_props.file_name]);
    % delta_mean_rfr = double(mean(delta_rfr,4,'omitnan'));
    % delta_ffnn = ncread([delta_filepath '/delta_gobai-ffnn-' param_props.file_name '.nc'],['delta_' param_props.file_name]);
    % delta_mean_ffnn = double(mean(delta_ffnn,4,'omitnan'));
    % delta_gbm = ncread([delta_filepath '/delta_gobai-gbm-' param_props.file_name '.nc'],['delta_' param_props.file_name]);
    % delta_mean_gbm = double(mean(delta_gbm,4,'omitnan'));
    vol = weights3d(lon,lat,depth);
    vol(isnan(delta_mean)) = NaN;
    delta_wtd_mean(:,:,m) = sum(delta_mean.*vol,3,'omitnan')./sum(vol,3,'omitnan');
    % delta_wtd_mean_rfr(:,:,m) = sum(delta_mean_rfr.*vol,3,'omitnan')./sum(vol,3,'omitnan');
    % delta_wtd_mean_ffnn(:,:,m) = sum(delta_mean_ffnn.*vol,3,'omitnan')./sum(vol,3,'omitnan');
    % delta_wtd_mean_gbm(:,:,m) = sum(delta_mean_gbm.*vol,3,'omitnan')./sum(vol,3,'omitnan');

    % load parameter from gobai
    var = ncread(gobai_filepath,param_props.file_name);
    var_mean = double(mean(var,4,'omitnan'));
    vol = weights3d(lon,lat,depth);
    vol(isnan(delta_mean)) = NaN;
    idx_depth = find(abs(depth-500)==min(abs(depth-500)));
    var_wtd_mean(:,:,m) = sum(var_mean(:,:,1:idx_depth).*vol(:,:,1:idx_depth),3,'omitnan')./...
        sum(vol(:,:,1:idx_depth),3,'omitnan');
    
    % plot AVG differences
    figure; hold on;
    worldmap([-90 90],[20 380]);
    set(gca,'fontsize',12);
    pcolorm(lat,[lon;lon(end)+1],[delta_wtd_mean(:,:,m);delta_wtd_mean(end,:,m)]');
    title(model_types{m});
    plot_land('map');
    c=colorbar;
    caxis([-25 25]);
    colormap(cmocean('balance'));
    c.Label.String = ['Avg. \Delta[O_{2}]_{(GOBAI - ' model_types{m} ')}'];
    mlabel off;
    plabel off;
    if ~isfolder([param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr'])
        mkdir([param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr']);
    end
    export_fig(gcf,[param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr' ...
        '/delta.png'],'-transparent');
    close;

    % % plot RFR differences
    % figure; hold on;
    % worldmap([-90 90],[20 380]);
    % set(gca,'fontsize',12);
    % pcolorm(lat,[lon;lon(end)+1],[delta_wtd_mean_rfr(:,:,m);delta_wtd_mean_rfr(end,:,m)]');
    % title(model_types{m});
    % plot_land('map');
    % c=colorbar;
    % caxis([-25 25]);
    % colormap(cmocean('balance'));
    % c.Label.String = ['Avg. \Delta[O_{2}]_{(GOBAI - RFR - ' model_types{m} ')}'];
    % mlabel off;
    % plabel off;
    % if ~isfolder([param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr'])
    %     mkdir([param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr']);
    % end
    % export_fig(gcf,[param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr' ...
    %     '/delta-rfr.png'],'-transparent');
    % close;
    % 
    % % plot FFNN differences
    % figure; hold on;
    % worldmap([-90 90],[20 380]);
    % set(gca,'fontsize',12);
    % pcolorm(lat,[lon;lon(end)+1],[delta_wtd_mean_ffnn(:,:,m);delta_wtd_mean_ffnn(end,:,m)]');
    % title(model_types{m});
    % plot_land('map');
    % c=colorbar;
    % caxis([-25 25]);
    % colormap(cmocean('balance'));
    % c.Label.String = ['Avg. \Delta[O_{2}]_{(GOBAI - FFNN - ' model_types{m} ')}'];
    % mlabel off;
    % plabel off;
    % if ~isfolder([param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr'])
    %     mkdir([param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr']);
    % end
    % export_fig(gcf,[param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr' ...
    %     '/delta-ffnn.png'],'-transparent');
    % close;
    % 
    % % plot GBM differences
    % figure; hold on;
    % worldmap([-90 90],[20 380]);
    % set(gca,'fontsize',12);
    % pcolorm(lat,[lon;lon(end)+1],[delta_wtd_mean_gbm(:,:,m);delta_wtd_mean_gbm(end,:,m)]');
    % title(model_types{m});
    % plot_land('map');
    % c=colorbar;
    % caxis([-25 25]);
    % colormap(cmocean('balance'));
    % c.Label.String = ['Avg. \Delta[O_{2}]_{(GOBAI - GBM - ' model_types{m} ')}'];
    % mlabel off;
    % plabel off;
    % if ~isfolder([param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr'])
    %     mkdir([param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr']);
    % end
    % export_fig(gcf,[param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr' ...
    %     '/delta-gbm.png'],'-transparent');
    % close;

    % plot mean to 500m
    figure(5); hold on;
    worldmap([-90 90],[20 380]);
    set(gca,'fontsize',12);
    pcolorm(lat,[lon;lon(end)+1],[var_wtd_mean(:,:,m);var_wtd_mean(end,:,m)]');
    title(['GOBAI-O_{2(' model_types{m} ')}']);
    plot_land('map');
    c=colorbar;
    %caxis([-25 ]);
    colormap(cmocean('ice'));
    c.Label.String = 'Avg. [O_{2}] (0-500m, \mumol kg^{-1})';
    mlabel off;
    plabel off;
    if ~isfolder([param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr'])
        mkdir([param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr']);
    end
    export_fig(gcf,[param_props.dir_name '/Figures/' model_types{m} '/' realizations{m} '_gr' ...
        '/' param_props.fig_name '.png'],'-transparent');
    close;

    %% display statistics
    % global monthly
    disp([model_types{m} ' Avg. Diff. (Global Monthly Means) = ' ...
        num2str(mean(gobai_mean-cmip_mean,1,'omitnan')) ' +/- ' ...
        num2str(std(gobai_mean-cmip_mean,1,'omitnan')) ' umol/kg']);
    % grid-cell-level
    disp([model_types{m} ' Avg. Diff. (Grid-cell-level Values) = ' ...
        num2str(mean(delta(:),1,'omitnan')) ' +/- ' ...
        num2str(std(delta(:),1,'omitnan')) 'umol/kg']);
    % cmip trend
    idx = ~isnan(cmip_mean);
    [~,yr_cmip,x_cmip,err_cmip] = leastsq2(time(idx),cmip_mean(idx),time(1),2,[365.2424 365.2424/2]);
    [~,acor,lag,dof] = autocov(time(idx),yr_cmip(idx),365.2424*5); % Determined lagged autocorrelation of residuals
    % figure; plot(lag,acor); 
    % figure; plot(time(idx),yf_cmip(idx),time(idx),cmip_mean(idx));
    edof = dof-6; if edof < 1; edof = 1; end
    trend = x_cmip(2)*365.2424*10; % scale to decadal
    uncer = err_cmip(2)*365.2424*10*(sqrt(length(time))/sqrt(edof))*2; % scale to decadal, by eDOF, and to 95%
    disp([model_types{m} ' Trend  = ' num2str(trend) ' +/- ' num2str(uncer) ' umol/kg/year']);
    % gobai trend
    idx = ~isnan(gobai_mean);
    [~,yr_gobai,x_gobai,err_gobai] = leastsq2(time(idx),gobai_mean(idx),time(1),2,[365.2424 365.2424/2]);
    [~,acor,lag,dof] = autocov(time(idx),yr_gobai(idx),365.2424*5); % Determined lagged autocorrelation of residuals
    edof = dof-6; if edof < 1; edof = 1; end
    trend = x_gobai(2)*365.2424*10; % scale to decadal
    uncer = err_gobai(2)*365.2424*10*(sqrt(length(time(idx)))/sqrt(edof))*2; % scale to decadal, by eDOF, and to 95%
    disp([model_types{m} ' Trend  = ' num2str(trend) ' +/- ' num2str(uncer) ' umol/kg/year']);
    % trend diff
    disp([model_types{m} ' Trend Diff. = ' num2str((x_gobai(2)-x_cmip(2))*365.2424*10) ...
        ' umol/kg/year']);

    %% display algorith-specific statistics
    % % global monthly (rfr)
    % disp([model_types{m} ' RFR Avg. Diff. (Global Monthly Means) = ' ...
    %     num2str(mean(rfr_mean-cmip_mean,1,'omitnan')) ' +/- ' ...
    %     num2str(std(rfr_mean-cmip_mean,1,'omitnan')) ' umol/kg']);
    % % global monthly (ffnn)
    % disp([model_types{m} ' FFNN Avg. Diff. (Global Monthly Means) = ' ...
    %     num2str(mean(ffnn_mean-cmip_mean,1,'omitnan')) ' +/- ' ...
    %     num2str(std(ffnn_mean-cmip_mean,1,'omitnan')) ' umol/kg']);
    % % global monthly (gbm)
    % disp([model_types{m} ' GBM Avg. Diff. (Global Monthly Means) = ' ...
    %     num2str(mean(gbm_mean-cmip_mean,1,'omitnan')) ' +/- ' ...
    %     num2str(std(gbm_mean-cmip_mean,1,'omitnan')) ' umol/kg']);
    % global monthly (ffnn-1)
    disp([model_types{m} ' FFNN Avg. Diff. (Global Monthly Means) = ' ...
        num2str(mean(ffnn_mean_1-cmip_mean,1,'omitnan')) ' +/- ' ...
        num2str(std(ffnn_mean_1-cmip_mean,1,'omitnan')) ' umol/kg']);
    % global monthly (ffnn-2)
    disp([model_types{m} ' FFNN Avg. Diff. (Global Monthly Means) = ' ...
        num2str(mean(ffnn_mean_2-cmip_mean,1,'omitnan')) ' +/- ' ...
        num2str(std(ffnn_mean_2-cmip_mean,1,'omitnan')) ' umol/kg']);
    % global monthly (ffnn-3)
    disp([model_types{m} ' FFNN Avg. Diff. (Global Monthly Means) = ' ...
        num2str(mean(ffnn_mean_3-cmip_mean,1,'omitnan')) ' +/- ' ...
        num2str(std(ffnn_mean_3-cmip_mean,1,'omitnan')) ' umol/kg']);

end

figure(1); hold on;
%legend(ts_plot,model_types,'location','northoutside','numcolumns',5,'FontSize',12);
export_fig(gcf,[param_props.dir_name '/Figures/osse_timeseries_residuals.png'],'-transparent'); close;

figure(2); hold on; close;

figure(3); hold on; close;

figure(4); hold on;
legend(prof_plot,model_types,'location','northoutside','FontSize',12);
export_fig(gcf,[param_props.dir_name '/Figures/osse_profile_delta.png'],'-transparent'); close;

% plot ensemble mean differences
figure('visible','on');
worldmap([-90 90],[20 380]);
set(gca,'fontsize',12);
pcolorm(lat,[lon;lon(end)+1],[mean(delta_wtd_mean,3,'omitnan');mean(delta_wtd_mean(end,:,:),3,'omitnan')]');
title('Ensemble Average');
plot_land('map');
c=colorbar;
caxis([-25 25]);
colormap(cmocean('balance'));
c.Label.String = ['Avg. \Delta[O_{2}]_{(GOBAI - ' model_types{m} ')}'];
mlabel off;
plabel off;
export_fig(gcf,[param_props.dir_name '/Figures/ensemble_mean_delta.png'],'-transparent');
close

% plot ensemble variability
figure('visible','on');
worldmap([-90 90],[20 380]);
set(gca,'fontsize',12);
pcolorm(lat,[lon;lon(end)+1],[std(delta_wtd_mean,[],3,'omitnan');std(delta_wtd_mean(end,:,:),[],3,'omitnan')]');
title('Ensemble Variability');
plot_land('map');
c=colorbar;
caxis([0 20]);
colormap(cmocean('tempo'));
c.Label.String = ['\Delta[O_{2}]_{(GOBAI - ESM)} Var.'];
mlabel off;
plabel off;
export_fig(gcf,[param_props.dir_name '/Figures/ensemble_mean_variability.png'],'-transparent');
close

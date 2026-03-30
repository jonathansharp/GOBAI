% compare_osse
%
% DESCRIPTION:
% This function compares reconstructed fields from GOBAI OSSEs
% to full model gridded fields.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/7/2024

function compare_osse(param_props,fpaths,base_grid,file_date,...
    float_file_ext,num_clusters,start_year,snap_date,train_ratio,...
    val_ratio,test_ratio,rlz,float_ext,glodap_ext,ctd_ext,...
    numWorkers_predict,varargin)

%% process date
date_str = num2str(snap_date);
% regions
regions = {'global' 'atlantic' 'pacific' 'indian' 'med' 'southern'};

%% process optional input arguments
% pre-allocate
coverage = '';
% process inputs
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'coverage')
        coverage = ['_' num2str(varargin{i+1}) '_coverage'];
    end
end

%% combined filepath
% osse
osse_filepath = ...
    [fpaths.param_path 'GOBAI/' base_grid '/FFNN/c' ...
        num2str(num_clusters) '_' file_date float_file_ext '/train' ...
        num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
        num2str(100*test_ratio) '/' float_ext glodap_ext ctd_ext coverage ...
        '/gobai-' param_props.file_name '.nc'];
% cmip
path2 = ['_Omon_' base_grid '_'];
path3 = ['_' rlz '_gr'];
cmip_filepath = [fpaths.model_path base_grid '/combined/regridded/' param_props.file_name path2 ...
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
cmip_dens_filepath = [fpaths.model_path base_grid '/combined/regridded/' ...
    'dens' path2 'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
% delta
delta_filepath = [fpaths.param_path 'GOBAI/' base_grid '/DELTA/c' ...
    num2str(num_clusters) '_' file_date float_file_ext '/' ...
    float_ext glodap_ext ctd_ext];

%% load dimensions
lat = ncread(osse_filepath,'lat');
lon = ncread(osse_filepath,'lon');
depth = ncread(osse_filepath,'depth');
dpth_bnds = ncread(osse_filepath,'dpth_bnds');
time = ncread(osse_filepath,'time');

%% load and evaluate monthly reconstructed gobai in relation to cmip
for r = 1:length(regions)
    gobai_mean.(regions{r}) = nan(length(time),1);
    cmip_mean.(regions{r}) = nan(length(time),1);
    gobai_inv.(regions{r}) = nan(length(time),1);
    cmip_inv.(regions{r}) = nan(length(time),1);
    gobai_mean_both.(regions{r}) = nan(length(time),1);
    cmip_mean_both.(regions{r}) = nan(length(time),1);
    gobai_inv_both.(regions{r}) = nan(length(time),1);
    cmip_inv_both.(regions{r}) = nan(length(time),1);
end
gobai_depth_mean = nan(length(time),length(depth));
cmip_depth_mean = nan(length(time),length(depth));

% %% load T/S timeframe
% rg_time = [pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'];

% column mean
gobai_col_mean = nan(length(lon),length(lat),length(time));
cmip_col_mean = nan(length(lon),length(lat),length(time));

% create folder for figures
if ~isfolder([param_props.dir_name '/Figures/' base_grid '/' rlz '_gr'])
    mkdir([param_props.dir_name '/Figures/' base_grid '/' rlz '_gr']);
end

%% loop through each month
for m = 1:length(time)

    %% load monthly output
    gobai = ncread(osse_filepath,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);
    cmip = ncread(cmip_filepath,param_props.file_name,[1 1 1 m],[Inf Inf Inf 1]);
    cmip_dens = ncread(cmip_dens_filepath,'dens',[1 1 1 m],[Inf Inf Inf 1]);

    %% apply RG mask
    % % load RG grid
    % lat_rg = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Latitude');
    % lon_rg = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Longitude');
    % temp_rg = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'],...
    %     'Temperature',[1 1 1 m],[Inf Inf 1 1]);
    % mask_rg = (isnan(temp_rg));
    % % reorder and pad mask
    % mask_cmip = [mask_rg(341:end,:,:);mask_rg(1:340,:,:)];
    % mask_cmip = [true(length(lon),25,length(depth)),...
    %     repmat(mask_cmip,1,1,length(depth)),true(length(lon),10,length(depth))];
    % % apply mask to monthly output
    % gobai(mask_cmip) = NaN;
    % cmip(mask_cmip) = NaN;

    %% determine differences between masked grids
    delta = gobai-cmip;

    %% calculate means and inventories
    vol = single(calculate_volume(lat,lon,depth,dpth_bnds));
    kg = cmip_dens.*vol; % kilograms in each cell
    for r = 1:length(regions)
        if strcmp(regions{r},'global')
            idx_gobai = ~isnan(gobai);
            idx_cmip = ~isnan(cmip);
        elseif strcmp(regions{r},'atlantic')
            reg_idx = ncread('RECCAP2_region_masks_all_v20221025.nc','atlantic');
            idx_gobai = ~isnan(gobai) & reg_idx >= 1 & reg_idx <= 5;
            idx_cmip = ~isnan(cmip) & reg_idx >= 1 & reg_idx <= 5;

        elseif strcmp(regions{r},'med')
            reg_idx = ncread('RECCAP2_region_masks_all_v20221025.nc','atlantic');
            idx_gobai = ~isnan(gobai) & reg_idx == 6;
            idx_cmip = ~isnan(cmip) & reg_idx == 6;
        else
            reg_idx = ncread('RECCAP2_region_masks_all_v20221025.nc',regions{r});
            idx_gobai = ~isnan(gobai) & reg_idx > 0;
            idx_cmip = ~isnan(cmip) & reg_idx > 0;
        end
        idx_both = idx_gobai & idx_cmip;
        % % Eliminate bottom depth for MPI-ESM1-2-LR
        % if strcmp(base_grid,'MPI-ESM1-2-LR'); idx_gobai(:,:,end) = false; end
        gobai_mean.(regions{r})(m) = ... % umol/kg
            sum(gobai(idx_gobai).*vol(idx_gobai))./sum(vol(idx_gobai));
        cmip_mean.(regions{r})(m) = ... % umol/kg
            sum(cmip(idx_cmip).*vol(idx_cmip))./sum(vol(idx_cmip));
        gobai_mean_both.(regions{r})(m) = ... % umol/kg
            sum(gobai(idx_both).*vol(idx_both))./sum(vol(idx_both));
        cmip_mean_both.(regions{r})(m) = ... % umol/kg
            sum(cmip(idx_both).*vol(idx_both))./sum(vol(idx_both));
        
        gobai_inv.(regions{r})(m) = ... % Pmol
            sum(gobai(idx_gobai).*kg(idx_gobai))./(10^21);
        cmip_inv.(regions{r})(m) = ... % Pmol
            sum(cmip(idx_cmip).*kg(idx_cmip))./(10^21);
        gobai_inv_both.(regions{r})(m) = ... % Pmol
            sum(gobai(idx_both).*kg(idx_both))./(10^21);
        cmip_inv_both.(regions{r})(m) = ... % Pmol
            sum(cmip(idx_both).*kg(idx_both))./(10^21);

    end

    %% calculate column means
    % gobai (mol/m2)
    
    % replicate depth bin heights
    dpth_bins_3d = repmat(permute(diff(dpth_bnds),[3 2 1]),length(lon),length(lat),1);
    % Eliminate bottom depth for MPI-ESM1-2-LR
    % if strcmp(base_grid,'MPI-ESM1-2-LR'); idx_gobai(:,:,end) = false; end
    gobai_col_mean_temp = sum(((gobai./10^6).*1026).*dpth_bins_3d,3,'omitnan'); % density estimate for now...
    gobai_col_mean_temp(gobai_col_mean_temp==0) = NaN;
    gobai_col_mean(:,:,m) = gobai_col_mean_temp;
    % cmip (mol/m2)
    cmip_col_mean_temp = sum(((cmip./10^6).*1026).*dpth_bins_3d,3,'omitnan'); % density estimate for now...
    cmip_col_mean_temp(cmip_col_mean_temp==0) = NaN;
    cmip_col_mean(:,:,m) = cmip_col_mean_temp;

    %% calculate global means by depth
    idx_gobai = ~isnan(gobai);
    idx_cmip = ~isnan(cmip);
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

%% calculate column mean trends
gobai_trend = nan(length(lon),length(lat));
cmip_trend = nan(length(lon),length(lat));
for x = 1:length(lon)
    for y = 1:length(lat)
        if any(~isnan(gobai_col_mean(x,y,:)))
            % gobai trend
            gobai_col_mean_temp = squeeze(gobai_col_mean(x,y,:));
            tr = polyfit(time,gobai_col_mean_temp,1);
            gobai_trend(x,y) = tr(1).*365.25.*10;
        end
        if any(~isnan(cmip_col_mean(x,y,:)))
            % cmip trend
            cmip_col_mean_temp = squeeze(cmip_col_mean(x,y,:));
            tr = polyfit(time,cmip_col_mean_temp,1);
            cmip_trend(x,y) = tr(1).*365.25.*10;
        end
    end
end

%% decompose gobai reconstruction
% % define path
% gobai_decomp_path = [fpaths.param_path 'GOBAI/' base_grid '/FFNN/c' ...
%         num2str(num_clusters) '_' file_date float_file_ext '/train' ...
%         num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
%         num2str(100*test_ratio) '/' float_ext glodap_ext ctd_ext];
% % initiate parallel pool
% tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;
% parfor z = 1%:length(depth)
%     % load gobai for depth level
%     gobai = squeeze(ncread(osse_filepath,param_props.file_name,...
%         [1 1 1 1],[Inf Inf z Inf]));
%     gobai_init = nan(size(gobai,1),size(gobai,2));
%     gobai_seas = nan(size(gobai));
%     gobai_tr = nan(size(gobai));
%     gobai_resid = nan(size(gobai));
%     for x = 1:length(lon)
%         for y = 1:length(lat)
%             if any(~isnan(gobai(x,y,:)))
%                 gobai_temp = squeeze(gobai(x,y,:));
%                 % fit model
%                 periods = [365.2425 365.2425/2]; % annual and semi-annual
%                 [yf,yr,coeff] = ... % fit model
%                     leastsq2(time,gobai_temp,time(1),2,periods);
%                 gobai_init(x,y) = coeff(1);
%                 gobai_tr(x,y,:) = coeff(2).*(time-time(1));
%                 gobai_seas(x,y,:) = ...
%                     coeff(3).*cos((2*pi/periods(1))*(time-time(1))) + ...
%                     coeff(4).*sin((2*pi/periods(1))*(time-time(1))) + ...
%                     coeff(5).*cos((2*pi/periods(2))*(time-time(1))) + ...
%                     coeff(6).*sin((2*pi/periods(2))*(time-time(1)));
%                 gobai_resid(x,y,:) = yr;
%             end
%         end
%     end
%     % save temporary decomposition file
%     parsave([gobai_decomp_path '/gobai_decomp_' num2str(depth(z))],...
%         gobai_init,'gobai_init',gobai_tr,'gobai_tr',...
%         gobai_seas,'gobai_seas',gobai_resid,'gobai_resid');
% end
% 
% 
% % % assemble gobai decomposition in NetCDF
% % for z = 1:length(depth)
% %     load([gobai_decomp_path '/gobai_decomp_' num2str(depth(z))],'gobai_init','gobai_tr',...
% %         'gobai_seas','gobai_resid');
% %     if z ==1
% %         ncsave_4d([gobai_decomp_path '/gobai_decomp.nc'],...
% %             {'lon' lon 'longitude' 'degrees east'},...
% %             {'lat' lat 'latitude' 'degrees north'},...
% %             {'depth' depth(z) 'depth' 'meters'},...
% %             {'time' time 'time' 'days since 0000-01-01'},...
% %             {var_name var_temp long_name units});
% %     else
% %         ncwrite(nc_filepath,'time',time(t),t);
% %         ncwrite(nc_filepath,var_name,var_temp,[1 1 1 t]);
% %     end
% %     nccreate();
% % end
% 
% % end parallel session
% delete(gcp('nocreate'));
% keyboard

%% decompose cmip reconstruction
% % define path
% cmip_decomp_path = [fpaths.model_path base_grid '/combined/regridded'];
% % initiate parallel pool
% tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;
% parfor z = 1%:length(depth)
%     % load cmip for depth level
%     cmip = squeeze(ncread(cmip_filepath,param_props.file_name,...
%         [1 1 1 1],[Inf Inf z Inf]));
%     cmip_init = nan(size(cmip,1),size(cmip,2));
%     cmip_seas = nan(size(cmip));
%     cmip_tr = nan(size(cmip));
%     cmip_resid = nan(size(cmip));
%     for x = 1:length(lon)
%         for y = 1:length(lat)
%             if any(~isnan(cmip(x,y,:)))
%                 cmip_temp = squeeze(cmip(x,y,:));
%                 % fit model
%                 periods = [365.2425 365.2425/2]; % annual and semi-annual
%                 [yf,yr,coeff] = ... % fit model
%                     leastsq2(time,cmip_temp,time(1),2,periods);
%                 cmip_init(x,y) = coeff(1);
%                 cmip_tr(x,y,:) = coeff(2).*(time-time(1));
%                 cmip_seas(x,y,:) = ...
%                     coeff(3).*cos((2*pi/periods(1))*(time-time(1))) + ...
%                     coeff(4).*sin((2*pi/periods(1))*(time-time(1))) + ...
%                     coeff(5).*cos((2*pi/periods(2))*(time-time(1))) + ...
%                     coeff(6).*sin((2*pi/periods(2))*(time-time(1)));
%                 cmip_resid(x,y,:) = yr;
%             end
%         end
%     end
%     % save temporary decomposition file
%     parsave([cmip_decomp_path '/cmip_decomp_' num2str(depth(z))],...
%         cmip_init,'cmip_init',cmip_tr,'cmip_tr',...
%         cmip_seas,'cmip_seas',cmip_resid,'cmip_resid');
% end
% 
% % end parallel session
% delete(gcp('nocreate'));
% keyboard

%% plot timeseries
for r = 1:length(regions)
    % gobai inventory decomposition
    periods = [365.2425 365.2425/2]; % annual and semi-annual
    [yf_gobai,yr_gobai,coeff_gobai] = ... % fit model
        leastsq2(time,gobai_inv_both.(regions{r}),time(1),2,periods);
    gobai_mean = mean(yf_gobai);
    gobai_tr = coeff_gobai(2).*(time-time(1))-mean(coeff_gobai(2).*(time-time(1)));
    gobai_seas = ...
        coeff_gobai(3).*cos((2*pi/periods(1))*(time-time(1))) + ...
        coeff_gobai(4).*sin((2*pi/periods(1))*(time-time(1))) + ...
        coeff_gobai(5).*cos((2*pi/periods(2))*(time-time(1))) + ...
        coeff_gobai(6).*sin((2*pi/periods(2))*(time-time(1)));
    gobai_resid = yr_gobai;
    % cmip inventory decomposition
    periods = [365.2425 365.2425/2]; % annual and semi-annual
    [yf_cmip,yr_cmip,coeff_cmip] = ... % fit model
        leastsq2(time,cmip_inv_both.(regions{r}),time(1),2,periods);
    cmip_mean = mean(yf_cmip);
    cmip_tr = coeff_cmip(2).*(time-time(1))-mean(coeff_cmip(2).*(time-time(1)));
    cmip_seas = ...
        coeff_cmip(3).*cos((2*pi/periods(1))*(time-time(1))) + ...
        coeff_cmip(4).*sin((2*pi/periods(1))*(time-time(1))) + ...
        coeff_cmip(5).*cos((2*pi/periods(2))*(time-time(1))) + ...
        coeff_cmip(6).*sin((2*pi/periods(2))*(time-time(1)));
    cmip_resid = yr_cmip;
    % for everything
    corr_all = corr(gobai_inv_both.(regions{r}),cmip_inv_both.(regions{r}));
    rmsd_all = sqrt(mean((gobai_inv_both.(regions{r})-cmip_inv_both.(regions{r})).^2));
    std_gobai_all = std(gobai_inv_both.(regions{r}));
    std_cmip_all = std(cmip_inv_both.(regions{r}));
    figure('Position',[100 100 800 400]); hold on;
    title([regions{r} ', corr. coeff. = ' num2str(round(corr_all,3)) ', '...
        'RMSD = ' num2str(round(rmsd_all,3)) ' Pmol']);
    plot(datenum(1950,1,1)+time,cmip_inv_both.(regions{r}),...
        datenum(1950,1,1)+time,gobai_inv_both.(regions{r}),'LineWidth',2);
    ylabel([param_props.label ' inventory (Pmol)']);
    datetick('x');
    export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
            '/total_comp_' regions{r} '_' float_ext glodap_ext ctd_ext '.png']); close;
    % for seasonal
    corr_seas = corr(gobai_seas,cmip_seas);
    rmsd_seas = sqrt(mean((gobai_seas-cmip_seas).^2));
    std_gobai_seas = std(gobai_seas);
    std_cmip_seas = std(cmip_seas);
    figure('Position',[100 100 800 400]); hold on;
    title([regions{r} ', corr. coeff. = ' num2str(round(corr_seas,3)) ', '...
        'RMSD = ' num2str(round(rmsd_seas,3)) ' Pmol']);
    plot(datenum(1950,1,1)+time,cmip_seas,datenum(1950,1,1)+time,gobai_seas,'LineWidth',2);
    ylabel([param_props.label ' inventory (Pmol)']);
    datetick('x');
    legend({base_grid ['GOBAI-' param_props.dir_name '_{(' base_grid ')}']});
    export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
            '/seas_comp_' regions{r} '_' float_ext glodap_ext ctd_ext '.png']); close;
    % for trend
    corr_tr = corr(gobai_tr,cmip_tr);
    rmsd_tr = sqrt(mean((gobai_tr-cmip_tr).^2));
    std_gobai_tr = std(gobai_tr);
    std_cmip_tr = std(cmip_tr);
    figure('Position',[100 100 800 400]); hold on;
    title([regions{r} ', corr. coeff. = ' num2str(round(corr_tr,3)) ', '...
        'RMSD = ' num2str(round(rmsd_tr,3)) ' Pmol']);
    plot(datenum(1950,1,1)+time,cmip_tr,datenum(1950,1,1)+time,gobai_tr,'LineWidth',2);
    ylabel([param_props.label ' inventory (Pmol)']);
    datetick('x');
    legend({base_grid ['GOBAI-' param_props.dir_name '_{(' base_grid ')}']});
    export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
            '/trend_comp_' regions{r} '_' float_ext glodap_ext ctd_ext '.png']); close;
    % for residual
    corr_resid = corr(gobai_resid,cmip_resid);
    rmsd_resid = sqrt(mean((gobai_resid-cmip_resid).^2));
    std_gobai_resid = std(gobai_resid);
    std_cmip_resid = std(cmip_resid);
    figure('Position',[100 100 800 400]); hold on;
    title([regions{r} ', corr. coeff. = ' num2str(round(corr_resid,3)) ', '...
        'RMSD = ' num2str(round(rmsd_resid,3)) ' Pmol']);
    plot(datenum(1950,1,1)+time,cmip_resid,datenum(1950,1,1)+time,gobai_resid,'LineWidth',2);
    ylabel([param_props.label ' inventory (Pmol)']);
    datetick('x');
    legend({base_grid ['GOBAI-' param_props.dir_name '_{(' base_grid ')}']});
    export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
            '/resid_comp_' regions{r} '_' float_ext glodap_ext ctd_ext '.png']); close;
end

%% plot taylor diagram
% theta_seas = acos(corr_seas);  % angle from correlation
% r_seas = std_gobai_seas/std_cmip_seas;   % normalized standard deviation
% theta_tr = acos(corr_tr);  % angle from correlation
% r_tr = std_gobai_tr/std_cmip_tr;   % normalized standard deviation
% theta_resid = acos(corr_resid);  % angle from correlation
% r_resid = std_gobai_resid/std_cmip_resid;   % normalized standard deviation
% % draw correlation grid
% polaraxes; hold on;
% corr_ticks = [0.2 0.4 0.6 0.8 0.9 0.95 1];
% for c = corr_ticks
%     t = acos(c);
%     polarplot([t t],[0 2],'k:');
% end
% % centered RMS contours
% E = 0.5:0.5:2;
% theta_grid = linspace(0,pi/2,200);
% for e = E
%     r_rms = cos(theta_grid) + sqrt(e^2 - sin(theta_grid).^2);
%     polarplot(theta_grid,r_rms,'--','LineWidth',1)
% end
% % draw standard deviation circles
% rmax = 2;
% th = linspace(0,pi/2,200);
% for s = 0.5:0.5:rmax
%     polarplot(th,s*ones(size(th)),'k:')
% end
% % plot reference point
% polarplot(0,1,'kp','MarkerFaceColor','k','MarkerSize',12);
% % plot results
% polarplot(theta_seas,r_seas,'o','MarkerSize',8);
% polarplot(theta_tr,r_tr,'^','MarkerSize',8);
% polarplot(theta_resid,r_resid,'s','MarkerSize',8);
% % axis limits
% rlim([0 rmax]);
% thetalim([0 90]);
% % legend
% legend(['Seasonal' 'Trend' 'Residual']);

%% plot global gobai trends
figure('visible','on');
worldmap('world');
set(gca,'fontsize',12);
title(['Reconstructed ' base_grid]);
mlabel off; plabel off;
pcolorm(lat,[lon;lon(end)+1],[gobai_trend;gobai_trend(end,:)]');
plot_land('map');
c=colorbar('Location','southoutside');
clim([-param_props.edges(end)/50 param_props.edges(end)/50]);
c.Label.String = [' Column ' param_props.label ...
    ' Trend (mol m^{-2} dec^{-1})'];
colormap(cmocean('balance'));
export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
        '/gobai_col_mean_trend_' float_ext glodap_ext ctd_ext '.png'],'-transparent'); close;
close

%% plot global cmip trends
figure('visible','on');
worldmap('world');
set(gca,'fontsize',12);
title(base_grid);
mlabel off; plabel off;
pcolorm(lat,[lon;lon(end)+1],[cmip_trend;cmip_trend(end,:)]');
plot_land('map');
c=colorbar('Location','southoutside');
clim([-param_props.edges(end)/50 param_props.edges(end)/50]);
c.Label.String = [' Column ' param_props.label ...
    ' Trend (mol m^{-2} dec^{-1})'];
colormap(cmocean('balance'));
export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
        '/model_col_mean_trend_' float_ext glodap_ext ctd_ext '.png'],'-transparent'); close;
close

%% save global means
if ~isfolder([param_props.dir_name '/Data/' base_grid]); mkdir([param_props.dir_name '/Data/' base_grid '/' rlz '_gr']); end
save([param_props.dir_name '/Data/' base_grid '/' rlz '_gr/statistics_' ...
    float_ext glodap_ext ctd_ext '.mat'],...
    'gobai_mean','gobai_inv','gobai_depth_mean','cmip_mean','cmip_inv','cmip_depth_mean');

%% plot global mean timeseries
% figure('visible','on');
% plot(datenum(1950,1,1)+time,cmip_mean,datenum(1950,1,1)+time,gobai_mean,'LineWidth',2);
% legend({base_grid ['GOBAI-' param_props.dir_name '_{(' base_grid ')}']});
% datetick('x','keeplimits');
% if ~isfolder([param_props.dir_name '/Figures/' base_grid '/' rlz '_gr'])
%     mkdir([param_props.dir_name '/Figures/' base_grid '/' rlz '_gr']);
% end
% export_fig(gcf,[param_props.dir_name '/Figures/' base_grid '/' rlz '_gr' ...
%         '/timeseries_' float_ext glodap_ext ctd_ext '.png'],'-transparent'); close;

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
        '/profile_' float_ext glodap_ext ctd_ext '.png'],'-transparent'); close;

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
        '/profile_delta_' float_ext glodap_ext ctd_ext '.png'],'-transparent'); close;

end
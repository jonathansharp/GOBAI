function gobai_global_timeseries(base_grid,param_props,fpaths,num_clusters,file_date,...
    float_file_ext,train_ratio,val_ratio,test_ratio,flt,gld,ctd,start_year,end_year)

%% define dataset extensions
if flt == 1; float_ext = 'f'; else float_ext = ''; end
if gld == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if ctd == 1; ctd_ext = 'w'; else ctd_ext = ''; end

%% define gobai file name
path_ext = ['GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) ...
        '_' file_date float_file_ext '/train' num2str(100*train_ratio) ...
        '_val' num2str(100*val_ratio) '_test' num2str(100*test_ratio) ...
        '/' float_ext glodap_ext ctd_ext '/'];
gobai_alg_dir = [fpaths.param_path path_ext];
filename = [gobai_alg_dir 'gobai-' param_props.file_name '.nc'];

%% download dimensions
GOBAI.lon = ncread(filename,'lon');
GOBAI.lat = ncread(filename,'lat');
GOBAI.pres = ncread(filename,'pres');
GOBAI.time = ncread(filename,'time');

%% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);

% define mask
if ~isfile([gobai_alg_dir 'mask.mat'])
    mask = true(size(GOBAI.vol));
    for t = 1:length(GOBAI.time)
        tic
        gobai_tmp = ncread(filename,'o2',[1 1 1 t],[Inf Inf Inf 1]);
        mask(isnan(gobai_tmp)) = false;
        toc
    end
    save([gobai_alg_dir 'mask.mat'],'mask');
else
    load([gobai_alg_dir 'mask.mat'],'mask');
end
GOBAI.vol(~mask) = NaN;
% download oxygen and calculate global mean and inventory
if ~isfile([gobai_alg_dir 'global_mean.mat']) || ...
        ~isfile([gobai_alg_dir 'global_inv.mat'])
    global_mean = nan(size(GOBAI.time));
    global_inv = nan(size(GOBAI.time));
    TS = load_RFROM_dim(fpaths.temp_path,'v2.2',start_year,end_year);
    cnt = 1;
    for m = 1:length(TS.months)
        % determine number of weeks in RFROM file
        nc_atts = ncinfo([fpaths.temp_path 'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
            gobai_tmp = ncread(filename,'o2',[1 1 1 cnt],[Inf Inf Inf 1]);
            temp_cns_tmp = ncread([fpaths.temp_path 'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_' ...
                num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                    'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
            sal_abs_tmp = ncread([fpaths.sal_path 'RFROM_SAL_v2.2/RFROMV22_SAL_STABLE_' ...
                num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                    'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
            dens_tmp = gsw_rho(sal_abs_tmp,temp_cns_tmp,...
                repmat(permute(GOBAI.pres,[3 2 1]),length(GOBAI.lon),length(GOBAI.lat),1));
            % mask out values
            gobai_tmp(~mask) = NaN; dens_tmp(~mask) = NaN;
            % inventory in Pmol (umol/kg * kg/m3 * m3)
            global_inv(cnt) = sum(gobai_tmp(:).*dens_tmp(:).*GOBAI.vol(:),'omitnan')./10^6./10^15;
            % average concentration in umol/kg
            global_mean(cnt) = sum(gobai_tmp(:).*...
                GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
            % increase counter
            cnt = cnt + 1;
        end
    end
    % save global inventory and mean
    save([gobai_alg_dir 'global_inv'],'global_inv');
    save([gobai_alg_dir 'global_mean'],'global_mean');
else
    % load global inventory and mean
    load([gobai_alg_dir 'global_inv'],'global_inv');
    load([gobai_alg_dir 'global_mean'],'global_mean');
end
% plot global inventory
h=figure; hold on; set(gcf,'Position',[100 100 800 400]);
plot(datenum(1950,0,0)+double(GOBAI.time),global_inv,'k','LineWidth',2);
global_inv_smoothed = smooth(global_inv,52);
plot(datenum(1950,0,0)+double(GOBAI.time(27:end-26)),...
    global_inv_smoothed(27:end-26),'k','LineWidth',1,'LineStyle',':');
datetick('x'); ylabel('0-2000 dbar Oxygen Inventory (Pmol)');
export_fig(h,[gobai_alg_dir 'global_inv.png'],'-transparent'); close
% plot global mean
h=figure; hold on; set(gcf,'Position',[100 100 800 400]);
plot(datenum(1950,0,0)+double(GOBAI.time),global_mean,'k','LineWidth',2);
global_mean_smoothed = smooth(global_mean,52);
plot(datenum(1950,0,0)+double(GOBAI.time(27:end-26)),...
    global_mean_smoothed(27:end-26),'k','LineWidth',1,'LineStyle',':');
datetick('x'); ylabel('0-2000 dbar Mean [O_{2}] (Pmol)');
export_fig(h,[gobai_alg_dir 'global_mean.png'],'-transparent'); close

end
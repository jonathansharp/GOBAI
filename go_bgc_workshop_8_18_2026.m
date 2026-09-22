tic

% set lat/lon limits
latlim=[30 40];
lonlim=[140 150];

% select pH profiles
[so_ph_floats,so_ph_float_profs] = select_profiles(lonlim,latlim,[],[],...
    'sensor','PH_IN_SITU_TOTAL',...
    'mode','D'); % include only pH floats

% display the number of matching floats and profiles
disp(['# of matching profiles: ' num2str(sum(cellfun('length',...
    so_ph_float_profs)))]);
disp(['# of matching floats: ' num2str(length(so_ph_floats))]);

% get and save data
[data,mdata] = load_float_data(so_ph_floats,... % specify WMO number
    {'PSAL','TEMP','PH_IN_SITU_TOTAL'}); % specify variables
% save('go_bgc_ws_data','data','mdata','-v7.3');
% load('go_bgc_ws_data');

% assemble into vectors
float_nums = fieldnames(data);
vars = fieldnames(data.(float_nums{1}));
for v = 1:length(vars)
    data_array.(vars{v}) = [];
end
for f = 1:length(float_nums)
    vars = fieldnames(data.(float_nums{f}));
    idx = (data.(float_nums{f}).PH_IN_SITU_TOTAL_ADJUSTED_QC == 1 | ...
        data.(float_nums{f}).PH_IN_SITU_TOTAL_ADJUSTED_QC == 2) & ...
        data.(float_nums{f}).PRES_ADJUSTED > 0 & data.(float_nums{f}).PRES_ADJUSTED <= 100;
    for v = 1:length(vars)
        temp_var = data.(float_nums{f}).(vars{v})(idx);
        data_array.(vars{v}) = [data_array.(vars{v});temp_var];
    end
end
clear data

% plot time histogram
figure('units','inches','Position',[1 1 13.33 5]);
histogram(data_array.TIME);
datetick('x');
exportgraphics(gcf,'go_bgc_ws/time_histogram.png'); close

% plot pH histogram
figure('units','inches','Position',[1 1 13.33 5]);
histogram(data_array.PH_IN_SITU_TOTAL_ADJUSTED);
exportgraphics(gcf,'go_bgc_ws/histogram.png'); close

% plot scatter
figure('units','inches','Position',[1 1 13.33 7.5]);
worldmap(latlim,lonlim); hold on;
idx = data_array.PRES_ADJUSTED > 45 & data_array.PRES_ADJUSTED < 55;
scatterm(data_array.LATITUDE(idx),data_array.LONGITUDE(idx),100,...
    data_array.PH_IN_SITU_TOTAL_ADJUSTED(idx),...
    'filled','MarkerEdgeColor','k','Linewidth',2);
plot_land('map');
colorbar; caxis([7.9 8.1]);
exportgraphics(gcf,'go_bgc_ws/pH_scatter.png'); close

% load temp and sal
temp_file = '/raid/Data/RFROM/RFROM_TEMP_v2.2-2025/';
sal_file = '/raid/Data/RFROM/RFROM_SAL_v2.2-2025/';
lon = ncread([temp_file 'RFROMV22_TEMP_STABLE_1993_01.nc'],'longitude');
lon_idx = find(lon >= lonlim(1) & lon <= lonlim(2));
lat = ncread([temp_file 'RFROMV22_TEMP_STABLE_1993_01.nc'],'latitude');
lat_idx = find(lat >= latlim(1) & lat <= latlim(2));
pres = ncread([temp_file 'RFROMV22_TEMP_STABLE_1993_01.nc'],'mean_pressure');
pres_idx = find(pres >= 0 & pres <= 100);
temp = []; sal = []; time = [];
for y = 2021
    for m = 1:12
        temp_temp = ncread([temp_file 'RFROMV22_TEMP_STABLE_' num2str(y) ...
            '_' sprintf('%02d',m) '.nc'],'ocean_temperature',...
            [lon_idx(1) lat_idx(1) pres_idx(1) 1],...
            [length(lon_idx) length(lat_idx) length(pres_idx) Inf]);
        temp = cat(4,temp,temp_temp);
        sal_temp = ncread([sal_file 'RFROMV22_SAL_STABLE_' num2str(y) ...
            '_' sprintf('%02d',m) '.nc'],'ocean_salinity',...
            [lon_idx(1) lat_idx(1) pres_idx(1) 1],...
            [length(lon_idx) length(lat_idx) length(pres_idx) Inf]);
        sal = cat(4,sal,sal_temp);
        time = [time;ncread([sal_file 'RFROMV22_SAL_STABLE_' num2str(y) ...
            '_' sprintf('%02d',m) '.nc'],'time')];
    end
end
clear temp_temp sal_temp

% plot temp
figure('units','inches','Position',[1 1 13.33 7.5]);
worldmap(latlim,lonlim); hold on;
contourfm(double(lat(lat_idx)),double(lon(lon_idx)),double(temp(:,:,6,30)'));
plot_land('map');
colorbar; colormap(cmocean('thermal'));
exportgraphics(gcf,'go_bgc_ws/temp.png'); close

% plot sal
figure('units','inches','Position',[1 1 13.33 7.5]);
worldmap(latlim,lonlim); hold on;
contourfm(double(lat(lat_idx)),double(lon(lon_idx)),double(sal(:,:,6,30)'));
plot_land('map');
colorbar; colormap(cmocean('haline'));
exportgraphics(gcf,'go_bgc_ws/sal.png'); close

%% train rf
numTrees = 100;
minleaf = 1;
X = [data_array.TEMP_ADJUSTED data_array.PSAL_ADJUSTED ...
    data_array.PRES_ADJUSTED];
Y = data_array.PH_IN_SITU_TOTAL_ADJUSTED;
rfModel = TreeBagger(numTrees, X, Y, 'Method', 'regression','MinLeafSize',minleaf);

% plot pH map
figure('units','inches','Position',[1 1 13.33 7.5]);
worldmap(latlim,lonlim); hold on;
cns_temp_temp = temp(:,:,6,30);
abs_sal_temp = sal(:,:,6,30);
pres_temp = repmat(100,40,40);
[lon_temp,lat_temp] = ndgrid(lon(lon_idx),lat(lon_idx));
temp_temp = gsw_t_from_CT(abs_sal_temp,cns_temp_temp,pres(6));
sal_temp = gsw_SP_from_SA(abs_sal_temp,pres(6),lon_temp,lat_temp);
time_temp = repmat(double(time(30)+datenum(1950,1,1)),40,40);
% pH_map = predict(rfModel,[temp_temp(:),sal_temp(:),pres_temp(:),time_temp(:)]);
pH_map = predict(rfModel,[temp_temp(:),sal_temp(:),pres_temp(:)]);
pH_map = reshape(pH_map,40,40);
idx = isnan(temp(:,:,6,30)) & isnan(sal(:,:,6,30));
pH_map(idx) = NaN;
contourfm(double(lat(lat_idx)),double(lon(lon_idx)),pH_map',7.9:0.02:8.15);
colorbar; caxis([7.9 8.1]);
data_array.DATE = datevec(data_array.TIME);
idx = data_array.PRES_ADJUSTED > 45 & data_array.PRES_ADJUSTED < 55 & ...
    data_array.DATE(:,2) == 7;
scatterm(data_array.LATITUDE(idx),data_array.LONGITUDE(idx),100,...
    data_array.PH_IN_SITU_TOTAL_ADJUSTED(idx),...
    'filled','MarkerEdgeColor','k','Linewidth',2);
plot_land('map');
colorbar;
exportgraphics(gcf,'go_bgc_ws/pH_map_no_t.png'); close

% =========================================================================
% K-Fold Cross-Validation for Random Forest Regression
% =========================================================================

% 1. Setup Parameters
K = 5;                  % Number of folds

% Note: Replace X and Y with your actual data matrices
% X: [N x P] Feature matrix (N samples, P features)
% Y: [N x 1] Target vector

% 2. Create K-Fold Partition
rng(42); % Set random seed for reproducibility
cv = cvpartition(length(Y), 'KFold', K);

% Preallocate arrays to store fold metrics
rmse_folds = zeros(K, 1);
mae_folds  = zeros(K, 1);
r2_folds   = zeros(K, 1);

% Preallocate vector for out-of-fold predictions
y_pred_cv = zeros(size(Y));

% 3. Run K-Fold Evaluation Loop
fprintf('Starting %d-Fold Cross-Validation...\n', K);

for k = 1:K
    % Get train and validation logical indices for current fold
    trainIdx = cv.training(k);
    valIdx   = cv.test(k);

    % Split data
    X_train = X(trainIdx, :);
    y_train = Y(trainIdx);
    X_val   = X(valIdx, :);
    y_val   = Y(valIdx);

    % Train Random Forest Regressor
    % 'Method', 'Bag' creates a Random Forest ensemble of decision trees
    rfModel = TreeBagger(numTrees, X_train, y_train, ...
        'Method', 'regression','MinLeafSize',minleaf);

    % Predict on validation set
    y_pred = predict(rfModel, X_val);
    y_pred_cv(valIdx) = y_pred; % Store out-of-fold predictions

    % Compute Metrics for this Fold
    rmse_folds(k) = sqrt(mean((y_val - y_pred).^2));
    mae_folds(k)  = mean(abs(y_val - y_pred));

    ss_tot = sum((y_val - mean(y_val)).^2);
    ss_res = sum((y_val - y_pred).^2);
    r2_folds(k)   = 1 - (ss_res / ss_tot);

    fprintf('  Fold %d/%d - RMSE: %.4f | MAE: %.4f | R^2: %.4f\n', ...
        k, K, rmse_folds(k), mae_folds(k), r2_folds(k));
end

% 4. Summary Results Across Folds
fprintf('\n================ CV SUMMARY ================\n');
fprintf('Mean RMSE : %.4f (± %.4f)\n', mean(rmse_folds), std(rmse_folds));
fprintf('Mean MAE  : %.4f (± %.4f)\n', mean(mae_folds), std(mae_folds));
fprintf('Mean R^2   : %.4f (± %.4f)\n', mean(r2_folds), std(r2_folds));

% 5. Visualization: Actual vs. Predicted (Out-of-Fold)
figure;
scatter(Y, y_pred_cv, 25, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
plot([min(Y), max(Y)], [min(Y), max(Y)], 'r--', 'LineWidth', 1.5);
xlabel('Actual pH');
ylabel('Predicted pH (Out-of-Fold)');
title(sprintf('%d-Fold CV: Actual vs. Predicted (R^2 = %.3f)', K, mean(r2_folds)));
grid on;
exportgraphics(gcf,'go_bgc_ws/pH_rmse_no_t.png'); close

% plot timeseries

%%
toc
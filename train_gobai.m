% train_gobai
%
% DESCRIPTION:
% This function trains algorithms in GMM clusters for GOBAI.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 3/10/2026

function train_gobai(alg_type,param_props,base_grid,file_date,...
    float_file_ext,num_clusters,variables,thresh,numWorkers_train,...
    snap_date,flt,gld,ctd,varargin)

%% define dataset extensions
if flt == 1; float_ext = 'f'; else float_ext = ''; end
if gld == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if ctd == 1; ctd_ext = 'w'; else ctd_ext = ''; end

%% set defaults and process optional input arguments
num_folds = 1;
glodap_only = 0;
data_per = 1;
coverage = '';
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'num_folds')
        num_folds = varargin{i+1};
    elseif strcmpi(varargin{i}, 'glodap_only')
        glodap_only = varargin{i+1};
    elseif strcmpi(varargin{i}, 'reduce_data')
        data_per = varargin{i+1};
    elseif strcmpi(varargin{i}, 'coverage')
        coverage = ['_' num2str(varargin{i+1}) '_coverage'];
    end
end

%% process necessary input arguments for model parameters
% pre-allocate
train_ratio = NaN;
val_ratio = NaN;
test_ratio = NaN;
numtrees = NaN;
minLeafSize = NaN;
numstumps = NaN;
numbins = NaN;
% process based on algorithm type
if strcmp(alg_type,'FFNN')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'train_ratio')
            train_ratio = varargin{i+1};
        elseif strcmpi(varargin{i}, 'val_ratio')
            val_ratio = varargin{i+1};
        elseif strcmpi(varargin{i}, 'test_ratio')
            test_ratio = varargin{i+1};
        end
    end
elseif strcmp(alg_type,'RFR')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'numtrees')
            numtrees = varargin{i+1};
        elseif strcmpi(varargin{i}, 'minLeafSize')
            minLeafSize = varargin{i+1};
        end
    end
elseif strcmp(alg_type,'GBM')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'numstumps')
            numstumps = varargin{i+1};
        elseif strcmpi(varargin{i}, 'numbins')
            numbins = varargin{i+1};
        end
    end
else
    disp('"alg_type" must be "FFNN", "RFR", or "GBM"')
end

%% process date
date_str = num2str(snap_date);

%% directory base
if strcmp(alg_type,'FFNN')
    dir_base = create_dir_base(alg_type,{num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
elseif strcmp(alg_type,'RFR')
    dir_base = create_dir_base(alg_type,{num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
elseif strcmp(alg_type,'GBM')
    dir_base = create_dir_base(alg_type,{num_clusters;file_date;...
        float_file_ext;numstumps;numbins});
end

%% load data
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    load([param_props.dir_name '/Data/processed_all_' ...
        param_props.file_name '_data_' float_ext glodap_ext ctd_ext ...
        '_' file_date float_file_ext '.mat'],'all_data');
else % for cmip models
    load([param_props.dir_name '/Data/' base_grid '_' ...
        param_props.file_name '_data_' float_ext glodap_ext ctd_ext ...
        '_' file_date float_file_ext coverage '.mat'],'all_data');
end

%% load data clusters
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    load([param_props.dir_name '/Data/GMM/' ...
        'all_data_clusters_' num2str(num_clusters) '_' ...
        float_ext glodap_ext ctd_ext '_' file_date float_file_ext ...
        '.mat'],'all_data_clusters');
else
    load([param_props.dir_name '/Data/GMM_' base_grid '/' ...
        'all_data_clusters_' num2str(num_clusters) '_' ...
        float_ext glodap_ext ctd_ext '_' ...
        file_date float_file_ext coverage '.mat'],'all_data_clusters');
end

%% load data cluster indices (for k-fold testing)
if num_folds > 1
    load([param_props.dir_name '/Data/k_fold_data_indices_' num2str(num_clusters) ...
        '_' num2str(num_folds) '_' float_ext glodap_ext ctd_ext '_' ...
        file_date float_file_ext '.mat'],'num_folds','train_idx','test_idx');
else
    % set NaNs so parloops run
    train_idx = NaN; test_idx = NaN;
end

%% remove float data for GLODAP only test
if glodap_only
    glodap_idx = all_data.platform < 10^6;
    vars = fieldnames(all_data);
    for v = 1:length(vars)
        all_data.(vars{v}) = all_data.(vars{v})(glodap_idx);
    end
end

%% create directory and file names
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    alg_dir = [param_props.dir_name '/Models/' dir_base];
else
    alg_dir = [param_props.dir_name '/Models/' base_grid '/' dir_base];
end
alg_fnames = cell(num_folds,num_clusters);
for f = 1:num_folds
    for c = 1:num_clusters
        if num_folds > 1
            alg_fnames(f,c) = {[alg_type '_' param_props.file_name '_C' ...
                num2str(c) '_F' num2str(f) '_test_' float_ext glodap_ext ctd_ext]};
        else
            alg_fnames(c) = {[alg_type '_' param_props.file_name '_C' ...
                num2str(c) '_' float_ext glodap_ext ctd_ext]};
        end
    end
end

%% variables for k-fold test
if num_folds > 1
    kfold_dir = [param_props.dir_name '/KFold/' alg_type '/c' num2str(num_clusters) '_' file_date float_file_ext];
    fig_dir = [param_props.dir_name '/Figures/KFold/' alg_type '/c' num2str(num_clusters) '_' file_date float_file_ext];
    if strcmp(alg_type,'RFR')
        kfold_name = ['RFR_output_tr' num2str(numtrees) '_lf' num2str(minLeafSize)];
        fig_name_1 = ['k_fold_comparison_tr' num2str(numtrees) '_lf' ...
            num2str(minLeafSize) '_' float_ext glodap_ext ctd_ext '.png'];
        fig_name_2 = ['k_fold_spatial_comparison_tr' num2str(numtrees) ...
            '_lf' num2str(minLeafSize) '_' float_ext glodap_ext ctd_ext '.png'];
    elseif strcmp(alg_type,'FFNN')
        kfold_name = ['FFNN_output_train' num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' num2str(100*val_ratio)];
        fig_name_1 = ['k_fold_comparison_train' num2str(100*train_ratio) ...
            '_val' num2str(100*val_ratio) '_test' num2str(100*val_ratio) ...
            '_' float_ext glodap_ext ctd_ext '.png'];
        fig_name_2 = ['k_fold_spatial_comparison_train' num2str(100*train_ratio) ...
            '_val' num2str(100*val_ratio) '_test' num2str(100*val_ratio) ...
            '_' float_ext glodap_ext ctd_ext '.png'];
    elseif strcmp(alg_type,'GBM')
        kfold_name = ['GBM_output_tr' num2str(numstumps) '_bin' num2str(numbins)];
        fig_name_1 = ['k_fold_comparison_tr' num2str(numstumps) '_bin' ...
            num2str(numbins) '_' float_ext glodap_ext ctd_ext '.png'];
        fig_name_2 = ['k_fold_spatial_comparison_tr' num2str(numstumps) ...
            '_bin' num2str(numbins) '_' float_ext glodap_ext ctd_ext '.png'];
    end
end

%% fit algorithms

% start timing training
tStart = tic;

% parameter that matches fold number with cluster number
folds = repelem(1:num_folds,1,num_clusters)';
clusters = repmat(1:num_clusters,1,num_folds)';

% fit models
% if num_folds == 1 & strcmp(alg_type,'FFNN')
    % set up parallel pool
    tic; parpool(numWorkers_train); fprintf('Pool initiation: '); toc;
    for cnt = 1:num_folds*num_clusters
        train_models(param_props,num_folds,...
            alg_dir,alg_fnames,variables,all_data,all_data_clusters,...
            train_idx,test_idx,data_per,alg_type,train_ratio,test_ratio,val_ratio,...
            numtrees,minLeafSize,numstumps,numbins,thresh,'yes',...
            folds(cnt),clusters(cnt),base_grid,coverage);
    end
% else
%     % set up parallel pool
%     tic; parpool(numWorkers_train); fprintf('Pool initiation: '); toc;
%     parfor cnt = 1:num_folds*num_clusters
%         train_models(param_props,num_folds,...
%             alg_dir,alg_fnames,variables,all_data,all_data_clusters,...
%             train_idx,test_idx,data_per,alg_type,train_ratio,test_ratio,val_ratio,...
%             numtrees,minLeafSize,numstumps,numbins,thresh,'no',...
%             folds(cnt),clusters(cnt));
%     end
% end

% end parallel session
delete(gcp('nocreate'));

% stop timing full script
if num_folds > 1
    fprintf([alg_type ' k-Fold Training for: ']);
else
    fprintf([alg_type ' Training for: ']);
end

% print elapsed time in minutes
tElapsed = toc(tStart);
disp(['Elapsed time is ' num2str(tElapsed/60) ' minutes.'])

%% calculate statistics for k-fold test
if num_folds > 1

% calculate weighted average over each cluster using probabilities
for f = 1:num_folds
    % pre-allocate probabilities
    probs_matrix = [];
    for c = 1:num_clusters
        % load output
        load([alg_dir '/' alg_fnames{f,c}],'output','output_ESPER','obs_index_test')
        alg_output.(['f' num2str(f)])(:,c) = output;
        esper_output.(['f' num2str(f)])(:,c) = output_ESPER;
        % assemble matrix of probabilities greater than the threshold (5%)
        probs_array = all_data_clusters.(['c' num2str(c)])(obs_index_test);
        probs_array(probs_array < thresh) = NaN;
        probs_matrix = [probs_matrix,probs_array];
    end
    alg_output.(['f' num2str(f) '_mean']) = ...
        double(sum(alg_output.(['f' num2str(f)]).*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
    esper_output.(['f' num2str(f) '_mean']) = ...
        double(sum(esper_output.(['f' num2str(f)]).*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
end
% clean up
clear f c
% aggregate output from all folds
alg_output.(['k_fold_test_' param_props.file_name]) = nan(size(all_data.(param_props.file_name)));
esper_output.(['k_fold_test_' param_props.file_name]) = nan(size(all_data.(param_props.file_name)));
for f = 1:num_folds
    % load output index
    load([alg_dir '/' alg_fnames{f,1}],'obs_index_test')
    % aggregate
    alg_output.(['k_fold_test_' param_props.file_name])(obs_index_test) = ...
        alg_output.(['f' num2str(f) '_mean']);
    esper_output.(['k_fold_test_' param_props.file_name])(obs_index_test) = ...
        esper_output.(['f' num2str(f) '_mean']);
end
% compare k-fold output to data
alg_output.k_fold_delta = alg_output.(['k_fold_test_' param_props.file_name]) - all_data.(param_props.file_name);
esper_output.k_fold_delta = esper_output.(['k_fold_test_' param_props.file_name]) - all_data.(param_props.file_name);
% calculate error stats for GOBAI algorithm
alg_mean_err = mean(alg_output.k_fold_delta,'omitnan');
alg_med_err = median(alg_output.k_fold_delta,'omitnan');
alg_rmse = sqrt(mean(alg_output.k_fold_delta.^2,'omitnan'));
alg_med_abs_err = median(abs(alg_output.k_fold_delta),'omitnan');
% calculate error stats for ESPER
esper_mean_err = mean(esper_output.k_fold_delta,'omitnan');
esper_med_err = median(esper_output.k_fold_delta,'omitnan');
esper_rmse = sqrt(mean(esper_output.k_fold_delta.^2,'omitnan'));
esper_med_abs_err = median(abs(esper_output.k_fold_delta),'omitnan');
% print error stats for GOBAI algorithm
fprintf([num2str(num_clusters) ' Clusters:\n']);
fprintf(['Mean Error (GOBAI) = ' num2str(alg_mean_err) ' ' param_props.units '\n']);
fprintf(['Median Error (GOBAI) = ' num2str(alg_med_err) ' ' param_props.units '\n']);
fprintf(['RMSE (GOBAI) = ' num2str(alg_rmse) ' ' param_props.units '\n']);
fprintf(['Median Abs. Error (GOBAI) = ' num2str(alg_med_abs_err) ' ' param_props.units '\n']);
% print error stats for ESPER
fprintf(['Mean Error (ESPER-NN) = ' num2str(esper_mean_err) ' ' param_props.units '\n']);
fprintf(['Median Error (ESPER-NN) = ' num2str(esper_med_err) ' ' param_props.units '\n']);
fprintf(['RMSE (ESPER-NN) = ' num2str(esper_rmse) ' ' param_props.units '\n']);
fprintf(['Median Abs. Error (ESPER-NN) = ' num2str(esper_med_abs_err) ' ' param_props.units '\n']);
% save predicted data
if ~isfolder([pwd '/' kfold_dir]); mkdir(kfold_dir); end
save([kfold_dir '/' kfold_name],'alg_output','alg_rmse',...
    'alg_med_err','alg_mean_err','esper_output','esper_rmse',...
    'esper_med_err','esper_mean_err','-v7.3');
clear alg_output alg_rmse alg_med_err alg_mean_err
clear esper_output esper_rmse esper_med_err esper_mean_err

%% plot histogram of errors for GOBAI algorithm
load([kfold_dir '/' kfold_name],'alg_output','alg_rmse');
figure('visible','on'); hold on;
set(gca,'fontsize',12);
set(gcf,'position',[100 100 600 400]);
[counts,bin_centers] = hist3([all_data.(param_props.file_name) alg_output.(['k_fold_test_' param_props.file_name])],...
    'Edges',{param_props.edges param_props.edges});
h=pcolor(bin_centers{1},bin_centers{2},counts');
plot([param_props.edges(1) param_props.edges(end)],[param_props.edges(1) param_props.edges(end)],'k--');
set(h,'EdgeColor','none');
xlim([param_props.edges(1) param_props.edges(end)]);
ylim([param_props.edges(1) param_props.edges(end)]);
xlabel(['Measured ' param_props.label ' ' param_props.units]);
ylabel([alg_type ' ' param_props.label ' ' param_props.units]);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
clim([1e0 1e5]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text(param_props.edges(1)+(2/5)*(param_props.edges(end)-param_props.edges(1)),...
    param_props.edges(1)+(1/10)*(param_props.edges(end)-param_props.edges(1)),...
    ['RMSE = ' num2str(round(alg_rmse,param_props.dec_points)) ' ' param_props.units],'fontsize',12);
if ~isfolder([pwd '/' fig_dir]); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/' fig_name_1]);
% clean up
clear counts bin_centers h p myColorMap
close

%% plot histogram of errors for ESPER
load([kfold_dir '/' kfold_name],'esper_output','esper_rmse');
figure('visible','on'); hold on;
set(gca,'fontsize',12);
set(gcf,'position',[100 100 600 400]);
[counts,bin_centers] = hist3([all_data.(param_props.file_name) esper_output.(['k_fold_test_' param_props.file_name])],...
    'Edges',{param_props.edges param_props.edges});
h=pcolor(bin_centers{1},bin_centers{2},counts');
plot([param_props.edges(1) param_props.edges(end)],[param_props.edges(1) param_props.edges(end)],'k--');
set(h,'EdgeColor','none');
xlim([param_props.edges(1) param_props.edges(end)]);
ylim([param_props.edges(1) param_props.edges(end)]);
xlabel(['Measured ' param_props.label ' ' param_props.units]);
ylabel(['ESPER-NN ' param_props.label ' ' param_props.units]);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
clim([1e0 1e5]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text(param_props.edges(1)+(2/5)*(param_props.edges(end)-param_props.edges(1)),...
    param_props.edges(1)+(1/10)*(param_props.edges(end)-param_props.edges(1)),...
    ['RMSE = ' num2str(round(esper_rmse,param_props.dec_points)) ' ' param_props.units],'fontsize',12);
if ~isfolder([pwd '/' fig_dir]); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/ESPER_' fig_name_1]);
% clean up
clear counts bin_centers h p myColorMap
close

%% plot gridded errors
% determine bin number of each test data point on 1 degree grid
lon_edges = -180:180; lon = -179.5:179.5;
lat_edges = -90:90; lat = -89.5:89.5;
[~,~,Xnum] = histcounts(all_data.longitude,lon_edges);
[~,~,Ynum] = histcounts(all_data.latitude,lat_edges);
% accumulate 3D grid of test data point errors
subs = [Xnum, Ynum];
idx_subs = any(subs==0,2);
sz = [length(lon),length(lat)];
alg_output.k_fold_delta_spatial = accumarray(subs(~idx_subs,:),...
    abs(alg_output.k_fold_delta(~idx_subs)),sz,@nanmean);
clear subs sz
% plot map
figure; hold on
worldmap([-90 90],[20 380]);
setm(gca,'mapprojection','robinson');
set(gcf,'units','inches','position',[0 5 20 10]);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(lat,[lon lon(end)+1],[alg_output.k_fold_delta_spatial ...
    alg_output.k_fold_delta_spatial(:,end)]');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
cmap = cmocean('amp'); cmap(1,:) = 1; colormap(cmap);
clim([0 (2/50)*(param_props.edges(end)-param_props.edges(1))]);
c=colorbar('location','southoutside');
c.Label.String = ['Average Absolute \Delta' param_props.label];
c.FontSize = 22;
c.TickLength = 0;
mlabel off; plabel off;
if ~isfolder([pwd '/' fig_dir]); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/' fig_name_2]);
% clean up
clear land cmap c
close

end

end

%% function to train ML models
function train_models(param_props,num_folds,...
    alg_dir,alg_fnames,variables,all_data,all_data_clusters,...
    train_idx,test_idx,data_per,alg_type,train_ratio,test_ratio,val_ratio,...
    numtrees,minLeafSize,numstumps,numbins,thresh,par_use,f,c,base_grid,...
    coverage)

%% define index for observations
if num_folds > 1
    obs_index_train = train_idx.(['f' num2str(f)]);
else
    obs_index_train = true(size(all_data.temperature));
end
% reduce training data volume if applicable:
% number of observations marked to use for training
num_obs = length(obs_index_train);
% unique random numbers (same for each cluster) equal to observational index
rng(f); numbers = randperm(num_obs);
% adjust observation index to fraction (i.e., 'data_per') its current sum
obs_index_train(numbers > (data_per.*num_obs)) = false;

% exclude Mediterranean cluster data before 2000
mean_lat = mean(all_data.latitude(all_data_clusters.(['c' num2str(c)]) > thresh));
mean_lon = mean(all_data.longitude(all_data_clusters.(['c' num2str(c)]) > thresh));
% if strcmp(base_grid,'RG')
if mean_lat > 32 && mean_lat < 41 && mean_lon > 10 && mean_lon < 30 % check for Med cluster
    in_med_new_idx = (all_data.longitude > -7 & all_data.longitude < 36) & ...
        (all_data.latitude > 30 & all_data.latitude < 45) & all_data.year > 2000;
    disp(['EXCLUDING YEAR FOR MED CLUSTER, #' num2str(c)])
    obs_index_train(all_data_clusters.(['c' num2str(c)]) > thresh & ...
        ~in_med_new_idx) = false;
    variables(strcmp(variables,'year')) = [];
    % elseif c == 5
    % figure; scatter(all_data.year(all_data_clusters.clusters == 5 & in_med_new_idx),...
    %      all_data.o2(all_data_clusters.clusters == 5 & in_med_new_idx));
    % 
    % figure; scatter(all_data.longitude(all_data_clusters.clusters == 9 & in_med_new_idx),...
    %      all_data.latitude(all_data_clusters.clusters == 9 & in_med_new_idx)); plot_land('xy');
    % 
    % figure; scatter(all_data.longitude(all_data_clusters.clusters == 5),...
    %      all_data.latitude(all_data_clusters.clusters == 5),5,...
    %      all_data.o2(all_data_clusters.clusters == 5));
    %     plot_land('xy'); colorbar;
end
% end

%% define observations index for k-fold testing
if num_folds > 1
    obs_index_test = test_idx.(['f' num2str(f)]);
    % reduce data volume if applicable:
    % number of observations marked to use for training
    num_obs = length(obs_index_test);
    % unique random numbers (same for each cluster) equal to observational index
    rng(f); numbers = randperm(num_obs)';
    % adjust observation index to fraction (i.e., 'data_per') its current sum
    obs_index_test(numbers > (data_per.*num_obs)) = false;
end

% check for data in cluster
if any(all_data_clusters.(['c' num2str(c)]) > thresh)

    % start timing fit
    tic
    
    %% fit model for each cluster
    if strcmp(alg_type,'FFNN')
        % define model parameters and train FFNN
        nodes1 = [10 15 20]; nodes2 = [20 15 10];
        alg = fit_FFNN(param_props.file_name,all_data,all_data_clusters.(['c' num2str(c)]),...
            obs_index_train,variables,nodes1,nodes2,train_ratio,val_ratio,test_ratio,...
            thresh,par_use);
    elseif strcmp(alg_type,'RFR')
        % define model parameters and train RFR
        NumPredictors = ceil(sqrt(length(variables)));
        alg = fit_RFR(param_props.file_name,all_data,all_data_clusters.(['c' num2str(c)]),...
            obs_index_train,variables,numtrees,minLeafSize,NumPredictors,0,thresh);
        alg = compact(alg); % convert RFR to compact
    elseif strcmp(alg_type,'GBM')
        % train GBM
        alg = fit_GBM(param_props.file_name,all_data,all_data_clusters.(['c' num2str(c)]),...
            obs_index_train,variables,numstumps,numbins,thresh);
    end
    
    %% stop timing fit
    if num_folds > 1
        fprintf(['Train ' alg_type ' for - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
    else
        fprintf(['Train ' alg_type ' for - Cluster #' num2str(c) ': ']);
    end
    
    toc

    %% for k-fold test
    if num_folds > 1

        % start timing predictions
        tic

        %% predict data for each cluster
        try
            if strcmp(alg_type,'FFNN')
                output = ...
                    run_FFNN(alg,all_data,all_data_clusters.(['c' num2str(c)]),...
                    obs_index_test,variables,thresh);
            elseif strcmp(alg_type,'RFR')
                output = ...
                    run_RFR(alg,all_data,all_data_clusters.(['c' num2str(c)]),...
                    obs_index_test,variables,thresh);
            elseif strcmp(alg_type,'GBM')
                output = ...
                    run_GBM(alg,all_data,all_data_clusters.(['c' num2str(c)]),...
                    obs_index_test,variables,thresh);
            end
            % estimate with ESPER neural network for oxygen
            if strcmp(param_props.file_name,'o2')
                idx_test = obs_index_test & ...
                    all_data_clusters.(['c' num2str(c)]) > thresh;
                ESPER_calc = ESPER_NN(7,[all_data.longitude(idx_test),...
                                           all_data.latitude(idx_test),...
                                           all_data.depth(idx_test)],...
                                           [all_data.salinity(idx_test),...
                                           all_data.temperature(idx_test)],...
                                           [1 2],'Equations',8);
                output_ESPER = nan(size(idx_test));
                output_ESPER(idx_test) = ESPER_calc.oxygen;
                output_ESPER = output_ESPER(obs_index_test);
            elseif strcmp(param_props.file_name,'no3')
                idx_test = obs_index_test & ...
                    all_data_clusters.(['c' num2str(c)]) > thresh;
                ESPER_calc = ESPER_NN(5,[all_data.longitude(idx_test),...
                                           all_data.latitude(idx_test),...
                                           all_data.depth(idx_test)],...
                                           [all_data.salinity(idx_test),...
                                           all_data.temperature(idx_test),...
                                           all_data.o2(idx_test)],...
                                           [1 2 6],'Equations',7);
                output_ESPER = nan(size(idx_test));
                output_ESPER(idx_test) = ESPER_calc.nitrate;
                output_ESPER = output_ESPER(obs_index_test);
            else
                output_ESPER = nan(sum(obs_index_test),1);
            end
        catch
            % if there aren't enough test data points, use NaNs
            output = nan(sum(obs_index_test),1);
            output_ESPER = nan(sum(obs_index_test),1);
        end

        % stop timing predictions
        fprintf(['Run ' alg_type ' for - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc

    end

% add NaNs to structures if no data in cluster 
else

    % print information
    if num_folds > 1
        fprintf(['Train ' alg_type ' for - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
        disp(' ');
    else
        fprintf(['Train ' alg_type ' for - Cluster #' num2str(c) ': N/A']);
        disp(' ');
    end

    % start timing predictions
    if num_folds > 1

        % add NaNs if no model
        % output = nan(sum(test_idx.(['f' num2str(f)])),num_clusters);

        % print information
        fprintf(['Run ' alg_type ' for - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
        disp(' ');

    end

end

% save test model and output for each cluster
if ~isfolder([pwd '/' alg_dir]); mkdir(alg_dir); end
if any(all_data_clusters.(['c' num2str(c)]) > thresh)
    if num_folds > 1
        parsave([alg_dir '/' alg_fnames{f,c} '.mat'],alg,alg_type,output,'output',output_ESPER,'output_ESPER',obs_index_test,'obs_index_test');
    else
        % genFunction(FFNN.n1,['functions/' 'FFNN1_o2_C13_rg.m'])
        parsave([alg_dir '/' alg_fnames{f,c} coverage '.mat'],alg,alg_type);
    end
else
    if num_folds > 1
        parsave([alg_dir '/' alg_fnames{f,c} '.mat'],output,'output',output_ESPER,'output_ESPER',obs_index_test,'obs_index_test');
    end
end

end

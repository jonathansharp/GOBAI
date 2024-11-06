% train_gbm
%
% DESCRIPTION:
% This function trains gradient boosted machines in GMM clusters for GOBAI.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/5/2024

function train_gbm(param,base_grid,file_date,...
    float_file_ext,num_clusters,variables,...
    numstumps,numbins,thresh,numWorkers_train,varargin)

%% set defaults and process optional input arguments
num_folds = 1;
glodap_only = 0;
data_per = 1;
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'num_folds')
        num_folds = varargin{i+1};
    elseif strcmpi(varargin{i}, 'glodap_only')
        glodap_only = varargin{i+1};
    elseif strcmpi(varargin{i}, 'reduce_data')
        data_per = varargin{i+1};
    end
end

%% directory base
dir_base = create_dir_base('GBM',{base_grid;num_clusters;file_date;...
    float_file_ext;numstumps;numbins});

%% process parameter name
[param1,param2,param3,~,~,~,~,param_edges] = param_name(param);

%% load data
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    load([param1 '/Data/processed_all_' param '_data_' file_date float_file_ext '.mat'],'all_data');
else
    load([param1 '/Data/' base_grid '_' param '_data_' file_date float_file_ext '.mat'],'all_data');
end

%% load data clusters
load([param1 '/Data/all_data_clusters_' base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');

%% load data cluster indices (for k-fold testing)
if num_folds > 1
    load([param1 '/Data/k_fold_data_indices_'  base_grid '_' num2str(num_clusters) ...
        '_' num2str(num_folds) '_' file_date float_file_ext '.mat'],...
        'num_folds','train_idx','test_idx');
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
gbm_dir = [param1 '/Models/' dir_base];
gbm_fnames = cell(num_folds,num_clusters,1);
for f = 1:num_folds
    for c = 1:num_clusters
        if num_folds > 1
            gbm_fnames(c) = {['GBM_' param2 '_C' num2str(c) '_F' num2str(f) '_test']};
        else
            gbm_fnames(c) = {['GBM_' param2 '_C' num2str(c)]};
        end
    end
end

%% variables for k-fold test
if num_folds > 1
    kfold_dir = [param1 '/KFold/GBM/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
    kfold_name = ['GBM_output_tr' num2str(numstumps) '_bin' num2str(numbins)];
    fig_dir = [param1 '/Figures/KFold/GBM/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
    fig_name_1 = ['k_fold_comparison_tr' num2str(numstumps) '_bin' num2str(numbins) '.png'];
    fig_name_2 = ['k_fold_spatial_comparison_tr' num2str(numstumps) '_bin' num2str(numbins) '.png'];
end

%% fit GBMs using all data

% start timing training
tStart = tic;

% set up parallel pool
% tic; parpool(num_clusters); fprintf('Pool initiation:'); toc;

% define model parameters
NumPredictors = ceil(sqrt(length(variables)));

% fit models for each fold (if applicable)
for f = 1:num_folds

    % fit models for each cluster
    for c = 1:num_clusters
    
        % check for data in cluster
        if any(all_data_clusters.clusters == c)
        
            % start timing fit
            tic
            
            %% define index for observations
            if num_folds > 1
                obs_index = train_idx.(['f' num2str(f)]);
            else
                obs_index = true(size(all_data.temperature));
            end
            
            %% reduce data volume if applicable
            % number of observations marked to use for training
            num_obs = length(obs_index);
            % unique random numbers equal to observational index
            numbers = randperm(num_obs);
            % adjust observation index to fraction (i.e., 'data_per') its current sum
            obs_index(numbers > (data_per.*num_obs)) = false;
            
            %% fit model for each cluster
            GBM = ...
                fit_GBM(param2,all_data,all_data_clusters.(['c' num2str(c)]),...
                obs_index,variables,numstumps,numbins,thresh);
            
            %% stop timing fit
            if num_folds > 1
                fprintf(['Train GBM for ' base_grid ' - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
            else
                fprintf(['Train GBM for ' base_grid ' - Cluster #' num2str(c) ': ']);
            end

            toc

            %% for k-fold test
            if num_folds > 1

                % start timing predictions
                tic

                %% define index for observations
                obs_index = test_idx.(['f' num2str(f)]);

                %% reduce data volume if applicable
                num_obs = sum(obs_index);
                numbers = randperm(length(obs_index));
                obs_index = numbers < (data_per.*num_obs);

                %% predict data for each cluster
                output = ...
                    run_GBM(GBM,all_data,all_data_clusters.(['c' num2str(c)]),...
                    obs_index,variables,thresh);
        
                % stop timing predictions
                fprintf(['Run GBM for ' base_grid ' - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
                toc

            end

        % add NaNs to structures if no data in cluster 
        else

            % print information
            if num_folds > 1
                fprintf(['Train GBM for ' base_grid ' - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
                disp
            else
                fprintf(['Train GBM for ' base_grid ' - Cluster #' num2str(c) ': N/A']);
                disp('');
            end

            % start timing predictions
            if num_folds > 1

                % add NaNs if no model
                output = nan(sum(test_idx.(['f' num2str(f)])),num_clusters);

                % print information
                fprintf(['Run GBM for ' base_grid ' - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
                disp('');

            end

        end

        % save test model and output for each cluster
        if ~isfolder([pwd '/' gbm_dir]); mkdir(gbm_dir); end
        if any(all_data_clusters.clusters == c)
            if num_folds > 1
                parsave([gbm_dir '/' gbm_fnames{f,c}],GBM,'GBM',output,'output');
            else
                parsave([gbm_dir '/' gbm_fnames{f,c}],GBM,'GBM');
            end
        else
            if num_folds > 1
                parsave([gbm_dir '/' gbm_fnames{f,c}],output,'output');
            end
        end

    end

end

% end parallel session
delete(gcp('nocreate'));

% stop timing full script
if num_folds > 1
    fprintf(['GBM k-Fold Training for ' base_grid ': ']);
else
    fprintf(['GBM Training for ' base_grid ': ']);
end

% print elapsed time in minutes
tElapsed = toc(tStart);
disp(['Elapsed time is ' num2str(min(tElapsed)) ' minutes.'])

%% calculate statistics for k-fold test
if num_folds > 1

% calculate weighted average over each cluster using probabilities
for f = 1:num_folds
    % pre-allocate probabilities
    probs_matrix = [];
    for c = 1:num_clusters
        % assemble matrix of probabilities greater than the threshold (5%)
        probs_array = all_data_clusters.(['c' num2str(c)])(test_idx.(['f' num2str(f)]));
        probs_array(probs_array < thresh) = NaN;
        probs_matrix = [probs_matrix,probs_array];
        clear probs_array
        % load output
        load([gbm_dir '/' gbm_fnames{f,c}],'output')
        gbm_output.(['f' num2str(f)])(:,c) = output;
        clear output
    end
    gbm_output.(['f' num2str(f) '_mean']) = ...
        double(sum(gbm_output.(['f' num2str(f)]).*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
end
% clean up
clear f c
% aggregate output from all folds
gbm_output.(['k_fold_test_' param2]) = nan(size(all_data.(param2)));
for f = 1:num_folds
    gbm_output.(['k_fold_test_' param2])(test_idx.(['f' num2str(f)])) = ...
        gbm_output.(['f' num2str(f) '_mean']);
end
% compare k-fold output to data
gbm_output.k_fold_delta = gbm_output.(['k_fold_test_' param2]) - all_data.(param2);
% calculate error stats
gbm_mean_err = mean(gbm_output.k_fold_delta);
gbm_med_err = median(gbm_output.k_fold_delta);
gbm_rmse = sqrt(mean(gbm_output.k_fold_delta.^2));
% save predicted data
if ~isfolder([pwd '/' kfold_dir]); mkdir(kfold_dir); end
save([kfold_dir '/' kfold_name],'gbm_output','gbm_rmse',...
    'gbm_med_err','gbm_mean_err','-v7.3');
clear gbm_output gbm_rmse gbm_med_err gbm_mean_err

%% plot histogram of errors
load([kfold_dir '/' kfold_name],'gbm_output','gbm_rmse');
figure('visible','off'); hold on;
set(gca,'fontsize',12);
set(gcf,'position',[100 100 600 400]);
[counts,bin_centers] = hist3([all_data.(param2) gbm_output.(['k_fold_test_' param2])],...
    'Edges',{param_edges param_edges});
h=pcolor(bin_centers{1},bin_centers{2},counts');
plot([param_edges(1) param_edges(end)],[param_edges(1) param_edges(end)],'k--');
set(h,'EdgeColor','none');
xlim([param_edges(1) param_edges(end)]); ylim([param_edges(1) param_edges(end)]);
xlabel(['Measured ' param2 ' (\mumol kg^{-1})']);
ylabel(['GBM ' param2 ' (\mumol kg^{-1})']);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
caxis([1e0 1e5]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text((3/5)*param_edges(end),(1/10)*param_edges(end),...
    ['RMSE = ' num2str(round(gbm_rmse,1)) '\mumol kg^{-1}'],'fontsize',12);
if ~isfolder([pwd '/' fig_dir]); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/' fig_name_1]);
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
gbm_output.k_fold_delta_spatial = accumarray(subs(~idx_subs,:),...
    abs(gbm_output.k_fold_delta(~idx_subs)),sz,@nanmean);
clear subs sz
% plot map
figure; hold on
worldmap([-90 90],[20 380]);
setm(gca,'mapprojection','robinson');
set(gcf,'units','inches','position',[0 5 20 10]);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(lat,[lon lon(end)+1],[gbm_output.k_fold_delta_spatial ...
    gbm_output.k_fold_delta_spatial(:,end)]');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
cmap = cmocean('amp'); cmap(1,:) = 1; colormap(cmap);
caxis([0 (2/50)*param_edges(end)]);
c=colorbar('location','southoutside');
c.Label.String = ['Average Absolute \Delta' param3];
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

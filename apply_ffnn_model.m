%% embedded function for processing 3D grids and applying FFNN models
function apply_ffnn_model(TS,num_clusters,ffnn_dir,ffnn_fnames,...
    base_grid,m,w,xdim,ydim,zdim,variables_TS,thresh,gobai_ffnn_dir)

    % convert to arrays
    TS_index = ~isnan(TS.temperature);
    vars = fieldnames(TS);
    for v = 1:length(vars)
        if ndims(TS.(vars{v})) == 3
            TS.([vars{v} '_array']) = TS.(vars{v})(TS_index);
            TS = rmfield(TS,vars{v});
        end
    end

    % replicate time variables as arrays
    vars = fieldnames(TS);
    for v = 1:length(vars)
        if length(TS.(vars{v})) == 1
            TS.([vars{v} '_array']) = repmat(TS.(vars{v}),size(TS.temperature_array));
            TS = rmfield(TS,vars{v});
        end
    end

    % calculate absolute salinity, conservative temperature, potential density
    TS.salinity_abs_array = gsw_SA_from_SP(TS.salinity_array,TS.pressure_array,...
        TS.longitude_array,TS.latitude_array);
    TS.temperature_cns_array = gsw_CT_from_t(TS.salinity_abs_array,...
        TS.temperature_array,TS.pressure_array);
    TS.sigma_array = gsw_sigma0(TS.salinity_abs_array,TS.temperature_cns_array);

    % pre-allocate
    gobai_matrix = single(nan(length(TS.temperature_array),num_clusters));
    probs_matrix = single(nan(length(TS.temperature_array),num_clusters));

    % apply models for each cluster
    for c = 1:num_clusters

      % check for existence of model
      if isfile([ffnn_dir '/' ffnn_fnames{c} '.mat'])

        % load GMM cluster probabilities for this cluster and month, and convert to array
        GMM_probs = ...
            load(['Data/GMM_' base_grid '_' num2str(num_clusters) '/c' num2str(c) ...
            '/m' num2str(m) '_w' num2str(w)],'GMM_cluster_probs');
        GMM_probs.probabilities_array = GMM_probs.GMM_cluster_probs(TS_index); % convert to array
        GMM_probs.probabilities_array(GMM_probs.probabilities_array < thresh) = NaN; % remove probabilities below thresh
        GMM_probs = rmfield(GMM_probs,{'GMM_cluster_probs'}); % remove 3D matrix
        probs_matrix(:,c) = GMM_probs.probabilities_array; % add to probability matrix

        % load model for this cluster
        alg = load([ffnn_dir '/' ffnn_fnames{c}],'FFNN');
    
        % predict data for each cluster
        gobai_matrix(:,c) = ...
            run_FFNN(alg.FFNN,TS,GMM_probs.probabilities_array,...
            true(size(TS.temperature_array)),variables_TS,thresh);

      end

    end

    % calculate weighted average over each cluster using probabilities
    gobai_array = ...
        double(sum(gobai_matrix.*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
    
    % convert back to 3D grid
    gobai_3d = nan(xdim,ydim,zdim);
    gobai_3d(TS_index) = gobai_array;
    
    % save monthly output
    if ~isfolder([pwd '/' gobai_ffnn_dir]); mkdir(gobai_ffnn_dir); end
%     nccreate([gobai_ffnn_dir 'm' num2str(m) '_w' num2str(w) '.nc'],...
%         gobai_3d,'gobai',w,'w',m,'m');
%     
    parsave([gobai_ffnn_dir 'm' num2str(m) '_w' num2str(w)],...
        gobai_3d,'gobai',w,'w',m,'m');

end

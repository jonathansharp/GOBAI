% combine_gobai
%
% DESCRIPTION:
% This function combines output from gobai gridded
% fields obtained via different machine learning algorithms.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/22/2023

%% create directory names
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
gobai_ffnn_dir = ... % FFNN
    ['Data/GOBAI/FFNN_c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*val_ratio) '/'];
gobai_rfr_dir = ... % RFR
    ['Data/GOBAI/RFR_c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numtrees) '_lf' num2str(minLeafSize) '/'];
gobai_gbm_dir = ... % GBM
    ['Data/GOBAI/GBM_c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numstumps) '/'];
gobai_dir = ... % final product
    ['Data/GOBAI/' num2str(num_clusters) '_' file_date float_file_ext '/'];

%% process time
start_year = 2004;
end_year = floor(snap_date/1e2);
end_month = mod(snap_date,1e2);
years = repelem(start_year:end_year,12)';
months = repmat(1:12,1,length(years)/12)';
years = years(1:end-(12-end_month));
months = months(1:end-(12-end_month));
clear start_year end_year end_month

%% average over each time window
parfor m = 1:length(years)
    
    % load monthly output (FFNN)
    gobai_3d_ffnn = load([gobai_ffnn_dir 'm' num2str(m)],'x');
    % load monthly output (RFR)
    gobai_3d_rfr = load([gobai_rfr_dir 'm' num2str(m)],'x');
    % load monthly output (GBM)
    gobai_3d_gbm = load([gobai_gbm_dir 'm' num2str(m)],'x');

    % average monthly outputs
    gobai_3d_ens = mean(cat(4,gobai_3d_ffnn,gobai_3d_rfr,gobai_3d_gbm));

    % save each timestep
    parsave_v1([gobai_dir 'm' num2str(m)],gobai_3d_ens);

end

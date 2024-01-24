% Runs all scripts to train and apply GOBAI-O2 GBM algorithms

addpath(genpath(pwd));
try
    train_gbm;
    mod_type = 'GBM';
    %plot_gobai_mean;
    plot_gobai_animation;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
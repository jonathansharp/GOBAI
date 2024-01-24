% Runs all scripts to train and apply GOBAI-O2 FFNN algorithms

addpath(genpath(pwd));
try
    train_ffnn;
    mod_type = 'FFNN';
    %plot_gobai_mean;
    plot_gobai_animation;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
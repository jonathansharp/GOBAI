% Runs all scripts to train and apply GOBAI-O2 RFR algorithms

addpath(genpath(pwd));
try
    train_rfr;
    mod_type = 'RFR';
    %plot_gobai_mean;
    plot_gobai_animation;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
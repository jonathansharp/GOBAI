% Runs all scripts to combine GOBAI-O2 FFNN gridded fields

try
    combine_gobai;
    mod_type = 'ENS';
    %plot_gobai_mean;
    plot_gobai_animation;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
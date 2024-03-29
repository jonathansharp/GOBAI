% Runs all scripts to create GOBAI-NO3

addpath(genpath(pwd));
try
    create_gobai_no3;
    plot_gobai_no3;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
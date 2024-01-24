% Runs all scripts to create GOBAI-O2

addpath(genpath(pwd));
try
    create_gobai_o2;
    plot_gobai_o2;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
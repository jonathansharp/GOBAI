% Runs all scripts to create GOBAI-O2

addpath(genpath(pwd));
try
    fit_full_models;
    create_gobai_o2;
    plot_gobai_o2;
    disp('success!');
catch ME
    disp(ME);
    disp(['Function: ' ME.stack(1).name]);
    disp(['Line: ' num2str(ME.stack(1).line)]);
end
exit
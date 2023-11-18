% Runs all scripts to train GOBAI-O2 FFNN algorithms

addpath(genpath(pwd));
try
    train_ffnn;
    disp('success!');
catch ME
    disp(ME);
    disp(['Function: ' ME.stack(1).name]);
    disp(['Line: ' num2str(ME.stack(1).line)]);
end
exit
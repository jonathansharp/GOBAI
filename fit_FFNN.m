function FFNN = fit_FFNN(target,data,probabilities,index,...
    variables,nodes1,nodes2,train_ratio,val_ratio,test_ratio,thresh,par_use)
    
% index for training based on input index and cluster probabilities
idx_train = index & probabilities > thresh;

% assemble matrix of predictors
predictors = nan(sum(idx_train),length(variables));
for v = 1:length(variables)
    predictors(:,v) = data.(variables{v})(idx_train);
end

% determine number of networks
numNets = length(nodes1);

% fit each network
for n = 1:numNets
    % training function
    % trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation
    trainFcn = 'trainbr';  % Baysian Regularization    
    % create a Fitting Network
    net = feedforwardnet([nodes1(n) nodes2(n)],trainFcn);
    % Set training parameter criteria
    net.trainParam.max_fail = 6; % default: 6
    net.trainParam.mu_max = 1e+10; % default: 1e+10
    net.trainParam.min_grad = 1e-7; % default: 1e-7
    net.trainParam.epochs = 1000; % default: 1000
    net.trainParam.showWindow = false;
    % Set performance parameters
    net.performParam.regularization = 0;
    % Setup Division of Data for Training, Validation, Testing
    net.divideParam.trainRatio = train_ratio;
    net.divideParam.valRatio = val_ratio;
    net.divideParam.testRatio = test_ratio;
    % Train the Network
    FFNN.(['n' num2str(n)]) = ...
        train(net,predictors',data.(target)(idx_train)','useParallel',par_use);
    % weird error was fixed by removing cpt_map toolbox from path 
end

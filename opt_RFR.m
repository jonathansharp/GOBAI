function results = opt_RFR(target,data,probabilities,index,variables,thresh)
    
% index for training based on input index and cluster probabilities
idx_train = index & probabilities > thresh;

% assemble matrix of predictors
predictors = nan(sum(idx_train),length(variables));
for v = 1:length(variables)
    predictors(:,v) = data.(variables{v})(idx_train);
end

% define hyperparameters
maxMinLS = 20;
minLS = optimizableVariable('minLS',[1,maxMinLS],'Type','integer');
numPTS = optimizableVariable('numPTS',[1,size(X,2)-1],'Type','integer');
hyperparametersRF = [minLS; numPTS];

% define objective function
function oobErr = oobErrRF(params,predictors,data,target,idx_train,numtrees)
%oobErrRF Trains random forest and estimates out-of-bag quantile error
%   oobErr trains a random forest of 300 regression trees using the
%   predictor data in X and the parameter specification in params, and then
%   returns the out-of-bag quantile error based on the median. X is a table
%   and params is an array of OptimizableVariable objects corresponding to
%   the minimum leaf size and number of predictors to sample at each node.
randomForest = TreeBagger(numtrees,predictors,data.(target)(idx_train),...
    'Method','regression','OOBPrediction','on','MinLeafSize',params.minLS,...
    'NumPredictorstoSample',params.numPTS);
oobErr = oobQuantileError(randomForest);
end

% optimize objective function
results = bayesopt(@(params)oobErrRF(params,predictors,data,target,idx_train,numtrees),...
    hyperparametersRF,'AcquisitionFunctionName','expected-improvement-plus','Verbose',0);

end

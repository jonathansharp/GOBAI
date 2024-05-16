function GBM = fit_GBM(target,data,probabilities,index,...
    variables,numstumps,numbins,thresh)
    
% index for training based on input index and cluster probabilities
idx_train = index & probabilities > thresh;

% assemble matrix of predictors
predictors = nan(sum(idx_train),length(variables));
for v = 1:length(variables)
    predictors(:,v) = data.(variables{v})(idx_train);
end

lrn_rate = 0.55; % optimized 

% perform random forest regression fit
GBM = fitrensemble(predictors,data.(target)(idx_train),...
    'NumLearningCycles',numstumps,'NumBins',numbins,'LearnRate',lrn_rate);

GBM = fitrensemble(predictors,data.(target)(idx_train),...
    'CrossVal','on');

GBM = fitrensemble(predictors,data.(target)(idx_train),...
    'Method','LSBoost','NumLearningCycles',50,'OptimizeHyperparameters','LearnRate');


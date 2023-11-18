function RFR = fit_RFR(target,data,probabilities,index,...
    variables,numtrees,minLeafSize,NumPredictors,oob,thresh)
    
% index for training based on input index and cluster probabilities
idx_train = index & probabilities > thresh;

% define whether to obtain predictor importance
if oob == 1; oob_info = 'on'; else oob_info = 'off'; end

% assemble matrix of predictors
predictors = nan(sum(idx_train),length(variables));
for v = 1:length(variables)
    predictors(:,v) = data.(variables{v})(idx_train);
end

% use parallel computing
paroptions = statset('UseParallel',true);

% perform random forest regression fit
RFR = ...
    TreeBagger(numtrees,predictors,data.(target)(idx_train),...
    'Method','r','Options',paroptions,...
    'OOBPrediction',oob_info,...
    'OOBPredictorImportance',oob_info,...
    'MinLeafSize',minLeafSize,...
    'NumPredictorsToSample',NumPredictors);

if oob == 1
    % Plot Out-of-Bag RMSE based on tree number
    figure;
    plot(sqrt(oobError(RFR)),'k','LineWidth',2);
    xlabel('Number of Grown Trees');
    ylabel('Out-of-Bag Root Mean Squared Error');
    close
    % Plot importance of each predictor
    figure;
    bar(model_oxy_RF.OOBPermutedPredictorDeltaError);
    xlabel('Predictor') ;
    ylabel('Out-of-Bag Feature Importance');
    xticklabels(OXY.train.Properties.VariableNames(all_oxy_RF(1:end-1)));
    exportgraphics(gcf,['Figures/RFR_importance_' type '_' basins{b} '.png']);
    close
end

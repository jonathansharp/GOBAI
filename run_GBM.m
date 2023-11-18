function  output = run_GBM(GBM,data,probabilities,index,variables,thresh)

% index for prediction based on input index and cluster probabilities
idx_test = index & probabilities > thresh;
test_sum = sum(idx_test);

% assemble matrix of predictors
predictors = nan(test_sum,length(variables));
for v = 1:length(variables)
    predictors(:,v) = data.(variables{v})(idx_test);
end

% compute predictions from GBM
output_temp = predict(GBM,predictors);

% expand output to match length of input data
output = nan(size(idx_test)); % pre-allocate long output array
output(idx_test) = output_temp; % add data to long output array
output = output(index); % condense long output array to match input index

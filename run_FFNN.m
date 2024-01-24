function  output = run_FFNN(FFNN,data,probabilities,index,variables,thresh)

% index for prediction based on input index anzid cluster probabilities
idx_test = index & probabilities > thresh;
test_sum = sum(idx_test);

% assemble matrix of predictors
predictors = nan(test_sum,length(variables));
for v = 1:length(variables)
    predictors(:,v) = data.(variables{v})(idx_test);
end

% determine number of networks
numNets = length(fieldnames(FFNN));

% make predictions using each network
output_temp = nan(test_sum,numNets);
for n = 1:numNets
    output_temp(:,n) = FFNN.(['n' num2str(n)])(predictors')';
end

% average across each network
output_mean = mean(output_temp,2);

% expand output to match length of input data
output = nan(size(idx_test)); % pre-allocate long output array
output(idx_test) = output_mean; % add data to long output array
output = output(index); % condense long output array to match input index

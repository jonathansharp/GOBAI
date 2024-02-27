function  output = run_RFR(RFR,data,probabilities,index,variables,thresh)

% index for prediction based on input index and cluster probabilities
idx_test = index & probabilities > thresh;
test_sum = sum(idx_test);

% assemble matrix of predictors
predictors = nan(test_sum,length(variables));
for v = 1:length(variables)
    predictors(:,v) = data.(variables{v})(idx_test);
end

% define size and number of chunks
chunkSize = 100000;
chunkNums = ceil(test_sum/chunkSize);
% pre-allocate output
output_temp = nan(chunkSize,chunkNums);
% pad predictors to be divisible by chunk size
predictors = [predictors;nan(chunkSize*chunkNums-test_sum,length(variables))];
%predictors = gpuArray(predictors);

% compute predictions, 100,000 data points at a time
for chunk = 1:chunkNums
    idx1 = (chunk-1)*chunkSize + 1;
    idx2 = (chunk-1)*chunkSize + chunkSize;
    output_temp(:,chunk) = predict(RFR,predictors(idx1:idx2,:));
end

% unwrap output and remove predictions from padded NaNs
output_temp = output_temp(:);
output_temp = output_temp(1:test_sum);

% expand output to match length of input data
output = nan(size(idx_test)); % pre-allocate long output array
output(idx_test) = output_temp; % add data to long output array
output = output(index); % condense long output array to match input index

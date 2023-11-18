function SVM = = opt_SVM(target,data,probabilities,index,...
    variables,thresh)
    
% index for training based on input index and cluster probabilities
idx_train = index & probabilities > thresh;

% assemble matrix of predictors
predictors = nan(sum(idx_train),length(variables));
for v = 1:length(variables)
    predictors(:,v) = data.(variables{v})(idx_train);
end

% perform random forest regression fit
lambda = logspace(-10,-1,10);
SVM = fitrlinear(predictors',data.(target)(idx_train)',...
    'ObservationsIn','columns','KFold',5,'Lambda',lambda,...
    'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
mse = kfoldLoss(SVM);
SVM = fitrlinear(predictors',data.(target)(idx_train)',...
    'ObservationsIn','columns','Lambda',lambda,...
    'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
numNZCoeff = sum(SVM.Beta~=0);

figure;
[h,hL1,hL2] = plotyy(log10(lambda),log10(mse),log10(lambda),log10(numNZCoeff));
hL1.Marker = 'o';
hL2.Marker = 'o';
ylabel(h(1),'log_{10} MSE');
ylabel(h(2),'log_{10} non-zero coeff freq');
xlabel('log_{10} lambda');
hold off
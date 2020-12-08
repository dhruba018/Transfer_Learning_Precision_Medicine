function [predY, nnWeights, nnIdx] = knnregress(trainX, trainY, testX, varargin)
% KNNREGRESS performs k-nearest neighbors (kNN) regression on given data.
    % NOTATION: function [predY, nnWeights] = KNNREGRESS(trainX, trainY, testX, K, dType)
    %   Inputs:
    %       trainX, trainY => Training data
    %       testX => Test data covariates
    %       K (optional) => Numbers of NNs to be used. Default = 10
    %       dType (optional) => Distance type to be used in determining NNs. Default = euclidean.
    %   Outputs:
    %       predY => Predicted kNN response
    %       nnWeights => kNN weights corresponding to each sample in testX. Size = nTest x K.
    %       nnIdx => Indices of K NNs for each sample in testX. Size = nTest x K.
    %
    % (c) 2018 S.R.Dhruba
    
    % Check optional parametes...
    if nargin < 3
        error('Insufficient inputs!!!')
    end
    K = 10 * (nargin < 4) + varargin{1} * (nargin >= 4);                                        % #neighbors 
    dType = 'euclidean';    if (nargin == 5),     dType = varargin{2};      end       % Distance type
    
    % Check input sizes...
    if size(trainX, 1) ~= size(trainY, 1)
        error('Training datasets must have the same size!!!')
    elseif size(trainX, 2) ~= size(testX, 2)
        error('Test data must have the same number of features as Training data!!!')
    elseif size(trainY, 2) > 1
        error('Only univariate regression is available... trainY must have a single column!!!')
    elseif K > size(trainX, 2)
        error('K must be less or equal than the number of features!!!')
    end
    
    % Form data matrices...
    nTrain = size(trainX, 1);       nTest = size(testX, 1);
    X = [trainX; testX];                                          % Complete covariate matrix
    dX = squareform(pdist(X, dType)');                 % Distance matrix
    
    % kNN responses for new data...
    dXtest = dX(nTrain+1:end, 1:nTrain);               % dist(testX, trainX)
    [~, nnIdx] = sort(dXtest, 2, 'ascend');             % Ordered neighboring samples in trainX
    nnWeights = zeros(nTest, K);                           % kNN weights
    predY = zeros(nTest, 1);                                  % kNN prediction
    for k = 1 : nTest
        dist_knn = dXtest(k, nnIdx(k, 1:K));             % Dist of testX(k, :) from K NNs in trainX
        nnWeights(k, :) = (1 ./ dist_knn) / sum((1 ./ dist_knn));       % kNN weights
        predY(k) = nnWeights(k, :) * trainY(nnIdx(k, 1:K));             % kNN fit
    end
    
end

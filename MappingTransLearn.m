function pPred = MappingTransLearn(pTrain, sTrain, pTest, sTest)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%      pPred = MappingTransLearn(pTrain, sTrain, pTest, sTest)
% Nonlinear Mapping based Transfer Learning Prediction models, where the
% input & output data in primary are mapped to the secondary space using
% polynomial regression mapping.
%
% INPUT ARGUMENTS:
%      pTrain:  n1 x (p+1) matrix of primary training data. First p columns
%                   contain training input Xp_train. Last column is the
%                   training response yp_train.
%      sTrain:  n1 x (p+1) matrix of secondary training data. First p
%                   columns contain training input Xs_train. Last column
%                   is the training response ys_train.
%      pTest:   n2 x p matrix of primary test data. Contains test input
%                   Xp_test.
%      sTest:   n2 x (p+1) matrix of secondary test data. First p columns
%                  contain test input Xs_test. Last column is the test
%                  response ys_test.
% 
% OUTPUT ARGUMENTS:
%      pPred:  n2 x 1 matrix of primary prediction output.
%
% (c) 2017 S. R. Dhruba
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pPred = MappingTransLearn(pTrain, sTrain, pTest, sTest)
% Inputs...
G11 = pTrain(:, 1:end-1);           G12 = pTest;                         % Primary X
G21 = sTrain(:, 1:end-1);           G22 = sTest(:, 1:end-1);         % Secondary X
d11 = pTrain(:, end);                                                             % Primary y
d21 = sTrain(:, end);                d22 = sTest(:, end);               % Secondary y
G1 = [G11; G12];                       G2 = [G21; G22];            d2 = [d21; d22];
[n, p] = size(G1);
nTrain = size(G11, 1);               nTest = size(G12, 1);


% % Input Mapping...
%       x2 = h(x1) = a0 + a1*x1 = X1*A
fprintf('>> Retrieving input mapping... ')
nX = 1;
H12 = zeros(nX+1, p);                                           % Mapping
for j = 1:p
    G11A = ones(nTrain, nX+1);                               % G1 = [1, g1]
    for k = 1:nX
        G11A(:, k+1) = G11(:, j).^k;
    end
    H12(:, j) = pinv(G11A) * G21(:, j);                     % A = pinv(G1)*g2
end
fprintf('Finished! \n')

% Map to secondary...
fprintf('>> Mapping test input... ')
G22m = zeros(nTest, p);                                       % Mapped gx
for j = 1:p
    G12A = ones(nTest, nX+1);
    for k = 1:nX
        G12A(:, k+1) = G12(:, j).^k;
    end
    G22m(:, j) = G12A * H12(:, j);                          % g2 = G1*A
end
fprintf('Finished! \n')


% % Response Mapping...
%       y1 = f(y2) = b0 + b1*y2 + b2*y2^2 = Y2*B
fprintf('>> Retrieving response mapping... ')
nY = 2;
D21 = ones(nTrain, nY+1);                                     % D2 = [1, d2, d2^2]
D22 = ones(nTest, nY+1);
for k = 1:nY
    D21(:, k+1) = d21.^k;
    D22(:, k+1) = d22.^k;
end
F21 = pinv(D21) * d11;                                          % B = pinv(D2)*d1
% d12m = D22*F21;                                                 % y1 = Y2*B
fprintf('Finished! \n')


% % BIAS CORRECTION...
rng(0)                                                                  % Seed
% Validation set...
fprintf('>> Retrieving Bias Model:Started... \n')
nVal = round(nTrain/2);
idx1 = sort(randperm(nTrain, nVal))';
idx2 = (1:nTrain)';             idx2(idx1) = [ ];
val_G11 = G11(idx1, :);        val_G12 = G11(idx2, :);
val_G21 = G21(idx1, :);       val_G22 = G21(idx2, :);
val_d11 = d11(idx1);            val_d12 = d11(idx2, :);
val_d21 = d21(idx1);           val_d22 = d21(idx2, :);

fprintf('>> \t Inner MP loop... \n')
% MP inner recursion for Bias model...
% Input map...
h12 = zeros(nX+1, p);
for j = 1:p
    val_G11A = ones(size(val_G11, 1), nX+1);           % G1 = [1, g1]
    for k = 1:nX
        val_G11A(:, k+1) = val_G11(:, j).^k;
    end
    h12(:, j) = pinv(val_G11A) * val_G21(:, j);          % A = pinv(G1)*g2
end

val_G21m = zeros(size(val_G12));
for j = 1:p
    val_G12A = ones(size(val_G12, 1), nX+1);
    for k = 1:nX
        val_G12A(:, k+1) = val_G12(:, j).^k;
    end
    val_G21m(:, j) = val_G12A * h12(:, j);               % g2 = G1*A
end

% Response map...
val_D21 = ones(size(val_d21, 1), nY+1);                % D2 = [1, d2, d2^2]
for k = 1:nY
    val_D21(:, k+1) = val_d21.^k;
end
f21 = pinv(val_D21) * val_d11;                              % B = pinv(D2)*d1


% % Prediction model training w/ secondary...
fprintf('>> \t Training Prediction Model... ')
n_tree = 100;       Train_X2 = G2;      Train_y2 = d2;
M = TreeBagger(n_tree, Train_X2, Train_y2, 'method', 'regression');
fprintf('Finished! \n')


% Prediction for validation set...
Valid_X21 = val_G21m;
Pred_y21 = predict(M, Valid_X21);
val_d22m = Pred_y21;

% Map to original...
val_D22m = ones(size(val_d22m, 1), nY+1);           % D2 = [1, d2, d2^2]
for k = 1:nY
    val_D22m(:, k+1) = val_d22m.^k;
end
val_d12m = val_D22m*f21;                                   % d1 = D2*B

% Learning Bias Model...
Train_X_bias = [val_G12, val_d12];     Train_y_bias = val_d12 - val_d12m;
M_bias = TreeBagger(n_tree, Train_X_bias, Train_y_bias, 'method', 'regression');
fprintf('>> Finished! \n')


% Mapped Prediction...
fprintf('>> Performing MP... \n')
Test_X21 = G22m;
Pred_y21 = predict(M, Test_X21);
d22m = Pred_y21;

% Map to Original Space...
fprintf('>> \t Mapping back to Primary Domain! \n')
D22m = ones(size(d22m, 1), nY+1);                        % D2 = [1, d2, d2^2]
for k = 1:nY
    D22m(:, k+1) = d22m.^k;
end
d12m = D22m*F21;                                               % d1 = D2*B

fprintf('>> \t Estimating Bias! \n')
% Estimating Bias...
Test_X_bias = [G12, d12m];
Pred_y_bias = predict(M_bias, Test_X_bias);
d12b = Pred_y_bias;                                            % Bias Estimate

% Bias corrected prediction...
fprintf('>> \t Performing Bias Corrected Prediction \n')
d12mb = d12m + d12b;
pPred = d12mb;
fprintf('>> Finished! \n')




function pPred = LatentPredTransLearn(pTrain, sTrain, pTest, sTest)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%      pPred = LatentPredTransLearn(pTrain, sTrain, pTest, sTest)
% Latent Variable Cost Optimization based Transfer Learning Prediction
% models. Includes three prediction models --
%      Latent Regression Prediction (LRP) 
%      Latent-Latent Prediction (LLP) 
%      Combined Latent Prediction (CLP) 
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
%      pPred:  n2 x 3 matrix of primary prediction output. First column
%                  contains the LRP prediction yp_pred_LRP, second column
%                  contains the LLP prediction, yp_pred_LLP & the third
%                  column contains the CLP prediction, yp_pred_CLP.
%
% (c) 2017 S. R. Dhruba
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pPred = LatentPredTransLearn(pTrain, sTrain, pTest, sTest)
% Inputs...
X11 = pTrain(:, 1:end-1);           X12 = pTest;                         % Primary X
X21 = sTrain(:, 1:end-1);           X22 = sTest(:, 1:end-1);        % Secondary X
y11 = pTrain(:, end);                                                             % Primary y
y21 = sTrain(:, end);                 y22 = sTest(:, end);              % Secondary y
X1 = [X11; X12];                        X2 = [X21; X22];
[n, p] = size(X1);
nTrain = size(X11, 1);                nTest = size(X12, 1);

% Cost optimization...
A = [ ];    b = [ ];    Aeq = [0, 1, 1];      beq = 1;
lb = [-1, 0, 0]';   ub = [1, 1, 1]';
options = optimoptions(@fmincon, 'display', 'off');

% Response... y11 & y21
fprintf('>> Cost optimization of response data... ')
J = @(c) OptFuncTL(y11, y21, c);        c0 = [-1, 0, 1]';
c_opt = fmincon(J, c0, A, b, Aeq, beq, lb, ub, [ ], options);
[~, w1, a1_opt, a2_opt] = OptFuncTL(y11, y21, c_opt);                 % w = c0 + c1*y1 + c2*y2
fprintf('Finished! \n')

% Predictors... X1 & X2
fprintf('>> Cost optimization of input data... ')
L_opt = zeros(3, p);                                               % Λ = [λ1, λ2, ... , λp]
V = zeros(n, p);
for k = 1:p
    x1k = X1(:, k);     x2k = X2(:, k);                          % x_{1k}, x_{2k}
    Jk = @(Lk) OptFuncTL(x1k, x2k, Lk);         Lk0 = [-1, 0, 1]';
    Lk_opt = fmincon(Jk, Lk0, A, b, Aeq, beq, lb, ub, [ ], options);       % λ_k
    L_opt(:, k) = Lk_opt;
    [~, vk] = OptFuncTL(x1k, x2k, Lk_opt);                % v_k = λk_0 + λk_1 x_{1k} + λk_2 x_{2k} 
    V(:, k) = vk;
end
V1 = V(1:nTrain, :);     V2 = V(nTrain+1:end, :);
fprintf('Finished! \n')

% LATENT REGRESSION PREDICTION (LRP)...
fprintf('>> Performing LRP... \n')
Y21 = [ones(nTrain, 1), y21];
Y22 = [ones(nTest, 1), y22];
b2 = pinv(Y21)*w1;                                              % w1 = Y21*b2
w2_LRP = Y22*b2;                                              % w2 = Y22*b2
W2_LRP = [ones(nTest, 1), w2_LRP];
y12_LRP = W2_LRP*a1_opt;                                 % Final prediction => y12 = a10 + a11*w2
y12_LRP(y12_LRP < 0) = 0;

% LATENT-LATENT PREDICTION (LLP)...
fprintf('>> Performing LLP... \n')
n_tree = 100;       Train_X = V1;      Train_y = w1;        Test_X = V2;
M = TreeBagger(n_tree, Train_X, Train_y, 'method', 'regression');
Pred_y = predict(M, Test_X);
w2_LLP = Pred_y;                                                % Predicted latent vector
W2_LLP = [ones(nTest, 1), w2_LLP];
y12_LLP = W2_LLP*a1_opt;                                  % Final prediction => y12 = a10 + a11*w2
y12_LLP(y12_LLP < 0) = 0;

% COMBINED LATENT PREDICTION (CLP)...
fprintf('>> Performing CLP... \n')
w2_CLP = (w2_LLP + w2_LRP)/2;
W2_CLP = [ones(size(w2_CLP)), w2_CLP];
y12_CLP = W2_CLP*a1_opt;                                 % Final prediction => y12 = a10 + a11*w2
y12_CLP(y12_CLP < 0) = 0;

% Output...
pPred = [y12_LRP, y12_LLP, y12_CLP];
fprintf('>>\n >> Finished!!! \n')
end


% Cost Function...
function [J, u, a1, a2] = OptFuncTL(v1, v2, w)
u = [ones(size(v1)), v1, v2]*w;                               % u = w0 + w1*v1 + w2*v2
U = [ones(size(u)), u];
a1 = pinv(U)*v1;
a2 = pinv(U)*v2;

v1_hat = U*a1;      v1_hat(v1_hat < 0) = 0;
v2_hat = U*a2;     v2_hat(v2_hat < 0) = 0;
J = (norm(v1 - v1_hat)^2 + norm(v2 - v2_hat)^2) / (corr(v1, u) + corr(v2, u));
end




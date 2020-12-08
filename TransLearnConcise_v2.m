clc;    clear;      close all

% CCLE...
[~, ~, gxdata2] = xlsread('Data\CCLE_gene_expression_Oct_13_SRD.xlsx');
data22 = cell2mat(gxdata2(2:end, 2:end));
gn2 = gxdata2(2:end, 1);
CL2 = gxdata2(1, 2:end)';
for k = 1:length(CL2)
    if isnumeric(CL2{k})
        CL2{k} = num2str(CL2{k});
    end
end

% GDSC v6...
[~, ~, gxdata3] = xlsread('Data\GDSC_v6_gene_expression_Oct_13_SRD.xlsx'); 
data33 = cell2mat(gxdata3(2:end, 2:end));
gn3 = gxdata3(2:end, 1);
CL3 = gxdata3(1, 2:end)';
for k = 1:length(CL3)
    if isnumeric(CL3{k})
        CL3{k} = num2str(CL3{k});
    end
end

% % COMMON SET
% Common genes...
cgn = [ ];
for i = 1:length(gn3)
    j = find(strcmpi(gn3(i), gn2));
    if ~isempty(j)
        cgn = [cgn; [i, j]];
    end
end

% Common CLs...
ccl = [ ];
for i = 1:length(CL3)
    j = find(strcmpi(CL3(i), CL2));
    if ~isempty(j)
        ccl = [ccl; [i*ones(size(j)), j]];
    end
end

% GX data for common genes & CLs...
ggn = gn3(cgn(:, 1));
gcl = CL3(ccl(:, 1));
gdata3 = data33(cgn(:, 1), ccl(:, 1))';
gdata2 = data22(cgn(:, 2), ccl(:, 2))';

% median(diag(corr(gdata3', gdata2', 'type', 'pearson')))          % GX correlation
% median(diag(corr(gdata3', gdata2', 'type', 'spearman')))

% % DRUG DATA
% GDSC drug AUC data...
[~, ~, drugdata3] = xlsread('GDSC_v6_Common_Drugs_AUC.xlsx');
DRUGDATA3 = 1 - cell2mat(drugdata3(2:end, 2:end));
DDATA3_CL = drugdata3(2:end, 1);
DCL3 = cell(size(DDATA3_CL));
for k = 1:length(DDATA3_CL)
    if ~isnumeric(DDATA3_CL{k})
        C = strsplit(DDATA3_CL{k}, '-');
        DCL3{k} = strcat(C{1:end});
    else
        DCL3{k} = num2str(DDATA3_CL{k});
    end
end
dList = drugdata3(1, 2:end)';

%% DRUGWISE CALCULATION
kd = 1                                                                 % Drug no
DDATA3 = DRUGDATA3(:, kd);

% CCLE drug AUC data...
sheet = ['Sheet', num2str(kd)];
[~, ~, drugdata2] = xlsread('CCLE_Drug_Sensitivity_Raziur_Extended.xlsx', sheet);
DDATA2 = cell2mat(drugdata2(2:end, 13))/8;      % AUC in column 13
DDATA2_CL = drugdata2(2:end, 1);
DCL2 = cell(size(DDATA2_CL));
for k = 1:length(DDATA2_CL)
    C = strsplit(DDATA2_CL{k}, '_');
    DCL2(k) = C(1);
end

% Common drug CLs...
dccl = [ ];
for i = 1:length(DCL3)
    j = find(strcmpi(DCL3(i), DCL2));
    if ~isempty(j)
        dccl = [dccl; [i*ones(size(j)), j]];
    end
end
dcl = DCL3(dccl(:, 1));
ddata3 = DDATA3(dccl(:, 1));
ddata2 = DDATA2(dccl(:, 2));

% Common GX & DR CLs...
gdcl = [ ];
for j = 1:length(dcl)
    i = find(strcmpi(dcl(j), gcl));
    if ~isempty(i)
        gdcl = [gdcl; [i, j*ones(size(i))]];
    end
end

% Common gx & dr data...
fdcl = gcl(gdcl(:, 1));
gdata33 = gdata3(gdcl(:, 1), :);        ddata33 = ddata3(gdcl(:, 2));
gdata22 = gdata2(gdcl(:, 1), :);        ddata22 = ddata2(gdcl(:, 2));

ind3 = find(~isnan(ddata33));
ind2 = find(~isnan(ddata22));
cind = [ ];
for i = 1:length(ind3)
    I = find(ind3(i) == ind2);
    cind = [cind; ind2(I)];
end
fCL = fdcl(cind);
nCL = length(cind);

% Available response CLs...
gx3 = gdata33(cind, :);         dr3 = ddata33(cind);
gx2 = gdata22(cind, :);         dr2 = ddata22(cind);

% Per column normalization...
ngx3 = gx3 - (ones(nCL, 1)*min(gx3, [ ], 1));   ngx3 = ngx3 ./ (ones(nCL, 1)*max(gx3, [ ], 1));
ngx2 = gx2 - (ones(nCL, 1)*min(gx2, [ ], 1));	ngx2 = ngx2 ./ (ones(nCL, 1)*max(gx2, [ ], 1));
% ngx3 = zscore(gx3);
% ngx2 = zscore(gx2);

% corr(dr3, dr2, 'type', 'pearson')                      % DR correlation
% corr(dr3, dr2, 'type', 'spearman')


% % FEATURE SELECTION
% gx3, dr3, gx2, dr2
K = 10;
rank3 = relieff(gx3, dr3, K)';
rank2 = relieff(gx2, dr2, K)';

Nx = 200;                                                            % Chosen features
IDX3 = rank3(1:Nx);
IDX2 = rank2(1:Nx);
for k = 1:Nx
    ind = find(IDX3(k) == IDX2);
    IDX2(ind) = [ ];
end
Nc = length(IDX3) - length(IDX2);                       % Common features
p = 2*Nx - Nc;                                                     % Total features
IDX = sort([IDX3; IDX2]);
Feat = [Nx, Nc, p]

% Feature selected sets...
X1 = gx3(:, IDX);       X2 = gx2(:, IDX);
% X1 = ngx3(:, IDX);     X2 = ngx2(:, IDX);              % Normalized data
y1 = dr3;                    y2 = dr2;

% % PLOT RESPONSE DISTRIBUTIONS...
idx = [1, 3, 6, 7, 9, 10, 13; 1:3, 5:8];
nBin = 50;

figure(11),     subplot(2,4, idx(2, idx(1, :) == kd))
histogram(y1, nBin, 'normalization', 'probability', 'FaceColor', 'm'),      hold on
histogram(y2, nBin, 'normalization', 'probability', 'FaceColor', 'y'),      hold off
title(dList{kd}, 'FontName', 'Book Antiqua', 'FontWeight', 'Bold', 'FontSize', 13)
xlabel('\bf\itx', 'FontName', 'Book Antiqua', 'FontSize', 11)
ylabel('\bf\itf(x)', 'FontName', 'Book Antiqua', 'FontSize', 11)

% f1 = histcounts(y1, nBin)'/numel(y1);
% v1 = linspace(min(y1), max(y1), nBin)';

% f2 = histcounts(y2, nBin)' / numel(y2);
% v2 = linspace(min(y2), max(y2), nBin)';

v = linspace(0, 1, nBin);
f1 = ksdensity(y1, v);      % f1 = f1/max(f1);
f2 = ksdensity(y2, v);     % f2 = f2/max(f2);

% figure
% plot(v, f1, 'r-.', v, f2, 'b-.', 'linewidth', 3)
% xlabel('$x$', 'interpreter', 'latex'),    ylabel('$f(x)$', 'interpreter', 'latex')
% title(['AUC Distributions for ', dList{kd}], 'interpreter', 'latex')
% hold on,
% histogram(y1, nBin, 'normalization', 'pdf'),    hold on
% histogram(y2, nBin, 'normalization', 'pdf')
% legend('\bfGDSC', '\bfCCLE')

% figure,
% histfit(y1, nBin, 'kernel'),    hold on
% histfit(y2, nBin, 'kernel')
% xlim([0, 1])

%% MAPPED PREDICTION - K-FOLD CROSS-VALIDATION...
% 1 = GDSC, 2 = CCLE
% y1x = 0.65*y1 + 0.35;
y1x = y1;

rng(0)
Ntimes = 1;
Metric = zeros(3, 4, Ntimes);
Rfold = zeros(Ntimes, 2);
for i = 1:Ntimes
    % rng(7)
    n_tree = 100;
    nSamp = 50;     nFold = round(nCL / nSamp);
    Idx = crossvalind('kfold', nCL, nFold);
    rFold = zeros(nFold, 2);      nSets = zeros(nFold, 2);
    r1 = zeros(nFold, 4);           rs1 = zeros(nFold, 4);     err1 = zeros(nFold, 4);
    tic
    % figure
    for k = 1:nFold
        % fprintf('Fold %d\n', k)
        idxTrain = (Idx == k);           idxTest = (Idx ~= k);         % Fold data Idx's
        X11 = X1(idxTrain, :);            X12 = X1(idxTest, :);
        X21 = X2(idxTrain, :);           X22 = X2(idxTest, :);
        y11 = y1x(idxTrain);                y12 = y1x(idxTest);
        y21 = y2(idxTrain);               y22 = y2(idxTest);
        
        rFold(k, :) = [corr(y11, y21, 'type', 'pearson'), corr(y11, y21, 'type', 'spearman')];
        
        % GDSC Mapped prediction...
        pTrain1 = [X11, y11];        pTest1 = X12;
        sTrain1 = [X21, y21];       sTest1 = [X22, y22];
        pPred1 = MappingTransLearn(pTrain1, sTrain1, pTest1, sTest1, n_tree);
        y12mb = pPred1;
        
        % Direct prediction...
        rng(0)                                          % RF seed for reproducibility & order independence
        Train_X1 = X11;      Train_y1 = y11;        Test_X1 = X12;
        model_1 = TreeBagger(n_tree, Train_X1, Train_y1, 'method', 'regression');
        Predict_y1 = predict(model_1, Test_X1);
        y12p = Predict_y1;
        
        % CCLE model prediction...
        rng(0)                                          % RF seed for reproducibility & order independence
        Test_X2 = X12;      Train_X2 = [X21; X22];      Train_y2 = [y21; y22];
        model_2 = TreeBagger(n_tree, Train_X2, Train_y2, 'method', 'regression');
        Predict_y2 = predict(model_2, Test_X2);
        y12c = Predict_y2;
        
        % Combined model prediction...
        rng(0)                                          % RF seed for reproducibility & order independence
        Test_X3 = X12;      Train_X3 = [X11; X22];      Train_y3 = [y11; y22];
        model_3 = TreeBagger(n_tree, Train_X3, Train_y3, 'method', 'regression');
        Predict_y3 = predict(model_3, Test_X3);
        y12cg = Predict_y3;
        
        nSets(k, :) = [numel(y11), numel(y12)];
        r1(k, :) = [corr(y12, y12mb, 'type', 'pearson'),  corr(y12, y12c, 'type', 'pearson'),...
                            corr(y12, y12cg, 'type', 'pearson'),  corr(y12, y12p, 'type', 'pearson')];
        rs1(k, :) = [corr(y12, y12mb, 'type', 'spearman'),  corr(y12, y12c, 'type', 'spearman'),...
                            corr(y12, y12cg, 'type', 'spearman'),  corr(y12, y12p, 'type', 'spearman')];
        err1(k, :) = sqrt([mean((y12 - y12mb).^2),  mean((y12 - y12c).^2),...
                                    mean((y12 - y12cg).^2),  mean((y12 - y12p).^2)]) / range(y12);
        
        % 	subplot(3, 2, k)
        %     scatter(y12, y12mb),    xlabel('observed'),     ylabel('predicted'),    hold on
        %     scatter(y12, y12cg),     xlabel('observed'),     ylabel('predicted'),    hold off;   box on
        %     legend('MP', 'CMP'),    title(['Fold ', num2str(k)])
    end
    toc
    
    Rfold(i, :) = mean(rFold, 1);
    Metric(:, :, i) = [mean(r1, 1); mean(rs1, 1); mean(err1, 1)];
end

% % Table...
fprintf('K-Fold Cross-validation Averaged over %d times...\n', Ntimes)
fprintf('\tK = %d, Mean nTrain = %0.2f, Mean nTest = %0.2f\n', nFold, mean(nSets, 1))
disp([corr(y1, y2, 'type', 'pearson'), corr(y1, y2, 'type', 'spearman'); mean(Rfold, 1)])
METRIC = mean(Metric, 3);
TBL = array2table(METRIC, 'variablenames', {'MP', 'CP', 'CMP', 'DP'},...
                                                'rownames', {'PCORR'; 'SCORR'; 'NRMSE'});
disp(TBL)


%% HISTOGRAM EQUALIZATION...
nBin = 100;
% y1
f1 = histcounts(y1, nBin)'/numel(y1);
F1 = cumsum(f1);
v1 = linspace(min(y1), max(y1), nBin)';

z2 = rand(size(y1));
f2 = histcounts(z2, nBin)' / numel(z2);
F2 = cumsum(f2);
v2 = linspace(min(z2), max(z2), nBin)';

Match = zeros(nBin, 1);
for k = 1:nBin
    [~, ind] = min(abs(F1(k) - F2));
    Match(k) = ind;
end
y1h = Match(round(y1 * (nBin - 1)) + 1) / (nBin - 1);
f1h = histcounts(y1h, nBin)' / numel(y1h);

figure
subplot(311),    histogram(y1, nBin, 'normalization', 'probability')
xlim([0, 1]),      legend('GDSC'),      title(['Drug ', num2str(kd)])
subplot(312),    histogram(y1h, nBin, 'normalization', 'probability')
xlim([0, 1]),      legend('GDSC-HE')
subplot(313),    histogram(y2, nBin, 'normalization', 'probability')
xlim([0, 1]),      legend('CCLE')

figure,     ksdensity(y1, linspace(0, 1, 100)),   xlim([0, 1])
hold on,    ksdensity(y1h, linspace(0, 1, 100)),   xlim([0, 1])
hold on,    ksdensity(y2, linspace(0, 1, 100)),   xlim([0, 1])
legend('GDSC', 'GDSC-HE', 'CCLE'),   title(['Drug ', num2str(kd)])



%% PLOTTING...
figure,     
subplot(121),   ksdensity(y12),      hold on,    ksdensity(y12mb)
hold on,           ksdensity(y12cg),  hold on,    ksdensity(y12p)
hold on,           ksdensity(y12c)
title('KS Density of y'),    legend('Actual', 'MP', 'CMP', 'DP', 'CP')

nBin = 20;
subplot(122),   histogram(y12, nBin),      hold on,    histogram(y12mb, nBin)
hold on,           histogram(y12cg, nBin),   hold on,    histogram(y12p, nBin)
hold on,           histogram(y12c, nBin)
title('Histogram of y'),    legend('Actual', 'MP', 'CMP', 'DP', 'CP')

VarFeat = [var(X2)', var(Train_X2)', var(Train_X3)', var(Train_X1)'];
VarTable = array2table([VarFeat; mean(VarFeat, 1)],...
                                        'rownames', [strtrim(cellstr(num2str((1:p)'))); 'Mean'],...
                                        'variablenames', {'MP', 'CP', 'CMP', 'DP'});

% Covariate dist...
covIdx = 1;
figure,     
subplot(121),   ksdensity(X2(:, covIdx)),              hold on,    ksdensity(Train_X2(:, covIdx))
hold on,           ksdensity(Train_X3(:, covIdx)),      hold on,    ksdensity(Train_X1(:, covIdx))
title(['KS Density of X(:, ', num2str(covIdx), ')']),    legend('MP', 'CP', 'CMP', 'DP')

nBin = 20; 
subplot(122),   histogram(X2(:, covIdx), nBin, 'facecolor', 'b')
hold on,    histogram(Train_X2(:, covIdx), nBin, 'facecolor', 'r')
hold on,	histogram(Train_X3(:, covIdx), nBin, 'facecolor', 'k')
hold on,    histogram(Train_X1(:, covIdx), nBin, 'facecolor', 'g')
title(['Histogram of X(:, ', num2str(covIdx), ')']),    legend('MP', 'CP', 'CMP', 'DP')


%% LATENT PREDICTION - K-FOLD CROSS-VALIDATION...
% 1 = GDSC, 2 = CCLE
rng(14)
Ntimes = 1;
Rfold = zeros(Ntimes, 2);
Metric1 = zeros(2, 4, Ntimes);      Metric2 = zeros(2, 4, Ntimes);
for i = 1:Ntimes
    % rng(7)
    n_tree = 100;
    nSamp = 50;     nFold = round(nCL / nSamp);
    Idx = crossvalind('kfold', nCL, nFold);
    rFold = zeros(nFold, 2);
    nSets1 = zeros(nFold, 2);	nSets2 = zeros(nFold, 2);
    r1 = zeros(nFold, 4);          err1 = zeros(nFold, 4);
    r2 = zeros(nFold, 4);          err2 = zeros(nFold, 4);
    tic
    % figure
    for k = 1:nFold
        % fprintf('Fold %d\n', k)
        idxTrain = (Idx == k);           idxTest = (Idx ~= k);         % Fold data Idx's
        X11 = X1(idxTrain, :);            X12 = X1(idxTest, :);
        X21 = X2(idxTrain, :);           X22 = X2(idxTest, :);
        y11 = y1(idxTrain);                y12 = y1(idxTest);
        y21 = y2(idxTrain);               y22 = y2(idxTest);
        
        rFold(k, :) = [corr(y11, y21, 'type', 'pearson'), corr(y11, y21, 'type', 'spearman')];
                
        % GDSC Latent prediction...
        pTrain1 = [X11, y11];        pTest1 = X12;
        sTrain1 = [X21, y21];       sTest1 = [X22, y22];
        pPred1 = LatentPredTransLearn(pTrain1, sTrain1, pTest1, sTest1);        % LP Methods
        y12w = pPred1(:, 1);         y12L = pPred1(:, 2);         y12e = pPred1(:, 3);
        
        % CCLE Latent prediction...
        pTrain2 = [X21, y21];        pTest2 = X22;
        sTrain2 = [X11, y11];         sTest2 = [X12, y12];
        pPred2 = LatentPredTransLearn(pTrain2, sTrain2, pTest2, sTest2);        % LP Methods
        y22w = pPred2(:, 1);         y22L = pPred2(:, 2);         y22e = pPred2(:, 3);
         
        % GDSC Direct prediction...
        rng(0)
        Train_X1 = X11;      Train_y1 = y11;        Test_X1 = X12;
        mdl_1 = TreeBagger(n_tree, Train_X1, Train_y1, 'method', 'regression');
        Predict_y1 = predict(mdl_1, Test_X1);
        y12p = Predict_y1;        y12p(y12p < 0) = 0;
        
        % CCLE Direct prediction...
        rng(0)
        Train_X2 = X21;      Train_y2 = y21;        Test_X2 = X22;
        mdl_2 = TreeBagger(n_tree, Train_X2, Train_y2, 'method', 'regression');
        Predict_y2 = predict(mdl_2, Test_X2);
        y22p = Predict_y2;        y22p(y22p < 0) = 0;
        
        % Metrics...
        nSets1(k, :) = [numel(y11), numel(y12)];
        r1(k, :) = [corr(y12, y12w), corr(y12, y12L), corr(y12, y12e), corr(y12, y12p)];
        err1(k, :) = sqrt([mean((y12 - y12w).^2),  mean((y12 - y12L).^2),...
                                    mean((y12 - y12e).^2),  mean((y12 - y12p).^2)]) / range(y12);
        
        nSets2(k, :) = [numel(y21), numel(y22)];
        r2(k, :) = [corr(y22, y22w), corr(y22, y22L), corr(y22, y22e), corr(y22, y22p)];
        err2(k, :) = sqrt([mean((y22 - y22w).^2),  mean((y22 - y22L).^2),...
                                    mean((y22 - y22e).^2),  mean((y22 - y22p).^2)]) / range(y22);
    end
    toc
    
    Rfold(i, :) = mean(rFold, 1);
    Metric1(:, :, i) = [mean(r1, 1); mean(err1, 1)];
    Metric2(:, :, i) = [mean(r2, 1); mean(err2, 1)];
end

% % Table...
fprintf('GDSC latent predictions using CCLE...\n')
fprintf('K-Fold Cross-validation Averaged over %d times...\n', Ntimes)
fprintf('\tK = %d, Mean nTrain = %0.2f, Mean nTest = %0.2f\n', nFold, mean(nSets1, 1))
disp([corr(y1, y2, 'type', 'pearson'), corr(y1, y2, 'type', 'spearman'); mean(Rfold, 1)])
METRIC1 = mean(Metric1, 3);
TBL1 = array2table(METRIC1, 'variablenames', {'LRP', 'LLP', 'CLP', 'DP'},...
                                                'rownames', {'PCORR'; 'NRMSE'});
disp(TBL1)

fprintf('CCLE latent predictions using GDSC...\n')
fprintf('K-Fold Cross-validation Averaged over %d times...\n', Ntimes)
fprintf('\tK = %d, Mean nTrain = %0.2f, Mean nTest = %0.2f\n', nFold, mean(nSets1, 1))
disp([corr(y1, y2, 'type', 'pearson'), corr(y1, y2, 'type', 'spearman'); mean(Rfold, 1)])
METRIC2 = mean(Metric2, 3);
TBL2 = array2table(METRIC2, 'variablenames', {'LRP', 'LLP', 'CLP', 'DP'},...
                                                'rownames', {'PCORR'; 'NRMSE'});
disp(TBL2)




%% RANDOM HOLDOUT SET
% y1x = 0.65*y1 + 0.35;
% rng(0)
nSamp = 50;                                                         % Known smaller subset
ind1 = sort(randperm(nCL, nSamp))';
ind2 = (1:nCL)';            ind2(ind1) = [ ];
X11 = X1(ind1, :);          X12 = X1(ind2, :);
X21 = X2(ind1, :);         X22 = X2(ind2, :);
y11 = y1(ind1);              y12 = y1(ind2);            Y1 = [ones(size(y1)), y1];
y21 = y2(ind1);             y22 = y2(ind2);           Y2 = [ones(size(y2)), y2];
[corr(y1, y2, 'type', 'pearson'), corr(y1, y2, 'type', 'spearman');...
corr(y11, y21, 'type', 'pearson'), corr(y11, y21, 'type', 'spearman')]

%% PREDICTION USING MP FUNCTION...
% GDSC...
pTrain1 = [X11, y11];        pTest1 = X12;
sTrain1 = [X21, y21];       sTest1 = [X22, y22];
pPred1 = MappingTransLearn(pTrain1, sTrain1, pTest1, sTest1);            % MP Method
y12mb = pPred1;

% % CCLE...
% pTrain2 = [X21, y21];        pTest2 = X22;
% sTrain2 = [X11, y11];         sTest2 = [X12, y12];
% pPred2 = MappingTransLearn(pTrain2, sTrain2, pTest2, sTest2);            % MP Method
% y22mb = pPred2;

% COMPARISON TABLE...
n_tree = 100;

tic
% CCLE model prediction...
% Using original data matrix...
rng(0)                  % RF seed for reproducibility & order independence
Test_X2 = X12;      Train_X2 = [X21; X22];      Train_y2 = [y21; y22];
model_2 = TreeBagger(n_tree, Train_X2, Train_y2, 'method', 'regression');
Predict_y2 = predict(model_2, Test_X2);
y12c = Predict_y2;
toc

tic
% Combined model prediction...
rng(0)                  % RF seed for reproducibility & order independence
Test_X3 = X12;      Train_X3 = [X11; X22];      Train_y3 = [y11; y22];
model_3 = TreeBagger(n_tree, Train_X3, Train_y3, 'method', 'regression');
Predict_y3 = predict(model_3, Test_X3);
y12cg = Predict_y3;
toc

tic
% Direct prediction...
% Trained w/ GDSC only...
rng(0)                  % RF seed for reproducibility & order independence
Train_X1 = X11;      Train_y1 = y11;        Test_X1 = X12;
model_1 = TreeBagger(n_tree, Train_X1, Train_y1, 'method', 'regression');
Predict_y1 = predict(model_1, Test_X1);
y12p = Predict_y1;
toc

% % Table...
METRIC = [corr(y12, y12mb, 'type', 'pearson'), corr(y12, y12c, 'type', 'pearson'),...
                  corr(y12, y12cg, 'type', 'pearson'), corr(y12, y12p, 'type', 'pearson');...
                  corr(y12, y12mb, 'type', 'spearman'), corr(y12, y12c, 'type', 'spearman'),...
                  corr(y12, y12cg, 'type', 'spearman'), corr(y12, y12p, 'type', 'spearman');...
                  sqrt(mse(y12, y12mb))/range(y12), sqrt(mse(y12, y12c))/range(y12),...
                  sqrt(mse(y12, y12cg))/range(y12), sqrt(mse(y12, y12p))/range(y12);...
                  mean(y12mb), mean(y12c), mean(y12cg), mean(y12p)];
TBL = array2table(METRIC', 'variablenames', {'r'; 'rs'; 'NRMSE'; 'MRV'},...
                                                'rownames', {'MP', 'CP', 'CMP', 'DP'});
disp(TBL)
fprintf('Actual MRV = %f\n', mean(y12))


%%
g = 50;
figure,   ksdensity(X1(:, g))
hold on,  ksdensity(X2(:, g))
hold on,  ksdensity(Train_X3(:, g))
legend('GDSC', 'CCLE', 'Appended')
title(['Gene ', num2str(g)])


%% PREDICTION USING LP FUNCTION...
% GDSC...
pTrain1 = [X11, y11];        pTest1 = X12;
sTrain1 = [X21, y21];       sTest1 = [X22, y22];
pPred1 = LatentPredTransLearn(pTrain1, sTrain1, pTest1, sTest1);        % LP Methods
y12w = pPred1(:, 1);         y12L = pPred1(:, 2);         y12e = pPred1(:, 3);

% CCLE...
pTrain2 = [X21, y21];        pTest2 = X22;
sTrain2 = [X11, y11];         sTest2 = [X12, y12];
pPred2 = LatentPredTransLearn(pTrain2, sTrain2, pTest2, sTest2);        % LP Methods
y22w = pPred2(:, 1);         y22L = pPred2(:, 2);         y22e = pPred2(:, 3);

% LP METHODS COMPARISON TABLES...
n_tree = 100;
tic
% Direct prediction...
% GDSC...
rng(0)
Train_X1 = X11;      Train_y1 = y11;        Test_X1 = X12;
mdl_1 = TreeBagger(n_tree, Train_X1, Train_y1, 'method', 'regression');
Predict_y1 = predict(mdl_1, Test_X1);
y12p = Predict_y1;        y12p(y12p < 0) = 0;
toc

tic
% CCLE...
rng(0)
Train_X2 = X21;      Train_y2 = y21;        Test_X2 = X22;
mdl_2 = TreeBagger(n_tree, Train_X2, Train_y2, 'method', 'regression');
Predict_y2 = predict(mdl_2, Test_X2);
y22p = Predict_y2;        y22p(y22p < 0) = 0;
toc

% c_opt
METRIC1 = [corr(y12, y12w), corr(y12, y12L), corr(y12, y12e), corr(y12, y12p);
                    sqrt(mse(y12, y12w))/range(y12), sqrt(mse(y12, y12L))/range(y12),... 
                    sqrt(mse(y12, y12e))/range(y12), sqrt(mse(y12, y12p))/range(y12)];
    
TBL_GDSC = array2table(METRIC1, 'rownames', {'CORR'; 'NRMSE';}, ...
                                                        'variablenames', {'LRP', 'LLP', 'CLP', 'DP'})


METRIC2 = [corr(y22, y22w), corr(y22, y22L), corr(y22, y22e), corr(y22, y22p);
                    sqrt(mse(y22, y22w))/range(y22), sqrt(mse(y22, y22L))/range(y22),... 
                    sqrt(mse(y22, y22e))/range(y22), sqrt(mse(y22, y22p))/range(y22)];
    
TBL_CCLE = array2table(METRIC2, 'rownames', {'CORR'; 'NRMSE';}, ...
                                                        'variablenames', {'LRP', 'LLP', 'CLP', 'DP'})




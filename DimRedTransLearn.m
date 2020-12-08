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


%% FEATURE SELECTION
K = 10;
rank3 = relieff(gx3, dr3, K)';
rank2 = relieff(gx2, dr2, K)';

Nx = 500;                                                            % Chosen features
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

xgx3 = ngx3(:, IDX);
xgx2 = ngx2(:, IDX);


%% PCA RECONSTRUCTION
N = 150;
% VAR = 0.9;
[e22, e33, e32] = PCAerror2(gx3, gx2, N);          % CCLE -> GDSC
[~, ~, e322] = PCAerror2(gx3, gx2, N, 1);            % (n-1)*p
result = [e22, e33, e32, e322]

% PCA after normalization...
[ne22, ne33, ne32] = PCAerror2(ngx3, ngx2, N);           % CCLE -> GDSC
[nE33, nE22, nE23] = PCAerror2(ngx2, ngx3, N);          % GDSC -> CCLE
result = [ne22, ne33, ne32]
RESULT = [nE33, nE22, nE23]

% PCA after feature selection..
VAR = 0.95;
[xe22, xe33, xe32, N2] = PCAerror2(xgx3, xgx2, VAR);           % CCLE -> GDSC
[xE33, xE22, xE23, N3] = PCAerror2(xgx2, xgx3, VAR);          % GDSC -> CCLE
result = [N2, xe22, xe33, xe32]
RESULT = [N3, xE33, xE22, xE23]

%% CCLE INTO 2 SETS
nSet = 20;                                                            % No of picked random subsets
rat = 0.85;
result2 = zeros(nSet, 4);
result3 = zeros(nSet, 4);
for k = 1:nSet
    % rng(0)
    L = round(nCL * rat);                                        % Random 'L' cell-lines
    is1 = sort(randperm(size(xgx2, 1), L))';
    is2 = (1:size(xgx2, 1))';       is2(is1) = [ ];
    
    cset1 = xgx2(is1, :);                                         % Training set => CCLE
    gset1 = xgx3(is1, :);                                         % Training set => GDSC
    cset2 = xgx2(is2, :);                                        % Test set => CCLE
    gset2 = xgx3(is2, :);                                        % Test set => GDSC
    
    VAR = 0.98;
    [pe11, pe22, pe21, N] = PCAerror2(cset2, cset1, VAR);       % CCLE reconstruction
    result2(k, :) = [N, pe11, pe22, pe21];
    [pe11, pe22, pe21, N] = PCAerror2(gset2, cset1, VAR);       % GDSC reconstruction
    result3(k, :) = [N, pe11, pe22, pe21];
end
cerr = mean(result2, 1)
gerr = mean(result3, 1)

%% WHOLE GX DATA
% gn3, CL3, data33', gn2, CL2, data22'
% DCL3, DDATA3, DCL2, DDATA2
gDATA3 = data33';
gDATA2 = data22';

% GDSC...
ccl3 = [ ];
for i = 1:length(CL3)
    j = find(strcmpi(CL3(i), DCL3));
    if ~isempty(j)
        ccl3 = [ccl3; [i*ones(size(j)), j]];
    end
end
CL33 = CL3(ccl3(:, 1));
gDATA33 = gDATA3(ccl3(:, 1), :);
dDATA33 = DDATA3(ccl3(:, 2));

% CCLE...
ccl2 = [ ];
for i = 1:length(CL2)
    j = find(strcmpi(CL2(i), DCL2));
    if ~isempty(j)
        ccl2 = [ccl2; [i*ones(size(j)), j]];
    end
end
CL22 = CL2(ccl2(:, 1));
gDATA22 = gDATA2(ccl2(:, 1), :);
dDATA22 = DDATA2(ccl2(:, 2));

% Available response CLs...
IND2 = find(~isnan(dDATA22));
IND3 = find(~isnan(dDATA33));
GX3 = gDATA33(IND3, :);         DR3 = dDATA33(IND3);
GX2 = gDATA22(IND2, :);         DR2 = dDATA22(IND2);

NGX3 = GX3 ./ (ones(size(GX3, 1), 1)*max(GX3, [ ], 1));               % Normalized column
NGX2 = GX2 ./ (ones(size(GX2, 1), 1)*max(GX2, [ ], 1));    

% Feature selection...
K = 10;
RANK3 = relieff(GX3, DR3, K)';
RANK2 = relieff(GX2, DR2, K)';

%% CCLE INTO 2 SETS
nSet = 20;                                                           % No of picked random subsets
RESULT = zeros(nSet, 4);
for k = 1:nSet
    % rng(0)
    L = 350;                                                          % Random 'L' cell-lines
    is1 = sort(randperm(size(NGX2, 1), L))';
    is2 = (1:size(NGX2, 1))';       is2(is1) = [ ];
    
    % P = size(NGX2, 2);
    P = 1000;                                                          % No of features
    gset1 = NGX2(is1, RANK2(1:P));
    gset2 = NGX2(is2, RANK2(1:P));
    
    VAR = 0.98;
    [PE11, PE22, PE21, N] = PCAerror2(gset2, gset1, VAR);
    RESULT(k, :) = [N, PE11, PE22, PE21];
end
mean(RESULT, 1)


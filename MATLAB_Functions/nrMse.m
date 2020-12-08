function nRMSE = nrMse(Xactual, Xestimate, varargin)
% function err = nrMse(Xactual, Xestimate, factor)
    % Returns normalized root mean square error for 2 vectors or matrices
    
    flag = isvector(Xactual);
    if (size(Xactual) ~= size(Xestimate))
        error('Inputs must have the same size!')
    end
    
    RMSE = sqrt(norm(Xactual - Xestimate, 'fro')^2 / numel(Xactual);
    
    factor = 'range';
    if nargin > 2
        factor = varargin{1};
    end
    
    switch factor
        case 'range'
            NF = @(x) max(x(:)) - min(x(:));
        case 'std'
            switch flag
                case 1
                    NF = @(x) std(x);
                case 0
                    NF = @(x) std2(x);
            end
        case 'mean'
            NF = @(x) mean(x(:));
    end
        
    nRMSE = RMSE / NF(Xactual);
    

end
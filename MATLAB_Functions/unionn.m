function [C, varargout] = unionn(varargin)
% UNIONN performs union of N input sets.
    % NOTATION: function [C, I1, I2,...] = UNIONN(A1, A2,... , setOrder)
    %   Inputs:
    %       A1, A2,... => Input sets
    %       setOrder (optional) => Order of the united set. Use either 'sorted' or 'stable'. 
    %                                           Default output is 'sorted'.
    %   Outputs:
    %       C => United set. Order defined by setOrder input.
    %       I1, I2,... => Indices of the elements contibuted by A1, A2,... in C. Order defined by setOrder input.
    %
    % (c) 2018 S.R.Dhruba
    
    % Check inputs...
    sortOrders = {'sorted', 'stable'};
    if nargin < 2
        error('Insufficient inputs!!! At least two input sets are required.')
    else
        orderFlag = 0;
        if and(ischar(class(varargin{end})), length(varargin{end}) == 6)
            orderFlag = sum(strcmpi(sortOrders, varargin{end}));
        end
        numSets = nargin - orderFlag;                      % #input sets
        inClasses = cellfun(@class, varargin(1 : numSets), 'uniformoutput', 0);
        switch (sum(cell2mat(regexpi(inClasses, inClasses(1)))) == numSets)     % Check if input classes are same
            case 1
                setOrder = varargin{end};
                A = varargin(1 : numSets);
            case 0
                error('All input sets must have the same class!!!')
        end
    end
    
    % Perform union & retrieve indices...
    C = union(A{1}, A{2});
    for k = 3 : numSets
        C = union(C, A{k});                                       % Get the overall union
    end
    C = union(A{1}, C, setOrder);                            % Get the ordered union set
    
    D = C;
    Aidx = cell(1, numSets);
    for k = 1 : numSets
        [~, ~, Aidx{k}] = intersect(D, A{k}, setOrder);  % Get the contributing element indices for each set
        D = setdiff(D, A{k}, 'stable');
    end
    varargout = Aidx;
    
end

